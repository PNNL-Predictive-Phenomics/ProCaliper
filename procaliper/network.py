from __future__ import annotations

import distanceclosure as dc
import networkx as nx
import numpy as np

import procaliper.protein_structure.distance as psd
from procaliper._protein import Protein


def contact_network(protein: Protein, max_dist_angstroms: float = 10.0) -> nx.Graph:
    """Constructs a contact network from a protein.

    Args:
        protein (Protein): Protein object.
        max_dist_angstroms (float, optional): Maximum distance between residues to be considered a contact. Defaults to 10.0.

    Returns:
        nx.Graph: Contact network.
    """
    return nx.from_numpy_array(
        psd.contact_map(protein.get_biopython_structure(), max_dist_angstroms)
    )


def distance_network(protein: Protein, max_dist_angstroms: float = 20) -> nx.Graph:
    """Constructs a distance network from a protein.

    Args:
        protein (Protein): Protein object.
        max_dist_angstroms (float, optional): Maximum distance between residues.
            Values greater than this will be set to np.inf. Defaults to 20.

    Returns:
        nx.Graph: Distance network.
    """
    g = nx.from_numpy_array(
        psd.distance_matrix(protein.get_biopython_structure(), max_dist_angstroms)
    )
    for u, v, d in list(g.edges(data=True)):
        if d["weight"] == np.inf:
            g.remove_edge(u, v)
        else:
            d["proximity"] = 1 / (d["weight"] + 1)
            d["d2"] = d["weight"] ** 2
    return g


def _region_label(site_list_zero_indexed: list[int], sequence: str) -> str:
    """WARNING: assumes sorted (ascending) list as input!"""
    s1 = sequence[site_list_zero_indexed[0]]
    s2 = sequence[site_list_zero_indexed[-1]]
    if len(site_list_zero_indexed) == 1:
        return f"{s1}{site_list_zero_indexed[0] + 1}"
    return f"{s1}{site_list_zero_indexed[0] + 1}..{s2}{site_list_zero_indexed[-1] + 1}"


def _region_type(region_name: str) -> str:
    if region_name.startswith("p_"):
        return "PTM site"
    elif region_name.startswith("b_"):
        return "binding site"
    elif region_name.startswith("a_"):
        return "active site"
    else:
        return "annotated region"


def regulatory_distance_network(protein: Protein) -> nx.Graph:
    """Constructs a regulatory region distance network from a protein.

    Distances are computed between PTM sites, annotated regions, binding sites, and active sites.

    Node labels will be 1-indexed and inclusive (e.g., `"K5..C7"` refers to residues 5, 6, and 7).
    The letter in front of the index refers to the first and last amino acid in the region.

    Args:
        protein (Protein): Protein object.

    Returns:
        nx.Graph: Distance network.
    """

    ptms = {f"p_{i}": [i] for i, x in enumerate(protein.site_annotations.ptm) if x}
    binding = {
        f"b_{i}": [i] for i, x in enumerate(protein.site_annotations.binding) if x
    }
    active = {f"a_{i}": [i] for i, x in enumerate(protein.site_annotations.active) if x}
    regions = protein.site_annotations.regions
    domains = protein.site_annotations.domains

    all_regs = {**ptms, **binding, **active, **regions, **domains}

    protein_residues = protein.get_biopython_residues()
    all_regs_residues = {
        k: [protein_residues[i] for i in v] for k, v in all_regs.items()
    }

    g = nx.Graph()
    for k, v in all_regs.items():
        g.add_node(
            k,
            label=_region_label(v, protein.data["sequence"]),
            region_type=_region_type(k),
            residues=v,
        )

    for k1, v1 in all_regs_residues.items():
        for k2, v2 in all_regs_residues.items():
            if k1 == k2:
                continue
            weight = psd.region_distance(v1, v2)
            g.add_edge(k1, k2, weight=weight, d2=weight**2, proximity=1 / (weight + 1))

    return g


def euclidean_backbone(g: nx.Graph) -> nx.Graph:
    """Returns the Euclidean backbone of a distance network.

    The Euclidean backbone of a weighted graph g is the smallest subgraph of g that contains
    all shortest paths where a path length is determined by the square root of the sum of the
    squared edge weights.

    This is useful for sparsifying a distance network without disconnecting it.

    Args:
        g (nx.Graph): Distance network. Edges must have an attribute "d2" representing the
            squared edge weight. This is computed by `distance_network` and
            `regulatory_distance_network` automatically.

    Returns:
        nx.Graph: Euclidean backbone.
    """
    return dc.backbone(g, weight="d2", kind="metric")
