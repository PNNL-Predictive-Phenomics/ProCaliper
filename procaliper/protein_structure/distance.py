from __future__ import annotations

from itertools import combinations, product
from typing import Any

import numpy as np
import numpy.typing as npt
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

"""Module for computing distances between residues and regions on a protein.
"""


def region_distance(
    region_1: list[Residue], region_2: list[Residue]
) -> np.floating[Any]:
    """Compute the distance between two regions of a protein, in Angstroms.

    Args:
        region_1 (list[Residue]): first region
        region_2 (list[Residue]): second region

    Returns:
        np.floating[Any]: minimum distance between the two regions
    """
    return min(residue_distance(r1, r2) for r1, r2 in product(region_1, region_2))


def residue_distance(
    r1: Residue,
    r2: Residue,
) -> np.floating[Any]:
    """Compute the distance between two residues, in Angstroms.

    Args:
        r1 (Residue): first residue
        r2 (Residue): second residue

    Returns:
        np.floating[Any]: distance between the two residues
    """
    dv = r1["CA"].coord - r2["CA"].coord
    return np.linalg.norm(dv)


def distance_matrix(
    structure: Structure, thresh: float = np.inf
) -> npt.NDArray[np.float64]:
    """Compute a distance matrix for a protein structure.

    Args:
        structure (Structure): protein structure
        thresh (float, optional): threshold for distance. Defaults to np.inf.
            Distances greater than this will be set to np.inf.

    Returns:
        npt.NDArray[np.float64]: distance matrix with shape nxn where n is the
            number of residues in the structure.
    """
    residues = [res for model in structure for chain in model for res in chain]
    residues = list(enumerate(residues))
    adj = np.ones((len(residues), len(residues))) * np.inf

    # a residue has zero distance to itself
    for i in range(len(residues)):
        adj[i, i] = 0

    for (row, r1), (col, r2) in combinations(residues, 2):
        dist = residue_distance(r1, r2)
        if dist <= thresh:
            adj[row, col] = dist
            adj[col, row] = adj[row, col]
    return adj


def proximity_matrix(
    structure: Structure, thresh: float = 0
) -> npt.NDArray[np.float64]:
    """Compute a proximity matrix for a protein structure.

    Args:
        structure (Structure): protein structure
        thresh (float, optional): threshold for proximity. Defaults to 0. Proximity
            less than this will be set to 0.

    Returns:
        npt.NDArray[np.float64]: proximity matrix with shape nxn where n is the
            number of residues in the structure.
    """
    residues = [res for model in structure for chain in model for res in chain]
    residues = list(enumerate(residues))
    adj = np.zeros((len(residues), len(residues)))

    # a residue has proximity 1 to itself
    for i in range(len(residues)):
        adj[i, i] = 1

    for (row, r1), (col, r2) in combinations(residues, 2):
        prox = 1 / (residue_distance(r1, r2) + 1)
        if prox >= thresh:
            adj[row, col] = prox
            adj[col, row] = adj[row, col]
    return adj


def contact_map(
    structure: Structure, max_dist_angsrtom: float = 10
) -> npt.NDArray[np.int8]:
    """A contact map for a protein structure.

    Args:
        structure (Structure): protein structure
        max_dist_angsrtom (float, optional): Largest distance to consider a contact,
            in Angstroms. Defaults to 10.

    Returns:
        npt.NDArray[np.float64]: contact map with shape nxn where n is the
            number of residues in the structure.
    """
    residues = [res for model in structure for chain in model for res in chain]
    residues = list(enumerate(residues))
    adj = np.zeros((len(residues), len(residues)), dtype=np.int8)

    # a residue has zero distance to itself
    for i in range(len(residues)):
        adj[i, i] = np.int8(1)

    for (row, r1), (col, r2) in combinations(residues, 2):
        dist = residue_distance(r1, r2)
        if dist <= max_dist_angsrtom:
            adj[row, col] = np.int8(1)
            adj[col, row] = np.int8(1)
    return adj
