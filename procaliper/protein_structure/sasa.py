from __future__ import annotations

from typing import TypedDict

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

from .hyperparameters import HYPERPARAMETERS

N_POINTS = int(
    HYPERPARAMETERS["n_points"]
)  # Dictates the accuracy of the ShrakeRupley calculation. Higher values are more accurate but slower to calculate. Default = 100


class SASAData(TypedDict):
    """Data class for holding SASA data from computed from a PDB file.

    Attributes:
        all_sasa_value (list[float]): The overall SASA value for each site
            (computed as sum of atom SASA values).
        atom_sasa_values (list[list[float]]): The SASA value for the each atom
            in each sites. Atoms are ordered from C-terminus to N-terminus
            according to standard pdb order. For example, in CYS, the last atom
            is always the SG sulfur.
        residue_number (list[int]): The residue number for the site.
        residue_name (list[str]): The residue name (three-letter amino acid
            abbreviation) for the sites.
    """

    all_sasa_value: list[float]
    atom_sasa_values: list[list[float]]
    residue_number: list[int]
    residue_name: list[str]


def calculate_sasa(pdb_filename: str) -> SASAData:
    """Compute the SASA values for all CYS sites in a PDB file.

    Uses the ShrakeRupley algorithm implemented in `Bio.PDB.SASA.ShrakeRupley`
    with a probe radius of 1.40.

    Args:
        pdb_filename (str): The path to the PDB file.

    Returns:
        SASAData: A data class for holding SASA data from computed from a PDB
            file."""
    p = PDBParser(QUIET=True)
    struct = p.get_structure("", pdb_filename)  # type: ignore

    sr = ShrakeRupley(probe_radius=1.40, n_points=N_POINTS, radii_dict=None)

    # Calc sasa values from Residues, then from atoms
    sr.compute(struct, level="A")  # type: ignore

    # Set up dict
    res = SASAData(
        {
            "all_sasa_value": [],
            "atom_sasa_values": [],
            "residue_number": [],
            "residue_name": [],
        }
    )

    # Fill dict with CYS sites
    for x in struct.child_list:  # type: ignore
        for y in x.child_list:  # type: ignore
            for z in y.child_list:  # type: ignore
                sv = [x.sasa for x in z.child_list]  # type: ignore
                res["all_sasa_value"].append(sum(sv))  # type: ignore
                res["atom_sasa_values"].append(sv)  # type: ignore
                res["residue_number"].append(int(z.id[1]))  # type: ignore
                res["residue_name"].append(z.resname)  # type: ignore

    return res
