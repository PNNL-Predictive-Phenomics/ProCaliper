from __future__ import annotations

from typing import TypedDict

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.SASA import ShrakeRupley

N_POINTS = 100
PROBE_RADIUS = 1.40


class SASAData(TypedDict):
    """Data class for holding SASA data from computed from a PDB file.

    Array index corresponds to residue number in the PDB. Note that Python
    arrays are 0-indexed and PDB files are 1-indexed, so Python index 0
    corresponds to residue 1.

    Attributes:
        all_sasa_value (list[float]): The overall SASA value for each site
            (computed as sum of atom SASA values).
        atom_sasa_values (list[list[float]]): The SASA value for the each atom
            in each sites. Atoms are ordered from C-terminus to N-terminus
            according to standard pdb order. For example, in CYS, the last atom
            is always the SG sulfur.
    """

    all_sasa_value: list[float]
    atom_sasa_values: list[list[float]]


def calculate_sasa(pdb_filename: str) -> SASAData:
    """Compute the SASA values for all CYS sites in a PDB file.

    Uses the ShrakeRupley algorithm implemented in `Bio.PDB.SASA.ShrakeRupley`
    with a probe radius of 1.40 and 100 points.

    Args:
        pdb_filename (str): The path to the PDB file.

    Returns:
        SASAData: A data class for holding SASA data from computed from a PDB
            file."""
    p = PDBParser(QUIET=True)
    struct = p.get_structure("", pdb_filename)  # type: ignore

    sr = ShrakeRupley(probe_radius=PROBE_RADIUS, n_points=N_POINTS, radii_dict=None)

    # Calc sasa values from Residues (from atoms)
    sr.compute(struct, level="R")  # type: ignore

    # Set up dict
    res = SASAData(
        {
            "all_sasa_value": [],
            "atom_sasa_values": [],
        }
    )

    # Fill dict with CYS sites
    for x in struct.child_list:  # type: ignore
        for y in x.child_list:  # type: ignore
            for z in y.child_list:  # type: ignore
                res["all_sasa_value"].append(z.sasa)  # type: ignore
                res["atom_sasa_values"].append([x.sasa for x in z.child_list])  # type: ignore

    return res
