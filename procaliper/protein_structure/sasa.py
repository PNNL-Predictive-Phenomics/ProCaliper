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
        entry (list[str]): An entry name corresponding the protein (typically will be UniProt ID).
        all_sasa_value (list[float]): The overall SASA value for all CYS sites.
        sg_sasa_value (list[float]): The SASA value for the CYS sites at SG atom.
        residue_id (list[int]): The residue ID for the CYS sites.
        residue_name (list[str]): The residue name for the CYS sites.
        b_factor (list[float]): The B factor for the CYS sites.
    """

    entry: list[str]
    all_sasa_value: list[float]
    sg_sasa_value: list[float]
    residue_id: list[int]
    residue_name: list[str]
    b_factor: list[float]


def calculate_sasa(pdb_filename: str, shortname: str) -> SASAData:
    """Compute the SASA values for all CYS sites in a PDB file.

    Uses the ShrakeRupley algorithm implemented in `Bio.PDB.SASA.ShrakeRupley`
    with a probe radius of 1.40.

    Args:
        pdb_filename (str): The path to the PDB file. shortname (str): The
            shortname of the protein (typically will be UniProt ID).

    Returns:
        SASAData: A data class for holding SASA data from computed from a PDB
            file."""
    p = PDBParser(QUIET=True)
    struct = p.get_structure(shortname, pdb_filename)  # type: ignore

    sr = ShrakeRupley(probe_radius=1.40, n_points=N_POINTS, radii_dict=None)

    # Calc sasa values from Residues, then from atoms
    sr.compute(struct, level="R")  # type: ignore
    sr.compute(struct, level="A")  # type: ignore

    # Set up dict
    res = SASAData(
        {
            "entry": [],
            "all_sasa_value": [],
            "sg_sasa_value": [],
            "residue_id": [],
            "residue_name": [],
            "b_factor": [],
        }
    )

    # Fill dict with CYS sites
    for x in struct.child_list:  # type: ignore
        for y in x.child_list:  # type: ignore
            for z in y.child_list:  # type: ignore
                if z.resname == "CYS":  # type: ignore
                    res["entry"].append(shortname)
                    res["all_sasa_value"].append(z.sasa)  # type: ignore
                    res["sg_sasa_value"].append(z.child_list[5].sasa)  # type: ignore
                    res["residue_id"].append(int(z.id[1]))  # type: ignore
                    res["residue_name"].append(z.resname)  # type: ignore
                    res["b_factor"].append(z.child_list[5].get_bfactor())  # type: ignore

    return res
