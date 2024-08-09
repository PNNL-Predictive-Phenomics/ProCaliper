from __future__ import annotations

from typing import TypedDict

from Bio.PDB import Entity, PDBParser
from Bio.PDB.SASA import ShrakeRupley

from .hyperparameters import HYPERPARAMETERS

N_POINTS = int(
    HYPERPARAMETERS["n_points"]
)  # Dictates the accuracy of the ShrakeRupley calculation. Higher values are more accurate but slower to calculate. Default = 100

"""
    returns a dict whith names "Entry", ,"all_sasa_value", "sg_sasa_value", "residue_id", 'residue_name'

"""


class SASAData(TypedDict):
    entry: list[str]
    all_sasa_value: list[float]
    sg_sasa_value: list[float]
    residue_id: list[int]
    residue_name: list[str]
    b_factor: list[float]


def calculate_sasa(pdb_filename: str, shortname: str):
    p = PDBParser(QUIET=True)
    struct = p.get_structure(shortname, pdb_filename)  # type: ignore
    if not isinstance(struct, Entity):
        raise TypeError("`struct` must be a Bio.PDB.Entity")
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
