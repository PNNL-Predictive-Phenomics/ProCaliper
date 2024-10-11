from __future__ import annotations

from typing import TypedDict, cast

from biopandas.pdb import PandasPdb
from openbabel import openbabel as ob
from openbabel import pybel

from .hyperparameters import HYPERPARAMETERS

# Removes annoying warning messages
pybel.ob.obErrorLog.SetOutputLevel(0)  # type: ignore

# METHOD_USED determines the method used for charge_calculations. examples ('qtpie', 'eem', 'gasteiger')
# For a full list reference https://open-babel.readthedocs.io/en/latest/Charges/charges.html
METHOD_USED = str(HYPERPARAMETERS.get("charge_method_used"))

"""
    Takes a pdb file and the method used ('qtpie', 'eem', etc) 
        and returns a dict of charge values for CYS sites.
"""


class ChargeData(TypedDict):
    """
    A data class for holding charge data from computed from a PDB file.

    Attributes:
        charges (list[list[float]]): The charge value for atoms in the residue,
            ordered from C-terminus to N-terminus according to standard pdb order.
            For example, in CYS, the last atom is always the SG sulfur.
        method (list[str]): The method used for the charge calculation.
        residue_number (list[int]): The residue number for the site.
        residue_name (list[str]): The residue name (three-letter amino acid
            abbreviation) for the sites.
    """

    charge: list[list[float]]
    charge_method: list[str]
    residue_number: list[int]
    residue_name: list[str]


def calculate_charge(pdb_filename: str) -> ChargeData:
    """Computes the charge of CYS sites in a PDB file.

    By default, the method used is 'gasteiger', but this is configurable in
    `hyperparameters.py`.

    Args:
        pdb_filename (str): The path to the PDB file. shortname (str): The
            shortname of the protein (typically will be UniProt ID).

    Raises:
        ValueError: If the charge method is not found.

    Returns:
        ChargeData: A data class for holding charge data from computed from a
            PDB file.
    """
    pbmol = next(pybel.readfile("pdb", pdb_filename))  # type: ignore
    mol = pbmol.OBMol  # type: ignore

    # Applies the model and computes charges.
    ob_charge_model = ob.OBChargeModel.FindType(METHOD_USED)  # type: ignore

    if not ob_charge_model:  # type: ignore
        raise ValueError("Charge method not found. Please check hyperparameters.py")
    ob_charge_model.ComputeCharges(mol)  # type: ignore

    charges = cast(list[float], ob_charge_model.GetPartialCharges())  # type: ignore

    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_filename)  # type: ignore

    # Set up dict
    res = ChargeData(
        {
            "charge": [],
            "charge_method": [],
            "residue_number": [],
            "residue_name": [],
        }
    )

    for res_num, residue in ppdb.df["ATOM"].groupby("residue_number"):
        res["charge"].append([charges[x - 1] for x in sorted(residue["atom_number"])])
        res["charge_method"].append(METHOD_USED)
        res["residue_number"].append(int(res_num))
        res["residue_name"].append(
            residue["residue_name"].iloc[0]
        )  # all residue names should be the same because these atoms are from the same residue

    return res
