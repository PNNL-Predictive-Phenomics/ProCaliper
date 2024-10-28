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
        all_charge_value (list[float]): The charge value for all CYS sites (summed over atoms).
        sg_charge_value (list[float]): The charge value for the CYS sites at SG atom.
        method (list[str]): The method used for the charge calculation.
        residue_id (list[int]): The residue ID for the CYS sites.
        residue_name (list[str]): The residue name for the CYS sites.
    """

    all_charge_value: list[float]
    sg_charge_value: list[float]
    method: list[str]
    residue_id: list[int]
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
            "all_charge_value": [],
            "sg_charge_value": [],
            "method": [],
            "residue_id": [],
            "residue_name": [],
        }
    )
    residue_charges: list[float] = []
    for x in range(len(ppdb.df["ATOM"])):  # type: ignore
        if ppdb.df["ATOM"]["atom_name"][x] == "N":  # type: ignore
            residue_charges = []

        # Add charge of atom to total for residue
        residue_charges.append(charges[x])

        # Adds data to dict when CYS site read is over
        if ppdb.df["ATOM"]["atom_name"][x] == "SG":  # type: ignore
            res["all_charge_value"].append(float(sum(residue_charges)))
            res["sg_charge_value"].append(float(residue_charges[-1]))
            res["method"].append(METHOD_USED)
            res["residue_id"].append(int(ppdb.df["ATOM"]["residue_number"][x]))  # type: ignore
            res["residue_name"].append(ppdb.df["ATOM"]["residue_name"][x])  # type: ignore

    return res
