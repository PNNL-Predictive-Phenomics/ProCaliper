from __future__ import annotations

from typing import TypedDict, cast

from biopandas.pdb import PandasPdb
from openbabel import openbabel as ob
from openbabel import pybel

# Removes annoying warning messages
pybel.ob.obErrorLog.SetOutputLevel(0)

"""
    Takes a pdb file and the method used ('qtpie', 'eem', etc) 
        and returns a dict of charge values for residue sites.
"""


class ChargeData(TypedDict):
    """
    A data class for holding charge data from computed from a PDB file.

    Array index corresponds to residue number in the PDB. Note that Python
    arrays are 0-indexed and PDB files are 1-indexed, so Python index 0
    corresponds to residue 1.

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


def calculate_charge(pdb_filename: str, method: str = "gasteiger") -> ChargeData:
    """Computes the charge of residue sites in a PDB file.

    By default, the method used is 'gasteiger', but this is configurable in
    `hyperparameters.py`.

    Args:
        pdb_filename (str): The path to the PDB file. shortname (str): The
            shortname of the protein (typically will be UniProt ID).
        method (str, optional): The method used for the charge calculation.
            Examples include 'qtpie', 'eem', 'gasteiger'. Defaults to
            'gasteiger'. For a full list reference
            https://open-babel.readthedocs.io/en/latest/Charges/charges.html


    Raises:
        ValueError: If the charge method is not found.

    Returns:
        ChargeData: A data class for holding charge data from computed from a
            PDB file.
    """
    pbmol = next(pybel.readfile("pdb", pdb_filename))
    mol = pbmol.OBMol

    # Applies the model and computes charges.
    ob_charge_model = ob.OBChargeModel.FindType(method)

    if not ob_charge_model:
        raise ValueError("Charge method not found. Please check hyperparameters.py")
    ob_charge_model.ComputeCharges(mol)

    charges = cast(list[float], ob_charge_model.GetPartialCharges())

    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_filename)

    # Set up dict
    res = ChargeData(
        {
            "charge": [],
            "charge_method": [],
        }
    )

    for _, residue in sorted(ppdb.df["ATOM"].groupby("residue_number")):
        res["charge"].append([charges[x - 1] for x in sorted(residue["atom_number"])])
        res["charge_method"].append(method)

    return res
