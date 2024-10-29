from __future__ import annotations

from biopandas.pdb import PandasPdb


def residue_pLDDT(pdb_filename: str) -> list[float]:
    """Extracts the pLDDT confidence for each residue in a PDB file.

    We assume that the pLDDT confidences are in the B-factor entries of the PDB
    file. If this information is provided at the atom-level, the maximimum value
    across the residue is used.

    Args:
        pdb_filename (str): The path to the PDB file.

    Returns:
        list[float]: The pLDDT confidence for each residue in the PDB file.
    """
    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_filename)  # type: ignore

    vals = []
    for _, res in ppdb.df["ATOM"].groupby("residue_number"):  # type: ignore
        vals.append(res["b_factor"].max())

    return vals
