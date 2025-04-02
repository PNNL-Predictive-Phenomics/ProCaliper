from __future__ import annotations

import logging
from typing import Literal, TypedDict

import propka
import propka.run
from biopandas.pdb import PandasPdb

"""
Module for computing titration data (e.g. the acid dissociation constant $pK_a$)
for a protein at the site-level.

This module uses propka to compute the titration data by default. If `pypka` (and
its non-free dependencies) are installed, pypka can be used instead. Similarly,
if `pkAI` is installed, `pkAI` can be used instead. These are provided by the 
`procaliper[pka]` extra.
"""

logging.getLogger("propka").setLevel(logging.WARNING)
NEUTRAL_PH = 7.0


class TitrationData(TypedDict):
    """Data class for titration data.

    Array index corresponds to residue number in the PDB. Note that Python
    arrays are 0-indexed and PDB files are 1-indexed, so Python index 0
    corresponds to residue 1.

    Attributes:
        pKa (list[float | None]): The pKa values for the titration data.
            Non-titratable sites are assigned `None` values. protonation_state
        (list[tuple[str, float]]): The expected protonation states for the
            titration data. The first element of the tuple is the state of the
            site and the second element is the average protonation of the site.
            Non-titratable sites are assigned `("undefined", nan)`.
    """

    pKa: list[float | None]
    protonation_state: list[tuple[str, float | str]]


def _state_from_pk(pk: float | None) -> tuple[str, float | str]:
    """This function is modified from PypKa.titration.getProtState."""
    state = "undefined"
    if pk is not None:
        average_prot = 10 ** (pk - NEUTRAL_PH) / (1 + 10 ** (pk - NEUTRAL_PH))
    else:
        return (state, float("nan"))

    if isinstance(average_prot, str):
        return state, average_prot

    if average_prot > 0.9:
        state = "protonated"
    elif average_prot < 0.1:
        state = "deprotonated"

    return state, average_prot


def calculate_titration_propka(pdb_filename: str) -> TitrationData:
    """Uses propka to calculate titration data for the protein.

    Args:
        pdb_filename (str): The path to the PDB file.

    Returns:
        TitrationData: The titration data for the protein.
    """
    mol = propka.run.single(pdb_filename, optargs=["--quiet"], write_pka=False)
    gs = mol.conformations["AVR"].groups

    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_filename)

    seq = {
        i: res["residue_name"].iloc[0]
        for i, res in ppdb.df["ATOM"].groupby("residue_number")
    }
    pks = {group.atom.res_num: group.pka_value for group in gs}
    sv = sorted(seq.items())
    return TitrationData(
        # pKa=[group.pka_value for group in gs],
        pKa=[pks[i] if i in pks else None for i, _ in sv],
        protonation_state=[_state_from_pk(pks[i] if i in pks else 0) for i, _ in sv],
    )


try:
    from pypka import Titration

    def calculate_titration_pypka(
        pdb_filename: str,
        cpu_limit: int | None = None,
        epsin: float = 15.0,
        ionic_strength: float = 0.1,
        periodic_boundary_dims: int = 0,
        sites: str | dict[str, tuple[str,]] = "all",
    ) -> TitrationData:
        """Uses pypka to calculate titration data for the protein.
        Arguments are identical to those of `pypka.Titration`; see `pypka.Titration` for more details.

        Args:
            pdb_filename (str): The path to the PDB file.
            cpu_limit (int, optional): The number of CPUs to use. Defaults to -1 (no limit).
            epsin (float, optional): The value of epsin. Defaults to 15.0.
            ionic_strength (float, optional): The value of ionic strength. Defaults to 0.1.
            periodic_boundary_dims (int, optional): The value of periodic boundary dimensions. Defaults to 0.
            sites (str | dict[str, tuple[str,]], optional): The value of sites. Defaults to "all".

        Returns:
            TitrationData: The titration data for the protein."""

        if cpu_limit is None:
            cpu_limit = -1

        titr_params = {
            "structure": pdb_filename,
            "ncpus": cpu_limit,
            "epsin": epsin,
            "ionicstr": ionic_strength,
            "pbc_dimensions": periodic_boundary_dims,
            "sites": sites,
        }

        titr = Titration(titr_params)
        # Note: iterating through titr iterates through its sites. Each site has
        # a pK value that defaults to `None` if we do not calculate pK for that
        # site. However, we will still loop over all sites in the pdb by iterating
        # through the titr object, even if we do not calculate pK.

        return TitrationData(
            pKs=[site.pK for site in titr],  # type: ignore
            protonation_state=[site.getProtState(NEUTRAL_PH)[0] for site in titr],
        )

except ImportError:

    def calculate_titration_pypka(
        pdb_filename: str,
        cpu_limit: int | None = None,
        epsin: float = 15.0,
        ionic_strength: float = 0.1,
        periodic_boundary_dims: int = 0,
        sites: str | dict[str, tuple[str,]] = "all",
    ) -> TitrationData:
        raise ImportError(
            "pypka not installed. Install with `pip install pypka` (note that libgfortran4 will be required.)"
        )


try:
    from pkai.pKAI import pKAI

    def calculate_titration_pkai(
        pdb_filename: str,
        model_name: Literal["pKAI", "pKAI+"] = "pKAI",
        device: Literal["cpu", "gpu"] = "cpu",
        threads: int | None = None,
    ) -> TitrationData:
        """Uses pkai to calculate titration data for the protein.

        This uses a deep-learning model to predict the titration values of the
        titratable protein sites.

        Args:
            pdb_filename (str): The path to the PDB file.
            model_name (Literal["pKAI", "pKAI+"], optional): The name of the deep learning
                model to use. Defaults to "pKAI".
            device (Literal["cpu", "gpu"], optional): The device to use. Defaults to "cpu".
            threads (int, optional): The number of threads to use. Defaults to None.

        Returns:
            TitrationData: The titration data for the protein."""
        predictions = pKAI(
            pdb_filename, model_name=model_name, device=device, threads=threads
        )

        # Note: pKAI returns only titratable residues, so we will need to loop carefully.
        # To do this, we will read the pdb file and loop through the residues there, comparing
        # their names and numbers to the ones returned by pKAI.

        ppdb = PandasPdb()
        ppdb.read_pdb(pdb_filename)

        seq = sorted({i for i, _ in ppdb.df["ATOM"].groupby("residue_number")})
        pred_dict = {resnumb: pk for _, resnumb, _, pk in predictions}

        pKs = [pred_dict[i] if i in pred_dict else None for i in seq]
        states = [_state_from_pk(pk) for pk in pKs]

        return TitrationData(
            pKa=pKs,
            protonation_state=states,
        )

except ImportError:

    def calculate_titration_pkai(
        pdb_filename: str,
        model_name: Literal["pKAI", "pKAI+"] = "pKAI",
        device: Literal["cpu", "gpu"] = "cpu",
        threads: int | None = None,
    ) -> TitrationData:
        raise ImportError(
            "pkai not installed. Install with `pip install pkai` (note that only Python 3.11 or lower is supported.)"
        )
