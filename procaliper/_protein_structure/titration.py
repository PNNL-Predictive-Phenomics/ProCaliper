from __future__ import annotations

from typing import TypedDict

from pypka import Titration

NEUTRAL_PH = 7.0


class TitrationData(TypedDict):
    residue_names: list[str]
    residue_numbers: list[int]
    pKs: list[float]
    states: list[tuple[str, float]]


def calculate_titration(
    pdb_filename: str,
    cpu_limit: int | None = None,
    epsin: float = 15.0,
    ionic_strength: float = 0.1,
    periodic_boundary_dims: int = 0,
    sites: str | dict[str, tuple[str,]] = "all",
) -> TitrationData:
    titr_params = {
        "structure": pdb_filename,
        "ncpus": cpu_limit,
        "epsin": epsin,
        "ionicstr": ionic_strength,
        "pbc_dimensions": periodic_boundary_dims,
        "sites": sites,
    }

    titr = Titration(titr_params)

    return TitrationData(
        residue_names=[site.res_name for site in titr],  # type: ignore
        residue_numbers=[site.res_number for site in titr],  # type: ignore
        pKs=[site.pK for site in titr],  # type: ignore
        states=[site.getProtState(NEUTRAL_PH)[0] for site in titr],  # type: ignore
    )
