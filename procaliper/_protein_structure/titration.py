from __future__ import annotations

from typing import Literal, TypedDict

from pkai.pKAI import pKAI

NEUTRAL_PH = 7.0


class TitrationData(TypedDict):
    residue_names: list[str]
    residue_numbers: list[int]
    pKs: list[float]
    states: list[tuple[str, float | str]]


def _state_from_pk(pk: float | None) -> tuple[str, float | str]:
    """This function is modified from PypKa.titration.getProtState.

    Parameters
    ----------
    pk : float | None
        _description_

    Returns
    -------
    tuple[str, float | str]
        _description_
    """
    state = "undefined"
    if pk is not None:
        average_prot = 10 ** (pk - NEUTRAL_PH) / (1 + 10 ** (pk - NEUTRAL_PH))
    else:
        return (state, "pk Not In Range")

    if isinstance(average_prot, str):
        return state, average_prot

    if average_prot > 0.9:
        state = "protonated"
    elif average_prot < 0.1:
        state = "deprotonated"

    return state, average_prot


def estimate_titration(
    pdb_filename: str,
    model_name: Literal["pKAI", "pKAI+"] = "pKAI",
    device: Literal["cpu", "gpu"] = "cpu",
    threads=None,
) -> TitrationData:
    predictions = pKAI(
        pdb_filename, model_name=model_name, device=device, threads=threads
    )

    residue_names = []
    residue_numbers = []
    pKs = []
    states = []
    for _, resnumb, resname, pk in predictions:  # we do not use the chain here
        residue_names.append(resname)
        residue_numbers.append(resnumb)
        pKs.append(pk)
        states.append(_state_from_pk(pk))

    return TitrationData(
        residue_names=residue_names,
        residue_numbers=residue_numbers,
        pKs=pKs,
        states=states,
    )


try:
    from pypka import Titration

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

except ImportError:

    def calculate_titration(
        pdb_filename: str,
        cpu_limit: int | None = None,
        epsin: float = 15.0,
        ionic_strength: float = 0.1,
        periodic_boundary_dims: int = 0,
        sites: str | dict[str, tuple[str,]] = "all",
    ) -> TitrationData:
        raise ImportError("pypka not installed. Install with `pip install pypka`")
