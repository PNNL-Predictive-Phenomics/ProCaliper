from __future__ import annotations

from typing import Callable

import matplotlib as plt

from .._protein import Protein

"""
Module for generating nglview widgets from proteins. Requires `nglview` be
installed. This is provided by the `procaliper[viz]` extra.
"""

try:
    import nglview
except ImportError:
    raise ImportError(
        "`nglview` is not installed. Install with `pip install nglview` or install procaliper with `pip install procaliper[viz]`"
    )


def protein_to_nglview(protein: Protein) -> nglview.NGLWidget:
    """Generates an nglview widget from a protein that has an associated PDB file.

    Must run `protein.fetch_pdb` first or specify an abosulute path to the PDB
    in `protein.pdb_location_absolute`.

    Args:
        protein (Protein): The protein object to visualize.

    Raises:
        ValueError: If the PDB location is not set.

    Returns:
        nglview.NGLWidget: an nglview widget
    """
    if not protein.pdb_location_absolute:
        raise ValueError("PDB location not set; use `fetch_pdb` first")
    return nglview.show_file(protein.pdb_location_absolute)


def _default_float_to_hex(x: float) -> str:
    return f"#{int((1-x)*255):02x}{int(1*255):02x}{int((1-x)*255):02x}"


def _default_float_to_hex_rb(x: float) -> str:
    if x < 0:
        x = -x
        return f"#{int((1-x)*255):02x}{int((1-x)*255):02x}{int(1*255):02x}"
    else:
        return f"#{int(1*255):02x}{int((1-x)*255):02x}{int((1-x)*255):02x}"


def _matplotlib_cmapper(cmap: str, two_sided: bool) -> Callable[[float], str]:
    if two_sided:
        return lambda x: plt.colors.to_hex(plt.colormaps.get_cmap(cmap)(x / 2 + 0.5))

    return lambda x: plt.colors.to_hex(plt.colormaps.get_cmap("viridis_r")(x))


def ngl_scheme(
    data: list[float],
    color_mapper: Callable[[float], str] | None = None,
    two_sided: bool = False,
    rescale: bool = True,
) -> list[tuple[str, str]]:
    """Converts a list of values to an nglview color scheme.

    Args:
        data (list[float]): The list of values to convert.
        color_mapper (Callable[[float], str] | str | None, optional): Function that
            converts a float to a hex color in the form `"#RRGGBB"`. If a string, it
            should be the name of a matplotlib colormap. If `None`, a default function
            is used that interpolates between white and green (one-sided) or red and blue (two-sided).
            Defaults to `None`.
        two_sided (bool, optional): Whether to use a two-sided color scheme. If
            `False`, we assume `data` only contains positive values. Defaults to
            `False`.
        rescale (bool, optional): Whether to rescale the values to be between
            0 and 1 (one-sided) or -1 and 1 (two-sided). Defaults to `True`.
            If `False`, the values are not rescaled, but are assumed to fall within
            [0, 1] or [-1, 1] depending on `two_sided`.

    Returns:
        list[tuple[str, str]]: A list of color and residue number tuples that
            are compatible with nglview.
    """
    if color_mapper is None:
        if two_sided:
            color_mapper = _default_float_to_hex_rb
        else:
            color_mapper = _default_float_to_hex

    if isinstance(color_mapper, str):
        color_mapper = _matplotlib_cmapper(color_mapper, two_sided)

    if rescale:
        maxx = max(data)
        scale = max(min(data), abs(maxx)) if two_sided else maxx

        if scale == 0:
            data_scaled = [0.0] * len(data)
        else:
            data_scaled = [x / maxx for x in data]
    else:
        data_scaled = data

    return [(color_mapper(x), f"{i+1}") for i, x in enumerate(data_scaled)]
