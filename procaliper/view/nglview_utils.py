from __future__ import annotations

from typing import Callable

from .._protein import Protein

try:
    import nglview
except ImportError:
    raise ImportError(
        "`nglview` is not installed. Install with `pip install nglview` or install procaliper with `pip install procaliper[viz]`"
    )


def protein_to_nglview(protein: Protein) -> nglview.NGLWidget:
    if not protein.pdb_location_absolute:
        raise ValueError("PDB location not set; use `fetch_pdb` first")
    return nglview.show_file(protein.pdb_location_absolute)


def _default_float_to_hex(x):
    return f"#{int((1-x)*255):02x}{int(1*255):02x}{int((1-x)*255):02x}"


def _default_float_to_hex_rb(x):
    if x < 0:
        x = -x
        return f"#{int(1*255):02x}{int((1-x)*255):02x}{int((1-x)*255):02x}"
    else:
        return f"#{int((1-x)*255):02x}{int((1-x)*255):02x}{int(1*255):02x}"


def ngl_scheme(
    data: list[float],
    float_to_hex: Callable[[float], str] | None = None,
    two_sided: bool = False,
) -> list[tuple[str, str]]:
    if float_to_hex is None:
        if two_sided:
            float_to_hex = _default_float_to_hex_rb
        else:
            float_to_hex = _default_float_to_hex

    maxx = max(data)
    scale = max(min(data), abs(maxx)) if two_sided else maxx

    if scale == 0:
        data_scaled = [0] * len(data)
    else:
        data_scaled = [x / maxx for x in data]

    return [(float_to_hex(x), f"{i+1}") for i, x in enumerate(data_scaled)]
