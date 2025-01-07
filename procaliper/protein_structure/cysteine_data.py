from __future__ import annotations

from typing import TypedDict, cast

import numpy as np
from biopandas.pdb import PandasPdb

"""
Computes spatial data for cysteine sites in a protein. Useful for studying
disulfide bonds.
"""


class CysteineData(TypedDict):
    """Data class for holding size data from computed from a PDB file.

    Non-CYS sites are assigned `None` values.

    Array index corresponds to residue number in the PDB. Note that Python
    arrays are 0-indexed and PDB files are 1-indexed, so Python index 0
    corresponds to residue 1.

    Attributes:
        cys_ratio (list[float | None]): The ratio of CYS sites to total sites.
        min_dist_to_closest_sulfur (list[float | None]): The minimum distance to the closest sulfur for each CYS site.
        sulfur_closeness_rating_scaled (list[float | None]): The sulfur closeness rating scaled for the CYS sites."""

    cys_ratio: list[float | None]
    min_dist_to_closest_sulfur: list[float | None]
    sulfur_closeness_rating_scaled: list[float | None]


def calculate_cysteine_data(pdb_filename: str) -> CysteineData:
    """Calculates spatial data for a protein from a PDB file.

    Args:
        pdb_filename (str): The path to the PDB file.

    Returns:
        CysteineData: A data class for holding size data from computed from a PDB file.
    """
    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_filename)

    res = CysteineData(
        {
            "cys_ratio": [],
            "min_dist_to_closest_sulfur": [],
            "sulfur_closeness_rating_scaled": [],
        }
    )

    total_residue = cast(int, max(ppdb.df["ATOM"]["residue_number"]))

    cys_positions: list[tuple[float, float, float]] = []
    for x in range(len(ppdb.df["ATOM"])):
        if ppdb.df["ATOM"]["residue_name"][x] == "CYS":
            if ppdb.df["ATOM"]["atom_name"][x] == "SG":
                cys_positions.append(
                    (
                        ppdb.df["ATOM"]["x_coord"][x],
                        ppdb.df["ATOM"]["y_coord"][x],
                        ppdb.df["ATOM"]["z_coord"][x],
                    )
                )
    total_cys_sites = len(cys_positions)

    cys_index = 0

    for _, grp in sorted(ppdb.df["ATOM"].groupby("residue_number")):
        if grp["residue_name"].max() == "CYS":
            sg_closeness_rating_scaled = 0
            x_p, y_p, z_p = cys_positions[cys_index]
            min_distance = 1000  # Initialize with a large number

            points_excluding_index = (
                cys_positions[:cys_index] + cys_positions[cys_index + 1 :]
            )
            for point in points_excluding_index:
                x_q, y_q, z_q = point
                distance = np.sqrt(
                    (x_p - x_q) ** 2 + (y_p - y_q) ** 2 + (z_p - z_q) ** 2
                )
                if distance < min_distance:
                    min_distance = distance
                sg_closeness_rating_scaled += 10 / ((distance + 1) ** 2)

            cys_index += 1

            res["cys_ratio"].append(float(total_cys_sites) / float(total_residue))
            res["min_dist_to_closest_sulfur"].append(min_distance)
            res["sulfur_closeness_rating_scaled"].append(sg_closeness_rating_scaled)
        else:
            res["cys_ratio"].append(None)
            res["min_dist_to_closest_sulfur"].append(None)
            res["sulfur_closeness_rating_scaled"].append(None)

    return res
