from __future__ import annotations

from typing import TypedDict, cast

import numpy as np
from biopandas.pdb import PandasPdb

"""
Calculates the total amount of cystein sites in a protein and the total amount of atoms in a protein
Appends that data to cystein site data.

"""


class SizeData(TypedDict):
    """Data class for holding size data from computed from a PDB file.

    Attributes:
        entry (list[str]): An entry name corresponding the protein (typically will be UniProt ID).
        cys_ratio (list[float]): The ratio of CYS sites to total sites.
        min_dist_to_closest_sulfur (list[float]): The minimum distance to the closest sulfur for each CYS site.
        sulfur_closeness_rating_scaled (list[float]): The sulfur closeness rating scaled for the CYS sites.
        pLDDT (list[float]): The pLDDT values for the CYS sites.
        residue_id (list[int]): The residue ID for the CYS sites.
        residue_name (list[str]): The residue name for the CYS sites."""

    entry: list[str]
    cys_ratio: list[float]
    min_dist_to_closest_sulfur: list[float]
    sulfur_closeness_rating_scaled: list[float]
    pLDDT: list[float]
    residue_id: list[int]
    residue_name: list[str]


def calculate_size(pdb_filename: str, shortname: str) -> SizeData:
    """Calculates spatial data for a protein from a PDB file.

    Args:
        pdb_filename (str): The path to the PDB file.
        shortname (str): The shortname of the protein (typically will be UniProt ID).

    Returns:
        SizeData: A data class for holding size data from computed from a PDB file.
    """
    ppdb = PandasPdb()
    ppdb.read_pdb(pdb_filename)  # type: ignore

    res = SizeData(
        {
            "entry": [],
            "cys_ratio": [],
            "min_dist_to_closest_sulfur": [],
            "sulfur_closeness_rating_scaled": [],
            "pLDDT": [],
            "residue_id": [],
            "residue_name": [],
        }
    )

    total_residue = cast(int, max(ppdb.df["ATOM"]["residue_number"]))  # type: ignore

    cys_nums: set[int] = set()
    for x in range(len(ppdb.df["ATOM"])):  # type: ignore
        if ppdb.df["ATOM"]["residue_name"][x] == "CYS":  # type: ignore
            cys_nums.add(int(ppdb.df["ATOM"]["residue_number"][x]))  # type: ignore

    total_cys_sites = len(cys_nums)

    all_locations: list[tuple[float, float, float]] = []

    for x in range(len(ppdb.df["ATOM"])):  # type: ignore
        if ppdb.df["ATOM"]["residue_name"][x] == "CYS":  # type: ignore
            if ppdb.df["ATOM"]["atom_name"][x] == "SG":  # type: ignore
                all_locations.append(
                    (
                        ppdb.df["ATOM"]["x_coord"][x],  # type: ignore
                        ppdb.df["ATOM"]["y_coord"][x],  # type: ignore
                        ppdb.df["ATOM"]["z_coord"][x],  # type: ignore
                    )
                )

    location_index = 0

    for x in range(len(ppdb.df["ATOM"])):  # type: ignore
        if ppdb.df["ATOM"]["residue_name"][x] == "CYS":  # type: ignore
            # Adds data to dict when CYS site read is over
            if ppdb.df["ATOM"]["atom_name"][x] == "SG":  # type: ignore
                sg_closeness_rating_scaled = 0

                x_p, y_p, z_p = all_locations[location_index]
                min_distance = 1000  # Initialize with a large number

                points_excluding_index = (
                    all_locations[:location_index] + all_locations[location_index + 1 :]
                )
                for point in points_excluding_index:
                    x_q, y_q, z_q = point
                    distance = np.sqrt(
                        (x_p - x_q) ** 2 + (y_p - y_q) ** 2 + (z_p - z_q) ** 2
                    )
                    if distance < min_distance:
                        min_distance = distance
                    sg_closeness_rating_scaled += 10 / ((distance + 1) ** 2)

                location_index += 1

                res["entry"].append(shortname)
                res["cys_ratio"].append(float(total_cys_sites) / float(total_residue))
                res["min_dist_to_closest_sulfur"].append(min_distance)
                res["sulfur_closeness_rating_scaled"].append(sg_closeness_rating_scaled)
                res["pLDDT"].append(ppdb.df["ATOM"]["b_factor"][x])  # type: ignore
                res["residue_id"].append(int(ppdb.df["ATOM"]["residue_number"][x]))  # type: ignore
                res["residue_name"].append(ppdb.df["ATOM"]["residue_name"][x])  # type: ignore

    return res
