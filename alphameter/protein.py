"""
Reference:
https://github.com/PhenoMeters/PTM_ML/tree/hunter/DatabaseScripts
"""

from __future__ import annotations

import os
import re
from typing import Any, cast

import requests

import alphameter._protein_structure as structure
from alphameter.type_aliases import AminoAcidLetter


class Protein:
    UNIPROT_SITE_PATTERNS = {
        "Active site": [(r"ACT_SITE (\d+);", False)],
        "Binding site": [
            (r"BINDING (\d+);", False),
            (r"BINDING (\d+)\.\.(\d+);", False),
        ],
        "DNA binding": [
            (r"DNA_BIND (\d+);", False),
            (r"DNA_BIND (\d+)\.\.(\d+);", False),
        ],
        "Disulfide bond": [
            (r"DISULFID (\d+);", False),
            (r"DISULFID (\d+)\.\.(\d+);", False),
        ],
        "Beta strand": [(r"STRAND (\d+);", True), (r"STRAND (\d+)\.\.(\d+);", True)],
        "Helix": [(r"HELIX (\d+);", True), (r"HELIX (\d+)\.\.(\d+);", True)],
        "Turn": [(r"TURN (\d+);", True), (r"TURN (\d+)\.\.(\d+);", True)],
    }

    def __init__(self) -> None:
        self.data: dict[str, Any] = {}
        self.pdb_location_relative: str | None = None
        self.pdb_location_absolute: str | None = None

        self.sasa_data: structure.sasa.SASAData | None = None
        self.charge_data: structure.charge.ChargeData | None = None
        self.size_data: structure.size.SizeData | None = None
        pass

    def _rectify_data_labels(self) -> None:
        """
        Standardize the features names in self.data
        """
        pass

    @classmethod
    def from_uniprot_row(cls, row: dict[str, Any]) -> Protein:
        p = cls()
        p.data["Sequence"] = row["Sequence"]

        for key, value in row.items():
            if key in cls.UNIPROT_SITE_PATTERNS:
                p.data[f"{key}_sites"] = p._extract_sites(
                    value,
                    cls.UNIPROT_SITE_PATTERNS[key],
                )
                p.data[f"{key}_cysteine_sites"] = [
                    site
                    for site in p.data[f"{key}_sites"]
                    if p._is_site_aa(site, aa="C")
                ]
            else:
                p.data[key] = value

        p._rectify_data_labels()
        return p

    def get_sasa(self) -> structure.sasa.SASAData:
        if self.sasa_data:
            return self.sasa_data

        if self.pdb_location_absolute:
            self.sasa_data = structure.sasa.calculate_sasa(
                self.pdb_location_absolute,
                self.data["Entry"],
            )
            return self.sasa_data
        else:
            raise ValueError(
                "SASA data not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def get_charge(self) -> structure.charge.ChargeData:
        if self.charge_data:
            return self.charge_data

        if self.pdb_location_absolute:
            self.charge_data = structure.charge.calculate_charge(
                self.pdb_location_absolute,
                self.data["Entry"],
            )
            return self.charge_data
        else:
            raise ValueError(
                "Charge data not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def get_size(self) -> structure.size.SizeData:
        if self.size_data:
            return self.size_data

        if self.pdb_location_absolute:
            self.size_data = structure.size.calculate_size(
                self.pdb_location_absolute,
                self.data["Entry"],
            )
            return self.size_data
        else:
            raise ValueError(
                "Size data not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def unravel_sites(
        self,
        selected_aas: None | set[AminoAcidLetter] = None,
        selected_keys: None | set[str] = None,
    ) -> list[dict[str, Any]]:
        res: list[dict[str, Any]] = []

        if not selected_keys:
            selected_keys = set(self.data.keys()) - {"Sequence"}

        site_keys = set(Protein.UNIPROT_SITE_PATTERNS.keys()) & selected_keys
        other_keys = selected_keys - site_keys

        for index, site in enumerate(self.data["Sequence"]):
            site_dict: dict[str, Any] = {k: self.data[k] for k in other_keys}
            site_dict["Letter"] = site
            site_dict["Position"] = index + 1
            if selected_aas and site not in selected_aas:
                continue

            for key in site_keys:
                site_dict[key] = index in self.data[f"{key}_sites"]

            res.append(site_dict)

        return res

    def fetch_pdb(self, save_path: str | None = None, url: str | None = None) -> None:
        if not url:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{self.data['Entry']}-F1-model_v4.pdb"
        if not save_path:
            save_path = f"{self.data['Entry']}.pdb"

        response = requests.get(url)

        if response.status_code != 200:
            raise Exception(f"Failed to fetch PDB: {response.status_code}")

        with open(save_path, "wb+") as f:
            f.write(response.content)

        self.pdb_location_relative = save_path
        self.pdb_location_absolute = os.path.abspath(save_path)

    def _extract_sites(
        self, site_description: str, patterns: list[tuple[str, bool]]
    ) -> list[int]:
        sites: list[int] = []
        if (
            str(site_description) == "nan"
        ):  # this will be the missing value default in pandas--is there a more elegant way to handle this?
            return sites
        for pattern, expand_range in patterns:
            print(site_description)
            matches = cast(list[str], re.findall(pattern, site_description))

            for match in matches:
                if isinstance(match, tuple):
                    start, end = int(match[0]), int(match[1])
                    if expand_range:
                        sites.extend(
                            range(start, end + 1)
                        )  # Include all values in the range from start to end
                    else:
                        sites.extend([start, end])  # Add start and end points
                else:
                    sites.append(int(match))
        return sites

    def _is_site_aa(self, site: int, aa: AminoAcidLetter = "C") -> bool:
        """_summary_

        Parameters
        ----------
        site : int
            Position of amino acid (1-indexed)
        aa : str, optional
            Amino acid code to test against by default "C" (cysteine)

        Returns
        -------
        bool
            True if the amino acid is a the position `site` is as specified.

        Raises
        ------
        ValueError
            If the protein does not have a defined sequence.
        """
        if "Sequence" not in self.data:
            raise ValueError("Sequence entry not found in data")

        sequence = self.data["Sequence"]

        return site <= len(sequence) and sequence[site - 1] == aa
