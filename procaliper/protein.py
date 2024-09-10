"""
Reference:
https://github.com/PhenoMeters/PTM_ML/tree/hunter/DatabaseScripts
"""

from __future__ import annotations

import os
import re
from typing import Any, cast

import pandas as pd
import requests
from UniProtMapper import ProtMapper

import procaliper._protein_structure as structure
from procaliper.type_aliases import AminoAcidLetter


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

    UNIPROT_SITE_PATTERNS_RECTIFIED = {
        "active_site": [(r"ACT_SITE (\d+);", False)],
        "binding_site": [
            (r"BINDING (\d+);", False),
            (r"BINDING (\d+)\.\.(\d+);", False),
        ],
        "dna_binding": [
            (r"DNA_BIND (\d+);", False),
            (r"DNA_BIND (\d+)\.\.(\d+);", False),
        ],
        "disulfide_bond": [
            (r"DISULFID (\d+);", False),
            (r"DISULFID (\d+)\.\.(\d+);", False),
        ],
        "beta_strand": [(r"STRAND (\d+);", True), (r"STRAND (\d+)\.\.(\d+);", True)],
        "helix": [(r"HELIX (\d+);", True), (r"HELIX (\d+)\.\.(\d+);", True)],
        "turn": [(r"TURN (\d+);", True), (r"TURN (\d+)\.\.(\d+);", True)],
    }

    UNIPROT_API_DEFAULT_FIELDS = [
        "id",
        "reviewed",
        "protein_name",
        "gene_names",
        "organism_name",
        "length",
        "sequence",
        "ft_act_site",
        "ft_binding",
        "ft_dna_bind",
        "ft_disulfid",
        "ft_strand",
        "ft_helix",
        "ft_turn",
    ]

    def __init__(self) -> None:
        self.data: dict[str, Any] = {}
        self.pdb_location_relative: str | None = None
        self.pdb_location_absolute: str | None = None

        self.sasa_data: structure.sasa.SASAData | None = None
        self.charge_data: structure.charge.ChargeData | None = None
        self.size_data: structure.size.SizeData | None = None
        self.titration_data: structure.titration.TitrationData | None = None
        pass

    def _rectify_data_labels(self) -> None:
        """
        Standardize the features names in self.data

        Replaces all spaces with underscores and lowercases the keys
        """
        for k in list(self.data.keys()):
            new_key = k.replace(" ", "_").lower()
            self.data[new_key] = self.data.pop(k)

    @classmethod
    def from_uniprot_row(cls, row: dict[str, Any]) -> Protein:
        p = cls()
        p.data["sequence"] = row["Sequence"]

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

    @classmethod
    def from_uniprot_id(
        cls, uniprot_id: str, fields: list[str] | None = None
    ) -> Protein:
        if not fields:
            fields = cls.UNIPROT_API_DEFAULT_FIELDS

        mapper = ProtMapper()

        result, error = mapper.get(ids=[uniprot_id], fields=fields)  # type: ignore
        if error:
            raise ValueError(f"Uniprot id {error} not retrieved")
        result.rename(columns={"From": "entry"}, inplace=True)
        if "Length" in result.columns:
            result["Length"] = pd.to_numeric(result["Length"])  # type: ignore
        return cls.from_uniprot_row(result.iloc[0].to_dict())  # type: ignore

    @classmethod
    def list_from_uniprot_ids(
        cls, uniprot_ids: list[str], fields: list[str] | None = None
    ) -> list[Protein]:
        if not fields:
            fields = cls.UNIPROT_API_DEFAULT_FIELDS

        mapper = ProtMapper()

        result, error = mapper.get(ids=uniprot_ids, fields=fields)  # type: ignore
        if error:
            raise ValueError(f"Uniprot id {error} not retrieved")
        result.rename(columns={"From": "entry"}, inplace=True)

        if "Length" in result.columns:
            result["Length"] = pd.to_numeric(result["Length"])  # type: ignore
        return [cls.from_uniprot_row(row.to_dict()) for _, row in result.iterrows()]  # type: ignore

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Protein):
            return False
        return (
            self.data == other.data
            and self.sasa_data == other.sasa_data
            and self.charge_data == other.charge_data
            and self.size_data == other.size_data
        )

    def get_sasa(self) -> structure.sasa.SASAData:
        if self.sasa_data:
            return self.sasa_data

        if self.pdb_location_absolute:
            self.sasa_data = structure.sasa.calculate_sasa(
                self.pdb_location_absolute,
                self.data["entry"],
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
                self.data["entry"],
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
                self.data["entry"],
            )
            return self.size_data
        else:
            raise ValueError(
                "Size data not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def get_titration(self) -> structure.titration.TitrationData:
        if self.titration_data:
            return self.titration_data

        if self.pdb_location_absolute:
            self.titration_data = structure.titration.calculate_titration(
                self.pdb_location_absolute,
            )
            return self.titration_data
        else:
            raise ValueError(
                "Titration data not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def get_titration_estimate(self) -> structure.titration.TitrationData:
        if self.titration_data:
            return self.titration_data

        if self.pdb_location_absolute:
            self.titration_data = structure.titration.estimate_titration(
                self.pdb_location_absolute,
            )
            return self.titration_data
        else:
            raise ValueError(
                "Titration data not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def unravel_sites(
        self,
        selected_aas: None | set[AminoAcidLetter] = None,
        selected_keys: None | set[str] = None,
    ) -> dict[str, list[Any]]:
        if not selected_keys:
            selected_keys = set(self.data.keys()) - {"sequence"}

        site_keys = set(Protein.UNIPROT_SITE_PATTERNS_RECTIFIED.keys()) & selected_keys
        other_keys = selected_keys - site_keys

        res: dict[str, list[Any]] = {k: [] for k in other_keys | site_keys}
        res["letter"] = []
        res["position"] = []

        for index, site in enumerate(self.data["sequence"]):
            site_dict: dict[str, Any] = {k: self.data[k] for k in other_keys}
            site_dict["letter"] = site
            site_dict["position"] = index + 1
            if selected_aas and site not in selected_aas:
                continue

            for key in site_keys:
                site_dict[key] = index in self.data[f"{key}_sites"]

            for k, v in site_dict.items():
                res[k].append(v)

        return res

    def fetch_pdb(self, save_path: str | None = None, url: str | None = None) -> None:
        if not url:
            url = f"https://alphafold.ebi.ac.uk/files/AF-{self.data['entry']}-F1-model_v4.pdb"
        if not save_path:
            save_path = f"{self.data['entry']}.pdb"

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
        if "sequence" not in self.data:
            raise ValueError("Sequence entry not found in data")

        sequence = self.data["sequence"]

        return site <= len(sequence) and sequence[site - 1] == aa
