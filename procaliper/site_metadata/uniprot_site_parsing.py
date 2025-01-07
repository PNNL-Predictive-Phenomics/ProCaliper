from __future__ import annotations

from typing import Any

"""
Module for parsing UniProt site annotations.
"""


class SiteAnnotations:
    fields_by_description_type = {
        "BINDING": ["ligand"],
        "ACT_SITE": ["note"],
        "DNA_BIND": [],
        "DISULFID": [],
        "HELIX": [],
        "TURN": [],
        "STRAND": [],
    }

    def __init__(self, sequence: str) -> None:
        self.residue_letter: list[str] = list(sequence)
        self.residue_number: list[int] = list(range(1, len(sequence) + 1))
        self.binding: list[bool] = [False] * len(sequence)
        self.active: list[bool] = [False] * len(sequence)
        self.dna_binding: list[bool] = [False] * len(sequence)
        self.disulfide_bond: list[bool] = [False] * len(sequence)
        self.helix: list[bool] = [False] * len(sequence)
        self.turn: list[bool] = [False] * len(sequence)
        self.beta_strand: list[bool] = [False] * len(sequence)

        self.binding_data: list[dict[str, str]] = [{} for _ in range(len(sequence))]
        self.active_data: list[dict[str, str]] = [{} for _ in range(len(sequence))]

    def table(self) -> dict[str, list[Any]]:
        tbl: dict[str, list[Any]] = {}

        tbl["residue_letter"] = self.residue_letter
        tbl["residue_number"] = self.residue_number
        tbl["binding"] = self.binding
        tbl["active"] = self.active
        tbl["dna_binding"] = self.dna_binding
        tbl["disulfide_bond"] = self.disulfide_bond
        tbl["helix"] = self.helix
        tbl["turn"] = self.turn
        tbl["beta_strand"] = self.beta_strand
        tbl["binding_data"] = self.binding_data
        tbl["active_data"] = self.active_data

        return tbl

    def __len__(self) -> int:
        return len(self.residue_letter)

    def _parse_description(
        self,
        description_type: str,
        description: str,
        extract_data: bool | None = None,
    ) -> tuple[list[bool], list[dict[str, str]] | None]:
        # example of descrition:
        # DISULFID 28..87; /evidence="ECO:0000255|PROSITE-ProRule:PRU00114"; DISULFID 105; /note="Interchain (with heavy chain)"

        site_matches = [False] * len(self)

        site_data: list[dict[str, str]] | None = None

        if extract_data is None:
            extract_data = bool(self.fields_by_description_type[description_type])
        if extract_data:
            site_data = [{} for _ in range(len(self))]

        if description_type not in self.fields_by_description_type:
            raise NotImplementedError(f"Unknown description type: {description_type}")
        if (
            not description or description != description
        ):  # not-equal check is for pandas nans
            return site_matches, site_data
        if description_type not in description:
            raise ValueError(
                f"{description_type} does not appear in the description: {description}"
            )

        stretches = description.split(description_type)

        # first stretch is always empty
        for stretch in stretches[1:]:
            fields = stretch.split(";")
            # first field is always site numbers
            se = fields[0].strip().split("..")
            start, end = len(self), len(self)
            if len(se) == 1:
                start, end = int(se[0]) - 1, int(se[0]) - 1  # uniprot 1-indexes sites
            elif len(se) == 2:
                start, end = int(se[0]) - 1, int(se[1]) - 1

            if start >= len(self) or end >= len(self) or start > end:
                raise ValueError(
                    f"Improperly formatted descritpion; site numbers not recognized: {stretch} in {description}"
                )

            field_sites = list(range(start, end + 1))
            for s in field_sites:
                site_matches[s] = True

            if len(fields) == 1 or site_data is None:
                continue

            for field in fields[1:]:
                field = field.strip()
                for field_id in self.fields_by_description_type[description_type]:
                    if not field.startswith(f"/{field_id}="):
                        continue
                    field_data = field.removeprefix(f"/{field_id}=")
                    for s in field_sites:
                        if field_id not in site_data[s]:
                            site_data[s][field_id] = field_data
                        else:
                            site_data[s][field_id] += "," + field_data

        return site_matches, site_data

    def extract_annotation(
        self,
        description_type: str,
        description: str,
        extract_data: bool | None = None,
    ) -> None:
        matches, data = self._parse_description(
            description_type, description, extract_data
        )
        if description_type == "ACT_SITE":
            self.active = matches
            if data:
                self.active_data = data
        elif description_type == "BINDING":
            self.binding = matches
            if data:
                self.binding_data = data
        elif description_type == "DNA_BIND":
            self.dna_binding = matches
        elif description_type == "DISULFID":
            self.disulfide_bond = matches
        elif description_type == "STRAND":
            self.beta_strand = matches
        elif description_type == "HELIX":
            self.helix = matches
        elif description_type == "TURN":
            self.turn = matches
        else:
            raise ValueError(f"Unrecognized description type {description_type}")
