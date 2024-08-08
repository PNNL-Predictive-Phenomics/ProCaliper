from __future__ import annotations

import re
from typing import Any, cast

# import pandas as pd


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
        pass

    def _rectify_data_labels(self) -> None:
        """
        Standardize the features names in self.data
        """
        pass

    @classmethod
    def from_uniprot_row(cls, row: dict[str, Any]) -> Protein:
        p = cls()

        for key in cls.UNIPROT_SITE_PATTERNS:
            p.data[f"{key}_sites"] = p._extract_sites(
                row[key],
                cls.UNIPROT_SITE_PATTERNS[key],
            )

        p._rectify_data_labels()
        return p

    def _extract_sites(
        self, site_description: str, patterns: list[tuple[str, bool]]
    ) -> list[int]:
        sites: list[int] = []
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
