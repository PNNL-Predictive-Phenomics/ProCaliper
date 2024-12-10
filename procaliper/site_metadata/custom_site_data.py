from __future__ import annotations

from typing import Any


class CustomSiteData:
    def __init__(self, residue_number: list[int], data: dict[str, list[Any]]) -> None:
        self.residue_number = residue_number
        for key, value in data.items():
            setattr(self, key, value)

        self.keys = {"residue_number"} | set(data.keys())

    @classmethod
    def from_dict(
        cls, data: dict[str, Any], residue_index="residue_number"
    ) -> CustomSiteData:
        if residue_index not in data:
            raise ValueError("CustomSiteData must have a residue_number key.")
        return cls(data[residue_index], data)

    def table(self) -> dict[str, list[Any]]:
        return {k: getattr(self, k) for k in self.keys}

    def add_residue_numbers(self, residue_number: list[int] | int) -> None:
        if isinstance(residue_number, int):
            self.residue_number = list(range(1, residue_number + 1))
        else:
            self.residue_number = residue_number

    def add_site_data(self, key: str, row: list[Any], overwrite: bool = False) -> None:
        if hasattr(self, key) and not overwrite:
            raise KeyError(
                f"CustomSiteData already has a {key} key and overwrite is False."
            )

        if len(row) != len(self.residue_number):
            raise ValueError(
                f"CustomSiteData has {len(self.residue_number)} residues, but {key} has {len(row)} values."
                " Perhaps you forgot to call add_residue_numbers?"
            )

        setattr(self, key, row)
        self.keys.add(key)
