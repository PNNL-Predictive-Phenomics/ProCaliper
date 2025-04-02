from __future__ import annotations

from typing import Any

"""
Module for storing custom residue-level data that can be loaded in from
user-provided tables.
"""


class CustomSiteData:
    """Class for storing custom site-level data."""

    def __init__(self, residue_number: list[int], data: dict[str, list[Any]]) -> None:
        self.residue_number = residue_number
        for key, value in data.items():
            setattr(self, key, value)

        self.keys = {"residue_number"} | set(data.keys())

    @classmethod
    def from_dict(
        cls,
        data: dict[str, list[Any]],
        residue_index_feature_name: str = "residue_number",
    ) -> CustomSiteData:
        """Create a CustomSiteData object from a dictionary of data.

        Args:
            data (dict[str, list[Any]]): Data dictionary indexed by feature
                name. Each value must be a list of the same length as the
                residue number feature. Must include a residue number key.
            residue_index_feature_name (str, optional): The name of the feature
                that contains the residue number. Defaults to "residue_number".

        Raises:
            ValueError: If the residue number feature is not in the data.

        Returns:
            CustomSiteData: A CustomSiteData object. that contains the data.
        """
        if residue_index_feature_name not in data:
            raise ValueError("CustomSiteData must have a residue_number key.")
        return cls(data[residue_index_feature_name], data)

    def table(self) -> dict[str, list[Any]]:
        """Return a dictionary of the data in the CustomSiteData object.

        Returns:
            dict[str, list[Any]]: A dictionary of the data in the CustomSiteData
                object.
        """
        return {k: getattr(self, k) for k in self.keys}

    def add_residue_numbers(self, residue_number: list[int] | int) -> None:
        """Specify the number of residues in the CustomSiteData object.

        Args:
            residue_number (list[int] | int): If an integer, the number of
                residues. If a list of integers, the list of residue numbers.
        """
        if isinstance(residue_number, int):
            self.residue_number = list(range(1, residue_number + 1))
        else:
            self.residue_number = residue_number

    def add_site_data(self, key: str, row: list[Any], overwrite: bool = False) -> None:
        """Add a site-level feature to the CustomSiteData object.

        Args:
            key (str): The name of the feature to add.
            row (list[Any]): The values for the feature.

            overwrite (bool, optional): Whether to overwrite an existing
                feature. Defaults to False.

        Raises:
            KeyError: If overwrite is False and the feature already exists.
            ValueError: If the number of values in the feature does not match
                the number of residues.
        """
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
