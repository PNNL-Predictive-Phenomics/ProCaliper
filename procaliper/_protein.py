"""
Contains the main class for holding single protein data.
"""

from __future__ import annotations

import os
import re
from itertools import chain
from typing import Any, cast

import pandas as pd
import requests
from UniProtMapper import ProtMapper

import procaliper.protein_structure as structure
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

        self.confidence_data: list[float] | None = None
        self.sasa_data: structure.sasa.SASAData | None = None
        self.charge_data: structure.charge.ChargeData | None = None
        self.cysteine_data: structure.cysteine_data.CysteineData | None = None
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
        """Create a new Protein object from a row from a Uniprot table

        Args:
            row (dict[str, Any]): Contains the data from the Uniprot table. Must
                have "Sequence" or "sequence" as a key.

        Raises:
            ValueError: If "Sequence" or "sequence" is not found in the row.

        Returns:
            Protein: A processed and standardized protein object.
        """
        p = cls()
        if "Sequence" in row:
            p.data["sequence"] = row["Sequence"]
        elif "sequence" in row:
            p.data["sequence"] = row["sequence"]
        else:
            raise ValueError(f"Sequence not found in row: {row}")

        for key, value in row.items():
            if key in cls.UNIPROT_SITE_PATTERNS:
                p.data[f"{key}_sites"] = p._extract_sites(
                    value,
                    cls.UNIPROT_SITE_PATTERNS[key],
                )
                # p.data[f"{key}_cysteine_sites"] = [
                #     site
                #     for site in p.data[f"{key}_sites"]
                #     if p._is_site_aa(site, aa="C")
                # ]
            else:
                p.data[key] = value

        p._rectify_data_labels()
        return p

    @classmethod
    def from_uniprot_id(
        cls,
        uniprot_id: str,
        fields: list[str] | None = None,
        from_db: str = "UniProtKB_AC-ID",
        to_db: str = "UniProtKB-Swiss-Prot",
    ) -> Protein:
        """Create a new Protein object from a Uniprot ID (fetches with Uniprot API)

        Args:
            uniprot_id (str): The Uniprot ID of the protein.
            fields (list[str] | None, optional): The fields to retrieve from
                Uniprot. If `None`, `Protein.UNIPROT_API_DEFAULT_FIELDS` is used.
            from_db (str, optional): The database to retrieve the ID from.
                Defaults to "UniProtKB_AC-ID".
            to_db (str, optional): The database to map to.
                Defaults to "UniProtKB-Swiss-Prot".

        Raises:
            ValueError: If we cannot retrieve the Uniprot ID.

        Returns:
            Protein: A processed and standardized protein object.
        """

        if not fields:
            fields = cls.UNIPROT_API_DEFAULT_FIELDS

        mapper = ProtMapper()

        result, error = mapper.get(  # type: ignore
            ids=[uniprot_id], fields=fields, from_db=from_db, to_db=to_db
        )
        if error:
            raise ValueError(f"Uniprot id not retrieved: {error}")
        result.rename(columns={"From": "entry"}, inplace=True)
        if "Length" in result.columns:
            result["Length"] = pd.to_numeric(result["Length"])  # type: ignore
        return cls.from_uniprot_row(result.iloc[0].to_dict())  # type: ignore

    @classmethod
    def list_from_uniprot_ids(
        cls,
        uniprot_ids: list[str],
        fields: list[str] | None = None,
        from_db: str = "UniProtKB_AC-ID",
        to_db: str = "UniProtKB-Swiss-Prot",
    ) -> list[Protein]:
        """Create a list of Protein objects from a list of Uniprot IDs (fetches with Uniprot API)

        Args:
            uniprot_ids (list[str]): The Uniprot IDs of the proteins.
            fields (list[str] | None, optional): The fields to retrieve from
                Uniprot. If `None`, `Protein.UNIPROT_API_DEFAULT_FIELDS` is used.
            from_db (str, optional): The database to retrieve the IDs from.
                Defaults to "UniProtKB_AC-ID".
            to_db (str, optional): The database to map to.
                Defaults to "UniProtKB-Swiss-Prot".

        Raises:
            ValueError: If we cannot retrieve the Uniprot IDs.

        Returns:
            list[Protein]: A list of processed and standardized protein objects.
        """
        if not fields:
            fields = cls.UNIPROT_API_DEFAULT_FIELDS

        mapper = ProtMapper()

        result, error = mapper.get(  # type: ignore
            ids=uniprot_ids, fields=fields, from_db=from_db, to_db=to_db
        )
        if error:
            raise ValueError(f"Uniprot id not retrieved: {error}")
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
            and self.cysteine_data == other.cysteine_data
        )

    def residue_data_frame(self) -> pd.DataFrame:
        d = dict(
            chain(
                self.get_charge().items(),
                self.get_sasa().items(),
                self.get_cysteine_data().items(),
                self.get_titration().items(),
            )
        )
        d["pLDDT"] = self.get_confidence()

        return pd.DataFrame(d)

    def get_confidence(self) -> list[float]:
        """Fetches precomputed confidence data from pdb file.

        Must run `self.fetch_pdb` first or specify an abosulute path to the PDB
        file in `self.pdb_location_absolute`.

        Raises:
            ValueError: If `confidence_data` is not already stored and
                `pdb_location_absolute` is not set.

        Returns:
            list[float]: A list of confidence values for each residue.
        """
        if self.confidence_data:
            return self.confidence_data

        if self.pdb_location_absolute:
            self.confidence_data = structure.confidence.residue_pLDDT(
                self.pdb_location_absolute,
            )
            return self.confidence_data
        else:
            raise ValueError(
                "Confidence data not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def get_sasa(self) -> structure.sasa.SASAData:
        """Fetches precomputed SASA data for the protein, or computes it.

        Must run `self.fetch_pdb` first or specify an abosulute path to the PDB
        file in `self.pdb_location_absolute`.

        Raises:
            ValueError: If `sasa_data` is not already stored and
                `pdb_location_absolute` is not set.

        Returns:
            structure.sasa.SASAData: A :class:`protein_structure.sasa.SASAData`
                object containing the SASA values for residues and atoms.
        """
        if self.sasa_data:
            return self.sasa_data

        if self.pdb_location_absolute:
            self.sasa_data = structure.sasa.calculate_sasa(
                self.pdb_location_absolute,
            )
            return self.sasa_data
        else:
            raise ValueError(
                "SASA data not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def get_charge(self, method="gasteiger") -> structure.charge.ChargeData:
        """Fetches precomputed charge data for the protein, or computes it.

        Must run `self.fetch_pdb` first or specify an abosulute path to the PDB
        file in `self.pdb_location_absolute`.

        Args:
            method (str, optional): The method used for the charge calculation.
                Examples include 'qtpie', 'eem', 'gasteiger'. Defaults to
                'gasteiger'. For a full list reference
                https://open-babel.readthedocs.io/en/latest/Charges/charges.html

        Raises:
            ValueError: If `charge_data` is not already stored and
                `pdb_location_absolute` is not set.

        Returns:
            structure.charge.ChargeData: A :class:`protein_structure.charge.ChargeData`
                object containing the charge values for residues and atoms.
        """
        if self.charge_data:
            if self.charge_data["charge_method"]:
                if self.charge_data["charge_method"][0] == method:
                    return self.charge_data

        if self.pdb_location_absolute:
            self.charge_data = structure.charge.calculate_charge(
                self.pdb_location_absolute,
                method=method,
            )

            self.last_charge_method = method

            return self.charge_data
        else:
            raise ValueError(
                "Charge data for specified method not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def get_cysteine_data(self) -> structure.cysteine_data.CysteineData:
        """Fetches precomputed size data for the protein, or computes it.

        Must run `self.fetch_pdb` first or specify an abosulute path to the PDB
        file in `self.pdb_location_absolute`.

        Raises:
            ValueError: If `cysteine_data` is not already stored and
                `pdb_location_absolute` is not set.

        Returns:
            structure.size.CysteineData: A :class:`protein_structure.size.CysteineData`
                object containing the size values for cystein sites.
        """
        if self.cysteine_data:
            return self.cysteine_data

        if self.pdb_location_absolute:
            self.cysteine_data = structure.cysteine_data.calculate_cysteine_data(
                self.pdb_location_absolute,
            )
            return self.cysteine_data
        else:
            raise ValueError(
                "Size data not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def get_titration(self) -> structure.titration.TitrationData:
        """Runs the default titration calculation for the protein.

        Equivalent to running `self.get_titration_from_propka`.

        Must run `self.fetch_pdb` first or specify an abosulute path to the PDB
        file in `self.pdb_location_absolute`.

        Raises:
            ValueError: If `titration_data` is not already stored and
                `pdb_location_absolute` is not set.

        Returns:
            structure.titration.TitrationData: A
                :class:`protein_structure.titration.TitrationData` object containing
                the titration values for residues.
        """
        return self.get_titration_from_propka()

    def get_titration_from_propka(self) -> structure.titration.TitrationData:
        """Fetches precomputed titration data for the protein, or computes it.

        Uses :func:`protein_structure.titration.calculate_titration_propka` if
        `self.titration_data` is not already stored.

        Must run `self.fetch_pdb` first or specify an abosulute path to the PDB
        file in `self.pdb_location_absolute`.

        Raises:
            ValueError: If `titration_data` is not already stored and
                `pdb_location_absolute` is not set.

        Returns:
            structure.titration.TitrationData: A
                :class:`protein_structure.titration.TitrationData` object containing
                the titration values for residues."""
        if self.titration_data:
            return self.titration_data

        if self.pdb_location_absolute:
            self.titration_data = structure.titration.calculate_titration_propka(
                self.pdb_location_absolute,
            )
            return self.titration_data
        else:
            raise ValueError(
                "Titration data not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def get_titration_from_pypka(self) -> structure.titration.TitrationData:
        """Fetches precomputed titration data for the protein, or computes it.

        Uses :func:`protein_structure.titration.calculate_titration_pypka` if
        `self.titration_data` is not already stored. Requires pypka to be
        installed, which has dependencies that are not FOSS. Please be sure to
        verify that you are legally allowed to use pypka.

        Must run `self.fetch_pdb` first or specify an abosulute path to the PDB
        file in `self.pdb_location_absolute`.

        Raises:
            ValueError: If `titration_data` is not already stored and
                `pdb_location_absolute` is not set. ImportError: If pypka is not
                installed.

        Returns:
            structure.titration.TitrationData: A
                :class:`protein_structure.titration.TitrationData` object containing
                the titration values forresidues."""

        if self.titration_data:
            return self.titration_data

        if self.pdb_location_absolute:
            self.titration_data = structure.titration.calculate_titration_pypka(
                self.pdb_location_absolute,
            )
            return self.titration_data
        else:
            raise ValueError(
                "Titration data not stored, and PDB location not set; use `fetch_pdb` first"
            )

    def get_titration_from_pkai(self) -> structure.titration.TitrationData:
        """Fetches precomputed titration data for the protein, or computes it.

        Uses :func:`protein_structure.titration.calculate_titration_pkai` if
        `self.titration_data` is not already stored. Requires pkai to be
        installed. Note that this method is a deep-learning model, not a
        physics-based calculation.

        Must run `self.fetch_pdb` first or specify an abosulute path to the PDB
        file in `self.pdb_location_absolute`.

        Raises:
            ValueError: If `titration_data` is not already stored and
                `pdb_location_absolute` is not set.

        Returns: structure.titration.TitrationData: A
            :class:`protein_structure.titration.TitrationData` object containing
                the titration values for residues."""
        if self.titration_data:
            return self.titration_data

        if self.pdb_location_absolute:
            self.titration_data = structure.titration.calculate_titration_pkai(
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
        """Split the protein into individual sites, recording values for each.

        Args:
            selected_aas: A set of amino acids letters to include in the output.
                If `None` (default), all amino acids will be included.
            selected_keys: A set of keys to include in the output. Possible keys
                can be viewed in `Protein.UNIPROT_SITE_PATTERNS_RECTIFIED`. If
                `None` (default), all keys for which this `Protein` object has
                data will be included.

        Returns:
            dict[str, list[Any]]: A dictionary mapping keys to lists of values.
                Each list is a parallel array of the same length as the protein
                sequence (after filtering out non-selected amino acids)."""
        if not selected_keys:
            selected_keys = set(self.data.keys()) - {"sequence"}

        site_keys = (
            set(x + "_sites" for x in Protein.UNIPROT_SITE_PATTERNS_RECTIFIED.keys())
            & selected_keys
        )
        other_keys = selected_keys - site_keys

        res: dict[str, list[Any]] = {k: [] for k in other_keys}
        for k in site_keys:
            res["is_" + k.removesuffix("_sites")] = []  # type: ignore
        res["residue_letter"] = []
        res["residue_number"] = []

        for index, site in enumerate(self.data["sequence"]):
            if selected_aas and site not in selected_aas:
                continue
            res["residue_letter"].append(site)
            res["residue_number"].append(index + 1)
            for k in site_keys:
                res["is_" + k.removesuffix("_sites")].append(  # type: ignore
                    (index + 1) in self.data[k]
                )
            for k in other_keys:
                res[k].append(self.data[k])  # will be the same for all sites

        return res

    def fetch_pdb(self, save_path: str | None = None, url: str | None = None) -> None:
        """Fetches the PDB file for the protein (from the AlphaFold database by default).

        Args:
            save_path (str | None, optional): The path to save the PDB file to.
                Defaults to `None`.
            url (str | None, optional): The URL to fetch the PDB file from.
                Defaults to `None`, in which case the AlphaFold database is used.

        Raises:
            Exception: If the response status code is not 200, meaning we could
                not fetch the PDB from the database."""
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

    def register_local_pdb(self, path_to_pdb_file: str | None = None) -> None:
        """Sets pdb file for protein object using local pdb file.

        Args:
            path_to_pdb_file (str | None, optional): Path to local PDB file.
                Defaults to `None`, in which case it assumes a file with 'entry'.pdb."""
        if not path_to_pdb_file:
            path_to_pdb_file = f"{self.data['entry']}.pdb"
        self.pdb_location_relative = path_to_pdb_file
        self.pdb_location_absolute = os.path.abspath(path_to_pdb_file)

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
        if "sequence" not in self.data:
            raise ValueError("Sequence entry not found in data")

        sequence = self.data["sequence"]

        return site <= len(sequence) and sequence[site - 1] == aa
