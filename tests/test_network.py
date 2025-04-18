from procaliper import Protein
from procaliper.network import regulatory_distance_network

INCOMPLETE_PDB_PATH = "tests/test_data/1nhz.pdb"
INCOMPLETE_PDB_UNIPROT_ID = "P04150"

COMPLETE_PDB_PATH = "tests/test_data/hsp90.pdb"
COMPLETE_PDB_UNIPROT_ID = "P07900"


def test_regulatory_distance_network() -> None:
    hsp90 = Protein.from_uniprot_id(COMPLETE_PDB_UNIPROT_ID)
    hsp90.register_local_pdb(path_to_pdb_file=COMPLETE_PDB_PATH)

    reg = regulatory_distance_network(hsp90)
    assert reg is not None


def test_regulatory_distance_network_incomplete() -> None:
    nhz = Protein.from_uniprot_id(INCOMPLETE_PDB_UNIPROT_ID)
    nhz.register_local_pdb(path_to_pdb_file=INCOMPLETE_PDB_PATH)

    reg = regulatory_distance_network(nhz)
    assert reg is not None
