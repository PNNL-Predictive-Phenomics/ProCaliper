from procaliper import Protein

TEST_PDB_PATH = "tests/test_data/1nhz.pdb"
TEST_PDB_UNIPROT_ID = "P04150"

TEST_PDB_SEQ_LEN = 777

TEST_PDB_OBS_RES = 239

TEST_PDB_FIRST_RES = 530
TEST_PDB_LAST_RES = 776


def test_pdb_indexing() -> None:
    protein = Protein.from_uniprot_id(TEST_PDB_UNIPROT_ID)
    protein.register_local_pdb(path_to_pdb_file=TEST_PDB_PATH)

    assert protein.structure_index is not None
    assert protein.sequence_position_to_structure_index is not None

    assert protein.structure_index[0] == TEST_PDB_FIRST_RES
    assert protein.structure_index[-1] == TEST_PDB_LAST_RES
    assert len(protein.structure_index) == TEST_PDB_OBS_RES

    assert protein.sequence_position_to_structure_index[TEST_PDB_FIRST_RES] == 0
    assert (
        protein.sequence_position_to_structure_index[TEST_PDB_LAST_RES]
        == TEST_PDB_OBS_RES - 1
    )
    assert len(protein.sequence_position_to_structure_index) == TEST_PDB_OBS_RES

    assert len(protein.data["sequence"]) == TEST_PDB_SEQ_LEN

    assert len(protein.get_charge()["charge"]) == TEST_PDB_OBS_RES
    assert len(protein.get_sasa()["all_sasa_value"]) == TEST_PDB_OBS_RES
    assert (
        len(protein.get_cysteine_data()["min_dist_to_closest_sulfur"])
        == TEST_PDB_OBS_RES
    )
    assert len(protein.get_titration()["pKa"]) == TEST_PDB_OBS_RES
