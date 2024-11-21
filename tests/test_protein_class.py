import pandas as pd

from procaliper import Protein

TEST_DATA_PATH = (
    "tests/test_data/uniprotkb_Human_AND_model_organism_9606_2024_08_07.tsv"
)


def test_read_uniprot_row() -> None:
    COMPARISON_ENTRY_1 = "A0A0B4J2F0"
    COMPARISON_SEQUENCE_1 = "MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS"

    COMPARISON_ENTRY_2 = "A0A0K2S4Q6"
    COMPARISON_DISULFIDE_2 = list(range(43, 111 + 1))  # 43..111, expand = True
    # STRAND 79..81; /evidence="ECO:0007829|PDB:7EMF"	HELIX 83..86; /evidence="ECO:0007829|PDB:7EMF"; HELIX 90..96; /evidence="ECO:0007829|PDB:7EMF"; HELIX 112..116; /evidence="ECO:0007829|PDB:7EMF"; HELIX 128..138; /evidence="ECO:0007829|PDB:7EMF"; HELIX 147..152; /evidence="ECO:0007829|PDB:7EMF"
    COMPARISON_ENTRY_3 = "A0JLT2"
    COMPARISON_STRAND = [79, 80, 81]
    COMPARISON_HELIX = (
        [83, 84, 85, 86]
        + [90, 91, 92, 93, 94, 95, 96]
        + [112, 113, 114, 115, 116]
        + [128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138]
        + [147, 148, 149, 150, 151, 152]
    )
    COMPARISON_TURN = [117, 118, 119]

    df = pd.read_csv(  # type: ignore
        TEST_DATA_PATH,
        sep="\t",
        nrows=20,
    )

    assert df is not None

    for _, row in df.iterrows():  # type: ignore
        protein = Protein.from_uniprot_row(row)  # type: ignore

        assert protein is not None

        if protein.data["entry"] == COMPARISON_ENTRY_1:
            assert protein.data["sequence"] == COMPARISON_SEQUENCE_1

        if protein.data["entry"] == COMPARISON_ENTRY_2:
            assert protein.site_annotations is not None
            dbonds = protein.site_annotations.disulfide_bond
            assert dbonds is not None
            left = [i + 1 for i, v in enumerate(dbonds) if v]
            assert left == COMPARISON_DISULFIDE_2

        if protein.data["entry"] == COMPARISON_ENTRY_3:
            assert protein.site_annotations is not None
            strand = protein.site_annotations.beta_strand
            assert strand is not None
            left = [i + 1 for i, v in enumerate(strand) if v]
            assert left == COMPARISON_STRAND

            helix = protein.site_annotations.helix
            assert helix is not None
            left = [i + 1 for i, v in enumerate(helix) if v]
            assert left == COMPARISON_HELIX

            turn = protein.site_annotations.turn
            assert turn is not None
            left = [i + 1 for i, v in enumerate(turn) if v]
            assert left == COMPARISON_TURN


def test_unravel():
    TEST_HEADER = "Entry	Reviewed	Entry Name	Protein names	Gene Names	Organism	Length	Sequence	Active site	Binding site	DNA binding	Disulfide bond	Beta strand	Helix	Turn"
    TEST_ROW = "A0A0B4J2F0	reviewed	PIOS1_HUMAN	Protein PIGBOS1 (PIGB opposite strand protein 1)	PIGBOS1	Homo sapiens (Human)	54	MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS							"

    row_dict = {k: v for k, v in zip(TEST_HEADER.split("\t"), TEST_ROW.split("\t"))}
    protein = Protein.from_uniprot_row(row_dict)
    print(protein.data.keys())
    unravelled = protein.unravel_sites(
        selected_aas={"M"},
        selected_keys={"entry", "turn", "residue_letter", "residue_number"},
    )

    expected = {
        "entry": ["A0A0B4J2F0", "A0A0B4J2F0"],
        "turn": [False, False],
        "residue_letter": ["M", "M"],
        "residue_number": [1, 43],
    }

    assert unravelled == expected


def test_fetch_pdb():
    TEST_HEADER = "Entry	Reviewed	Entry Name	Protein names	Gene Names	Organism	Length	Sequence	Active site	Binding site	DNA binding	Disulfide bond	Beta strand	Helix	Turn"
    TEST_ROW = "A0A0B4J2F0	reviewed	PIOS1_HUMAN	Protein PIGBOS1 (PIGB opposite strand protein 1)	PIGBOS1	Homo sapiens (Human)	54	MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS							"

    row_dict = {k: v for k, v in zip(TEST_HEADER.split("\t"), TEST_ROW.split("\t"))}
    protein = Protein.from_uniprot_row(row_dict)
    protein.fetch_pdb(save_path="tests/test_data/outputs/test_pdb.pdb")


def test_structure_run_only():
    TEST_HEADER = "Entry	Reviewed	Entry Name	Protein names	Gene Names	Organism	Length	Sequence	Active site	Binding site	DNA binding	Disulfide bond	Beta strand	Helix	Turn"
    TEST_ROW = "A0A0B4J2F0	reviewed	PIOS1_HUMAN	Protein PIGBOS1 (PIGB opposite strand protein 1)	PIGBOS1	Homo sapiens (Human)	54	MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS							"

    row_dict = {k: v for k, v in zip(TEST_HEADER.split("\t"), TEST_ROW.split("\t"))}
    protein = Protein.from_uniprot_row(row_dict)
    protein.fetch_pdb(save_path="tests/test_data/outputs/test_pdb.pdb")

    protein.get_confidence()

    protein.get_charge()
    protein.get_sasa()
    protein.get_cysteine_data()

    protein.get_titration()
    protein.get_titration_from_propka()

    try:
        protein.get_titration_from_pypka()  # optional dependency
    except ImportError:
        pass
    try:
        protein.get_titration_from_pkai()  # optional dependency
    except ImportError:
        pass


def test_residue_data_frame_run_only():
    TEST_HEADER = "Entry	Reviewed	Entry Name	Protein names	Gene Names	Organism	Length	Sequence	Active site	Binding site	DNA binding	Disulfide bond	Beta strand	Helix	Turn"
    TEST_ROW = "A0A0B4J2F0	reviewed	PIOS1_HUMAN	Protein PIGBOS1 (PIGB opposite strand protein 1)	PIGBOS1	Homo sapiens (Human)	54	MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS							"

    row_dict = {k: v for k, v in zip(TEST_HEADER.split("\t"), TEST_ROW.split("\t"))}
    protein = Protein.from_uniprot_row(row_dict)
    protein.fetch_pdb(save_path="tests/test_data/outputs/test_pdb.pdb")

    protein.residue_data_frame()


def test_uniprot_api():
    df = pd.read_csv(  # type: ignore
        TEST_DATA_PATH,
        sep="\t",
        nrows=5,
    )

    ids: list[str] = df["Entry"].to_list()

    print(Protein.from_uniprot_id(ids[0]).data)
    print(Protein.from_uniprot_row(df.iloc[0].to_dict()).data)  # type: ignore

    assert Protein.from_uniprot_id(ids[0]) == Protein.from_uniprot_row(
        df.iloc[0].to_dict()  # type: ignore
    )

    assert Protein.list_from_uniprot_ids(ids) == [
        Protein.from_uniprot_row(row.to_dict())  # type: ignore
        for _, row in df.iterrows()  # type: ignore
    ]
