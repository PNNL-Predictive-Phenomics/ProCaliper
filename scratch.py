from __future__ import annotations

import propka
import propka.run

from procaliper._protein_structure.titration import _state_from_pk
from procaliper.protein import Protein

TEST_HEADER = "Entry	Reviewed	Entry Name	Protein names	Gene Names	Organism	Length	Sequence	Active site	Binding site	DNA binding	Disulfide bond	Beta strand	Helix	Turn"
TEST_ROW = "A0A0B4J2F0	reviewed	PIOS1_HUMAN	Protein PIGBOS1 (PIGB opposite strand protein 1)	PIGBOS1	Homo sapiens (Human)	54	MFRRLTFAQLLFATVLGIAGGVYIFQPVFEQYAKDQKELKEKMQLVQESEEKKS							"

row_dict = {k: v for k, v in zip(TEST_HEADER.split("\t"), TEST_ROW.split("\t"))}
protein = Protein.from_uniprot_row(row_dict)
protein.fetch_pdb(save_path="tests/test_data/outputs/test_pdb.pdb")
protein.get_titration_from_pypka()


print("\n" * 10)
# print("-" * 100)
if protein.titration_data:
    for num, res, pka, stt in zip(
        protein.titration_data["residue_numbers"],
        protein.titration_data["residue_names"],
        protein.titration_data["pKs"],
        protein.titration_data["states"],
    ):
        print(num, res, pka, stt)


TEST_PATH = "tests/test_data/outputs/test_pdb.pdb"

mol = propka.run.single(TEST_PATH, optargs=["--quiet"], write_pka=False)
print("\n" * 2)
for group in mol.conformations["AVR"].groups:
    print(
        group.atom.res_num,
        group.atom.res_name,
        group.pka_value,
        _state_from_pk(group.pka_value),
    )

# for atom in mol.conformations["1A"].atoms:
#     if atom.group is not None:
#         print(atom.res_num, atom.group.residue_type, atom.group.pka_value)
