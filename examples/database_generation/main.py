import pandas as pd

import procaliper as am

df = pd.read_csv("examples/database_generation/sample_uniprot.tsv", sep="\t")  # type: ignore

sdf = None

for _, row in df.iterrows():  # type: ignore
    protein = am.Protein.from_uniprot_row(row)  # type: ignore

    protein.fetch_pdb(
        save_path=f"examples/database_generation/pdb_files/{protein.data['entry']}.pdb"
    )

    sasa = pd.DataFrame(protein.get_sasa())
    charge = pd.DataFrame(protein.get_charge())
    titr = pd.DataFrame(protein.get_titration())

    sasa.drop(columns=["residue_name", "residue_number"], inplace=True)
    charge.drop(columns=["residue_name", "residue_number"], inplace=True)
    titr.drop(columns=["residue_name", "residue_number"], inplace=True)

    conf = pd.DataFrame(protein.get_confidence(), columns=["pLDDT"])

    unravelled = pd.DataFrame(protein.unravel_sites())  # all data, all sites

    row_df = pd.concat([unravelled, sasa, charge, conf, titr], axis=1)

    if sdf is None:
        sdf = row_df
    else:
        sdf = pd.concat([sdf, row_df])

assert sdf is not None

print(sdf)
print(sdf.columns)

print("Number of binding sites in sample:", sdf["binding"].sum())
