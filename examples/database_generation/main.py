import pandas as pd

import procaliper

df = pd.read_csv("examples/database_generation/sample_uniprot.tsv", sep="\t")

sdf = None

for _, row in df.iterrows():
    protein = procaliper.Protein.from_uniprot_row(row)  # type: ignore

    protein.fetch_pdb(
        save_path=f"examples/database_generation/pdb_files/{protein.data['entry']}.pdb"
    )

    sasa = pd.DataFrame(protein.get_sasa())
    charge = pd.DataFrame(protein.get_charge())
    titr = pd.DataFrame(protein.get_titration())

    conf = pd.DataFrame(protein.get_confidence(), columns=["pLDDT"])

    unravelled = pd.DataFrame(protein.unravel_sites())  # all data, all sites

    row_df = pd.concat([unravelled, sasa, charge, conf, titr], axis=1)

    if sdf is None:
        sdf = row_df
    else:
        sdf = pd.concat([sdf, row_df])

assert sdf is not None

first_columns = ["entry", "residue_letter", "residue_number"]
new_column_order = first_columns + list(sdf.columns.difference(first_columns))

sdf = sdf.reindex(columns=new_column_order)

print(sdf)
print(sdf.columns)

print("Number of binding sites in sample:", sdf["binding"].sum())
