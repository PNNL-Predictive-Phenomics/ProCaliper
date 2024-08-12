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
    size = pd.DataFrame(protein.get_size())

    # entry info will be in unravelled
    sasa.drop(columns=["entry"], inplace=True)
    charge.drop(columns=["entry"], inplace=True)
    size.drop(columns=["entry"], inplace=True)

    unravelled = pd.DataFrame(protein.unravel_sites(selected_aas={"C"}))

    row_df = pd.concat([unravelled, sasa, charge, size], axis=1)

    if sdf is None:
        sdf = row_df
    else:
        sdf = pd.concat([sdf, row_df])

assert sdf is not None

print(sdf)
print(sdf.columns)
