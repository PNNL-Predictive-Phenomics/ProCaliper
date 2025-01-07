import pandas as pd

import procaliper

"""
This script is an example that creates a pandas dataframe from a UniProt table.

In the input table, each row represents a UniProt entry (i.e., a protein).
"""


if __name__ == "__main__":
    # Read in the UniProt table as a pandas dataframe. This type of file can be
    # downloaded from the UniProt website using the search bar.
    df = pd.read_csv("examples/database_generation/sample_uniprot.tsv", sep="\t")

    sdf = None  # this will store the final dataframe

    # Iterate over the rows in the UniProt table. In principle, we don't need to
    # use pandas for this, but it is easier and doesn't slow things down too
    # much, even for human-proteome scale datasets (~20K rows).
    for _, row in df.iterrows():
        # This is the primary way we read in tabular data into ProCaliper.
        # Note that `procaliper` expects `row` to be a dictionary, but that a
        # pandas Series can be used as well (as here). Your type checker may
        # complain about this.
        protein = procaliper.Protein.from_uniprot_row(row)  # type: ignore

        # Download and save the structure. If the structure has already been
        # downloaded, use `protein.register_local_pdb` instead.
        protein.fetch_pdb(
            save_path=f"examples/database_generation/pdb_files/{protein.data['entry']}.pdb"
        )
        # Get a table of site-level metadata extracted from the pdb and store it in a dataframe.
        unravelled = pd.DataFrame(protein.unravel_sites())  # all data, all sites

        # Compute structure features using default arguements. See the documentation
        # for more details.
        sasa = pd.DataFrame(protein.get_sasa())
        charge = pd.DataFrame(protein.get_charge())
        titr = pd.DataFrame(protein.get_titration())

        # pLDDT is extracted directly from the PDB file using the
        # `get_confidence` method.
        conf = pd.DataFrame(protein.get_confidence(), columns=["pLDDT"])

        # Concatenate the metadata and structure features into a single dataframe
        row_df = pd.concat([unravelled, sasa, charge, conf, titr], axis=1)

        # Append the dataframe to the final dataframe
        if sdf is None:
            sdf = row_df
        else:
            sdf = pd.concat([sdf, row_df])

    assert sdf is not None

    # Reorder the columns
    first_columns = ["entry", "residue_letter", "residue_number"]
    new_column_order = first_columns + list(sdf.columns.difference(first_columns))
    sdf = sdf.reindex(columns=new_column_order)

    # Print the dataframe and its columns
    print(sdf)
    print(sdf.columns)

    # A very simple example of how to use the dataframe
    print(
        "Number of binding sites in sample (across all proteins):",
        sdf["binding"].sum(),
    )
