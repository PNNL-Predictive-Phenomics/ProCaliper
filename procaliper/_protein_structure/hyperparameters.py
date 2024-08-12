HYPERPARAMETERS = {
    # Hyperparameters
    # for more information visit tutorial.md
    # General Parameters
    "path_to_dataframe": "input/10400_active_binding.csv",
    "dataframe_uniprot_col_name": "Entry",  # Column name of dataframe where uniprot ids are stored
    "dataframe_site_number_col_name": "site_number",  # Col name of dataframe where cys (residue) numbers are stored
    "path_to_save": "output/10400_features.csv",  # Where the data calculated from the pdbs will be saved to. Must end in .csv
    "path_to_pdb_db": "/Users/uffo617/Documents/PTM_ML/DatabaseScripts/1_CreateDatabase/1_GetPDBs/output/10400pdb",  # Path to location of pdb folder. pdb filenames must include their uniprot id.
    "files_per_save": 1,  # Indicates the number of files processed before each save. Recomended 1-10
    "output_file_name": "output/10400_active_binding_w_features.csv",
    # SASA Parameters
    "n_points": 100,  # Determines the accuracy of the sasa calculation. Higher values are more accurate but take more time to calculate. Recommended 100-1000
    # Charge Parameters
    "charge_method_used": "gasteiger",  # method of charge calculation, for all options check  https://open-babel.readthedocs.io/en/latest/Charges/charges.html
}
