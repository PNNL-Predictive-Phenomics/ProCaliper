![PyPI](https://img.shields.io/pypi/v/procaliper.svg)
[![Versions](https://img.shields.io/badge/Tested%20On%20Python-3.9,%203.10,%203.11,%203.12-blue.svg)](https://shields.io/)

# procaliper


A tool for fetching and organizing protein data in order to perform structure calculations.

# Installation

To install from [PyPI](https://pypi.org/project/procaliper/):

```shell
pip install procaliper
```
## extras
During installation, `procaliper` will also install `openbabel-wheel`. If you require specialized features of `openbabel` that are not available in this precompiled version, please run `pip uninstall openbabel-wheel` after `procaliper` is installed and provide your own version of the `openbabel` python library.

Optional feature dependencies can be installed as follows:

```shell
pip install procaliper[nglview,pka]
``` 

The `nglview` extra provides the ability to visualize protein structures in a graphical notebook environment.

The `pka` extra provides additional methods for computing disassociation constants (typically denoted $pK_a$). Note that installing this extra requires obtaining  a `DelPhi` license. Furthermore, these additional methods require an older version of `numpy` (version `1.26.4`) and `python` version between `3.9` and `3.11` to properly function. If the `pka` extra is not installed, `procaliper` will use [propka](https://github.com/jensengroup/propka) for $pK_a$ calculation.

# Basic Usage
A basic example is provided here. See the examples folder for further examples.

```python
import procaliper as pc

# Read in or download protein structure
protein = pc.Protein.from_uniprot_id("A0A0B4J2F0") # create a protein object from UniProt metadata using a UniProt ID
# alternatively, you can read in a protein from a UniProt Table using `Protein.from_uniprot_row`

# Download & save protein structure from AlphaFold
protein.fetch_pdb(save_path="A0A0B4J2F0.pdb") 
# Alternatively, `protein.register_local_pdb(file)` can be used to specify a previously downloaded pdb file

# Compute structure features
sasa = protein.get_sasa() # compute site-level SASA values
charge = protein.get_charge() # compute site-level charges
titr = protein.get_titration() # compute site-level titration (pKa) data
# The results from the above calculations are autmatically stored in the `protein` object.

# Get a table of site-level data
site_data = protein.unravel_sites() # returns a dictionary of lists; readable, e.g., by `pandas`
```

# Contributing

Contributions are welcome!

See [Contributing](https://github.com/PhenoMeters/ProCaliper/blob/main/CONTRIBUTING.md) for detailed instructions.

# License

This code is released under GPL v3. See our [License](https://github.com/PhenoMeters/ProCaliper/blob/main/LICENSE) for detailed information.

By default, `procaliper` does not have any dependencies that are not freely licensed. However, additional features can be installed that rely on restricted software. That is, please note that some *optional* dependencies are not FOSS. Specifically, the `pka` extra requires a `DelPhi` license. We encourage caution when using software that is not free and open source, especially for contributions to the scientific literature.

# Credits

This work is supported by the Predictive Phenomics Initiative at Pacific Northwest National Laboratory.