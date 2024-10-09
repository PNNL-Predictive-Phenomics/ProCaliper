"""Top-level package for procaliper."""

__author__ = """AlphaMeter"""
__email__ = "song.feng@pnnl.gov"
__version__ = "0.1.0"

import procaliper.protein_structure as protein_structure
from procaliper._protein import Protein as Protein

__all__ = ["Protein", "protein_structure"]
