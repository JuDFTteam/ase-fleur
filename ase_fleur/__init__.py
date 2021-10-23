# -*- coding: utf-8 -*-
"""
ASE package for interfacing with the FLEUR code
"""
from ase.calculators.calculator import register_calculator_class
from ase_fleur.calculator import Fleur, FleurProfile

# register calculator class (while there is no entry point in ase)
register_calculator_class("fleur", Fleur)

__all__ = ("Fleur", "FleurProfile")
__version__ = "0.0.1"
