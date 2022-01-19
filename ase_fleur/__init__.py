# -*- coding: utf-8 -*-
"""
ASE package for interfacing with the FLEUR code
"""
from ase.calculators.calculator import register_calculator_class
from ase_fleur.calculator import Fleur, FleurProfile

# register calculator class (while there is no entry point in ase)
register_calculator_class("fleur", Fleur)

__all__ = ("Fleur", "FleurProfile")

__copyright__ = "Copyright (c), Forschungszentrum Jülich GmbH, IAS-1/PGI-1, Germany. All rights reserved."
__license__ = "MIT license, see LICENSE file."
__version__ = "0.0.1"
__authors__ = "The JuDFT team"
