# -*- coding: utf-8 -*-
"""
Tests of the fleur calculator class
"""
from ase.build import bulk
import pytest


def verify(calc):
    assert calc.get_fermi_level() is not None
    assert calc.get_ibz_k_points() is not None
    assert calc.get_eigenvalues(spin=0, kpt=0) is not None
    assert calc.get_number_of_spins() is not None
    assert calc.get_k_point_weights() is not None


@pytest.mark.calculator("fleur")
def test_main(factory):
    """
    Basic test of Fleur calculator
    """
    atoms = bulk("Si")
    atoms.calc = factory.calc()
    atoms.get_potential_energy()
    verify(atoms.calc)
