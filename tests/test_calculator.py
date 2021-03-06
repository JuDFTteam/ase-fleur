# -*- coding: utf-8 -*-
"""
Tests of the fleur calculator class
"""
from ase.build import bulk
import pytest
import numpy as np


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


@pytest.mark.calculator("fleur")
def test_force(factory):
    """
    Basic test of Fleur calculator
    """
    atoms = bulk("Si")
    atoms.calc = factory.calc()
    assert atoms.get_forces() == pytest.approx(np.array([[0, 0, 0], [0, 0, 0]]))
    verify(atoms.calc)


@pytest.mark.calculator("fleur")
def test_version(factory):
    """
    Test of version parsing
    """
    assert factory.factory.version() == "6.0"
