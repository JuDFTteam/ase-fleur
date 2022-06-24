# -*- coding: utf-8 -*-
"""
Tests of the fleur calculator class
"""
from ase_fleur.calculator import FleurProfile, Fleur
from ase.test.factories import factory
from ase.build import bulk
import pytest


@factory("fleur")
class FleurFactory:
    """
    Factory for use in ase tests of the Calculator class
    """

    def __init__(self, executable, inpgen_executable):
        self.executable = executable
        self.inpgen_executable = inpgen_executable

    def _profile(self):
        return FleurProfile([self.executable], [self.inpgen_executable])

    def version(self):
        self._profile().version()

    def calc(self, **kwargs):
        return Fleur(profile=self._profile(), **kwargs)

    @classmethod
    def fromconfig(cls, config):
        return cls(config.executables["fleur"], config.executables["fleur_inpgen"])


def verify(calc):
    assert calc.get_fermi_level() is not None
    assert calc.get_ibz_k_points() is not None
    assert calc.get_eigenvalues(spin=0, kpt=0) is not None
    assert calc.get_number_of_spins() is not None
    assert calc.get_k_point_weights() is not None


@pytest.mark.calculator
def test_main(fleur_factory):
    """
    Basic test of Fleur calculator
    """
    atoms = bulk("Si")
    atoms.calc = fleur_factory.calc()
    atoms.get_potential_energy()
    verify(atoms.calc)
