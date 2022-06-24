# -*- coding: utf-8 -*-
"""
Test configuration
"""
from ase_fleur.calculator import FleurProfile, Fleur
from ase.test.factories import factory as factory_dec, Factories, CalculatorInputs

import pytest


@factory_dec("fleur")
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


def pytest_addoption(parser):
    """
    Add cmdline options
    """
    parser.addoption("--calculator", action="store_true", help="Run the calculator tests")


def pytest_configure(config):
    """
    Add markers
    """
    config.addinivalue_line("markers", "calculator: test running the FleurCalculator")


def pytest_collection_modifyitems(session, config, items):
    """
    Deselect the calculator tests if not specified
    """
    if not config.getoption("--calculator"):
        for item in items.copy():
            if "calculator" in item.keywords:
                items.remove(item)

        config.hook.pytest_deselected(items=items)


def pytest_generate_tests(metafunc):
    """
    Parametrize calculator factories
    """
    from ase.test.factories import parametrize_calculator_tests

    parametrize_calculator_tests(metafunc)


@pytest.fixture(scope="session", name="factories")
def factories_fixture(pytestconfig):
    if pytestconfig.getoption("--calculator"):
        return Factories(["fleur"])
    return Factories([])


@pytest.fixture(name="factory")
def factory_fixture(request, factories):
    name, kwargs = request.param
    if not factories.installed(name):
        pytest.skip(f"Not installed: {name}")
    if not factories.enabled(name):
        pytest.skip(f"Not enabled: {name}")
    factory = factories[name]
    return CalculatorInputs(factory, kwargs)
