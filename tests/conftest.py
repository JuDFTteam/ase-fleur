# -*- coding: utf-8 -*-
"""
Test configuration
"""
from ase_fleur.calculator import FleurProfile, Fleur
from ase.test.factories import factory as factory_dec, Factories, CalculatorInputs
from ase.utils import workdir

import pytest
from pathlib import Path


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


@pytest.fixture(scope="session", autouse=True)
def sessionlevel_testing_path():
    # We cd into a tempdir so tests and fixtures won't create files
    # elsewhere (e.g. in the unsuspecting user's directory).
    #
    # However we regard it as an error if the tests leave files there,
    # because they can access each others' files and hence are not
    # independent.  Therefore we want them to explicitly use the
    # "testdir" fixture which ensures that each has a clean directory.
    #
    # To prevent tests from writing files, we chmod the directory.
    # But if the tests are killed, we cannot clean it up and it will
    # disturb other pytest runs if we use the pytest tempdir factory.
    #
    # So we use the tempfile module for this temporary directory.
    import tempfile

    with tempfile.TemporaryDirectory(prefix="ase-test-workdir-") as tempdir:
        path = Path(tempdir)
        path.chmod(0o555)
        with workdir(path):
            yield path
        path.chmod(0o755)
