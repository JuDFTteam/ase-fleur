# -*- coding: utf-8 -*-
"""
Test configuration
"""
import ase_fleur  # Make sure the calculator class is registered
from ase.test.factories import Factories

import pytest

version = ase_fleur.__version__


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


@pytest.fixture(scope="session")
def generate_factories(pytestconfig):
    if pytestconfig.getoption("--calculator"):
        return Factories(["fleur"])
    return None
