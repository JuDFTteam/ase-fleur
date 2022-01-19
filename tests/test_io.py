# -*- coding: utf-8 -*-
"""
Tests of the io formats
"""
import pytest
from pathlib import Path
import numpy as np

from ase.io import read, write
from ase.build import bulk
from ase.calculators.calculator import compare_atoms

TEST_FILES_DIR = Path(__file__).resolve().parent / ".." / "test-files"


def test_fleur_inpgen_input_roundtrip(tmp_path):
    m1 = bulk("Ti")
    write(tmp_path / "fleur_inpgen", images=m1, format="fleur-inpgen")
    m2 = read(tmp_path / "fleur_inpgen", format="fleur-inpgen")

    # (How many decimals?)
    assert not compare_atoms(m1, m2, tol=1e-7)


def test_read_fleur_xml():
    atoms = read(TEST_FILES_DIR / "inp.xml")

    assert all(atoms.symbols == "Si")
    assert atoms.positions == pytest.approx(
        np.array([[1.28265213, 1.28265213, 1.28265213], [-1.28265213, -1.28265213, -1.28265213]])
    )
    assert all(atoms.pbc)
    assert atoms.cell[:] == pytest.approx(
        np.array([[0.0, 5.13060853, 5.13060853], [5.13060853, 0.0, 5.13060853], [5.13060853, 5.13060853, 0.0]])
    )


def test_read_fleur_outxml():
    atoms = read(TEST_FILES_DIR / "out.xml")

    assert all(atoms.symbols == "Si")
    assert atoms.positions == pytest.approx(
        np.array([[1.28265213, 1.28265213, 1.28265213], [-1.28265213, -1.28265213, -1.28265213]])
    )
    assert all(atoms.pbc)
    assert atoms.cell[:] == pytest.approx(
        np.array([[0.0, 5.13060853, 5.13060853], [5.13060853, 0.0, 5.13060853], [5.13060853, 5.13060853, 0.0]])
    )

    assert dict(atoms.calc.properties()) == {
        "free_energy": pytest.approx(-15784.360931872383),
        "natoms": 2,
        "charges": pytest.approx(np.array([12.2309952, 12.2309952])),
        "fermi_level": pytest.approx(0.1848170588),
        "energy": pytest.approx(-15784.360931872383),
    }
