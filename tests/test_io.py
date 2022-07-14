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

    param = 2.715
    assert all(atoms.symbols == "Si")
    assert atoms.positions == pytest.approx(
        np.array([[param / 4, param / 4, param / 4], [-param / 4, -param / 4, -param / 4]])
    )
    assert all(atoms.pbc)
    assert atoms.cell[:] == pytest.approx(np.array([[0.0, param, param], [param, 0.0, param], [param, param, 0.0]]))


def test_read_fleur_outxml():
    atoms = read(TEST_FILES_DIR / "out.xml")

    param = 2.715
    assert all(atoms.symbols == "Si")
    assert atoms.positions == pytest.approx(
        np.array([[param / 4, param / 4, param / 4], [-param / 4, -param / 4, -param / 4]])
    )
    assert all(atoms.pbc)
    assert atoms.cell[:] == pytest.approx(np.array([[0.0, param, param], [param, 0.0, param], [param, param, 0.0]]))

    assert dict(atoms.calc.properties()) == {
        "free_energy": pytest.approx(-15784.360931872383),
        "natoms": 2,
        "charges": pytest.approx(np.array([12.2309952, 12.2309952])),
        "fermi_level": pytest.approx(0.1848170588),
        "energy": pytest.approx(-15784.360931872383),
        "ibz_kpoints": pytest.approx(np.array([[0.1530809, 0.1530809, 0.1530809], [0.4592427, 0.1530809, 0.1530809]])),
        "kpoint_weights": pytest.approx(np.array([2 / 8, 6 / 8])),
        "nkpts": 2,
        "nbands": 19,
        "nspins": 1,
        "eigenvalues": pytest.approx(
            np.array(
                [
                    [
                        [
                            -4.70104467,
                            -4.70082255,
                            -3.08535988,
                            -3.08535988,
                            -3.08532042,
                            -3.08477431,
                            -3.08471498,
                            -3.08471498,
                            -0.19766336,
                            0.06831202,
                            0.18481706,
                            0.18481706,
                            0.28734007,
                            0.34124866,
                            0.34124866,
                            0.46788854,
                            0.47699595,
                            0.47699595,
                            0.54370476,
                        ],
                        [
                            -4.70098049,
                            -4.70088247,
                            -3.08541659,
                            -3.08540615,
                            -3.08517064,
                            -3.08490304,
                            -3.08466068,
                            -3.08465487,
                            -0.12710011,
                            -0.02193719,
                            0.08264381,
                            0.13155591,
                            0.26834718,
                            0.37615645,
                            0.42415907,
                            0.43011481,
                            0.58858562,
                            0.60013418,
                            0.65977508,
                        ],
                    ]
                ]
            )
        ),
    }


def test_read_fleur_outxml_magnetic():
    atoms = read(TEST_FILES_DIR / "out_magnetic.xml")

    param = 1.804494289
    assert all(atoms.symbols == "Fe2")
    assert atoms.positions == pytest.approx(np.array([[0, 0, 0]]))
    assert all(atoms.pbc)
    assert atoms.cell[:] == pytest.approx(np.array([[0, param, param], [param, 0, param], [param, param, 0]]))

    assert dict(atoms.calc.properties()) == {
        "free_energy": pytest.approx(-34642.48614156948),
        "natoms": 1,
        "charges": pytest.approx(np.array([25.0630064])),
        "fermi_level": pytest.approx(0.322012778),
        "energy": pytest.approx(-34642.48614156948),
        "ibz_kpoints": pytest.approx(np.array([[0.0, 0.0, 0.0]])),
        "kpoint_weights": pytest.approx(np.array([1])),
        "magmom": pytest.approx(3.9362688000000006),
        "magmoms": pytest.approx(np.array([3.93482625])),
        "nkpts": 1,
        "nbands": 14,
        "nspins": 1,  # 1 because the system is non-collinear
        "eigenvalues": pytest.approx(
            np.array(
                [
                    [
                        [
                            0.03188352,
                            0.07041465,
                            0.19753007,
                            0.19830856,
                            0.23418157,
                            0.24844237,
                            0.26124592,
                            0.32201278,
                            0.34943278,
                            0.34977853,
                            0.39773002,
                            0.43696835,
                            1.07972073,
                            1.09726426,
                        ]
                    ]
                ]
            )
        ),
    }
