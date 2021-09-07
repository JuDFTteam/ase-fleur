# -*- coding: utf-8 -*-
"""Reads Fleur files.

Read structures from inpgen input files or inp.xml files

Built for fleur versions after 0.27

The package masci-tools (version 0.4.11 or greater) is required by this module:

    masci-tools - https://pypi.org/project/masci-tools/

"""
from ase.atoms import Atoms
from ase.utils import writer
from ase.utils.plugins import ExternalIOFormat
from masci_tools.io.fleur_inpgen import write_inpgen_file, read_inpgen_file
from masci_tools.io.io_fleurxml import load_inpxml
from masci_tools.util.xml.xml_getters import get_structure_data


inpgen_format = ExternalIOFormat(
    desc="FLEUR inpgen input file",
    code="1F",
    module="ase_fleur.io",
    magic=[b"*inpgen *", b"*input generator *"],
    glob="inp_*",
)

xml_format = ExternalIOFormat(
    desc="FLEUR XML input file",
    code="1F",
    module="ase_fleur.io",
    magic=b"*fleurInputVersion *",
)


def read_fleur_inpgen(fileobj, index=-1):
    """Reads structure from fleur inpgen file.

    Parameters
    ----------
    fileobj: file object
        File handle from which data should be read.

    Other parameters
    ----------------
    index: integer -1
        Not used in this implementation.
    """

    # Last parameter is lapw parameters and is not used here
    cell, atoms, pbc, _ = read_inpgen_file(fileobj, convert_to_angstroem=False)
    positions, symbols, _ = zip(*atoms)

    return Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell)


def read_fleur_xml(fileobj, index=-1):
    """Reads structure from fleur xml file.

    Parameters
    ----------
    fileobj: file object
        File handle from which data should be read.

    Other parameters
    ----------------
    index: integer -1
        Not used in this implementation.
    """

    xmltree, schema_dict = load_inpxml(fileobj)
    atoms, cell, pbc = get_structure_data(xmltree, schema_dict, site_namedtuple=True, convert_to_angstroem=False)
    positions, symbols, _ = zip(*atoms)

    return Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell)


@writer
def write_fleur_inpgen(fileobj, atoms, parameters=None, **kwargs):
    """writes fleur input structure in inpgen input file

    Parameters
    ----------
    filename : str
        Name of file to which data should be written.
    images : Atom Object or List of Atoms objects
        This function will write the first Atoms object to file.

    Returns
    -------
    """

    atom_sites = [(atom.position, atom.symbol, atom.symbol) for atom in atoms]

    write_inpgen_file(
        atoms.cell,
        atom_sites,
        pbc=atoms.pbc,
        input_params=parameters,
        file=fileobj,
        convert_from_angstroem=False,
        **kwargs,
    )
