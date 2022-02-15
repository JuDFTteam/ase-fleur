# -*- coding: utf-8 -*-
"""Reads Fleur files.

Read structures from inpgen input files or inp.xml files

Built for fleur versions after 0.27

The package masci-tools (version 0.4.11 or greater) is required by this module:

    masci-tools - https://pypi.org/project/masci-tools/

"""
import io

import numpy as np

from ase.atoms import Atoms
from ase.calculators.singlepoint import SinglePointDFTCalculator
from ase.utils import writer
from ase.utils.plugins import ExternalIOFormat

from masci_tools.io.fleur_inpgen import write_inpgen_file, read_inpgen_file
from masci_tools.io.io_fleurxml import load_inpxml, load_outxml
from masci_tools.io.common_functions import convert_to_pystd
from masci_tools.util.parse_tasks_decorators import conversion_function
from masci_tools.util.xml.xml_getters import get_structure_data, get_kpoints_data
from masci_tools.util.schema_dict_util import eval_simple_xpath, get_number_of_nodes
from masci_tools.util.parse_utils import Conversion
from masci_tools.io.parsers.fleur import outxml_parser


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
    magic=b"*fleurInputVersion*",
)

outxml_format = ExternalIOFormat(
    desc="FLEUR output XML input file",
    code="1F",
    module="ase_fleur.io",
    magic=b"*fleurOutputVersion*",
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
    try:
        xmltree, schema_dict = load_outxml(fileobj)
    except (ValueError, IndexError):
        if isinstance(fileobj, io.IOBase):
            fileobj.seek(0)
        xmltree, schema_dict = load_inpxml(fileobj)

    atoms, cell, pbc = get_structure_data(xmltree, schema_dict, site_namedtuple=True, convert_to_angstroem=False)
    positions, symbols, _ = zip(*atoms)

    return Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell)


OUTXML_ADDITIONAL_TASKS = {
    "total_energy_ase": {
        "_minimal": True,
        "_conversions": [
            Conversion(
                name="convert_htr_to_ev",
                kwargs={
                    "name": "total_energy",
                    "converted_name": "total_energy_ev",
                },
            )
        ],
        "total_energy": {"parse_type": "attrib", "path_spec": {"name": "value", "tag_name": "totalEnergy"}},
    },
    "atom_charges": {
        "_conversions": ["calculate_total_charge_atoms"],
        "parsed_atom_charges": {
            "parse_type": "attrib",
            "path_spec": {"name": "total", "tag_name": "mtCharge", "contains": "valence"},
        },
        "parsed_corestates": {
            "parse_type": "allAttribs",
            "path_spec": {"name": "coreStates"},
            "kwargs": {
                "subtags": True,
            },
            "flat": False,
        },
    },
}
MAGMOMS_KEY = "magnetic_moments"
MAGMOM_KEY = "total_magnetic_moment_cell"
FORCES_KEY = "force_atoms"


@conversion_function
def calculate_total_charge_atoms(out_dict, logger):
    """
    Calculate the the total charge per atom

    :param out_dict: dict with the already parsed information
    """
    mt_charges = out_dict.pop("parsed_atom_charges", None)
    if mt_charges is None:
        if logger is not None:
            logger.warning("calculate_total_charge_atoms got None for Muffin-tin charges")
        return out_dict

    corestates = out_dict.pop("parsed_corestates", None)
    if corestates is None:
        if logger is not None:
            logger.warning("calculate_total_charge_atoms got None for corestates")
        return out_dict

    spins = max(corestates["spin"]) if isinstance(corestates["spin"], list) else corestates["spin"]

    # 1. Calculate the core charges per atomtype
    if not isinstance(corestates["state"], list):
        corestates["state"] = [corestates["state"]]
    if spins == 2:
        corestates["state"] = corestates["state"][: len(corestates["state"]) // 2]

    corecharges = []
    for states in corestates["state"]:
        weights = states["weight"]
        if isinstance(weights, list):
            weights = sum(weights)
        corecharges.append(weights)
    corecharges = np.array(corecharges)

    # 2. Calculate the mt charges
    mt_charges = mt_charges[0]
    if not isinstance(mt_charges, list):
        mt_charges = [mt_charges]

    if spins == 2 and len(mt_charges) % 2 != 0:
        if logger is not None:
            logger.warning("calculate_total_charge_atoms got spins=2 and odd number of mt charges")
        else:
            raise ValueError("calculate_total_charge_atoms got spins=2 and odd number of mt charges")
        return out_dict

    if spins == 2:
        mt_charges = [
            charge_up + charge_dn
            for charge_up, charge_dn in zip(mt_charges[: len(mt_charges) // 2], mt_charges[len(mt_charges) // 2 :])
        ]
    mt_charges = np.array(mt_charges)

    atom_charges = convert_to_pystd(mt_charges + corecharges)
    out_dict.setdefault("atom_charges", []).append(atom_charges)

    return out_dict


def read_fleur_outxml(fileobj, index=-1):
    """Reads structure and results from fleur out.xml file.

    Parameters
    ----------
    fileobj: file object
        File handle from which data should be read.

    Other parameters
    ----------------
    index: integer -1
        Not used in this implementation.
    """
    structure = read_fleur_xml(fileobj)

    if isinstance(fileobj, io.IOBase):
        fileobj.seek(0)
    xmltree, schema_dict = load_outxml(fileobj)
    kpoints, weights, cell, pbc = get_kpoints_data(xmltree, schema_dict, only_used=True, convert_to_angstroem=False)

    parser_warnings = {}
    results = outxml_parser(xmltree, parser_info_out=parser_warnings, additional_tasks=OUTXML_ADDITIONAL_TASKS)

    # The outxml_parser can only easily supply the charges per atom type for now
    # so we need to blow this up to per atom
    equivalent_atoms = []
    atom_groups = eval_simple_xpath(xmltree, schema_dict, "atomGroup", list_return=True)
    for group in atom_groups:
        if results["fleur_modes"]["film"]:
            equivalent_atoms.append(get_number_of_nodes(group, schema_dict, "filmpos"))
        else:
            equivalent_atoms.append(get_number_of_nodes(group, schema_dict, "relpos"))

    atomtype_charges = results.get("atom_charges")
    atom_charges = []
    if atomtype_charges:
        for charge, n_equiv in zip(atomtype_charges, equivalent_atoms):
            atom_charges.extend([charge] * n_equiv)
    atom_charges = np.array(atom_charges)

    results_dict = {
        "efermi": results.get("fermi_energy"),
        "free_energy": results.get("energy"),
        "energy": results.get("total_energy_ev"),
        "charges": atom_charges,
    }

    if MAGMOM_KEY in results:
        results_dict["magmoms"] = np.array(results[MAGMOMS_KEY])
        results_dict["magmom"] = results[MAGMOM_KEY]

    if FORCES_KEY in results:
        results_dict["forces"] = np.array([force for _, force in results[FORCES_KEY]])

    calc = SinglePointDFTCalculator(structure, **results_dict)
    structure.calc = calc
    return structure


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
