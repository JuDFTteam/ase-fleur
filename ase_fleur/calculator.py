# -*- coding: utf-8 -*-
"""
This module defines a calculator for the Fleur code starting from version v27
"""
from pathlib import Path
import warnings
import re

from masci_tools.io.fleurxmlmodifier import FleurXMLModifier
from masci_tools.io.parsers.fleur import outxml_parser

from ase.calculators.genericfileio import GenericFileIOCalculator, CalculatorTemplate

from ase_fleur.io import write_fleur_inpgen, read_fleur_outxml


class FleurProfile:
    """
    Profile for executing the Fleur code

    :param argv: arguments for the Fleur code
    :param inpgen_argv: arguments for the input generator for the Fleur code
    """

    def __init__(self, argv, inpgen_argv):
        self.argv = argv
        self.inpgen_argv = inpgen_argv

    def version(self) -> str:
        """
        Return the version string of the fleur code in this profile
        """
        from subprocess import check_output
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            with open(Path(td) / "err", "w", encoding="utf8") as err:
                out = check_output(self.argv + ["-info"], stderr=err, cwd=td).decode("utf-8")
        m = re.findall(r"^\s*MaX\-Release (.*)\(www\.max\-centre\.eu\)", out, flags=re.MULTILINE)
        if not m:
            raise ValueError(f"Could not retrieve version from output: {out}")
        return m[0].strip()

    def run(self, directory, outputfile, error_file):
        """
        Run Fleur in the given directory

        :param directory: path to the execution directory
        :param outputfile: path to the file for the stdout output
        :param error_file: path to the file for the stderr output
        """
        from subprocess import check_call

        with open(outputfile, "w", encoding="utf8") as fd:
            with open(error_file, "w", encoding="utf8") as ferr:
                check_call(self.argv, stdout=fd, stderr=ferr, cwd=directory)
        with open(outputfile) as f:
            print(f.read())

    def run_inpgen(self, directory, inputfile, outputfile, error_file):
        """
        Run inpgen in the given directory

        :param directory: path to the execution directory
        :param inputfile: path to the input file for the inpgen
        :param outputfile: path to the file for the stdout output
        :param error_file: path to the file for the stderr output
        """
        from subprocess import check_call

        with open(outputfile, "w", encoding="utf8") as fd:
            with open(error_file, "w", encoding="utf8") as ferr:
                check_call(self.inpgen_argv + ["-f", str(inputfile)], stdout=fd, stderr=ferr, cwd=directory)


class FleurTemplate(CalculatorTemplate):
    """
    Template defining a Fleur Calculation
    """

    def __init__(self, *, inpgen_profile):
        super().__init__(
            name="fleur",
            implemented_properties=("energy", "forces", "magmom", "magmoms", "efermi", "free_energy", "charges"),
        )
        self.stdout_file = "fleur.log"
        self.error_file = "error.log"
        self.output_file = "out.xml"
        self.inpgen_profile = inpgen_profile
        self.max_runs = 3
        self.iter_per_run = 30
        self.density_converged = 1e-6
        self.force_convergence = {"force_converged": 0.002, "qfix": 2, "forcealpha": 1.0, "forcemix": "straight"}

    def write_input(self, directory, atoms, parameters, properties):
        """
        Create Fleur inp.xml file from atoms object by calling the
        Fleur inpgen

        :param directory: path to the calculation directory
        :param atoms: ase.Atoms object to use
        :param parameters: Dict with inpgen parameters
                           Changes to be done after the inpgen was run
                           can be specified in the entry inpxml_changes
        """

        # Sketch
        # 1. Create inpgen input using the fleur IO format
        directory = Path(directory)
        directory.mkdir(exist_ok=True, parents=True)
        parameters = dict(parameters)
        inp_changes = parameters.pop("inpxml_changes", [])
        if "title" not in parameters:
            parameters["title"] = "Fleur inpgen input generated from ASE"
        else:
            if all(s not in parameters["title"] for s in ("inpgen", "input generator")):
                warnings.warn("inpgen or inputgenerator has to appear in the inpgen file title" "Added to the end")
                parameters["title"] += " (inpgen)"

        inputfile = directory / "fleur.in"
        write_fleur_inpgen(inputfile, atoms, parameters=parameters)

        # 2. Run inpgen
        self.execute_inpgen(directory, self.inpgen_profile, inputfile)

        # 3. Modify inp.xml according to set parameters
        fm = FleurXMLModifier()
        fm.set_inpchanges({"itmax": self.iter_per_run, "mindistance": self.density_converged})

        if "forces" in properties:
            fm.set_inpchanges(
                {
                    "force_converged": self.force_convergence["force_converged"],
                    "l_f": True,
                    "qfix": self.force_convergence["qfix"],
                    "forcealpha": self.force_convergence["forcealpha"],
                    "forcemix": self.force_convergence["forcemix"],
                }
            )

        # 4. ggf. make custom modifications using the FleurXMLmodifier
        if inp_changes:
            fm.add_task_list(inp_changes)

        xmltree, _ = fm.modify_xmlfile(directory / "inp.xml")
        xmltree.write(directory / "inp.xml", encoding="utf-8", pretty_print=True)

    def execute(self, directory, profile) -> None:
        """
        Execute Fleur multiple times until the calculation is either converged
        or a maximum number of iterations is reached

        :param directory: Path to the calculation directory
        :param profile: FleurProfile to use
        """
        converged = False
        run = 1
        while not converged and run <= self.max_runs:
            profile.run(directory, self.stdout_file, self.error_file)
            converged = self.check_convergence(directory)

    def execute_inpgen(self, directory, profile, inputfile) -> None:
        """
        Execute Fleur inpgen to create the Fleur inp.xml

        :param directory: Path to the calculation directory
        :param profile: FleurProfile to use
        :param inputfile: Path to the inputfile to use
        """
        profile.run_inpgen(directory, inputfile, self.stdout_file, self.error_file)

    def check_convergence(self, directory):
        """
        Check if the calculation is converged

        :param directory: Path to the calculation directory
        """
        fleur_results = outxml_parser(directory / self.output_file)

        MAGNETIC_DISTANCE_KEY = "overall_density_convergence"
        DISTANCE_KEY = "density_convergence"

        distance = fleur_results.get(MAGNETIC_DISTANCE_KEY, fleur_results.get(DISTANCE_KEY))
        if distance is None:
            raise ValueError("Could not find charge density distance in output file")

        return distance < self.density_converged

    def read_results(self, directory):
        """
        Read the calculation results from the produced out.xml

        :param directory: Path to the calculation directory
        """
        atoms = read_fleur_outxml(directory / self.output_file)
        return dict(atoms.calc.properties())


class Fleur(GenericFileIOCalculator):
    """
    Ase Calculator for FLEUR calculations
    """

    def __init__(self, *, profile=None, directory=".", **kwargs):

        if profile is None:
            profile = FleurProfile(["fleur"], ["inpgen"])

        super().__init__(
            template=FleurTemplate(inpgen_profile=profile), profile=profile, directory=directory, parameters=kwargs
        )
