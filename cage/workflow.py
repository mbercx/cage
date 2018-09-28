# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

from cage.cli.commands.setup import optimize, chainsetup, spheresetup
from cage.core import Cage

from fireworks import Firework, LaunchPad, PyTask, FWorker, \
    Workflow, ScriptTask

"""
Workflow setup for the cage calculations.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Jun 2018"

# Set up the Launchpad for the workflows
LAUNCHPAD = LaunchPad(host="ds221271.mlab.com", port=21271, name="cage",
                      username="mbercx", password="yosemite1")

RUN_NWCHEM_SCRIPT = "/g/g91/bercx1/local/scripts/job_workflow.sh"
RUN_NWCHEM_COMMAND = "srun -N1 -n36 /g/g91/bercx1/nwchem/nwchem-6.6/bin/LINUX64/nwchem"


def landscape_workflow(filename, cation, facets, operation, end_radii, nradii,
                       adensity):
    """
    Workflow to calculate the 2D landscape along a List of facets.

    Currently based upon the directory structure created by the setup
    script. This works fine, but perhaps isn't the cleanest solution.

    Args:
        filename:
        cation:
        facets:
        operation:
        end_radii:
        nradii:
        adensity:

    Returns:

    """
    # Load the cage molecule from the filename provided
    c = Cage.from_file(filename)
    molecule = c.composition.reduced_formula

    # Set up the calculation directories and input
    chainsetup(filename, cation, facets, operation, end_radii, nradii,
               adensity)
    chain_dir = "chain" + "_" + operation

    wf_list = []

    # For each edge
    for edge in [d for d in os.listdir(chain_dir) if "edge" in d]:

        fws_list = []

        # Set up the list of FireTasks, i.e. energy calculations:
        for geo in [d for d in os.listdir(os.path.join(chain_dir, edge))
                    if "geo" in d]:

            dir = os.path.abspath(os.path.join(chain_dir, edge, geo))
            nwchem_command = RUN_NWCHEM_COMMAND + " " \
                             + os.path.join(dir, "input") + " > " \
                             + os.path.join(dir, "result.out")

            fws_list.append(
                Firework(tasks= [ScriptTask.from_str(nwchem_command)],
                         name= operation + " " + geo)
            )

        workflow_name = "Landscape: " + cation + " on " + molecule + " (" + \
                        edge + ")"
        wf_list.append(
            Workflow(fireworks=fws_list, name=workflow_name)
        )

    # Set the workflows up to run from the launchpad
    LAUNCHPAD.bulk_add_wfs(wf_list)

def sphere_workflow(filename, cation, radius, density):
    """
    Set up a spherical landscape calculation of a specified radius.

    Args:
        filename:
        cation:
        radius:
        density:

    Returns:

    """

    c = Cage.from_file(filename)
    molecule = c.composition.reduced_formula

    # Set up the calculation directories and input
    sphere_dir = spheresetup(filename=filename,
                             cation=cation,
                             radius=radius,
                             density=density)

    fws_list = []

    # Set up the list of FireTasks, i.e. energy calculations:
    for geo in [d for d in os.listdir(sphere_dir) if "geo" in d]:

        dir = os.path.abspath(os.path.join(sphere_dir, geo))

        nwchem_command = RUN_NWCHEM_COMMAND + " " \
                         + os.path.join(dir, "input") + " > " \
                         + os.path.join(dir, "result.out")

        fws_list.append(
            Firework(tasks= [ScriptTask.from_str(nwchem_command)],
                     name= "energy " + geo)
        )

    workflow_name = "Landscape: " + cation + " on " + molecule + " (" + \
                     sphere_dir + ")"

    wf_list = [Workflow(fireworks=fws_list, name=workflow_name), ]

    LAUNCHPAD.bulk_add_wfs(wf_list)


def optimize_workflow(filename, charge=0):
    """
    Workflow for the optimization of a molecule

    Args:
        filename:
        charge:

    Returns:

    """
    # Load the cage molecule from the filename provided
    molecule = Cage.from_file(filename)

    # Set up optimization calculation
    optimize(filename, charge)

    optimize_dir = os.path.join(os.getcwd(), "optimize")

    optimize_command = RUN_NWCHEM_COMMAND + " " \
                       + os.path.join(optimize_dir, "input") + " > " \
                       + os.path.join(optimize_dir, "result.out")

    fw = Firework(tasks=[ScriptTask.from_str(optimize_command)],
                  name="Run Nwchem")

    LAUNCHPAD.add_wf(
        Workflow(fireworks=[fw],
                 name="Optimize " + molecule.composition.reduced_formula)
    )


def test_workflow():
    current_dir = os.getcwd()

    fw = Firework(ScriptTask.from_str("echo 'This worked!' >> " +
                                      current_dir + "/test"))

    LAUNCHPAD.add_wf(
        Workflow(fireworks=[fw],
                 name="test workflow")
    )
