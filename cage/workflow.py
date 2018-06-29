# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import subprocess
import shlex

from cage.cli.commands.setup import chainsetup
from cage.core import Cage

from custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, \
    UnconvergedErrorHandler
from custodian.vasp.jobs import VaspJob

from fireworks import Firework, LaunchPad, PyTask, FWorker, \
    Workflow
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

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

def run_nwchem(directory):
    """
    Method that simply runs VASP in the directory that is specified. Mainly
    used to set up a PyTask that uses outputs from PyTasks in previous
    FireWorks to run VASP in the appropriate directory.

    Args:
        directory: Absolute path to the directory in which VASP should be run.

    Returns:
        None

    """

    print(directory)

    command = RUN_NWCHEM_COMMAND + " " + os.path.abspath(directory) + \
              "/input > " + os.path.abspath(directory) + "/result.out"
    subprocess.Popen(command)

def landscape_workflow(filename, cation, facets, operation, end_radii, nradii,
               adensity):
    """

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

    # Set up the calculation directories and input
    chainsetup(filename, cation, facets, operation, end_radii, nradii,
               adensity)

    fw_list = []

    # For each edge
    for edge in [d for d in os.listdir("chain") if "edge" in d]:

        task_list = []

        # Set up the list of FireTasks, i.e. energy calculations:
        for geo in [d for d in os.listdir(os.path.join("chain", edge))
                    if "geo" in d]:

            task_list.append(
                PyTask(func="cage.workflow.run_nwchem",
                       args=[os.path.abspath(os.path.join("chain", edge,
                                                            geo))])
            )

        fw_list.append(
            Firework(tasks=task_list, name="Landscape: " + edge)
        )

    c = Cage.from_file(filename)
    molecule = c.composition.reduced_formula

    LAUNCHPAD.add_wf(
        Workflow(fireworks=fw_list,
                 name="Landscape: " + cation + " on " + molecule)
    )


