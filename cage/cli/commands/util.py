
import os

import pymatgen as pmg
import pymatgen.io.nwchem as nwchem

from cage.landscape import LandscapeAnalyzer
from cage import Cage
from json import JSONDecodeError

"""
Small utility scripts that support the use of the cage package.

"""

# Elements to ignore for the surface facet determination
IGNORE = (pmg.Element('Li'), pmg.Element('Na'), pmg.Element('Mg'),
          pmg.Element('H'), pmg.Element('I'), pmg.Element('Br'),
          pmg.Element('Cl'), pmg.Element('F'))

OUTPUT_FILE = "result.out"
JOB_SCRIPT = "job_nwchem.sh"


def geo(output_file):
    """
    Quickly write out the initial and the final configuration of a nwchem
    optimization to xyz files.

    Args:
        output_file (str): Output file of the nwchem calculation.
    """
    data = nwchem.NwOutput(output_file).data[-1]
    data['molecules'][0].to(fmt='xyz', filename='initial_geo.xyz')
    data['molecules'][-1].to(fmt='xyz', filename='final_geo.xyz')


def energy(output_file):
    data = nwchem.NwOutput(output_file).data[-1]
    print("Total energy = " + str(data["energies"][-1]))


def check_calculation(output):
    """
    Script that checks if a calculation has completed successfully from the
    output file(s).

    Args:
        output (str): Output file or directory

    """

    # TODO Make this method recursively check subdirectories

    if os.path.isdir(output):

        dir_list = [directory for directory in os.listdir(output)
                    if os.path.isdir(directory)]

        for directory in dir_list:

            file = os.path.join(directory, OUTPUT_FILE)

            try:
                out = nwchem.NwOutput(file, fmt='json')
            except JSONDecodeError:
                try:
                    out = nwchem.NwOutput(file)
                except:
                    raise IOError('File not found.')
            except FileNotFoundError:
                print("No output file found in " + directory)

            try:
                error = False
                for data in out.data:
                    if data['has_error']:
                        error = True

                print('File: ' + os.path.abspath(file))
                if out.data[-1]['task_time'] != 0:
                    print('Calculation completed in ' + str(
                        out.data[-1]['task_time']) + 's')
                else:
                    print(
                        'No timing information found. Calculation might not '
                        'have completed successfully.')

                print('Calculation has error: ' + str(error))

            except NameError:
                print("No data found in file!")

    else:
        try:
            out = nwchem.NwOutput(output, fmt='json')
        except JSONDecodeError:
            try:
                out = nwchem.NwOutput(output)
            except:
                raise IOError('File not found.')

        try:
            error = False
            for data in out.data:
                if data['has_error']:
                    error = True

            print('File: ' + os.path.abspath(output))
            if out.data[-1]['task_time'] != 0:
                print('Calculation completed in ' + str(
                    out.data[-1]['task_time']) + 's')
            else:
                print('No timing information found. Calculation might not '
                      'have completed successfully.')

            print('Calculation has error: ' + str(error))

        except NameError:
            print("No data found in file!")


def gather_landscape(directory):
    """
    Gather the results from a landscape calculation.

    Args:
        directory (str): Directory which contains all the geometry directories
            that describe the landscape.
    """

    lands_analyzer = LandscapeAnalyzer.from_data(data_dir=directory)
    lands_analyzer.to(os.path.join(directory, "landscape.json"))


def process_output(output):
    """
    Process the results in an output file or all subdirectories in a directory.

    Args:
        output (str): File or directory that contains the output. If a
            directory is provided, the output has to be in the subdirectories
            of that directory.
    """
    if os.path.isdir(output):

        dir_list = [directory for directory in os.listdir(output)
                    if os.path.isdir(directory)]

        for directory in dir_list:

            print("Processing output in " +
                  os.path.join(directory, OUTPUT_FILE) +
                  "...")
            out = nwchem.NwOutput(os.path.join(directory, OUTPUT_FILE))

            try:
                error = False
                for output in out.data:
                    if output['has_error']:
                        error = True

                if error:
                    print("File: " + os.path.join(directory, OUTPUT_FILE) +
                          " contains errors!")

                elif out.data[-1]['task_time'] == 0:
                    print('No timing information found in ' +
                          os.path.join(directory, OUTPUT_FILE) + ".")

                else:
                    out.to_file(os.path.join(directory, 'data.json'))

            except NameError:

                print("No data found in file. ")

            except IndexError:

                print("Data is empty!")

    else:

        output = os.path.abspath(output)
        print('Processing output in ' + output)

        try:
            out = nwchem.NwOutput(output)
        except:
            raise IOError('Could not find proper nwchem output file.')

        try:
            error = False
            for output in out.data:
                if output['has_error']:
                    error = True

            if error:
                print("File: " + output + " contains errors!")

            elif out.data[-1]['task_time'] == 0:
                print('No timing information found in ' + output + ".")

            else:
                out.to_file(os.path.join(os.path.dirname(output),
                                         'data.json'))

        except NameError:

            print("No data found in file. ")

        except IndexError:

            print("Data is empty!")

        out.to_file(os.path.join(os.path.dirname(output), 'data.json'))


def search_and_reboot(directory):
    """
    Look to a directory for calculations in its subdirectories, and reboot them
    if they did not complete successfully.

    Args:
        directory (str):
    """
    dir_list = [subdir for subdir in os.listdir(directory)
                if os.path.isdir(subdir)]

    for directory in dir_list:

        print("Checking output in " + os.path.join(directory, OUTPUT_FILE) +
              "...")
        output = nwchem.NwOutput(os.path.join(directory, OUTPUT_FILE))

        try:
            error = False
            for data in output.data:
                if data['has_error']:
                    error = True

            if error:
                print("File: " + os.path.join(directory, OUTPUT_FILE) +
                      " contains errors! Simply rebooting is probably not "
                      "sufficient.")

            if output.data[-1]['task_time'] == 0:
                print('No timing information found in ' +
                      os.path.join(directory, OUTPUT_FILE) +
                      '. Rebooting calculation...')
                os.system(
                    "sh -c 'cd " + directory + " && msub " + JOB_SCRIPT + " '")

        except NameError:

            print("No data found in file. Rebooting calculation...")
            os.system("sh -c 'cd " + directory + " && msub " + JOB_SCRIPT
                      + " '")

        except IndexError:

            print("Data is empty! Rebooting Calculation...")
            os.system("sh -c 'cd " + directory + " && msub " + JOB_SCRIPT
                      + " '")


def visualize_facets(filename):
    """
    Visualize the facets of a molecule based on a structure file.

    Args:
        filename (str): Structure file of the molecule.
    """

    filename = str(os.path.abspath(filename))

    try:
        # Load the POSCAR into a Cage
        anion = Cage.from_poscar(filename)
    except ValueError:
        # If that fails, try other file formats supported by pymatgen
        anion = Cage.from_file(filename)

    anion.find_surface_facets(ignore=IGNORE)

    facet_filename = "".join(filename.split("/")[-1].split(".")[0:-1])\
                     + ".vesta"

    anion.visualize_facets(facet_filename, ignore=IGNORE)
