
import os

import pymatgen.io.nwchem as nwchem

from json import JSONDecodeError

"""
Small utility scripts that support the use of the cage package.

"""


def geo(output_file):
    """
    Quickly write out the initial and the final configuration of a nwchem
    optimization to xyz files.

    Args:
        output_file (str): Output file of the nwchem calculation.
    """
    data = nwchem.NwOutput(output_file).data[-1]
    data['molecules'][0].to(fmt='xyz', filename='initial_mol.xyz')
    data['molecules'][-1].to(fmt='xyz', filename='final_mol.xyz')


def check_calculation(output_file):
    """
    Script that checks if a calculation has completed successfully from the
    ouput file.
    """
    # TODO Add method of extracting data more quickly
    #TODO Make method independent of the name of the output file.

    if os.path.isdir(output_file):

        dir_list = [directory for directory in os.listdir(output_file)
                    if os.path.isdir(directory)]

        for directory in dir_list:

            file = os.path.join(directory, 'result.out')

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
            out = nwchem.NwOutput(output_file, fmt='json')
        except JSONDecodeError:
            try:
                out = nwchem.NwOutput(output_file)
            except:
                raise IOError('File not found.')

        try:
            error = False
            for data in out.data:
                if data['has_error']:
                    error = True

            print('File: ' + os.path.abspath(output_file))
            if out.data[-1]['task_time'] != 0:
                print('Calculation completed in ' + str(
                    out.data[-1]['task_time']) + 's')
            else:
                print('No timing information found. Calculation might not '
                      'have completed successfully.')

            print('Calculation has error: ' + str(error))

        except NameError:
            print("No data found in file!")


def process_output(output_file):
    """
    Process the results in an output file or all subdirectories in a directory.

    Args:
        output_file:

    Returns:

    """
    if os.path.isdir(output_file):

        dir_list = [directory for directory in os.listdir(output_file)
                    if os.path.isdir(directory)]

        for directory in dir_list:

            print("Processing output in " +
                  os.path.join(directory, 'result.out') +
                  "...")
            output = nwchem.NwOutput(os.path.join(directory, 'result.out'))

            try:
                error = False
                for data in output.data:
                    if data['has_error']:
                        error = True

                if error:
                    print("File: " + os.path.join(directory, 'result.out') +
                          " contains errors!")

                elif output.data[-1]['task_time'] == 0:
                    print('No timing information found in ' +
                          os.path.join(directory, 'result.out') + ".")

                else:
                    output.to_file(os.path.join(directory, 'data.json'))

            except NameError:

                print("No data found in file. ")

            except IndexError:

                print("Data is empty!")

    else:

        output_file = os.path.abspath(output_file)
        print('Processing output in ' + output_file)

        try:
            output = nwchem.NwOutput(output_file)
        except:
            raise IOError('Could not find proper nwchem output file.')

        try:
            error = False
            for data in output.data:
                if data['has_error']:
                    error = True

            if error:
                print("File: " + output_file + " contains errors!")

            elif output.data[-1]['task_time'] == 0:
                print('No timing information found in ' + output_file + ".")

            else:
                output.to_file(os.path.join(os.path.dirname(output_file),
                                            'data.json'))

        except NameError:

            print("No data found in file. ")

        except IndexError:

            print("Data is empty!")

        output.to_file(os.path.join(os.path.dirname(output_file), 'data.json'))


def search_and_reboot(dir_name):
    """

    Args:
        dir_name:

    Returns:

    """
    dir_list = [directory for directory in os.listdir(dir_name)
                if os.path.isdir(directory)]

    for dir_name in dir_list:

        print("Checking output in " + os.path.join(dir_name, 'result.out') +
              "...")
        output = nwchem.NwOutput(os.path.join(dir_name, 'result.out'))

        try:
            error = False
            for data in output.data:
                if data['has_error']:
                    error = True

            if error:
                print("File: " + os.path.join(dir_name, 'result.out') +
                      " contains errors! Simply rebooting is probably not "
                      "sufficient.")

            if output.data[-1]['task_time'] == 0:
                print('No timing information found in ' +
                      os.path.join(dir_name, 'result.out') +
                      '. Rebooting calculation...')
                os.system(
                    "sh -c 'cd " + dir_name + " && msub job_nwchem.sh '")

        except NameError:

            print("No data found in file. Rebooting calculation...")
            os.system("sh -c 'cd " + dir_name + " && msub job_nwchem.sh '")

        except IndexError:

            print("Data is empty! Rebooting Calculation...")
            os.system("sh -c 'cd " + dir_name + " && msub job_nwchem.sh '")
