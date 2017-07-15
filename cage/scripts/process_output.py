# Encoding: utf-8

import sys
import os

import pymatgen.io.nwchem as nwchem

if os.path.isdir(sys.argv[1]):

    dirname = sys.argv[1]

    dir_list = [directory for directory in os.listdir(dirname)
                if os.path.isdir(directory)]

    for directory in dir_list:

        print("Processing output in " + os.path.join(directory, 'result.out') +
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

    output_file = sys.argv[1]
    dir_name = os.path.dirname(output_file)
    print('Processing output in ' + output_file)

    try:
        output = nwchem.NwOutput(output_file)
    except:
        raise IOError('Could not find proper nwchem output file.')

    output.to_file(os.path.join(dir_name, 'data.json'))



