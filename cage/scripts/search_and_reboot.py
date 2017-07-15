# Encoding: utf-8

import sys
import os
import pymatgen.io.nwchem as nw

dirname = sys.argv[1]

dir_list = [directory for directory in os.listdir(dirname)
            if os.path.isdir(directory)]

for directory in dir_list:

    print("Checking output in " + os.path.join(directory, 'result.out') +
          "...")
    output = nw.NwOutput(os.path.join(directory, 'result.out'))

    try:
        error = False
        for data in output.data:
            if data['has_error']:
                error = True

        if error:
            print("File: " + os.path.join(directory, 'result.out') +
                  " contains errors! Simply rebooting is probably not "
                  "sufficient.")

        if output.data[-1]['task_time'] == 0:

            print('No timing information found in ' +
                  os.path.join(directory, 'result.out') +
                  '. Rebooting calculation...')
            os.system("sh -c 'cd " + directory + " && msub job_nwchem.sh '")

    except NameError:

        print("No data found in file. Rebooting calculation...")
        os.system("sh -c 'cd " + directory + " && msub job_nwchem.sh '")

    except IndexError:

        print("Data is empty! Rebooting Calculation...")
        os.system("sh -c 'cd " + directory + " && msub job_nwchem.sh '")
