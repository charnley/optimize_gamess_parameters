#!/usr/bin/env python

import numpy as np
import re

import gamess
import mopac

def get_coordinates_xyz(filename):
    """
    Get coordinates from filename and return a vectorset with all the
    coordinates, in XYZ format.

    Parameters
    ----------
    filename : string
        Filename to read

    Returns
    -------
    atoms : list
        List of atomic types
    V : array
        (N,3) where N is number of atoms

    """

    f = open(filename, 'r')
    V = list()
    atoms = list()
    n_atoms = 0

    # Read the first line to obtain the number of atoms to read
    try:
        n_atoms = int(f.readline())
    except ValueError:
        print("Could not obtain the number of atoms in the .xyz file. "+filename)
        return None

    # Skip the title line
    f.readline()

    # Use the number of atoms to not read beyond the end of a file
    for lines_read, line in enumerate(f):

        if lines_read == n_atoms:
            break

        atom = re.findall(r'[a-zA-Z]+', line)[0]
        # atom = atom.upper()

        numbers = re.findall(r'[-]?\d+\.\d*(?:[Ee][-\+]\d+)?', line)
        numbers = [float(number) for number in numbers]

        # The numbers are not valid unless we obtain exacly three
        if len(numbers) == 3:
            V.append(np.array(numbers))
            atoms.append(atom)
        else:
            exit("Reading the .xyz file failed in line {0}. Please check the format.".format(lines_read + 2))

    f.close()
    atoms = np.array(atoms)
    V = np.array(V)
    return atoms, V



def get_xyzs(name_list, opt_dir):

    if not opt_dir[-1] == "/": opt_dir += "/"

    xyz_list = {}

    for name in name_list:
        atoms, coord = get_coordinates_xyz(opt_dir + name + ".xyz")
        xyz_list[name] = [atoms, coord]

    return xyz_list



def get_relative_jobs(filename, solvents, absjobs):

    absdb = {}
    for job in absjobs:

        name = job['file']
        solvent_name = job['solvent']
        jobid = job['id']

        print name

        if name not in absdb.keys(): absdb[name] = {}

        absdb[name][solvent_name] = jobid


    f = open(filename)

    header = f.next()
    header = header.split()
    n_cols = len(header)

    no_idx = header.index("No.")
    fileid_idx = header.index("FileHandle")
    solvent_idx = header.index("Solvent")
    type_idx = header.index("type")
    charge_idx = header.index("Charge")
    deltagsolv_idx = header.index("DeltaGsolv")
    eps_idx = header.index("eps")

    jobs = []
    molecules = {}

    for line in f:
        line = line.split()

        if line[type_idx] == "abs": continue

        job = {}
        job['id'] = int(line[no_idx])
        job['file'] = line[fileid_idx]
        job['solvent'] = line[solvent_idx].split("-")
        job['delta_g_solv'] = float(line[deltagsolv_idx])
        job['eps'] = float(line[eps_idx])
        job['charge'] = int(line[charge_idx])

        jobs.append(job)

        print job['solvent'], job['file'], job['delta_g_solv'], job['id']

        if job['file'] not in absdb.keys():
            print job['file'], "no way"
            continue

        print absdb[job['file']]


        quit()

    quit()

    return jobs



def get_full_database(filename):

    f = open(filename)

    header = f.next()
    header = header.split()
    n_cols = len(header)

    no_idx = header.index("No.")
    fileid_idx = header.index("FileHandle")
    solvent_idx = header.index("Solvent")
    type_idx = header.index("type")
    charge_idx = header.index("Charge")
    deltagsolv_idx = header.index("DeltaGsolv")
    eps_idx = header.index("eps")

    molecules = {}
    abs_jobs = []
    rel_jobs = []
    molecule_names = []
    solvents = {}

    for line in f:
        line = line.split()

        isabs = line[type_idx] == "abs"

        job = {}
        job['id'] = int(line[no_idx])
        job['file'] = line[fileid_idx]
        job['solvent'] = line[solvent_idx]
        job['delta_g_solv'] = float(line[deltagsolv_idx])
        job['eps'] = float(line[eps_idx])
        job['charge'] = int(line[charge_idx])

        # molecules[line[fileid_idx]] = int(line[charge_idx])

        molecule_names.append(line[fileid_idx])

        if isabs:
            abs_jobs.append(job)
            # eps db
            solvents[job['eps']] = job['solvent']
        else:
            rel_jobs.append(job)


    molecule_names = np.unique(molecule_names)

    return solvents, molecule_names, rel_jobs, abs_jobs


def print_result(jobs):

    print "id, eps, charge, free"
    for job in jobs:
        print ", ".join([str(job['id']), str(job['eps']), str(job['charge']), job['file'], str(job['delta_g_solv'])])



    return


def main():

    import argparse
    import sys

    description = ""

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options]',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-q', '--program_type', action='store', default="gamess", help='Input file program type')

    parser.add_argument('-p', '--parameters', action='store', nargs='+', default='', help='CSV file containing Radii information')

    parser.add_argument('-e', '--header', action='store', default='', help='GAMESS header file')

    parser.add_argument('-m', '--mnsoldb', action='store', default=False, help='MNSOL database file')

    parser.add_argument('-f', '--xyz_folder', action='store', default=False, help='Folder of gas optimized XYZ files')

    parser.add_argument('-x', '--exam_results', action='store_true', default=False, help='Print the results')

    args = parser.parse_args()


    # Get header
    if args.header:
        with open(args.header, 'r') as f:
            header = f.read()


    # Read radii parameters
    parameters = None
    if args.parameters:

        parameters = np.zeros((107))

        for parameter_file in args.parameters:
            with open(parameter_file) as f:
                parameter_content = f.read()
                parameter_content = parameter_content.split("\n")
                parameter_content = list(filter(None, parameter_content))

                if ":" in parameter_content[0]:
                    # Dictionary
                    for line in parameter_content:
                        line = line.split(":")
                        idx = int(line[0])
                        val = float(line[1])

                        parameters[idx-1] = val

                else:
                    # just a plain list
                    for i, line in enumerate(parameter_content):
                        parameters[i] = float(line)


    # Load mnsol database and XYZ structures
    solvents, molecule_names, rel_jobs, abs_jobs = get_full_database(args.mnsoldb)

    # relative_jobs = get_relative_jobs(args.mnsoldb, solvents, jobs)

    if args.exam_results:

        # keys = solvents.keys()
        # keys.sort()
        #
        # for key in keys:
        #     print key, solvents[key]

        print_result(jobs)


        quit()


    molecule_xyzs = get_xyzs(molecule_names, args.xyz_folder)

    if args.program_type == "gamess":
        make_input_file = gamess.make_input_file
        ext = "inp"
    elif args.program_type == "mopac":
        make_input_file = mopac.make_input_file
        ext = "mop"

    for job in abs_jobs:

        continue

        inpfile = str(job['id']) + "." + ext
        atoms, coord = molecule_xyzs[job['file']]
        charge = job['charge']
        eps = job['eps']

        print inpfile, eps, job['file']

        header_i = header.replace('REPLACEEPS', str(eps))

        with open(inpfile, 'w') as f:
            f.write(make_input_file(atoms, coord, header_i, charge, parameters))

    for job in rel_jobs:

        inpfile = str(job['id']) + ".inp"
        atoms, coord = molecule_xyzs[job['file']]
        charge = job['charge']
        eps = job['eps']

        print inpfile, eps, job['file']

        header_i = header.replace('REPLACEEPS', str(eps))

        # with open(inpfile, 'w') as f:
        #     f.write(make_input_file(atoms, coord, header_i, charge, parameters))



if __name__ == "__main__":
    main()

