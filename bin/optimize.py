#!/usr/bin/env python

import numpy as np
import scipy.optimize as opt
from copy import deepcopy
import numpy.linalg as linalg

import gamess as gms
import nodes
from shell import shell

from time import time

n_atom_types = 107


__XYZDIR__ = "/home/charnley/dev/2017-mnsol-data/jobs/xyz/" # local
__XYZDIR__ = "/home/charnley/dev/2017-mnsol/jobs/xyz/" # sunray
__XYZDIR__ = "/home/charnley/dev/2017-mnsol/jobs/xyz-gamess-pm6/" # sunray, pm6 gas


def rmsd(V, W):
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += (v - w)**2.0
    return np.sqrt(rmsd/N)


def main():

    import argparse
    import sys
    import os

    description = ""

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options]',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-p', '--parameters', action='store', nargs='+', default='', help='')
    parser.add_argument('-l', '--xyz_list', action='store', default='', help='list of molecules')
    parser.add_argument('-d', '--reference', action='store', default='', help='Reference list')
    parser.add_argument('-e', '--header', action='store', default='', help='Header file')
    parser.add_argument('-m', '--mask', action='store', default=None, nargs='+', type=str, help='Only use atoms from this subset')

    parser.add_argument('-w', '--double_workers', action='store_true', default=False, help='Double the amount of workers (usefull for SQM methods)')

    args = parser.parse_args()

    try:
        os.environ["SLURM_JOB_ID"]
        nodelist, workers, jobid = nodes.get_slurm_information()
        workers = 8
        chunksize = 90

        print "Submitted to:", nodelist, "on", workers, "workers. id:", jobid

    except KeyError:
        nodelist = ["node634", "node678", "node637", "node662"]
        workers = 8
        chunksize = 90
        print "Could not find SLURM enviroment", nodelist, "workers=",workers, "chunksize=",chunksize

    if args.double_workers:
        workers *= 2

    slave_master = nodes.hostname()

    # define gamess header
    if args.header:
        header_type = args.header.split("/")
        header_type = header_type[-1]

        print "header:", header_type

        __XYZDIR__ = "/home/charnley/dev/2017-pcm-parameters/xyz/"+header_type+"/"
        with open(args.header, 'r') as myfile:
                header = myfile.read()
    else:
        header = """

 $system
   mwords=250
 $end

 $basis
    gbasis=PM6
 $end

 $contrl
    scftyp=RHF
    icharg=0
    runtyp=energy
 $end

 $scf
    ! less output
    npunch=1
 $end

 $pcm
    solvnt=WATER
    smd=.t.
 $end

"""

    if args.reference:
        expe_values = args.reference
    else:
        # Default value for sunray
        expe_values = "/home/charnley/dev/2017-pcm-parameters/optimize_gamess_parameters/reference_water.csv"

    exam_subset = args.xyz_list # list of moleculenames

    # load in all charges
    # load in experimental values (dict)
    charges = {}
    references = {}
    f = open(expe_values, 'r')
    for line in f:
        line = line.split(";")
        charge = int(line[-2])
        solvation_energy = float(line[-1])
        name = line[1]

        charges[name] = charge
        references[name] = solvation_energy

    f.close()

    # load in all XYZ files
    unique_atoms = []
    molecules = []
    molecules_coordinates = []
    molecules_atoms = []
    molecules_charges = []
    molecules_references = []
    f = open(exam_subset, 'r')
    for line in f:
        line = line.replace("\n", "")

        atoms, coordinates = gms.get_coordinates_xyz(__XYZDIR__ + line + '.xyz')

        # mask. sort out if masked is used
        if args.mask:
            if not set(np.unique(atoms)).issubset(args.mask):
                continue

        molecules.append(line)
        molecules_coordinates.append(coordinates)
        molecules_atoms.append(atoms)
        molecules_charges.append(charges[line])
        molecules_references.append(references[line])

        atoms = np.unique(atoms)
        atoms = [str(x) for x in atoms] # fixes weird type error
        unique_atoms += atoms

    f.close()


    molecules = np.array(molecules)
    molecules_charges = np.array(molecules_charges)
    molecules_references = np.array(molecules_references)

    if args.mask:
        print
        print "using subset left from mask", args.mask
        print molecules
        print

    # find unique atoms
    unique_atoms = np.unique(unique_atoms)
    unique_parameters = [gms.get_atom(atom) for atom in unique_atoms]
    unique_parameters = np.array(unique_parameters)
    unique_parameters.sort()
    unique_parameters -= 1 # convert from atom to index
    n_parameters = len(unique_parameters)

    print "Optimizing atom types:", unique_parameters+1

    # start parameters
    parameters = np.zeros((n_atom_types))

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

    else:
        parameters += 2.30 # Klamt default parameter


    # optimize parameters
    parameters_optimize = parameters[unique_parameters]
    print "Starting from:", parameters_optimize



    def energy_rmsd(parameter_view, verbose=True):

        if verbose:
            start_time = time()

        # local = deepcopy(global) # 107
        # local[view] = parameter_optimize # 7
        # local med optimize # 107

        parameter_view = [round(x, 2) for x in parameter_view]

        parameters_local = np.zeros((n_atom_types))
        parameters_local[unique_parameters] = parameter_view

        # energies = gms.get_energies(molecules_atoms,
        #                    molecules_coordinates,
        #                    molecules_charges,
        #                    header,
        #                    parameters_local, workers=30)

        energies = gms.get_energies_nodes(molecules_atoms,
                        molecules_coordinates,
                        molecules_charges,
                        header,
                        parameters_local,
                        workers=workers,
                        chunksize=chunksize,
                        node_list=nodelist,
                        hostname=slave_master)

        find_nan = np.argwhere(np.isnan(energies))

        energies = np.array(energies)

        if len(find_nan) != 0 and verbose:
            print find_nan
            for nanidx in find_nan:
                print "nan:", molecules[nanidx]
                energies[nanidx] = molecules_references[nanidx]*1.5

        ermsd = rmsd(molecules_references, energies)

        if verbose:
            endtime = time() - start_time
            if ermsd > 20.0:
                print "ENERGIES:", list(energies)
                print "REFERENC:", list(molecules_references)

            print "PARAM:", list(parameter_view), "RMSD:", ermsd, "TIME:", endtime

        return ermsd


    def energy_jacobian(parameters, dx=0.001):

        gradient = np.zeros((n_parameters))

        for i, p in enumerate(parameters):

            cparameters = deepcopy(parameters)

            dp = cparameters[i] * dx

            cparameters[i] += dp
            rmsd_high = energy_rmsd(cparameters)

            cparameters[i] -= (2.0 * dp)
            rmsd_low = energy_rmsd(cparameters)

            dr = (rmsd_high - rmsd_low) / (2.0 * dp)
            gradient[i] = dr


        norm = linalg.norm(gradient)
        print "PARAM:", list(parameters), "NORM:", norm

        return gradient

    def energy_jacobian_nodes(parameters, dx=0.001):


        return


    # print energy_rmsd(parameters_optimize)
    # print energy_jacobian(parameters_optimize)
    #
    # All radii should be positive
    bounds = [(0.5, None) for _ in range(n_parameters)]

    # Methods
    # method_name = "L-BFGS-B"
    method_name = "Powell"
    method_name = "Nelder-Mead"

    function_name = energy_rmsd
    # jacobian_function = energy_jacobian

    parameters_initial = deepcopy(parameters_optimize)

    # For nelder-mead
    print opt.minimize(function_name, parameters_initial,
                 method=method_name,
                 options={"maxiter": 300, "disp": True})

    # for powell

    # for LBFGS
    # print opt.minimize(function_name, parameters_initial,
    #              # jac=jacobian_function,
    #              method=method_name,
    #              bounds=bounds,
    #              options={"maxiter": 300, "disp": True, "eps": 0.01})

    # print unique_parameters+1

if __name__ == "__main__":
    main()

