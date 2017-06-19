#!/usr/bin/env python

import sys
import numpy as np
import scipy.optimize as opt
from copy import deepcopy
import numpy.linalg as linalg

import gamess as gms

n_atom_types = 107


__XYZDIR__ = "/home/charnley/dev/2017-mnsol-data/jobs/xyz/" # local
__XYZDIR__ = "/home/charnley/dev/2017-mnsol/jobs/xyz/" # sunray


def rmsd(V, W):
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += (v - w)**2.0
    return np.sqrt(rmsd/N)


def main():

    args = sys.argv[1:]

    expe_values = args[0] # list of solvation free energies
    exam_subset = args[1] # list of moleculenames


    # TODO load in all charges
    # TODO load in experimental values (dict)
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

    # TODO load in all XYZ files

    unique_atoms = []

    molecules = []
    molecules_coordinates = []
    molecules_atoms = []
    molecules_charges = []
    molecules_references = []
    f = open(exam_subset, 'r')
    for line in f:
        line = line.replace("\n", "")
        molecules.append(line)

        atoms, coordinates = gms.get_coordinates_xyz(__XYZDIR__ + line + '.xyz')

        molecules.append(line)
        molecules_coordinates.append(coordinates)
        molecules_atoms.append(atoms)
        molecules_charges.append(charges[line])
        molecules_references.append(references[line])

        unique_atoms += np.unique(atoms)

    f.close()

    # find unique atoms
    unique_atoms = np.unique(unique_atoms)
    unique_parameters = [gms.get_atom(atom) for atom in unique_atoms]
    unique_parameters = np.array(unique_parameters)
    unique_parameters.sort()
    unique_parameters -= 1 # convert from atom to index
    n_parameters = len(unique_parameters)

    # start parameters
    parameters = np.zeros((n_atom_types))
    parameters += 2.223 # Klamt default parameter

    # optimize parameters
    parameters_optimize = parameters[unique_parameters]

    # define gamess header
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

    ! diis=.t.
 $end

 $pcm
    solvnt=WATER
 $end


 $tescav
    ! mthall=1
    ! ntsall=60
 $end

"""


    def energy_rmsd(parameter_view):

        # local = deepcopy(global) # 107
        # local[view] = parameter_optimize # 7
        # local med optimize # 107

        parameters_local = np.zeros((n_atom_types))
        parameters_local[unique_parameters] = parameter_view

        energies = gms.get_energies(molecules_atoms,
                           molecules_coordinates,
                           molecules_charges,
                           header,
                           parameters_local, workers=30)

        energies = np.array(energies)
        ermsd = rmsd(molecules_references, energies)

        print "PARAM:", list(parameter_view), "RMSD:", ermsd

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


    # print energy_rmsd(parameters_optimize)
    # print energy_jacobian(parameters_optimize)

    # All radii should be positive
    bounds = [(0.5, None) for _ in range(n_parameters)]

    # Methods
    method_name = "L-BFGS-B"
    # method_name = "Nelder-Mead"

    # method_name = "L-BFGS-B"
    function_name = energy_rmsd
    jacobian_function = energy_jacobian

    parameters_initial = deepcopy(parameters_optimize)

    # print opt.minimize(function_name, parameters_initial,
    #              method=method_name,
    #              bounds=bounds,
    #              options={"maxiter": 300, "disp": True, "eps": 0.001})

    print opt.minimize(function_name, parameters_initial,
                 jac=jacobian_function,
                 method=method_name,
                 bounds=bounds,
                 options={"maxiter": 300, "disp": True, "eps": 0.001})

    print unique_parameters+1

if __name__ == "__main__":
    main()

