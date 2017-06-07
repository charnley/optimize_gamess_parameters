#!/usr/bin/env python

import sys
import numpy as np
import scipy.optimize as opt

import gamess as gms

n_atom_types = 107


__XYZDIR__ = "/home/charnley/dev/2017-mnsol-data/jobs/xyz/"


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
        molecules_coordinates.append(coordinates)
        molecules_atoms.append(atoms)
        molecules_charges.append(charges[line])
        molecules_references.append(references[line])

    f.close()

    # start parameters
    parameters = np.zeros((n_atom_types))
    parameters += 2.223 # Klamt default parameter

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


    print gms.get_energies(molecules_atoms,
                           molecules_coordinates,
                           molecules_charges,
                           header,
                           parameters, workers=1)

    # All radii should be positive
    bounds = [(0.0, None) for _ in range(n_atom_types)]

    # Methods
    # method_name = "L-BFGS-B"
    # method_name = "Nelder-Mead"

    method_name = "L-BFGS-B"

    # opt.minimize(function_name, parameters,
    #              jac=jacobian_function,
    #              method=method_name,
    #              bounds=bounds)

    # minimize(cost.cost, parameters, jac=cost.jacobian, method="L-BFGS-B",
    #         options={"maxiter": 1000, "disp": True}, bounds=bounds)


if __name__ == "__main__":
    main()

