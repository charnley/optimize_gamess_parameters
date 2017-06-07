#!/usr/bin/env python

from subprocess import Popen, PIPE

import multiprocessing as mp
import numpy.linalg as linalg
import numpy as np
import re
import os

from numpy.linalg import norm

__GAMESS__ = "/home/charnley/opt/gamess-mndod/rungms"
__GAMESS__ = "/home/charnley/dev/2016-mndod/gamess/git-mndod/rungms"
__SCARTCH__ = "scr"
__GMSSCR__ = "/home/charnley/scr"


global __ATOM_LIST
__ATOM_LIST__ = [ x.strip() for x in ['h ','he', \
      'li','be','b ','c ','n ','o ','f ','ne', \
      'na','mg','al','si','p ','s ','cl','ar', \
      'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', \
      'zn','ga','ge','as','se','br','kr', \
      'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', \
      'cd','in','sn','sb','te','i ','xe', \
      'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', \
      'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', \
      'au','hg','tl','pb','bi','po','at','rn', \
      'fr','ra','ac','th','pa','u ','np','pu'] ]


def shell(cmd, shell=False):

    if shell:
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    else:
        cmd = cmd.split()
        p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)

    output, err = p.communicate()
    return output


def get_atom(atom):
    global __ATOM_LIST__
    atom = atom.lower()
    return __ATOM_LIST__.index(atom) + 1


def delete_scr(filename):
    filename = filename.split(".")[0]
    filename += "*"
    os.chdir(__GMSSCR__)
    os.remove(filename)


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
        exit("Could not obtain the number of atoms in the .xyz file.")

    # Skip the title line
    f.readline()

    # Use the number of atoms to not read beyond the end of a file
    for lines_read, line in enumerate(f):

        if lines_read == n_atoms:
            break

        atom = re.findall(r'[a-zA-Z]+', line)[0]
        atom = atom.upper()

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


def generate_rin(atom_list, radius_list):

    header = ""
    fmt = "   rin({}) = {}"

    for i, atom in enumerate(atom_list):
        idx = get_atom(atom) - 1
        header += fmt.format(i+1, radius_list[idx]) + "\n"

    return header


def genereate_data(atom_list, coordinates):

    header = ""

    return header


def get_energy(atom_list, coordinates, header, charge, parameters, filename="test"):

    coordinate_line = "{:2s}     {:3d}.0      {:13.10f}    {:13.10f}   {:13.10f}\n"

    # write the input file
    f = open(filename + ".inp", 'w')

    # fix charge
    header = header.replace("icharg=0", "icharg="+str(charge))

    # energy model
    f.write(header + "\n")

    # radii parameters
    f.write(" $pcmcav\n")
    f.write(" \n\n")
    f.write(generate_rin(atom_list, parameters))
    f.write(" $end\n\n")

    # Coordinates
    f.write(" $data\n")
    f.write("title\nC1\n")
    for atom, coord in zip(atom_list, coordinates):
        f.write(coordinate_line.format(atom, get_atom(atom), *coord))
    f.write(" $end\n")
    f.close()

    # run the calculation
    cmd = __GAMESS__ + " " + filename + ".inp" + " | grep -A1 \"LEC+CAV+DIS+REP\" "
    out = shell(cmd , shell=True)
    numbers = re.findall(r'[-]?\d+\.\d*(?:[Ee][-\+]\d+)?', out)

    # clean scratch
    # TODO should be more general
    shell("rm ~/scr/"+filename+"*", shell=True)

    if not len(numbers) > 0:
        return float("nan")

    # clear gms tmpfiles
    os.remove(filename+".inp")

    # get the energy
    energy = float(numbers[-1])

    return energy


def split_jobs(jobs, workers):
    """ Split job array into MP """
    return [jobs[i::workers] for i in xrange(workers)]


def energy_worker(i, ids, matoms, mcoordinates, mcharges, parameters, header, energies):

    for ii, atoms, coordinates, charge in zip(ids, matoms, mcoordinates, mcharges):
        energies[ii] = get_energy(atoms, coordinates, header, charge, parameters, filename=str(ii))
        print energies[ii]

    return


def get_energies(molecules_atoms, molecules_coordinates, molecules_charges, header, parameters, workers=1):

    n_molecules = len(molecules_atoms)

    energies = mp.Array("d", [0.0 for _ in xrange(n_molecules)])

    molecules_id = split_jobs(range(n_molecules), workers)
    molecules_atoms = split_jobs(molecules_atoms, workers)
    molecules_coordinates = split_jobs(molecules_coordinates, workers)
    molecules_charges = split_jobs(molecules_charges, workers)

    os.chdir(__SCARTCH__)

    processes = [mp.Process(target=energy_worker,
        args=(i, molecules_id[i],
                 molecules_atoms[i],
                 molecules_coordinates[i],
                 molecules_charges[i],
                 parameters, header, energies)) for i in xrange(workers)]

    for p in processes: p.start()
    for p in processes: p.join()

    print list(energies)

    return energies


def get_jacobian(parameters, dx=0.001, workers=1):

    grad = []

    for i, p in enumerate(parameters):
        print p

    norm = linalg.norm(gradient)

    return gradient


def numerical_jacobian(F, r, h=0.001):
    """
    Create Jacobian for function F, using numerical differentiation
    """
    n = len(r)
    J = [[0.0]*n for i in range(n)]
    for i in range(0, n):
        r_a = copy.deepcopy(r)
        r_b = copy.deepcopy(r)
        r_a[i] = r_a[i] + h
        r_b[i] = r_b[i] - h
        f_a = F(r_a)
        f_b = F(r_b)
        for j in range(0, n):
            J[j][i] = (f_a[j] - f_b[j])/(2.0*h)
    return J



if __name__ == "__main__":

    import sys

    args = sys.argv[1:]

    xyz_file = args[0]
    parameter_file = args[1]

    atom_list, coordinates = get_coordinates_xyz(xyz_file)

    parameters = []
    f = open(parameter_file)
    for line in f:
        line = float(line)
        parameters.append(line)
    f.close()

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
    charge = 0
    os.chdir(__SCARTCH__)

    print get_energy(atom_list, coordinates, header, charge, parameters)

