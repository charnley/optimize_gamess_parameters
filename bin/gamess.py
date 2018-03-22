#!/usr/bin/env python

from subprocess import Popen, PIPE

import multiprocessing as mp
import numpy.linalg as linalg
import numpy as np
import re
import os
from time import time

from numpy.linalg import norm

import socket
import nodes

__GAMESS__ = "/home/charnley/opt/gamess-mndod/rungms"
__GAMESS__ = "/opt/gamess/mndod/rungms"
__GAMESS__ = "/home/charnley/opt/gamess/mndod-fast/rungms"
# __GAMESS__ = "/scratch/666/mndod-fast/rungms"

__SCRATCH__ = "scr"
__GMSSCR__ = "/home/charnley/scr"

__XYZDIR__ = "/home/charnley/dev/2017-mnsol/jobs/xyz/" # sunray
__XYZDIR__ = "/home/charnley/dev/2017-mnsol/jobs/xyz-gamess-hf/" # hf gas
__XYZDIR__ = "/home/charnley/dev/2017-mnsol/jobs/xyz-gamess-pm6/" # hf gas

# TODO Change XYZDIR for header type!


global __DFTBHUB__
__DFTBHUB__ = {
        "br": -0.0573,
        "c": -0.1492,
        "ca": -0.0340,
        "cl": -0.0697,
        "f": -0.1623,
        "h": -0.1857,
        "i": -0.0433,
        "k": -0.0339,
        "mg": -0.02,
        "n": -0.1535,
        "na": -0.0454,
        "o": -0.1575,
        "p": -0.14,
        "s": -0.11,
        "zn": -0.03 }



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


def generate_dftb_header(atoms):

    header = ""

    # unique, but in the right order
    indexes = np.unique(atoms, return_index=True)[1]
    atoms = [atoms[index] for index in sorted(indexes)]

    # Generate hubder
    line = ""
    # line = " $dftb hubder(1)="
    # hubs = []
    # for atom in atoms:
    #     hubs.append(str(__DFTBHUB__[atom.lower()]))
        # line = line+str(__DFTBHUB__[atom])+","

    # line += ",".join(hubs) + " $end\n"

    __dftbdir__ = "/home/andersx/dftb/3ob-3-1/"
    dftbsk = ""
    dftbsk += " $dftbsk\n"

    for atom1 in atoms:
        for atom2 in atoms:
            dftbsk += atom1 + " " + atom2 + " \"" + __dftbdir__ + atom1 + "-" + atom2 + '.skf\"\n'

    dftbsk += " $end\n"

    header += "\n" + line + "\n" + dftbsk

    return header


def generate_rin(atom_list, radius_list):

    header = ""
    fmt = "   rin({}) = {}"

    if type(radius_list) == type(dict()):

        for i, atom in enumerate(atom_list):
            idx = get_atom(atom)
            header += fmt.format(i+1, radius_list[idx]) + "\n"

    else: # type list
        for i, atom in enumerate(atom_list):
            idx = get_atom(atom) - 1
            header += fmt.format(i+1, radius_list[idx]) + "\n"

    return header


def genereate_data(atom_list, coordinates):

    header = ""

    return header


def make_input_file(atoms, coordinates, header, charge, parameters):

    out = ""
    coordinate_line = "{:2s}     {:3d}.0      {:13.10f}    {:13.10f}   {:13.10f}\n"

    # fix charge
    header = header.replace("icharg=0", "icharg="+str(charge))

    if "dftb" in header:
        header += generate_dftb_header(atoms)

    # energy model
    out += header +"\n"

    if parameters is not None:
        # radii parameters
        out += " $pcmcav\n"
        out += " \n\n"
        out += generate_rin(atoms, parameters)
        out += " $end\n\n"

    # Coordinates
    out += " $data\n"
    out += "title\nC1\n"
    for atom, coord in zip(atoms, coordinates):
        out += coordinate_line.format(atom, get_atom(atom), *coord)
    out += " $end\n"

    return out


def get_energy(atom_list, coordinates, header, charge, parameters, filename="test", smd=True):

    # Generate GAMESS header
    filecontent = make_input_file(atom_list, coordinates, header, charge, parameters)

    # write the input file
    f = open(filename + ".inp", 'w')
    f.write(filecontent)
    f.close()

    # TODO Check for uncov PCM

    # run the calculation
    if smd:
        cmd = __GAMESS__ + " " + filename + ".inp" + r' | grep "FREE ENERGY OF SOLVATION\|FINAL R" '

    else:
        cmd = __GAMESS__ + " " + filename + ".inp" + r' | grep -A1 "LEC+CAV+DIS+REP\|FINAL R" '

    out = shell(cmd , shell=True)
    numbers = re.findall(r'[-]?\d+\.\d*(?:[Ee][-\+]\d+)?', out)

    # get the energy
    if not len(numbers) > 0:
        energy = float("nan")
        hostname = socket.gethostname()
        shell('cp '+filename+'.inp /home/charnley/dev/2017-pcm-parameters/fails/'+hostname+'_'+filename+'.inp', shell=True)

    else:
        if smd:
            gas = float(numbers[0])
            pcm = float(numbers[1])

            # Test if solvent scf converged
            if pcm == 0.0:
                return float("nan")

            energy = float(numbers[-2])
        else:
            energy = float(numbers[-1])

    if abs(energy) > 500.0:
        hostname = socket.gethostname()
        shell('cp '+filename+'.inp /home/charnley/dev/2017-pcm-parameters/fails_conv/'+hostname+'_'+filename+'.inp', shell=True)


    # clean scratch
    # TODO should be more general
    shell("rm "+filename+".*", shell=True)

    # clear gms tmpfiles
    # os.remove(filename+".inp")

    return energy


def split_jobs(jobs, workers):
    """ Split job array into MP """
    return [jobs[i::workers] for i in xrange(workers)]


def energy_worker(i, ids, matoms, mcoordinates, mcharges, parameters, header, energies):

    for ii, atoms, coordinates, charge in zip(ids, matoms, mcoordinates, mcharges):
        start_time = time()
        energies[ii] = get_energy(atoms, coordinates, header, charge, parameters, filename=str(ii))
        print ii, energies[ii], time()-start_time

    return


def get_energies(molecules_atoms, molecules_coordinates, molecules_charges, header, parameters, workers=1):

    n_molecules = len(molecules_atoms)

    energies = mp.Array("d", [0.0 for _ in xrange(n_molecules)])

    molecules_id = split_jobs(range(n_molecules), workers)
    molecules_atoms = split_jobs(molecules_atoms, workers)
    molecules_coordinates = split_jobs(molecules_coordinates, workers)
    molecules_charges = split_jobs(molecules_charges, workers)

    pwd = os.getcwd()

    os.chdir(__SCRATCH__)

    processes = [mp.Process(target=energy_worker,
        args=(i, molecules_id[i],
                 molecules_atoms[i],
                 molecules_coordinates[i],
                 molecules_charges[i],
                 parameters, header, energies)) for i in xrange(workers)]

    for p in processes: p.start()
    for p in processes: p.join()

    os.chdir(pwd)

    return energies


def get_energies_nodes(molecules_atoms,
                       molecules_coordinates,
                       molecules_charges,
                       header,
                       parameters,
                       workers=1,
                       chunksize = 50,
                       node_list=["localhost"],
                       hostname="sunray"):


    PORTNUM = 5000
    AUTHKEY = "anotherkeykey"

    manager = nodes.make_server_manager(PORTNUM, AUTHKEY)

    n_molecules = len(molecules_atoms)
    energies = np.zeros(n_molecules)
    n_nodes = len(node_list)

    if chunksize*n_nodes > n_molecules:
        chunksize = int(n_molecules/n_nodes)

    shared_jobs = manager.get_job_queue()
    shared_results = manager.get_result_queue()

    for i in range(0, n_molecules, chunksize):

        job_molecules_id = list(range(n_molecules))[i:i + chunksize]
        job_molecules_atoms = molecules_atoms[i:i + chunksize]
        job_molecules_coordinates = molecules_coordinates[i:i + chunksize]
        job_molecules_charges = molecules_charges[i:i + chunksize]

        shared_jobs.put([job_molecules_id,
                         job_molecules_atoms,
                         job_molecules_coordinates,
                         job_molecules_charges,
                         header,
                         parameters])


    # Start worker nodes
    #TODO Change gamess to more general
    slaves = nodes.make_slaves(node_list, hostname, PORTNUM, AUTHKEY, workers=workers,
            gamess_tar="/home/charnley/opt/gamess/github.tar.gz")

    # for slave in slaves:
        # out = slave.stdout.readline()
        # out, err = slave.communicate()
        # print out
        # print err

    numresults = 0
    resultdict = {}

    while numresults < n_molecules:
        outdict = shared_results.get()
        resultdict.update(outdict)
        numresults += len(outdict)

    nodes.terminate(slaves)

    for i in resultdict:
        energies[i] = resultdict[i]

    return energies



def get_energies_parameters_nodes(molecules_atoms, molecules_coordinates, molecules_charges,
                header, parameters_list,
                parameter_index_list=None,
                workers=1,
                chunksize=10,
                node_list=[""]):


    n_molecules = len(molecules_atoms)
    n_parameters = len(parameters_list)

    chunksize = min(chunksize, n_molecules)

    print "making zeroes"

    energies_list_list = np.zeros((n_molecules, n_parameters))
    energies_list_list = energies_list_list.flatten()
    n_jobs = len(energies_list_list)

    print "making server"

    PORTNUM = 5000
    AUTHKEY = "haxorboy"
    hostname = nodes.hostname()

    manager = nodes.make_server_manager(PORTNUM, AUTHKEY)
    shared_jobs = manager.get_job_queue()
    shared_results = manager.get_result_queue()

    jobs_idx = list(range(n_jobs))

    print "n_jobs", n_jobs

    print "par * mol", n_parameters, n_molecules, n_parameters * n_molecules

    for j in range(n_parameters):

        parameters = parameters_list[j]

        for i in range(0, n_molecules, chunksize):

            # TODO raisecondiition
            # should be
            # to = min(frm+chunksize, (j+1)*n_molcules)
            frm = j*n_molecules + i
            to = frm + chunksize

            job_idx_list = jobs_idx[frm:to]

            job_molecules_atoms = molecules_atoms[i:i + chunksize]
            job_molecules_coordinates = molecules_coordinates[i:i + chunksize]
            job_molecules_charges = molecules_charges[i:i + chunksize]

            shared_jobs.put([job_idx_list,
                            job_molecules_atoms,
                            job_molecules_coordinates,
                            job_molecules_charges,
                            header,
                            parameters])


    print "Generating slaves"

    # Start worker nodes
    slaves = nodes.make_slaves(node_list, hostname, PORTNUM, AUTHKEY, workers=workers)

    numresults = 0
    resultdict = {}

    print "Waiting for results"

    while numresults < n_jobs:
        outdict = shared_results.get()
        resultdict.update(outdict)
        numresults += len(outdict)

    nodes.terminate(slaves)

    for i in resultdict:
        energies_list_list[i] = resultdict[i]

    energies_list_list = energies_list_list.reshape((n_parameters, n_molecules))

    return energies_list_list


def main():

    import argparse
    import sys

    description = ""

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options]',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-p', '--parameters', action='store', nargs='+', default='', help='')
    parser.add_argument('-e', '--header', action='store', default='', help='')

    parser.add_argument('-x', '--xyz_file', action='store', default='', help='')
    parser.add_argument('-c', '--charge', type=int, action='store', default=0, help='')

    parser.add_argument('-l', '--xyz_list', action='store', default='', help='')

    parser.add_argument('-s', '--scan', action='store_true', default=False, help='')
    parser.add_argument('-m', '--exam', action='store_true', default=False, help='')

    parser.add_argument('-w', '--double_workers', action='store_true', default=False, help='Double the amount of workers (usefull for SQM methods)')

    args = parser.parse_args()


    # Node setup
    try:
        nodelist, workers, jobid = nodes.get_slurm_information()

        if args.double_workers:
            workers *= 2

    except KeyError:
        nodelist=["node634", "node678", "node637", "node662"]
        workers=8
        job=""

    # for sqm
    # workers = 16
    chunksize = 90

    # Read parameters for solvent radii
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


    # define gamess header
    if args.header:
        header_type = args.header.split("/")
        header_type = header_type[-1]
        __XYZDIR__ = "/home/charnley/dev/2017-pcm-parameters/xyz/"+header_type+"/"

        with open(args.header, 'r') as myfile:
                header = myfile.read()

    else:
        header = """ \

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
    mxts=15000

    smd=.t.
 $end

 $tescav
    mthall=4
    ntsall=60
 $end


"""


    # generate input file for molecule using specific parameters
    if args.xyz_file:
        atoms, coordinates = get_coordinates_xyz(args.xyz_file)

        # test dict
        # parameters = {}
        # parameters[1] = 2.0
        # parameters[6] = 5.0
        # parameters[7] = 0.9

        print make_input_file(atoms, coordinates, header, args.charge, parameters)
        quit()


    if args.xyz_list:

        charge_file = "/home/charnley/dev/2017-pcm-parameters/lists/mnsol_water.charge"
        moldb = {}
        f = open(charge_file, 'r')
        for line in f:
            line = line.split(',')
            charge = int(line[1])
            name = line[0]

            moldb[name] = {}
            moldb[name]['charge'] = charge
        f.close()

        # Read the molecules from molecule list
        molecules = [] # names
        molecules_coordinates = [] # xyz
        molecules_atoms = [] # atom type
        molecules_charges = [] # charges

        f = open(args.xyz_list, 'r')
        for line in f:

            line = line.replace("\n", "")
            line = line.split()

            out = get_coordinates_xyz(__XYZDIR__ + line[0] + '.xyz')

            if not out:
                continue

            atoms, coordinates = out

            molecules.append(line[0])
            molecules_coordinates.append(coordinates)
            molecules_atoms.append(atoms)
            molecules_charges.append(moldb[line[0]]['charge'])

        f.close()


    if args.xyz_list and args.exam:

        for molecule, atoms, charges, coordinates in zip(molecules, molecules_atoms, molecules_charges, molecules_coordinates):

            print molecule, charge

            f = open(molecule+'.inp', 'w')
            f.write(make_input_file(atoms, coordinates, header, charges, parameters))
            f.close()

        quit()


    if args.xyz_list and args.scan:

        # Create parameter list
        parameters_list = []
        parameter_range = np.arange(0.5, 3.1, 0.1)
        parameters_scan_list = [1, 6]
        parameters = np.zeros(107)

        print "Generating parameter scan for", parameters_scan_list

        for pi in parameter_range:
            for pj in parameter_range:
                par = {}
                par[1] = pi
                par[6] = pj
                parameters_list.append(par)

        # number of molecules in list
        n_molecules = len(molecules)

        print "submitting", len(parameters)

        energies = get_energies_parameters_nodes(molecules_atoms, molecules_coordinates, molecules_charges,
                header, parameters_list,
                parameter_index_list=parameters_scan_list,
                workers=workers,
                chunksize=chunksize,
                node_list=nodelist)

        print "molecules", molecules

        for i, par in enumerate(parameters_list):
            print par[1], par[6], "#", list(energies[i])

        quit()




    if args.xyz_list and not args.scan:

        # number of molecules in list
        n_molecules = len(molecules)

        # workers per node
        workers = 16

        # split job list in chunks of
        chunksize = 30

        start_time = time()


        energies = get_energies_nodes(molecules_atoms,
                        molecules_coordinates,
                        molecules_charges,
                        header,
                        parameters,
                        workers=workers,
                        chunksize=chunksize,
                        node_list=nodelist,
                        hostname=nodes.hostname())

        # energies = get_energies_nodes(molecules_atoms,
        #                 molecules_coordinates,
        #                 molecules_charges,
        #                 header,
        #                 parameters,
        #                 workers=workers,
        #                 chunksize=chunksize,
        #                 # node_list=["node634"])
        #                 # node_list=["node634", "node678"])

        endtime = time() - start_time

        for i in xrange(n_molecules):
            print molecules[i], energies[i]

        print "time:", endtime


if __name__ == "__main__":


    main()






    # legacy code

    quit()

    import sys

    args = sys.argv[1:]

    if len(args) == 0:
        print "test usage:"
        print "gamess.py xyz_list parameter_file"
        quit()

    molecules_file = args[0]
    parameter_file = args[1]

    # Read parameters for solvent radii
    parameters = []
    f = open(parameter_file)
    for line in f:
        line = float(line)
        parameters.append(line)
    f.close()

    # Read the molecules from molecule list
    molecules = [] # names
    molecules_coordinates = [] # xyz
    molecules_atoms = [] # atom type
    molecules_charges = [] # charges

    f = open(molecules_file, 'r')
    for line in f:

        line = line.replace("\n", "")
        line = line.split()

        atoms, coordinates = get_coordinates_xyz(__XYZDIR__ + line[0] + '.xyz')

        molecules.append(line[0])
        molecules_coordinates.append(coordinates)
        molecules_atoms.append(atoms)
        molecules_charges.append(line[1])

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

    ! Calculate the cavitation energy
    ICAV=1

    ! Calculation of dispersion and repulsion free energy
    IDISP=1

 $end

"""


    # charge = 0
    # os.chdir(__)
    # print get_energy(molecules_atoms, molecules_coordinates, header, molecules_charges, parameters)

    # energies = get_energies(molecules_atoms, molecules_coordinates, molecules_charges, header, parameters, workers=1)


    for chunksize in range(10, 210, 10):

        start_time = time()

        energies = get_energies_nodes(molecules_atoms,
                        molecules_coordinates,
                        molecules_charges,
                        header,
                        parameters,
                        workers=16,
                        chunksize=chunksize,
                        # node_list=["node634"])
                        # node_list=["node634", "node678"])
                        node_list=["node634", "node678", "node637", "node662"])

        endtime = time() - start_time

        print chunksize, endtime


