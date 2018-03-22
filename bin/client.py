#!/usr/bin/env python

import argparse
import sys
import os

import socket
import time
import Queue as queue

import numpy as np

import multiprocessing

from multiprocessing.managers import SyncManager

import gamess
from shell import shell

def make_client_manager(ip, port, authkey):

    class ServerQueueManager(SyncManager):
        pass

    ServerQueueManager.register('get_job_queue')
    ServerQueueManager.register('get_result_queue')

    manager = ServerQueueManager(address=(ip, port), authkey=authkey)
    manager.connect()

    print 'Client connected to %s:%s' % (ip, port)

    return manager


def factorize_naive(n):
    """ A naive factorization method. Take integer 'n', return list of
        factors.
    """
    if n < 2:
        return []
    factors = []
    p = 2

    while True:
        if n == 1:
            return factors

        r = n % p
        if r == 0:
            factors.append(p)
            n = n / p
        elif p * p >= n:
            factors.append(n)
            return factors
        elif p > 2:
            # Advance in steps of 2 over odd numbers
            p += 2
        else:
            # If p == 2, get to 3
            p += 1
    assert False, "unreachable"


def test_client(job_queue, result_queue, workers=1):

    print "Starting client"

    while True:
        try:
            job = job_queue.get_nowait()
            outdict = {n: factorize_naive(n) for n in job}
            result_queue.put(outdict)

        except queue.Empty:
            break

    print "Stopping client"

    return


def start_client(job_queue, result_queue, stay_awake=False, workers=1):

    while True:
        try:

            print "getting job"
            job = job_queue.get_nowait()

            job_molecules_id,\
            job_molecules_atoms,\
            job_molecules_coordinates,\
            job_molecules_charges, \
            header, \
            parameters = job

            N = len(job_molecules_id)

            energies = np.zeros(N)

            energies = gamess.get_energies(job_molecules_atoms,
                         job_molecules_coordinates,
                         job_molecules_charges,
                         header,
                         parameters,
                         workers=workers)

            outdict = {key: energy for key, energy in zip(job_molecules_id, energies)}

            # time.sleep(1)

            result_queue.put(outdict)

        except queue.Empty:

            print "Empty queue"

            return

    return


def setup_client(directory, gamessball="/home/charnley/opt/gamess/mndod-fast.tar.gz"):
    """
    Setup client
    """

    tarfile = gamessball.split("/")[-1]
    gmsexe = tarfile.split('.')[0] + "/rungms "
    rungms = directory + "/" + gmsexe

    if not directory == "" or not directory == ".":
        os.chdir(directory)

    # check if gamess is already copied
    if not os.path.exists(tarfile):
        print "GAMESS does not exists"
        print shell('cp '+gamessball+' .')
        print shell('tar -xvf '+tarfile)

    # check if gamess is correc tar
    else:
        print "gamess exists"
        md5_here = shell("md5sum " + tarfile).split()[0]
        md5_there = shell("md5sum " + gamessball).split()[0]

        if md5_here != md5_there:
            print "wrong version"
            print shell('cp '+gamessball+' .')
            print shell('tar -xvf '+tarfile)


    gamess.__GAMESS__ = rungms

    return

def clean_client():


    return

def main():

    description = """"""

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options]',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-w', '--workers', action='store', type=int, default=1, metavar='N', help='Number of CPU\'s to use on this slave')
    parser.add_argument('-i', '--ip', action='store', default="127.0.0.1", help='IP address of the slave manager. Must be hostname or IP. (default: 127.0.0.1)')
    parser.add_argument('-p', '--port', action='store', default=5000, type=int, help='Port address of the slave manager.')
    parser.add_argument('-a', '--auth', action='store', default='', help='Slave manager authentication key.')

    parser.add_argument('-s', '--stay_awake', action='store_true', help='Stay awake when queue is empty (default=False)')
    parser.add_argument('-d', '--dir', action='store', metavar='DIR', default='', help='Working directory (default=Empty)')

    parser.add_argument('-t', '--test', action='store_true', help='Test the client setup (from nodes.py)')

    parser.add_argument('-g', '--gamess', action='store', default='/home/charnley/opt/gamess/mndod-fast.tar.gz', help='tarball for gamess executable')

    args = parser.parse_args()

    hostname = socket.gethostname()
    workers = args.workers
    IP = args.ip
    PORTNUM = args.port
    AUTHKEY = args.auth

    manager = make_client_manager(IP, PORTNUM, AUTHKEY)
    job_queue = manager.get_job_queue()
    result_queue = manager.get_result_queue()

    print "client here"

    if args.test:
        print "Testing client"
        test_client(job_queue, result_queue, workers=workers)
        return


    # change directory and setup gamess locally
    setup_client(args.dir, gamessball=args.gamess)
    print shell('pwd')

    gamess.__SCRATCH__ = args.dir
    gamess.__SCARTCH__ = "scr"

    start_client(job_queue, result_queue, stay_awake=args.stay_awake, workers=workers)

    return


if __name__ == "__main__":

    main()

