#!/usr/bin/env python

# from multiprocessing.managers import BaseManager
from multiprocessing.managers import SyncManager

import subprocess
import os

import Queue
import time
import socket

import numpy as np

from shell import shell

# TODO Get hostnames
# TODO Get how many cpu's per hostnames

__here__ = os.path.abspath(os.path.dirname(__file__))

def make_server_manager(port, authkey, verbose=False):

    job_queue = Queue.Queue()
    result_queue = Queue.Queue()

    class JobQueueManager(SyncManager):
        pass

    JobQueueManager.register('get_job_queue', callable=lambda: job_queue)
    JobQueueManager.register('get_result_queue', callable=lambda: result_queue)

    manager = JobQueueManager(address=('', port), authkey=authkey)

    manager.start()

    if verbose:
        print 'Server started at port %s' % port

    return manager


def make_slave(slave, host, port, authkey,
               directory="/scratch",
               workers=1,
               debug=False,
               verbose=False,
               gamess_tar=''):

    """

    return: slave, a process

    """

    client = __here__
    client += "/client.py"
    client += " --ip " + host
    client += " --port " + str(port)
    client += " --auth " + authkey
    client += " --dir " + directory

    client += " --workers " + str(workers)

    if gamess_tar != '':
        client += " --gamess " + gamess_tar

    if debug:
        client += " --test"

    pipe = True
    if pipe:
        client += " > /home/charnley/dev/2017-pcm-parameters/node_info/"+slave+".log " + " 2> /home/charnley/dev/2017-pcm-parameters/node_info/"+slave+".log"

    cmd = []
    cmd.append("ssh")
    cmd.append(slave)
    cmd.append('"' + client + '"')

    cmd = " ".join(cmd)

    if debug:
        print cmd

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, shell=True)

    return proc


def make_slaves(slaves, hostname, port, authkey,
                directory="/scratch",
                workers=1,
                debug=False,
                gamess_tar=''):

    procs = []

    for node in slaves:
        procs.append(make_slave(node,
                                hostname,
                                port,
                                authkey,
                                directory=directory,
                                workers=workers,
                                debug=debug,
                                gamess_tar=gamess_tar))

    return procs


def terminate(slaves):
    """ Terminate slave process if they are still on
    """

    [proc.terminate() for proc in slaves]

    return


def get_slurm_information():

    # Is this run in a SLURM Enviroment?
    slurm_id = os.environ["SLURM_JOB_ID"]
    slurm_nodes = os.environ["SLURM_JOB_NODELIST"]
    slurm_nodes = shell('scontrol show hostname', shell=True)
    slurm_workers = os.environ["SLURM_JOB_CPUS_PER_NODE"]

    # clean up nodes
    slurm_nodes = slurm_nodes.split("\n")
    slurm_nodes = [node.strip() for node in slurm_nodes]
    slurm_nodes = list(filter(None, slurm_nodes))

    # Just need the ncpus
    slurm_workers = slurm_workers.split('(')[0]
    slurm_workers = int(slurm_workers)

    return slurm_nodes, slurm_workers, slurm_id


def check_for_slurm():
    return os.environ["SLURM_JOB_ID"] != ""


def hostname():
    return socket.gethostname()


if __name__ == "__main__":

    """ Test the network communication is working
    """

    PORTNUM = 5000
    AUTHKEY = "haxorboy"
    hostname = socket.gethostname()
    node_list = ["node634", "node678"]

    node_list, workers, jobid = get_slurm_information()

    print "removing host from node_list", hostname, node_list
    node_list.remove(hostname)
    print "Testing on", node_list

    manager = make_server_manager(PORTNUM, AUTHKEY, verbose=True)
    shared_jobs = manager.get_job_queue()
    shared_results = manager.get_result_queue()

    N = 100000
    nums = np.arange(N)

    chunksize = 100

    for i in range(0, len(nums), chunksize):
        shared_jobs.put(nums[i:i + chunksize])

    # Start worker nodes
    slaves = make_slaves(node_list, hostname, PORTNUM, AUTHKEY, debug=True, directory="/scratch/"+jobid)

    # Make everything sync up
    time.sleep(5)

    # What did the slaves say?
    for slave in slaves:
        out = slave.stdout.readline()
        # out, err = slave.communicate()
        print out

    # Wait until all results are ready in shared_result_q
    numresults = 0
    resultdict = {}
    while numresults < N:
        print "Waiting for results"
        outdict = shared_results.get()
        resultdict.update(outdict)
        numresults += len(outdict)
        print "Total results now", numresults

    # Kill slaves
    print "Killing slaves"
    [proc.terminate() for proc in slaves]

    # print resultdict
    print "Total results gatherd:", numresults

    time.sleep(2)
    manager.shutdown()

    print "Done testing"

