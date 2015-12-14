import os
import xml.etree.ElementTree as ET

import numpy as np


def parseTimingData(timingfile):
    """Parses the timing data of one rank.

    Returns a dictionary that maps timer names to measurements."""
    with open(timingfile, "r") as file:
        pairs = [line.split(" ") for line in file.readlines()]

        return {timer: float(time) for (timer, time) in pairs}


class Job:
    @staticmethod
    def valid(path):
        """A job is valid if it has timing data."""
        timingdir = os.path.join(path, "timing")

        return os.path.isdir(timingdir) and len(os.listdir(timingdir)) > 0

    @staticmethod
    def parse(path):
        scenario = ET.parse(os.path.join(path, "scenario.xml"))
        parallel = scenario.find("parallel")
        nx, ny, nz = [int(parallel.attrib["numProcessors" + d])
                      for d in ["X", "Y", "Z"]]
        rank = nx * ny * nz
        timingpath = os.path.join(path, "timing")
        timingdata = [
            parseTimingData(os.path.join(timingpath, "rank-{}".format(r)))
            for r in range(rank)
        ]

        return Job(rank, nx, ny, nz, path, scenario, timingdata)

    def __init__(self, rank, nx, ny, nz, path, scenario, timingdata):
        self.rank = rank
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.path = path
        self.scenario = scenario
        self.timingdata = timingdata


def parseExperiment(directory):
    """Parses the output of a run of the `scale` script as a list of jobs."""
    jobs = [os.path.join(directory, j) for j in os.listdir(directory)]
    data = [Job.parse(j) for j in jobs if Job.valid(j)]

    return data


def parseExperiments(directories):
    """Parse and merge the results of multiple runs of scaling experiment.

    Returns a tuple of measurement names and two nested dicts that map ranks and
    processor configurations to evaluations, where an evaluation is a tuple of
    sample size, array of means and array of variances. The two arrays are in
    the same order as the keys. So if "total" is at index 4 in the list of keys,
    the mean total is at index 4 of the means array.

    """
    experiments = [parseExperiment(d) for d in directories]
    data = {}

    for e in experiments:
        for job in e:
            if not job.rank in data:
                data[job.rank] = {}

            nconf = (job.nx, job.ny, job.nz)
            if not nconf in data[job.rank]:
                data[job.rank][nconf] = []

            data[job.rank][nconf] += job.timingdata

    firstRank = list(data.values())[0]
    firstMeasurements = list(firstRank.values())[0]

    # Assume that all measurements have the same keys
    keys = sorted(list(firstMeasurements[0].keys()))
    npdata = {}

    for rank, confs in data.items():
        npdata[rank] = {}

        for conf, measurements in confs.items():
            n = len(measurements)
            m = len(keys)
            npmeasurements = np.empty((m, n))

            for i in range(n):
                for j in range(m):
                    npmeasurements[j, i] = measurements[i][keys[j]]

            npdata[rank][conf] = (n, np.mean(npmeasurements,
                                             axis=1), np.var(npmeasurements,
                                                             axis=1))

    return (keys, npdata)
