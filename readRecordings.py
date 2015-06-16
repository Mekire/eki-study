# Author: Erin C. McKiernan
# Code adapted from Neo IO documentation https://pythonhosted.org/neo/io.html 

from __future__ import print_function
import os
from neo import io
from matplotlib import pyplot as plt


DEBUG = False


ABF_DEFAULTS = {"xmin" : 160, "xmax" : 280,
                "ymin" :-160, "ymax" :  20,
                "figsize" : None}

SMR_DEFAULTS = {"xmin" : 163, "xmax" : 276,
                "ymin" : -70, "ymax" :  0,
                "figsize" : (18,5)}


def _display(bl, **kwargs):
    """
    Display function for both readABF and readSMR.
    """
    if DEBUG:
        print(bl.segments)
        print(bl.segments[0].analogsignals)
        print(bl.segments[0].eventarrays)
    for seg in bl.segments:
        siglist = seg.analogsignals
        timepoints = siglist[0].times
        plt.figure(figsize=kwargs["figsize"])
        plt.plot(timepoints, siglist[0])
        plt.plot(timepoints, siglist[1])
        plt.xlim(kwargs["xmin"], kwargs["xmax"])
        plt.ylim(kwargs["ymin"], kwargs["ymax"])
        plt.xlabel("Time (seconds)")
        plt.ylabel("Voltage (mV)")
        plt.show()


def readABF(recname="ABF-files/09721000.abf", dirname=".", **range_kwargs):
    """
    Reading and plotting recordings in ABF (Axon) format.
    See ABF_DEFAULTS for legal keyword arguments.
    """
    kwargs = ABF_DEFAULTS.copy()
    kwargs.update(range_kwargs)
    filename = os.path.join(dirname, recname)
    r = io.AxonIO(filename)
    bl = r.read_block(lazy=False, cascade=True)
    _display(bl, **kwargs)


def readSMR(recname="SMR-files/09706000.SMR", dirname=".", **range_kwargs):
    """
    Reading and plotting recordings in SMR (Spike2) format
    See SMR_DEFAULTS for legal keyword arguments.
    """
    kwargs = SMR_DEFAULTS.copy()
    kwargs.update(range_kwargs)
    filename = os.path.join(dirname, recname)
    r = io.Spike2IO(filename)
    bl = r.read(lazy=False, cascade=True)[0]
    _display(bl, **kwargs)
