"""
Implements an API for interacting with data.
Functionality mirrors that which was demonstrated in the iPython notebook.
"""


from __future__ import print_function
import scipy as sc
import burstAnalysis as burst


KINDS = ("EKI", "WT")

MEASURES = ("burstDur", "cycleDur", "dutyCycle", "qI")

MEASURES_FULL = {"burstDur"  : "Burst Duration",
                 "cycleDur"  : "Cycle Duration",
                 "dutyCycle" : "Duty Cycle",
                 "qI"        : "Quiescence Interval"}


def print_data(data):
    """
    Print all data from data dictionary.
    """
    for kind in KINDS:
        print("\n{} data:".format(kind))
        for data_set in data[kind]:
            print(data_set, end="\n\n")


def print_specific(measures_data, specific):
    """
    Print the data of a specific measure for both EKI and WT.
    """
    for kind in KINDS:
        print("\n{} {}:".format(kind, MEASURES_FULL[specific]))
        for measure in measures_data[kind][specific]:
            print(measure, end="\n\n")


def print_all_specific(measures_data):
    """
    Print the data for all measures for both EKI and WT.
    """
    for measure in MEASURES:
        print_specific(measures_data, measure)


def load_data(csv_files):
    """
    Loads csv data into a dictionary of the form:
        {"EKI" : eki_data, "WT" : wt_data}

    csv_files should be a dictionary of the form:
        {"EKI" : "eki_filename.csv", "WT" : "wt_filename.csv"}
    """
    return {k : burst.getCSVData(v) for k,v in csv_files.items()}


def calc_measure_data(data):
    """
    Given a data dictionary created by load_data, create a dictionary of
    all measures for both EKI and WT.
    """
    return {k : burst.calcBurstMeasures(v) for k,v in data.items()}


def calc_all_CDF(measures_data, binEdges):
    """
    Calculate CDFs for all measures and place in a dictionary.

    This function can be used individually, but if the intent is to calculate
    quants and quartiles, use calc_quants_quarts as the function is called
    there for you. 
    """
    cdfs = {'EKI': {}, 'WT': {}}
    for kind in KINDS:
        for measure in MEASURES:
            samples = measures_data[kind][measure]
            cdfs[kind][measure] = burst.calcCDFs(sampleList=samples,
                                                 binEdges=binEdges)
    return cdfs


def get_mean_CDF(cdfs, specific):
    """
    Get the means of the CDFs for one specific measure.
    Called in calc_quants_quarts.
    """
    return {kind : cdfs[kind][specific].mean(0) for kind in KINDS}


def calc_quants_quarts(measures_data, specific, binSize=1e-3, maxPt=30.0):
    """
    Given a dictionary of all measure data and a specific measure to
    investigate, finds the quants and quarts of that measure.

    Return value is a tuple of the quant and quart dictionaries.
    """
    binEdges = sc.arange(0, maxPt, binSize)
    cdfs = calc_all_CDF(measures_data, binEdges)
    means = get_mean_CDF(cdfs, specific)
    quants = {k : sc.interpolate.interp1d(means[k],binEdges) for k in KINDS}
    quartiles = {k : quants[k](sc.arange(0.25,1.00,0.25)) for k in KINDS}
    return quants, quartiles


def get_min_max_quants(quants):
    """
    Finds both the min and max quants and returns them as a tuple.
    """
    min_quants = {k : float(quants[k](0.01)) for k in KINDS}
    max_quants = {k : float(quants[k](1.0)) for k in KINDS}
    return min_quants, max_quants


def print_quants_quartiles(quartiles, min_quants, max_quants):
    """
    Print function for the quants and quartiles.
    """
    fields = "min & q1 & q2 & q3 & max \\"
    template = "{0:.2f} & {2:.2f} & {3:.2f} & {4:.2f} & {1:.2f} \\"
    for k in KINDS:
        filled = template.format(min_quants[k], max_quants[k], *quartiles[k])
        print("\n".join([k, fields, filled]), end="\n\n")
