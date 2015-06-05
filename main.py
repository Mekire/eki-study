#!/usr/bin/env python
from __future__ import print_function
import readRecordings
import interface


def main():
    """
    Generates and displays data in a similar format to the iPython notebook.
    """
    readRecordings.readSMR()
    readRecordings.readABF()
    csvFiles = {"EKI":"MN1-Ib_EKI.csv", "WT":"MN1-Ib_WT.csv"}
    data = interface.load_data(csvFiles)
    bMeasures = interface.calc_measure_data(data)
    quants, quartiles = interface.calc_quants_quarts(bMeasures, "burstDur")
    min_quant, max_quant = interface.get_min_max_quants(quants)

    interface.print_data(data)
    interface.print_specific(bMeasures, "burstDur")
    interface.print_all_specific(bMeasures)
    interface.print_quants_quartiles(quartiles, min_quant, max_quant)
    

if __name__ == "__main__":
    main()


