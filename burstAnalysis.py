# Author: Erin C. McKiernan
# Last modified: 1 June 2015


# import packages
from __future__ import division
import csv
import scipy as sc
import scipy.stats
import scipy.interpolate
import numpy as np
import pylab as gr


DEBUG = False

KINDS = ("EKI", "WT")

MEASURES = ("burstDur", "cycleDur", "dutyCycle", "qI")

MEASURES_FULL = {"burstDur"  : "Burst Duration",
                 "cycleDur"  : "Cycle Duration",
                 "dutyCycle" : "Duty Cycle",
                 "qI"        : "Quiescence Interval"}

PLOT_MEASURE_INFO = {"burstDur"  : {"xlabel"  : "Burst Duration (secs)",
                                    "ylabel"  : "Relative Frequency",
                                    "maxTime" : 30.0,
                                    "bin"     : 1.0},
                     "cycleDur"  : {"xlabel"  : "Cycle Duration (secs)",
                                    "ylabel"  : "",
                                    "maxTime" : 30.0,
                                    "bin"     : 1.0},
                     "dutyCycle" : {"xlabel"  : "Duty Cycle",
                                    "ylabel"  : "Relative Frequency",
                                    "maxTime" : 1,
                                    "bin"     : 0.04},
                     "qI"        : {"xlabel"  : "Quiescence Interval (secs)",
                                    "ylabel"  : "",
                                    "maxTime" : 30.0,
                                    "bin"     : 1.0}}

 
def getCSVData(CSVfile, nHeaderLines=1):
    """
    Extract data from the csv files

    If there are more header rows, change the 1 to the number of header rows.
    """
    data = []
    with open(CSVfile,'rb') as csvf:
        csvdata = csv.reader(csvf, delimiter=',', quotechar='|')  
        for nrows,row in enumerate(csvdata):
            if nrows >= nHeaderLines:
                floatStrings = " ".join(row[1:]).split() 
                data.append(map(float, floatStrings))
    return data


def testMonotonicIncrease(aRRay):
    """
    Test if the averages of CDFs are monotonically increasing.
    Note: CDFs should always be increasing - this is used for debugging.

    Untested since change. -Sean
    """
    x = not len(gr.find(gr.diff(aRRay)<0)) > 0
    description = 'non-decreasing' if x else 'NOT monotonic'
    print('The array is {}.'.format(description))
    return x
    

def calcBurstMeasures(data):
    """
    Calculate the burst measures from one data set. The input argument
    is a list whose elements are realizations of a bursting process
    containing the beginning and end of each burst as they occurred.

    Burst duration - time elapsed from start to end of one burst
    Cycle duration - time elapsed from start of one burst to start of next 
    Duty cycle - burst duration divided by cycle duration 
    Quiescence interval - time elapsed from end of burst to start of next

    Example:
    dataEKI = getCSVData(csvFiles['EKI'])
    bMeasures = calcBurstMeasures(dataEKI)
    """
    measures = {measure : [] for measure in MEASURES}
    for datum in data:
        burstStart = sc.array(datum[0::2])
        burstStop = sc.array(datum[1::2])
        burstDur = burstStop-burstStart
        cycleDur = sc.diff(burstStart)
        measures["burstDur"].append(burstDur)
        measures["cycleDur"].append(cycleDur)
        measures["dutyCycle"].append(burstDur[:-1]/cycleDur)
        measures["qI"].append(burstStart[1:]-burstStop[:-1])
        if DEBUG:
            for attr in MEASURES:
                print("{}:\n{}".format(MEASURES_FULL[attr],measures[attr][-1]))
            print("\n")
    return measures


def calcBins():
    """
    Set intervals and get bins for each measure.
    """
    range_args = {'burstDur' : (1.0, 50.0, 1.0), 'cycleDur' : (1.0, 50.0, 1.0),
                  'dutyCycle' : (0.0, 1.0, 0.04), 'qI' : (1.0, 50.0, 1.0)}
    return {measure : sc.arange(*range_args[measure]) for measure in MEASURES}


def calcCDFs(sampleList, binEdges):
    """
    Calculate cumulative distribution functions (CDFs).
    
    Example:
    dataEKI = getCSVData(csvFiles['EKI'])
    bMeasures = calcBurstMeasures(dataEKI)
    cdfs = calcCDFs(sampleList=bMeasures['burstDur'],binEdges=bins['burstDur'])
    """
    cdfs = []
    for i,sample in enumerate(sampleList):
        nPts = sc.float64(len(sample))
        cdfs.append(sc.zeros(len(binEdges)))
        for j,point in enumerate(binEdges):
            cdfs[i][j] = len(sc.where(sample<=point)[0])/nPts
    return sc.array(cdfs)
            

def graphCDFsOneSample(cdfs, binEdges, ax, grpColor='r',
                       indivMark=':', strLabel='avg'):
    """
    Plot cumulative distribution functions (CDFs).

    Example usage:
    fig = gr.figure(num=figNum,figsize=(9,10))
    gr.ioff();
    ax1 = fig.add_subplot(111)
    cdfs = calcCDFs(sampleList=bMeasures['burstDur'],
                    binEdges=bins['burstDur'])
    graphCDFsOneSample(cdfs, binEdges=bE, ax=ax1, grpColor='r',
                       indivMark=':',strLabel='avg')
    gr.ion()
    gr.draw()
    """
    cdfAvg = cdfs.mean(0)
    for sample in cdfs:
        ax.plot(binEdges, sample, grpColor+indivMark, lw=1)        
    ax.plot(binEdges, cdfAvg, grpColor, alpha=0.6, lw=3, label=strLabel)
    ax.legend(loc='lower right')    
    return ax


def graphMultiCDF(bins, bMeasures):
    """
    Multipanel figure showing different comparisons between cdfs. 
    Assumes you get 4 lists of pairs of cdfs, 1 for each measure.
    """
    colors = {'EKI' : 'r', 'WT' : 'b'}
    cdfs = {'EKI': {}, 'WT': {}}
    for kind in cdfs:
        for attr in MEASURES:
            cdfs[kind][attr] = calcCDFs(sampleList=bMeasures[kind][attr],
                                        binEdges=bins[attr])
    f4 = gr.figure(figsize=(15,15))
    gr.ioff()
    rows, cols = 2, 2
    for i,attr in enumerate(MEASURES, start=1):
        ax = f4.add_subplot(rows, cols, i)
        for kind in cdfs:
            myCDF = cdfs[kind][attr]
            myBins = bins[attr]
            graphCDFsOneSample(cdfs=myCDF, binEdges=myBins, ax=ax,
                               grpColor=colors[kind], strLabel=kind)
            ax.set_title(MEASURES_FULL[attr])
            if attr != 'dutyCycle':
                ax.set_xlabel('secs')
                ax.set_xlim(0, 30)
    gr.ion()
    gr.draw()


def calcAvgRelFreq(Xbins, Xmeasures):
    """
    Calculate average relative frequencies.
    """
    rf = sc.zeros((len(Xmeasures),len(Xbins)-1))
    for n,measure in enumerate(Xmeasures):
        rf[n,:] = sc.histogram(measure, bins=Xbins)[0]/len(measure)
    return {"relFreqs" : rf, "avgRelFreq" : sc.mean(rf, 0)}


def getRelFreqs(bins, csvFiles): 
    """
    Set intervals and get bins to calculate relative frequencies. 
    """
    burstRFs = {"EKI" : {}, "WT" : {}}
    for k,v in csvFiles.items():
        print(v)
        burstMeasure = calcBurstMeasures(getCSVData(v,1))
        for attr in MEASURES:
            burstRFs[k][attr] = calcAvgRelFreq(bins[attr], burstMeasure[attr])
    return burstRFs


def pdfComparison(csvFiles, label1='WT', label2='EKI', maxRF=0.33, alpha1=1,
                  alpha2=1, skinfactor=0.3):
    """
    Plots histograms for group comparisons.
    Specifics for individual plot appearance are contained in the dictionary
    PLOT_MEASURE_INFO found at the top of this module.
    """
    bins = calcBins()
    burstRFs = getRelFreqs(bins, csvFiles)
    rows, cols = 2, 2
    ax = []
    f0 = gr.figure(figsize=(13,13))
    gr.ioff()
    for i,attr in enumerate(MEASURES, start=1):
        plot = f0.add_subplot(rows, cols, i)
        measure = PLOT_MEASURE_INFO[attr]
        binw = measure["bin"]
        shift = binw*(1-skinfactor)/2.0
        plot.bar(bins[attr][:-1],
                 burstRFs[label1][attr]['avgRelFreq'],
                 width=binw, color='w', alpha=alpha1, label=label1)
        plot.bar(bins[attr][:-1]+shift,
                 burstRFs[label2][attr]['avgRelFreq'],
                 width=binw*skinfactor, color='k', alpha=alpha2, label=label2)
        plot.legend(shadow=True, fancybox=True, ncol=1, prop={'size':16})
        plot.set_xlim(0, measure["maxTime"])
        plot.set_ylim(0, maxRF)
        plot.set_xlabel(measure["xlabel"])
        plot.set_ylabel(measure["ylabel"])
    f0.text(0.05,0.95,'A',size=20)
    f0.text(0.5,0.95,'B',size=20)
    f0.text(0.05,0.5,'C',size=20)
    f0.text(0.5,0.5,'D',size=20)
    f0.subplots_adjust(left=0.12, bottom=0.1, right=0.95,
                       top=0.95, wspace=0.2, hspace=0.2)
    gr.ion()
    gr.draw()
