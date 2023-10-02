#!/usr/bin/env python
'''
    File name: OBS-HailDataPreprocessor.py
    Author: Andreas Prein
    E-mail: prein@ucar.edu
    Date created: 27.03.2017
    Date last modified: 27.03.2017

    ##############################################################
    Purpos:

    Contains all statistic functions used in the syntethic hail model

'''

from pylab import *
def fnConditional_PDFs(rgrObs,         # gridded observed hail
                       rgrLonO,                   #  longitude grid of obs
                       rgrLatO,                    #  latitude grid of obs
                       rgrVars,                    # hail model imput variables
                       rgrLonV,                    # hail model imput variables longitudes
                       rgrLatV,                     # hail model imput variables latitudes
                       rgsVars,                     # list of variable names
                       iBox= None):            # grid cells added/substracted for spatial averaging

    # This function calculates PDFs of variables conditioned on the
    # location and time of observed hail

    import datetime
    import glob
    from netCDF4 import Dataset
    import sys, traceback
    import dateutil.parser as dparser
    import string
    # from ipdb import set_trace as stop
    import numpy as np
    import os
    import subprocess
    import pandas as pd
    from scipy import stats
    import copy
    import pylab as plt
    import shlex
    import bisect
    import csv
    # from shapely.geometry import Point
    # from shapely.geometry.polygon import Polygon
    from matplotlib import path
    import os.path
    import matplotlib.gridspec as gridspec
    
    from matplotlib.patches import Polygon
    # from mpl_toolkits.basemap import Basemap, cm
    import string
    from datetime import timedelta
    from scipy import interpolate
    import scipy

    # grid cells added/substracted for spatial averaging
    if iBox == None:
        iBox=2
        
    rgsVars=list(rgsVars)
    # for these variables the maximum value in the search box is calculated
    rgsMaxVars=['CAPEmax',
#                 'CIN',
#                 'SRH01',
                'SRH03',
#                 'VS01',
                'VS03',
#                 'VS06',
#                 'VS012'
               ]

    # loop over all observed hail reports
    rgiHailObs=np.array(np.where(rgrObs == 1))
    rgrConditVars=np.zeros((rgrVars.shape[1], rgiHailObs.shape[1])); rgrConditVars[:]=np.nan
    for ih in range(rgiHailObs.shape[1]):
        iTT=rgiHailObs[0][ih]
        rLatO=rgrLatO[rgiHailObs[1][ih],rgiHailObs[2][ih]]
        rLonO=rgrLonO[rgiHailObs[1][ih],rgiHailObs[2][ih]]
        rgrVloc=np.where((rgrLonV == rLonO) & (rgrLatV == rLatO))
        try:
            iLaAct=rgrVloc[0][0]
        except:
            stop()
        iLoAct=rgrVloc[1][0]
        # loop over variables
        try:
            iTTmax=np.where(np.max(rgrVars[iTT:iTT+2,rgsVars.index('CAPEmax'),iLaAct-iBox:iLaAct+iBox+1, iLoAct-iBox:iLoAct+iBox+1]) == rgrVars[iTT:iTT+2,rgsVars.index('CAPEmax'),iLaAct-iBox:iLaAct+iBox+1, iLoAct-iBox:iLoAct+iBox+1])[0][0]
        except:
            continue
        iTT=iTTmax+iTT
        for va in range(len(rgsVars)):
            if rgsVars[va] in rgsMaxVars:
                rgrConditVars[rgsVars.index(rgsVars[va]),ih]=np.nanmax(rgrVars[iTT,rgsVars.index(rgsVars[va]),iLaAct-iBox:iLaAct+iBox+1, iLoAct-iBox:iLoAct+iBox+1])
            else:
                rgrConditVars[rgsVars.index(rgsVars[va]),ih]=np.nanmean(rgrVars[iTT,rgsVars.index(rgsVars[va]),iLaAct-iBox:iLaAct+iBox+1, iLoAct-iBox:iLoAct+iBox+1])


    if np.sum(np.isnan(rgrConditVars)) > 0:
        rgiNaN=np.where(np.isnan(rgrConditVars[0,:]) == 1)
        rgrConditVars=np.delete(rgrConditVars,rgiNaN,1)

    # calculate the PDF of the samples
    rgrPercentilesHailVar=np.zeros((rgrVars.shape[1],101)); rgrPercentilesHailVar[:]=np.nan
    for va in range(rgrVars.shape[1]):
        rgrPercentilesHailVar[va,:]=scipy.stats.scoreatpercentile(rgrConditVars[va,:], range(101))

    return rgrConditVars, rgrPercentilesHailVar





# __________________________________________________________________
# __________________________________________________________________
# function that transferes values to probabilities

def transfere(x,rgrTH):
    import numpy as np
    # x ... variable that should be transfered to probabilities
    # rgrTH ... threshold for this variable

    # this function sets the probability to zero for values below the 1st
    # threshold, to 0 to 1 between the 1st and 2nd threshold, to 1 between
    # the 2nd and 3rd threshold and to 0 above the 3rd threshold

    rMeanTH1=np.mean([rgrTH[1],rgrTH[0]]); rTHtrans1=(rgrTH[1]-rgrTH[0])*0.3
    rMeanTH2=np.mean([rgrTH[3],rgrTH[2]]); rTHtrans2=(rgrTH[3]-rgrTH[2])*0.3
    rTHMedian=np.mean([rgrTH[1],rgrTH[2]])
    y=np.piecewise(x, [(x < rTHMedian),\
                           (x >= rTHMedian)],\
                           [lambda x: 0.5+0.5*np.tanh((x-rMeanTH1)/rTHtrans1),\
                                lambda x: 1-(0.5+0.5*np.tanh((x-rMeanTH2)/rTHtrans2))])
    # plt.plot(x.flatten(),y.flatten(),'ro')
    # plt.show()
    # stop()
    return y
