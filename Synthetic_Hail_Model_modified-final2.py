#!/usr/bin/env python
# coding: utf-8

'''
    File name: SyntheticHailModel.py
    Author: Boris Blanc
    E-mail: boris.blanc@env.ethz.ch
    Date created: 24.08.2017
    Date last modified: 30.07.2023

    ##############################################################
    Purpos:

    Evaluates the conditions under which hail forms.

    1) hail observations are imported

    2) ERA-Interim conditions at the day and location where hail occured
    are read

    3) conditions for non-hail days are calculated

    4) investigation if hail days are significantly different from non hail
    days in several U.S. subregions

    5) trends in days with hail conditions are calculated

'''

from dateutil import rrule
import datetime
from datetime import timedelta
import glob
from netCDF4 import Dataset
import sys, traceback
import dateutil.parser as dparser
import string
import numpy as np
import numpy.ma as ma
import os
# from mpl_toolkits import basemap
# import ESMF
import pickle
import subprocess
import pandas as pd
from scipy import stats
import copy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import pylab as plt
import random
import scipy.ndimage as ndimage
import matplotlib.gridspec as gridspec
# from mpl_toolkits.basemap import Basemap, cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from pylab import *
import string
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import shapefile
# import shapely.geometry
import shapefile
import math
from scipy.stats.kde import gaussian_kde
from math import radians, cos, sin, asin, sqrt
# from shapely.geometry import Polygon, Point
from scipy.interpolate import interp1d
import csv
import os.path
import matplotlib.gridspec as gridspec
import scipy
import matplotlib.path as mplPath
import calendar 
from calendar import monthrange
from calendar import isleap
from pdb import set_trace as stop

from numpy import linspace, meshgrid
from scipy.interpolate import griddata

import os
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

from cartopy import config
import cartopy.crs as ccrs
import geopandas as gpd
import cartopy.feature as cfeature

from matplotlib.cm import ScalarMappable

import random
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as mcolors
import matplotlib.patches as mpatch
from matplotlib.ticker import LogLocator
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.colors import ListedColormap
from netCDF4 import Dataset
import h5py


def grid(x, y, z, resX=20, resY=20):
    "Convert 3 column data to matplotlib grid"
    xi = linspace(min(x), max(x), resX)
    yi = linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi, interp='linear')
    X, Y = meshgrid(xi, yi)
    return X, Y, Z


########################################
#                            Settings

iRadius = int(sys.argv[1]) # grid cells that are blended out around hail observations (4)
qT=int(sys.argv[2]) #4
qD=int(sys.argv[3]) #26
YYYY=int(sys.argv[4]) #26

dStartDay=datetime.datetime(YYYY, 1, 1,0,0)
dStopDay=datetime.datetime(YYYY, 12, 31,23,59)
#generate time vectors
rgdTimeDD = pd.date_range(dStartDay, end=dStopDay, freq='d')
rgdTime1H = pd.date_range(dStartDay, end=dStopDay, freq='1h')

rgdFullTime=pd.date_range(datetime.datetime(1959, 1, 1, 0),
                          end=datetime.datetime(2022, 12, 31, 23), freq='1h')
rgdFullTime_day=pd.date_range(datetime.datetime(1959, 1, 1, 0),
                          end=datetime.datetime(2022, 12, 31, 23), freq='d')

iMonths=np.unique(rgdTimeDD.month)
iYears=np.unique(rgdTimeDD.year)

rgsLableABC=string.ascii_lowercase+string.ascii_uppercase

sSaveDataDir='/glade/scratch/bblanc/ERA5_hail_model/Hail_Probabilities_final/'

sERAconstantFields='/glade/scratch/bblanc/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc'

rgsModelVars=['CAPEmax',
              'FLH',
              'VS03',
              'SRH03']

rgsVarCombinations=['CAPEmax-FLH']
rgsFinVars=['CAPEmax-FLH','CAPEmax','SRH03','VS03']

rgsVarsAct=(np.array([rgsModelVars[ll].split('-') for ll in range(len(rgsModelVars))]))
rgsERAVariables = np.array([item for sublist in rgsVarsAct for item in sublist])


rgsVarMaxLim=['FLH']
rgsVarMinLim=['CAPEmax','SRH03','VS03','FLH']

rMinHailSize= 2.5 # this is the maximum diameter of a hail stone in cm

# # provide radio sounding data if desired
# sSoundDataFile='/glade/p_old/work/prein/observations/RadioSoundings_Uni-Wyoming/data/'
# rgsSoundStat=['72469']

# iCorrelPred=20 # 20 # how many 2D predictors should be used
# iPosComb=scipy.math.factorial(len(rgsERAVariables))/(scipy.math.factorial(2) * scipy.math.factorial(len(rgsERAVariables)-2)) # max combinations of variables


iPercentile=np.array([qT,qD+qT])  # percentiles that are excluded from the distribution
Qtail=iPercentile[0]
QDelta=iPercentile[1]-Qtail
rgrPercSteps=np.linspace(Qtail,QDelta+Qtail,5)
rgrProbSteps=np.append(rgrPercSteps[0]-(rgrPercSteps[1]-rgrPercSteps[0]),rgrPercSteps)
rgrProb=0.5+0.5*np.tanh((rgrProbSteps-np.mean([rgrProbSteps[0],rgrProbSteps[-1]]))/((rgrProbSteps[-1]-rgrProbSteps[0])*0.3))[1:]
rgrProb[-1]=1


setup_string = 'search-rad-'+str(iRadius)+'_qT-'+str(qT)+'_qD-'+str(qD)


#-----------------------------------------------------------------------
# Read the observations from NOAA
#-----------------------------------------------------------------------

sSaveFolder="/glade/work/bblanc/HailObs/SPC_data/ncdf_files/"
sDataSet = 'NOAA'


#RawData = None
hailsize = ['Hail', 'HailSize']
rgrNOAAObs=np.zeros((2, len(rgdTimeDD), 130, 300))

for ii in range(len(hailsize)):
    dd=0
    for yy in range(len(iYears)):
        #print('Loadind NOAA data for year '+str(iYears[yy]))
        sFileName=sSaveFolder+'SPC-Hail-StormReports_gridded-75km_'+str(iYears[yy])+'.nc'

        # read in the variables
        ncid=Dataset(sFileName, mode='r')

        rgrLatGrid=np.squeeze(ncid.variables['lat'][:,0])
        rgrLonGrid=np.squeeze(ncid.variables['lon'][0,:])
        rgrLonNOAA=np.asarray(([rgrLonGrid,]*rgrLatGrid.shape[0]))
        rgrLatNOAA=np.asarray(([rgrLatGrid,]*rgrLonGrid.shape[0])).transpose()
        yearlength=365 + calendar.isleap(iYears[yy])
        
        rgrNOAAObs[ii,dd:dd+yearlength,:,:]=np.squeeze(ncid.variables[hailsize[ii]][:])
        
        dd=dd+yearlength
        ncid.close()

#Get monthly and yearly hail observations
daily_NOAAobs = np.zeros((len(iYears), 12, rgrLonNOAA.shape[0], rgrLonNOAA.shape[1]))
for yy in range(len(iYears)):
    for mm in range(12):
        daily_NOAAobs[yy,mm,:,:] = np.sum(rgrNOAAObs[0,(rgdTimeDD.year == iYears[yy]) & (rgdTimeDD.month == (mm+1)),:,:], axis=0)
        
monthly_NOAAobs=np.mean(daily_NOAAobs, axis=0)
yearly_NOAAobs=np.sum(daily_NOAAobs, axis=1)

#------------------------------------------------#
# Save Latitude, Longitude and Height in arrays
#------------------------------------------------#

rgsERAdata='/glade/scratch/bblanc/ERA5_hail_model/ERA5-hailpredictors/'
sERAconstantFields='/glade/scratch/bblanc/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc'
ncid=Dataset(sERAconstantFields, mode='r')
rgrLat=np.squeeze(ncid.variables['latitude'][:])
rgrLon=np.squeeze(ncid.variables['longitude'][:])
rgrHeight=(np.squeeze(ncid.variables['Z'][:]))/9.81
ncid.close()

rgiSize=rgrHeight.shape

# all lon grids have to be -180 to 180 degree
rgrLonERA=np.asarray(([rgrLon,]*rgrLat.shape[0]))
rgrLatERA=np.asarray(([rgrLat,]*rgrLon.shape[0])).transpose()
rgi180=np.where(rgrLonERA > 180)
rgrLonERA[rgi180]=rgrLonERA[rgi180]-360.

# check if we only want to read a subregion
rgiDomain=None
if rgiDomain != None:
    iloW=rgiDomain[3]
    iloE=rgiDomain[1]
    ilaN=rgiDomain[0]
    ilaS=rgiDomain[2]
else:
    iloW=rgrLonERA.shape[1]
    iloE=0
    ilaN=0
    ilaS=rgrLonERA.shape[0]
    

# ___________________________________________________________________________
#   Load the hail predictors
# ___________________________________________________________________________


rgrERAVarall=np.zeros((len(rgdTimeDD), len(rgsERAVariables), rgrLatERA.shape[0],rgrLatERA.shape[1]))
iyear=0
#loop over years 
for yy in range(len(iYears)):
    iday=0
    outfile = '/glade/scratch/bblanc/ERA5_hail_model/ERA5_yearly_hailpreedictors/ERA5_dailymax_' + str(iYears[yy])+ '.npz'
    yearlength=365 + isleap(iYears[yy])
    
    if os.path.isfile(outfile) == False:
        rgrERAVarsyy=np.zeros((yearlength, len(rgsERAVariables), rgrLatERA.shape[0],rgrLatERA.shape[1]))

        #Loop over months
        for mm in range(len(iMonths)):
            # __________________________________________________________________
            # start reading the ERA data
            dDate=rgdTimeDD[0]
            monthlength=int(monthrange(iYears[yy],mm+1)[1])

            rgrERAdataAll=np.zeros((monthlength*24, len(rgsERAVariables), rgrLatERA.shape[0],rgrLatERA.shape[1])); rgrERAdataAll[:]=np.nan

            print (str(iYears[yy])+'    Loading ERA-5 from year '+str(iYears) +', days:'+ str(iday)+':'+str(monthlength+iday))

            # loop over variables
            sFileName=rgsERAdata +str(iYears[yy])+str("%02d" % (mm+1))+'_ERA-5_HailPredictors_newSRH03.nc'
            ncid=Dataset(sFileName, mode='r')

            for va in range(len(rgsERAVariables)):
                rgrDataTMP=np.squeeze(ncid.variables[rgsERAVariables[va]][:])

                try:
                    rgrERAdataAll[:,va,:,:]=rgrDataTMP
                except:
                    rgiSizeTMP=rgrDataTMP.shape
                    if (rgiSizeTMP[1]-rgrERAdataAll.shape[2] < 2) & (rgiSizeTMP[2]-rgrERAdataAll.shape[3] < 2):
                        rgrERAdataAll[:,va,:rgiSizeTMP[1],:rgiSizeTMP[2]]=rgrDataTMP
                    else:
                        print ('    DATASETS DO NOT HAVE THE SAME DIMENSIONS!')
                        stop()
            ncid.close()    

            #reshaping the monthly data into daily data to be able to extract the daily max CAPE
            dayrgrERAdataAll=rgrERAdataAll.reshape(int(monthlength),24,rgrERAdataAll.shape[1],rgrERAdataAll.shape[2],rgrERAdataAll.shape[3])
            timemax=np.argmax(dayrgrERAdataAll[:,:,0,:,:],axis=1)
            timemax2=np.zeros((dayrgrERAdataAll.shape[0],dayrgrERAdataAll.shape[1],dayrgrERAdataAll.shape[2],dayrgrERAdataAll.shape[3],dayrgrERAdataAll.shape[4]))
            timemax2[:,:,:,:,:]=timemax[:,None,None,:]
            timemax2=timemax2.astype(int)

            #Take the predictors at the time of max CAPE for each day
            rgrERAVars=np.take_along_axis(dayrgrERAdataAll,timemax2,axis=1)
            rgrERAVars=rgrERAVars[:,1,:,:,:]
            
            #save monthly data in a yearly array
            rgrERAVarsyy[iday:iday+monthlength,:,:,:]=np.abs(rgrERAVars)

            iday=monthlength+iday
        
        np.savez(outfile,
                rgrERAVarsyy = rgrERAVarsyy,
                rgrLat = rgrLat,
                rgrLon = rgrLon,
                rgrHeight = rgrHeight,
                rgdTimeDD = rgdTimeDD)
        
        #save the yearly data in an array contaning the whole time period
        rgrERAVarall[iyear:iyear+yearlength,:,:,:]=np.abs(rgrERAVarsyy)
        
        iyear=iyear+yearlength
    
    else:
        print ('The file already exists')
        data_tmp = np.load(outfile)
        rgrERAVarall[iyear:iyear+yearlength,:,:,:] = data_tmp['rgrERAVarsyy']
        rgrLat = data_tmp['rgrLat']
        rgrLon = data_tmp['rgrLon']
        rgrHeight = data_tmp['rgrHeight']
        #rgdTimeDD = pd.to_datetime(data_tmp['rgdTimeDD'])
        iyear=iyear+yearlength
    
    
# ___________________________________________________________________________
#   Load the saved CONDITIONED HAIL PDFs FOR SELECTED REGIONS
# ___________________________________________________________________________

sSavePercent=sSaveDataDir+'HailConditions_SRs_newSRH03_1990-2020_search-rad-4_qT-4_qD-26.pkl'

iFile=open(sSavePercent,'rb')
dfSave = pickle.load(iFile)
grFullHailData=dfSave['grFullHailData']
grConPercentilesSR=dfSave['grConPercentilesSR']
rgsERAVariables=dfSave['rgsERAVariables']
rgrCorrels=dfSave['rgrCorrels']
iFile.close()
   

# ___________________________________________________________________________
# ___________________________________________________________________________
#   CALCULATE AND SAVE THE HAIL PROBABILITY FROM ERA

# Train model over CONUS

# # define hail conditions based on PDFs
rgrPercentilesHailVarUS=grConPercentilesSR['CONUS']  

rgrHailSaveDir=sSaveDataDir+'HailProbabilities-ERA-5_4_qT-4_qD-26_CONUS_test_absSRH/ff'
sTargetSRUS='CONUS'
grThresholdsUS={}
rgiVarLim=np.where(np.array(['-' in rgsModelVars[ii] for ii in range(len(rgsModelVars))]) == False)[0]
for va in range(len(rgiVarLim)):
    sVarAct=rgsModelVars[rgiVarLim[va]]
    rgrV1=rgrPercentilesHailVarUS[(rgsERAVariables == sVarAct),:]
    q1=iPercentile[0]; q2=iPercentile[1]; q3=100-q2; q4=100-q1
    if (sVarAct in rgsVarMaxLim) & (sVarAct in rgsVarMinLim):
        # variable has uppe and lower limit
        grThresholdsUS[sVarAct]=[np.percentile(rgrV1,q1),np.percentile(rgrV1,q2),np.percentile(rgrV1,q3),np.percentile(rgrV1,q4)]
    elif  (sVarAct in rgsVarMaxLim):
        # variable has uppe limit
        try:
            grThresholdsUS[sVarAct]=[-99999.,-99999.,np.percentile(rgrV1,q3),np.percentile(rgrV1,q4)]
        except:
            stop()
    elif  (sVarAct in rgsVarMinLim):
        # variable has lower limit
        grThresholdsUS[sVarAct]=[np.percentile(rgrV1,q1),np.percentile(rgrV1,q2),999999., 999999.]
        
from IO_Functions import fnHailProbability
outfile=rgrHailSaveDir
rgrTotalHailProbability=fnHailProbability(rgrERAVarall[:,:,:,:],
                  rgrLatERA,
                  rgrLonERA,
                  grThresholdsUS,
                  rgrCorrels['CONUS'],
                  rgdTimeDD,
                  rgrHailSaveDir,
                  rgsERAVariables,
                  rgsFinVars,
                  rgrProb,
                  exp_name = setup_string)
    
        
print('done')
