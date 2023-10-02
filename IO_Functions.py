#!/usr/bin/env python
'''
    File name: OBS-HailDataPreprocessor.py
    Author: Andreas Prein
    E-mail: prein@ucar.edu
    Date created: 27.03.2017
    Date last modified: 27.03.2017

    ##############################################################
    Purpos:

    Reads in raw hail report from different sources

    Current datasets are:
    1) NOAA storm events database (http://www.spc.noaa.gov/wcm/index.html#data)
    2) BoM severe storm archive (http://www.bom.gov.au/australia/stormarchive/storm.php?stormType=hail)
    3) ESSL European severe weather database (data not public - on request)


'''



    
from pylab import *
def fnHailObs_to_HailGridded(sDataSet,   # possibilities are (NOAA, BoM, or ESSL)
                             dStartDate,
                             dStopDate,
                             sRawData=None,
                             sSaveFolder=None,
                             rMinHailSize=None):
    import datetime
    import glob
    from netCDF4 import Dataset
    import sys, traceback
    import dateutil.parser as dparser
    import string
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
    # from ipdb import set_trace as stop
    from matplotlib.patches import Polygon
    # from mpl_toolkits.basemap import Basemap, cm
    import string
    from datetime import timedelta
    # from jdcal import gcal2jd, jd2gcal
    from scipy import interpolate
    
    if sSaveFolder == None:
        sSaveFolder="/glade/work/bblanc/HailObs/"
    if rMinHailSize == None:
        print (" ")
        print ('    ________________________')
        print ("    The default hail diamter is 2.5 cm")
        rMinHailSize=2.5   # in cm 

    #############################################
    #                      read in all the hail reports
    #############################################

    if sDataSet == 'NOAA':
        
        if (sRawData == None):
            sRawData="/glade/work/bblanc/HailObs/NOAA/StormEvents_details-ftp_v1.0_d2010_c20220425.csv"
        with open(sRawData, 'rt') as f:
            reader = csv.reader(f)
            rgrHailDataRaw = list(reader)
        # generate time vector
        rgdTime=[datetime.datetime(int(float(rgrHailDataRaw[ii][1])), \
                                       int(float(rgrHailDataRaw[ii][2])), int(float(rgrHailDataRaw[ii][3])), \
                                       int(float(rgrHailDataRaw[ii][5][0:2])), int(float(rgrHailDataRaw[ii][5][3:5])), int(float(rgrHailDataRaw[ii][5][6:8]))) \
                     for ii in range(len(rgrHailDataRaw))]
        rgiYearsAll=np.array(([rgdTime[ii].year for ii in range(len(rgdTime))]))
        rgsState=np.array([rgrHailDataRaw[ii][7] for ii in range(len(rgrHailDataRaw))])
        rgiState=np.array([int(rgrHailDataRaw[ii][8]) for ii in range(len(rgrHailDataRaw))])
        rgrHailSize=np.array([float(rgrHailDataRaw[ii][10]) for ii in range(len(rgrHailDataRaw))])*2.5
        rgrHailLat=np.array([float(rgrHailDataRaw[ii][15]) for ii in range(len(rgrHailDataRaw))])
        rgrHailLon=np.array([float(rgrHailDataRaw[ii][16]) for ii in range(len(rgrHailDataRaw))])
        # create grid for hail mapping
        # load grid from netcdf
        ncid=Dataset('/glade/work/bblanc/HailObs/NOAA/HailReports-ERA5/NOAA-Hail-StormReports_gridded-75km_2010.nc', mode='r')
        rgrLatGrid=np.squeeze(ncid.variables['latitude'][:,0])
        rgrLonGrid=np.squeeze(ncid.variables['longitude'][0,:])
        ncid.close()
        # rGridSpacing=0.75 #1. #0.75  # degrees
        # rgrLonGrid=np.array(np.arange(-125.25,-63,rGridSpacing))
        # rgrLatGrid=np.array(np.arange(20.25,60,rGridSpacing))
        rgrLonGrid2D=np.asarray(([rgrLonGrid,]*rgrLatGrid.shape[0]))
        rgrLatGrid2D=np.asarray(([rgrLatGrid,]*rgrLonGrid.shape[0])).transpose()
        rgrHailDensity=np.zeros((12, len(rgrLonGrid),len(rgrLatGrid))); rgrHailDensity[:]=0
        
        for ii in range(len(hailsize)):
            print('Loadind NOAA data')
            sFileName=sSaveFolder+'NOAA-Hail-StormReports_gridded-75km_'+str(dStartDay.year)+'.nc'

            # read in the variables
            ncid=Dataset(sFileName, mode='r')
            HailObs[ii,:,:,:]=np.squeeze(ncid.variables[hailsize[ii]][:])

        ncid.close()
        # clean up
        os.system("rm "+sFileName)
        
    elif sDataSet == 'BoM':
        if (sRawData == None):
            sRawData="/glade/work/bblanc/HailObs/BoM/BoM_hail_1950-2023.csv"
        with open(sRawData, 'rb') as f:
            reader = csv.reader(f)
            rgrHailDataRaw = list(reader)
        # generate time vector
        rgdTime=[datetime.datetime(int(rgrHailDataRaw[ii][1][0:4]), \
                                       int(rgrHailDataRaw[ii][1][5:7]), int(rgrHailDataRaw[ii][1][8:10]), \
                                       int(rgrHailDataRaw[ii][1][11:13]), int(rgrHailDataRaw[ii][1][14:16]), int(rgrHailDataRaw[ii][1][17:19])) \
                     for ii in range(1,len(rgrHailDataRaw))]
        rgiYearsAll=np.array(([int(rgdTime[ii].year) for ii in range(len(rgdTime))]))
        rgsState=np.array([rgrHailDataRaw[ii][5] for ii in range(1,len(rgrHailDataRaw))])
        rgsNearestTown=np.array([rgrHailDataRaw[ii][4] for ii in range(1,len(rgrHailDataRaw))])
        rgrHailSize=np.array([float(rgrHailDataRaw[ii][6]) for ii in range(1,len(rgrHailDataRaw))])
        rgrHailLat=np.array([float(rgrHailDataRaw[ii][2]) for ii in range(1,len(rgrHailDataRaw))])
        rgrHailLon=np.array([float(rgrHailDataRaw[ii][3]) for ii in range(1,len(rgrHailDataRaw))])
        # create grid for annual hail count plot
        ncid=Dataset('/glade/work/bblanc/HailObs/BoM/HailReports-ERA5/BoM-Hail-StormReports_gridded-75km_2010.nc', mode='r')
        rgrLatGrid=np.squeeze(ncid.variables['latitude'][:,0])
        rgrLonGrid=np.squeeze(ncid.variables['longitude'][0,:])
        ncid.close()
        rgrLonGrid2D=np.asarray(([rgrLonGrid,]*rgrLatGrid.shape[0]))
        rgrLatGrid2D=np.asarray(([rgrLatGrid,]*rgrLonGrid.shape[0])).transpose()
        rgrHailDensity=np.zeros((12, len(rgrLonGrid),len(rgrLatGrid))); rgrHailDensity[:]=0
        # rGridSpacing=0.75 #1. #0.75  # degrees
        # rgrLonGrid=np.array(np.arange(111,156,rGridSpacing))
        # rgrLatGrid=np.array(np.arange(-46.5,-9,rGridSpacing))
        # rgrLonGrid2D=np.asarray(([rgrLonGrid,]*rgrLatGrid.shape[0]))
        # rgrLatGrid2D=np.asarray(([rgrLatGrid,]*rgrLonGrid.shape[0])).transpose()
        # rgrHailDensity=np.zeros((12, len(rgrLonGrid),len(rgrLatGrid))); rgrHailDensity[:]=0
        
    
    elif sDataSet == 'ESSL':
        if (sRawData == None):
            sRawData="/glade/u/home/bblanc/SyntheticHailModel/observations/ESSL-Hail/original/"
        # read the 3 hail observation files with different quallity controle levels
        rgsObsFiles=os.listdir(sRawData)
        for fi in range(len(rgsObsFiles)):
            with open(sRawData+rgsObsFiles[fi], 'rb') as f:
                reader = csv.reader(f)
                if fi == 0:
                    rgrHailDataRaw = np.array(list(reader))
                else:
                    rgrHailDataTMP=list(reader)
                    rgrUnpack = np.chararray((len(rgrHailDataTMP), 95), itemsize=100); rgrUnpack[:]=''
                    for ii in range(len(rgrHailDataTMP)):
                        try:
                            rgrUnpack[ii,:]=rgrHailDataTMP[ii]
                        except:
                            print ('        skip '+str(ii))
                    rgrHailDataRaw=np.append(rgrHailDataRaw,rgrUnpack[1:,:], axis=0)
        rgrObsHeader=rgrHailDataRaw[0,:]; rgrObsHeader=[rgrObsHeader[ii] for ii in range(len(rgrObsHeader))]
        rgrHailDataRaw=np.delete(rgrHailDataRaw, 0, 0)
        # remove non data records
        rgiNonData=np.where(rgrHailDataRaw[:,rgrObsHeader.index('TIME_EVENT')] != '')[0]
        rgrHailDataRaw=rgrHailDataRaw[rgiNonData,:]
        # hail diamater must be >= rMinHailSize
        rgiHailSize=np.where(rgrHailDataRaw[:,rgrObsHeader.index('MAX_HAIL_DIAMETER')] != '')[0]
        rgrHailDataRaw=rgrHailDataRaw[rgiHailSize,:]
        rgiHailSize=np.where(rgrHailDataRaw[:,rgrObsHeader.index('MAX_HAIL_DIAMETER')].astype('float') >= rMinHailSize)[0]
        rgrHailDataRaw=rgrHailDataRaw[rgiHailSize,:]
        # generate time vector
        rgdTime=[datetime.datetime(int(rgrHailDataRaw[ii,rgrObsHeader.index('TIME_EVENT')][0:4]), \
                                       int(rgrHailDataRaw[ii,rgrObsHeader.index('TIME_EVENT')][5:7]), int(rgrHailDataRaw[ii,rgrObsHeader.index('TIME_EVENT')][8:10]), \
                                       int(rgrHailDataRaw[ii,rgrObsHeader.index('TIME_EVENT')][11:13]), int(rgrHailDataRaw[ii,rgrObsHeader.index('TIME_EVENT')][14:16]), int(rgrHailDataRaw[ii,rgrObsHeader.index('TIME_EVENT')][17:19])) \
                     for ii in range(1,len(rgrHailDataRaw))]
        # get only hail reports that occured after 1979
        rgiYearsAll=np.array(([int(rgdTime[ii].year) for ii in range(len(rgdTime))]))
        startYear=1979
        rgiYY=np.where(rgiYearsAll >= startYear)[0]
        rgrHailDataRaw=rgrHailDataRaw[rgiYY,:]
        rgdTime=np.array(rgdTime)[rgiYY]

        rgrHailLat=rgrHailDataRaw[:,rgrObsHeader.index('LATITUDE')].astype('float')
        rgrHailLon=rgrHailDataRaw[:,rgrObsHeader.index('LONGITUDE')].astype('float')
        rgsState=rgrHailDataRaw[:,rgrObsHeader.index('COUNTRY')]
        rgrHailSize=rgrHailDataRaw[:,rgrObsHeader.index('MAX_HAIL_DIAMETER')].astype('float')
        rgiYears=np.unique([rgdTime[ii].year for ii in range(len(rgdTime))])

        # create grid for annual hail count plot
        ncid=Dataset('/glade/u/home/bblanc/SyntheticHailModel/observations/ESSL-Hail/gridded/ESSL-Hail-StormReports_gridded-75km_2009.nc', mode='r')
        rgrLatGrid=np.squeeze(ncid.variables['lat'][:,0])
        rgrLonGrid=np.squeeze(ncid.variables['lon'][0,:])
        ncid.close()
        rgrLonGrid2D=np.asarray(([rgrLonGrid,]*rgrLatGrid.shape[0]))
        rgrLatGrid2D=np.asarray(([rgrLatGrid,]*rgrLonGrid.shape[0])).transpose()
        rgrHailDensity=np.zeros((12, len(rgrLonGrid),len(rgrLatGrid))); rgrHailDensity[:]=0
        # rGridSpacing=0.75 #1. #0.75  # degrees
        # rgrLonGrid=np.array(np.arange(-11.25,45,rGridSpacing))
        # rgrLatGrid=np.array(np.arange(30,71.25,rGridSpacing))
        # rgrLonGrid2D=np.asarray(([rgrLonGrid,]*rgrLatGrid.shape[0]))
        # rgrLatGrid2D=np.asarray(([rgrLatGrid,]*rgrLonGrid.shape[0])).transpose()
        # rgrHailDensity=np.zeros((12, len(rgrLonGrid),len(rgrLatGrid))); rgrHailDensity[:]=0


    # use only reports that are within start & end date
    startYear=dStartDate.year
    rgi1959=((rgiYearsAll <= dStopDate.year) & (rgiYearsAll >= startYear) & (rgrHailSize >= rMinHailSize))
    rgdTime=np.array(rgdTime)[rgi1959]
    rgiYearsAll=rgiYearsAll[rgi1959]
    rgsState=rgsState[rgi1959]
    rgrHailSize=rgrHailSize[rgi1959]
    rgrHailLat=rgrHailLat[rgi1959]
    rgrHailLon=rgrHailLon[rgi1959]
    rgiYears=np.unique([rgdTime[ii].year for ii in range(len(rgdTime))])

    # correct time to UTC
    for he in range(len(rgrHailLat)):
        iUTC=np.round(24./360.*(rgrHailLon[he]))
        rgdTime[he]=rgdTime[he]+timedelta(hours=iUTC)


    #############################################
    # write netcdfs of gridded hail reports in annual files
    rgdTimeDays = pd.date_range(dStartDate, end=dStopDate, freq='d')
    rgrOBShail=np.zeros((len(rgdTimeDays),rgrLonGrid2D.shape[0],rgrLonGrid2D.shape[1])); rgrOBShail[:]=np.nan
    for yy in range(len(rgiYears)):
        sFile=sSaveFolder+'/'+sDataSet+'-Hail-StormReports_gridded-31km_'+str(rgiYears[yy])+'.nc_Uncompressed'
        sFileFin=sSaveFolder+'/'+sDataSet+'-Hail-StormReports_gridded-31km_'+str(rgiYears[yy])+'.nc'
        rgdDays=np.where(rgdTimeDays.year == rgiYears[yy])[0]
        if os.path.isfile(sFileFin) != 1:
            print ('    Process Year: '+str(rgiYears[yy]))
            # create matrix
            rgiYear=np.where(np.array([rgdTime[ii].year for ii in range(len(rgrHailLon))]) == rgiYears[yy])[0]
            rgdDaysAct=np.where(rgdTimeDays.year == rgiYears[yy])[0]
            rgrHailGridAct=np.zeros((len(rgdDaysAct), len(rgrLatGrid),len(rgrLonGrid))); rgrHailGridAct[:]=0
            # loop over hail reports and fill in 1 into matrix
            for hh in range(len(rgiYear)):
                iLa=np.where(np.min(np.abs(rgrLatGrid-rgrHailLat[rgiYear[hh]])) == np.abs(rgrLatGrid-rgrHailLat[rgiYear[hh]]))[0][0]
                iLo=np.where(np.min(np.abs(rgrLonGrid-rgrHailLon[rgiYear[hh]])) == np.abs(rgrLonGrid-rgrHailLon[rgiYear[hh]]))[0][0]
                # distance = (rgrLatGrid2D-rgrHailLat[rgiYear[hh]])**2 + (rgrLonGrid2D-rgrHailLon[rgiYear[hh]])**2
                # distance1D = distance.ravel()
                # distanceSORT=distance1D.argsort()[0]
                # rgrMinDist=distance1D[distanceSORT]
                # iLo=np.where(rgrMinDist == distance)[0][0]
                # iLa=np.where(rgrMinDist == distance)[1][0]
                # stop()
                iTime=rgdTime[rgiYear[hh]]
                iTime=np.where((rgdTimeDays[rgdDaysAct].year==iTime.year) & (rgdTimeDays[rgdDaysAct].month==iTime.month) & (rgdTimeDays[rgdDaysAct].day==iTime.day) == 1)[0][0]
                try:
                    rgrHailGridAct[iTime,iLa,iLo]=1# +rgrHailGridAct[iTime,iLa,iLo]
                except:
                    stop(rgiYear)

            # save data to matrix
            rgrOBShail[rgdDays,:,:]=rgrHailGridAct

            # ________________________________________________________________________
            # write the netcdf
            try:
                root_grp = Dataset(sFile, 'w', format='NETCDF4')
            except:
                stop()
            # dimensions
            root_grp.createDimension('time', None)
            root_grp.createDimension('longitude', len(rgrLonGrid))
            root_grp.createDimension('latitude', len(rgrLatGrid))
            # variables
            temp = root_grp.createVariable('Hail', 'f4', ('time','rlat','rlon',),fill_value=-99999)
            lat = root_grp.createVariable('latitude', 'f4', ('latitude','longitude',))
            lon = root_grp.createVariable('longitude', 'f4', ('latitude','longitude',))
            rlat = root_grp.createVariable('latitude', 'f4', ('latitude',))
            rlon = root_grp.createVariable('longitude', 'f4', ('longitude',))
            time = root_grp.createVariable('time', 'f8', ('time',))
            # Variable Attributes
            temp.units = "HailLarger1inch"
            temp.long_name = "HailLarger1inch"
            temp.coordinates = "lon lat"
            temp.cell_methods = "time: sum"
            # temp._FillValue = prism_nodata
            # temp._MissingValue = prism_nodata

            time.calendar = "gregorian"
            time.units = "days since "+str(startYear)+"-1-1 12:00:00"
            time.standard_name = "time"
            time.long_name = "time"
            time.axis = "T"

            lon.standard_name = "longitude"
            lon.long_name = "longitude"
            lon.units = "degrees_east"

            lat.standard_name = "latitude"
            lat.long_name = "latitude"
            lat.units = "degrees_north"

            rlon.standard_name = "grid_longitude"
            rlon.units = "degrees"
            rlon.axis = "X"

            rlat.standard_name = "grid_latitude"
            rlat.units = "degrees"
            rlat.axis = "Y"

            # write data to netcdf
            rlat[:]=rgrLatGrid
            lat[:]=rgrLatGrid2D
            rlon[:]=rgrLonGrid
            lon[:]=rgrLonGrid2D
            temp[:]=rgrHailGridAct
            time[:]=rgdDaysAct
            root_grp.close()

            # compress the netcdf file
            subprocess.Popen("nccopy -k 4 -d 1 -s "+sFile+' '+sFileFin, shell=True)
            import time
            time.sleep(10)
            subprocess.Popen("rm  "+sFile, shell=True)
            # nccopy -k 4 -d 1 -s infile outfile
        else:
            print ('    Read File: '+sFileFin)
            try:
                ncid=Dataset(sFileFin, mode='r')
            except:
                stop()
            rgrOBShail[rgdDays,:,:]=np.squeeze(ncid.variables['Hail'][:])
            ncid.close()

    return rgrOBShail, rgrLonGrid2D, rgrLatGrid2D





# _____________________________________________________________________
# _____________________________________________________________________
# _____________________________________________________________________
from pylab import *
def fnRead_ERA_Data(dStartDay,
                    dStopDay,
                    rgsERAVariables,
                    rgiDomain=None,
                    rgsERAdata=None):

    import datetime
    import glob
    from netCDF4 import Dataset
    import sys, traceback
    import dateutil.parser as dparser
    import string
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
    # from ipdb import set_trace as stop
    from matplotlib.patches import Polygon
    # from mpl_toolkits.basemap import Basemap, cm
    import string
    from datetime import timedelta
    # from jdcal import gcal2jd, jd2gcal
    from scipy import interpolate
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.gridspec as gridspec

    # read in the constant fields
    if rgsERAdata != None:
        rgsERAdata=rgsERAdata
    else:
        rgsERAdata='/glade/scratch/bblanc/ERA5_hail_model/'
    sERAconstantFields='/glade/scratch/bblanc/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc'
    ncid=Dataset(sERAconstantFields, mode='r')
    rgrLat=np.squeeze(ncid.variables['latitude'][:])
    rgrLon=np.squeeze(ncid.variables['longitude'][:])
    rgrHeight=(np.squeeze(ncid.variables['Z'][:]))/9.81
    ncid.close()

    rgiSize=rgrHeight.shape
    # make coordinates 2D
    rgrLonGrid2D=np.asarray(([rgrLon,]*rgrLat.shape[0]))
    rgrLatGrid2D=np.asarray(([rgrLat,]*rgrLon.shape[0])).transpose()
    # all lon grids have to be -180 to 180 degree
    rgi180=np.where(rgrLonGrid2D > 180)
    rgrLonGrid2D[rgi180]=rgrLonGrid2D[rgi180]-360.

    # generate the time vector
    rgdTime1H = pd.date_range(dStartDay, end=dStopDay, freq='1h')
    iYears=np.unique(rgdTime1H.year)
    iMonth=np.unique(rgdTime1H.month) 
    # check if we only want to read a subregion
    if rgiDomain != None:
        iloW=rgiDomain[3]
        iloE=rgiDomain[1]
        ilaN=rgiDomain[0]
        ilaS=rgiDomain[2]
    else:
        iloW=rgrLonGrid2D.shape[1]
        iloE=0
        ilaN=0
        ilaS=rgrLonGrid2D.shape[0]

    # __________________________________________________________________
    # start reading the ERA data
    rgrERAdataAll=np.zeros((len(iMonth), len(rgsERAVariables),ilaS-ilaN,iloW-iloE)); rgrERAdataAll[:]=np.nan

    print ('    Loading ERA-5 ') 
    rgdTimeDD = pd.date_range(dStartDay, end=dStopDay, freq='d')
    dDate=rgdTimeDD[1]
    # loop over variables
    sFileName=rgsERAdata +str(dDate.year)+str("%02d" % (dDate.month))+'_ERA-5_HailPredictors_from-RDA'+'.nc'
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

    rgrLonGrid2D=rgrLonGrid2D[ilaN:ilaS,iloE:iloW]
    rgrLatGrid2D=rgrLatGrid2D[ilaN:ilaS,iloE:iloW]

    return rgrERAdataAll, rgrLonGrid2D, rgrLatGrid2D









# _________________________________________________________________
# _________________________________________________________________

def fnHailProbability(rgrData,
                      rgrLat,
                      rgrLon,
                      grThresholds,
                      Correl_Paths,
                      rgdTime,
                      sSaveDataDir,
                      rgsDataVariables,
                      rgsModelVars,
                      rgrProb,
                      exp_name):
    import datetime
    import glob
    from netCDF4 import Dataset
    import sys, traceback
    import dateutil.parser as dparser
    import string
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
    # from ipdb import set_trace as stop
    from matplotlib.patches import Polygon
    # from mpl_toolkits.basemap import Basemap, cm
    import string
    from datetime import timedelta
    # from jdcal import gcal2jd, jd2gcal
    from scipy import interpolate
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.gridspec as gridspec
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.gridspec as gridspec
    from pdb import set_trace as stop

    rgsMaxVars=['CAPEmax',
                'SRH03',
                'VS03'
               ]
    rgrHailProb=np.zeros((rgrData.shape[0],len(rgsModelVars),rgrData.shape[2],rgrData.shape[3])); rgrHailProb[:]=0
    rgiSizeTMP=rgrData[:,0,:,:].shape
    from Statistic_Functions import transfere
    ii=0
    #iterate over the single variables
    for va in range(len(rgsMaxVars)):
        sVarTMP=rgsMaxVars[va]
        
        rgrHailProb[:,rgsModelVars.index(sVarTMP),:,:]=transfere(rgrData[:,rgsDataVariables.tolist().index(sVarTMP),:,:], grThresholds[sVarTMP])
        print ('        Hail Probability for '+sVarTMP)
        
    # do the combined variables
    for va in range(len(Correl_Paths.keys())):
        
        sVarTMP=list(Correl_Paths.keys())[va]
        print ('        Hail Probability for '+sVarTMP)
        rgsVarsAct=sVarTMP.split('-')
        xx=np.where(np.array(rgsDataVariables) == rgsVarsAct[0])[0][0]
        yy=np.where(np.array(rgsDataVariables) == rgsVarsAct[1])[0][0]
        rgsSteps=list(Correl_Paths[sVarTMP].keys())
        for pe in range(len(rgsSteps)):
            rgrPath = Correl_Paths[sVarTMP][str(rgsSteps[pe])]

            # plt.plot(rgrPath.vertices[:,0],rgrPath.vertices[:,1]); plt.show()

            # rgrPointsL=np.append(rgrXVals[:,None],rgrYVals[:,0,pe][:,None], axis=1)
            # rgrPointsH=np.append(rgrXVals[:,None],rgrYVals[:,1,pe][:,None], axis=1)
            # rgrPoints=np.append(rgrPointsH,rgrPointsL[::-1,:],axis=0)
            # rgrPath = mplPath.Path(rgrPoints)
            try:
                rgiEnvTMP=(rgrPath.contains_points(np.array([rgrData[:,xx,:,:].ravel(),rgrData[:,yy,:,:].ravel()]).transpose()))
            except:
                stop()
            rgiEnvTMP=rgiEnvTMP.reshape(rgrData.shape[0], rgrData.shape[2], rgrData.shape[3])
            TMP=rgrHailProb[:,rgsModelVars.index(sVarTMP),:,:]
            TMP[rgiEnvTMP]=rgrProb[pe]
            rgrHailProb[:,rgsModelVars.index(sVarTMP),:,:]=TMP
            # plt.contourf(rgrHailProb[0,rgsModelVars.index(sVarTMP),:,:], levels=np.linspace(0,1,10),extemd='both'); plt.show()

    rgrTotalHailProbability=np.prod(rgrHailProb, axis=1)
    # set nan to 0
    rgrTotalHailProbability[np.isnan(rgrTotalHailProbability)]=0

    # ________________________________________________________________________
    # write out annual NetCDFs
    rgrTimeNetCDF=pd.date_range(datetime.datetime(1959, 1, 1,12), end=rgdTime[-1], freq='d')
    rgiTime=np.array(range(len(rgrTimeNetCDF)))
    rgiShape=rgrHailProb.shape
    nodata = -9999
    iYears=np.unique(rgdTime.year)

    for yy in range(len(iYears)):
        rgrDaysAct=np.where(rgdTime.year ==iYears[yy])
        rgiTimeAct=rgiTime[np.where(rgrTimeNetCDF.year ==iYears[yy])]

        sFile=sSaveDataDir+'/'+str(iYears[yy])+'_HailProbabilities_ERA-5_hourly.nc-Uncompressed'
        sFileFin=sSaveDataDir+'/'+str(iYears[yy])+'_HailProbabilities_ERA-5_hourly.nc'
        if os.path.isdir(sSaveDataDir) != 1:
            subprocess.call(["mkdir","-p",sSaveDataDir])
        root_grp = Dataset(sFileFin, 'w', format='NETCDF4')
        # dimensions
        root_grp.createDimension('time', None)
        root_grp.createDimension('rlon', rgiShape[3])
        root_grp.createDimension('rlat', rgiShape[2])
        root_grp.createDimension('vars', rgiShape[1])
        # variables
        # CAPE = root_grp.createVariable('CAPE', 'f4', ('time','rlat','rlon',),fill_value=nodata)
        # VS0_3 = root_grp.createVariable('VS0_3', 'f4', ('time','rlat','rlon',),fill_value=nodata)
        # LRml = root_grp.createVariable('LRml', 'f4', ('time','rlat','rlon',),fill_value=nodata)
        # ThettaEpbl_MRml = root_grp.createVariable('ThettaEpbl-MRml', 'f4', ('time','rlat','rlon',),fill_value=nodata)
        # TCW_MRpbl = root_grp.createVariable('TCW-MRpbl', 'f4', ('time','rlat','rlon',),fill_value=nodata)
        # CAPE_CIN = root_grp.createVariable('CAPE-CIN', 'f4', ('time','rlat','rlon',),fill_value=nodata)
        # FreezingLH_MRpbl = root_grp.createVariable('FreezingLH-MRpbl', 'f4', ('time','rlat','rlon',),fill_value=nodata)
        # RHbML_RHpbl = root_grp.createVariable('RHbML-RHpbl', 'f4', ('time','rlat','rlon',),fill_value=nodata)
        HailProb= root_grp.createVariable('HailProb', 'f4', ('time','rlat','rlon',),fill_value=nodata)
        HailProbVars= root_grp.createVariable('HailProbVars', 'f4', ('time','vars','rlat','rlon',),fill_value=nodata)
        lat = root_grp.createVariable('lat', 'f4', ('rlat','rlon',))
        lon = root_grp.createVariable('lon', 'f4', ('rlat','rlon',))
        time = root_grp.createVariable('time', 'f8', ('time',))
        # Variable Attributes

        HailProbVars.long_name = 'convective available potential energy'
        HailProbVars.coordinates = "lon lat"
        HailProbVars.longname = ', '.join(rgsModelVars)
        HailProbVars.unit = "%"

        # CAPE.long_name = 'convective available potential energy'
        # CAPE.coordinates = "lon lat"
        # CAPE.cell_methods = "time: max"

        # VS0_3.long_name = 'vector shear between the surface and 6 km'
        # VS0_3.coordinates = "lon lat"
        # VS0_3.cell_methods = "time: max"

        # LRml.long_name = 'mid level lapse rate'
        # LRml.coordinates = "lon lat"
        # LRml.cell_methods = "time: mean"

        # ThettaEpbl_MRml.long_name = 'relation between Theta-E and mid level mixing ratio'
        # ThettaEpbl_MRml.coordinates = "lon lat"
        # ThettaEpbl_MRml.cell_methods = "time: mean"

        # TCW_MRpbl.long_name = 'relation between precipitable water and boundary layer mixing ratio'
        # TCW_MRpbl.coordinates = "lon lat"
        # TCW_MRpbl.cell_methods = "time: mean"

        # CAPE_CIN.long_name = 'relation between CAPE and CIN'
        # CAPE_CIN.coordinates = "lon lat"
        # CAPE_CIN.cell_methods = "time: max"

        # FreezingLH_MRpbl.long_name = 'relation between freezing level height and boundary layer mixing ratio'
        # FreezingLH_MRpbl.coordinates = "lon lat"
        # FreezingLH_MRpbl.cell_methods = "time: max"

        # RHbML_RHpbl.long_name = 'Relation between relative humidity below the melting level and the PBL'
        # RHbML_RHpbl.coordinates = "lon lat"
        # RHbML_RHpbl.cell_methods = "time: mean"

        HailProb.units = "%"
        HailProb.long_name = 'Total hail probability'
        HailProb.coordinates = "lon lat"
        HailProb.cell_methods = "time: mean"

        time.calendar = "gregorian"
        time.units = "days since 1959-1-1 12:00:00"
        time.standard_name = "time"
        time.long_name = "time"
        time.axis = "T"

        lon.standard_name = "longitude"
        lon.long_name = "longitude"
        lon.units = "degrees_east"

        lat.standard_name = "latitude"
        lat.long_name = "latitude"
        lat.units = "degrees_north"

        # write data to netcdf
        lat[:]=rgrLat
        lon[:]=rgrLon
        HailProbVars[:]=np.squeeze(rgrHailProb[rgrDaysAct,:,:,:])
        # CAPE[:]=np.squeeze(rgrHailProb[rgrDaysAct,grThresholds.keys().index('CAPE'),:,:])
        # VS0_3[:]=np.squeeze(rgrHailProb[rgrDaysAct,rgsModelVars.index('VS0_3'),:,:])
        # LRml[:]=np.squeeze(rgrHailProb[rgrDaysAct,rgsModelVars.index('LRml'),:,:])
        # ThettaEpbl_MRml[:]=np.squeeze(rgrHailProb[rgrDaysAct,rgsModelVars.index('ThettaEpbl-MRml'),:,:])
        # TCW_MRpbl[:]=np.squeeze(rgrHailProb[rgrDaysAct,rgsModelVars.index('TCW-MRpbl'),:,:])
        # CAPE_CIN[:]=np.squeeze(rgrHailProb[rgrDaysAct,rgsModelVars.index('CAPE-CIN'),:,:])
        # FreezingLH_MRpbl[:]=np.squeeze(rgrHailProb[rgrDaysAct,rgsModelVars.index('FreezingLH-MRpbl'),:,:])
        # RHbML_RHpbl[:]=np.squeeze(rgrHailProb[rgrDaysAct,rgsModelVars.index('RHbML-RHpbl'),:,:])
        HailProb[:]=np.squeeze(rgrTotalHailProbability[rgrDaysAct,:,:])
        time[:]=rgiTimeAct
        root_grp.close()

        # # compress the netcdf file
        # subprocess.Popen("nccopy -k 4 -d 1 -s "+sFile+' '+sFileFin, shell=True)
        # import time
        # time.sleep(10)
        # subprocess.Popen("rm  "+sFile, shell=True)
        # # nccopy -k 4 -d 1 -s infile outfile
        
    return rgrTotalHailProbability





# _____________________________________________________________________
# _____________________________________________________________________
# _____________________________________________________________________

def fn_Read_Soundings(sROdata,
                      rgiStatNr,
                      rgdTime,
                      rgsVar):

    rgrRSdataAll=np.zeros((len(rgdTime), len(rgiStatNr), len(rgsVar))); rgrRSdataAll[:]=np.nan

    # generate time vector
    dStartDate=rgdTime[0]; dStartDate = dStartDate.replace(hour=0)
    dStopDate=rgdTime[-1]; dStopDate = dStopDate.replace(hour=23)
    rgdTime = pd.date_range(dStartDate, end=dStopDate, freq='m')
    rgdTimeRS = pd.date_range(dStartDate, end=dStopDate, freq='12H')

    dcRSdata={}
    # loop over stations
    for st in range(len(rgiStatNr)):
        sStat=str(rgiStatNr[st])
        # matrix that hold the data
        rgrData=np.zeros((len(rgdTimeRS),11,200)); rgrData[:]=np.nan
        rgrCIN=np.zeros((len(rgdTimeRS))); rgrCIN[:]=np.nan
        rgrCAPE=np.copy(rgrCIN)
        rgrPWT=np.copy(rgrCIN)

        # loop over the demanded years
        rgsFileName=glob.glob(sROdata+sStat+'*')
        if len(rgsFileName) != 0:
            iTime = 0
            for yy in range(len(rgsFileName)):
                sActFileName=sorted(rgsFileName)[yy]
                iSkip=0
                if sActFileName != '':
                    ncid = Dataset(sActFileName, mode='r') # open the netcdf
                    ncid.variables['time'][:]
                    rgdTimeAct=ncid.variables['time'][:]
                    sDaysSince=ncid.variables['time'].units
                    sDaysSince = dparser.parse(sDaysSince,fuzzy=True)
                    sTimestep=ncid.variables['time'].units.split(' ', 1)[0]
                    if iTime == 0:
                        sVariabbles=ncid.variables['data'].long_name.split(', ')
                        sUnits=ncid.variables['data'].units.split(', ')
                        # cobine variable with units
                        for va in range(len(sVariabbles)):
                            sVariabbles[va]=sVariabbles[va]+' ['+sUnits[va]+']'
                        dcMetaDat={'ID':getattr(ncid,'stationID'),
                                   'lat':getattr(ncid,'stationLat'),
                                   'lon':getattr(ncid,'stationLon'),
                                   'elev':getattr(ncid,'stationHeight')}
                        iTime = 1
                    ncid.close()
                    if str(sTimestep) == 'days':
                        rDaySlizes=1.
                    elif str(sTimestep) == 'hours':
                        rDaySlizes=24.
                    elif str(sTimestep) == 'minutes':
                        rDaySlizes=24.*60.
                    elif str(sTimestep) == 'seconds':
                        rDaySlizes=24.*60.*60.
                    else:
                        print ('        The timestep in the netcdf file is not implemented!')
                        exit
                    dFileStart = sDaysSince + datetime.timedelta(days=rgdTimeAct[0]/rDaySlizes)
                    dFileEnd = sDaysSince + datetime.timedelta(days=rgdTimeAct[-1]/rDaySlizes)
                    # since this is not the actual frequency of the data we now derive it
                    rFrequency=float(rgdTimeAct.shape[0])/float((dFileEnd-dFileStart).days+1) #[slices per day]
                    # quickfix which removes the time-zone information
                    dFileStart=dFileStart.replace(tzinfo=None)
                    dFileEnd=dFileEnd.replace(tzinfo=None)
                    if dStartDate >= dFileStart and dStartDate <= dFileEnd:
                        # the time period is covered but we need an offset
                        idelta = dStartDate-dFileStart
                        rOffset=int(idelta.days*rFrequency+idelta.seconds/(60.*60.*24.)*rFrequency)
                    elif dStartDate <= dFileStart:
                        # we can read from the beginning of the file
                        rOffset=int(0)
                    else:
                        iSkip=1

                    if dStopDate >= dFileEnd:
                        # we can read untill the end of the file
                        rCount=int(len(rgdTimeAct)-1.)
                    elif dStopDate < dFileEnd and dStopDate > dFileStart:
                        # we stop before the end of the file
                        idelta = dStopDate-dFileStart
                        rCount=int(idelta.days*rFrequency+idelta.seconds/(60.*60.*24.)*rFrequency)
                    else:
                        iSkip=1

                    if iSkip == 0:
                        # read in the data
                        print ('        '+sActFileName)
                        ncid = Dataset(sActFileName, mode='r') # open the netcdf
                        rgrDataAct = np.array(ncid.variables['data'][rOffset:rCount+1,:,:])
                        rgrCINact= ncid.variables['CIN'][rOffset:rCount+1]
                        rgrCAPEact= ncid.variables['CAPE'][rOffset:rCount+1]
                        rgrPWTact= ncid.variables['PRECW'][rOffset:rCount+1]
                        # Check for missing values
                        vars=ncid.variables
                        grData=vars['data']
                        dcVarAtt=grData.ncattrs()
                        if '_missingvalue' in dcVarAtt:
                            rNaNval=ncid.variables['data']._MissingValue
                        elif '_FillValue' in dcVarAtt:
                            rNaNval=ncid.variables['data']._FillValue
                        else:
                            rNaNval='none'
                        ncid.close()
                        if rNaNval != 'none':
                            rgrDataAct[np.array(rgrDataAct)==np.array(rNaNval)] = np.nan
                        iStaTime=np.where(dFileStart.year == rgdTimeRS.year)[0]
                        rgrData[iStaTime,:,:]=rgrDataAct
                        rgrCIN[iStaTime]=rgrCINact
                        rgrCAPE[iStaTime]=rgrCAPEact
                        rgrPWT[iStaTime]=rgrPWTact
                    # exit for loop if we have all the data we needed
                    if dStopDate < dFileEnd:
                        break

            #  store the data in pandas pannel
            from HelperFunctions import TemporalFrequency
            try:
                sTempFreq=TemporalFrequency(rFrequency)
            except:
                stop()
            rgiSize=rgrData.shape
            try:
                dfData = pd.Panel(rgrData, items=rgdTimeRS,  major_axis=sVariabbles, minor_axis=range(0,rgiSize[2]))
            except:
                stop()
            dfCIN=pd.Series(rgrCIN, index=rgdTimeRS)
            dfCAPE=pd.Series(rgrCAPE, index=rgdTimeRS)
            dfPWT=pd.Series(rgrPWT, index=rgdTimeRS)
            dcRSdata[sStat]=[dfData,dfCIN,dfCAPE,dfPWT,dcMetaDat]


    ######################################
    # derive the hail data for all stations
    rgsVariables=['CAPE','RHbf','RHml','VS0_1','VS0_3','VS0_6','VS6_12','PW','FLH']
    rgrHailPara=np.zeros((len(rgdTimeRS), len(dcRSdata.keys()),len(rgsVariables))); rgrHailPara[:]=np.nan
    rgrLon=np.zeros((len(dcRSdata.keys()))); rgrLon[:]=np.nan
    rgrLat=np.copy(rgrLon)
    rgrElev=np.copy(rgrLon)
    rgsStatID=['']*len(dcRSdata.keys())
    rgsStatNr=['']*len(dcRSdata.keys())
    for st in range(len(dcRSdata.keys())):
         # derive station methadata
         rgrLon[st]=dcRSdata[dcRSdata.keys()[st]][4]['lon']
         rgrLat[st]=dcRSdata[dcRSdata.keys()[st]][4]['lat']
         rgrElev[st]=dcRSdata[dcRSdata.keys()[st]][4]['elev']
         rgsStatID[st]=dcRSdata[dcRSdata.keys()[st]][4]['ID']
         rgsStatNr[st]=dcRSdata.keys()[st]
         # derive heights
         rgrDataAct=dcRSdata[dcRSdata.keys()[st]][0]
         rgrHeight=np.array(rgrDataAct.major_xs('HGHT [m]'))-dcRSdata[dcRSdata.keys()[st]][4]['elev'] # to height above surface
         rgrPRES=np.array(rgrDataAct.major_xs('PRES [hPa]'))
         rgrTEMP=np.array(rgrDataAct.major_xs('TEMP [C]'))
         rgrMR=np.array(rgrDataAct.major_xs('MIXR [g/kg]'))/1000.  # convert to kg/kg
         rgrWD=np.array(rgrDataAct.major_xs('DRCT [deg]'))
         rgrWS=np.array(rgrDataAct.major_xs('SKNT [knot]'))*0.514444 # convert to m/s
         # calculate U and V component
         rgrU = -rgrWS * np.sin(rgrWD * np.pi/180.0)
         rgrV = -rgrWS * np.cos(rgrWD * np.pi/180.0)

         # calcualte satturation mixing ratio
         from thermodynamics import VaporPressure, MixRatio
         rgrSVP=VaporPressure(rgrTEMP, phase="liquid")
         rgrMRs=MixRatio(rgrSVP,rgrPRES*100.)

         for va in range(len(rgsVariables)):
             if rgsVariables[va] == 'CAPE':
                 # derive cape instead of taking it from the sounding data
                 for ii in range(rgrPRES.shape[1]):
                     rgiReal=(np.where(~np.isnan(rgrTEMP[:,ii]) == True))[0]
                     try:
                         CAPE,CIN=wrf.cape_3d(rgrPRES[rgiReal,ii], (rgrTEMP[rgiReal,ii]+273.15), rgrMR[rgiReal,ii], rgrHeight[rgiReal,ii], dcRSdata[dcRSdata.keys()[st]][4]['elev'] , rgrPRES[rgiReal[0],ii], ter_follow=False, missing=-99999., meta=False)
                         rgrHailPara[ii,st,va]=np.max(CAPE)
                     except:
                         # take CAPE from sounding report
                         rgrHailPara[ii,st,va]=dcRSdata[dcRSdata.keys()[st]][2][ii]
                 # # Read in CAPE provided with the sounding
                 # rgrHailPara[:,st,va]=dcRSdata[dcRSdata.keys()[st]][2]
             if rgsVariables[va] == 'RHbf':
                  for tt in range(len(rgdTimeRS)):
                       if np.sum(~np.isnan(rgrTEMP[:,tt])) > 0:
                            iFL=np.where((np.nanmin(abs(rgrTEMP[:,tt])) == abs(rgrTEMP[:,tt])))[0][0]
                            rgrMRbf=rgrMR[0:iFL,tt]
                            rgrMRsBF=rgrMRs[0:iFL,tt]
                            # wheighted according to height
                            rgrHGTcopy=np.copy(rgrHeight[:,tt])
                            rgrHGTcopy=rgrHGTcopy[1:]-rgrHGTcopy[:-1]
                            rgrHGTcopy=rgrHGTcopy[0:iFL]
                            rgrHGTsum=np.nansum(rgrHGTcopy, axis=0)
                            rgrHailPara[tt,st,va]=((np.nansum(rgrMRbf*rgrHGTcopy,axis=0)/rgrHGTsum)/(np.nansum(rgrMRsBF*rgrHGTcopy, axis=0)/rgrHGTsum))
             if rgsVariables[va] == 'RHml':
                  for tt in range(len(rgdTimeRS)):
                      if np.sum(~np.isnan(rgrPRES[:,tt])) > 0:
                          i700=np.where((np.nanmin(abs(rgrPRES[:,tt]-700)) == abs(rgrPRES[:,tt]-700)))[0][0]
                          i500=np.where((np.nanmin(abs(rgrPRES[:,tt]-500)) == abs(rgrPRES[:,tt]-500)))[0][0]
                          rgrMRbf=rgrMR[i700:i500,tt]
                          rgrMRsBF=rgrMRs[i700:i500,tt]
                          # wheighted according to height
                          rgrHGTcopy=np.copy(rgrHeight[:,tt])
                          rgrHGTcopy=rgrHGTcopy[1:]-rgrHGTcopy[:-1]
                          rgrHGTcopy=rgrHGTcopy[i700:i500]
                          rgrHGTsum=np.nansum(rgrHGTcopy, axis=0)
                          rgrHailPara[tt,st,va]=((np.nansum(rgrMRbf*rgrHGTcopy,axis=0)/rgrHGTsum)/(np.nansum(rgrMRsBF*rgrHGTcopy, axis=0)/rgrHGTsum))
             if rgsVariables[va] == 'VS0_1':
                  for tt in range(len(rgdTimeRS)):
                       if np.sum(~np.isnan(rgrHeight[:,tt])) > 0:
                            try:
                                 i1000=np.where((np.nanmin(abs(rgrHeight[:,tt]-1000)) == abs(rgrHeight[:,tt]-1000)))[0][0]
                                 rgrHailPara[tt,st,va]=((rgrU[i1000,tt]-rgrU[:,tt][~np.isnan(rgrU[:,tt])][0])**2+\
                                                            (rgrV[i1000,tt]-rgrV[:,tt][~np.isnan(rgrV[:,tt])][0])**2)**0.5
                            except:
                                 continue
             if rgsVariables[va] == 'VS0_3':
                  for tt in range(len(rgdTimeRS)):
                       if np.sum(~np.isnan(rgrHeight[:,tt])) > 0:
                            try:
                                 i3000=np.where((np.nanmin(abs(rgrHeight[:,tt]-3000)) == abs(rgrHeight[:,tt]-3000)))[0][0]
                                 rgrHailPara[tt,st,va]=((rgrU[i3000,tt]-rgrU[:,tt][~np.isnan(rgrU[:,tt])][0])**2+\
                                                            (rgrV[i3000,tt]-rgrV[:,tt][~np.isnan(rgrV[:,tt])][0])**2)**0.5
                            except:
                                 continue
             if rgsVariables[va] == 'VS0_6':
                  for tt in range(len(rgdTimeRS)):
                       if np.sum(~np.isnan(rgrHeight[:,tt])) > 0:
                            try:
                                 i6000=np.where((np.nanmin(abs(rgrHeight[:,tt]-6000)) == abs(rgrHeight[:,tt]-6000)))[0][0]
                                 rgrHailPara[tt,st,va]=((rgrU[i6000,tt]-rgrU[:,tt][~np.isnan(rgrU[:,tt])][0])**2+\
                                                            (rgrV[i6000,tt]-rgrV[:,tt][~np.isnan(rgrV[:,tt])][0])**2)**0.5
                            except:
                                 continue
             if rgsVariables[va] == 'VS6_12':
                  for tt in range(len(rgdTimeRS)):
                       if np.sum(~np.isnan(rgrHeight[:,tt])) > 0:
                            try:
                                 i6000=np.where((np.nanmin(abs(rgrHeight[:,tt]-6000)) == abs(rgrHeight[:,tt]-6000)))[0][0]
                                 i12000=np.where((np.nanmin(abs(rgrHeight[:,tt]-12000)) == abs(rgrHeight[:,tt]-12000)))[0][0]
                                 rgrHailPara[tt,st,va]=((rgrU[i12000,tt]-rgrU[i6000,tt])**2+\
                                                            (rgrV[i12000,tt]-rgrV[i6000,tt])**2)**0.5
                            except:
                                 continue
             if rgsVariables[va] == 'PW':
                  rgrHailPara[:,st,va]=np.array(dcRSdata[dcRSdata.keys()[st]][3])
             if rgsVariables[va] == 'FLH':
                 for tt in range(len(rgdTimeRS)):
                     if np.sum(~np.isnan(rgrTEMP[:,tt])) > 0:
                         rgiFL=np.where((np.nanmin(abs(rgrTEMP[:,tt])) == abs(rgrTEMP[:,tt])))[0][0]
                         rgrHailPara[tt,st,va]=rgrHeight[rgiFL,tt]
    # reform for daily statisitcs calculation
    rgrHailPara=np.reshape(rgrHailPara, (len(rgdTimeRS)/2,2, len(dcRSdata.keys()),len(rgsVariables)))

    # remove -99999 values
    rgrHailPara[(rgrHailPara==-99999)]=np.nan
    # for CAPE calculate the maximum with an overlap of +-1/2 day
    rgrCAPEmax=np.zeros((rgrRSdataAll.shape[0], rgrRSdataAll.shape[1])); rgrCAPEmax[:]=np.nan
    for dd in range(rgrRSdataAll.shape[0]):
        ddb=dd-1
        dde=dd+1
        if dd == 0:
            ddb=0
        elif dd == rgrRSdataAll.shape[0]-1:
            dde=rgrRSdataAll.shape[0]-1
        rgrCAPEmax[dd,:]=np.max(np.append(rgrHailPara[dd,:,:,rgsVariables.index('CAPE')], \
                                              [rgrHailPara[ddb,1,:,rgsVariables.index('CAPE')][0],rgrHailPara[dde,0,:,rgsVariables.index('CAPE')][0]]))
    rgrRSdataAll[:,:,rgsVar.index('CAPE')]=rgrCAPEmax
    # rgrRSdataAll[:,:,rgsVar.index('CAPE')]=np.max(rgrHailPara[:,:,:,rgsVariables.index('CAPE')], axis=1)
    rgrRSdataAll[:,:,rgsVar.index('RHbML')]=np.mean(rgrHailPara[:,:,:,rgsVariables.index('RHbf')], axis=1)*100.
    rgrRSdataAll[:,:,rgsVar.index('RHml')]=np.mean(rgrHailPara[:,:,:,rgsVariables.index('RHml')], axis=1)*100.
    rgrRSdataAll[:,:,rgsVar.index('VS0_1')]=np.max(rgrHailPara[:,:,:,rgsVariables.index('VS0_1')], axis=1)
    rgrVS0_3=np.max(rgrHailPara[:,:,:,rgsVariables.index('VS0_3')], axis=1); rgrVS0_3[(rgrVS0_3==-99999)]=np.nan
    rgrRSdataAll[:,:,rgsVar.index('VS0_6')]=np.max(rgrHailPara[:,:,:,rgsVariables.index('VS0_6')], axis=1)
    rgrRSdataAll[:,:,rgsVar.index('VS6_12')]=np.max(rgrHailPara[:,:,:,rgsVariables.index('VS6_12')], axis=1)
    rgrRSdataAll[:,:,rgsVar.index('TCW')]=np.mean(rgrHailPara[:,:,:,rgsVariables.index('PW')], axis=1)
    rgrRSdataAll[:,:,rgsVar.index('FreezingLH')]=np.mean(rgrHailPara[:,:,:,rgsVariables.index('FLH')], axis=1)

    rgrRSdataAll[(rgrRSdataAll==-99999)]=np.nan

    return rgrRSdataAll, rgrLat, rgrLon

    stop()
