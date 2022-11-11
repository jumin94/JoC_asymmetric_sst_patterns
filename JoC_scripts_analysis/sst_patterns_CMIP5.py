#Evaluate model ensemble SST changes

#------------------------------------------IMPORTS---------------------------------------------
import numpy as np
import matplotlib.path as mpath
import os
import glob
import pandas as pd
import xarray as xr
import netCDF4
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.util import add_cyclic_point
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 0.5  # previous pdf hatch linewidth
import cartopy.util as cutil
import logging
import os, fnmatch
import warnings
warnings.filterwarnings("ignore")
import utilities.funciones as funciones
import utilities.sst_analysis as sst_analysis 

#---------------------------------------------LOAD DATA--------------------------------
#Data paths
path = '/storage/silver/acrcc/co157815/fromjasmin/cmip5'
path_rean = '/datos/ERA5/mon'

#Load model data
models = ['ACCESS1-0','ACCESS1-3','bcc-csm1-1','bcc-csm1-1-m','CanESM2','CCSM4','CESM1-CAM5','CMCC-CESM','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M']

var = 'tos'
scenarios = ['historical','rcp85']
variables = ['tos']
dato = funciones.cargo_todo(scenarios,models,path,var)

#Evaluate SST changes for all models
path_out = '/home/users/co157815/JoC_paper/JoC_results'
os.chdir(path_out)
os.getcwd()
#Create folder to save results
os.makedirs('sst_changes',exist_ok=True)
ssts = sst_analysis.changes_list(dato,scenarios,models,'DJF')
#Evaluate ensemble mean, zonally symmetric component and asymmetric component
asym = (sst_analysis.components(models[0],ssts)[2] ) * 0 #asym = (sst_analysis.components(models[0],ssts)[2] ) * 0 
cont = 0
for model in models:
    asym = asym + sst_analysis.components(model,ssts)[2] #asym = asym + sst_analysis.components(model,ssts)[2]
    cont += 1
    
asym = asym / cont
asym.to_netcdf(path_out+'/sst_changes/asymmetric_SST_MEM_CMIP5.nc')              

total = ( sst_analysis.components(models[0],ssts)[2] ) * 0 
cont = 0
for model in models:
    total = total + sst_analysis.components(model,ssts)[0]
    cont += 1
    
total = total / cont
total.to_netcdf(path_out+'/sst_changes/total_SST_MEM_CMIP5.nc')              

zonal = ( sst_analysis.components(models[0],ssts)[2] ) * 0 
cont = 0
for model in models:
    zonal = zonal + sst_analysis.components(model,ssts)[1]
    cont += 1
    
zonal = zonal / cont
zonal.to_netcdf(path_out+'/sst_changes/zonally_symmetric_SST_MEM_CMIP5.nc')              


#Evaluate standard deviation of the above fields
from utilities.funciones import estadisticos
ens_mean,ens_std = estadisticos(models,ssts)
ens_std[0].to_netcdf(path_out+'/sst_changes/total_SST_std_MEM_CMIP5.nc')              
ens_std[1].to_netcdf(path_out+'/sst_changes/zonally_symmetric_SST_std_MEM_CMIP5.nc')  
ens_std[2].to_netcdf(path_out+'/sst_changes/asymmetric_SST_std_MEM_CMIP5.nc') 
