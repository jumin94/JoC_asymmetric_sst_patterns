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
import utilities.funciones
import utilities.sst_analysis 

#---------------------------------------------LOAD DATA--------------------------------
#Data paths
path = '/datos/julia.mindlin/CMIP6_ensambles/preprocesados'
path_rean = '/datos/ERA5/mon'

#Load model data
models = ['ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CAMS-CSM1-0',
       'CanESM5', 'CESM2', 'CESM2-WACCM', 'CMCC-CM2-SR5', 'CNRM-CM6-1',
       'CNRM-ESM2-1', 'EC-Earth3', 'FGOALS-g3', 'HadGEM3-GC31-LL',
       'HadGEM3-GC31-MM', 'IITM-ESM', 'INM-CM4-8', 'INM-CM5-0',
       'KACE-1-0-G', 'MIROC6', 'MIROC-ES2L', 'MPI-ESM1-2-HR',
       'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NESM3', 'NorESM2-LM', 'NorESM2-MM',
       'TaiESM1', 'UKESM1-0-LL'
]

var = 'mon/tos'
scenarios = ['historical','ssp585']
variables = ['tos']
dato = funciones.cargo_todo(scenarios,models,path,var)

#Evaluate SST changes for all models
path_out = '/home/julia.mindlin/Tesis/JoC_paper/JoC_results'
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
asym.to_netcdf(path_out+'/sst_changes/asymmetric_SST_MEM_CMIP6.nc')              

total = ( sst_analysis.components(models[0],ssts)[2] ) * 0 
cont = 0
for model in models:
    total = total + sst_analysis.components(model,ssts)[0]
    cont += 1
    
total = total / cont
total.to_netcdf(path_out+'/sst_changes/total_SST_MEM_CMIP6.nc')              

zonal = ( sst_analysis.components(models[0],ssts)[2] ) * 0 
cont = 0
for model in models:
    zonal = zonal + sst_analysis.components(model,ssts)[1]
    cont += 1
    
zonal = zonal / cont
zonal.to_netcdf(path_out+'/sst_changes/zonally_symmetric_SST_MEM_CMIP6.nc')              


#Evaluate standard deviation of the above fields
from utilities.funciones import estadisticos
ens_mean,ens_std = estadisticos(models,ssts)
ens_std[0].to_netcdf(path_out+'/sst_changes/total_SST_std_MEM_CMIP6.nc')              
ens_std[1].to_netcdf(path_out+'/sst_changes/zonally_symmetric_SST_std_MEM_CMIP6.nc')  
ens_std[2].to_netcdf(path_out+'/sst_changes/asymmetric_SST_std_MEM_CMIP6.nc') 
