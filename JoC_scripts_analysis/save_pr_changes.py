#Evaluate precipitation sensitivity to remote drivers 

#---------------------------------------------------IMPORTS---------------------------------------------------
import numpy as np
import pandas as pd
import xarray as xr
import os, fnmatch
import glob
import utilities.open_data as od
import utilities.regression as reg_am
import utilities.csv2nc
import warnings
warnings.filterwarnings("ignore")

#---------------------------------------------------LOAD DATA--------------------------------------------------------
#Define path, model list and experiments
models = ['ACCESS1-0','ACCESS1-3','bcc-csm1-1','bcc-csm1-1-m','CanESM2','CCSM4','CMCC-CESM','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M']

experiments = ['historical','rcp85']

path_data = '/storage/silver/acrcc/co157815/fromjasmin/cmip5'
path_results = '/home/users/co157815/JoC_paper/JoC_results'
os.chdir(path_data)
os.getcwd()

#Create dictionary with data files
var = 'pr'
dato = od.cargo_todo(experiments,models,path_data,var)

#Index
gloW  = pd.read_csv(path_results+'/indices_CMIP5/GW_VB_index50_70.csv')
gw_index = gloW.iloc[:,4].values
#--------------------------------------------------CALCULATIONS--------------------------------------------------
#Create regression class
reg = reg_am.across_models()
var = 'pr'
#Generate regression data
reg.regression_data(dato,experiments,models,gw_index,var)
#Create folder to save data if necessary
os.chdir(path_results)
os.getcwd()
os.makedirs('pr_changes_CMIP5',exist_ok=True)
os.chdir(path_results+'/pr_changes_CMIP5')
os.getcwd()
os.makedirs(var,exist_ok=True)
#Save
for m in range(len(models)):
    reg.psl_change[m].to_netcdf(path_results+'/pr_changes_CMIP5'+'/pr_change_'+models[m]+'_1940-1969_2070-2099.nc')
