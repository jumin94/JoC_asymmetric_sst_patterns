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
var = 'tos'
dato = od.cargo_todo(experiments,models,path_data,var)

#Load remote driver change in each model (indices)
#Open indices - los mismos de siempre
gloW  = pd.read_csv(path_results+'/indices_CMIP5/GW_VB_index50_70.csv')
gw_index = gloW.iloc[:,4].values
ta = pd.read_csv(path_results+'/indices_CMIP5/TW_index.csv')
vb = pd.read_csv(path_results+'/indices_CMIP5/VB_indexGWscaled50_70.csv')
sst_C_std = pd.read_csv(path_results+'/indices_CMIP5/C_std_asym_index_DJF.csv')
sst_E_std = pd.read_csv(path_results+'/indices_CMIP5/E_std_asym_index_DJF.csv')
#Select values
TA = ta.iloc[:,4].values
VB = vb.iloc[:,5].values
SST_C_std = sst_C_std.iloc[:,2].values
SST_E_std = sst_E_std.iloc[:,2].values
#--------------------------------------------------CALCULATIONS--------------------------------------------------
#Create regression class
reg = reg_am.across_models()
var = 'tos'
#Generate regression data
reg.regression_data(dato,experiments,models,gw_index,var)
#Create folder to save data if necessary
os.chdir(path_results)
os.getcwd()
os.makedirs('sensitivity_maps_CMIP5',exist_ok=True)
os.chdir(path_results+'/sensitivity_maps_CMIP5')
os.getcwd()
os.makedirs(var,exist_ok=True)
#Perform regression
path_maps = path_results+'/sensitivity_maps_CMIP5/'+var
indices = [TA,VB,SST_C_std,SST_E_std] 
indices_names = ['TA','VB','C_index','E_index'] 
reg.perform_regression(indices,indices_names,gw_index,path_maps)
#Create .nc files from .csv generated by the regression function
file_list = utilities.csv2nc.csv_to_nc(path_maps)
