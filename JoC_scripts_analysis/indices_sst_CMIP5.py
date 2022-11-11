#Indices ------------------------------------------------------------------------------------------

#This code evaluates all the remote driver indices needed for a storyline analysis of the SH atmospheric ciruclation in CMIP6. This is: 

#Asymmetric SST anomalies evaluated as Central and Eastern Pacific Warming patterns 

#Climatological change is defined as the is ensemble mean SSP5-8.8 2070-2099 vs. ensemble mean historical simulations 1940-1969 except for VB delay where 1950-1979 is used. 

#Author: Julia Mindlin

#Global warming
import numpy as np
import pandas as pd
import xarray as xr
import utilities.index_module as im
from utilities.index_module import cargo_todo
import os

models = ['ACCESS1-0','ACCESS1-3','bcc-csm1-1','bcc-csm1-1-m','CanESM2','CCSM4','CESM1-CAM5','CMCC-CESM','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M']

experiments = ['historical','rcp85']

path_data = '/storage/silver/acrcc/co157815/fromjasmin/cmip5' 
path_results = '/home/users/co157815/JoC_paper/JoC_results'
os.chdir(path_data)
os.getcwd()
os.makedirs('indices_CMIP5',exist_ok=True)


#Asymmetric Sea Surface temperature change
#SST changes
var = 'tos'
variables = ['tos']
dato = cargo_todo(path_data,experiments,models,var)
seasons = ['DJF'] 

print(dato)

#Eastern STD 
box = [0,-10,260,290]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    #components, 0: full, 1: symmetric, 2: asymmetric
    D_SST_E_std = im.sst_index_asym(models,box,ssts)
    D_SST_E_std = {season: D_SST_E_std}
    D_SST_E_std = pd.DataFrame(D_SST_E_std)
    D_SST_E_std.insert(0,"Modelo", models,True)
    D_SST_E_std.to_csv(path_results+'/indices_CMIP5/E_std_asym_index_'+season+'.csv',float_format='%g')

#Central STD
box = [5,-5,180,250]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    D_SST_C_std = im.sst_index_asym(models,box,ssts)
    D_SST_C_std = {season: D_SST_C_std}
    D_SST_C_std = pd.DataFrame(D_SST_C_std)
    D_SST_C_std.insert(0,"Modelo", models,True)
    D_SST_C_std.to_csv(path_results+'/indices_CMIP5/C_std_asym_index_'+season+'.csv',float_format='%g')


#Niño1.2
box = [0,-10,270,280]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    #components, 0: full, 1: symmetric, 2: asymmetric
    D_SST_12 = im.sst_index_asym(models,box,ssts)
    D_SST_12 = {season: D_SST_12}
    D_SST_12 = pd.DataFrame(D_SST_12)
    D_SST_12.insert(0,"Modelo", models,True)
    D_SST_12.to_csv(path_results+'/indices_CMIP5/Nino12_asym_index_'+season+'.csv',float_format='%g')

#Niño4
box = [5,-5,160,210]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    D_SST_4 = im.sst_index_asym(models,box,ssts)
    D_SST_4 = {season: D_SST_4}
    D_SST_4 = pd.DataFrame(D_SST_4)
    D_SST_4.insert(0,"Modelo", models,True)
    D_SST_4.to_csv(path_results+'/indices_CMIP5/Nino4_asym_index_'+season+'.csv',float_format='%g')

#C index
D_SST_C = 1.7*D_SST_4.iloc[:,1] - 0.1*D_SST_12.iloc[:,1]
D_SST_C.to_csv(path_results+'/indices_CMIP5/C_asym_index_'+season+'.csv',float_format='%g')
#E index
D_SST_E = D_SST_12.iloc[:,1] - 0.5*D_SST_4.iloc[:,1] 
D_SST_E.to_csv(path_results+'/indices_CMIP5/E_asym_index_'+season+'.csv',float_format='%g')


#Niño1.2
box = [0,-10,270,280]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    #components, 0: full, 1: symmetric, 2: asymmetric
    D_SST_12 = im.sst_index_full(models,box,ssts)
    D_SST_12 = {season: D_SST_12}
    D_SST_12 = pd.DataFrame(D_SST_12)
    D_SST_12.insert(0,"Modelo", models,True)
    D_SST_12.to_csv(path_results+'/indices_CMIP5/Nino12_full_index_'+season+'.csv',float_format='%g')

#Niño4
box = [5,-5,160,210]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    D_SST_4 = im.sst_index_full(models,box,ssts)
    D_SST_4 = {season: D_SST_4}
    D_SST_4 = pd.DataFrame(D_SST_4)
    D_SST_4.insert(0,"Modelo", models,True)
    D_SST_4.to_csv(path_results+'/indices_CMIP5/Nino4_full_index_'+season+'.csv',float_format='%g')

#C index
D_SST_C = 1.7*D_SST_4.iloc[:,1] - 0.1*D_SST_12.iloc[:,1]
D_SST_C.to_csv(path_results+'/indices_CMIP5/C_full_index_'+season+'.csv',float_format='%g')
#E index
D_SST_E = D_SST_12.iloc[:,1] - 0.5*D_SST_4.iloc[:,1] 
D_SST_E.to_csv(path_results+'/indices_CMIP5/E_full_index_'+season+'.csv',float_format='%g')
