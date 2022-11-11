#Indices ------------------------------------------------------------------------------------------

#This code evaluates all the remote driver indices needed for a storyline analysis of the SH atmospheric ciruclation in CMIP6. This is: 

#Global Warming - Global mean warming from surface temperature (tas) area weighted 
#Tropical Warming - Tropical upper tropospheric warming zonally averaged and averaged between 15N and 15S. 
#Vortex breakdown delay with respect to 1940 - 1969 evaluated with a linear regression including EESCs
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

models = ['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CAMS-CSM1-0','CanESM5','CESM2_','CESM2-WACCM','CMCC-CM2-SR5','CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','FGOALS-g3','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IITM-ESM','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC6','MIROC-ES2L','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NESM3','NorESM2-LM','NorESM2-MM','TaiESM1','UKESM1-0-LL']

experiments = ['historical','ssp585']

path_data = '/datos/julia.mindlin/CMIP6_ensambles/preprocesados' 
path_results = '/home/julia.mindlin/Tesis/JoC_paper/JoC_results'
os.chdir(path_data)
os.getcwd()
os.makedirs('indices_CMIP6',exist_ok=True)


#Asymmetric Sea Surface temperature change
#SST changes
var = 'mon/tos'
variables = ['tos']
dato = cargo_todo(ruta,experiments,models,var)
seasons = ['DJF'] 


#Eastern STD 
box = [0,-10,260,290]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    #components, 0: full, 1: symmetric, 2: asymmetric
    D_SST_E_std = im.sst_index_asym(models,box,ssts)
    D_SST_E_std = {season: D_SST_E_std}
    D_SST_E_std = pd.DataFrame(D_SST_E_std)
    D_SST_E_std.insert(0,"Modelo", models,True)
    D_SST_E_std.to_csv(path_results+'/indices_CMIP6/E_std_asym_index_'+season+'.csv',float_format='%g')

#Central STD
box = [5,-5,180,250]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    D_SST_C_std = im.sst_index_asym(models,box,ssts)
    D_SST_C_std = {season: D_SST_C_std}
    D_SST_C_std = pd.DataFrame(D_SST_C_std)
    D_SST_C_std.insert(0,"Modelo", models,True)
    D_SST_C_std.to_csv(path_results+'/indices_CMIP6/C_std_asym_index_'+season+'.csv',float_format='%g')


#Niño1.2
box = [0,-10,270,280]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    #components, 0: full, 1: symmetric, 2: asymmetric
    D_SST_12 = im.sst_index_asym(models,box,ssts)
    D_SST_12 = {season: D_SST_12}
    D_SST_12 = pd.DataFrame(D_SST_12)
    D_SST_12.insert(0,"Modelo", models,True)
    D_SST_12.to_csv(path_results+'/indices_CMIP6/Nino12_asym_index_'+season+'.csv',float_format='%g')

#Niño4
box = [5,-5,160,210]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    D_SST_4 = im.sst_index_asym(models,box,ssts)
    D_SST_4 = {season: D_SST_4}
    D_SST_4 = pd.DataFrame(D_SST_4)
    D_SST_4.insert(0,"Modelo", models,True)
    D_SST_4.to_csv(path_results+'/indices_CMIP6/Nino4_asym_index_'+season+'.csv',float_format='%g')

#C index
D_SST_C = 1.7*D_SST_4.iloc[:,1] - 0.1*D_SST_12.iloc[:,1]
D_SST_C.to_csv(path_results+'/indices_CMIP6'/C_asym_index_'+season+'.csv',float_format='%g')
#E index
D_SST_E = D_SST_12.iloc[:,1] - 0.5*D_SST_4.iloc[:,1] 
D_SST_E.to_csv(path_results+'/indices_CMIP6/E_asym_index_'+season+'.csv',float_format='%g')


#Niño1.2
box = [0,-10,270,280]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    #components, 0: full, 1: symmetric, 2: asymmetric
    D_SST_12 = im.sst_index_full(models,box,ssts)
    D_SST_12 = {season: D_SST_12}
    D_SST_12 = pd.DataFrame(D_SST_12)
    D_SST_12.insert(0,"Modelo", models,True)
    D_SST_12.to_csv(path_results+'/indices_CMIP6/Nino12_full_index_'+season+'.csv',float_format='%g')

#Niño4
box = [5,-5,160,210]
for season in seasons:
    ssts = im.changes_list(dato,experiments,models,season)
    D_SST_4 = im.sst_index_full(models,box,ssts)
    D_SST_4 = {season: D_SST_4}
    D_SST_4 = pd.DataFrame(D_SST_4)
    D_SST_4.insert(0,"Modelo", models,True)
    D_SST_4.to_csv(path_results+'/indices_CMIP6/Nino4_full_index_'+season+'.csv',float_format='%g')

#C index
D_SST_C = 1.7*D_SST_4.iloc[:,1] - 0.1*D_SST_12.iloc[:,1]
D_SST_C.to_csv(path_results+'/indices_CMIP6/C_full_index_'+season+'.csv',float_format='%g')
#E index
D_SST_E = D_SST_12.iloc[:,1] - 0.5*D_SST_4.iloc[:,1] 
D_SST_E.to_csv(path_results+'/indices_CMIP6/E_full_index_'+season+'.csv',float_format='%g')
