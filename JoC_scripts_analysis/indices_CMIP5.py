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

#Global warming
var = 'mon/tas'
variables = ['tas']
dato = cargo_todo(path_data,experiments,models,var)
seasons = ['DJF'] 

for season in seasons:
    GW = im.global_warming(dato,season,models,experiments)
    GW = {season: DT}
    GW = pd.DataFrame(GW)
    GW.insert(0,"Modelo", models,True)
    GW.to_csv(path_results+'/indices_CMIP6/GW_index_'+season+'.csv',float_format='%g')
    
#Tropical warming
var = 'mon/ta'
variables = ['ta']
ruta = path
dato = cargo_todo(ruta,scenarios,models,var)
seasons = ['DJF'] #['MAM','JJA','SON']

for season in seasons:
    DT = im.tropical_warming(dato,season,models,experiments) 
    TA = {season: DT}
    TA = pd.DataFrame(TA)
    TA.insert(0,"Modelo", models,True)
    TA.to_csv(path_results+'/indices_CMIP6/TA_index_'+season+'.csv',float_format='%g')

    
#Vortex Breakdown-----------------------------------------------------------------------
### TO BE COPIED FROM PREVIOUS CHAPTER

