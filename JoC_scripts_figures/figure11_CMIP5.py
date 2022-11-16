#Figure 11
#Scatter plots

## ---------------------------------------------------IMPORTS-----------------------------------------------------
import numpy as np
import pandas as pd
import xarray as xr
import os, fnmatch
import glob
import utilities.open_data as od
import utilities.regression as reg_am
import utilities.csv2nc
import utilities.plot_sensitivity_maps
import regionmask
import matplotlib.pyplot as plt 

#---------------------------------------------------FUNCIONES----------------------------------------------------
def std(dato):
    return (dato - np.mean(dato))/np.std(dato)

#---------------------------------------------------OPEN DATA----------------------------------------------------
models = ['ACCESS1-0','ACCESS1-3','bcc-csm1-1','bcc-csm1-1-m','CanESM2','CCSM4','CMCC-CESM','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','HadGEM2-CC','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M']

experiments = ['historical','rcp85']

path_data = '/storage/silver/acrcc/co157815/fromjasmin/cmip5'
path_results = '/home/users/co157815/JoC_paper/JoC_results'
os.chdir(path_data)
os.getcwd()

#Create dictionary
gloW  = pd.read_csv(path_results+'/indices_CMIP5/GW_VB_index50_70.csv')
gw_DJF = gloW.iloc[:,4].values
TA = pd.read_csv(path_results+'/indices_CMIP5/TW_index.csv')
VB = pd.read_csv(path_results+'/indices_CMIP5/VB_indexGWscaled50_70.csv')
SST_C_std = pd.read_csv(path_results+'/indices_CMIP5/C_std_asym_index_DJF.csv')
SST_E_std = pd.read_csv(path_results+'/indices_CMIP5/E_std_asym_index_DJF.csv')
#Select values
tw_DJF = TA.iloc[:,4].values
vb = VB.iloc[:,5].values
c_index_DJF = SST_C_std.iloc[:,2].values
e_index_DJF = SST_E_std.iloc[:,2].values

#---------------------------------------------------------------CODE---------------------------------------------
#Create regression class
reg = reg_am.across_models()
var = 'pr'
#Generate regression data
reg.regression_data(dato,scenarios,models,gw_index,var)
pr_change = reg.psl_change

#Define masks
GW = pr_change[0]
#Box from paper 2020
mask_mindlin = GW.where(GW.lon<305).where(GW.lon>292).where(GW.lat>-42).where(GW.lat<-27) / GW.where(GW.lon<305).where(GW.lon>292).where(GW.lat>-42).where(GW.lat<-27)
#Box from IPCC
import regionmask
mask_IPCC = regionmask.defined_regions.ar6.land.mask(GW).where(regionmask.defined_regions.ar6.land.mask(GW) == 14)/14
#Box from paper Lean
#38.75∘S–26.25∘S, 66.25–61.25∘W
mask_diaz = GW.where(GW.lon<298.75).where(GW.lon>293.75).where(GW.lat>-38.75).where(GW.lat<-26.25) / GW.where(GW.lon<298.75).where(GW.lon>293.75).where(GW.lat>-38.75).where(GW.lat<-26.25)
#Box from Clementine
#32°S:25°S,60°W:50°W
mask_junquas = GW.where(GW.lon<310).where(GW.lon>300).where(GW.lat>-32).where(GW.lat<-25) / GW.where(GW.lon<310).where(GW.lon>300).where(GW.lat>-32).where(GW.lat<-25)
#Box from Clementine
#32°S:25°S,60°W:50°W  
mask_zilli = GW.where(GW.lon<320).where(GW.lon>311).where(GW.lat>-25).where(GW.lat<-20) / GW.where(GW.lon<310).where(GW.lon>300).where(GW.lat>-32).where(GW.lat<-25)
#Box from Paula
#32°S:25°S,60°W:50°W  
mask_gonzalez = GW.where(GW.lon<315).where(GW.lon>295).where(GW.lat>-40).where(GW.lat<-20) / GW.where(GW.lon<295).where(GW.lon>315).where(GW.lat>-40).where(GW.lat<-20)

#Boxplot precip changes in each box for each model
SESA_changes = {}
SESA_changes['IPCC'] = []
SESA_changes['mindlin'] = []
SESA_changes['diaz'] = []
SESA_changes['junquas'] = []
SESA_changes['zilli'] = []
SESA_changes['gonzalez'] = []
for i in range(len(models)):
    A = mask_IPCC.values*pr_change[i].values*86400*gw_DJF[i]
    IPCC_mean = A[~np.isnan(A)].mean() 
    SESA_changes['IPCC'].append(IPCC_mean)
    SESA_changes['mindlin'].append((pr_change[i].sel(lon=slice(292,305)).sel(lat=slice(-27,-42))*86400*gw_DJF[i]).mean(dim='lon').mean(dim='lat').values)
    SESA_changes['diaz'].append((pr_change[i].sel(lon=slice(293.75,298.75)).sel(lat=slice(-26.25,-38.75))*86400*gw_DJF[i]).mean(dim='lon').mean(dim='lat').values)
    SESA_changes['junquas'].append((pr_change[i].sel(lon=slice(300,310)).sel(lat=slice(-25,-32))*86400*gw_DJF[i]).mean(dim='lon').mean(dim='lat').values)
    SESA_changes['zilli'].append((pr_change[i].sel(lon=slice(311,320)).sel(lat=slice(-20,-25))*86400*gw_DJF[i]).mean(dim='lon').mean(dim='lat').values)
    SESA_changes['gonzalez'].append((pr_change[i].sel(lon=slice(295,315)).sel(lat=slice(-20,-40))*86400*gw_DJF[i]).mean(dim='lon').mean(dim='lat').values)
    
#Open sensitiviy maps
GlobalWarming = xr.open_dataset(path_results+'/sensitivity_maps_CMIP5/Aij.nc')*86400
TropicalWarming = xr.open_dataset(path_results+'/sensitivity_maps_CMIP5/TAij.nc')*86400
VorBreak_GW = xr.open_dataset(path_results+'/sensitivity_maps_CMIP5/VBij.nc')*86400
SST_C = xr.open_dataset(path_results+'/sensitivity_maps_CMIP5/SST_1ij.nc')*86400
SST_E = xr.open_dataset(path_results+'/sensitivity_maps_CMIP5/SST_2ij.nc')*86400

GlobalWarmingp = xr.open_dataset(path_maps+'/sensitivity_maps_CMIP5/Apij.nc')
TropicalWarmingp = xr.open_dataset(path_maps+'/sensitivity_maps_CMIP5/TApij.nc')
VorBreak_GWp = xr.open_dataset(path_maps+'/sensitivity_maps_CMIP5/VBpij.nc')
SST_Cp = xr.open_dataset(path_maps+'/sensitivity_maps_CMIP5/SST_1pij.nc')
SST_Ep = xr.open_dataset(path_maps+'/sensitivity_maps_CMIP5/SST_2pij.nc')

#Storyline construction
storyline_construction = []
ts = np.arange(-1.5,1.5,.1)
for k in range(len(ts)):
   storyline_construction.append(GlobalWarming.coef + (SST_C.coef)*ts[k]) #+ TropicalWarming.coef*ts[k] + VorBreak_GW.coef*ts[k]

SESA_storylines = {}
SESA_storylines['IPCC'] = []
SESA_storylines['mindlin'] = []
SESA_storylines['diaz'] = []
SESA_storylines['junquas'] = []
SESA_storylines['gonzalez'] = []
SESA_storylines['zilli'] = []
for i in range(len(ts)):
    #gw_DJF_modeled = np.random.uniform(low=np.min(gw_DJF)+.5, high=np.max(gw_DJF)-.5, size=28)
    A = mask_IPCC.values*storyline_construction[i].values #*gw_DJF_modeled[i]
    IPCC_mean = A[~np.isnan(A)].mean() 
    SESA_storylines['IPCC'].append(IPCC_mean)
    SESA_storylines['mindlin'].append((storyline_construction[i].sel(lon=slice(292,305)).sel(lat=slice(-42,-27))).mean(dim='lon').mean(dim='lat').values) #*gw_DJF_modeled[i]
    SESA_storylines['diaz'].append((storyline_construction[i].sel(lon=slice(293.75,298.75)).sel(lat=slice(-38.75,-26.25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines['junquas'].append((storyline_construction[i].sel(lon=slice(300,310)).sel(lat=slice(-32,-25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines['zilli'].append((storyline_construction[i].sel(lon=slice(311,320)).sel(lat=slice(-25,-20))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines['gonzalez'].append((storyline_construction[i].sel(lon=slice(295,315)).sel(lat=slice(-40,-20))).mean(dim='lon').mean(dim='lat').values)
    
storyline_reconstruction = []
for k in range(len(models)):
   storyline_reconstruction.append(GlobalWarming.coef + (SST_E.coef)*std(e_index_DJF)[k] + (SST_C.coef)*std(c_index_DJF)[k] + (TropicalWarming.coef)*std(tw_DJF)[k] + VorBreak_GW.coef*std(vb)[k]) #+ TropicalWarming.coef*ts[k] + VorBreak_GW.coef*ts[k]

SESA_storylines_recon = {}
SESA_storylines_recon['IPCC'] = []
SESA_storylines_recon['mindlin'] = []
SESA_storylines_recon['diaz'] = []
SESA_storylines_recon['junquas'] = []
SESA_storylines_recon['gonzalez'] = []
SESA_storylines_recon['zilli'] = []
for i in range(len(models)):
    #gw_DJF_modeled = np.random.uniform(low=np.min(gw_DJF)+.5, high=np.max(gw_DJF)-.5, size=28)
    A = mask_IPCC.values*storyline_reconstruction[i].values #*gw_DJF_modeled[i]
    IPCC_mean = A[~np.isnan(A)].mean() 
    SESA_storylines_recon['IPCC'].append(IPCC_mean)
    SESA_storylines_recon['mindlin'].append((storyline_reconstruction[i].sel(lon=slice(292,305)).sel(lat=slice(-42,-27))).mean(dim='lon').mean(dim='lat').values) #*gw_DJF_modeled[i]
    SESA_storylines_recon['diaz'].append((storyline_reconstruction[i].sel(lon=slice(293.75,298.75)).sel(lat=slice(-38.75,-26.25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_recon['junquas'].append((storyline_reconstruction[i].sel(lon=slice(300,310)).sel(lat=slice(-32,-25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_recon['zilli'].append((storyline_reconstruction[i].sel(lon=slice(311,320)).sel(lat=slice(-25,-20))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_recon['gonzalez'].append((storyline_reconstruction[i].sel(lon=slice(295,315)).sel(lat=slice(-40,-20))).mean(dim='lon').mean(dim='lat').values)
    
storyline_construction = []
ts = np.arange(-1.5,1.5,.1)
for k in range(len(ts)):
   storyline_construction.append(GlobalWarming.coef + (SST_E.coef)*ts[k]) #+ TropicalWarming.coef*ts[k] + VorBreak_GW.coef*ts[k]

SESA_storylines = {}
SESA_storylines['IPCC'] = []
SESA_storylines['mindlin'] = []
SESA_storylines['diaz'] = []
SESA_storylines['junquas'] = []
SESA_storylines['gonzalez'] = []
SESA_storylines['zilli'] = []
for i in range(len(ts)):
    #gw_DJF_modeled = np.random.uniform(low=np.min(gw_DJF)+.5, high=np.max(gw_DJF)-.5, size=28)
    A = mask_IPCC.values*storyline_construction[i].values #*gw_DJF_modeled[i]
    IPCC_mean = A[~np.isnan(A)].mean() 
    SESA_storylines['IPCC'].append(IPCC_mean)
    SESA_storylines['mindlin'].append((storyline_construction[i].sel(lon=slice(292,305)).sel(lat=slice(-42,-27))).mean(dim='lon').mean(dim='lat').values) #*gw_DJF_modeled[i]
    SESA_storylines['diaz'].append((storyline_construction[i].sel(lon=slice(293.75,298.75)).sel(lat=slice(-38.75,-26.25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines['junquas'].append((storyline_construction[i].sel(lon=slice(300,310)).sel(lat=slice(-32,-25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines['zilli'].append((storyline_construction[i].sel(lon=slice(311,320)).sel(lat=slice(-25,-20))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines['gonzalez'].append((storyline_construction[i].sel(lon=slice(295,315)).sel(lat=slice(-40,-20))).mean(dim='lon').mean(dim='lat').values)
    
    
storyline_reconstruction = []
for k in range(len(models)):
   storyline_reconstruction.append(GlobalWarming.coef + (SST_E.coef)*std(e_index_DJF)[k]  + (TropicalWarming.coef)*std(tw_DJF)[k] + VorBreak_GW.coef*std(vb)[k]) #+ TropicalWarming.coef*ts[k] + VorBreak_GW.coef*ts[k]

SESA_storylines_recon = {}
SESA_storylines_recon['IPCC'] = []
SESA_storylines_recon['mindlin'] = []
SESA_storylines_recon['diaz'] = []
SESA_storylines_recon['junquas'] = []
SESA_storylines_recon['gonzalez'] = []
SESA_storylines_recon['zilli'] = []
for i in range(len(models)):
    #gw_DJF_modeled = np.random.uniform(low=np.min(gw_DJF)+.5, high=np.max(gw_DJF)-.5, size=28)
    A = mask_IPCC.values*storyline_reconstruction[i].values #*gw_DJF_modeled[i]
    IPCC_mean = A[~np.isnan(A)].mean() 
    SESA_storylines_recon['IPCC'].append(IPCC_mean)
    SESA_storylines_recon['mindlin'].append((storyline_reconstruction[i].sel(lon=slice(292,305)).sel(lat=slice(-42,-27))).mean(dim='lon').mean(dim='lat').values) #*gw_DJF_modeled[i]
    SESA_storylines_recon['diaz'].append((storyline_reconstruction[i].sel(lon=slice(293.75,298.75)).sel(lat=slice(-38.75,-26.25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_recon['junquas'].append((storyline_reconstruction[i].sel(lon=slice(300,310)).sel(lat=slice(-32,-25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_recon['zilli'].append((storyline_reconstruction[i].sel(lon=slice(311,320)).sel(lat=slice(-25,-20))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_recon['gonzalez'].append((storyline_reconstruction[i].sel(lon=slice(295,315)).sel(lat=slice(-40,-20))).mean(dim='lon').mean(dim='lat').values)
    
    
storyline_construction = []
ts = np.arange(-1.5,1.5,.1)
for k in range(len(ts)):
   storyline_construction.append(GlobalWarming.coef + (SST_C.coef)*ts[k]) #+ TropicalWarming.coef*ts[k] + VorBreak_GW.coef*ts[k]

SESA_storylines = {}
SESA_storylines['IPCC'] = []
SESA_storylines['mindlin'] = []
SESA_storylines['diaz'] = []
SESA_storylines['junquas'] = []
SESA_storylines['gonzalez'] = []
SESA_storylines['zilli'] = []
for i in range(len(ts)):
    #gw_DJF_modeled = np.random.uniform(low=np.min(gw_DJF)+.5, high=np.max(gw_DJF)-.5, size=28)
    A = mask_IPCC.values*storyline_construction[i].values #*gw_DJF_modeled[i]
    IPCC_mean = A[~np.isnan(A)].mean() 
    SESA_storylines['IPCC'].append(IPCC_mean)
    SESA_storylines['mindlin'].append((storyline_construction[i].sel(lon=slice(292,305)).sel(lat=slice(-42,-27))).mean(dim='lon').mean(dim='lat').values) #*gw_DJF_modeled[i]
    SESA_storylines['diaz'].append((storyline_construction[i].sel(lon=slice(293.75,298.75)).sel(lat=slice(-38.75,-26.25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines['junquas'].append((storyline_construction[i].sel(lon=slice(300,310)).sel(lat=slice(-32,-25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines['zilli'].append((storyline_construction[i].sel(lon=slice(311,320)).sel(lat=slice(-25,-20))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines['gonzalez'].append((storyline_construction[i].sel(lon=slice(295,315)).sel(lat=slice(-40,-20))).mean(dim='lon').mean(dim='lat').values)
    
storyline_construction_E = []
ts = np.arange(-1.5,1.5,.1)
for k in range(len(ts)):
   storyline_construction_E.append(GlobalWarming.coef + (SST_E.coef)*ts[k]) #+ TropicalWarming.coef*ts[k] + VorBreak_GW.coef*ts[k]

SESA_storylines_E = {}
SESA_storylines_E['IPCC'] = []
SESA_storylines_E['mindlin'] = []
SESA_storylines_E['diaz'] = []
SESA_storylines_E['junquas'] = []
SESA_storylines_E['gonzalez'] = []
SESA_storylines_E['zilli'] = []
for i in range(len(ts)):
    #gw_DJF_modeled = np.random.uniform(low=np.min(gw_DJF)+.5, high=np.max(gw_DJF)-.5, size=28)
    A = mask_IPCC.values*storyline_construction_E[i].values #*gw_DJF_modeled[i]
    IPCC_mean = A[~np.isnan(A)].mean() 
    SESA_storylines_E['IPCC'].append(IPCC_mean)
    SESA_storylines_E['mindlin'].append((storyline_construction_E[i].sel(lon=slice(292,305)).sel(lat=slice(-42,-27))).mean(dim='lon').mean(dim='lat').values) #*gw_DJF_modeled[i]
    SESA_storylines_E['diaz'].append((storyline_construction_E[i].sel(lon=slice(293.75,298.75)).sel(lat=slice(-38.75,-26.25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_E['junquas'].append((storyline_construction_E[i].sel(lon=slice(300,310)).sel(lat=slice(-32,-25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_E['zilli'].append((storyline_construction_E[i].sel(lon=slice(311,320)).sel(lat=slice(-25,-20))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_E['gonzalez'].append((storyline_construction_E[i].sel(lon=slice(295,315)).sel(lat=slice(-40,-20))).mean(dim='lon').mean(dim='lat').values)
    
storyline_reconstruction = []
for k in range(len(models)):
   storyline_reconstruction.append(GlobalWarming.coef + (SST_C.coef)*std(c_index_DJF)[k] + (SST_E.coef)*std(e_index_DJF)[k]  - (TropicalWarming.coef)*std(tw_DJF)[k] + VorBreak_GW.coef*std(vb)[k]) #+ TropicalWarming.coef*ts[k] + VorBreak_GW.coef*ts[k]

SESA_storylines_recon = {}
SESA_storylines_recon['IPCC'] = []
SESA_storylines_recon['mindlin'] = []
SESA_storylines_recon['diaz'] = []
SESA_storylines_recon['junquas'] = []
SESA_storylines_recon['gonzalez'] = []
SESA_storylines_recon['zilli'] = []
for i in range(len(models)):
    #gw_DJF_modeled = np.random.uniform(low=np.min(gw_DJF)+.5, high=np.max(gw_DJF)-.5, size=28)
    A = mask_IPCC.values*storyline_reconstruction[i].values #*gw_DJF_modeled[i]
    IPCC_mean = A[~np.isnan(A)].mean() 
    SESA_storylines_recon['IPCC'].append(IPCC_mean)
    SESA_storylines_recon['mindlin'].append((storyline_reconstruction[i].sel(lon=slice(292,305)).sel(lat=slice(-42,-27))).mean(dim='lon').mean(dim='lat').values) #*gw_DJF_modeled[i]
    SESA_storylines_recon['diaz'].append((storyline_reconstruction[i].sel(lon=slice(293.75,298.75)).sel(lat=slice(-38.75,-26.25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_recon['junquas'].append((storyline_reconstruction[i].sel(lon=slice(300,310)).sel(lat=slice(-32,-25))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_recon['zilli'].append((storyline_reconstruction[i].sel(lon=slice(311,320)).sel(lat=slice(-25,-20))).mean(dim='lon').mean(dim='lat').values)
    SESA_storylines_recon['gonzalez'].append((storyline_reconstruction[i].sel(lon=slice(295,315)).sel(lat=slice(-40,-20))).mean(dim='lon').mean(dim='lat').values)
    
storyline_reconstruction_brasil_C = []
for k in range(len(models)):
    storyline_reconstruction_brasil_C.append(GlobalWarming.coef + (SST_C.coef)*std(c_index_DJF)[k] - (SST_E.coef)*std(e_index_DJF)[k]  + (TropicalWarming.coef)*std(tw_DJF)[k] - VorBreak_GW.coef*std(vb)[k]) #+ TropicalWarming.coef*ts[k] + VorBreak_GW.coef*ts[k]

SESA_storylines_recon_brasil_C = {}
SESA_storylines_recon_brasil_C['zilli'] = []
for i in range(len(models)):
    SESA_storylines_recon_brasil_C['zilli'].append((storyline_reconstruction_brasil_C[i].sel(lon=slice(311,320)).sel(lat=slice(-25,-20))).mean(dim='lon').mean(dim='lat').values)

    
storyline_reconstruction_brasil_E = []
for k in range(len(models)):
    storyline_reconstruction_brasil_E.append(GlobalWarming.coef - (SST_C.coef)*std(c_index_DJF)[k] + (SST_E.coef)*std(e_index_DJF)[k]  - (TropicalWarming.coef)*std(tw_DJF)[k] + VorBreak_GW.coef*std(vb)[k]) #+ TropicalWarming.coef*ts[k] + VorBreak_GW.coef*ts[k]

SESA_storylines_recon_brasil_E = {}
SESA_storylines_recon_brasil_E['zilli'] = []
for i in range(len(models)):
    SESA_storylines_recon_brasil_E['zilli'].append((storyline_reconstruction_brasil_E[i].sel(lon=slice(311,320)).sel(lat=slice(-25,-20))).mean(dim='lon').mean(dim='lat').values)

    
def figure():
    fig, ax = plt.subplots(3,2, sharex='col', sharey='row',figsize=(20,15))
    ax[0,0].plot(ts,SESA_storylines['IPCC'],color='k')
    ax[0,0].scatter(std(c_index_DJF),SESA_storylines_recon['IPCC'],color='red',marker='D')
    ax[0,0].set_title('Region A',fontsize=20)
    ax[0,0].hlines(0,-1.7,3,color='k',linewidth=1,linestyle='dashed')
    ax[0,0].grid(True)
    ax[0,0].set_ylabel(r'$\Delta$ pr [mm day$^{-1}$]',fontsize=24)
    ax[0,0].tick_params(axis='both', which='both', labelsize=18,
                   bottom=False, top=False, labelbottom=False,
                   left=False, right=False, labelleft=True)
    ax[0,0].set_ylim(-.2,.4)


    ax[1,0].plot(ts,SESA_storylines['junquas'],color='k')
    ax[1,0].scatter(std(c_index_DJF),SESA_storylines_recon['junquas'],color='red',marker='D')
    ax[1,0].hlines(0,-1.7,2.5,color='k',linewidth=1,linestyle='dashed')
    #ax[1].scatter(,color='red',marker='D',linewidth=2.5)
    ax[1,0].set_title('Region B',fontsize=20)
    ax[1,0].grid(True)
    ax[1,0].tick_params(axis='both', which='both', labelsize=18,
                   bottom=False, top=False, labelbottom=False,
                   left=False, right=False, labelleft=True)
    ax[1,0].set_ylabel(r'$\Delta$ pr [mm day$^{-1}$]',fontsize=24)
    ax[1,0].set_ylim(-.2,.4)

    ax[2,0].plot(ts,SESA_storylines['zilli'],color='k')
    ax[2,0].scatter(std(c_index_DJF),SESA_storylines_recon_brasil_C['zilli'],color='red',marker='D')
    ax[2,0].hlines(0,-1.7,2.5,color='k',linewidth=1,linestyle='dashed')
    ax[2,0].set_title('Region F',fontsize=20)
    ax[2,0].grid(True, 'major', 'both', ls='--', lw=.5, c='k', alpha=.7)
    ax[2,0].tick_params(axis='both', which='both', labelsize=18,
                   bottom=False, top=False, labelbottom=True,
                   left=False, right=False, labelleft=True)
    ax[2,0].set_ylabel(r'$\Delta$ pr [mm day$^{-1}$]',fontsize=24)
    ax[2,0].set_xlabel('storyline coefficient (CP\') ',fontsize=24)
    ax[2,0].set_xlim(-1.7,2.5)
    ax[2,0].set_ylim(-.2,.4)
    #ax[0,0].legend(loc='upper left',ncol=4)

    ax[0,1].plot(ts,SESA_storylines_E['IPCC'],color='k')
    ax[0,1].scatter(std(e_index_DJF),SESA_storylines_recon['IPCC'],color='red',marker='D')
    ax[0,1].set_title('Region A',fontsize=20)
    ax[0,1].hlines(0,-1.7,3,color='k',linewidth=1,linestyle='dashed')
    ax[0,1].grid(True)
    ax[0,1].set_ylabel(r'$\Delta$ pr [mm day$^{-1}$]',fontsize=24)
    ax[0,1].tick_params(axis='both', which='both', labelsize=18,
                   bottom=False, top=False, labelbottom=False,
                   left=False, right=False, labelleft=True)
    ax[0,1].set_ylim(-.2,.4)


    ax[1,1].plot(ts,SESA_storylines_E['junquas'],color='k')
    ax[1,1].scatter(std(e_index_DJF),SESA_storylines_recon['junquas'],color='red',marker='D')
    ax[1,1].hlines(0,-1.7,2.5,color='k',linewidth=1,linestyle='dashed')
    #ax[1].scatter(,color='red',marker='D',linewidth=2.5)
    ax[1,1].set_title('Region B',fontsize=20)
    ax[1,1].grid(True)
    ax[1,1].tick_params(axis='both', which='both', labelsize=18,
                   bottom=False, top=False, labelbottom=False,
                   left=False, right=False, labelleft=True)
    ax[1,1].set_ylabel(r'$\Delta$ pr [mm day$^{-1}$]',fontsize=24)
    ax[1,1].set_ylim(-.2,.4)

    ax[2,1].plot(ts,SESA_storylines_E['zilli'],color='k')
    ax[2,1].scatter(std(e_index_DJF),SESA_storylines_recon_brasil_E['zilli'],color='red',marker='D')
    ax[2,1].hlines(0,-1.7,2.5,color='k',linewidth=1,linestyle='dashed')
    ax[2,1].set_title('Region F',fontsize=20)
    ax[2,1].grid(True, 'major', 'both', ls='--', lw=.5, c='k', alpha=.7)
    ax[2,1].tick_params(axis='both', which='both', labelsize=18,
                   bottom=False, top=False, labelbottom=True,
                   left=False, right=False, labelleft=True)
    ax[2,1].set_ylabel(r'$\Delta$ pr [mm day$^{-1}$]',fontsize=24)
    ax[2,1].set_xlabel('storyline coefficient  (EP\') ',fontsize=24)
    ax[2,1].set_xlim(-1.7,2.5)
    ax[2,1].set_ylim(-.2,.4)
    ax[2,1].legend(loc='upper left',ncol=4)
    return fig

path_fig = '/home/users/co157815/JoC_paper/JoC_scripts_figures'
figure11 = figure()
figure11.save_fig(path_fig+'/figure11_scatterplots_storyline_vs_pacific_indices_CMIP5.nc')
