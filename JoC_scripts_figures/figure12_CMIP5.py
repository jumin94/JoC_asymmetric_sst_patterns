#Figure 12
#Series per degree of warming in regions A, B, F and two storyline coefficients

#---------------------------------------------------IMPORTS------------------------------------------------
## Imports
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
import matplotlib as mpl
from cartopy.util import add_cyclic_point
import matplotlib.pyplot as plt
from global_land_mask import globe

#---------------------------------------------------FUNCTIONS--------------------------------------------------
def std(index):
    return (index - np.mean(index))/np.std(index)

#---------------------------------------------------CODE------------------------------------------------------------
#Open data
ruta = '/datos/julia.mindlin/CMIP6_ensambles/preprocesados' 
models = [
    'ACCESS-CM2', 'ACCESS-ESM1-5', 'BCC-CSM2-MR', 'CAMS-CSM1-0',
    'CanESM5', 'CESM2_', 'CESM2-WACCM','CMCC-CM2-SR5','CNRM-CM6-1',
    'CNRM-ESM2-1','EC-Earth3', 'FGOALS-g3', 'HadGEM3-GC31-LL','HadGEM3-GC31-MM',
    'IITM-ESM','INM-CM4-8','INM-CM5-0','KACE-1-0-G',
    'MIROC6','MIROC-ES2L', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR',
    'MRI-ESM2-0', 'NESM3', 'NorESM2-LM', 'NorESM2-MM', 'TaiESM1','UKESM1-0-LL'
    ]

scenarios = ['historical','ssp585']
os.chdir(ruta)
os.getcwd()

#Create dictionary
var = 'mon/pr'
dato = od.cargo_todo(scenarios,models,ruta,var)
path_results = '/home/julia.mindlin/Tesis/JoC_paper/JoC_results'
#Open indices for DJF season
season_name = 'DJF'
GW  = pd.read_csv(path_results+'/indices_CMIP6/GW_index_'+season_name+'.csv')
gw_DJF = GW.iloc[:,2].values
E_index  = pd.read_csv(path_results+'/indices_CMIP6/E_std_asym_index_'+season_name+'.csv')
e_index_DJF = E_index.iloc[:,2]
C_index = pd.read_csv(path_results+'/indices_CMIP6/C_std_asym_index_'+season_name+'.csv')
c_index_DJF = C_index.iloc[:,2]
TW = pd.read_csv(path_results+'/indices_CMIP6/TW_index_'+season_name+'.csv')
tw_DJF = TW.iloc[:,2]
VB = pd.read_csv(path_results+'/indices_CMIP6/VB_regresion_coef_all_models.csv')
vb =  VB.iloc[:,2]
tw_DJF = tw_DJF/gw_DJF

#Create regression class
reg = reg_am.across_models()
var = 'pr'
#Generate regression data (precip changes)
reg.regression_data(dato,scenarios,models,gw_DJF,var)
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


#Calculo los cambios promedidados en cada caja
#Boxplot precip changes in each box for each model and each storyline
SESA_changes = {}
SESA_changes['IPCC'] = []
SESA_changes['mindlin'] = []
SESA_changes['diaz'] = []
SESA_changes['junquas'] = []
SESA_changes['zilli'] = []
SESA_changes['gonzalez'] = []
for i in range(len(models)):
    A = mask_IPCC.values*pr_change[i].values*86400
    IPCC_mean = A[~np.isnan(A)].mean() 
    SESA_changes['IPCC'].append(IPCC_mean)
    SESA_changes['mindlin'].append((pr_change[i].sel(lon=slice(292,305)).sel(lat=slice(-27,-42))*86400).mean(dim='lon').mean(dim='lat').values)
    SESA_changes['diaz'].append((pr_change[i].sel(lon=slice(293.75,298.75)).sel(lat=slice(-26.25,-38.75))*86400).mean(dim='lon').mean(dim='lat').values)
    SESA_changes['junquas'].append((pr_change[i].sel(lon=slice(300,310)).sel(lat=slice(-25,-32))*86400).mean(dim='lon').mean(dim='lat').values)
    SESA_changes['zilli'].append((pr_change[i].sel(lon=slice(311,320)).sel(lat=slice(-20,-25))*86400).mean(dim='lon').mean(dim='lat').values)
    SESA_changes['gonzalez'].append((pr_change[i].sel(lon=slice(295,315)).sel(lat=slice(-20,-40))*86400).mean(dim='lon').mean(dim='lat').values)
    

#Load sensitivity maps
GlobalWarming = xr.open_dataset(path_results+'/sensitivity_maps_CMIP6/Aij.nc')*86400
TropicalWarming = xr.open_dataset(path_results+'/sensitivity_maps_CMIP6/TAij.nc')*86400
VorBreak_GW = xr.open_dataset(path_results+'/sensitivity_maps_CMIP6/VBij.nc')*86400
SST_C = xr.open_dataset(path_results+'/sensitivity_maps_CMIP6/SST_1ij.nc')*86400
SST_E = xr.open_dataset(path_results+'/sensitivity_maps_CMIP6/SST_2ij.nc')*86400

GlobalWarmingp = xr.open_dataset(path_maps+'/sensitivity_maps_CMIP6/Apij.nc')
TropicalWarmingp = xr.open_dataset(path_maps+'/sensitivity_maps_CMIP6/TApij.nc')
VorBreak_GWp = xr.open_dataset(path_maps+'/sensitivity_maps_CMIP6/VBpij.nc')
SST_Cp = xr.open_dataset(path_maps+'/sensitivity_maps_CMIP6/SST_1pij.nc')
SST_Ep = xr.open_dataset(path_maps+'/sensitivity_maps_CMIP6/SST_2pij.nc')

#Mask IPCC box
A = mask_IPCC.values*GlobalWarming.coef.values
GW_IPCC_mean = A[~np.isnan(A)].mean() 
A = mask_IPCC.values*TropicalWarming.coef.values
TW_IPCC_mean = A[~np.isnan(A)].mean() 
A = mask_IPCC.values*SST_C.coef.values
SST_C_IPCC_mean = A[~np.isnan(A)].mean() 
A = mask_IPCC.values*SST_E.coef.values
SST_E_IPCC_mean = A[~np.isnan(A)].mean() 
A = mask_IPCC.values*VorBreak_GW.coef.values
VB_IPCC_mean = A[~np.isnan(A)].mean() 

#Observed precipitation to evaluate observed trends
path_GPCP = '/datos'
precip_gpcp = xr.open_dataset(path_GPCP+'/precip.mon.mean.GPCP.nc')

lat = precip_gpcp.lat.values
lon = precip_gpcp.lon.values - 180
lon_grid, lat_grid = np.meshgrid(lon,lat)
globe_land_mask = globe.is_land(lat_grid, lon_grid)

def land_masked_field(dato):
    masked = dato*globe_land_mask
    return masked.where(masked != 0)

precip_DJF_1989_2018 = (precip_gpcp.sel(time=slice('1989','2018')).precip/30).groupby('time.season').mean(dim='time').sel(season='DJF')
precip_DJF_1989_2018.attrs = precip_gpcp.precip.attrs
precip_DJF_1950_1979 = (precip_gpcp.sel(time=slice('1950','1979')).precip/30).groupby('time.season').mean(dim='time').sel(season='DJF')
precip_DJF_1950_1979.attrs = precip_gpcp.precip.attrs


#------------------------------------------------------FIGURE---------------------------------------------------------------------
def figure():
    fig, ax = plt.subplots(3, sharex='col', sharey='row',figsize=(12,15))

    GW = np.array([0,1,2,3,4,5,6])
    ts = 1.22
    ax[0].plot(GW,GW*GW_IPCC_mean*ts,color='red')
    ax[0].set_title('Region A',fontsize=20)
    for m in range(len(models)):
        GW_model = np.linspace(0,gw_DJF[m],10)
        ax[0].plot(GW_model,GW_model*SESA_changes['IPCC'][m],color='grey')
    ax[0].plot(GW_model,GW_model*SESA_changes['IPCC'][m],color='grey',label='CMIP6 models')
    observed_rate = (land_masked_field(precip_DJF_1989_2018).sel(lon=slice(290,326)).sel(lat=slice(-37,-20)).mean(dim='lon').mean(dim='lat') - land_masked_field(precip_DJF_1950_1979).sel(lon=slice(290,326)).sel(lat=slice(-37,-20)).mean(dim='lon').mean(dim='lat')).values/.8
    ax[0].plot(GW,GW*observed_rate,color='k',linestyle='dashed',linewidth=2.5,label='extrapolated observed trend (GPCP)')
    ax[0].hlines(0,0,7,color='k',linewidth=0.5,linestyle='dashed')
    ax[0].plot(GW,GW*GW_IPCC_mean,color='k',linewidth=2.5,label='CMIP6 ensemble mean')
    ax[0].plot(GW,GW*(-TW_IPCC_mean*ts + GW_IPCC_mean),color='red',linewidth=1.5,label='wet storyline ts = 1.22')
    ax[0].plot(GW,GW*((VB_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='red',linewidth=1.5)
    ax[0].plot(GW,GW*((VB_IPCC_mean + SST_C_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='red',linewidth=1.5)
    ax[0].plot(GW,GW*((VB_IPCC_mean + SST_C_IPCC_mean + SST_E_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='red',linewidth=1.5)
    ts = -1.22
    ax[0].plot(GW,GW*(-TW_IPCC_mean*ts + GW_IPCC_mean),color='blue',linewidth=1.5,label = 'dry storyline ts = -1.22')
    ax[0].plot(GW,GW*((VB_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='blue',linewidth=1.5)
    ax[0].plot(GW,GW*((VB_IPCC_mean + SST_C_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='blue',linewidth=1.5)
    ax[0].plot(GW,GW*((VB_IPCC_mean + SST_C_IPCC_mean + SST_E_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='blue',linewidth=1.5)
    #Scatter
    ts = 1.22
    ax[0].scatter(GW,GW*(-TW_IPCC_mean*ts + GW_IPCC_mean),color='k',marker='x',linewidth=1.5,label = 'TW')
    ax[0].scatter(GW,GW*((VB_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='k',marker='o',linewidth=1.5,label = 'TW,VB')
    ax[0].scatter(GW,GW*((VB_IPCC_mean + SST_C_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='k',marker='v',linewidth=1.5,label = 'TW,VB,CP')
    ax[0].scatter(GW,GW*((VB_IPCC_mean + SST_C_IPCC_mean + SST_E_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='k',marker='D',linewidth=1.5,label = 'TW,VB,CP,EP')
    ax[0].scatter(GW,GW*(-TW_IPCC_mean*ts + GW_IPCC_mean),color='red',marker='x',linewidth=1.5)
    ax[0].scatter(GW,GW*((VB_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='red',marker='o',linewidth=1.5)
    ax[0].scatter(GW,GW*((VB_IPCC_mean + SST_C_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='red',marker='v',linewidth=1.5)
    ax[0].scatter(GW,GW*((VB_IPCC_mean + SST_C_IPCC_mean + SST_E_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='red',marker='D',linewidth=1.5)
    ts = -1.22
    ax[0].scatter(GW,GW*(-TW_IPCC_mean*ts + GW_IPCC_mean),color='blue',marker='x',linewidth=1.5)
    ax[0].scatter(GW,GW*((VB_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='blue',marker='o',linewidth=1.5)
    ax[0].scatter(GW,GW*((VB_IPCC_mean + SST_C_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='blue',marker='v',linewidth=1.5)
    ax[0].scatter(GW,GW*((VB_IPCC_mean + SST_C_IPCC_mean + SST_E_IPCC_mean - TW_IPCC_mean)*ts + GW_IPCC_mean),color='blue',marker='D',linewidth=1.5)
    ax[0].set_ylim(-1,2)
    ax[0].grid(True)
    ax[0].set_ylabel(r'$\Delta$ pr [mm day$^{-1}$]',fontsize=24)
    ax[0].tick_params(axis='both', which='both', labelsize=18,
                   bottom=False, top=False, labelbottom=False,
                   left=False, right=False, labelleft=True)

    box = [300,310,-32,-25]
    for m in range(len(models)):
        GW_model = np.linspace(0,gw_DJF[m],10)
        ax[1].plot(GW_model,GW_model*SESA_changes['junquas'][m],color='grey')
    ax[1].hlines(0,0,7,color='k',linewidth=1,linestyle='dashed')
    observed_rate = (precip_DJF_1989_2018.sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat') - precip_DJF_1950_1979.sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat')).values/.8
    ax[1].plot(GW,GW*observed_rate,color='k',linestyle='dashed',linewidth=2.5)
    ts = 1.22
    sl0 = GlobalWarming.coef
    sl1 = - TropicalWarming.coef
    sl2 = VorBreak_GW.coef - TropicalWarming.coef 
    sl3 = VorBreak_GW.coef + SST_C.coef -TropicalWarming.coef 
    sl4 = VorBreak_GW.coef + SST_C.coef + SST_E.coef - TropicalWarming.coef
    ax[1].plot(GW,GW*sl0.sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='k',linewidth=2.5)
    ax[1].plot(GW,GW*(sl1*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',linewidth=1.5)
    ax[1].plot(GW,GW*((sl2*ts + sl0)).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',linewidth=1.5)
    ax[1].plot(GW,GW*((sl3*ts + sl0) ).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',linewidth=1.5)
    ax[1].plot(GW,GW*((sl4*ts + sl0)).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',linewidth=1.5)
    ts = -1.22
    ax[1].plot(GW,GW*(sl1*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',linewidth=1.5)
    ax[1].plot(GW,GW*(sl2*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',linewidth=1.5)
    ax[1].plot(GW,GW*(sl3*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',linewidth=1.5)
    ax[1].plot(GW,GW*(sl4*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',linewidth=1.5)
    #Scatter
    ts = 1.22
    ax[1].scatter(GW,GW*(sl1*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',marker='x',linewidth=1.5)
    ax[1].scatter(GW,GW*(sl2*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',marker='o',linewidth=1.5)
    ax[1].scatter(GW,GW*(sl3*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',marker='v',linewidth=1.5)
    ax[1].scatter(GW,GW*(sl4*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',marker='D',linewidth=1.5)
    ts = -1.22
    ax[1].scatter(GW,GW*(sl1*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',marker='x',linewidth=1.5)
    ax[1].scatter(GW,GW*(sl2*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',marker='o',linewidth=1.5)
    ax[1].scatter(GW,GW*(sl3*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',marker='v',linewidth=1.5)
    ax[1].scatter(GW,GW*(sl4*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',marker='D',linewidth=1.5)

    ax[1].set_title('Region B',fontsize=20)
    ax[1].grid(True)
    ax[1].tick_params(axis='both', which='both', labelsize=18,
                   bottom=False, top=False, labelbottom=False,
                   left=False, right=False, labelleft=True)
    ax[1].set_ylabel(r'$\Delta$ pr [mm day$^{-1}$]',fontsize=24)

    box = [311,320,-25,-20]
    for m in range(len(models)):
        GW_model = np.linspace(0,gw_DJF[m],10)
        plt.plot(GW_model,GW_model*SESA_changes['zilli'][m],color='grey')
    ax[2].hlines(0,0,7,color='k',linewidth=1,linestyle='dashed')
    observed_rate = (precip_DJF_1989_2018.sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat') - precip_DJF_1950_1979.sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat')).values/.8
    ax[2].plot(GW,GW*observed_rate,color='k',linestyle='dashed',linewidth=2.5)
    sl0 = GlobalWarming.coef
    sl1 = - TropicalWarming.coef
    sl2 = VorBreak_GW.coef - TropicalWarming.coef 
    sl3 = VorBreak_GW.coef - SST_C.coef -TropicalWarming.coef 
    sl4 = VorBreak_GW.coef - SST_C.coef + SST_E.coef - TropicalWarming.coef
    ax[2].plot(GW,GW*sl0.sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='k',linewidth=2.5)
    ts = 1.22
    ax[2].plot(GW,GW*(sl1*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',linewidth=1.5)
    ax[2].plot(GW,GW*(sl2*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',linewidth=1.5)
    ax[2].plot(GW,GW*(sl3*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',linewidth=1.5)
    ax[2].plot(GW,GW*(sl4*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',linewidth=1.5)
    ax[2].scatter(GW,GW*(sl1*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',marker='x',linewidth=1.5)
    ax[2].scatter(GW,GW*(sl2*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',marker='o',linewidth=1.5)
    ax[2].scatter(GW,GW*(sl3*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',marker='v',linewidth=1.5)
    ax[2].scatter(GW,GW*(sl4*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='red',marker='D',linewidth=1.5)
    ts = -1.22
    ax[2].plot(GW,GW*(sl1*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',linewidth=1.5)
    ax[2].plot(GW,GW*(sl2*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',linewidth=1.5)
    ax[2].plot(GW,GW*(sl3*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',linewidth=1.5)
    ax[2].plot(GW,GW*(sl4*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',linewidth=1.5)
    ax[2].scatter(GW,GW*(sl1*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',marker='x',linewidth=1.5)
    ax[2].scatter(GW,GW*(sl2*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',marker='o',linewidth=1.5)
    ax[2].scatter(GW,GW*(sl3*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',marker='v',linewidth=1.5)
    ax[2].scatter(GW,GW*(sl4*ts + sl0).sel(lon=slice(box[0],box[1])).sel(lat=slice(box[2],box[3])).mean(dim='lon').mean(dim='lat').values,color='blue',marker='D',linewidth=1.5)
    ax[2].set_title('Region F',fontsize=20)
    ax[2].grid(True, 'major', 'both', ls='--', lw=.5, c='k', alpha=.7)


    ax[2].tick_params(axis='both', which='both', labelsize=18,
                   bottom=False, top=False, labelbottom=True,
                   left=False, right=False, labelleft=True)
    ax[2].set_ylabel(r'$\Delta$ pr [mm day$^{-1}$]',fontsize=24)
    ax[2].set_xlabel(r'$\Delta$ T [K]',fontsize=24)
    ax[2].set_xlim(0,7)
    ax[0].legend(loc='upper left',ncol=4)
    plt.show()
    return fig

path_fig = '/home/julia.mindlin/Tesis/JoC_paper/JoC_scripts_figures'
figure12 = figure()
figure12.save_fig(path_fig+'/figure12_storylines_per_degree_warming_PR.nc')

