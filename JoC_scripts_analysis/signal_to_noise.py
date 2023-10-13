import numpy as np
import pandas as pd
import xarray as xr
import os, fnmatch
import glob
import utilities.open_data as od
import utilities.regression as reg_am
import utilities.csv2nc
import utilities.plot_sensitivity_maps
import utilities.plot_elipse as plt_el
import utilities.plot_large_pp_map
import regionmask

#Open data--------------------------------------------
ruta = '/datos/julia.mindlin/CMIP6_ensambles/preprocesados' #Dropbox/DATOS_CMIP6' #/historical/mon/tas/past'
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

def cargo_todo_pi(models,ruta,var):
    os.chdir(ruta)
    os.getcwd()
    dic = {}
    dic['piControl'] = {}
    for scenario in dic.keys():
        listOfFiles = os.listdir(ruta+'/piControl/'+var)
        for model in models:
            dic['piControl'][model] = []
            pattern1 = "*"+model+"*piControl*"+"*T42*"
            for entry in listOfFiles:
                if fnmatch.fnmatch(entry,pattern1):
                    dato = xr.open_dataset(ruta+'/'+scenario+'/'+var+'/'+entry)
                    dic['piControl'][model].append(dato)
    return dic

def cross_year_season(month,season):
    #Season is a list with two values, begining and endig season
    return (month >= season[0]) | (month <= season[1])

#Create dictionary
var = 'mon/pr'
dato = od.cargo_todo(scenarios,models,'/pikachu'+ruta,var)
dato_picontrol = cargo_todo_pi(models,'/pikachu'+ruta,var)

#Open GW index (I have this index everywhere)
path_indices = '/home/julia.mindlin/Tesis/BSC/indices'
#Open indices for DJF season
season_name = 'DJF'
GW  = pd.read_csv(path_indices+'/GW_index_'+season_name+'.csv')
gw_DJF = GW.iloc[:,2].values

#Create regression class
reg = reg_am.across_models()
var = 'pr'
#Generate changes
reg.regression_data(dato,scenarios,models,gw_DJF,var)
pr_change_dic = reg.psl_change 
beta = {}
for i in range(len(models)):
    beta[models[i]] = pr_change_dic[i]*gw_DJF[i]

del pr_change_dic
#Evaluate standard deviation from piControl experiments in 1 year periods and standard deviation of no overlapping 30 year periods

# 1 year std for summer precip (interannual internal variability)
climate_mean_30year = {}
variability_mean_1year = {}
for model in models:
    pr = dato_picontrol['piControl'][model][0]
    pr_DJF = pr.sel(time=cross_year_season(pr['time.month'],[12,2]))
    #Saco los primeros dos meses (enero y febrero) y el ultimo mes (diciembre)
    #Calculo la media estacional de tres meses centrada 
    seasonal_means = pr_DJF.pr[2:-1,:,:].rolling(time=3,center=True).mean().dropna('time') 
    #Me quedo con un valor cada tres, que corresponde a las medias centradas en J
    seasonal_means_DJF = seasonal_means[::3,:,:]
    #Climate mean standard deviation
    climate_mean_30year[model] = seasonal_means_DJF.rolling(time=30,center=True).mean().dropna('time')[::30].std(dim='time')
    #Interannual variability mean value for 30 year values	
    standard_deviation_DJF = seasonal_means_DJF.rolling(time=30,center=True).std()
    non_overlap_std_DJF = standard_deviation_DJF.dropna('time')[::30]
    variability_mean_1year[model] = non_overlap_std_DJF.mean(dim='time')

#Evaluate signal to noise in each model
gamma = {}
for model in range(len(models)):
    gamma[model] = beta[model]/variability_mean_1year[model]
    
for i in range(len
for model in models:
    f[model] = climate_mean_30year[model]/variability_mean_1year[model]
               
               
import csv2nc
def create_x(i,j,pr):
        x = np.array([])
        for k in range(len(pr)):
            aux = pr[models[k]]
            x = np.append(x,aux[i-1,j-1].values)
            x_mean = np.mean(x); x_anom = (x - x_mean)
            x_new = x
        return x_new

#Evaluate median values for variability in 30-year and 1-year (this gives the correction f)
Aij = pd.DataFrame(columns=['a','lat','lon'])
lat = climate_mean_30year[models[0]].lat; lon = climate_mean_30year[models[0]].lon
for i in range(len(lat)):
    for j in range(len(lon)):
        x = create_x(i,j,climate_mean_30year)
        if np.any(np.isnan(x)):
            a = pd.DataFrame({'a':[np.nan],'lat':[lat[i-1].values.tolist()],'lon':[lon[j-1].values.tolist()]})
            Aij = pd.concat([Aij,a],axis=0)
            x = np.array([])
        else:
            a = pd.DataFrame({'a':[np.median(x)],'lat':[lat[i-1].values.tolist()],'lon':[lon[j-1].values.tolist()]})
            Aij = pd.concat([Aij,a],axis=0)
        
A = {'coef':Aij.iloc[:,0],'lat':Aij.iloc[:,1],'lon':Aij.iloc[:,2]}
Aaij = pd.DataFrame(A).fillna(0)
Aaij.to_csv('/home/julia.mindlin/Tesis/Capitulo3/scripts/CMIP6_storylines/DJF/signal_noise_data/MEM_median_30year_variance.csv', float_format='%g')
               
Aij = pd.DataFrame(columns=['a','lat','lon'])
lat = variability_mean_1year[models[0]].lat; lon = variability_mean_1year[models[0]].lon
for i in range(len(lat)):
    for j in range(len(lon)):
        x = create_x(i,j,variability_mean_1year)
        if np.any(np.isnan(x)):
            a = pd.DataFrame({'a':[np.nan],'lat':[lat[i-1].values.tolist()],'lon':[lon[j-1].values.tolist()]})
            Aij = pd.concat([Aij,a],axis=0)
            x = np.array([])
        else:
            a = pd.DataFrame({'a':[np.median(x)],'lat':[lat[i-1].values.tolist()],'lon':[lon[j-1].values.tolist()]})
            Aij = pd.concat([Aij,a],axis=0)
        
A = {'coef':Aij.iloc[:,0],'lat':Aij.iloc[:,1],'lon':Aij.iloc[:,2]}
Aaij = pd.DataFrame(A).fillna(0)
Aaij.to_csv('/home/julia.mindlin/Tesis/Capitulo3/scripts/CMIP6_storylines/DJF/signal_noise_data/MEM_median_1year_variance.csv', float_format='%g')

csv2nc.csv_to_nc('/home/julia.mindlin/Tesis/Capitulo3/scripts/CMIP6_storylines/DJF/signal_noise_data')
MEM_median_1year = xr.open_dataset('/home/julia.mindlin/Tesis/Capitulo3/scripts/CMIP6_storylines/DJF/signal_noise_data/MEM_median_1year_variance.nc')
MEM_median_30year = xr.open_dataset('/home/julia.mindlin/Tesis/Capitulo3/scripts/CMIP6_storylines/DJF/signal_noise_data/MEM_median_30year_variance.nc')
f = MEM_median_30year**2/MEM_median_1year**2
               

#finally, we have a map for gamma_forced
suma = beta[models[0]]**2 / variability_mean_1year[models[0]]**2
for model in models:
    suma = suma.fillna(0) + (beta[model]**2 / variability_mean_1year[model]**2)
    
suma = suma / len(models)
gamma_forced = f.copy()
gamma_forced.coef.values = np.sqrt(suma.values - 2*f.coef.values)
               
#Evaluo donde mas del 90% de los modelos tiene el mismo signo, positivo o negativo
positives = np.sign(beta[models[0]]).where(np.sign(beta[models[0]])>0).fillna(0)
negatives = np.sign(beta[models[0]]).where(np.sign(beta[models[0]])<0).fillna(0)
for model in models[1:]:
    positives = positives + np.sign(beta[model]).where(np.sign(beta[model])>0).fillna(0)
    negatives = negatives + np.sign(beta[model]).where(np.sign(beta[model])<0).fillna(0)
    
#Evaluo donde hay entre 1 y 5 negativos
negatives_123456 = negatives.where(negatives == -5).fillna(0)/-5 + negatives.where(negatives == -4).fillna(0)/-4 + negatives.where(negatives == -3).fillna(0)/-3 + negatives.where(negatives == -2).fillna(0)/-2 + negatives.where(negatives == -6).fillna(0)/-6 + negatives.where(negatives == -1).fillna(0)/-1
               
               
#Guardo todo esto 
path = '/home/julia.mindlin/Tesis/Capitulo3/scripts/CMIP6_storylines/DJF/signal_noise_data'
gamma_forced.to_netcdf(path+'/gamma_forced_CMIP6.nc')
positives.to_netcdf(path+'/number_of_models_positive_trend_CMIP6.nc')
negatives.to_netcdf(path+'/number_of_models_negative_trend_CMIP6.nc')              
               
