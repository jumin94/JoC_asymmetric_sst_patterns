#Python imports
import numpy as np
import pandas as pd
import xarray as xr
import os, fnmatch
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from numpy import linalg as la
from sklearn.linear_model import LinearRegression
import warnings
warnings.filterwarnings("ignore")
#My imports
import utilities.open_data as od
import utilities.regression as reg_am
import utilities.csv2nc
import utilities.plot_sensitivity_maps

#Functions
def std(x):
    return (x - np.mean(x))/np.std(x)


#Open data and define models
models = ['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CAMS-CSM1-0','CanESM5','CESM2_','CESM2-WACCM','CMCC-CM2-SR5',
          'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','FGOALS-g3','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IITM-ESM','INM-CM4-8',
          'INM-CM5-0','KACE-1-0-G','MIROC6','MIROC-ES2L','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0',
          'NESM3','NorESM2-LM','NorESM2-MM','TaiESM1','UKESM1-0-LL']

#Open indices
path_indices = '/home/julia.mindlin/Tesis/indices'
#Open indices for DJF season
season_name = 'DJF'
GW  = pd.read_csv(path_indices+'/GW_index_'+season_name+'.csv')
gw_DJF = GW.iloc[:,2].values
E_index  = pd.read_csv(path_indices+'/E_std_asym_index_'+season_name+'.csv')
e_index_DJF = E_index.iloc[:,2]
C_index = pd.read_csv(path_indices+'/C_std_asym_index_'+season_name+'.csv')
c_index_DJF = C_index.iloc[:,2]
TW = pd.read_csv(path_indices+'/TW_index_'+season_name+'.csv')
tw_DJF = TW.iloc[:,2]
VB = pd.read_csv(path_indices+'/VB_regresion_coef_all_models.csv')
vb_DJF =  VB.iloc[:,2]
#Open index errors for DJF season
GW  = pd.read_csv(path_indices+'/GW_index_errors_'+season_name+'.csv')
gw_index_errors_DJF = GW.iloc[:,3].values
E_index  = pd.read_csv(path_indices+'/E_std_asym_index_errors_'+season_name+'.csv')
e_index_errors_DJF = E_index.iloc[:,3]
C_index = pd.read_csv(path_indices+'/C_std_asym_index_errors_'+season_name+'.csv')
c_index_errors_DJF = C_index.iloc[:,3]
TW = pd.read_csv(path_indices+'/TW_index_errors_'+season_name+'.csv')
tw_errors_DJF = TW.iloc[:,3]
VB = pd.read_csv(path_indices+'/VB_regresion_coef_all_models.csv')
vb_errors =  VB.iloc[:,4]

#Open indices, evaluate mean and standard deviation
tw_DJF = tw_DJF/gw_DJF #Scale by global warming
std_tw = np.std(tw_DJF.values); mean_tw = np.mean(tw_DJF)
std_sv = np.std(vb_DJF.values); mean_sv = np.mean(vb_DJF.values) #scaled by global warming already
std_cp = np.std(c_index_DJF); mean_cp = np.mean(c_index_DJF)
std_ep = np.std(e_index_DJF); mean_ep = np.mean(e_index_DJF)
#Evaluate correlations between indices (this can be done much faster with other functions but just to be clear)
corr_tw_sv = np.corrcoef(std(tw_DJF),std(vb_DJF))[0,1]
corr_tw_cp = np.corrcoef(std(tw_DJF),std(c_index_DJF))[0,1]
corr_tw_ep = np.corrcoef(std(tw_DJF),std(e_index_DJF))[0,1]
corr_sv_cp = np.corrcoef(std(vb_DJF),std(c_index_DJF))[0,1]
corr_sv_ep = np.corrcoef(std(vb_DJF),std(e_index_DJF))[0,1]
corr_cp_ep = np.corrcoef(std(c_index_DJF),std(e_index_DJF))[0,1]


path_maps = '/home/julia.mindlin/Tesis/sensitivity_maps/DJF/pr'
#Create plots
GlobalWarming = xr.open_dataset(path_maps+'/Aij.nc')*86400
TropicalWarming = xr.open_dataset(path_maps+'/TAij.nc')*86400
VorBreak_GW = xr.open_dataset(path_maps+'/VBij.nc')*86400
SST_C = xr.open_dataset(path_maps+'/SST_1ij.nc')*86400
SST_E = xr.open_dataset(path_maps+'/SST_2ij.nc')*86400

GlobalWarmingp = xr.open_dataset(path_maps+'/Apij.nc')
TropicalWarmingp = xr.open_dataset(path_maps+'/TApij.nc')
VorBreak_GWp = xr.open_dataset(path_maps+'/VBpij.nc')
SST_Cp = xr.open_dataset(path_maps+'/SST_1pij.nc')
SST_Ep = xr.open_dataset(path_maps+'/SST_2pij.nc')

#Generate storylines
import scipy.stats as ss
#Simulate storylines with the same correlations that in the observed indices
# We need to recalculate the covariance matrix
# using the estimated paramaters
cov_matrix = [[std_tw**2, corr_tw_sv*std_tw*std_sv,corr_tw_cp*std_tw*std_cp,corr_tw_ep*std_tw*std_ep], 
              [corr_tw_sv*std_sv*std_tw, std_sv**2,corr_sv_cp*std_sv*std_cp,corr_sv_ep*std_sv*std_ep],
              [corr_tw_cp*std_cp*std_tw, corr_sv_cp*std_sv*std_cp,std_cp**2,corr_cp_ep*std_cp*std_ep],
              [corr_tw_ep*std_ep*std_tw, corr_sv_ep*std_sv*std_ep,corr_cp_ep*std_cp*std_ep,std_ep**2]]

simulated_remote_drivers = ss.multivariate_normal(
    mean=[mean_tw, mean_sv,mean_cp,mean_ep],
    cov=cov_matrix).rvs(60000, random_state=1)

import scipy.stats as ss
#Simulate storylines with the same correlations that in the observed indices
# We need to recalculate the covariance matrix
# using the estimated paramaters
cov_matrix = [[1**2, corr_tw_sv*std_tw*std_sv,corr_tw_cp*std_tw*std_cp,corr_tw_ep*std_tw*std_ep], 
              [corr_tw_sv*std_sv*std_tw, 1**2,corr_sv_cp*std_sv*std_cp,corr_sv_ep*std_sv*std_ep],
              [corr_tw_cp*std_cp*std_tw, corr_sv_cp*std_sv*std_cp,1**2,corr_cp_ep*std_cp*std_ep],
              [corr_tw_ep*std_ep*std_tw, corr_sv_ep*std_sv*std_ep,corr_cp_ep*std_cp*std_ep,1**2]]

#The remote drivers will be simulated from a multivariate normal
simulated_remote_drivers = ss.multivariate_normal(
    mean=[0, 0, 0, 0],
    cov=cov_matrix).rvs(100000, random_state=1)

#The simulated remote driver values come in tandem, each set of four values (one per remote driver) corresponds to one sample
#in this multivariate normal
storyline_dyn = []
ts = simulated_remote_drivers
#The for loop will compute the precipitation value under each storyline
for i in range(len(ts[:,0])):
    storyline_dyn.append(GlobalWarming.coef + TropicalWarming.coef*ts[i,0] + VorBreak_GW.coef*ts[i,1] + SST_C.coef*ts[i,2] + SST_E.coef*ts[i,3])

#Evaluate the average precipitation change under each storyline for each region in Figure 1
SESA_storylines_dyn = {}
SESA_storylines_dyn['mindlin'] = []
SESA_storylines_dyn['diaz'] = []
SESA_storylines_dyn['junquas'] = []
SESA_storylines_dyn['zilli'] = []
SESA_storylines_dyn['gonzalez'] = []
for i in range(20000):
    SESA_storylines_dyn['mindlin'].append((storyline_dyn[i].sel(lon=slice(292,305)).sel(lat=slice(-42,-27))).mean(dim='lon').mean(dim='lat').values*5)
    SESA_storylines_dyn['diaz'].append((storyline_dyn[i].sel(lon=slice(293.75,298.75)).sel(lat=slice(-38.75,-26.25))).mean(dim='lon').mean(dim='lat').values*5)
    SESA_storylines_dyn['junquas'].append((storyline_dyn[i].sel(lon=slice(300,310)).sel(lat=slice(-32,-25))).mean(dim='lon').mean(dim='lat').values*5)
    SESA_storylines_dyn['zilli'].append((storyline_dyn[i].sel(lon=slice(311,320)).sel(lat=slice(-25,-20))).mean(dim='lon').mean(dim='lat').values*5)
    SESA_storylines_dyn['gonzalez'].append((storyline_dyn[i].sel(lon=slice(295,315)).sel(lat=slice(-40,-20))).mean(dim='lon').mean(dim='lat').values*5)
    
#Sort the precipitation values 
precip_storyline_sorted = pd.Series(SESA_storylines_dyn['mindlin']).sort_values(ascending=False)

#Define a dictionary with the precipitation extremes 
def storyline_extremes(region = 'mindlin'):
    precip_storyline_sorted = pd.Series(SESA_storylines_dyn[region]).sort_values(ascending=False)
    precip_high_extreme = precip_storyline_sorted[int(len(precip_storyline_sorted)*(99000/100000)):]
    precip_low_extreme = precip_storyline_sorted[:int(len(precip_storyline_sorted)*(1000/100000))]
    TW_index_MC = pd.Series(ts[:,0]) #pd.Series(ts[:][0])
    VB_index_MC = pd.Series(ts[:,1])#pd.Series(ts[:][1])
    CP_index_MC = pd.Series(ts[:,2])#pd.Series(ts[:][2])
    EP_index_MC = pd.Series(ts[:,3])#pd.Series(ts[:][3])

    dic = {}
    dic['high'] = {}; dic['low'] = {}
    dic['high']['tw'] = TW_index_MC[precip_high_extreme.index]
    dic['low']['tw'] = TW_index_MC[precip_low_extreme.index]
    dic['high']['vb'] = VB_index_MC[precip_high_extreme.index]
    dic['low']['vb'] = VB_index_MC[precip_low_extreme.index]    
    dic['high']['ep'] = EP_index_MC[precip_high_extreme.index]
    dic['low']['ep'] = EP_index_MC[precip_low_extreme.index]
    dic['high']['cp'] = CP_index_MC[precip_high_extreme.index]
    dic['low']['cp'] = CP_index_MC[precip_low_extreme.index]
    return dic

from sklearn import linear_model

def figure12(scatters,tw,vb,ep,cp,tw_err,vb_err,ep_err,cp_err,TW_VB_label=[r'Tropical Warming ($t_s$)',r'Vortex Breakdown delay ($t_s$)'],EP_CP_label =[r'Eastern Pacific Warming ($t_s$)',r'Central Pacific Warming ($t_s$)']):
    
    #Define means
    mean_tw = np.mean(tw)
    mean_vb = np.mean(vb)
    mean_ep = np.mean(ep)
    mean_cp = np.mean(cp)    

    #Plot-----------------------------------------------------------------------
    markers = ['<','<','v','*','D','x','x','p','+','+','d','8','X','X','^','d','d','1','2','>','>','D','D','s','.','P', 'P', '3','4']
    #markers = ['o','o','v','8','*','D','D','^','.','.','h','H','<','<','v','8','8','D','X','x','x','p','p','+','_','d', 'd', '>','X','s']
    #fig = plt.figure()
    fig, ax = plt.subplots(2,2,figsize=(14,12),dpi=300)
    for px, py, t, l, k, f in zip(tw, vb, markers, models, tw_err, vb_err):
        ax[0,0].errorbar(px, py, xerr=k, yerr=f,fmt='',ecolor='gray', elinewidth=0.5,capthick=0.001)
        ax[0,0].scatter(px, py, marker=t,label=l)
        ax[1,0].errorbar(px, py, xerr=k, yerr=f,fmt='',ecolor='gray', elinewidth=0.5,capthick=0.001)
        ax[1,0].scatter(px, py, marker=t,label=l)

    for px, py, t, l, k, f in zip(ep, cp, markers, models, ep_err, cp_err):
        ax[0,1].errorbar(px, py, xerr=k, yerr=f,fmt='',ecolor='gray', elinewidth=0.5,capthick=0.001)
        ax[0,1].scatter(px, py, marker=t,label=l)
        ax[1,1].errorbar(px, py, xerr=k, yerr=f,fmt='',ecolor='gray', elinewidth=0.5,capthick=0.001)
        ax[1,1].scatter(px, py, marker=t,label=l)

    #box = ax.get_position()
    #ax[0].set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
    confidence_ellipse(tw, vb, ax[0,0], 'no',edgecolor='red',label='80 $\%$ confidence region')
    confidence_ellipse(tw, vb, ax[0,0], 'no',chi_squared=4.6,edgecolor='k',linestyle='--',alpha=0.5,label='$\pm$ 10 $\%$ confidence regions')
    confidence_ellipse(tw, vb, ax[0,0], 'no',chi_squared=2.4,edgecolor='k',linestyle='--',alpha=0.5)
    confidence_ellipse(ep, cp, ax[0,1], ' ',edgecolor='red',label='80 $\%$ confidence region')
    confidence_ellipse(ep, cp, ax[0,1], ' ',chi_squared=4.6,edgecolor='k',linestyle='--',alpha=0.5,label='$\pm$ 10 $\%$ confidence regions')
    confidence_ellipse(ep, cp, ax[0,1], ' ',chi_squared=2.4,edgecolor='k',linestyle='--',alpha=0.5)
    confidence_ellipse(tw, vb, ax[1,0], 'no',edgecolor='red',label='80 $\%$ confidence region')
    confidence_ellipse(tw, vb, ax[1,0], 'no',chi_squared=4.6,edgecolor='k',linestyle='--',alpha=0.5,label='$\pm$ 10 $\%$ confidence regions')
    confidence_ellipse(tw, vb, ax[1,0], 'no',chi_squared=2.4,edgecolor='k',linestyle='--',alpha=0.5)
    confidence_ellipse(ep, cp, ax[1,1], ' ',edgecolor='red',label='80 $\%$ confidence region')
    confidence_ellipse(ep, cp, ax[1,1], ' ',chi_squared=4.6,edgecolor='k',linestyle='--',alpha=0.5,label='$\pm$ 10 $\%$ confidence regions')
    confidence_ellipse(ep, cp, ax[1,1], ' ',chi_squared=2.4,edgecolor='k',linestyle='--',alpha=0.5)
    ax[0,0].axvline(mean_tw, c='grey', lw=1)
    ax[0,0].axhline(mean_vb, c='grey', lw=1)    
    ax[0,0].grid()
    ax[0,0].tick_params(labelsize=18)
    ax[0,1].axvline(mean_ep, c='grey', lw=1)
    ax[0,1].axhline(mean_cp, c='grey', lw=1)    
    ax[0,1].grid()
    ax[0,1].tick_params(labelsize=18)
    ax[1,0].axvline(mean_tw, c='grey', lw=1)
    ax[1,0].axhline(mean_vb, c='grey', lw=1)    
    ax[1,0].grid()
    ax[1,0].tick_params(labelsize=18)
    ax[1,1].axvline(mean_ep, c='grey', lw=1)
    ax[1,1].axhline(mean_cp, c='grey', lw=1)    
    ax[1,1].grid()
    ax[1,1].tick_params(labelsize=18)    
    #scatter plots region1
    TW_low_tercile = scatters['region1']['low']['tw']; VB_low_tercile = scatters['region1']['low']['vb']
    TW_high_tercile = scatters['region1']['high']['tw']; VB_high_tercile = scatters['region1']['high']['vb']
    EP_low_tercile = scatters['region1']['low']['ep']; CP_low_tercile = scatters['region1']['low']['cp']
    EP_high_tercile = scatters['region1']['high']['ep']; CP_high_tercile = scatters['region1']['high']['cp']    
    ax[0,0].scatter(TW_low_tercile,VB_low_tercile,color='blue',alpha=0.2,label='drying trend ')
    ax[0,0].scatter(np.mean(TW_low_tercile),np.mean(VB_low_tercile),marker='*',s=200,color='blue')
    ax[0,0].scatter(TW_high_tercile,VB_high_tercile,color='red',alpha=0.2,label='wetting trend ')
    ax[0,0].scatter(np.mean(TW_high_tercile),np.mean(VB_high_tercile),marker='*',s=200,color='red')
    ax[0,0].set_xlabel(TW_VB_label[0],fontsize=18)
    ax[0,0].set_ylabel(TW_VB_label[1],fontsize=18)
    ax[0,0].set_title('a',loc='left',fontsize=30)
    ax[0,1].scatter(EP_low_tercile,CP_low_tercile,color='blue',alpha=0.2)   
    ax[0,1].scatter(EP_high_tercile,CP_high_tercile,color='red',alpha=0.2)
    ax[0,1].scatter(np.mean(EP_low_tercile),np.mean(CP_low_tercile),marker='*',s=200,color='blue')
    ax[0,1].scatter(np.mean(EP_high_tercile),np.mean(CP_high_tercile),marker='*',s=200,color='red')
    ax[0,1].set_xlabel(EP_CP_label[0],fontsize=18)
    ax[0,1].set_ylabel(EP_CP_label[1],fontsize=18)
    ax[0,1].set_title('b',loc='left',fontsize=30)
    region_1_high = [np.mean(TW_high_tercile),np.mean(VB_high_tercile),np.mean(EP_high_tercile),np.mean(CP_high_tercile)]
    region_1_low = [np.mean(TW_low_tercile),np.mean(VB_low_tercile),np.mean(EP_low_tercile),np.mean(CP_low_tercile)]
    
    #scatter plots region2
    TW_low_tercile = scatters['region2']['low']['tw']; VB_low_tercile = scatters['region2']['low']['vb']
    TW_high_tercile = scatters['region2']['high']['tw']; VB_high_tercile = scatters['region2']['high']['vb']
    EP_low_tercile = scatters['region2']['low']['ep']; CP_low_tercile = scatters['region2']['low']['cp']
    EP_high_tercile = scatters['region2']['high']['ep']; CP_high_tercile = scatters['region2']['high']['cp']    
    ax[1,0].scatter(TW_low_tercile,VB_low_tercile,color='blue',alpha=0.2,label='wet storylines')
    ax[1,0].scatter(np.mean(TW_low_tercile),np.mean(VB_low_tercile),marker='*',s=200,color='blue',label='mean wet storyline ')
    ax[1,0].scatter(TW_high_tercile,VB_high_tercile,color='red',alpha=0.2,label='dry storylines')
    ax[1,0].scatter(np.mean(TW_high_tercile),np.mean(VB_high_tercile),marker='*',s=200,color='red',label='mean dry storyline ')
    ax[1,0].set_xlabel(TW_VB_label[0],fontsize=18)
    ax[1,0].set_ylabel(TW_VB_label[1],fontsize=18)
    ax[1,0].set_title('c',loc='left',fontsize=30)
    #ax[1,0].legend()
    ax[1,1].scatter(EP_low_tercile,CP_low_tercile,color='blue',alpha=0.2)   
    ax[1,1].scatter(EP_high_tercile,CP_high_tercile,color='red',alpha=0.2)
    ax[1,1].scatter(np.mean(EP_low_tercile),np.mean(CP_low_tercile),marker='*',s=200,color='blue')
    ax[1,1].scatter(np.mean(EP_high_tercile),np.mean(CP_high_tercile),marker='*',s=200,color='red')
    ax[1,1].set_xlabel(EP_CP_label[0],fontsize=18)
    ax[1,1].set_ylabel(EP_CP_label[1],fontsize=18)
    ax[1,1].set_title('d',loc='left',fontsize=30)
    ax[1,1].set_xlim(-4,4)
    ax[1,1].set_ylim(-4,4)
    ax[0,0].set_xlim(-4,4)
    ax[0,0].set_ylim(-4,4)
    ax[0,1].set_xlim(-4,4)
    ax[0,1].set_ylim(-4,4)
    ax[1,0].set_xlim(-4,4)
    ax[1,0].set_ylim(-4,4)
    lgd = ax[1,0].legend(loc='upper center', bbox_to_anchor=(1, -0.2), shadow=True, ncol=5)
    plt.subplots_adjust(bottom=0.05)
    region_2_high = [np.mean(TW_high_tercile),np.mean(VB_high_tercile),np.mean(EP_high_tercile),np.mean(CP_high_tercile)]
    region_2_low = [np.mean(TW_low_tercile),np.mean(VB_low_tercile),np.mean(EP_low_tercile),np.mean(CP_low_tercile)]
    return fig, region_1_high,region_1_low,region_2_high,region_2_low

    
ep = e_index_DJF.values; cp = c_index_DJF.values
ep_err = e_index_errors_DJF.values; cp_err = c_index_errors_DJF.values
tw = tw_DJF.values; vb = vb_DJF.values
tw_err = tw_errors_DJF /gw_DJF; vb_err = vb_errors.values /gw_DJF

scatters = {}
scatters['region1'] = storyline_extremes('mindlin'); scatters['region2'] = storyline_extremes('zilli'); 
fig, region_1_high,region_1_low,region_2_high,region_2_low = figure12(scatters,std(tw),std(vb),std(ep),std(cp),tw_err,vb_err,ep_err,cp_err)

