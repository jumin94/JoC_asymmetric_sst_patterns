#Figure 4

#--------------------------------------------------IMPORTS----------------------------------------
import numpy as np
import pandas as pd
import xarray as xr
import os, fnmatch
import glob
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from numpy import linalg as la
from sklearn.linear_model import LinearRegression
from sklearn import linear_model
import warnings
warnings.filterwarnings("ignore")

#--------------------------------------------------FUNCTIONS-----------------------------------------
def confidence_ellipse(x ,y, ax, corr,chi_squared=3.21, facecolor='none',**kwargs):
 if x.size != y.size:
  raise ValueError('x and y must be the same size')

 cov = np.cov(x,y)
 pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
 eigval, eigvec = la.eig(cov)
 largest_eigval = np.argmax(eigval)
 largest_eigvec = eigvec[:,largest_eigval]
 smallest_eigval = np.argmin(eigval)
 smallest_eigvec = eigvec[:,smallest_eigval]
 lamda1 = np.max(eigval)
 lamda2 = np.min(eigval)

 scale_x = np.sqrt(lamda1)
 scale_y = np.sqrt(lamda2)
 if corr == 'no':
    angle = 90.0 #np.arctan(smallest_eigvec[0]/smallest_eigvec[1])*180/np.pi
 else:
    angle = np.arctan(smallest_eigvec[0]/smallest_eigvec[1])*180/np.pi

 # Using a special case to obtain the eigenvalues of this
 # two-dimensionl dataset. Calculating standard deviations

 ell_radius_x = scale_x*np.sqrt(chi_squared)
 ell_radius_y = scale_y*np.sqrt(chi_squared)
 ellipse = Ellipse((0, 0), width=ell_radius_x * 2,height=ell_radius_y * 2, angle = -angle, facecolor=facecolor,**kwargs)

 # Calculating x mean
 mean_x = np.mean(x)
 # calculating y mean
 mean_y = np.mean(y)

 transf = transforms.Affine2D() \
     .translate(mean_x, mean_y)

 ellipse.set_transform(transf + ax.transData)
 return ax.add_patch(ellipse), print(angle), ellipse

def plot_ellipse(x,y,xerr,yerr,title='October - November',corr='no',x_label='Eastern Pacific Warming [K K$^{-1}$]',y_label='Central Pacific Warming [K K$^{-1}$]',fig_name='ellipse'):
    #Compute regression y on x
    x1 = x.reshape(-1, 1)
    y1 = y.reshape(-1, 1)
    linear_regressor = LinearRegression()  # create object for the class
    reg = linear_regressor.fit(x1, y1)  # perform linear regression
    X_pred = np.linspace(np.min(x)-1, np.max(x)+0.5, 31)
    X_pred = X_pred.reshape(-1, 1)
    Y_pred = linear_regressor.predict(X_pred)  # make predictions
    c = reg.coef_

    #Compute regression x on y
    reg2 = linear_regressor.fit(y1, x1)  # perform linear regression
    Y_pred2 = np.linspace(np.min(y), np.max(y), 31)
    Y_pred2 = Y_pred2.reshape(-1, 1)
    X_pred2 = linear_regressor.predict(Y_pred2)  # make predictions
    c2 = reg2.coef_

    #Define limits
    min_x = np.min(x) - 0.2*np.abs(np.max(x) - np.min(x))
    max_x = np.max(x) + 0.2*np.abs(np.max(x) - np.min(x))
    max_y = np.max(y) + 0.2*np.abs(np.max(y) - np.min(y))
    max_y = np.min(y) - 0.2*np.abs(np.max(y) - np.min(y))
    mean_x = np.mean(x)
    mean_y = np.mean(y)

    #Calcular las rectas x = y, x = -y
    Sx = np.std(x)
    Sy = np.std(y)
    S_ratio = Sy/Sx
    YeqX = S_ratio*X_pred - S_ratio*mean_x + mean_y
    YeqMinsX = S_ratio*mean_x + mean_y - S_ratio*X_pred


    #Plot-----------------------------------------------------------------------
    markers = ['<','<','v','*','D','x','x','p','+','+','d','8','X','X','^','d','d','1','2','>','>','D','D','s','.','P', 'P', '3','4']
    fig = plt.figure(dpi=300)
    fig, ax = plt.subplots()
    for px, py, t, l, k, f in zip(x, y, markers, models, xerr, yerr):
       ax.errorbar(px, py, xerr=k, yerr=f,fmt='',ecolor='gray', elinewidth=0.5,capthick=0.001)
       ax.scatter(px, py, marker=t,label=l)

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
    confidence_ellipse(x, y, ax, corr,edgecolor='red',label='80 $\%$ confidence region')
    confidence_ellipse(x, y, ax, corr,chi_squared=4.6,edgecolor='k',linestyle='--',alpha=0.5,label='$\pm$ 10 $\%$ confidence regions')
    confidence_ellipse(x, y, ax, corr,chi_squared=2.4,edgecolor='k',linestyle='--',alpha=0.5)
    ax.axvline(mean_x, c='grey', lw=1)
    ax.axhline(mean_y, c='grey', lw=1)
    #ax.plot(X_pred,YeqX ,color='black')
    #ax.plot(X_pred,YeqMinsX ,color='grey')
    ax.grid()
    ax.tick_params(labelsize=18)
    if corr == 'yes':
        r = np.corrcoef(x,y)[0,1]; chi = (1.26**2)*2
        ts1 = np.sqrt(((1-r**2)/(2*(1-r)))*chi)
        ts2 = np.sqrt(((1-r**2)/(2*(1+r)))*chi)
        story_x1 = [mean_x + ts1*np.std(x)]
        story_x2 = [mean_x - ts1*np.std(x)]
        story_y_red1 = [mean_y + ts1*np.std(y)]
        story_y_red2 =[mean_y - ts1*np.std(y)]
        story_x3 = [mean_x - ts2*np.std(x)]
        story_x4 = [mean_x + ts2*np.std(x)]
        story_y_red3 = [mean_y + ts2*np.std(y)]
        story_y_red4 =[mean_y - ts2*np.std(y)]
    else:
        story_x = [mean_x + 1.26*np.std(x),mean_x - 1.26*np.std(x)]
        story_y_red = [mean_y + 1.26*np.std(y),mean_y + 1.26*np.std(y)]
        story_y_blue =[mean_y - 1.26*np.std(y),mean_y - 1.26*np.std(y)]
    lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), shadow=True, ncol=5)
    plt.subplots_adjust(bottom=0.05)
    plt.xlabel(x_label,fontsize=18)
    plt.ylabel(y_label,fontsize=18)
    plt.text(mean_x + 1.8*np.std(x),mean_y + 2*np.std(y),'R='+str(round(np.corrcoef(x,y)[0,1],3)),fontsize=15)
    plt.title(title)
    plt.savefig(path_plots+'/'+fig_name+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')

    
#--------------------------------------------------CODE------------------------------------------------------
#Open index data
#Open data and define models
models = ['ACCESS-CM2','ACCESS-ESM1-5','BCC-CSM2-MR','CAMS-CSM1-0','CanESM5','CESM2_','CESM2-WACCM','CMCC-CM2-SR5',
          'CNRM-CM6-1','CNRM-ESM2-1','EC-Earth3','FGOALS-g3','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IITM-ESM',
          'INM-CM4-8', 'INM-CM5-0','KACE-1-0-G','MIROC6','MIROC-ES2L','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0',
          'NESM3','NorESM2-LM','NorESM2-MM','TaiESM1','UKESM1-0-LL']


path_indices = '/home/julia.mindlin/Tesis/BSC/indices'
#Open indices for MA season
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
vb =  VB.iloc[:,2]
#Open index errors for ON season
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

tw_DJF = tw_DJF/gw_DJF
    
x = e_index_DJF.values 
y = c_index_DJF.values
xerr = e_index_errors_DJF.values
yerr = c_index_errors_DJF.values

path_plots = '/home/julia.mindlin/Tesis/JoC_paper'
plot_ellipse(x,y,xerr,yerr,' ','yes',fig_name='C_E_ellipse')

x = tw_DJF.values
y = vb.values
xerr = tw_errors_DJF/gw_DJF
yerr = vb_errors.values /gw_DJF

plot_ellipse(x,y,xerr,yerr,' ','no','Tropical Warming [K K$^{-1}$]','Vortex Breakdown delay [K K$^{-1}$]',fig_name='TW_VB_ellipse')