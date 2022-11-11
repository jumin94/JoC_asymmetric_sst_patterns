#Figure 1 

#--------------------------------------------------IMPORTS--------------------------------------
import cartopy
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.util import add_cyclic_point
import cartopy.util as cutil
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import matplotlib.ticker as mticker
from matplotlib.ticker import MaxNLocator as  MaxNLocator
from matplotlib.colors import BoundaryNorm as BoundaryNorm
import matplotlib.ticker as mticker
import matplotlib
from matplotlib.pyplot import cm
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 0.5  # previous pdf hatch linewidth


#-------------------------------------------------FUNCTIONS------------------------------------------------
def add_box(ax,box,texto,color='black',text_color='black'):
    x_start, x_end = box[0], box[1]
    y_start, y_end = box[2], box[3]
    margin = 0.0007
    margin_fractions = np.array([margin, 1.0 - margin])
    x_lower, x_upper = x_start + (x_end - x_start) * margin_fractions
    y_lower, y_upper = y_start + (y_end - y_start) * margin_fractions
    box_x_points = x_lower + (x_upper - x_lower) * np.array([0, 1, 1, 0, 0])
    box_y_points = y_lower + (y_upper - y_lower) * np.array([0, 0, 1, 1, 0])
    ax.plot(box_x_points, box_y_points, transform=ccrs.PlateCarree(),linewidth=1.5, color=color, linestyle='-')
    #ax.text(x_start + (x_end - x_start)*0.4, y_start + (y_end - y_start)*0.4, texto,transform=ccrs.PlateCarree( ),color=text_color,fontsize=20)
      
def axis_setup(ax,title,extent,i=0):
    lat_ticks=[-55,-45,-35,-25,-15,-5];lon_ticks=[-70,-60,-50,-40]
    data_crs = ccrs.PlateCarree()
    if i == 0:
        ax.text(321,-22, 'A',color='black',transform=data_crs,fontsize=18)
        ax.text(305 - 2.8, -42 + 2.5, 'B',transform=data_crs,color='black',fontsize=18)
        ax.text(298.75 - 2.5, -38.75 + 1, 'C',transform=data_crs,color='black',fontsize=18)
        ax.text(310 - 2.5, -32 + 2, 'D',transform=data_crs,color='black',fontsize=18)
        ax.text(315 - 2.8, -40 + 2.5, 'E',transform=data_crs,color='black',fontsize=18)
        ax.text(311 + 7 ,-25 + 2.7,'F',transform=data_crs,color='black',fontsize=18)
        levels = [positives.min(),25,positives.max()]
        ax.contourf(positives.lon, positives.lat, positives,levels, transform=data_crs,levels=levels, hatches=[" ", "."], alpha=0)
        levels = [negatives.min(),-25,negatives.max()]
        ax.contourf(negatives.lon, negatives.lat, negatives,levels, transform=data_crs,levels=levels, hatches=[".", " "], alpha=0)
        levels = [gamma_forced.coef.min(),.8,gamma_forced.coef.max()]
        ax.contourf(gamma_forced.lon, gamma_forced.lat, gamma_forced.coef.values,levels, transform=data_crs,levels=levels, hatches=[" ", "o"], alpha=0)
        cnt = ax.contour(mask_IPCC.lon,mask_IPCC.lat,mask_IPCC.fillna(0),levels=[0,1],colors='k')
        box = [305,292,-42,-27] #Mindlin
        add_box(ax,box,' ','black','black') 
        box = [298.75,293.75,-38.75,-26.25]
        add_box(ax,box,' ','black','black') #Diaz
        box = [310,300,-32,-25] #Junquas
        add_box(ax,box,' ','black','black')
        box = [295,315,-40,-20] #Gonzalez
        add_box(ax,box,' ','black','black') 
        box = [311,320,-25,-20] #Zilli
        add_box(ax,box,' ','black','black') 
    elif i == 1:
        levels = [positives_model.min(),2,positives_model.max()]
        ax.contourf(positives_model.lon, positives_model.lat, positives_model,levels, transform=data_crs,levels=levels, hatches=[" "," ."], alpha=0)
        levels = [negatives_model.min(),-2,negatives_model.max()]
        ax.contourf(negatives_model.lon, negatives_model.lat, negatives_model,levels, transform=data_crs,levels=levels, hatches=[" ." ""], alpha=0)
        #levels = [gamma_forced.coef.min(),.8,gamma_forced.coef.max()]
        #ax.contourf(gamma_forced.lon, gamma_forced.lat, gamma_forced.coef.values,levels, transform=data_crs,levels=levels, hatches=[" ", "o"], alpha=0)
    plt.title(title,fontsize=18)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.2, color='gray', alpha=0.01, linestyle='-')
    gl.xlabels_top = False; gl.ylabels_left = True; gl.ylabels_right = False;gl.xlines = False; gl.yline = False
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 20, 'color': 'black'}; gl.ylabel_style = {'size': 20, 'color': 'black'}
    ax.add_feature(cartopy.feature.COASTLINE,alpha=.5)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    ax.gridlines(crs=data_crs, linewidth=0.3, linestyle='-')
    ax.set_extent(extent, ccrs.PlateCarree())
    return ax


def plot_figure(mem_in,model_in,path_fig,title,extent=[-180, 180, -90, 0]):
    
    lat_mem = mem_in.lat
    lat_model = model_in.lat
    lon = np.arange(0,357.188,2.81)
    mem, lon_c = add_cyclic_point(mem_in,lon)
    lon = np.arange(0,357.188,2.81)
    model, lon_c = add_cyclic_point(model_in,lon)
    
    cmapPr = mpl.colors.ListedColormap(['#754308', '#955910', '#b37625', '#ca9849', '#ddbe78', '#ebd6a2',
                                        '#f6eac9', '#f5f1e6','white','white','lightcyan','paleturquoise','turquoise',
                                        'mediumaquamarine','mediumseagreen','seagreen','green','darkgreen'])
                                                                                                            
    cmapPr.set_over('darkslategrey')
    cmapPr.set_under('#543005')
    
    #SESA
    fig = plt.figure(figsize=(10, 10),dpi=300,constrained_layout=True)
    projection = ccrs.PlateCarree(central_longitude=0)
    data_crs = ccrs.PlateCarree()
    
    ax1 = plt.subplot(1,2,1,projection=projection)
    clevels = np.arange(-0.25,0.3,0.05)
    im=ax1.contourf(lon_c, lat_mem, mem,clevels,transform=data_crs,cmap='BrBG',extend='both')
    ax1 = axis_setup(ax1,title[0],extent,0)
    ax2 = plt.subplot(1,2,2,projection=projection)
    clevels = np.arange(-0.25,0.3,0.05)
    im=ax2.contourf(lon_c, lat_model, model,clevels,transform=data_crs,cmap='BrBG',extend='both')
    ax2 = axis_setup(ax2,title[1],extent,1)
    
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    colorbar_axes = fig.add_axes([left+0.48, bottom-0.05, 0.03, height*1.2])
    cbar = fig.colorbar(im, colorbar_axes, orientation='vertical')
    cbar.set_label('mm day$^{-1}$',fontsize=22) #rotation = radianes
    cbar.ax.tick_params(axis='both',labelsize=22)

    #plt.savefig(path_fig+'/signal_noise_'+model+'.png')
    plt.clf
    
    return fig


def change_variability(model,n):
    experiments = ['historical','ssp585']
    #Evaluo el cambio climatologico
    periods = [[1940,1969],[2070,2099]]
    pr_hist = dic_ens_members[experiments[0]][model][n-1].sel(time=slice(str(periods[0][0]),str(periods[0][1])))
    pr_ssp = dic_ens_members[experiments[1]][model][n-1].sel(time=slice(str(periods[1][0]),str(periods[1][1])))
    pr_hist_DJF = pr_hist.sel(time=cross_year_season(pr_hist['time.month'],[12,2]))
    pr_ssp_DJF = pr_ssp.sel(time=cross_year_season(pr_ssp['time.month'],[12,2]))
    mean_change = pr_ssp_DJF.mean(dim='time').pr - pr_hist_DJF.mean(dim='time').pr
    #Calculo la media estacional de tres meses centrada 
    seasonal_means = dic_ens_members[experiments[0]][model][n-1].pr[2:-1,:,:].sel(time=slice('1850','2014')).rolling(time=3,center=True).mean().dropna('time') 
    #Me quedo con un valor cada tres, que corresponde a las medias centradas en J
    seasonal_means_DJF = seasonal_means[::3,:,:]
    #Climate mean standard deviation
    climate_mean_30year = seasonal_means_DJF.rolling(time=30,center=True).mean().dropna('time')[::30].std(dim='time')
    #Interannual variability mean value for 30 year values	
    standard_deviation_DJF = seasonal_means_DJF.rolling(time=30,center=True).std()
    non_overlap_std_DJF = standard_deviation_DJF.dropna('time')[::30]
    variability_mean_1year = non_overlap_std_DJF.mean(dim='time')    
    return mean_change,climate_mean_30year

def cargo_todos_ens_members(models,experiments,ruta,var):
    os.chdir(ruta)
    os.getcwd()
    dic = {}
    for experiment in experiments:
        dic[experiment] = {}
        for model in models:
            listOfFiles = os.listdir(ruta+'/'+experiment+'/'+var+'/all_ensemble_members')
            dic[experiment][model] = []
            pattern1 = "*"+model+"*"+experiment+"*"
            for entry in listOfFiles:
                if fnmatch.fnmatch(entry,pattern1):
                    dato = xr.open_dataset(ruta+'/'+experiment+'/'+var+'/all_ensemble_members'+'/'+entry)
                    dic[experiment][model].append(dato)
    return dic

#------------------------------------------------------CODE-----------------------------------------
#Open model data
path_indices = '/home/julia.mindlin/Tesis/BSC/indices'
season_name = 'DJF'
gw = pd.read_csv(path_indices+'/GW_index_'+season_name+'.csv')

#Open a map to use as grid reference
path_maps = '/home/julia.mindlin/Tesis/Capitulo3/scripts/CMIP6_storylines/DJF/sensitivity_maps/pr/C_E_std_max'
GlobalWarming = xr.open_dataset(path_maps+'/Aij.nc')*86400
GW = GlobalWarming

path_fig = '/home/julia.mindlin/Tesis/JoC_paper'
#Define regional boxes
mask_mindlin = GW.where(GW.lon<305).where(GW.lon>292).where(GW.lat>-42).where(GW.lat<-27) / GW.where(GW.lon<305).where(GW.lon>292).where(GW.lat>-42).where(GW.lat<-27)
#mask_mindlin = mask_mindlin.reindex(lat=list(reversed(mask_mindlin.lat)))
#Box from IPCC
import regionmask
mask_IPCC = regionmask.defined_regions.ar6.land.mask(GW).where(regionmask.defined_regions.ar6.land.mask(GW) == 14)/14
#mask_IPCC = mask_IPCC.reindex(lat=list(reversed(mask_IPCC.lat)))
#Box from paper Vera & Diaz
#38.75∘S–26.25∘S, 66.25–61.25∘W
mask_diaz = GW.where(GW.lon<298.75).where(GW.lon>293.75).where(GW.lat>-38.75).where(GW.lat<-26.25) / GW.where(GW.lon<298.75).where(GW.lon>293.75).where(GW.lat>-38.75).where(GW.lat<-26.25)
#mask_diaz = mask_diaz.reindex(lat=list(reversed(mask_diaz.lat)))
#Box from Junquas
#32°S:25°S,60°W:50°W
mask_junquas = GW.where(GW.lon<310).where(GW.lon>300).where(GW.lat>-32).where(GW.lat<-25) / GW.where(GW.lon<310).where(GW.lon>300).where(GW.lat>-32).where(GW.lat<-25)
#mask_junquas = mask_junquas.reindex(lat=list(reversed(mask_junquas.lat)))


#Open all ensemble members for the models that I will use as example
#Figure 1 in the introduction is to motivate that there are models which show a drying trend in SESA
#I load all the ensemble members of the models that have this behavior
models_fig1 = ['NESM3','MPI-ESM1-2-LR'] 
experiments = ['historical','ssp585']
os.chdir(ruta)
os.getcwd()

ruta = '/datos/julia.mindlin/CMIP6_ensambles/preprocesados'
var = 'mon/pr'
dic_ens_members = cargo_todos_ens_members(models_fig1,experiments,ruta,var)

#Evaluate models signal to noise
signal_noise = {}
change = {}
for model in models_fig1:
    signal_noise[model] = []
    change[model] = []
    for i in range(len(dic_ens_members[experiments[1]][model])):
        CH,VAR = change_variability(model,i+1)
        SN = CH/VAR
        signal_noise[model].append(SN)
        change[model].append(CH)

import matplotlib.pyplot as plt
for model in models_fig1:
    for i in range(len(signal_noise[model])):
        fig = plt.figure()
        (signal_noise[model][i]).sel(lat = slice(-50,-20)).sel(lon=slice(290,310)).where([(signal_noise[model][i]<=-1) & (change[model][i]<0)][0]).plot()
        
        
dic_masks = {}
import matplotlib.pyplot as plt
for model in models_fig1:
    for i in range(len(signal_noise[model])):
        fig = plt.figure()
        (signal_noise[model][i]).sel(lat = slice(-50,-20)).sel(lon=slice(290,310)).where([(signal_noise[model][i]<=-1) & (change[model][i]<0)][0]).plot()
        
#Genero figura
model = models_fig1[1]
mask_pos = [((signal_noise[model][0]>=1) & (change[model][0]>0))]
mask_neg = [((signal_noise[model][0]<=-1) & (change[model][0]<0))]
mask_positive = (signal_noise[model][0]).where(mask_pos[0])/(signal_noise[model][0]).where(mask_pos[0]).fillna(0)
mask_negative = (signal_noise[model][0]).where(mask_neg[0])/(signal_noise[model][0]).where(mask_neg[0]).fillna(0)
for i in range(len(signal_noise[model])-1):
    mask_pos = [((signal_noise[model][i+1]>=1) & (change[model][i+1]>0))]
    mask_neg = [((signal_noise[model][i+1]<=-1) & (change[model][i+1]<0))]
    mask_positive += (signal_noise[model][i+1]).where(mask_pos[0])/(signal_noise[model][i+1]).where(mask_pos[0]).fillna(0)
    mask_negative += (signal_noise[model][i+1]).where(mask_neg[0])/(signal_noise[model][i+1]).where(mask_neg[0]).fillna(0)

mean_signal_noise = signal_noise[model][0]
for i in range(len(signal_noise[model])-1):
    mean_signal_noise = signal_noise[model][i]
    
mean_signal_noise = mean_signal_noise / len(signal_noise[model])

dic_masks[model] = {}
dic_masks[model]['positive_change_signal_noise'] = mask_positive
dic_masks[model]['negative_change_signal_noise'] = mask_negative

#-----------------------------------------------FIGURE---------------------------------------------
model = models_fig1[1] #MPI
positives_model = (signal_noise[model][8]).where([(signal_noise[model][8]>=1) & (change[model][8]>0)][0]) #Ensemble member 9
negatives_model = (signal_noise[model][8]).where([(signal_noise[model][8]<=-1) & (change[model][8]<0)][0]) #Ensemble member 9

path = '/home/julia.mindlin/Tesis/JoC_paper/JoC_results'
gamma_forced = xr.open_dataset(path+'/gamma_forced_CMIP6.nc')
positives = xr.open_dataset(path+'/number_of_models_positive_trend_CMIP6.nc')
negatives = xr.open_dataset(path+'/number_of_models_negative_trend_CMIP6.nc')      

fig = plot_figure(GlobalWarming.coef,beta[model]*86400/gw.iloc[:,2].values[23],path_fig,['Multimodel Ensemble Mean',model],extent=[-73,-34,-57,0])
fig.savefig(path_fig+'figure1_signal_noise_+'+model+'.png')