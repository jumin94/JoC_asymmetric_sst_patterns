#Figure 8
#Precipitation sensitiviy maps

#--------------------------------------------------IMPORTS_------------------------------------
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
import warnings
warnings.filterwarnings("ignore")


#--------------------------------------------------FUNCTIONS---------------------------------------------------
def add_box(ax,box,texto,color='black',text_color='black'):
    x_start, x_end = box[0], box[1]
    y_start, y_end = box[2], box[3]
    margin = 0.0007
    margin_fractions = np.array([margin, 1.0 - margin])
    x_lower, x_upper = x_start + (x_end - x_start) * margin_fractions
    y_lower, y_upper = y_start + (y_end - y_start) * margin_fractions
    box_x_points = x_lower + (x_upper - x_lower) * np.array([0, 1, 1, 0, 0])
    box_y_points = y_lower + (y_upper - y_lower) * np.array([0, 0, 1, 1, 0])
    ax.plot(box_x_points, box_y_points, transform=ccrs.PlateCarree(),linewidth=1.5, color=color, linestyle='--')
      
def axis_setup(ax,title,extent):
    lat_ticks=[-55,-45,-35,-25,-15,-5];lon_ticks=[-70,-60,-50,-40]
    data_crs = ccrs.PlateCarree()
    plt.title(title,fontsize=18)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.2, color='gray', alpha=0.01, linestyle='-')
    gl.xlabels_top = False; gl.ylabels_left = True; gl.ylabels_right = False;gl.xlines = False; gl.yline = False
    gl.xlocator = mticker.FixedLocator(lon_ticks)
    gl.ylocator = mticker.FixedLocator(lat_ticks)
    gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 13, 'color': 'black'}; gl.ylabel_style = {'size': 13, 'color': 'black'}
    ax.add_feature(cartopy.feature.COASTLINE,alpha=.8)
    ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.8)
    ax.gridlines(crs=data_crs, linewidth=0.3, linestyle='-')
    ax.set_extent(extent, ccrs.PlateCarree())

    return ax
    
def plot_figure(data,datap,path_fig,title,cs_coords,extent=[-180, 180, -90, 0]):
    
    lat = data[0].lat
    lon = np.arange(0,357.188,2.81)
    tw, lon_c = add_cyclic_point(data[0],lon)
    lon = np.arange(0,357.188,2.81)
    vb, lon_c = add_cyclic_point(data[1],lon)
    lon = np.arange(0,357.188,2.81)
    cp, lon_c = add_cyclic_point(data[2],lon)
    lon = np.arange(0,357.188,2.81)
    ep, lon_c = add_cyclic_point(data[3],lon)
    lon = np.arange(0,357.188,2.81)
    twp, lon_c = add_cyclic_point(datap[0],lon)
    lon = np.arange(0,357.188,2.81)
    vbp, lon_c = add_cyclic_point(datap[1],lon)
    lon = np.arange(0,357.188,2.81)
    cpp, lon_c = add_cyclic_point(datap[2],lon)
    lon = np.arange(0,357.188,2.81)
    epp, lon_c = add_cyclic_point(datap[3],lon)

    cmapPr = mpl.colors.ListedColormap(['#754308', '#955910', '#b37625', '#ca9849', '#ddbe78', '#ebd6a2',
                                        '#f6eac9', '#f5f1e6','white','white','lightcyan','paleturquoise','turquoise',
                                        'mediumaquamarine','mediumseagreen','seagreen','green','darkgreen'])
                                                                                                            
    cmapPr.set_over('darkslategrey')
    cmapPr.set_under('#543005')
    
    #SoutherHemisphere Stereographic
    fig = plt.figure(figsize=(10, 15),dpi=300,constrained_layout=True)
    projection = ccrs.PlateCarree(central_longitude=0)
    data_crs = ccrs.PlateCarree()
    
    ax1 = plt.subplot(2,2,1,projection=projection);ax1 = axis_setup(ax1,title[0],extent)
    ax2 = plt.subplot(2,2,2,projection=projection);ax2 = axis_setup(ax2,title[1],extent)
    ax3 = plt.subplot(2,2,3,projection=projection);ax3 = axis_setup(ax3,title[2],extent)
    ax4 = plt.subplot(2,2,4,projection=projection);ax4 = axis_setup(ax4,title[3],extent)
    clevels = np.arange(-.18,.2,0.02)
    im=ax1.contourf(lon_c, lat, tw,clevels,transform=data_crs,cmap='BrBG',extend='both')
    ax1.contour(lon_c, lat,tw,levels=[0], transform=data_crs,cmap='Greys_r')
    #Transecta
    #ax1.plot(cs_coords[0],cs_coords[1],transform=data_crs,linewidth=2,color='r')
    #ax1.text(cs_coords[0][0]-4,cs_coords[1][0]-3,'SW',transform=data_crs,fontsize=19,color='k')
    #ax1.text(cs_coords[0][-1],cs_coords[1][-1]+1.2,'NE',transform=data_crs,fontsize=19,color='k')
    ax1.add_feature(cartopy.feature.COASTLINE,alpha=.5)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)    
    im=ax2.contourf(lon_c, lat, vb,clevels,transform=data_crs,cmap='BrBG',extend='both')
    ax2.contour(lon_c, lat,vb,levels=[0], transform=data_crs,cmap='Greys_r')
    ax2.add_feature(cartopy.feature.COASTLINE,alpha=.5)
    ax2.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)    
    im=ax3.contourf(lon_c, lat, cp,clevels,transform=data_crs,cmap='BrBG',extend='both')
    ax3.contour(lon_c, lat,cp,levels=[0], transform=data_crs,cmap='Greys_r')
    ax3.add_feature(cartopy.feature.COASTLINE,alpha=.5)
    ax3.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)    
    im=ax4.contourf(lon_c, lat, ep,clevels,transform=data_crs,cmap='BrBG',extend='both')
    ax4.contour(lon_c, lat,ep,levels=[0], transform=data_crs,cmap='Greys_r')
    ax4.add_feature(cartopy.feature.COASTLINE,alpha=.5)
    ax4.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    plevels = [twp.min(),0.2,twp.max()]
    ax1.contourf(lon_c, lat,twp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax2.contourf(lon_c, lat,twp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax3.contourf(lon_c, lat,twp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax4.contourf(lon_c, lat,twp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    plevels = [vbp.min(),0.2,vbp.max()]
    ax1.contourf(lon_c, lat,vbp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax2.contourf(lon_c, lat,vbp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax3.contourf(lon_c, lat,vbp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax4.contourf(lon_c, lat,vbp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    plevels = [cpp.min(),0.2,cpp.max()]
    ax1.contourf(lon_c, lat,cpp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax2.contourf(lon_c, lat,cpp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax3.contourf(lon_c, lat,cpp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax4.contourf(lon_c, lat,cpp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    plevels = [epp.min(),0.2,epp.max()]
    ax1.contourf(lon_c, lat,epp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax2.contourf(lon_c, lat,epp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax3.contourf(lon_c, lat,epp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    ax4.contourf(lon_c, lat,epp, transform=data_crs,levels=plevels, hatches=["...", ""],alpha=0.0005)
    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    colorbar_axes = fig.add_axes([left+0.4, bottom, 0.03, height*2.2])
    cbar = fig.colorbar(im, colorbar_axes, orientation='vertical')
    cbar.set_label('mm day$^{-1}$ K$^{-1}$',fontsize=18) #rotation = radianes
    cbar.ax.tick_params(axis='both',labelsize=18)
    plt.subplots_adjust(hspace=0.15,wspace=0.3)
    plt.clf
    
    return fig

def transecta(lons,lats,rds):
    rd_cs = []
    for rd in rds:
        transecta=[]
        for k in range(len(lons)):
            transecta.append(rd.sel(lat=lats[k],method='nearest').sel(lon=lons[k],method='nearest'))
        rd_cs.append(np.array(transecta))
    return rd_cs

def cs_plot(lons,css,names,title):
    fig = plt.figure(figsize=(15,10)); ax = plt.subplot(111)
    for cs,name in zip(css,names):
        ax.plot(cs,'-',linewidth=3,label = name)
    ax.axhline(0,0,10,color='k',linestyle='--')
    ax.set_ylabel('precipitation change [mm/day]',fontsize=22)
    ax.set_xlabel('AB transect',fontsize=22)
    ax.set_yticks([-.1,0,.1,.2,.3]); ax.set_yticklabels([-.1,0,.1,.2,.3],fontsize=22)
    #ax.set_xticks(np.arange(0,len(lons))[::3]);ax.set_xticklabels(xticks[::3],fontsize=12)
    ax.tick_params(axis='both',labelsize=22)
    plt.title(title)
    plt.legend()
    return fig

def fdr(pvals,alpha):
    p = pvals.stack(aux=('lat','lon'))
    p_sorted = p.sortby(p)
    p_fdr = np.arange(1,len(p_sorted.aux)+1,1)/len(p_sorted.aux)*alpha
    auxiliar = p_sorted.copy()
    auxiliar.values= p_fdr
    final = (p_sorted/ auxiliar)
    final = final.unstack().where(final.unstack() < 1)
    return final 


#--------------------------------------------------OPEN DATA--------------------------------------------
import matplotlib as mpl
from cartopy.util import add_cyclic_point
import matplotlib.pyplot as plt
import xarray as xr 
import numpy as np 

path_maps = '/home/julia.mindlin/Tesis/JoC_paper/JoC_results/sensitivity_maps_CMIP5/ua'
#path_maps = '/home/julia.mindlin/Tesis/Capitulo3/scripts/CMIP6_storylines/DJF/sensitivity_maps/pr/C_E_std_max'
#Create sensitivity maps
path_fig = '/home/julia.mindlin/Tesis/JoC_paper/JoC_figures'

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

#Define masks
GlobalWarming = xr.open_dataset(path_maps+'/Aij.nc')*86400
GW = GlobalWarming
mask_mindlin = GW.where(GW.lon<305).where(GW.lon>292).where(GW.lat>-42).where(GW.lat<-27) / GW.where(GW.lon<305).where(GW.lon>292).where(GW.lat>-42).where(GW.lat<-27)
#mask_mindlin = mask_mindlin.reindex(lat=list(reversed(mask_mindlin.lat)))
#Box from IPCC
import regionmask
mask_IPCC = regionmask.defined_regions.ar6.land.mask(GW).where(regionmask.defined_regions.ar6.land.mask(GW) == 14)/14
#mask_IPCC = mask_IPCC.reindex(lat=list(reversed(mask_IPCC.lat)))
#Box from paper Lean
#38.75∘S–26.25∘S, 66.25–61.25∘W
mask_diaz = GW.where(GW.lon<298.75).where(GW.lon>293.75).where(GW.lat>-38.75).where(GW.lat<-26.25) / GW.where(GW.lon<298.75).where(GW.lon>293.75).where(GW.lat>-38.75).where(GW.lat<-26.25)
#mask_diaz = mask_diaz.reindex(lat=list(reversed(mask_diaz.lat)))
#Box from Clementine
#32°S:25°S,60°W:50°W
mask_junquas = GW.where(GW.lon<310).where(GW.lon>300).where(GW.lat>-32).where(GW.lat<-25) / GW.where(GW.lon<310).where(GW.lon>300).where(GW.lat>-32).where(GW.lat<-25)
#mask_junquas = mask_junquas.reindex(lat=list(reversed(mask_junquas.lat)))
#Box for San Paulo
#32°S:25°S,60°W:50°W
mask_junquas = GW.where(GW.lon<310).where(GW.lon>300).where(GW.lat>-32).where(GW.lat<-25) / GW.where(GW.lon<310).where(GW.lon>300).where(GW.lat>-32).where(GW.lat<-25)
#mask_junquas = mask_junquas.reindex(lat=list(reversed(mask_junquas.lat)))

#Define plotting data
rds = [TropicalWarming.coef,VorBreak_GW.coef,SST_C.coef,SST_E.coef]
rds_pval = [TropicalWarmingp.coef,VorBreak_GWp.coef,SST_Cp.coef,SST_Ep.coef]
#rds_pval = [Fp.coef,Fp.coef,Fp.coef,Fp.coef]
names  = ['(a) Tropical Warming','(b) Vortex Breakdown delay','(c) Central Pacific Warming','(d) Eastern Pacific Warming']
#Cross section
lons = np.linspace(290,322,10)
lats = np.linspace(-38,-18,10)
cs_coords = [lons,lats]
#Plot
fig = plot_figure(rds,rds_pval,path_fig,names,cs_coords,extent=[-73,-34,-57,0])
fig.savefig(path_fig+'/precip_rds_CMIP5.png')