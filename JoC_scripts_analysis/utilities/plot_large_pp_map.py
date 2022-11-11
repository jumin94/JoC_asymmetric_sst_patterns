# Subplot mean precipitation changes
import cartopy.crs as ccrs
import numpy as np
import cartopy.feature
from cartopy.util import add_cyclic_point
import matplotlib.path as mpath
import os
import glob
import pandas as pd
import xarray as xr
import netCDF4
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 0.5  # previous pdf hatch linewidth
import cartopy.util as cutil


def add_box(box,texto,color='black',text_color='black'):
    x_start, x_end = box[0], box[1]
    y_start, y_end = box[2], box[3]
    margin = 0.0007
    margin_fractions = np.array([margin, 1.0 - margin])
    x_lower, x_upper = x_start + (x_end - x_start) * margin_fractions
    y_lower, y_upper = y_start + (y_end - y_start) * margin_fractions
    box_x_points = x_lower + (x_upper - x_lower) * np.array([0, 1, 1, 0, 0])
    box_y_points = y_lower + (y_upper - y_lower) * np.array([0, 0, 1, 1, 0])
    plt.plot(box_x_points, box_y_points, transform=ccrs.PlateCarree(),linewidth=1.5, color=color, linestyle='-')
    plt.text(x_start + (x_end - x_start)*0.2, y_start + (y_end - y_start)*0.6, texto,transform=ccrs.PlateCarree( ),color=text_color)

boxes_defauls =  [[160,210,-5,5],[190,240,-5,5],[210,270,-5,5],[270,280,-10,0]]
boxes_names_default = ['Niño 4','Niño 3.4','Niño 3','Niño 1+2']
def big_map(maps,k,index_name,levels=np.arange(-.2,.25,.05),box=boxes_defauls,box_names=boxes_names_default):
    cmapPr = mpl.colors.ListedColormap(['sienna','darkgoldenrod','burlywood','wheat',
                                        'moccasin','white','white','paleturquoise','mediumaquamarine',
                                        'mediumseagreen','seagreen','darkgreen'])

    cmapPr.set_over('darkslategrey')
    cmapPr.set_under('saddlebrown')

    #k = 0
    lon = np.arange(0, 360, 2.8125)
    #index_name = 'Global Warming'
    cyclic_data, cyclic_lons = add_cyclic_point(maps[k]*86400, coord=lon)
    #cyclic_pval, lon_c_pval = add_cyclic_point(maps_pval[k].coef,lon)

    fig = plt.figure(figsize=(8, 4),dpi=300,constrained_layout=True)
    data_crs = ccrs.PlateCarree()
    proj = ccrs.PlateCarree(180)
    ax1 = plt.subplot(1,1,1,projection=proj)
    im1=ax1.contourf(cyclic_lons, maps[k].lat,cyclic_data,levels,transform=data_crs,cmap=cmapPr,extend='both')
    #levels = [cyclic_pval.min(),0.1,cyclic_pval.max()]
    #ax1.contourf(lon_c_pval, maps[k].lat,cyclic_pval,levels, transform=data_crs,levels=levels, hatches=["...", ""], alpha=0.0005)
    ax1.set_title(index_name,fontsize=14)
    ax1.add_feature(cartopy.feature.COASTLINE,alpha=.5)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    ax1.gridlines(crs=data_crs, linewidth=0.3, linestyle='-')
    ax1.set_extent([-60, 180, -70, 40], ccrs.PlateCarree(central_longitude=180))
    #Add gridlines
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.2, color='gray', alpha=0.01, linestyle='-')
    #Add niño boxes
    for i in range(len(box)):
        add_box(box[i],box_names[i],'black','red')

    plt1_ax = plt.gca()
    left, bottom, width, height = plt1_ax.get_position().bounds
    colorbar_axes = fig.add_axes([left + 0.9, bottom,0.02, height])
    cbar = fig.colorbar(im1, colorbar_axes, orientation='vertical')
    cbar.set_label(r'K/K',fontsize=14) #rotation= radianes
    cbar.ax.tick_params(axis='both',labelsize=14)
    