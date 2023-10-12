#Figure2
path_results = '/home/julia.mindlin/Tesis/JoC_paper/JoC_results'

total = xr.open_dataset(path_results+'/sst_changes/total_SST_MEM_CMIP5.nc')
total_std = xr.open_dataset(path_results+'/sst_changes/total_SST_std_MEM_CMIP5.nc')
zonal = xr.open_dataset(path_results+'/sst_changes/asymmetric_SST_MEM_CMIP5.nc')
zonal_std = xr.open_dataset(path_results+'/sst_changes/asymmetric_SST_std_MEM_CMIP5.nc')
asym = xr.open_dataset(path_results+'/sst_changes/zonally_symmetric_SST_MEM_CMIP5.nc')
asym_std = xr.open_dataset(path_results+'/sst_changes/zonally_symmetric_SST_std_MEM_CMIP5.nc')                

field = [total,total_std,zonal,zonal_std,asym,asym_std]
levels = [np.arange(0,2,.2),np.arange(0,.5,.05),np.arange(1,1.5,.05),
          np.arange(0,.4,.05),np.arange(-.5,.55,.05),np.arange(0,.4,.05)]

cmap = ['OrRd','OrRd','OrRd','OrRd','RdBu_r','OrRd']
name = ['a) sea surface temperature (SST) change - MEM mean',
        'd) sea surface temperature (SST) change - MEM std',
        'b) zonally symmetric SST change - MEM mean',
        'e) zonally symmetric SST change - MEM std',
        'c) zonally asymmetric SST change -  MEM mean',
        'f) zonally asymmetric SST change -  MEM std']
        
adjust = [0.045,-0.02,-0.085,0.045,-0.02,-0.085]

fig = plt.figure(figsize=(10,8),dpi=300,constrained_layout=True)
data_crs = ccrs.PlateCarree()
proj = ccrs.PlateCarree(180)
for i in range(6):
    ax1 = plt.subplot(3,2,i+1,projection=proj)
    lon = np.arange(0, 360, 2.8125)
    cyclic_data, cyclic_lons = add_cyclic_point(field[i], coord=lon)
    im1=ax1.contourf(cyclic_lons, asym.lat,cyclic_data,levels[i],transform=data_crs,cmap=cmap[i],extend='both')
    ax1.set_title(name[i],fontsize=12)
    ax1.add_feature(cartopy.feature.COASTLINE)
    ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    ax1.gridlines(crs=data_crs, linewidth=0.3, linestyle='-')
    ax1.set_xticks(np.linspace(0, 360,  7), crs=data_crs)
    ax1.set_yticks(np.arange(-90, 50,  10), crs=data_crs)      
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax1.xaxis.set_major_formatter(lon_formatter)
    ax1.yaxis.set_major_formatter(lat_formatter)

    ax1.set_extent([-60, 120, -40, 40], ccrs.PlateCarree(central_longitude=180))
    #Add niño boxes
    # Niño 4
    box = [160,210,-5,5]
    #add_box(box,'','green','black')
    box = [190,240,-5,5]
    #add_box(box,'','green','black')
    box = [210,270,-5,5]
    #add_box(box,'','green','black')
    box = [270,280,-10,0]
    #add_box(box,'','green','black')
    box = [180,250,5,-5]
    add_box(box,'CP','black','black',x_middle=0.4,y_middle=0.7)
    box = [260,280,0,-10]
    add_box(box,'EP','black','black',x_middle=0.08,y_middle=0.8)

    left, bottom, width, height = ax1.get_position().bounds
    colorbar_axes = fig.add_axes([left + 0.34, bottom+0.01 ,0.01, height*0.9])
    cbar = fig.colorbar(im1, colorbar_axes, orientation='vertical')
    cbar.set_label(r'K/K',fontsize=10) #rotation= radianes
    cbar.ax.tick_params(axis='both',labelsize=10)
    plt.subplots_adjust(wspace=0.4)
    
    plt.savefig('/home/julia.mindlin/Tesis/Capitulo3/scripts/CMIP6_storylines/DJF/paper_figures/sst_scaled_option_all.pdf',dpi=300,bbox_inches="tight")
    plt.savefig('/home/julia.mindlin/Tesis/Capitulo3/scripts/CMIP6_storylines/DJF/paper_figures/sst_scaled_option_all.png',dpi=300,bbox_inches="tight")
    