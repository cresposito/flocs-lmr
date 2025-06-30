import os

import numpy as np
import pandas as pd
import geopandas as gpd
import contextily as cx
from shapely.geometry import Point
from shapely.geometry import LineString
import matplotlib.pyplot as plt
import importlib

import sedProfiles as sp

# %%
baseDataSetFolder=r'C:\Users\cesposito\The Water Institute of the Gulf\P-00703_NSF_Caltech - General\Data\\West Bay Report'
adcp_folder = os.path.join(baseDataSetFolder,'\ADCP Data')

# filename, and sheet for the detailed grain size data.
fn_gs = os.path.join(baseDataSetFolder,r'Grain Size Data\WB_Grain Size Compilation.xlsx')
sh_gs = 'WestBay_first_4_tab'
###THERE IS A SECOND TAB THAT NEEDS TO BE ADDED###
# sh_gs = 'WestBay_last_4_tab'

#parameters to include from grain size file header
params_fn_gs = ['File ID','Sample ID']     

# filename and sheet for station locations
fn_loc = os.path.join(baseDataSetFolder,r'Sample Locations\WestBay_Locations.xlsx')
sh_loc = 'IsoSamplingLocs'

# filename and sheet with detailed information about all of the sampling (ADCP, and isokinetic)
fn_smp = os.path.join(baseDataSetFolder,r'ADCP Data\WestBay_StationaryADCP&Isokinetic.xlsx')
sh_smp = 'WestBay'

# buffer for selecting relevant ADCP points, in meters
buffer = 100

# %%
# read the detailed grain size data
df_gs = sp.read_sed_samples(fn_gs,sh_gs,params_fn_gs)

# read location information
headers = ['station','lat','lon']
loc = pd.read_excel(fn_loc,sheet_name=sh_loc, header = None, names = headers)
loc = loc.drop_duplicates() #some are duplicated
geometry = [Point(xy) for xy in zip(loc['lon'], loc['lat'])]
gdf = gpd.GeoDataFrame(loc, geometry=geometry)
gdf.set_crs(epsg=4326, inplace=True)
loc = gdf


# read sampling information
smp = pd.read_excel(fn_smp,sheet_name = sh_smp)


# %%
# organize data into profiles
profiles_WB = sp.profiles_from_sampling_info(smp,loc)

# add grain size data
profiles_WB = sp.profiles_add_gs(profiles_WB,df_gs,sample_id_style=2)

# add adcp data
# ens_start = profiles[1]['adcp']['ADCP Ensemble Start'][0]
# ens_end = profiles[1]['adcp']['ADCP Ensemble End'][0]
# fn = sp.find_file(adcp_folder,profiles[1]['adcp']['Stationary ADCP File Name'][0])
# V_mag, t = sp.read_adcp_section(fn,ens_start,ens_end,profiles[1]['Total Depth (m)'])

# profiles_WB = sp.read_adcp_sections(profiles_WB,adcp_folder)

#%% read west bay data

#preread all adcp files
adcp_fns = []
for fns in smp['Stationary ADCP File Name']:
    #list of adcp file names. if there is a file with an r suffix, change it to t
    adcp_fns = adcp_fns + fns.replace(' ','').replace('r.','t.').split(';')    

adcp_fns = list(set(adcp_fns))

pth = r'C:\Users\cesposito\THE WATER INSTITUTE OF THE GULF\P-00703_NSF_Caltech - General\Data\West Bay Report\ADCP Data\Trip 1c 5-6 May 2009\ASCII'
A={}
for fn in adcp_fns:
    print(f'reading {fn}')
    A[fn]=sp.rdi_readin_adcp_VariableBins(os.path.join(pth,fn))
    
#%%
importlib.reload(sp)

profiles_WB = sp.attach_adcp_to_profiles(profiles_WB,A,buffer)

sp.profile_map_figure(profiles_WB, A)

       

# %% 
# calculate u*

profiles_WB = sp.vel_profiles(profiles_WB)

# %%
# make figures

# %% this will read west bay ADCP data
pth = r'C:\Users\cesposito\THE WATER INSTITUTE OF THE GULF\P-00703_NSF_Caltech - General\Data\West Bay Report\ADCP Data\Trip 1c 5-6 May 2009\ASCII'
fn_adcp = os.path.join(pth,'WBAY_166t.000')
import sedProfiles as sp
B = sp.rdi_readin_adcp_VariableBins(fn_adcp)

station='R-4.5-A'


#%%

# Reproject to Web Mercator (required for contextily basemaps)
loc_web = loc.to_crs(epsg=3857)
lw = loc_web[loc_web['station']==station]
lw = lw.drop_duplicates()

buffer = 100    #meters
lwb = lw.buffer(buffer)
lwb_gdf = gpd.GeoDataFrame(geometry=lwb,crs = lw.crs)
clip_polygon = lwb_gdf.geometry.iloc[0]


points = [Point(xy) for xy in zip(B['pos']['lon'], B['pos']['lat'])]  # x=lon, y=lat

# Create GeoDataFrame and set WGS84 CRS
points_gdf = gpd.GeoDataFrame(index=range(len(points)), geometry=points, crs='EPSG:4326')


# Reproject to match loc_web (EPSG:3857)
points_web = points_gdf.to_crs(epsg=3857)


# Boolean mask: which points are inside the polygon
inside_mask = points_web.geometry.within(clip_polygon)

# Get indices or values
inside_indices = points_web[inside_mask].index.to_list()


# Plot
fig, ax = plt.subplots()
lwb_gdf.plot(ax=ax,facecolor='none',edgecolor='purple')
loc_web.plot(ax=ax, markersize=10, color='red')
points_web.iloc[inside_indices].plot(ax=ax,color='blue',linewidth=2)
lw.plot(ax=ax,color='purple')

# Add basemap
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=12)

plt.show()
