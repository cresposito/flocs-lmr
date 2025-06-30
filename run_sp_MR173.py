import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as dt
import geopandas as gpd
import contextily as cx
from shapely.geometry import Point
import pickle


import sedProfiles as sp

# %%
baseDataSetFolder=r'C:\Users\cesposito\The Water Institute of the Gulf\P-00703_NSF_Caltech - General\Data\MR173 Head of Passes'
adcp_folder = os.path.join(baseDataSetFolder,'StationaryADCP')

# filename, and sheet for the detailed grain size data.
fn_gs = os.path.join(baseDataSetFolder,r'Isokinetic Suspended Sediment\MR173_Grain Size Compilation.xlsx')
sh_gs = 'MR173'

#parameters to include from grain size file header
params_fn_gs = ['File ID','Sample ID']     

# filename and sheet for station locations
fn_loc = os.path.join(baseDataSetFolder,r'Isokinetic Suspended Sediment\MR173_Suspended Sediment Concentrations.xlsx')
sh_loc = 'IsoSamplingLocs'

# filename and sheet with detailed information about all of the sampling (ADCP, and isokinetic)
fn_smp = os.path.join(baseDataSetFolder,r'StationaryADCP\MR173_StationaryADCP&Isokinetic.xlsx')
sh_smp = 'ADCP Vertical '


# %%
# read the detailed grain size data
df_gs = sp.read_sed_samples(fn_gs,sh_gs,params_fn_gs)

# read location information
headers = ['station','lat','lon']
loc = pd.read_excel(fn_loc,sheet_name=sh_loc, header = None, names = headers)

# read sampling information
smp = pd.read_excel(fn_smp,sheet_name = sh_smp)


# %%
# organize data into profiles
profiles_MR173 = sp.profiles_from_sampling_info(smp,loc)

# add grain size data
profiles_MR173 = sp.profiles_add_gs(profiles_MR173,df_gs,sample_id_style=1)

# add adcp data
# ens_start = profiles[1]['adcp']['ADCP Ensemble Start'][0]
# ens_end = profiles[1]['adcp']['ADCP Ensemble End'][0]
# fn = sp.find_file(adcp_folder,profiles[1]['adcp']['Stationary ADCP File Name'][0])
# V_mag, t = sp.read_adcp_section(fn,ens_start,ens_end,profiles[1]['Total Depth (m)'])

profiles_MR173 = sp.read_adcp_sections(profiles_MR173,adcp_folder)
profiles_MR173 = sp.vel_profiles(profiles_MR173)

# fn_MR173_pkl = r'C:\TEMP\MR173.pkl'
# with open(fn_MR173_pkl, 'wb') as f:
#     pickle.dump(profiles_MR173,f)

# %% map figure of station locations

fig, ax = plt.subplots()

Lat = []
Lon = []
Station = []
for key in profiles_MR173.keys():
    Lat.append(profiles_MR173[key]['Lat'])
    Lon.append(profiles_MR173[key]['Lon'])
    Station.append(profiles_MR173[key]['Station'])

# Create geometry from lon/lat
geometry = [Point(lon, lat) for lon, lat in zip(Lon, Lat)]

# Create GeoDataFrame with just Station and geometry
profile_gdf = gpd.GeoDataFrame({'Station': Station}, geometry=geometry, crs='EPSG:4326')
profile_gdf = profile_gdf.to_crs(epsg=3857)

profile_gdf.plot(ax=ax)
for x, y, label in zip(profile_gdf.geometry.x, profile_gdf.geometry.y, profile_gdf['Station']):
    ax.text(x, y, label, fontsize=8, ha='left', va='bottom')

# Add basemap
cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=12)

# plt.show()

#%% velocity profiles for all profiles
import sedProfiles as sp

for key in profiles_MR173.keys():
    
    fn_save = profiles_MR173[key]['Station']+'_'+str(profiles_MR173[key]['Date'])+'.png'
    illegal_chars = {'<', '>', ':', '"', '/', '\\', '|', '?', '*'}
    fn_save = ''.join(c if c not in illegal_chars else '_' for c in fn_save)
    fn_save = fn_save+'.png'
    
    pth = r'C:\Users\cesposito\THE WATER INSTITUTE OF THE GULF\P-00703_NSF_Caltech - General\Data\_ProcessedData'
    
    fig, ax = sp.velocity_profile_figure(profiles_MR173[key],fn_save=os.path.join(pth,fn_save))
    plt.close(fig)

#%%

#%% grain size profiles