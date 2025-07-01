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

# data set
data_set = 'WB'


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

# add location
for key in profiles_WB.keys():
    profiles_WB[key]['data_set'] = data_set

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
      

# %% 
# calculate u*

profiles_WB = sp.vel_profiles(profiles_WB)

# %% map figure

pth = r'C:\Users\cesposito\THE WATER INSTITUTE OF THE GULF\P-00703_NSF_Caltech - General\Data\_ProcessedData'
fn_save = 'WB_Map.png'
sp.profile_map_figure(profiles_WB, A,fn_save = os.path.join(pth,fn_save))


#%% velocity profiles for all profiles

for key in profiles_WB.keys():
    
    fn_save = (
        'VEL_'
        +profiles_WB[key]['data_set']
        +'_'
        +profiles_WB[key]['Station']
        +'_'
        +str(profiles_WB[key]['Date'])
        +'.png'
        )
        
    illegal_chars = {'<', '>', ':', '"', '/', '\\', '|', '?', '*'}
    fn_save = ''.join(c if c not in illegal_chars else '_' for c in fn_save)
    
    pth = r'C:\Users\cesposito\THE WATER INSTITUTE OF THE GULF\P-00703_NSF_Caltech - General\Data\_ProcessedData'
    
    profile_fig_out = sp.profile_figure(profiles_WB[key],fn_save=os.path.join(pth,fn_save))
    if profile_fig_out == None:
        pass
    else:
        fig, ax = profile_fig_out   
        plt.close(fig)

