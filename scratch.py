import os

import numpy as np
import pandas as pd

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
profiles = sp.profiles_from_sampling_info(smp,loc)

# add grain size data
profiles = sp.profiles_add_gs(profiles,df_gs)

# add adcp data
# ens_start = profiles[1]['adcp']['ADCP Ensemble Start'][0]
# ens_end = profiles[1]['adcp']['ADCP Ensemble End'][0]
# fn = sp.find_file(adcp_folder,profiles[1]['adcp']['Stationary ADCP File Name'][0])
# V_mag, t = sp.read_adcp_section(fn,ens_start,ens_end,profiles[1]['Total Depth (m)'])

profiles = sp.read_adcp_sections(profiles,adcp_folder)

# %% 
# calculate u*

# %%
# make figures

# %% 
# add other parameters
        
        
    


