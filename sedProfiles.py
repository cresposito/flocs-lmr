import dolfyn
import numpy as np
from dolfyn.adp import api
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as dt
import os
import pandas as pd
from sklearn import linear_model

def read_sed_samples(filename, sheet, parameters_to_include):
                     
# =============================================================================
    
    # reads an excel sheet with individual sediment grain size measurements.
    # sheet should be patterned on the MR-173 format

    # FILENAME:  string, path of the excel file
    # SHEET:     string, name of the sheet to read
    # PARAMETERS_TO_INCLUDE   list of strings, indicating the parameters from the headers to include. All other parameters in thea headers will be ignored
# =============================================================================
    
    # Read the Excel file
    df = pd.read_excel(filename, header=None)
    
    #this nifty line removes trailing :'s from any strings in the first column
    df[0] = df[0].apply(lambda x: x.rstrip(':') if isinstance(x, str) else x)
    
    # The first column contains parameter names
    parameter_names = df.iloc[:, 0]
    # parameter_names_str = parameter_names.astype(str)
    # parameter_names_str = parameter_names_str.str.replace(":", "", regex=False)
    # include_mask = parameter_names_str.isin(parameters_to_include)
    include_mask = parameter_names.isin(parameters_to_include)


    # Filter the parameter names
    filtered_parameters = parameter_names[include_mask].values
    
    #identify the location of the grain size data
    idx = df.index[df.iloc[:, 0] == 'Channel Diameter (Lower)']
    start_pos = df.index.get_loc(idx[0]+3)
    bin_lower = df.iloc[start_pos:, 0].values
    units = str(df.iloc[idx[0]+1,0])
    
    # Create one dictionary per data column (excluding the first)
    dicts = {}
    for col in df.columns[1:]:
        data = df.iloc[:, col][include_mask].values
        col_dict = dict(zip(filtered_parameters, data))
        col_dict['units'] = units
        col_dict['bin lower'] = bin_lower
        col_dict['gs'] = df.iloc[start_pos:, col].values
        dicts[col] = col_dict
    
    df_gs = pd.DataFrame.from_dict(dicts,orient='index')
    return df_gs

# def get_dir_size(path='.'):
#     ...:     total = 0
#     ...:     with os.scandir(path) as it:
#     ...:         for entry in it:
#     ...:             if entry.is_file():
#     ...:                 total += entry.stat().st_size
#     ...:             elif entry.is_dir():
#     ...:                 total += get_dir_size(entry.path)
    # ...:     return total):

def profiles_from_sampling_info(smp,loc):
    # organize data into profiles
    
    #SMP: dataframe with sampling information
    #LOC: dataframe with station location information

    #profiles will be defined by a unique combination of Station and Date
    unique_pairs = smp[['Station', 'Date']].drop_duplicates()
    unique_list = list(unique_pairs.itertuples(index=False, name=None))

    profiles={}
    for ctr, st_da in enumerate(unique_list):
        profiles[ctr]={} #create nested dictionary    
        
        #all rows in the sampling dataframe with this station and date. 
        #These should constitute a distinct profiles
        smp_prf = smp[(smp['Station'] == st_da[0]) & (smp['Date'] == st_da[1])].reset_index(drop=True)

        #pull data that applies to the entire profile from the first row
        profiles[ctr]['Station'] = smp_prf.loc[0]['Station']
        profiles[ctr]['Date'] = smp_prf.loc[0]['Date']
        profiles[ctr]['Total Depth (ft)'] = smp_prf.loc[0]['Total Depth (ft)']
        profiles[ctr]['Total Depth (m)'] = profiles[ctr]['Total Depth (ft)']*0.3048
        profiles[ctr]['Lat'] = loc[loc['station']==profiles[ctr]['Station']]['lat'].values    #for lat and lon, match the station name
        profiles[ctr]['Lon'] = loc[loc['station']==profiles[ctr]['Station']]['lon'].values

        profiles[ctr]['gs']={}
        profiles[ctr]['adcp']={}
        profiles[ctr]['bed']={}

        #all other data is obtained from each sample
        profiles[ctr]['gs']['Isokinetic Sample Depth'] = smp_prf['Isokinetic Sample Depth'].values
        profiles[ctr]['gs']['Depth Increment'] = smp_prf['Depth Increment'].values
        profiles[ctr]['gs']['Replicate Indicator'] = smp_prf['Replicate Indicator'].values
        profiles[ctr]['gs']['Date & Time'] = smp_prf['Date & Time'].values
        
        profiles[ctr]['adcp']['ADCP Ensemble Start'] = smp_prf['ADCP Ensemble Start'].values
        profiles[ctr]['adcp']['ADCP Ensemble End'] = smp_prf['ADCP Ensemble End'].values
        profiles[ctr]['adcp']['Stationary ADCP File Name'] = smp_prf['Stationary ADCP File Name'].values
        profiles[ctr]['adcp']['Date & Time'] = smp_prf['Date & Time'].values

    return profiles

def profiles_add_gs(profiles,df_gs):
    # adds grain size data to the profiles
    
    #PROFILES: dictionary of profiles, as read in by profiles_from_sampling_info
    #DF_GS:    dataframe with grain size measurements, as read in by read_sed_samples
    
    for key in profiles.keys():
        gs_rows = []    #this will store grain size data
        for ii, di in enumerate(profiles[key]['gs']['Depth Increment']):
            #create the sample name that will be used to identify this sample
            ss = profiles[key]['Station'] + '_' + str(di)
            if pd.notna(profiles[key]['gs']['Replicate Indicator'][ii]):
                ss = ss + profiles[key]['gs']['Replicate Indicator'][ii]
                
            aa = df_gs[df_gs['Sample ID']==ss]
            
            #ensure that aa has exactly one row
            #this means that the expected sample name is found exactly once in the 
            #detailed grain size data
            if len(aa)==1:
                print(ss+ ' FOUND')
                
                #simply store the dataframe with grain size data
                
                #profiles[key]['gs']['gs_'+str(ii)] = aa
                profiles[key]['gs']['Sample ID'] = aa.iloc[0]['Sample ID']
                profiles[key]['gs']['units'] = aa.iloc[0]['units']            
                profiles[key]['gs']['bin lower'] = aa.iloc[0]['bin lower']
                
                gs_rows.append(aa.iloc[0]['gs'])           
                
            else:
                print(ss+ ' NOT FOUND')
                
        profiles[key]['gs']['gs'] = np.vstack(gs_rows)
        
    return profiles

def read_adcp_sections(profiles,adcp_folder):
    
    for key in profiles.keys():
        print('------------------')
        print(profiles[key]['Station'])
        
        profiles[key]['adcp']['t']=[]
        profiles[key]['adcp']['z']=[]
        profiles[key]['adcp']['V_mag']=[]
        
        for i in range(len(profiles[key]['adcp']['Stationary ADCP File Name'])):      
            adcp_path = find_file(adcp_folder,profiles[key]['adcp']['Stationary ADCP File Name'][i])
            ens_start = profiles[key]['adcp']['ADCP Ensemble Start'][i]
            ens_end = profiles[key]['adcp']['ADCP Ensemble End'][i]
            
            try:
                t,z,V_mag = read_adcp_section(adcp_path,ens_start,ens_end,profiles[key]['Total Depth (m)'])
                
                profiles[key]['adcp']['t'].append(t)
                profiles[key]['adcp']['z'].append(z)
                profiles[key]['adcp']['V_mag'].append(V_mag)
                
                print('.....adcp data successfully read')
                
            except Exception as e:
                
                print('......failure reading adcp data')
                                    
            
    return profiles


            
def read_adcp_section(adcp_path,ens_start,ens_end,station_depth):
    
    ds=dolfyn.read(adcp_path)
    
    i1 = np.where(ds.number.values == ens_start)[0].item()
    i2 = np.where(ds.number.values == ens_end)[0].item()
    
    t = dolfyn.time.dt642date(ds.time)
    
    # bin centers
    bin_centers=ds["range"]
    
    # Calculate the mean distance to the bottom for each time
    mean_dist = ds["dist_bt"].mean(dim=["beam"])
    
    # Create a boolean mask where the condition is true
    mask = mean_dist > bin_centers
    
    V_masked=ds["vel"].where(mask, np.nan)
    V_mag=np.linalg.norm(V_masked,axis=0)
    
    #height above bed
    z=station_depth-ds['range']
    
    # clip for ensembles
    V_mag = V_mag[:, i1:i2]
    t = t[i1:i2]
    
    return t,z,V_mag
    
def find_file(start_dir, filename):
    #find a file named filename. Assumes that there is only one file within start_dir with this name.
    # Iterate through all files and directories in the start directory
    for root, dirs, files in os.walk(start_dir):
        # Check if the file we are looking for exists in the current directory
        if filename in files:
            # Return the full path to the file
            return os.path.join(root, filename)
    
    # If the file is not found in any subdirectory
    return None