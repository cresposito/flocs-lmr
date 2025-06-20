import dolfyn
import numpy as np
from dolfyn.adp import api
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as dt
import os
import pandas as pd
from sklearn import linear_model
from datetime import datetime
from collections import defaultdict


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

def rdi_readin_adcp_VariableBins(fn_adcp):
      
    xx = {}
    xx['top'] = {}
        
    with open(fn_adcp,'r') as fid:
        tt = rdi_topheader(fid)
        for k,v in tt.items():
            xx['top'][k] = v

        first_ensemble = True
        
        ctr = 0
        while True:
            ctr = ctr+1
            tline = fid.readline()
            if not tline: break
        
            aa = rdi_header(fid,tline)
            bb = rdi_read_adcp_ens(fid, aa['num_bins'])
                    
            if first_ensemble:
                recursive_init(xx, aa)                 # for nested scalar/list data
                recursive_init(xx, bb, stack_arrays=True)  # for nested arrays (e.g. bin data)
                first_ensemble = False
            else:
                recursive_append(xx, aa)
                recursive_append(xx, bb, stack_arrays=True)
       
    return xx

def recursive_init(xx, data, stack_arrays=False):
    for k, v in data.items():
        if isinstance(v, dict):
            xx[k] = {}
            recursive_init(xx[k], v, stack_arrays)
        elif isinstance(v, np.ndarray) and stack_arrays:
            # Initialize with transposed array (assumes 2D)
            xx[k] = v.T
        else:
            xx[k] = [v]  # Wrap scalar/list for append

def recursive_append(xx, data, stack_arrays=False):
    for k, v in data.items():
        if isinstance(v, dict):
            if k not in xx:
                xx[k] = {}
            recursive_append(xx[k], v, stack_arrays)
        elif isinstance(v, np.ndarray) and stack_arrays:
            # Stack arrays along the 2nd axis (columns = time/ensembles)
            xx[k] = np.column_stack([xx[k], v.T])
        else:
            xx[k].append(v)
            
def rdi_read_adcp_ens(fid, num_bins):
    # Read all lines at once
    lines = [fid.readline().strip().split() for _ in range(num_bins)]
    
    # Convert to float array
    bin_array = np.array(lines, dtype=float)  # shape (num_bins, 13)
    
    # clean for bad values here
    #####################
    #####################
   
    out = {}
    out['vel'] = {}
    out['bks'] = {}
    
    out['bin_dpt'] = bin_array[:, 0]
    out['bin_dsc'] = bin_array[:, 12]
    
    out['vel']['vel_mag'] = bin_array[:, 1]
    out['vel']['dir'] = bin_array[:, 2]
    out['vel']['evel'] = bin_array[:, 3]
    out['vel']['nvel'] = bin_array[:, 4]
    out['vel']['zvel'] = bin_array[:, 5]
    out['vel']['ervel'] = bin_array[:, 6]
    
    out['bks']['b1'] = bin_array[:, 7]
    out['bks']['b2'] = bin_array[:, 8]
    out['bks']['b3'] = bin_array[:, 9]
    out['bks']['b4'] = bin_array[:, 10]

    return out 
    

def rdi_topheader(fid):
    out = {}

    # First two lines are notes
    out['note1'] = fid.readline().strip()
    out['note2'] = fid.readline().strip()

    # Third line: read 7 float values
    line = fid.readline()
    parts = line.strip().split()

    out['dcl'] = int(parts[0])  #depth cell length (cm)
    out['bat'] = int(parts[1])  #blank after transmit (cm)
    out['dcn'] = int(parts[2])  #adcp depth cfom config. node (cm)
    out['ndc'] = int(parts[3])  #num depth cells
    out['ppe'] = int(parts[4])  #pings per ensemble 
    out['tpe'] = int(parts[5])  #time per ensemble (hundredths of seconds)
    out['prm'] = int(parts[6])  #profiling mode

    return out

def rdi_header(fid,tline):
    out = {}
    out['dpth']={}
    out['pos']={}
    out['disc']={}
    out['misc']={}
    out['bt_vel']={}
    out['dist']={}
    
    # Line 1: passed as tline (not read from file)
    aa = [float(x) for x in tline.strip().split()]
    year = aa[0] + 2000 if aa[0] < 100 else aa[0]  # Adjust 2-digit year
    dt = datetime(int(year), int(aa[1]), int(aa[2]), int(aa[3]), int(aa[4]), int(aa[5]),int(aa[6]*10)) #year, month, day, hours, minutes, seconds, microseconds (convert from 1/100 seconds as given in ascii output)
    
    
    # clean for bad values here
    #####################
    #####################
    
    
    
    #cast into dictionary here
    
    out['time'] = dt              # year, month, day, hour, minute, second, hundredths
    out['pitrl'] = aa[9:11]            # pitch & roll (degrees)
    
    out['misc']['corhead'] = aa[11]            # average adcp heading
    out['misc']['temp'] = aa[12]               # instrument temperature (degrees C)

    # Line 2
    tline = fid.readline()
    bb = [float(x) for x in tline.strip().split()]

    out['bt_vel']['bt_evel'] = bb[0]                # bottom tracking: east, north, vertical, error
    out['bt_vel']['bt_nvel'] = bb[1]                # bottom tracking: east, north, vertical, error
    out['bt_vel']['bt_zvel'] = bb[2]                # bottom tracking: east, north, vertical, error
    out['bt_vel']['bt_etvel'] = bb[3]                # bottom tracking: east, north, vertical, error
    
    
    out['dpth']['b1'] = bb[8]             # depth beams 1
    out['dpth']['b2'] = bb[9]             # depth beams 2
    out['dpth']['b3'] = bb[10]             # depth beams 3
    out['dpth']['b4'] = bb[11]             # depth beams 4
    out['dpth']['mean_depth'] = np.nanmean([bb[8], bb[9], bb[10], bb[11]])

    # Line 3
    tline = fid.readline()
    cc = [float(x) for x in tline.strip().split()]
    out['dist']['total'] = cc[0]            # total distance
    out['dist']['time_s'] = cc[1]            # elapsed time
    out['dist']['north'] = cc[2]            # distance north
    out['dist']['east'] = cc[3]            # distance east
    out['dist']['made_good'] = cc[4]            # distance made good

    # Line 4
    tline = fid.readline()
    dd = [float(x) for x in tline.strip().split()]
    out['pos']['lat'] = dd[0]            #Latitude
    out['pos']['lon'] = dd[1]            #Longitude
    out['pos']['Evel'] = dd[2]            #gga or vtg East velocity in m/s or ft/s
    out['pos']['Nvel'] = dd[3]            #gga or vtg North velocity in m/s or ft/s

    # Line 5
    tline = fid.readline()
    ee = [float(x) for x in tline.strip().split()]
    out['disc']['Q_middle'] = ee[0]
    out['disc']['Q_top'] = ee[1]
    out['disc']['Q_bot'] = ee[2]
    out['disc']['Q_start_shore']=ee[3]
    out['disc']['Q_end_shore']=ee[5]
    out['disc']['dist_start_shore']=ee[4]
    out['disc']['dist_end_shore']=ee[6]
    out['disc']['dpt_top_middle']=ee[7]
    out['disc']['dpt_bot_middle']=ee[8]

    # Line 6
    tline = fid.readline()
    ff = tline.strip().split()
    out['num_bins'] = int(ff[0])                # Only the first value is used for number of bins
    out['unit'] = ff[1]
    
    out['misc']['vel_ref'] = ff[2]
    out['misc']['intensity_unit'] = ff[3]
    out['misc']['intensity_scale_factor'] = float(ff[4])
    out['misc']['sound_absorption_factor'] = float(ff[5])

    return out