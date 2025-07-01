import dolfyn
import numpy as np
from dolfyn.adp import api
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as dt
import os
import pandas as pd
import geopandas as gpd
import contextily as cx
from sklearn import linear_model
from datetime import datetime
from collections import defaultdict
import pdb
from shapely.geometry import Point
from matplotlib.gridspec import GridSpec



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
    unique_pairs = unique_pairs.dropna()
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

def profiles_add_gs(profiles,df_gs,sample_id_style):
    # adds grain size data to the profiles
    
    #PROFILES: dictionary of profiles, as read in by profiles_from_sampling_info
    #DF_GS:    dataframe with grain size measurements, as read in by read_sed_samples
    
    for key in profiles.keys():
        gs_rows = []    #this will store grain size data
        
        #create the sample names that will be used to find the samples in this profile in the grain size data spreadsheet
        if sample_id_style ==1:
            #used for MR173
            station = profiles[key]['Station']
            dpth_inc = profiles[key]['gs']['Depth Increment']
            sample_id = [station + '_' + str(di) for di in dpth_inc]
            
        elif sample_id_style ==2:
            #used for West Bay
            arr = profiles[key]['gs']['Isokinetic Sample Depth']
            dpt = [f'{x:.1f}' if x % 1 else f'{int(x)}' for x in arr]

            sample_id = []            
            for dp in dpt:
                sample_id.append(profiles[key]['Station'] + ' ' + str(dp) + 'ft')        
        else:
            raise ValueError('invalid sample_id_style')
            
        #replicate indicator style should be the same for all input data sets
        for ii, ri in enumerate(profiles[key]['gs']['Replicate Indicator']):
            #the replicate indicator will be a string if it is defined
            if type(ri)==str:
                sample_id[ii]=sample_id[ii]+ri

        
        for ii, sid in enumerate(sample_id):
            # ss = profiles[key]['Station'] + '_' + str(sid)
            # if pd.notna(profiles[key]['gs']['Replicate Indicator'][ii]):
                # ss = ss + profiles[key]['gs']['Replicate Indicator'][ii]
                           
            aa = df_gs[df_gs['Sample ID']==sid]
            
            #ensure that aa has exactly one row
            #this means that the expected sample name is found exactly once in the 
            #detailed grain size data
            if len(aa)==1:
                print(sid+ ' FOUND')
                
                #simply store the dataframe with grain size data
                
                #profiles[key]['gs']['gs_'+str(ii)] = aa
                profiles[key]['gs']['Sample ID'] = aa.iloc[0]['Sample ID']
                profiles[key]['gs']['units'] = aa.iloc[0]['units']            
                profiles[key]['gs']['bin lower'] = aa.iloc[0]['bin lower']
                
                gs_rows.append(aa.iloc[0]['gs'])           
                
            else:
                print(sid+ ' NOT FOUND')
        
        profiles[key]['gs']['gs'] = np.vstack(gs_rows)
        
    return profiles

def read_adcp_sections(profiles,adcp_folder):
    
    for key in profiles.keys():
        print('------------------')
        print(profiles[key]['Station'])
        
        profiles[key]['adcp']['t']=[]
        profiles[key]['adcp']['bin_dpth']=[]
        profiles[key]['adcp']['bed_dpth']=[]
        profiles[key]['adcp']['z']=[]
        profiles[key]['adcp']['V_mag']=[]
        
        for i in range(len(profiles[key]['adcp']['Stationary ADCP File Name'])):      
            adcp_path = find_file(adcp_folder,profiles[key]['adcp']['Stationary ADCP File Name'][i])
            ens_start = profiles[key]['adcp']['ADCP Ensemble Start'][i]
            ens_end = profiles[key]['adcp']['ADCP Ensemble End'][i]
            
            try:
                t,bin_dpth, bed_dpth,z,V_mag = read_adcp_section(adcp_path,ens_start,ens_end)
                
                profiles[key]['adcp']['t'].append(t)
                profiles[key]['adcp']['bin_dpth'].append(bin_dpth)
                profiles[key]['adcp']['bed_dpth'].append(bed_dpth)
                profiles[key]['adcp']['z'].append(z)
                profiles[key]['adcp']['V_mag'].append(V_mag)
                
                print('.....adcp data successfully read')
                
            except Exception as e:
                
                print('......failure reading adcp data')
                                    
            
    return profiles


            
def read_adcp_section(adcp_path,ens_start,ens_end):
    
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
    
    #bin depths
    bin_dpth = ds['range']
    
    #bed depths
    bed_dpth = ds['dist_bt'].mean(axis=0).values[i1:i2]
    
    #height above bed
    dd=[]
    for d in bed_dpth:
        dd.append(d-bin_dpth)
    
    z=np.stack(dd,axis=1)
    
    # clip for ensembles
    V_mag = V_mag[:, i1:i2]
    t = t[i1:i2]
    
    return t,bin_dpth, bed_dpth,z,V_mag
    
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

def rdi_readin_adcp_VariableBins(fn_adcp,points_gdf = True):
    #fn_adcp:       filename with WinRiver classic ASCII output
    #points_gdf:    if true, include a geodataframe in the position field
      
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
            
            ###make a z field
            bb['z'] = aa['dpth']['mean_depth']-bb['bin_dpth']  #height above the bed added to header
                    
            if first_ensemble:
                recursive_init(xx, aa)                 # for nested scalar/list data
                recursive_init(xx, bb, stack_arrays=True)  # for nested arrays (e.g. bin data)
                first_ensemble = False
            else:
                recursive_append(xx, aa)
                recursive_append(xx, bb, stack_arrays=True)
    
    if points_gdf:
        points = [Point(xy) for xy in zip(xx['pos']['lon'], xx['pos']['lat'])]  # x=lon, y=lat
        
        # Create GeoDataFrame and set WGS84 CRS
        points_gdf = gpd.GeoDataFrame(index=range(len(points)), geometry=points, crs='EPSG:4326')

        #points stored in both wgs84(lat/lon) and in web mercator
        xx['pos']['points_wgs84']=points_gdf
        xx['pos']['points_web']=points_gdf.to_crs(epsg=3857)        
        
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
    
    # clean for bad values here
    #####################
    #####################

    return out

            
def rdi_read_adcp_ens(fid, num_bins):
    # Read all lines at once
    lines = [fid.readline().strip().split() for _ in range(num_bins)]
    
    # Convert to float array
    bin_array = np.array(lines, dtype=float)  # shape (num_bins, 13)
      
    out = {}
    out['vel'] = {}
    out['bks'] = {}
    
    out['bin_dpth'] = bin_array[:, 0]
    out['bin_disc'] = bin_array[:, 12]
    
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

    ####
    # clean for bad values here
    out['bin_disc'][out['bin_disc']==2147483647]=np.nan
    
    ii_bad_vel = out['vel']['vel_mag']==-32768

    for key in out['vel']:
        out['vel'][key][ii_bad_vel] = np.nan

    for key in out['bks']:
        out['bks'][key][ii_bad_vel] = np.nan
        

    return out 

def ft_to_SI(xx):
    
    if not all(u == 'ft' for u in xx['unit']):
        raise ValueError('all these should be ft...add other functionality as needed')
        
    
        

def attach_adcp_to_profiles(profiles,A,buffer):
    #profiles:      grain size profiles
    #A              dictionary of ADCP profiles as read in by rdi_readin_adcp_VariableBins
    #buffer         buffer to associate ADCP transect data with profile station points
    
    for key in profiles.keys():
        print(key)
        prf = profiles[key]
        
        #compute clipping polygon
        pt=Point(prf['Lon'],prf['Lat'])
        gdf = gpd.GeoDataFrame(geometry=[pt],crs='EPSG:4326')
        gdf = gdf.to_crs(epsg=3857)
        gdf_buff = gpd.GeoDataFrame(geometry=gdf.buffer(buffer),crs = gdf.crs)
        
        #get list of all adcp transects associated with this profile
        adcp_fns = []
        for lne in prf['adcp']['Stationary ADCP File Name']:
            #list of adcp file names. if there is a file with an r suffix, change it to t
            adcp_fns = adcp_fns + lne.replace(' ','').replace('r.','t.').split(';')  
            
        adcp_fns = list(set(adcp_fns))    

        #for each adcp transect in the adcp key
        prf['adcp']['V_mag']=[]
        prf['adcp']['Stationary ADCP File Name']=[]
        prf['adcp']['t']=[]
        prf['adcp']['z']=[]
        prf['adcp']['bed_dpth']=[]
        prf['adcp']['bin_dpth']=[]
        prf['adcp']['pos']=[]
        prf['station_pt']=gdf
        prf['buffer_poly']=gdf_buff
        
        for fn in adcp_fns:
            # print(fn)
            inside_mask = A[fn]['pos']['points_web'].geometry.within(gdf_buff.geometry.iloc[0])
            inside_indices = A[fn]['pos']['points_web'][inside_mask].index.to_list()
            # print(inside_indices)
            
            prf['adcp']['V_mag'].append(A[fn]['vel']['vel_mag'][:,inside_indices])
            prf['adcp']['z'].append(A[fn]['z'][:,inside_indices])
            prf['adcp']['Stationary ADCP File Name'].append(fn)
            prf['adcp']['t'].append(np.array(A[fn]['time'])[inside_indices])
            prf['adcp']['pos'].append(A[fn]['pos']['points_web'].iloc[inside_indices])
            prf['adcp']['bed_dpth'].append(np.array(A[fn]['dpth']['mean_depth'])[inside_indices])

                        
            #check that all columns are the same. if they are, set the bin depth column and the z column. if not, throw an error
            aa = A[fn]['bin_dpth']
            if np.all(aa == aa[:,0][:,None]):
                prf['adcp']['bin_dpth'].append(aa[:,0])
            else:
                raise ValueError('code does not handle variable bins')
                
        profiles[key]=prf
        
    return profiles


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




def vel_profiles(profiles):
    
    #PROFILES:  dict of sediment grain size profiles data

    for k1 in profiles.keys():
        ustar_list =[]
        z0_list = []
        ufit_list = []
        zfit_list =[]
        UU = []
        ZZ = []
        
        for i in range(len(profiles[k1]['adcp']['V_mag'])):
            if profiles[k1]['adcp']['V_mag'][i].shape[1]>0:
                ustar, z0, UU, ZZ = vel_profile(
                    V_mag = profiles[k1]['adcp']['V_mag'][i],
                    z = profiles[k1]['adcp']['z'][i]
                    )
    
                k=0.41            
                zfit = np.arange(0.1,profiles[k1]['adcp']['z'][i].max(),0.1)
                ufit = (ustar/k)*np.log(zfit/z0)            
                
                ustar_list.append(ustar)
                z0_list.append(z0)
                ufit_list.append(ufit)
                zfit_list.append(zfit)

        if len(ufit_list)>0:        #profiles with not valid fits will not get the ['fit'] field
            profiles[k1]['adcp']['fit']={}
            profiles[k1]['adcp']['fit']['ustar'] = ustar_list
            profiles[k1]['adcp']['fit']['z0'] = z0_list
            profiles[k1]['adcp']['fit']['ufit'] = ufit_list
            profiles[k1]['adcp']['fit']['zfit'] = zfit_list
            profiles[k1]['adcp']['fit']['UU'] = UU
            profiles[k1]['adcp']['fit']['ZZ'] = ZZ
        
    return profiles

        
def vel_profile(V_mag,z):
    # computes the u* and z0 for a section of adcp data
    # V_MAG     mxn numpy array of velocity data. 
    # Z         (m,) data array of z levels (height above bed) where the velocity data were collected
    
    ii = ~np.isnan(V_mag)
        
    UU = V_mag[ii]
    ZZ = z[ii]
    
    #keep only elements where both UU and ZZ are positive
    jj = (UU >= 0) & (ZZ >= 0)
    
    UU=UU[jj]
    ZZ=ZZ[jj]
    
    
    lnZ=np.log(ZZ)
       
    lnZ=lnZ.reshape(-1,1)
    regr=linear_model.LinearRegression()
    
    try:
        regr.fit(lnZ,UU)
        k=0.41
        ustar=regr.coef_*k
        z0=np.exp(-1*k*regr.intercept_/ustar)

    except ValueError as e:
        # Handle the ValueError (or any exception) here
        print(f"An error occurred: {e}")
        pdb.set_trace()    
    
    return ustar, z0, UU, ZZ

def profile_figure(profile,fn_save = None):

    station = profile['Station']
    date =  str(profile['Date'])
    dpth = profile['Total Depth (m)']
    adcp = profile['adcp']
    
    if 'fit' in profile['adcp']:
        adcp_fit = profile['adcp']['fit']
    else:
        print(f'no adcp fits in profile {station}')
        return
    
    num_fits = len(adcp['fit']['ufit'])
    num_rows = num_fits+1   #the last row is for the grain size profile
    
    # fig, axs = plt.subplots(num_rows,2,width_ratios=[3,1],figsize=(10,2+2*num_rows))
    
    fig = plt.figure(layout = 'constrained',figsize=(10,2+2*num_rows))
    fig.suptitle(f'{station} \n {date} \n Total Depth (m): {dpth:.2f}')
    
    gspec = GridSpec(num_rows,3,figure=fig)
    axs = np.empty((num_rows,3),dtype=object)
    
    for i in range(num_fits):
    
        t = adcp['t'][i]
        bin_dpth = adcp['bin_dpth'][i]
        bed_dpth = adcp['bed_dpth'][i]
        V_mag = adcp['V_mag'][i]

        ufit_list = adcp_fit['ufit'][i]
        zfit_list = adcp_fit['zfit'][i]
        ustar = adcp_fit['ustar'][i]
        z0 = adcp_fit['z0'][i]
        
        axs[i,0] = fig.add_subplot(gspec[i,:-1])
        axs[i,1] = fig.add_subplot(gspec[i,-1])
        
        axs[i,0].pcolormesh(
            t,
            bin_dpth,
            V_mag,
            cmap='Blues',
            shading='nearest'
            )
        
        axs[i,0].plot(t, bed_dpth, linestyle=':', linewidth=2, color='black')
        axs[i,0].yaxis.set_inverted(True)
        axs[i,0].set_ylim(bottom=bed_dpth.max()*1.15)
        axs[i,0].xaxis.set_major_formatter(dt.DateFormatter('%H:%M:%S'))
        axs[i,0].xaxis.set_major_locator(dt.SecondLocator(bysecond=range(0, 60, 15)))  # Tick every 15 seconds starting on the minute
        
        axs[i,0].set_ylabel('depth (m)')

        axs[i,1].plot(adcp_fit['UU'],adcp_fit['ZZ'],'.',color='lightgrey',zorder=0)
        axs[i,1].plot(ufit_list,zfit_list,label = f'u* = {ustar[0]:.2f} \nz0 = {z0[0]:.2f}')
        axs[i,1].plot([0,2], [0, 0], linestyle=':', linewidth=2, color='black')
        axs[i,1].set_ylim(bottom=-1)
        axs[i,1].set_xlim([0,2])
        axs[i,1].legend()
        axs[i,1].set_ylabel('height above \nbed (m)',labelpad=0)
        
    
    ax_bot = fig.add_subplot(gspec[-1,:])        
    plot_gs_profile(profile,ax_bot)
        
    if fn_save != None:
        fig.savefig(fn_save)
        
    return fig, axs

def plot_gs_profile(profile, ax):

    for i in range(profile['gs']['gs'].shape[0]):
        ax.plot(profile['gs']['bin lower'],profile['gs']['gs'][i,:]
                ,label = profile['gs']['Isokinetic Sample Depth'][i])

    ax.set_xscale('log')
    ax.legend(title = 'sample depth (m)')
    ax.set_xlabel('grain size (um)')
    ax.set_ylabel('volume % \n (right?)')
        
    return ax

def profile_map_figure(profiles,A,fn_save = None):
    #profiles:      grain size profiles
    #A              dictionary of ADCP profiles as read in by rdi_readin_adcp_VariableBins

    fig, ax = plt.subplots()

    for key in profiles.keys():
        profiles[key]['buffer_poly'].plot(ax=ax,facecolor='none',edgecolor='purple')
        profiles[key]['station_pt'].plot(ax=ax, markersize=10, color='red')
        
        for ii in range(len(profiles[key]['adcp']['pos'])):
            profiles[key]['adcp']['pos'][ii].plot(ax=ax,color='black')
        
    cx.add_basemap(ax, source=cx.providers.CartoDB.Positron, zoom=12)
    
    if fn_save != None:
        fig.savefig(fn_save)
        
    return fig, ax