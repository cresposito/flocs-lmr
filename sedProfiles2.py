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
    

def readInventory(inventoryFile,datafolder):
    df=pd.read_excel(inventoryFile)

    #add a new column with a unique identifier for each profile 
    # (assumes one profile per station)
    profiles=df['Station'].unique()
    df['ProfileNum']=None

    for i, st in enumerate(profiles):
        df.loc[df['Station'] == st, 'ProfileNum'] = i+1

    #reorder the columns so that 'ProfileNum' appears immediately after 'Station'
    cols=df.columns.tolist()
    index_station = cols.index('Station')
    index_profile = cols.index('ProfileNum')
    removed_element = cols.pop(index_profile)
    cols.insert(index_station + 1, removed_element)

    # rearrange dataframe columns
    df=df[cols]

    #add a column with the full ADCP path
    df['filepath']=None
    indices=df.index.tolist()
    for ii in indices:
        df.at[ii,'filepath']=find_file(datafolder, df.at[ii,'Stationary ADCP File Name'])

    #add a column for total depth in m
    df['Total Depth (m)']=df['Total Depth (ft)']*0.3048

    return df

def MR173_Plots(df,savePath):
    #produces figures for MR173
    
    stations=df['Station'].unique().tolist()

    # del(stations[0:17])   #temporary, for restarting

    stations = [sta for sta in stations if not pd.isna(sta)]    #filter out the nans.  not sure why they're there



    ctr=0
    for sta in stations:
        ctr=ctr+1
        # sta='SWPD-2'

        # indices with this station
        indices = df[df['Station'] == sta].index.tolist()
        print(sta)
        print(indices)

        try:
            df = ADCP_MR173_Plot1(df, indices, savePath)
        except:
            print('ERROR IN STATION ' + sta)
        
        
        # if ctr>2: break   #use for development, to run just once


def ADCP_MR173_Plot1(df, indices,savePath):

    fig = plt.figure(figsize=(10,2+2*len(indices)))
    gs = gridspec.GridSpec(len(indices)+1, 4)
    gs.update(wspace=0.25, hspace=0.5) # set the spacing between axes.

    stationDepth=df.iloc[indices[0]]['Total Depth (m)']
    ylims=[0-0.1*stationDepth,1.1*stationDepth]

    formatted_depth = f"{df.iloc[indices[0]]['Total Depth (m)']:.1f}"
    FigTitle = (df.iloc[indices[0]]['Station'] + "\n" +
                "Q Belle Chasse: -999" + "\n" +
                str(df.iloc[indices[0]]["Date"].date()) + "\n" +
                "bed depth (m): " + formatted_depth)
    
    fig.suptitle(FigTitle)


    for ii in range(0,len(indices)):
        adcp_path=df.iloc[indices[ii]]['filepath']    
        ds=dolfyn.read(adcp_path)

        ens=[int(df.iloc[indices[ii]]['ADCP Ensemble Start']),int(df.iloc[indices[ii]]['ADCP Ensemble End'])]        
        ens_ind=[]
        ens_ind.append(np.where(ds.number.values == ens[0])[0].item())
        ens_ind.append(np.where(ds.number.values == ens[1])[0].item())

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
        z=stationDepth-ds['range']

        # clip for ensembles
        V_mag = V_mag[:, ens_ind[0]:ens_ind[1]]
        t = t[ens_ind[0]:ens_ind[1]]

        if ii==0:
            V_mag_total = V_mag
            t_total = t
        else:
            V_mag_total = np.concatenate((V_mag_total,V_mag),axis=1)
            t_total = t_total + t

        # velocity pcolor
        ax = fig.add_subplot(gs[ii,0:3])
        try:
            pcm = ax.pcolormesh(
                t, 
                z, 
                V_mag, 
                cmap='Blues', shading='nearest'
                )
        except:
            pass
        ax.set_ylim(ylims)
        pcm.set_clim([0, 1.5])
        ax.plot([t[0], t[-1]], [0, 0], linestyle=':', linewidth=2, color='black')
        ax.plot([t[0], t[-1]], [stationDepth, stationDepth], linestyle=':', linewidth=2, color='blue')
        ax.xaxis.set_major_formatter(dt.DateFormatter('%H:%M:%S'))
        ax.xaxis.set_major_locator(dt.SecondLocator(bysecond=range(0, 60, 15)))  # Tick every 15 seconds starting on the minute
        plt.xticks(rotation=25)

        # velocity profile
        ax = fig.add_subplot(gs[ii,3])
        ax.plot(np.nanmean(V_mag,axis=1),z)
        ax.plot([0, 9999], [0, 0], linestyle=':', linewidth=2, color='black')
        ax.plot([0, 9999], [stationDepth, stationDepth], linestyle=':', linewidth=2, color='black')
        ax.set_ylim(ylims)
        ax.set_xlim([0, 1.5])

        
        u=np.nanmean(V_mag,axis=1)
        lnZ=np.log(z)
        lnZ=lnZ.values
        
        jj=~np.isnan(u)
        # kk=~np.isnan(np.squeeze(lnZ,axis=1))
        kk=~np.isnan(lnZ)
        ee=jj&kk       
        u=u[ee]
        lnZ=lnZ[ee]

        lnZ=lnZ.reshape(-1,1)
        regr=linear_model.LinearRegression()
        
        try:
            regr.fit(lnZ,u)
        except ValueError as e:
            # Handle the ValueError (or any exception) here
            print(f"An error occurred: {e}")
            pass  # Optionally, handle the error gracefully or do nothing

        k=0.41
        ustar=regr.coef_*k
        z0=np.exp(-1*k*regr.intercept_/ustar)

        zz=np.arange(0.1,8,0.1)
        uu=(ustar/k)*np.log(zz/z0)
        ax.plot(uu,zz)
# 
        # eq_str='u(z) = '
        ustar_str = f'{ustar.item():.3f}'
        z0_str = f"{z0.item():.3f}"
        ax.text(0.1, 0.9, 'u*:' + ustar_str + '\nz_0:' + z0_str, transform=ax.transAxes, ha='left', va='top', fontsize=8)

        #add u* and z0 to the dataframe (not for final combined panel)
        if not 'u*' in df.columns:
            df['u*']=None
            df['z0']=None
        
        df.at[indices[ii],'u*']=ustar
        df.at[indices[ii],'z0']=z0
    
    # add the total figure
        # velocity pcolor
    ax = fig.add_subplot(gs[ii+1,0:3])
    pcm = ax.pcolormesh(
        range(0,len(t_total)), 
        z, 
        V_mag_total, 
        cmap='Blues', shading='nearest'
        )
    ax.set_ylim(ylims)
    pcm.set_clim([0, 1.5])
    ax.plot([0, len(t_total)], [0, 0], linestyle=':', linewidth=2, color='black')
    ax.plot([0, len(t_total)], [stationDepth, stationDepth], linestyle=':', linewidth=2, color='blue')
    # ax.xaxis.set_major_formatter(dt.DateFormatter('%H:%M:%S'))
    # ax.xaxis.set_major_locator(dt.SecondLocator(bysecond=range(0, 60, 15)))  # Tick every 15 seconds starting on the minute
    # plt.xticks(rotation=25)

    # velocity profile
    ax = fig.add_subplot(gs[ii+1,3])
    ax.plot(np.nanmean(V_mag_total,axis=1),z)
    ax.plot([0, 9999], [0, 0], linestyle=':', linewidth=2, color='black')
    ax.plot([0, 9999], [stationDepth, stationDepth], linestyle=':', linewidth=2, color='black')
    ax.set_ylim(ylims)
    ax.set_xlim([0, 1.5])

    
    u=np.nanmean(V_mag_total,axis=1)
    lnZ=np.log(z)
    lnZ=lnZ.values

    jj=~np.isnan(u)
    # kk=~np.isnan(np.squeeze(lnZ,axis=1))
    kk=~np.isnan(lnZ)
    ee=jj&kk       
    u=u[ee]
    lnZ=lnZ[ee]

    lnZ=lnZ.reshape(-1,1)
    print(lnZ.shape)
    print(u.shape)
    regr=linear_model.LinearRegression()
    regr.fit(lnZ,u)

    k=0.41
    ustar=regr.coef_*k
    z0=np.exp(-1*k*regr.intercept_/ustar)

    zz=np.arange(0.1,8,0.1)
    uu=(ustar/k)*np.log(zz/z0)
    ax.plot(uu,zz)
# 
    # eq_str='u(z) = '
    ustar_str = f'{ustar.item():.3f}'
    z0_str = f"{z0.item():.3f}"
    ax.text(0.1, 0.9, 'u*:' + ustar_str + '\nz_0:' + z0_str, transform=ax.transAxes, ha='left', va='top', fontsize=8)

    # save the figure
    filename = "ADCP_"+df.iloc[indices[0]]['Station']+"_"+str(df.iloc[indices[0]]['Date & Time'])+".png"
    illegal_chars = {'<', '>', ':', '"', '/', '\\', '|', '?', '*'}
    sanitized_filename = ''.join(c if c not in illegal_chars else '_' for c in filename)

    fig.savefig(os.path.join(savePath, sanitized_filename))

    #save the dataframe. overwriting each time
    df.to_excel(os.path.join(savePath,'Inventory.xlsx'))

    return df


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