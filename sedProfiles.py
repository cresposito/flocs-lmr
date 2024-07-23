import dolfyn
import numpy as np
from dolfyn.adp import api
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as dt
import os
import pandas as pd


def readInventory(inventoryFile):
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

    return df

def read_ADCP_files_MR173(df,dataFolder):
    #reads in all ADCP data for the MR_173 project
    
    #make new column to contain adcp data
    df['ADCP V_mag']=None
    df['ADCP t']=None
    df['ADCP bins']=None

    for ind,row in df.iterrows():
        #locate ADCP file in dataFolder or its subfolders, and read selected ensembles
        # print(row)
        if isinstance(row['Stationary ADCP File Name'],str):
            adcp_path=find_file(dataFolder,row['Stationary ADCP File Name'])    
            ens=[int(row['ADCP Ensemble Start']),int(row['ADCP Ensemble End'])]        
            ens=ens

            V_mag, t, bin_centers = read_ADCP(adcp_path, ens)
            t = np.array(t)
            V_mag = np.array(V_mag)
            ds = np.array(bin_centers)
            
            df.loc[ind, 'ADCP V_mag'] = [V_mag]
            df.loc[ind, 'ADCP t'] = [[t]]
            df.loc[ind, 'ADCP bins'] = [[bin_centers]]

    return df


def read_ADCP(adcp_path, ens):
    ds=dolfyn.read(adcp_path)
    t = dolfyn.time.dt642date(ds.time)

    # Calculate the mean distance to the bottom for each time
    mean_dist = ds["dist_bt"].mean(dim=["beam"])

    # Create a boolean mask where the condition is true
    mask = mean_dist > ds["range"]

    V_masked=ds["vel"].where(mask, np.nan)
    V_mag=np.linalg.norm(V_masked,axis=0)

    # clip for ensembles and return
    V_mag = V_mag[:, ens[0]:ens[1]+1]
    t = t[ens[0]:ens[1]+1]

    return V_mag, t, ds["range"]

    

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