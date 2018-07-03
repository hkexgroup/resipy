# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 17:13:55 2017

@author: Sina
"""

import numpy as np
import os
import pandas as pd
import easygui as eg
#%%
datadirectory = 'F:\Sina\Google Drive\Python_code\With DCA\DCA data\Protocols\L3\ew' 
working = os.environ.get("WORKING_DIRECTORY", datadirectory) #add the folder with the data
os.chdir(working)
#%%
#b="\ ".strip()
#if datadirectory[-1:]==b:
#    fdir = datadirectory
#else:
#    fdir = datadirectory+b
#%% Functions
def commons_normal (data1,data2,num_data): #finding common measurements among time-lapse normal datasets (for time-lapse inversion)
    shape_array = pd.DataFrame(np.zeros([len(data1),2])).rename(columns={0:'length',1:'order'})
    for i in range(len(data1)):
        shape_array.iloc[i,:]= [len(data1[i]),i]
    common_index = shape_array.sort_values(by='length').reset_index(drop = True).iloc[0,1]
    common_final = []
    for i in range(num_data):
        temp_common = pd.merge(data2[i], data1[int(common_index)][['M', 'N', 'A', 'B']], how = 'inner').drop_duplicates(subset = ['M', 'N', 'A', 'B'])
        common_final.append(temp_common)
    min_num_array = pd.DataFrame(np.zeros([len(common_final),2])).rename(columns={0:'length',1:'order'})
    for i in range(len(common_final)):
        min_num_array.iloc[i,:]= [len(common_final[i]),i]
    min_nums = min_num_array.sort_values(by='length').reset_index(drop = True).iloc[0,0]
    return (common_final,min_nums)
#%% inputs
msg=['Select the file(s) of interest']
filestobeused = eg.fileopenbox(msg,multiple=True)
#%%
files=[]
names=[]
for infile in filestobeused:   # reading data
    name = os.path.basename(infile)[:-4]
#    data =  pd.read_csv(infile, sep = '\t', skiprows = 0, header =None)
    data =  pd.read_csv(infile, sep = '\t',skiprows = 1,header=None).rename(columns = {0:'Num',1:'M',2:'N',3:'A',4:'B',5:'R',6:'IP',7:'Rerr',8:'IPerr'})
    files.append(data)
    names.append(name)
#%%
common_normal,nums = commons_normal(files,files,len(files))
print(nums)
n,f=1,0
common_good =[]
for i in range(len(common_normal)):
    if len(common_normal[i])>=300:
        temp_common_good = common_normal[i]
        common_good.append(temp_common_good)
while n!=f:
    common_good,n = commons_normal(common_good,common_good,len(common_good))
    temp = pd.concat(common_good)
    f=len(temp)/len(common_good)
    if f>n-0.2 and f<n+0.2:
        break
#%%
for i in range(len(filestobeused)):
    common_good[i]['Num'] = np.arange(len(common_good[i]))+1
    save_file = open('protocol_%s.dat' % (names[i]), 'w')
    save_file.write(str(len(common_good[i])) + '\n')
    save_file.close()
    with open('protocol_%s.dat' % (names[i]), mode='a') as file:
        common_good[i].to_csv(file,sep='\t',header = None, index=None)