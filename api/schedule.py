# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 13:52:54 2018
generate a scheduling file for a dipole dipole array 
@author: jamyd91
"""
import numpy as np
from tqdm import tqdm 

def dpdp_schedule(e_num):
    
    e_idx = np.arange(e_num)+1 # individual electrode indexes
    
    pair_skip = np.arange(e_num-4) # number of electrodes to skip between each pair
    # | | --pair_skip-- | |
    
    e_skip = np.arange(int((e_num-4)/2)) # number of electrodes in between the charge or potential pair
    # |--e_skip--| --pair_skip-- |--e_skip--|
    
    layer_count = 0
    meas_count = 0
    cache = ''#'c1 c2 p1 p2\n'
    for es in tqdm(e_skip,desc='scheduler',ncols=100):
        for ps in pair_skip:
            layer_count += 1
            for i in e_idx:
                if i+ps+es+3>max(e_idx):
                    break
                meas_count += 1
                c1 = i
                c2 = i+es+1
                p1 = i+es+ps+2
                p2 = i+es+ps+3
                cache = cache + '{}\t{}\t{}\t{}\t{}\n'.format(meas_count,p1,p2,c1,c2)
    return str(meas_count)+'\n'+cache
                
        
            
#%% testing code (legacy)
    
# =============================================================================
# #%%convert res2dmod output into a protocal.dat
# import tkinter as tk
# from tkinter import filedialog
# import parsers as prs 
# 
# #import dat file
# print("please select the res2dmod mesh file you want to convert.\n")
# root=tk.Tk()
# root.withdraw()
# file_path=filedialog.askopenfilename(title='Select mesh file',filetypes=(("data files","*.dat"),("all files","*.*")))            
#        
# #%%
# elec,df = prs.res2invInputParser(file_path)
# prs.dataframe2dat(df)
# 
# 
# #%% test   genGeoFile_adv 
# import gmshWrap as gw
# 
# surf_x = [-1,10]
# surf_y = [0,0]
# elec_x = np.linspace(2,7,5)
# elec_y = [0]*5
# string1x = [1]*10
# string1y = np.linspace(-1,-10,10)
# string2x = [8]*10
# string2y = np.linspace(-1,-10,10)
# poly1x = [3,3,6,6]
# poly1y = [-2,-5,-5,-2]
# bound1x = [3.5,4.5,5.5]
# bound1y = [-6,-7,-8]
# 
# geom_input = {'surface': [surf_x,surf_y],
#               'electrode':[elec_x,elec_y],
#               'borehole1':[string1x,string1y],
#               'borehole2':[string2x,string2y],
#               'boundary1':[bound1x,bound1y],
#               'polygon1':[poly1x,poly1y]}
# 
# =============================================================================
#nodes,file_name = gw.genGeoFile_adv(geom_input, doi = 11)