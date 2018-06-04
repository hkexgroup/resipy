# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 10:29:45 2018
post processing errors in r2gui 
@author: jamyd91
"""
#import python standard libraries
import tkinter as tk
from tkinter import filedialog
#import conda libraries
import matplotlib.pyplot as plt
import numpy as np 

#return file path string 
def find_file_path():
    root=tk.Tk()
    root.withdraw()
    file_path=filedialog.askopenfilename(title='Select file',filetypes=(("data files","*.dat"),("all files","*.*")))
    return file_path

#%% open file 
def import_R2_err_dat(file_path):
    fh = open(file_path,'r')
    dump = fh.readlines()
    fh.close()
    line1 = dump[0].strip().split("  ")
    headers = [head for i,head in enumerate(line1) if head != ""]
    
    #% header dictionary
    error_dict = {"Normalised_Error":[],
                  "Observed_Rho":[],  
                  "Calculated_Rho":[],
                  "Original_Weight":[],
                  "Final_Weight":[],
                  "C":[],
                  "P+":[],
                  "P-":[],
                  "C+":[],
                  "C-":[]}
    
    for i in range(1,len(dump)):
        line_info = dump[i].split()
        error_dict["Normalised_Error"].append(float(line_info[0]))
        error_dict["Observed_Rho"].append(float(line_info[1]))
        error_dict["Calculated_Rho"].append(float(line_info[2]))
        error_dict["Original_Weight"].append(float(line_info[3]))
        error_dict["Final_Weight"].append(float(line_info[4]))
        error_dict["C"].append(int(line_info[5]))
        error_dict["P+"].append(int(line_info[6]))
        error_dict["P-"].append(int(line_info[7]))
        error_dict["C+"].append(int(line_info[8]))
        error_dict["C-"].append(int(line_info[9]))
        
    return error_dict

def disp_R2_errors(error_dict):    
    # plot normalised errors
    measurement_no = np.arange(1,len(error_dict["Normalised_Error"])+1)
    #make figure
    plt.figure()    
    plt.scatter(measurement_no,error_dict["Normalised_Error"])
    plt.ylabel("normalised error")
    plt.xlabel("measurement number")
    #add diagnositic lines
    y_pos_limit = (3,3)
    y_neg_limit = (-3,-3)
    baseline = (0,0)
    plt.plot((1,measurement_no[-1]),y_pos_limit,'r--')
    plt.plot((1,measurement_no[-1]),y_neg_limit,'r--')
    plt.plot((1,measurement_no[-1]),baseline,'k--')
    
    
    
#%% test block 
#open errors from r2 output file
file_path = find_file_path()

error_info = import_R2_err_dat(file_path)

disp_R2_errors(error_info)
