# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 10:29:45 2018
post processing errors in r2gui 
@author: jamyd91
"""
#import python standard libraries
#import tkinter as tk
#from tkinter import filedialog
#import conda libraries
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd

#return file path string 
#def find_file_path():
#    root=tk.Tk()
#    root.withdraw()
#    file_path=filedialog.askopenfilename(title='Select file',filetypes=(("data files","*.dat"),("all files","*.*")))
#    return file_path

#%% open file 
def import_R2_err_dat(file_path):
    """ Imports _err.dat file produced by R2.exe after inversion.
    
    Parameters
    ---------
    file_path: string
        Maps to the error.dat file.
    
    Returns
    ---------
    error_dict: dictionary
        A dictionary containing lists of each of the different parameters.
    """
    
    fh = open(file_path,'r')
    dump = fh.readlines()
    fh.close()
    line1 = dump[0].strip().split("  ")
    headers = [head for i,head in enumerate(line1) if head != ""]
    
    #% header dictionary
    if "Observed_Phase" in headers:
        error_dict = pd.read_csv(file_path, delim_whitespace=True)
        
    else:
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

def disp_R2_errors(error_dict, ax=None):
    """
    Displays the normalised errors associated with the inversion. 
    
    Parameters
    ---------
    error_dict: dictionary
        Varaible returned from import_R2_err_dat.
    ax: matplotlib axis handle, optional
        Resulting figure plots on to this axis.
    """
    # plot normalised errors
    measurement_no = np.arange(1,len(error_dict["Normalised_Error"])+1)
    #make figure
    if ax is None: 
        fig, ax = plt.subplots()    
    
    ax.scatter(measurement_no,error_dict["Normalised_Error"])
    ax.set_ylabel("Normalised Error")
    ax.set_xlabel("Measurement Number")
    #add diagnositic lines
    y_pos_limit = (3,3)
    y_neg_limit = (-3,-3)
    baseline = (0,0)
    ax.plot((1,measurement_no[-1]),y_pos_limit,'r--')
    ax.plot((1,measurement_no[-1]),y_neg_limit,'r--')
    ax.plot((1,measurement_no[-1]),baseline,'k--')
    
    
    
#%% test block 
#file_path = 'test/f001_err.dat'
#error_info2 = import_R2_err_dat(file_path)
#disp_R2_errors(error_info)