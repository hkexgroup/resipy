# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:49:22 2018
System check
@author: jamyd91
"""
import multiprocessing, platform, warnings, os # python standard libs
import re
from subprocess import PIPE, Popen

#%% ram amount check. 
#Now for complicated meshes we need alot more RAM. the below function is a os agnostic
#function for returning the amount of total ram. 
def systemCheck():
    print("________________System-Check__________________")
    #display processor info
    print("Processor info: %s"%platform.processor())
    num_threads = multiprocessing.cpu_count()
    print("Number of logical CPUs: %i"%num_threads)
    #this message will display if wine is not installed / detected
    helpful_msg ="""   
This version of pyR2 requires wine to run R2.exe, please consider installing
'wine is not an emulator' package @ https://www.winehq.org/. On linux wine can be found on
most reprositories (ubuntu/debian users can use "sudo apt install wine"). Wine acts as
a compatiblity layer between unix like OS systems (ie macOS and linux) and windows programs. 
    """
    msg_flag = False
    #check operating system 
    OpSys=platform.system()    
    if OpSys=='Darwin':
        print("Kernal type: macOS")
    else:
        print("Kernal type: %s"%OpSys)
    #check the amount of ram 
    if OpSys=="Linux":
        totalMemory = os.popen("free -m").readlines()[1].split()[1]
        #detect wine 
        is_wine = os.popen("wine --version").readlines()#[0].split()[0]
        if is_wine.find("wine") == -1:
            warnings.warn("Wine is not installed!", Warning)
            msg_flag = True
        else:
            wine_version = is_wine.split()[0].split('-')[1]
            print("Wine version = "+wine_version)
                          
    elif OpSys=="Windows":
        info = os.popen("systeminfo").readlines()
        for i,line in enumerate(info):
            if line.find("Total Physical Memory")!=-1:
                temp = line.split()[3]
                idx = temp.find(',')
                totalMemory = temp[0:idx]
                totalMemory += temp[idx+1:]
                
    elif OpSys=='Darwin':
        sysinfo = []
        info = Popen(['system_profiler SPHardwareDataType'], shell = True, stdout=PIPE, universal_newlines=True)
        for stdout_line in iter(info.stdout.readline, ''):
            sysinfo.append(stdout_line)
        memoryLine = [s for s in sysinfo if any(xs in s for xs in ['Memory'])] 
        totalMemory = re.findall('\\d+', memoryLine[0]) 
        totalMemory = int(totalMemory[0])*1000
        #detect wine 
        wineVersion = []
        is_wine = Popen(['/usr/local/bin/wine --version'], stdout=PIPE, shell=True, universal_newlines=True)
        for stdout_line in iter(is_wine.stdout.readline, ""):
            wineVersion.append(stdout_line)
        if not wineVersion:
            warnings.warn("Wine is not installed!", Warning)
            msg_flag = True
        else:
            wine_version = stdout_line.split()[0].split('-')[1]
            print("Wine version = "+wine_version)
#        totalMemory = os.popen("hwprefs memory_size").readlines()[0].split()[0]
#        totalMemory = int(totalMemory)*1000
#        #detect wine 
#        is_wine = os.popen("wine --version").readlines()#[0].split()[0]
#        if is_wine.find("wine") == -1:
#            warnings.warn("Wine is not installed!", Warning)
#            msg_flag = True
#        else:
#            wine_version = is_wine.split()[0].split('-')[1]
#            print("Wine version = "+wine_version)
        
    else:
        raise OSError("unrecognised/unsupported operating system")
        
    totalMemory = int(totalMemory)
    print("Total RAM available: %i Mb"%totalMemory)
    
    #print some warnings incase the user has a low end PC
    if totalMemory <= 4000:
        warnings.warn("The amount of RAM currently installed is low (<4Gb), complicated ERT problems may incur memory access voilations", Warning)
    if num_threads <=2:
        warnings.warn("Only one or two CPUs detected, multithreaded workflows will not perform well.", Warning)
    if msg_flag:
        print(helpful_msg)
    
    return {'memory':totalMemory,'core_count':num_threads,'OS':OpSys}

#%% 
info = systemCheck()