# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:49:22 2018
System check
@author: jamyd91
"""
import multiprocessing, platform, warnings, os # python standard libs

#%% ram amount check. 
#Now for complicated meshes we need alot more RAM. the below function is a os agnostic
#function for returning the amount of total ram. 
def systemCheck():
    print("________________System-Check__________________")
    #display processor info
    print("Processor info: %s"%platform.processor())
    num_threads = multiprocessing.cpu_count()
    print("Number of logical CPUs: %i"%num_threads)
    #check operating system 
    OpSys=platform.system()
    if OpSys=='Darwin':
        print("Kernal type: macOS")
    else:
        print("Kernal type: %s"%OpSys)
    #check the amount of ram 
    if OpSys=="Linux":
        totalMemory = os.popen("free -m").readlines()[1].split()[1]
        is_wine = os.popen("wine --version").readlines()[0].split()[0]
        if is_wine.find("wine") == -1:
            warnings.warn("It looks like wine is not installed")
        elif is_wine != "wine-3.0":
            warnings.warn("Using wine version 3.0 and later is reccomended")
            
    elif OpSys=="Windows":
        info = os.popen("systeminfo").readlines()
        for i,line in enumerate(info):
            if line.find("Total Physical Memory")!=-1:
                temp = line.split()[3]
                idx = temp.find(',')
                totalMemory = temp[0:idx]
                totalMemory += temp[idx+1:]
                
    elif OpSys=='Darwin':
        totalMemory = os.popen("hwprefs memory_size").readlines()[0].split()[0]
        totalMemory = int(totalMemory)*1000
        
    else:
        raise OSError("unrecognised/unsupported operating system")
        
    totalMemory = int(totalMemory)
    print("Total RAM available: %i Mb"%totalMemory)
    
    #print some warnings incase the user has a low end PC
    if totalMemory <= 4000:
        warnings.warn("The amount of RAM currently installed is quite low (<4Gb), complicated problems may incur memory access voilations")
    if num_threads <=2:
        warnings.warn("Your core count is quite small, multithreaded workflows will not perform well.")
    
    return {'memory':totalMemory,'core_count':num_threads,'OS':OpSys}

#%% 
systemCheck()