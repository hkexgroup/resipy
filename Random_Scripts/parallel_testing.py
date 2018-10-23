#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 12:20:01 2018
parellisation testing for batch inversion 
@author: jimmy
"""
import os,sys, platform, time #utility standard python packages 
pyR2_location = r".." #insert r2gui location here
sys.path.append(pyR2_location)
from api.R2 import R2 #import r2 class 
from api.meshTools import checkRAM
#multiprocessing lib
import multiprocessing as mp #python standard library for multiprocessing 
import shutil # python library - manages copying files etc
import random#these 2 modules will be used to create random folder names in the invdir 
import string
from subprocess import PIPE, Popen, call#python library to manage running os level commands 

#%%Return some information about the system 
print("Processor type: %s"%platform.processor())
num_threads = mp.cpu_count()
print("Number of logical processors: %i"%num_threads)
if num_threads<3:
    print("It looks like your computer only has 1 or 2 logical cores, consider getting a better PC you scrub")
totalMemory = checkRAM()
#print("Total RAM available: %iMb"%int(totalMemory))
if int(totalMemory)<4000:
    print("Your computer has a small amount of RAM, parallel processing is ill advised")

#%%create batch survey
k = R2()
test_dir = os.path.join(pyR2_location,'api/test/testTimelapse')
k.createBatchSurvey(test_dir)
k.dirname
k.invert()
#okay that works pretty well. 

#%% test making multiple directories 
#for now we will just use the files created by the previous block of code for a test case
letters = string.ascii_lowercase + string.ascii_uppercase
def parallel_R2(working_directory,time_step=0,output=None):
    """
    working_directory : str
    time_step : int 
    output : object (currently unimplimented)
        Pass a mp.Queue object here, useful for retrieving results 
    """
    if not os.path.isdir(working_directory):
        raise EnvironmentError("Supplied working directory apparently doesn't exist")
    #create random folder
    itr_dir_name = ''.join(random.choice(letters) for i in range(5))
    path = os.path.join(working_directory,itr_dir_name)
    while os.path.isdir(path):#this while loop adds in a fail safe incase the folder has already been created
        itr_dir_name = ''.join(random.choice(letters) for i in range(5))#make a new random folder
        path = os.path.join(working_directory,itr_dir_name)
    os.makedirs(path)
    print(path)
    #copy R2 and relevant input files to that random folder
    shutil.copyfile(os.path.join(working_directory,'R2.exe'),os.path.join(path,'R2.exe'))
    shutil.copyfile(os.path.join(working_directory,'protocol.dat'),os.path.join(path,'protocol.dat'))
    shutil.copyfile(os.path.join(working_directory,'R2.in'),os.path.join(path,'R2.in'))
    shutil.copyfile(os.path.join(working_directory,'res0.dat'),os.path.join(path,'res0.dat'))
    
    #run R2 
    if platform.system() == "Windows":#command line input will vary slighty by system 
        cmd_line = 'R2.exe'
    else:
        cmd_line = ['wine', 'R2.exe']
    
    os.chdir(path) #change to random folder   
    call(cmd_line)#run R2
    os.chdir(working_directory)
    #for arguments sake here I'm copying the output vtk file back in to the working directory once the 
    #inversion has finished. 
    result_name='time_step'+str(time_step)+'.vtk'
    shutil.copyfile(os.path.join(path,'f001_res.vtk'),os.path.join(working_directory,result_name))
    #Then delete the created folder otherwise things will get very messy!
    shutil.rmtree(path)
    return True # the function has completed 
    
#%%now actually put the above into practice 
t0 = time.time()
for i in range(6,10):
    parallel_R2(k.dirname,11) # first check the function works and do a computation in series 
print("time taken to do inversions in series %f s"%(time.time()-t0))

#multiprocessing using the pool method 
t1 = time.time()
pool = mp.Pool(processes=int(num_threads/2)) # now do computation in parallel? 
results = [pool.apply_async(parallel_R2, args=(k.dirname,x,)) for x in range(2)]
results = [p.get() for p in results]

print("time taken to do inversions in parallel %f s"%(time.time()-t1))

## Define an output queue
#t1 = time.time()
#output = mp.Queue()
#processes = [mp.Process(target=parallel_R2, args=(k.dirname, i, output)) for i in range(4)]
## Run processes
#for p in processes:
#    p.start()
#
## Exit the completed processes
#for p in processes:
#    p.join()
#print("time taken to do inversions in parallel %f s"%(time.time()-t1))    
    
    
    
    