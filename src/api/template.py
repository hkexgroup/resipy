# -*- coding: utf-8 -*-
#Template python scripts run outside of pyR2

parallelScript = r"""#import relevant modules
from multiprocessing import Pool
import subprocess
import shutil
import os

startupinfo = subprocess.STARTUPINFO()
startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW   

#directories in which to run Binley code 
dirs = {:s}

def execute(cmd):
    proc = subprocess.Popen(cmd, stdout = subprocess.PIPE, shell=False, universal_newlines=True, startupinfo=startupinfo)           
    for stdout_line in iter(proc.stdout.readline, ""):
        yield stdout_line
    proc.stdout.close()
    return_code = proc.wait()
    if return_code:
        print('error on return_code')        
      
def worker(dirname):
    exe_loc = {:s}
    os.chdir(dirname)
    num = int(dirname.split('\\')[-1]) # index of the survey
    cmd = '"'+exe_loc+'"'
    #print("\nRunning inversion number : %i"%num) # show that something is happening? 
    for text in execute(cmd): # EXCUTE THE BINLEY CODE
        print(text.rstrip())
    #now put files back into main directory as if it was a sequential inversion
    try:
        shutil.copy('f001_res.vtk',r'..\f%003i_res.vtk'%num)
        shutil.copy('f001_res.dat',r'..\f%003i_res.dat'%num)
        shutil.copy('f001_res.err',r'..\f%003i_res.err'%num)
    except FileNotFoundError:
        shutil.copy('f001.vtk',r'..\f%003i.vtk'%num)
        shutil.copy('f001.dat',r'..\f%003i.dat'%num)
        shutil.copy('f001.err',r'..\f%003i.err'%num)
    shutil.copy('electrodes.vtk',r'..\electrodes%003i.vtk'%num)
    shutil.copy('electrodes.dat',r'..\electrodes%003i.dat'%num)

if __name__ == '__main__':
    pool = Pool({:d})
    pool.map(worker, dirs)
    pool.close()  
    pool.join()
"""