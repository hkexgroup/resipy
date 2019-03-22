# -*- coding: utf-8 -*-
#Template python scripts run outside of pyR2

parallelScript = r"""#import relevant modules
from multiprocessing import Pool
import subprocess, shutil, os, sys

#directories in which to run Binley code 
dirs = {:s}

def execute(cmd):
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW  
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
    for text in execute(cmd): # EXCUTE THE BINLEY CODE
        print(text.rstrip())
    #now put files back into main directory as if it was a sequential inversion
    try:
        shutil.copy('f001_res.vtk',r'..\f%003i_res.vtk'%num)
        shutil.copy('f001_res.dat',r'..\f%003i_res.dat'%num)
        shutil.copy('f001_err.dat',r'..\f%003i_err.dat'%num)
    except FileNotFoundError:
        shutil.copy('f001.vtk',r'..\f%003i.vtk'%num)
        shutil.copy('f001.dat',r'..\f%003i.dat'%num)
        shutil.copy('f001.err',r'..\f%003i.err'%num)
    shutil.copy('electrodes.vtk',r'..\electrodes%003i.vtk'%num)
    shutil.copy('electrodes.dat',r'..\electrodes%003i.dat'%num)

print('----------- START OF PARALLISED INVERSION -------------')
print('--------------- PROCESSING %i DATASETS -----------------'%len(dirs))
if __name__ == '__main__':
    pool = Pool({:d})
    for i,_ in enumerate(pool.imap(worker, dirs)):
        print('\rCompleted job %i'%(i+1))
    pool.close()  
    pool.join()
"""

startAnmt = """from paraview.simple import * 
def start_cue(self):
	global annotations
	global maxIndex
	text_obj = Text()#make a text object
	annotations= []\n"""
    
endAnmt ="""	maxIndex = len(annotations)
def tick(self):
	global annotations
	global maxIndex
	index = int( self.GetClockTime() )
	if index >= maxIndex :
		 index = maxIndex - 1
	textSource = paraview.simple.FindSource('Text1')
	textSource.Text = annotations[index]
def end_cue(self): pass"""