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
    proc = subprocess.run(cmd, startupinfo=startupinfo)           
    if proc.returncode < 0 : # is negative if there was a problem (using run not Popen)
        print('error on return_code')     
      
def worker(dirname):
    exe_loc = {:s}
    os.chdir(dirname)
    num = int(dirname.split('\\')[-1]) # index of the survey
    cmd = '"'+exe_loc+'"'
    execute(cmd) # EXCUTE THE BINLEY CODE
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
    for f in os.listdir(dirname):
        if f.endswith('.out'):
            shutil.copy(f,r'..\inversion%003i.log'%num)
            
#def progress(iteration,total): # display a progress bar? not currently in use 
#    iteration += 1
#    barLength = 38 # Modify this to change the length of the progress bar
#    progress = iteration/total
#    if isinstance(progress, int):
#        progress = float(progress)
#    block = int(round(barLength*progress))
#    counter = ' %i/%i'%(iteration,total)
#    text = "\rCompleted: [%s]"%("#"*block + "-"*(barLength-block))+counter
#    sys.stdout.write(text)
#    sys.stdout.flush()

to_process=len(dirs)
print('\n----------- START OF PARALLISED INVERSION -------------',flush=True)
print('--------------- PROCESSING %i DATASETS -----------------\n'%to_process,flush=True)
if __name__ == '__main__':
    #progress(-1,to_process)
    pool = Pool({:d})
    for i,_ in enumerate(pool.imap(worker, dirs)):
        print('Completed job %i\n'%(i+1),flush=True)
        #progress(i,to_process)
    pool.close()  
    pool.join()
else:
    print("Parallel Script has not been run as main module, therefore parallel process has not been initiated")
print('\n',flush=True)
    
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