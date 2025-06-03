# -*- coding: utf-8 -*-
"""
This file is part of the ResIPy project (https://gitlab.com/hkex/resipy).
@licence: GPLv3
@author: ResIPy authors and contributors

The 'Project' class wraps all main interactions between R* executables
and other filtering or meshing part of the code. It's the entry point for
the user.
"""
ResIPy_version = '3.6.3' # ResIPy version (semantic versionning in use)

#import relevant modules
import os, sys, shutil, platform, warnings, time, glob # python standard libs
from subprocess import PIPE, call, Popen
import psutil
from copy import deepcopy
from threading import Thread

# used to download the binaries
import requests
import hashlib

import subprocess
import numpy as np # import default 3rd party libaries (can be downloaded from conda repositry, incl with winpython)
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.patches as mpatches
import matplotlib.path as mpath
from scipy.spatial import cKDTree

OS = platform.system()
sys.path.append(os.path.relpath('..'))

#import ResIPy resipy packages
from resipy.Survey import Survey, polyfit 
from resipy.parsers import geomParser
from resipy.r2in import write2in
import resipy.meshTools as mt
#from resipy.meshTools import cropSurface
from resipy.protocol import (dpdp1, dpdp2, wenner_alpha, wenner_beta, wenner,
                          wenner_gamma, schlum1, schlum2, multigrad)
from resipy.seqGen import Generator
from resipy.SelectPoints import SelectPoints
from resipy.saveData import (write2Res2DInv, write2csv, writeSrv)
from resipy.interpolation import rotGridData, invRotGridData, nearest3d

apiPath = os.path.abspath(os.path.join(os.path.abspath(__file__), '../'))
print('API path = ', apiPath)
print('ResIPy version = ',str(ResIPy_version))

warnings.simplefilter('default', category=DeprecationWarning) # this will show the deprecation warnings

'''NOTE
pre-processing and error models for unique, combined or multiple surveys:
    idea, using the value of the index argument to identify the scope of
    the function.
    index = -2 : apply/show data from combined survey (bigSurvey)
    index = -1 : apply to each datasets the same type of model
    index > 0 : apply an error model to the selected unique survey
'''
if int(np.version.version.split('.')[0]) >= 2: 
    # if numpy version > 2 then set print options to legacy, otherwise the gmsh wrapper will crash gmsh 
    np.set_printoptions(legacy='1.25') 

#%% check executables are here
def checkSHA1(fname):
    BUF_SIZE = 65536  # lets read stuff in 64kb chunks!
    sha1 = hashlib.sha1()
    with open(fname, 'rb') as f:
        while True:
            data = f.read(BUF_SIZE)
            if not data:
                break
            sha1.update(data)
    return sha1.hexdigest()

def checkExe(dirname):
    exes = ['cR2.exe','R3t.exe','cR3t.exe']#,'R2.exe','gmsh.exe']
    hashes = ['e35f0271439761726473fa2e696d63613226b2a5',
              'a3bf11e392f0628fe7241bfae401c32a10f116ed',
              '9337435f018264771470d5d4312908b0d1242af1',
              # '4aad36d5333ddf163c46bab9d3c2a799aa48716e',
              # '91bd6e5fcb01a11d241456479c203624d0e681ed'
              ]
    for i, exe in enumerate(exes):
        fname = os.path.join(dirname, exe)
        download = False
        if os.path.exists(fname) is not True:
            download = True
            print('{:s} not found, will download it...'.format(exe), end='', flush=True)
        else: # check if the file is up to date
            sha1 = checkSHA1(fname)
            if sha1 != hashes[i]:
                download = True
                print('{:s} needs to be updated...'.format(exe), end='', flush=True)
        if download:
            response = requests.get("https://gitlab.com/hkex/resipy/-/raw/master/src/resipy/exe/" + exe)
            with open(fname, 'wb') as f:
                f.write(response.content)
            print('done')
        else:
            print('{:s} found and up to date.'.format(exe))

# the below failed if no internet connection so let's put it in try/except                
try:
    checkExe(os.path.join(apiPath, 'exe'))
except Exception as e:
    pass

# little class for managing multiple processes (for parallel inversion)
class ProcsManagement(object): # little class to handle the kill
    def __init__(self, r2object):
        self.r2 = r2object
        self.killFlag = False
    def kill(self):
        self.killFlag = True
        print('killing...')
        self.r2.irunParallel2 = False # this will end the infinite loop
        procs = self.r2.procs # and kill the running processes
        for p in procs:
            p.terminate()
        print('all done')
        
#%% system check
def getSysStat():
    """Return processor speed and usage, and free RAM and usage. 
    
    Returns
    -------
    cpu_speed: float
        in Mhz.
    cpu_usage: float
        in percent. 
    ram_avail: float
        avialable memory.
    ram_usage: float
        in percent. 
    ram_total: float
        total memory in gb 
    """
    #processor info
    try: # for Apple silicon
        cpu_speed = psutil.cpu_freq()[0]
    except:
        cpu_speed = 0
    cpu_usage = psutil.cpu_percent()
        
    #check the amount of ram avialable 
    ram = psutil.virtual_memory()
    ram_avail = ram.available*9.31e-10
    ram_usage = ram.percent 
    ram_total = ram.total*9.31e-10
    return cpu_speed, cpu_usage, ram_avail, ram_usage, ram_total 

def getMacOSVersion(): # this is obsolete now
    OpSys=platform.system()    
    if OpSys=='Darwin':
        versionList = platform.mac_ver()[0].split('.')
        macVersion = float(versionList[0] + '.' + versionList[1]) # not getting patch version so xx.xx only
        if macVersion >= 10.15:
            return True
        else:
            return False

def whichWineMac():
    """
    Checks if 'wine' or 'wine64' should be used on macOS.

    Returns:
        str: The appropriate wine executable ('wine' or 'wine64'), or None if neither is found.
    """

    wine_paths = ['/usr/local/bin/wine', '/opt/homebrew/bin/wine', '/usr/bin/wine'] # common wine paths
    wine64_paths = ['/usr/local/bin/wine64', '/opt/homebrew/bin/wine64', '/usr/bin/wine64'] # common wine64 paths
    
    global wPath
    for path in wine_paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            try:
                # Basic check if wine is working
                subprocess.run([path, '--version'], capture_output=True, timeout=1)
                wPath = path
                return 'wine'
            except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
                pass # wine either timed out, was not found, or had an OS error.

    for path in wine64_paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            try:
                # Basic check if wine64 is working
                subprocess.run([path, '--version'], capture_output=True, timeout=1)
                wPath = path
                return 'wine64'
            except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
                pass # wine64 either timed out, was not found, or had an OS error.

    return None
        
def systemCheck(dump=print):
    """Performs a simple diagnostic of the system, no input commands needed. System
    info is printed to screen, number of CPUs, memory and OS. This check is 
    useful for parallel processing. 
    
    Parameters
    ----------
    dump : function
        stdout dump 
    
    Returns
    -------
    system_info: dict
        Dictionary keys refer information about the system 
    """
    dump("________________System-Check__________________")
    #check operating system 
    OpSys=platform.system()    
    if OpSys=='Darwin':
        dump("Kernel type: macOS")
    else:
        dump("Kernel type: %s"%OpSys)
    
    totalMemory = 0 # incase system can't figure it out!
    num_threads = 0
    
    #display processor info
    dump("Processor info: %s"%platform.processor())
    cpu_cores = psutil.cpu_count(logical=True)
    physical_cores = psutil.cpu_count(logical=False)
    if isinstance(physical_cores, int) is False:
        physical_cores = 1  # can be the case for raspberrypi

    try: # for Apple silicon
        max_freq = max(psutil.cpu_freq())
    except:
        max_freq = 0
    dump("%i Threads at <= %5.1f Mhz"%(cpu_cores,max_freq))
    dump("(%i)"%physical_cores)
    
    #check the amount of ram 
    ram = psutil.virtual_memory()
    totalMemory = ram.total*9.31e-10
    availMemory = ram.available*9.31e-10
    usage = ram.percent
    dump('Total memory = %3.1f Gb (usage = %3.1f)'%(totalMemory,usage))
    
    #wine check - this message will display if wine is not installed / detected
    helpful_msg ="""   
This version of ResIPy requires wine to run R2.exe, please consider installing
'wine is not an emulator' package @ https://www.winehq.org/. On linux wine can be found on
most reprositories (ubuntu/debian users can use "sudo apt install wine-stable"). Wine acts as
a compatiblity layer between unix like OS systems (ie macOS and linux) and windows programs. 
    """
    wineCheck = True
    msg_flag = False
    if OpSys=="Linux":
        #detect wine 
        p = Popen("wine --version", stdout=PIPE, shell=True)
        is_wine = str(p.stdout.readline())#[0].split()[0]
        if is_wine.find("wine") == -1:
            warnings.warn("Wine is not installed!", Warning)
            msg_flag = True
            wineCheck = False
        else:
            wine_version = is_wine.split()[0].split('-')[1]
            dump("Wine version = "+wine_version)
                          
    elif OpSys=="Windows":
        dump('Wine Version = Native Windows (N/A)')
                
    elif OpSys=='Darwin':
        if whichWineMac() is None:
            warnings.warn("Wine is not installed!", Warning)
            msg_flag = True
            wineCheck = False
        else:
            dump("Wine is installed")
        
        # try: # old method
        #     winetxt = 'wine'
        #     if getMacOSVersion(): # looks like new wine-stable is no longer wine64 and is just wine
        #         winetxt = 'wine64' 
        #     winePath = []
        #     wine_path = Popen(['which', winetxt], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
        #     for stdout_line in iter(wine_path.stdout.readline, ''):
        #         winePath.append(stdout_line)
        #     if winePath != []:
        #         is_wine = Popen(['%s' % (winePath[0].strip('\n')), '--version'], stdout=PIPE, shell = False, universal_newlines=True)
        #     else:
        #         global wPath
        #         try:
        #             is_wine = Popen(['/usr/local/bin/%s' % winetxt,'--version'], stdout=PIPE, shell = False, universal_newlines=True)
        #             wPath = '/usr/local/bin/'
        #         except:
        #             is_wine = Popen(['/opt/homebrew/bin/%s' % winetxt,'--version'], stdout=PIPE, shell = False, universal_newlines=True) # quick fix for M1 Macs
        #             wPath = '/opt/homebrew/bin/'
        #     wineVersion = []
        #     for stdout_line in iter(is_wine.stdout.readline, ""):
        #         wineVersion.append(stdout_line)
        #     wine_version = stdout_line.split()[0].split('-')[1]
        #     dump("Wine version = "+wine_version)
        # except:
        #     warnings.warn("Wine is not installed!", Warning)
        #     msg_flag = True
        #     wineCheck = False
        
    else:
        print(OpSys)
        raise OSError("unrecognised/unsupported operating system")
            
    if msg_flag:
        if dump != print:
            print(helpful_msg)
        else:
            dump(helpful_msg)
            
    #check graphics unit 
    dump("GPU info: "+str(mt.gpuinfo))
    
    return {'totalMemory':totalMemory,
            'availMemory':availMemory,
            'cpuCount':cpu_cores,
            'physicalCpuCount':physical_cores, 
            'maxFreq':max_freq,
            'OS':OpSys,
            'wineCheck':wineCheck,
            'GPU':mt.gpuinfo}

def donothing(x):
    pass
sysinfo = systemCheck(dump=donothing)

#%% useful function for directory management 
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)
    
    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)
    
    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

#%% other helper functions 
# distance matrix function for 2D (numpy based from https://stackoverflow.com/questions/22720864/efficiently-calculating-a-euclidean-distance-matrix-using-numpy)
def cdist(a):
    z = np.array([complex(x[0], x[1]) for x in a])
    return np.abs(z[...,np.newaxis]-z)

def addCustSeq(fname): # add custom sequence 
    seq = pd.read_csv(fname, header=0)
    if seq.shape[1] != 4:
        raise ValueError('The file should be a CSV file with headers with exactly 4 columns '
                         '(a, b, m, n) with electrode numbers.')
    else:
        return seq.values

#%% geometrical functions 
def bearing(dx,dy):
    if dx == 0 and dy == 0:
        raise ValueError('both dx and dy equal 0 - check no 2 electrodes occupy same xy coordinates')
    elif dx == 0 and dy > 0:
        return 0
    elif dx == 0 and dy < 0:
        return 180
    elif dx > 0 and dy == 0:
        return 90
    elif dx < 0 and dy == 0:
        return 270
    elif dx > 0 and dy > 0: 
        return np.rad2deg(np.arctan(dx/dy))
    elif dx > 0 and dy < 0: 
        return 180 + np.rad2deg(np.arctan(dx/dy))
    elif dx < 0 and dy < 0: 
        return 180 + np.rad2deg(np.arctan(dx/dy))
    elif dx < 0 and dy > 0: 
        return 360 + np.rad2deg(np.arctan(dx/dy))

# fit an angle to a line in the XY plane 
def fitXYangle(x,y):
    m,c = polyfit(x,y,1)
    xmdl = np.array([np.min(x),np.max(x)]) 
    ymdl = (m*xmdl) + c 
    dx = xmdl[0] - xmdl[1]
    dy = ymdl[0] - ymdl[1]
    a = bearing(dx,dy) + 90 
    if a > 360:
        a -= 360 
    return a 

# find min and max of an array with some padding
def findminmax(a, pad=1):
    a = np.array(a)
    delta = abs(np.max(a) - np.min(a))
    mina = np.min(a) - (delta*(pad/100))
    maxa = np.max(a) + (delta*(pad/100))
    
    if mina == maxa:
        mina -= 1
        maxa += 1
    return [mina, maxa]
    

#%% main Project class (called 'R2' in previous versions)
class Project(object): # Project master class instanciated by the GUI
    """Master class to handle all processing around the inversion codes.

    Parameters
    ----------
    dirname : str, optional
        Path of the working directory. Can also be set using `R2.setwd()`.
    typ : str, optional
        Either `R2` or `R3t` for 3D. Complex equivalents are `cR2` and `cR3t`.
        Automatically infered when creating the survey.
    """
    def __init__(self, dirname='', typ='R2'): # initiate R2 class
        self.apiPath = os.path.dirname(os.path.abspath(__file__)) # directory of the code
        if dirname == '':
            dirname = os.path.join(self.apiPath)
        else:
            dirname = os.path.abspath(dirname)
        print('Working directory is:', dirname)
        self.pinfo = {} # project information 
        self.setwd(dirname) # working directory (for the datas)
        self.elec = None # will be assigned when creating a survey
        self.surveys = [] # list of survey object
        self.surveysInfo = [] # info about surveys (date)
        self.bigSurvey = None # big combined survey created in the case of timelapse mode 
        self.mesh = None # mesh object (one per Project instance)
        self.meshParams = {} # mesh parameters passed to mesh creation scheme
        self.wholespace = False # flag for whole space problem 
        self.topo = pd.DataFrame(columns=['x','y','z']) # store additional topo points
        self.coordLocal = False # flag is local coordinate conversion required 
        self.coordParam = {'x0':None, 'y0':None, 'a':None} # coordinate conversion parameters 
        self.param = {} # dict configuration variables for inversion
        self.configFile = ''
        self.invLog = '' # to save inversion output - all R2.out files
        self.fwdLog = '' # to save forward modeling R2_forward.out files
        self.typ = typ # or cR2 or R3t, cR3
        self.err = False # if we want error in protocol.dat or not
        self.iBorehole = False # to tell the software to not plot pseudoSection
        self.iTimeLapse = False # to enable timelapse inversion
        self.iBatch = False # to enable batch inversion
        self.meshResults = [] # contains vtk mesh object of inverted section
        self.projs = [] # contains Instances of Project for pseudo 3D inversion from 2D lines
        self.projectPseudo3D = None # updates iteratively - for showing pseudo 3D inversion iterations, killing, etc.
        self.pseudo3DBreakFlag = False # flag to cancel inversions in a chain (pseudo3D only)
        self.pseudo3DSurvey = None # contains one survey instance with all 2D lines combined in a 3D grid
        self.pseudo3DMeshResult = None # contains pseudo 3D mesh result (for external use - e.g., ParaView)
        self.pseudo3DMeshResultList = None # contains pseudo 3D meshes (individual 2D lines in 3D space - for external use - e.g., ParaView)
        self.pseudo3Dfmd = None # pseudo3D project FMD (only used when user defines a FMD in createMultiMesh())
        self.sequence = None # quadrupoles sequence if forward model
        self.resist0 = None # initial resistivity
        self.iForward = False # if True, it will use the output of the forward
        # to run an inversion (and so need to reset the regions before this)
        self.fmd = None # depth of investigation below the surface [in survey units]
        self.proc = None # where the process to run R2/cR2 will be
        self.mproc = None # where the process for mesh building
        self.zlim = None # zlim to plot the mesh by default (from max(elec, topo) to min(doi, elec))
        self.trapeziod = None # trapeziod vertices of cropped 2D mesh (triangles removed from bottom corners)
        self.geom_input = {} # dictionnary used to create the mesh
        # attributes needed for independant error model for timelapse/batch inversion
        self.referenceMdl = False # is there a starting reference model already?
        self.fwdErrModel = False # is there is a impoforward modelling error already (due to the mesh)?
        self.custSeq = False # flag - True if of 3D custom sequence imported
        self.errTyp = 'global'# type of error model to be used in batch and timelapse surveys
        self.surfaceIdx = None # used to show plan view iterations of 3D inversions
        self.darkMode = False # If true, electrodes wil be plotted in white, else black
        self.iadvanced = True # If true, use the advanced mesh format for 3D mesh
        # summary information, shown went printing 
        self.pinfo = {'Data':False,
                      'Number of Surveys': 0,
                      'Number of electrodes':0,
                      'Borehole electrodes':False, 
                      'Remote electrodes':False, 
                      'Error model':False, 
                      'Mesh Type':'Not set',
                      'Number of Elements':0,
                      'Number of Nodes':0,
                      'Estimated RAM for forward solution (Gb)':0.0,
                      'Estimated RAM for inverse solution (Gb)':0.0,
                      'Convergence of inverse solution':False, 
                      'RMS of inverse solution':np.nan,
                      'RMS error as percentage estimate':np.nan, 
                      'Number of iterations':0,
                      'Median normalised error':np.nan,
                      'Working Directory':self.dirname}
        
        
    def summary(self): 
        text = '_____Project Summary_____\n'
        for key in self.pinfo.keys(): 
            text += key 
            text += ' : '
            v = self.pinfo[key] # dict value 
            if type(v) is int: 
                text += '{:d}'.format(v)
            elif type(v) is float:
                text += '{:6.4f}'.format(v)
            elif type(v) is bool: 
                if v:
                    text += 'Yes'
                else:
                    text += 'No'
            elif type(v) is str:
                text += '{:s}'.format(v)
            else:
                text += 'nan' 
            text += '\n'
            
        return text 
    
            
    def __str__(self):
        return self.summary() 
                    
    
    def setBorehole(self, val=False):
        """Set all surveys in borehole type if `True` is passed.
        """
        self.iBorehole = val
        for s in self.surveys:
            s.iBorehole = val
        if self.iBorehole: 
            self.pinfo['Borehole electrodes'] = True 


    def _num2elec(self, elec):
        """Return a formated dataframe with 'x','y','z','label','remote','buried'
        columns wether the input is a sparse matrix or dataframe.
        """
        if type(elec) == np.ndarray:
            if elec.shape[1] == 1:
                elec = pd.DataFrame(elec, columns=['x'])
            elif elec.shape[1] == 2:
                elec = pd.DataFrame(elec, columns=['x','z'])
            elif elec.shape[1] == 3:
                elec = pd.DataFrame(elec, columns=['x','y','z'])
            elif elec.shape[1] == 4:
                elec = pd.DataFrame(elec, columns=['x','y','z','buried'])
                elec.loc[:,'buried'] = elec['buried'].astype(bool)

        if 'line' in elec.columns and 'number' in elec.columns:                
            elecid  = elec.number.values # electrode number array 
            line = elec.line.values 
            label = ['1 1'] * len(elec)
            for i in range(len(elec)):
                label[i] = '%i %i'%(line[i],elecid[i])
            if 'label' in elec.columns:
                elec.drop(columns=['label'], inplace=True)
            elec.insert(0,'label',label)
        if 'y' not in elec.columns:
            elec.insert(len(elec.columns),'y',0)
            # elec['y'] = 0  # all on same line by default
        if 'z' not in elec.columns:
            elec.insert(len(elec.columns),'z',0)
            # elec['z'] = 0 # all flat topo by default
        if 'remote' not in elec.columns:
            elec.insert(len(elec.columns),'remote',False)
            # elec['remote'] = False # all non remote by default
        if 'buried' not in elec.columns:
            elec.insert(len(elec.columns),'buried',False)
            # elec['buried'] = False # all surface elec by default
        if 'label' not in elec.columns:
            label = (1 + np.arange(elec.shape[0])).astype(str) 
            elec.insert(0,'label',label)# all elec ordered and start at 1            
        elec = elec.astype({'x':float, 'y':float, 'z':float, 'buried':bool, 'remote':bool, 'label':str})
        
        return elec
    
    
    def _findRemote(self, elec):
        """Flag remote electrodes.
        
        Parameters
        ----------
        elec : dataframe
            contains the electrodes information
            
        Returns
        -------
        elec : dataframe
            remote flag added to electrodes
        """
        remote_flags = [-9999999, -999999, -99999,-9999,-999]
                    #9999999, 999999, 99999, 9999, 999] # values asssociated with remote electrodes
        iremote = np.in1d(elec['x'].values, remote_flags)
        iremote = np.isinf(elec[['x','y','z']].values).any(1) | iremote
        elec['remote'] = iremote
        if np.sum(iremote) > 0:
            print('Detected {:d} remote electrode.'.format(np.sum(iremote)))
            
        return elec
    

    def setCoordConv(self,flag=False,x0=None,y0=None,a=None):
        """
        Set that the project should be imported and exported according to a 
        coordinate conversion system. Generally UTM coordinates will cause 
        stability issues with mesh generation and the finite element solutions 
        used in the R* codes, hence the we use a local coordinate system. 
        Call this function before creating the mesh. 

        Parameters
        ----------
        flag : bool, optional
            Flag if coordinate conversion system is required. 
            The default is False. If flag is false but Project.coordLocal is 
            True then the function can be used to rvert electrode coordinates 
            back to their initial position.
        x0 : float, optional
            Reference X coordinate in UTM system. The default is None. If left 
            as none, then the minimum X coordinate of the electrodes will be 
            used. 
        y0 : float, optional
            Reference Y coordinate in UTM system. The default is None. If left 
            as none, then the minimum Y coordinate of the electrodes will be 
            used. 
        a : float, optional
            Rotation angle in degrees. The default is 0.
        """
        if flag and self.coordLocal:
            warnings.warn('Coordinate conversion has already been set! Are you sure want to try this again?')
            return 
        revert = False 
        if not flag and self.coordLocal:
            warnings.warn('reverting coordinates back to their original positions')
            revert = True 
            
        self.coordLocal = flag 
        
         #check coordParam is in the right format / data type 
        if self.coordParam is None:
            self.coordParam = {'x0':None, 'y0':None, 'a':None} 
        if len(self.coordParam.keys()) < 3:
            self.coordParam = {'x0':None, 'y0':None, 'a':None} 
            
        # update coordParam if conversion parameters are given 
        if x0 is not None:
            self.coordParam['x0'] = x0 
        if y0 is not None: 
            self.coordParam['y0'] = y0 
        if a is not None:
            self.coordParam['a'] = a 
            
        # if electrodes assigned to project, check if they need moving 
        if self.elec is not None: 
            # get some default conversion parameters if not provided 
            if self.coordParam['x0'] is None and x0 is None: 
                x0 = np.min(self.elec.x.values)
                self.coordParam['x0'] = x0 
                # print('Auto setting x0')
            elif self.coordParam['x0'] is not None: 
                x0 = self.coordParam['x0'] 
                
            if self.coordParam['y0'] is None and y0 is None: 
                y0 = np.min(self.elec.y.values)
                self.coordParam['y0'] = y0 
                # print('Auto setting y0')
            elif self.coordParam['y0'] is not None: 
                y0 = self.coordParam['y0'] 
                
            if self.coordParam['a'] is None and a is None:  
                # print('Auto fitting A')
                if '2' in self.typ and any(self.elec.y != self.elec.y[0]):
                    a = fitXYangle(self.elec.x.values, self.elec.y.values)
                else: 
                    a = 0 
                self.coordParam['a'] = a 
            elif self.coordParam['a'] is not None: 
                a = self.coordParam['a'] 
                
            if self.coordLocal: # move electrodes is local coordinates required 
                localx, localy = rotGridData(self.elec.x.values, self.elec.y.values, x0, y0, a)
                elec = self.elec.copy()
                elec.loc[:,'x'] = localx 
                elec.loc[:,'y'] = localy 
                self.elec = elec 
                for s in self.surveys: 
                    localx, localy = rotGridData(s.elec.x.values, s.elec.y.values, x0, y0, a)
                    s.elec.loc[:,'x'] = localx 
                    s.elec.loc[:,'y'] = localy 
                    s.computeK() 
            elif revert: # revert coordinates back using inverse rotation 
                utmx, utmy = invRotGridData(self.elec.x.values, self.elec.y.values, x0, y0, a)
                elec = self.elec.copy()
                elec.loc[:,'x'] = utmx
                elec.loc[:,'y'] = utmy 
                self.elec = elec 
                for s in self.surveys: 
                    utmx, utmy = invRotGridData(s.elec.x.values, s.elec.y.values, x0, y0, a)
                    s.elec.loc[:,'x'] = utmx 
                    s.elec.loc[:,'y'] = utmy 
                    # s.computeK() # doubt there would be a need to recompute k 
                
    
    def setElec(self, elec, elecList=None, _uiOverideCoordLocal=False):
        """Set electrodes. Automatically identified remote electrode.

        Parameters
        ----------
        elec : numpy array
            Array of NxM dimensions. N = number of electrodes, M = 2 for x,z or
            M = 3 if x,y,z coordinates are supplied.
        elecList : list, optional
            If not None then elec is ignored in favor of elecList. This option
            is to be used in the advanced use case where electrodes move which
            each survey. Each entry of the list is a numpy array the same format
            of 'elec', which is then assigned to each survey class.
        _uiOverideCoordLocal: bool, optional
            UI command only, set to to True when updating the project electrodes
            in the UI. 
        """
        if elecList is None: # same electrode set shared by all surveys (most common case)
            ok = False
            elec = self._num2elec(elec)
            if self.elec is not None: # electrode already inferred when parsing data
                if self.iForward: # in case of forward modelling, changing the number of electrodes is allowed
                    ok = True
                else: # check intersection of labels
                    if len(self.surveys) > 0:
                        s = self.surveys[0]
                        if s.hasElecString():
                            s1 = np.unique(elec['label'].values)
                            s2 = np.unique(s.df[['a','b','m','n']].values.flatten())
                        else:
                            s1 = np.unique(elec['label'].astype(float).astype(int).values)
                            s2 = np.unique(s.df[['a','b','m','n']].astype(float).astype(int).values)
                        x = np.intersect1d(s1, s2)
                        if len(x) == len(s2):
                            ok = True
                        else:
                            raise ValueError('The following electrode labels are missing'
                                  ' from the electrode declaration: ' + ', '.join(s2[~np.in1d(s2, x)]))
                    else:
                        if elec.shape[0] >= self.elec.shape[0]:
                            ok = True
                        else:
                            raise ValueError('The number of electrodes read ({:d}) is smaller'
                              ' than the number of electrode from data file ({:d})'.format(
                                  elec.shape[0], self.elec.shape[0]))                
            else:
                ok = True # first assignement of electrodes
            if ok:
                elec = self._findRemote(elec) # identification of remote electrode
                self.elec = elec 
                
                # for some reason these columns are sometimes missing in timelapse mode?? 
                if 'buried' not in elec.columns: 
                    elec['buried'] = False 
                if 'remote' not in elec.columns:
                    elec['remote'] = False 
                    
                for s in self.surveys:
                    s.elec = elec
            
            # check for shared electrodes (between lines) or duplicated
            self._checkElecDuplicate()
        
        else:
            #some error checking
            try:
                num_surveys = len(self.surveys) # number of surveys
                if len(elecList) != num_surveys:
                    raise ValueError("The number of electrode matrices must match the number of surveys")
            except AttributeError:
                raise AttributeError("No Survey attribute assocaited with R2 class, make sure you create a survey first")
            
            if self.elec is None:
                initElec = elecList[0]
                self.elec = self._num2elec(np.zeros((len(initElec),3)))
            for i, survey in enumerate(self.surveys):
                survey.elec = self._findRemote(self._num2elec(elecList[i])) # plus identification of remote electrode

        # do coordinate conversion if local grid required 
        if not _uiOverideCoordLocal: 
            if self.coordLocal and self.elec is not None:
                self.coordLocal = False # will be reset to true in the following function 
                self.setCoordConv(True)

        if len(self.surveys) > 0:
            self.computeFineMeshDepth()
            for s in self.surveys:
                s.setSeqIds()
                s.computeK() # recalculate K 
                
        self.pinfo['Number of electrodes'] = len(self.elec)
        if 'iremote' in self.elec.columns: 
            if any(elec['iremote']):
                self.pinfo['Remote electrodes'] = True 
                
    
    def mergeElec(self, dist=-1):
        """Merge electrodes that have less than a certain distance to eache other
        
        Parameters
        ----------
        dist : float, optional
            maximum distance of close electrodes in meters (electrodes that have a distance less than dist will be merged)
            -1 flag for auto calculation where dist = average(electrode spacing)/100
            
        Returns
        -------
        bool
            True: merge successful
            False: merge unsuccessful
        """
        
        from scipy.spatial import distance
        
        eleccoors = self.elec[['x','y','z']].values.copy()
        
        if dist == -1:
            dist = np.mean(distance.cdist(eleccoors, eleccoors))/100
        
        c = 0
        closepairsL = [] # empty list to get len(closepairsL) == 0 in loop's first step (i=0)
        for i in range(50):
            elecdists = distance.cdist(eleccoors, eleccoors)
            elecdists[np.triu_indices(len(elecdists))] = dist
            elecdists[elecdists == 0] = dist
            
            # indices of pairs with distance < dist
            closepairs = np.c_[np.where(elecdists < dist)[0], np.where(elecdists < dist)[1]]
            
            if len(closepairs) > len(closepairsL): 
                c += 1
            else:
                c = 0
            if c >= 4:
                break
            
            poses = np.zeros_like(closepairs[:,0])
            for i in range(len(poses)):
                poses[i] = closepairs[i,0] if closepairs[i,0] < closepairs[i,1] else closepairs[i,1]
            
            eleccoors[closepairs[:,0], :] = eleccoors[poses, :]
            closepairsL = closepairs.copy()
            elecdf_temp = self.elec.copy() # to make it compatible with the new self.setElec
            elecdf_temp[['x','y','z']] = eleccoors
            
            if len(closepairs) == 0:
                print('Merging close electrodes successful!')
                self.setElec(elecdf_temp)
                return True
        
        if len(closepairs) != 0: # solution did not converge
            print('Merging close electrodes unsuccessful - Choose a shorter merging distance (i.e., dist < {:.2f}m)'.format(dist))
            return False
    
    
    def _checkElecDuplicate(self):
        """Check for duplicated electrodes positions and merge them
        to the same label in this case. Useful for 3D survey from 2D lines,
        if multiple lines shared have electrodes in common.
        """
        # from: https://stackoverflow.com/questions/46629518/find-indices-of-duplicate-rows-in-pandas-dataframe
        elec = self.elec.copy()
        df = elec[['x','y','z']].copy()
        iunique = df.duplicated(keep=False)
        if np.sum(iunique) > 0:
            print('Merging electrodes positionned at the same location.')
            df = df[iunique]
            dups = df.groupby(list(df)).apply(lambda x: tuple(x.index)).tolist()
            
            for dup in dups:
                dico = {}
                label2keep = elec.loc[dup[0], 'label']
                for i, index in enumerate(dup[1:]):
                    label = elec.loc[index,'label'] # label of duplicated electrodes
                    dico[label] = label2keep
                    # check through elec for each survey                        
                    for survey in self.surveys:
                        if label2keep in survey.elec['label'].values:
                            # if the label2keep is present, just remove the duplicate
                             i2keep = survey.elec['label'] != label
                             survey.elec = survey.elec[i2keep].reset_index(drop=True)
                        else: # if not just renamed its label
                            survey.elec = survey.elec.replace(to_replace=dico)
                # and we replace (not remove) the label by the label of the first 
                # occurence of the duplicated electrode
                for survey in self.surveys:
                    survey.df = survey.df.replace(to_replace=dico)
                    
                # finally we remove duplicates from the Project instance itself
                if label2keep in self.elec['label'].values:
                    # if the label2keep is present, just remove the duplicate
                     i2keep = self.elec['label'] != label
                     self.elec = self.elec[i2keep].reset_index(drop=True)
                else: # if not just renamed its label
                    self.elec = self.elec.replace(to_replace=dico)
            
            for survey in self.surveys: 
                survey.setSeqIds() 
                
                
            
    def generateElec(self, nb=24, dx=0.5, dz=0, nline=1, lineSpacing=2):
        """Generate electrodes positions for 2D and 3D surveys.
        
        Parameters
        ----------
        nb : int, optional
            Number of electrodes per line. For 3D survey, if multiple lines,
            the total number of electrodes will be nb x nline.
        dx : float, optional
            Spacing in meters between electrodes in the X direction.
        dz : float, optional
            Increment in the Z direction (elevation) between consecutive electrodes.
        nline : int, optional
            Number of lines. 1 for 2D and multiple for 3D.
        lineSpacing : float, optional
            Spacing between lines (3D only).
        """
        # check type
        nb = int(nb)
        nline = int(nline)
        
        # electrode positions for one line
        elec = np.zeros((nb, 3))
        elec[:,0] = np.linspace(0, (nb-1)*dx, nb)
        elec[:,2] = np.linspace(0, (nb-1)*dz, nb)
        
        # specific 2D or 3D additions and building dfelec
        if (self.typ == 'R2') | (self.typ == 'cR2'):
            dfelec = pd.DataFrame(columns=['label','x','y','z','buried'])
            dfelec['label'] = (1 + np.arange(nb)).astype(int).astype(str)
            dfelec[['x','y','z']] = elec
            dfelec['buried'] = False
        else:
            grid = np.tile(elec.T, nline).T
            labels = []
            for i in range(nline):
                labels.append(['{:d} {:d}'.format(i+1, j+1) for j in range(nb)]) # string + elec number
                grid[i*nb:(i+1)*nb, 1] = (i+1)*lineSpacing # Y is line spacing
            dfelec = pd.DataFrame(columns=['label','x','y','z','buried'])
            dfelec['label'] = np.hstack(labels)
            dfelec['x'] = grid[:,0]
            dfelec['y'] = grid[:,1]
            dfelec['z'] = grid[:,2]
            dfelec['buried'] = False
            
        self.setElec(dfelec)
        
    def hasElecString(self):
        """Determine if a electrode strings are present in the electrode labels 

        Returns
        -------
        bool
            True if strings present in electrode label
        """
        df = self.elec
        if 'label' not in df.keys():
            return False
        else:
            label = df['label']
            for l in label:
                if len(l.split()) == 1:
                    return False
        return True
        
        
    def detectStrings(self, tolerance=15, max_itr=None):
        """Automatically detect electrode strings. If all electrodes lie on the same
        Y plane then the function assumes the problem is a borehole problem. 

        Parameters
        ----------
        tolerance : float, optional
            Maximum (+/-) bearing (ie directional angle) tolerance each subsiquent
            electrode may have. The default is 15.
        max_itr : int, optional
            Maximum number of searches that can be performed to find colinear
            nieghbouring electrodes. The default is number of electrodes plus 10.

        Raises
        ------
        ValueError
            if the change in x and y direction for 2 neighbouring electrodes is 
            the same. ie no 2 electrodes can occupy the same x y position in 
            this code. 
        ValueError
            if the change maxium number of searches is exceeded. 

        Returns
        -------
        list (of list)
            Each entry in the list corresponds to an electrode string, and is
            a list of integers which are the indices of the respective electrodes
            in self.elec
        """
        iremote = self.elec['remote'].values
        x = self.elec['x'].values[~iremote]
        y = self.elec['y'].values[~iremote]
        
        # deal with borehole case 
        if all(y[0] == y):
            print('Detecting strings for a borehole problem!')
            y = self.elec['z'].values[~iremote]
            
        # decide maximum number of iterations 
        if max_itr is None:
            max_itr = len(x)+10
        
        # mid point of survey 
        xm = np.mean(x)
        ym = np.mean(y)
        dist = np.sqrt((x-xm)**2 + (y-ym)**2) # distance to all other points in survey
        
        # find neighbours 
        points = np.array([x,y]).T
        tree = cKDTree(points)
        distm, niegh = tree.query(points,k=4)
        
        #start finding strings 
        string = [] # caches for saving found electrodes
        allocated = np.zeros(len(x),dtype=bool) 
        
        c = 0 
        while any(allocated==False):
            #color = tuple(np.random.rand(3)) # color for line 
            #find starting point
            si = np.argmax(dist) # start point index
            # ax.scatter(x[si],y[si],c='r')
            allocated[si] = True
            dist[si] = -1
            this_string = [si]
            #get start point neighbours 
            n = niegh[si,1:4]
            dxdy = points[n,:] - points[si,:] #deltas 
            d = distm[si,1:4]
            di = np.argmin(d) # index of min distance
            ni = n[di] # index of nearest neighbour
            this_string.append(ni) # save
            allocated[ni] = True
            dist[ni] = -1
            
            if (dxdy[di,0] == 0) and (dxdy[di,1] == 0): # fail safe. 
                # cant determine bearing direction in this case, therefore just 
                # assume all electrodes lie on one string
                this_string = [i+1 for i in range(len(x))]
                string.append(this_string)
                break
            
            #work out staring direction 
            theta = bearing(dxdy[di,0],dxdy[di,1])
            # ax.plot([x[si],x[ni]],[y[si],y[ni]],c=color)
            canidate = True
            
            while canidate:
                pi = ni # prevoius neighbour 
                n = niegh[ni,1:4]
                dxdy = points[n,:] - points[ni,:] 
                d = distm[ni,1:4]
                itr = len(n)
                thetav = np.zeros(itr)
                for i in range(len(n)):
                    thetav[i] = bearing(dxdy[i,0],dxdy[i,1])
                thtest = np.abs(thetav - theta) < tolerance # test if another point on desired bearing
                if all(thtest==False):
                    canidate = False # no more canidates 
                else:
                    di = np.argmin(d[thtest]) # index of min distance
                    ni = n[thtest][di] # index of next nearest neighbour
                    theta = bearing(x[ni]-x[pi],y[ni]-y[pi])
                    # ax.plot([x[pi],x[ni]],[y[pi],y[ni]],c=color)
                    allocated[ni] = True
                    this_string.append(ni)
                    dist[ni] = -1
                    #overwrite electrode string in data frame? #TODO
                    # label = self.elec.loc[ni,'label'].split() 
                c+=1
                if c>max_itr:
                    canidate=False
            if c>max_itr:
                raise ValueError('max number of iterations exceeded')
                break
            string.append(this_string)
            
        return string
    
    def elec2horidist(self):
        """
        Convert 2D xz data into a true horizontal distance. Assumes that survey 
        was done with a tape measure and the X distances are not true horizontal
        distance but distances measured along the ground. 

        """
        for s in self.surveys:
            s.elec2horidist() 
            
        # reset first survey electrodes to project electrodes 
        self.setElec(self.surveys[0].elec)


    def convertLocalGrid(self):
        """
        Converts UTM grid to local grid for stability in mesh generation and 
        the finite element solution. 
        
        Parameters
        ----------
        reset : bool, optional
            reset the grid to original
        """
        warnings.warn('This function is now redundant, use setCoordConv() which is more robust',
                      DeprecationWarning)
        for s in self.surveys:
            s.convertLocalGrid() 

        # reset first survey electrodes to project electrodes
        self.setElec(self.surveys[0].elec)


    def setwd(self, dirname):
        """Set the working directory.

        Parameters
        ----------
        dirname : str
            Path of the working directory.
        """
        wd = os.path.abspath(os.path.join(dirname, 'invdir'))
        if os.path.exists(dirname) is False:
            os.mkdir(dirname)
        if os.path.exists(wd):
            shutil.rmtree(wd)
            print('clearing dirname')
        os.mkdir(wd)
        self.dirname = wd
        self.pinfo['Working Directory'] = wd 


    def setTitle(self, linetitle):
        """Set the title of the survey name when inverting data. Input is a string.
        """
        if isinstance(linetitle, str):
            self.param['lineTitle'] = linetitle
        else:
            print("Cannot set Survey title as input is not a string")


    def saveProject(self, fname):
        """Save the current project will all dataset in custom 
        ResIPy format (.resipy) for future importation.
        """
        from zipfile import ZipFile, ZipInfo
        import json
        
        if fname[-7:] != '.resipy':
            fname = fname + '.resipy'
        
        # create save directory
        name = os.path.basename(fname)
        savedir = os.path.join(self.dirname, name)
        if os.path.exists(savedir):
            shutil.rmtree(savedir)
        os.mkdir(savedir)
        
        # add files to it
        if self.mesh is not None and self.pseudo3DSurvey is None:
            self.mesh.vtk(os.path.join(savedir, 'mesh.vtk'))
        self.elec.to_csv(os.path.join(savedir, 'elec.csv'), index=False)
        c = 0
        if self.iForward:
            c = 1
            dfseq = pd.DataFrame(self.sequence, columns=['a','b','m','n'])
            dfseq.to_csv(os.path.join(savedir, 'dfseq.csv'), index=False)
        for i, survey in enumerate(self.surveys):
            f = os.path.join(savedir, survey.name)
            survey.df.to_csv(f + '-df.csv', index=False)
            survey.elec.to_csv(f + '-elec.csv', index=False)
            if (c+i < len(self.meshResults)):
                self.meshResults[c + i].vtk(f + '.vtk')
        
        # batch/time-lapse
        if self.iBatch or self.iTimeLapse:
            f = os.path.join(savedir, 'bigSurvey')
            self.bigSurvey.df.to_csv(f + '-bigdf.csv', index=False)
            
        # pseudo 3D
        if self.pseudo3DSurvey is not None:
            f = os.path.join(savedir, 'pseudo3DSurvey')
            self.pseudo3DSurvey.df.to_csv(f + '-pseudo3Ddf.csv', index=False)
            self.pseudo3DSurvey.elec.to_csv(f + '-pseudo3Delec.csv', index=False)
            for proj in self.projs:
                if proj.mesh is not None:
                    name = proj.surveys[0].name
                    proj.mesh.vtk(os.path.join(savedir, '{}-ps3d.vtk'.format(name)))
        
        settings = {'surveysInfo': self.surveysInfo,
                    'topo': self.topo.to_dict(),
                    'coordLocal':self.coordLocal, 
                    'coordParam':self.coordParam, 
                    'typ': self.typ,
                    'err': self.err,
                    'iBorehole': self.iBorehole,
                    'iTimeLapse': self.iTimeLapse,
                    'iBatch': self.iBatch,
                    'sequence': self.sequence.tolist() if self.sequence is not None else None,
                    'resist0': list(self.resist0) if self.resist0 is not None else None,
                    'iForward': self.iForward,
                    'fmd': self.fmd,
                    'zlim': self.zlim,
                    'geom_input': self.geom_input, # need to make list of array?
                    'referenceMdl': self.referenceMdl,
                    'fwdErrModel': self.fwdErrModel,
                    'custSeq': self.custSeq,
                    'errTyp': self.errTyp,
                    'surfaceIdx': list(self.surfaceIdx) if self.surfaceIdx is not None else None
                   }
        with open(os.path.join(savedir, 'settings.json'), 'w') as f:
            f.write(json.dumps(settings))
            
        # also save summary information 
        with open(os.path.join(savedir,'summary.json'), 'w') as f:
            f.write(json.dumps(self.pinfo))
        
        # param as numpy array
        keys = ['num_regions', 'res0File', 'num_xz_poly', 'a_wgt', 'b_wgt',
              'lineTitle', 'job_type', 'flux_type', 'singular_type',
              'res_matrix', 'scale', 'regions', 'patch_x', 'patch_z',
              'inverse_type', 'target_decrease', 'qual_ratio', 'data_type',
              'reg_mode', 'tolerance', 'max_iter', 'error_mod', 'alpha_aniso',
              'alpha_s', 'min_error', 'rho_min', 'rho_max', 'mesh_type']
        sparams = {}
        for key in keys:
            if key in self.param.keys():
                sparams[key] = self.param[key]
        if 'node_elec' in self.param:
            sparams['node_elec'] = [list(self.param['node_elec'][0]),
                                    [int(a) for a in self.param['node_elec'][1]]]
            # int64 not JSON serializable so we convert it to int (int32)
        if 'xz_poly_table' in self.param:
            sparams['xz_poly_table'] = self.param['xz_poly_table'].tolist()
        if 'xy_poly_table' in self.param:
            sparams['xy_poly_table'] = self.param['xy_poly_table'].tolist()
        with open(os.path.join(savedir, 'params.json'), 'w') as f:
            f.write(json.dumps(sparams))
        
        with open(os.path.join(savedir, 'invLog.log'), 'w') as f:
            f.write(self.invLog)
        
        with open(os.path.join(savedir, 'fwdLog.log'), 'w') as f:
            f.write(self.fwdLog)
        
        # zip the directory, move it and clean
        with ZipFile(fname, 'w') as fz:
            for file in os.listdir(savedir):
                fz.write(os.path.join(savedir, file), file)
        shutil.rmtree(savedir)
        
        # TODO maybe add a self.uiParams = {} for UI specific parameters?
     
        
    def loadProject(self, fname):
        """Load data from project file.
        
        Parameters
        ----------
        fname : str
            Path where the file will be saved.
        """
        from zipfile import ZipFile, ZipInfo
        import json
        
        # create save directory
        name = os.path.basename(fname).replace('.resipy','')
        savedir = os.path.join(self.dirname, name)
        if os.path.exists(savedir):
            shutil.rmtree(savedir)
        os.mkdir(savedir)
        
        # read in zip and extract in working directory
        with ZipFile(fname, 'r') as fz:
            fz.extractall(savedir)
            
        # read files an reconstruct Survey objects
        self.meshResults = []
        self.surveys = []
        dfnames = glob.glob(os.path.join(savedir,'*-df.csv'))
        elecnames = glob.glob(os.path.join(savedir,'*-elec.csv'))
        for dfname, elecname in zip(dfnames, elecnames): # better than 100000 loops!
            df = pd.read_csv(dfname)
            for c in ['a','b','m','n']:
                df[c] = df[c].astype(str)
            dfelec = pd.read_csv(elecname)
            dfelec['label'] = dfelec['label'].astype(str)
            surveyName = os.path.basename(dfname[:-7])
            self.surveys.append(Survey(df=df, elec=dfelec, name=surveyName))
            elec = dfelec[~dfelec['remote']][['x','y','z']].values
            if os.path.exists(os.path.join(savedir, surveyName + '.vtk')):
                mesh = mt.vtk_import(os.path.join(savedir, surveyName + '.vtk'))
                mesh.setElec(elec[:,0], elec[:,1], elec[:,2])
                self.meshResults.append(mesh)
            self.surveys[-1].setSeqIds()

        # for i in range(100000): # don't think that someone will store so many
        #     f = os.path.join(savedir, 'survey{:d}'.format(i))
        #     if os.path.exists(f + '-df.csv'):
        #         df = pd.read_csv(f + '-df.csv')
        #         for c in ['a','b','m','n']:
        #             df[c] = df[c].astype(str)
        #         dfelec = pd.read_csv(f + '-elec.csv')
        #         dfelec['label'] = dfelec['label'].astype(str)
        #         self.surveys.append(Survey(df=df, elec=dfelec)) 
        #         elec = dfelec[~dfelec['remote']][['x','y','z']].values
        #         if os.path.exists(f + '.vtk'):
        #             mesh = mt.vtk_import(f + '.vtk')
        #             mesh.setElec(elec[:,0], elec[:,1], elec[:,2])
        #             self.meshResults.append(mesh)
        
        # batch/time-lapse
        fbig = os.path.join(savedir, 'bigSurvey')
        if os.path.exists(fbig + '-bigdf.csv'):
            df = pd.read_csv(fbig + '-bigdf.csv')
            for c in ['a','b','m','n']:
                df[c] = df[c].astype(str)
            self.bigSurvey = Survey(df=self.surveys[0].df, elec=self.surveys[0].elec)
            self.bigSurvey.df = df.copy() 
            self.bigSurvey.dfOrigin = df.copy()

        self.elec = pd.read_csv(os.path.join(savedir, 'elec.csv'))
        self.elec['label'] = self.elec['label'].astype(str)
        if os.path.exists(os.path.join(savedir, 'mesh.vtk')):
            self.importMesh(os.path.join(savedir, 'mesh.vtk')) # better
#         self.mesh = mt.vtk_import(os.path.join(savedir, 'mesh.vtk'))
        
        # read flags and settings
        with open(os.path.join(savedir, 'settings.json'), 'r') as f:
            settings = json.load(f)
        self.surveysInfo = settings['surveysInfo']
        self.topo = pd.DataFrame(settings['topo'])
        self.typ = settings['typ']
        self.err = settings['err']
        self.iBorehole = settings['iBorehole']
        self.iTimeLapse = settings['iTimeLapse']
        self.iBatch = settings['iBatch']
        self.sequence = np.array(settings['sequence']) if settings['sequence'] is not None else None
        self.resist0 = settings['resist0']
        self.iForward = settings['iForward']
        if self.iForward:
            if os.path.exists(os.path.join(savedir, 'dfseq.csv')):
                self.importSequence(os.path.join(savedir, 'dfseq.csv'))
            shutil.copytree(savedir, os.path.join(self.dirname, 'fwd'))
            shutil.move(os.path.join(self.dirname, 'fwd', 'mesh.vtk'),# needed for inverting a fwd_only project after loading
                        os.path.join(self.dirname, 'fwd', 'forward_model.vtk'))
        self.fmd = settings['fmd']
        self.zlim = settings['zlim']
        if self.iForward and self.mesh is not None:
            self.meshResults = [self.mesh] + self.meshResults
        try:
            self.geom_input = settings['geom_input'] # might failed if numpy array inside
        except:
            print('could not import geom_input')
        self.referenceMdl = settings['referenceMdl']
        self.fwdErrModel = settings['fwdErrModel']
        self.custSeq = settings['custSeq']
        self.errTyp = settings['errTyp']
        self.surfaceIdx = settings['surfaceIdx']
        # check if local coordinate flag in keys (need to do this for backward compatibility)
        if 'coordLocal' in settings.keys(): 
            self.coordLocal = settings['coordLocal'] 
        else: 
            self.coordLocal = False 
            
        if 'coordParam' in settings.keys():
            self.coordParam = settings['coordParam']
        else: 
            self.coordParam = {'x0':None,'y0':None,'a':0}
            
        # read parameters
        with open(os.path.join(savedir, 'params.json'), 'r') as f:
            sparams = json.load(f)
        if 'xz_poly_table' in sparams:
            sparams['xz_poly_table'] = np.array(sparams['xz_poly_table'])
        if 'xy_poly_table' in sparams:
            sparams['xy_poly_table'] = np.array(sparams['xy_poly_table'])
        sparams['mesh'] = self.mesh
        if 'node_elec' in sparams:
            sparams['node_elec'][0] = np.array(sparams['node_elec'][0])
            sparams['node_elec'][1] = np.array(sparams['node_elec'][1]).astype(int)
        self.param = sparams
        
        # pseudo 3D - must be here to take params from self
        fpseudo3D = os.path.join(savedir, 'pseudo3DSurvey')
        if os.path.exists(fpseudo3D + '-pseudo3Ddf.csv'):
            df = pd.read_csv(fpseudo3D + '-pseudo3Ddf.csv')
            for c in ['a','b','m','n']:
                df[c] = df[c].astype(str)
            dfelec = pd.read_csv(fpseudo3D + '-pseudo3Delec.csv')
            dfelec['label'] = dfelec['label'].astype(str)
            self.pseudo3DSurvey = Survey(df=df, elec=dfelec)
            for s in self.surveys:
                directory = os.path.join(self.dirname, s.name)
                os.mkdir(directory) # making separate inversion diectories
                proj = self._createProjects4Pseudo3D(dirname=directory, invtyp=self.typ) # non-parallel meshing
                proj.createSurvey(fname=None, name=s.name, df=s.df, elec=s.elec)
                self.projs.append(proj) # appending projects list for later use of meshing and inversion
                if os.path.exists(os.path.join(savedir, '{}-ps3d.vtk'.format(s.name))):
                    proj.importMesh(os.path.join(savedir, '{}-ps3d.vtk'.format(s.name)))
            if self.meshResults != []:
                for proj, mresult in zip(self.projs, self.meshResults):
                    proj.meshResults.append(deepcopy(mresult))
            self._updatePseudo3DSurvey()
            self.mesh = self.projs[0].mesh # we don't want an empty self.mesh if there are meshes in self.projs
        
        # populating inversion/forward modeling logs 
        if os.path.exists(os.path.join(savedir, 'invLog.log')): # compatibility with old saves
            with open(os.path.join(savedir, 'invLog.log'), 'r') as f:
                self.invLog = f.read()
            
            with open(os.path.join(savedir, 'fwdLog.log'), 'r') as f:
                self.fwdLog = f.read()
                
        # read summary information (if it is there)
        if os.path.exists(os.path.join(savedir, 'summary.json')):
            with open(os.path.join(savedir, 'summary.json'), 'r') as f:
                self.pinfo = json.load(f)

    def createSurvey(self, fname='', ftype='Syscal', info={}, spacing=None, 
                     parser=None, debug=True, estMemory=True, **kwargs):
        """Read electrodes and quadrupoles data and return 
        a survey object.

        Parameters
        ----------
        fname : str
            Filename to be parsed.
        ftype : str, optional
            Type of file to be parsed. Either 'Syscal','ProtocolDC','ResInv',
            'BGS Prime', 'ProtocolIP', 'Sting', 'ABEM-Lund', 'Lippmann' or 'ARES'.
        info : dict, optional
            Dictionnary of info about the survey.
        spacing : float, optional
            Electrode spacing to be passed to the parser function.
        parser : function, optional
            A parser function to be passed to `Survey` constructor.
        debug : bool, optional
            If True, information about the reciprocal measurements, default 
            filtering, etc. will be displayed.
        estMemory: bool, optional
            If true, estimate the amount of RAM required to do the inversion. 
            Default is True. 
        **kwargs: Keyword arguments to be passed to Survey()
        
        """
        self.surveys.append(Survey(fname, ftype, spacing=spacing, parser=parser, debug=debug, **kwargs))
        self.surveysInfo.append(info)
        self.setBorehole(self.iBorehole)

        # if all survey have magErr and phiErr then put self.err = True
        ok = True
        for s in self.surveys:
            if 'magErr' not in s.df.columns:
                ok = False
        if ok is True:
            print('magErr/phiErr columns detected, will be used in protocol.dat')
            self.err = True

        # define electrode position according to first survey
        if len(self.surveys) == 1:
            self.elec = None
            self.setElec(self.surveys[0].elec)
        
        #do a check for reading in 3D protocol files for 2D projects 
        if self.typ[-1] == '2' and 'Protocol' in ftype: 
            for s in self.surveys:
                if len(s.df['a'].iloc[0].split()) == 2: 
                    s._rmLineNum() 
                    
        # flag that data has been added 
        self.pinfo['Data'] = True 
        self.pinfo['Number of Surveys'] = 1 
        
        if estMemory: 
            _ = self._estimateMemory(dump=donothing)
        
            
    def addData(self,index=0, **kwargs):
        """Adds data to the survey - used usually to add reciprocal datasets
        
        Parameters
        ----------
        index: int
            Survey index to add data to. 
        **kwargs: Keyword arguments to be passed to Survey.addData()
        """
        
        self.surveys[index].addData(**kwargs)


    def createBatchSurvey(self, dirname, ftype='Syscal', info={}, spacing=None,
                          parser=None, isurveys=[], dump=None, debug=False):
        """Read multiples files from a folders (sorted by alphabetical order).

        Parameters
        ----------
        dirname : str
            Directory with files to be parsed.
        ftype : str, optional
            Type of file to be parsed. Either 'Syscal','ProtocolDC','ResInv',
            'BGS Prime', 'ProtocolIP', 'Sting', 'ABEM-Lund', 'Lippmann' or 'ARES'.
        info : dict, optional
            Dictionnary of info about the survey.
        spacing : float, optional
            Electrode spacing to be passed to the parser function.
        parser : function, optional
            A parser function to be passed to `Survey` constructor.
        isurveys : list, optional
            List of surveys index that will be used for error modelling and so
            reciprocal measurements. By default all surveys are used.
        dump : function, optional
            Function to dump the information message when importing the files.
        debug : bool, optional
            If True informations about reciprocal computation, default filtering
            and so on will be displayed.
        """
        self.createTimeLapseSurvey(dirname=dirname, ftype=ftype, info=info,
                                   spacing=spacing, isurveys=isurveys,
                                   parser=parser, dump=dump, debug=debug)
        self.iTimeLapse = False
        self.iBatch = True
        self.setBorehole(self.iBorehole)


    def createTimeLapseSurvey(self, dirname, ftype='Syscal', info={},
                              spacing=None, parser=None, isurveys=[],
                              dump=None, debug=False):
        """Read electrodes and quadrupoles data and return
        a survey object.

        Parameters
        ----------
        dirname : str or list of str
            Directory with files to be parsed or list of file to be parsed.
        ftype : str, optional
            Type of file to be parsed. Either 'Syscal','ProtocolDC','ResInv',
            'BGS Prime', 'ProtocolIP', 'Sting', 'ABEM-Lund', 'Lippmann' or 'ARES'.
        info : dict, optional
            Dictionnary of info about the survey. Put inverse_type = 1 to allow
            for a changing number of measurements between surveys.
        spacing : float, optional
            Electrode spacing to be passed to the parser function.
        parser : function, optional
            A parser function to be passed to `Survey` constructor.
        isurveys : list, optional
            List of surveys index that will be used for error modelling and so
            reciprocal measurements. By default all surveys are used.
        dump : function, optional
            Function to dump information message when importing the files.
        debug : bool, optional
            If True informations about reciprocal computation, default filtering
            and so on will be displayed.
        """
        if dump is None:
            def dump(x):
                print(x, end="")
        self.iTimeLapse = True
        self.iTimeLapseReciprocal = [] # true if survey has reciprocal
        self.surveys = [] # flush other survey
        if isinstance(dirname, list): # it's a list of filename
            if len(dirname) < 2:
                raise ValueError('at least two files needed for timelapse inversion')
                return
            files = dirname
        else: # it's a directory and we import all the files inside
            if os.path.isdir(dirname):
                files = [os.path.join(dirname, f) for f in np.sort(os.listdir(dirname)) if f[0] != '.']
                # this filter out hidden file as well
            else:
                raise ValueError('dirname should be a directory path or a list of filenames')


        for i, f in enumerate(files):
            self.createSurvey(f, ftype=ftype, parser=parser, spacing=spacing, 
                              debug=debug, estMemory=False)
            haveReciprocal = all(self.surveys[-1].df['irecip'].values == 0)
            self.iTimeLapseReciprocal.append(haveReciprocal)
            dump('\r{:d}/{:d} imported'.format(i+1, len(files)))
            # all surveys are imported whatever their length, they will be matched
            # later if reg_mode == 2 (difference inversion)
        dump('\n')
        self.iTimeLapseReciprocal = np.array(self.iTimeLapseReciprocal)
        # elec and borehole flags assign when first call to R2.createSurvey()
        
        # create bigSurvey (useful if we want to fit a single error model
        # based on the combined data of all the surveys)
        self.bigSurvey = Survey(files[0], ftype=ftype, spacing=spacing, debug=False, parser=parser)
        # then override the df
        if len(isurveys) == 0: # assume all surveys would be use for error modelling
            isurveys = np.ones(len(self.surveys), dtype=bool)
        isurveys = np.where(isurveys)[0] # convert to indices
        df = self.bigSurvey.df.copy()
        c = 0
        for i in isurveys[1:]:
            df2 = self.surveys[i].df
            ipos = df2['irecip'].values > 0
            ineg = df2['irecip'].values < 0
            df2.loc[ipos, 'irecip'] = df2[ipos]['irecip'] + c
            df2.loc[ineg, 'irecip'] = df2[ineg]['irecip'] - c
            #df = df.append(df2, sort=True) # sort to silence the future warning if concatenation axis is not aligned
            df = pd.concat([df, df2], ignore_index=True)
            c = c + df2.shape[0]
        self.bigSurvey.df = df.copy() # override it
        self.bigSurvey.dfOrigin = df.copy()
        
        # flag that data has been added 
        self.pinfo['Data'] = True 
        self.pinfo['Number of Surveys'] = len(self.surveys) 
        
        _ = self._estimateMemory(dump=donothing)


    def create3DSurvey(self, fname, lineSpacing=1, zigzag=False, ftype='Syscal',
                       name=None, parser=None, dump=None):
        """Create a 3D survey based on 2D regularly spaced surveys.
        
        Note: If one electrode has similar position between multiple lines (shared
        electrode), ResIPy will keep online one position of the electrode
        and share its label accross the lines.

        Parameters
        ----------
        fname : list of str
            List of 2D filenames in the right order for the grid or directory
            name (the files will be sorted alphabetically in this last case).
        lineSpacing : float, optional
            Spacing in meter between each line.
        zigzag : bool, optional
            If `True` then one survey out of two will be flipped.
            #TODO not implemented yet
        ftype : str, optional
            Type of file to be parsed. Either 'Syscal','ProtocolDC','ResInv',
            'BGS Prime', 'ProtocolIP', 'Sting', 'ABEM-Lund', 'Lippmann' or 'ARES'.
        name : str, optional
            Name of the merged 3D survey.
        """
        if isinstance(fname, list): # it's a list of filename
            fnames = fname
        else: # it's a directory and we import all the files inside
            if os.path.isdir(fname):
                fnames = [os.path.join(fname, f) for f in np.sort(os.listdir(fname)) if f[0] != '.']
                # this filter out hidden file as well
            else:
                raise ValueError('fname should be a directory path or a list of filenames')
        if dump is None:
            def dump(x):
                print(x)
                
        surveys = []
        for fname in fnames:
            surveys.append(Survey(fname, ftype=ftype, parser=parser,compRecip=False,debug=False))
            dump('Imported datafile: %s'%fname)
        survey0 = surveys[0]
        
        # check this is a regular grid (actually unregular grid works too
        # as we create one mesh for all)
        # nelec = survey0.elec.shape[0]
        # for s in surveys:
        #     if s.elec.shape[0] != nelec:
        #         raise ValueError('Survey {:s} has {:d} electrodes while the first survey has {:d}.'
        #                          'All surveys should have the same number of electrodes.'.format(s.name, s.elec.shape[0], nelec))
        
        # build global electrodes and merged dataframe
        elec = []
        dfs = []
        for i, s in enumerate(surveys):
            e = s.elec.copy()
            e['y'] = i*lineSpacing
            prefix = '{:d} '.format(i+1)
            e['label'] = prefix + e['label']
            elec.append(e)
            df = s.df.copy()
            df[['a','b','m','n']] = prefix + df[['a','b','m','n']]
            dfs.append(df)
        elec = pd.concat(elec, axis=0, sort=False).reset_index(drop=True)
        dfm = pd.concat(dfs, axis=0, sort=False).reset_index(drop=True)
        
        survey0.debug = True 
        survey0.elec = elec
        survey0.df = dfm
        survey0.dfOrigin = dfm # for raw phase plot
        survey0.dfReset = dfm # for reseting filters on res
        survey0.dfPhaseReset = dfm # for reseting filters on IP
        survey0.name = '3Dfrom2Dlines' if name is None else name
        survey0.setSeqIds()
        survey0.computeReciprocal() # compute reciprocals for survey 
        
        self.surveys = [survey0]
        self.elec = None
        self.setElec(elec)
        self.setBorehole(self.iBorehole)
        
        # flag that data has been added 
        self.pinfo['Data'] = True 
        self.pinfo['Number of Surveys'] = 1 
        _ = self._estimateMemory(dump=donothing)
        
    
    def createMergedSurveys(self, fname, ftype='Protocol DC', delimiter=',',
                            dump=None, debug=False):
        """Create one or multiple surveys by merging multiple files referenced 
        in singular human readable file with a specific format, survey details.
        See notes below. 

        Parameters
        ----------
        fname : str
            File path to .csv file. 
        ftype: str, optional
            File type, see ResIPy docs for types of file type avialable. Will be 
            overwritten if specified in the survey details file. 
        delimiter: str, optional
            Delimiter used in merge file. Defaults to ','. 
            
        Notes
        -----
        Format of survey details file should have up to 4 columns with these names. 
        fpath, sid (optional), and ftype (optional), string (optional). 
        
        fpath: str
            Path to survey file, best to use full system path 
        ftype: str, optional 
            File type, should correspond to the file 'ftype' for each file (they
            might be different for some reason, though best avoided). If not
            passed the file type can be set be the variable at the top of this
            function. Type of files avialable are: Either 'Syscal','ProtocolDC',
            'ResInv','BGS Prime', 'ProtocolIP', 'Sting', 'ABEM-Lund',
            'Lippmann' or 'ARES'.
        sid: int, optional 
            Survey index, used for making timelapse surveys from multiple files
        string: int, optional 
            Denotes the string number of the of the sruvey (useful for merging 
            3D surveys from 2D lines). 
        
        """

        finfo = pd.read_csv(fname, sep=delimiter) # file info dataframe 
        
        if 'fpath' not in finfo.columns: 
            msg = 'File paths are not defined in survey merge file, the file headers should read:\n'
            for s in ['fpath','ftype','sid','string']:
                msg += '%s '%s 
            msg += '\nNote that all but "ftype" columns are optional'
            raise Exception(msg)
            
        if 'ftype' not in finfo.columns: 
            finfo['ftype'] = ftype 
        
        if 'sid' not in finfo.columns: 
            finfo['sid'] = 0 
            
        if 'string' not in finfo.columns:
            finfo['string'] = 0 
            
        # get survey ids and unique survey indexes 
        # returned indexes are the first instance of each survey index 
        usid, _ = np.unique(finfo.sid.values, return_index=True)
        
        self.surveys = [] # flush other surveys
        
        # check paths
        paths = []
        for path in finfo['fpath']:
            fpath = path.strip()
            if os.path.isabs(fpath) is False:
                fpath = os.path.join(os.path.dirname(fname), fpath)
            paths.append(fpath)
        finfo['fpath'] = paths

        if len(usid) > 1: 
            self.iTimeLapse = True 
            
        c = 0 
        for i in usid: 
            self.createSurvey(finfo.fpath[i], ftype=finfo.ftype[i], 
                              debug=debug, estMemory=False, string=finfo.string[i])
            sidx = np.argwhere(finfo.sid.values == i).flatten().tolist() # survey index 

            _ = sidx.pop(0)
                
            for j in sidx: 
                self.addData(index=c, fname=finfo.fpath[j], ftype=finfo.ftype[j],
                             string=finfo.string[j])
                # add data for each survey 
            if c == 0:
                self.bigSurvey = Survey(finfo.fpath[i], ftype=finfo.ftype[i],
                                        string=finfo.string[i])
                for j in sidx: 
                    self.bigSurvey.addData(fname=finfo.fpath[j], ftype=finfo.ftype[j],
                                           string=finfo.string[j])
            c += 1 

        elecids = []
        for s in self.surveys:
            elecids += s.df[['a','b','m','n']].values.flatten().tolist()
            
        uelecids = sorted(np.unique(elecids)) # is this efficient enough for strings ?
        
        if self.elec is None or len(uelecids) != len(self.elec):
            xplaceholder = np.arange(len(uelecids),dtype=float)
            yplaceholder = np.zeros_like(xplaceholder)
            zplaceholder = np.zeros_like(xplaceholder)
            placeholder = [False]*len(uelecids)
            elec = {
                'label':uelecids,
                'x':xplaceholder, 
                'y':yplaceholder,
                'z':zplaceholder,
                'buried':placeholder,
                'remote':placeholder,
                }
            elec = pd.DataFrame(elec)
            self.elec = elec 
            for s in self.surveys:
                s.elec = elec.copy()
                
        self.mergeElec() 
        
        if self.iTimeLapse: # bit of clean up regarding timelapse case 
            self.iTimeLapseReciprocal = np.array([False]*len(self.surveys)) # true if survey has reciprocal
            for i in range(len(self.surveys)):
                haveReciprocal = all(self.surveys[i].df['irecip'].values == 0)
                self.iTimeLapseReciprocal[i] = haveReciprocal 
            
            # create bigSurvey (useful if we want to fit a single error model
            # based on the combined data of all the surveys)
            df = self.bigSurvey.df.copy()
            for i in range(1,len(self.surveys)):
                df2 = self.surveys[i].df
                ipos = df2['irecip'].values > 0
                ineg = df2['irecip'].values < 0
                df2.loc[ipos, 'irecip'] = df2[ipos]['irecip'] + c
                df2.loc[ineg, 'irecip'] = df2[ineg]['irecip'] - c
                df = pd.concat([df, df2], ignore_index=True)
                c = c + df2.shape[0]
            self.bigSurvey.df = df.copy() # override it
            self.bigSurvey.dfOrigin = df.copy()
            self.bigSurvey.ndata = df.shape[0]
        else:
            self.bigSurvey = None 
        
        # flag that data has been added 
        self.pinfo['Data'] = True 
        self.pinfo['Number of Surveys'] = len(self.surveys) 
        
        _ = self._estimateMemory(dump=donothing)


    def createPseudo3DSurvey(self, dirname, lineSpacing=1, ftype='Syscal', parser=None, **kwargs):
        """Create a pseudo 3D survey based on 2D surveys. Multiple 2D Projects to be turned into a single pseudo 3D survey.
            THIS WILL NEED CORRECT ELECTRODE LAYOUT - DONE IN self._updatePseudo3DSurvey()
        
        Parameters
        ----------
        dirname : list of str
            List of 2D filenames in the right order for the grid or directory
            name (the files will be sorted alphabetically in this last case).
        lineSpacing : float, optional
            Spacing in meter between each line.
        ftype : str, optional
            Type of the survey to choose which parser to use.
        kwargs : -
            Keyword arguments to be passed to Project.createBatchSurvey()
        """
        
        self.createBatchSurvey(dirname=dirname, ftype=ftype, parser=parser, **kwargs) # We need surveys in the master Project for data procesing (error modeling, etc.)
    
        elecList = []
        dfList = []
        for i, s in enumerate(self.surveys):
            directory = os.path.join(self.dirname, s.name)
            os.mkdir(directory) # making separate inversion diectories
            proj = self._createProjects4Pseudo3D(dirname=directory, invtyp=self.typ) # non-parallel meshing
            # proj.createSurvey(fname=None, name=s.name, df=s.df, elec=s.elec, **kwargs)
            proj.surveys = [s] # might be faster than importing survey again
            # this actually stores a reference to the same survey object, so it's all good
            # as filtering operation done in batch surveys will be directly done on the surveys of 
            # each project as it's the same Survey object (with two references)
            self.projs.append(proj) # appending projects list for later use of meshing and inversion
            e = s.elec.copy()
            self.projs[-1].elec = e 
            e['y'] = i*lineSpacing
            prefix = '{:d} '.format(i+1)
            e['label'] = prefix + e['label']
            elecList.append(e)
            df = s.df.copy()
            df[['a','b','m','n']] = prefix + df[['a','b','m','n']]
            dfList.append(df)
            
        elec = pd.concat(elecList, axis=0, sort=False).reset_index(drop=True)
        dfm = pd.concat(dfList, axis=0, sort=False).reset_index(drop=True)
        
        self.pseudo3DSurvey = Survey(fname=None, df=dfm, elec=elec)
        self.elec = None
        self.setElec(elec) # create initial electrodes df - to be populated later
        self.setBorehole(self.iBorehole)
        # flag that data has been added 
        self.pinfo['Data'] = True 
        self.pinfo['Number of Surveys'] = len(self.surveys)
        _ = self._estimateMemory(dump=donothing)
        
    
    
    def importPseudo3DElec(self, fname=''):
        """Import electrodes positions. The label columns should include line
        number separated by space (like in 3D):
            label,x,y,z
            1 3,0,0,0
            1 4,1,1,0
            1 5,1,2,1

        Parameters
        ----------
        fname : str
            Path of the CSV file containing the electrodes positions. It should contains 3 columns maximum with the X, Y, Z positions of the electrodes.
        """
        with open(fname, 'r') as f:
            try:
                float(f.readline().split(',')[0])
                header = None
            except Exception:
                header = 'infer'
        df = pd.read_csv(fname, header=header)
        
        if header is None:
            x = df.values[:,0]
            y = df.values[:,0]
            z = df.values[:,0]
            elec = pd.DataFrame({'x':x, 'y':y, 'z':z})
        else:
            elec = df
        self.setPseudo3DElec(elec)
    
    
    
    def setPseudo3DElec(self, elec):
        """Set pseudo 3D electrodes (with an electrode label as:
            <line number> <electrode number>).
    
        Parameters
        ----------
        elec : pd.DataFrame 
            Electrode data frame 
        """
        if 'line' in elec.columns and 'number' in elec.columns: 
            label = ['%i %i'%(elec.line[i], elec.number[i]) for i in range(len(elec))]
            elec['label'] = label 
        if 'remote' not in elec.columns: 
            elec['remote'] = False 
        if 'buried' not in elec.columns: 
            elec['buried'] = False 
        
        self.pseudo3DSurvey.elec = elec
        
        # take self.surveys information to inform all projects in self.projs
        self._updatePseudo3DSurvey()
        
        
        
    def _updatePseudo3DSurvey(self, elecList=None):
        """Update a pseudo 3D survey based on 2D surveys. 
            Cleaned data, updated electrodes will be inserted in each survey.
        
        Parameters
        ----------
        elecList : list of dataframes, optional
            List of electrodes dataframes - each df must have 2D like XYZ (rotated to have y=0).
        """
        if self.projs == []:
            raise ValueError('Survey needs to be created first! use Project.createPseudo3DSurvey()')
            
        elecList = self._create2DLines(elecList)

        for elecdf, proj, survey in zip(elecList, self.projs, self.surveys):
            survey.elec = elecdf.copy()
            proj.setElec(elecdf.copy())
            proj.surveys[0].df = survey.df.copy()
            proj.typ = self.typ    
            proj.err = self.err
            if self.pseudo3Dfmd is not None: # pseudo3D fmd has been forced elsewhere so needs to be added here
                proj.fmd = self.pseudo3Dfmd


    def split3DGrid(self, elec=None, changeLabel=True):
        """Split self.elec to available lines based on 'label' 
        
        
        Parameters
        ----------
        elec : dataframe, optional
            Contains the electrodes information. "label" column must be provided and
            have "<line number> <electrode number>" format.
        changeLable : bool, optional
            If True, the line number will be dropped from labels - Flase for GUI related surveys.
        
        Returns
        -------
        elecList : list of dataframes
            List of electrodes dataframes - each df can have a 3D like XYZ.
        """
        if elec is None:
            if self.pseudo3DSurvey is not None:
               elec =  self.pseudo3DSurvey.elec.copy()
            else:
               elec = self.elec.copy()
        elec[['lineNum', 'elecNum']] = elec['label'].str.split(expand=True)
        elecGroups = elec.groupby('lineNum', sort=False)#, as_index=False)
        elecdfs = [elecGroups.get_group(x) for x in elecGroups.groups]
        elecList = []
        for elecdfRaw in elecdfs:
            elecdf = elecdfRaw.copy().reset_index(drop=True) # void pandas setting with copy warning annoying error
            if 'line' in elecdf:
                elecdf.drop(columns=['line'], inplace=True)
            if 'number' in elecdf:
                elecdf.drop(columns=['number'], inplace=True)
            if changeLabel:
                elecdf['label'] = elecdf['elecNum'].values # it's 2D so let's get rid of line numbers in labels
            elecdf = elecdf.drop(['lineNum', 'elecNum'], axis=1)
            elecdf = self._findRemote(elecdf)
            elecList.append(elecdf)
        return elecList
    
    
    def _create2DLines(self, elecList=None):
        """Create a list of 2D electrode XYZ/topo where only x & z are variable and y=0.
            Simply, rotating an array of XY locations on x, y = 0 pivot to have all y values equal to zero
        
        Parameters
        ----------
        elecList : list of dataframes, optional
            List of electrodes dataframes - each df can have a 3D like XYZ.
        
        Returns
        -------
        elecList : list of dataframes
            List of electrodes dataframes - each df will have 2D like XYZ (rotated to have y=0).
        """
        
        if elecList is None:
            elecList = self.split3DGrid()
            
        for elecdf in elecList:
            # transforming line to start from x, y = 0
            ie = elecdf[~elecdf['remote']].index.values
            elecdf.loc[ie,'x'] = elecdf.loc[ie,'x'].values - elecdf.loc[np.argmin(elecdf.loc[ie,'x']),'x']
            if np.all(elecdf.loc[ie,'x'] == elecdf.loc[0,'x']): # line is vertical
                elecdf.loc[ie,'y'] = elecdf.loc[ie,'y'].values - elecdf.loc[np.argmin(elecdf.loc[ie,'y']),'y']
            else:
                elecdf.loc[ie,'y'] = elecdf.loc[ie,'y'].values - elecdf.loc[np.argmin(elecdf.loc[ie,'x']),'y']
            delx = elecdf.loc[ie,'x'].max() - elecdf.loc[ie,'x'].min()
            dely = elecdf.loc[ie,'y'].max() - elecdf.loc[ie,'y'].min()
            f = np.inf if delx == 0 else np.abs((dely)/(delx))
            rotangle = np.arctan(f)
            if elecdf.loc[np.argmin(elecdf.loc[ie,'x']),'y'] > elecdf.loc[np.argmax(elecdf.loc[ie,'x']),'y']: # CCW rotation needed
                rotangle *= -1
            # rotation     
            rotmat = np.array([[np.cos(rotangle), np.sin(rotangle)], 
                                [-np.sin(rotangle), np.cos(rotangle)]])
            xy = np.array([elecdf.loc[ie,'x'].values, elecdf.loc[ie,'y'].values])
            newmat = np.dot(rotmat, xy).T
            elecdf.loc[ie,'x'] = newmat[:,0].copy()
            elecdf['y'] = 0 # to make sure we don't end up with super small values
                        
        return elecList
    
    
    def createMultiMesh(self, runParallel=True, **kwargs):
        """Create multiple meshes from avalable Projects in self.projs.
        
        Parameters
        ----------
        runParallel : bool, optional
            if True, mesh generation will run in multiple threads.
        kwargs : -
            Keyword arguments to be passed to mesh generation schemes
        """
        # TODO the below doesn't work well, we need to use
        # 1) tempfile for creating temporary directory and then the 
        # write all files to those and then use a similar approach to 
        # Project.runParallel() with a while loop and limiting given ncores
        # however, it's fast enough like this as it's a sequence of 2D meshes
        for proj in self.projs:
            if runParallel:
                p = Thread(target=proj.createMesh, kwargs=kwargs)
                p.start()
                p.join()
            else:
                proj.createMesh(**kwargs)
                
        self.mesh = self.projs[0].mesh # just to have a populated mesh in master Project!
        if 'fmd' in kwargs:
            self.pseudo3Dfmd = kwargs['fmd']
            
        # check ram requirements 
        _ = self._estimateMemory(dump=donothing) 
   
    
    
    def showPseudo3DMesh(self, ax=None, color_map='Greys', meshList=None,
                         cropMesh=True, color_bar=False, returnMesh=False,
                         cropMaxDepth=False, clipCorners=False, pvshow=True, **kwargs):
        """Show 2D meshes in 3D view
        
        Parameters
        ----------
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.
        color_map : str, optional
            Name of the colormap to be used.
        meshList : list of Mesh classes
            If not None, pseudo 3D meshes will be plotted by default.
        cropMesh : bool, optional
            If True, 2D mesh will be bound to electrodes and zlim.
        color_bar : Boolean, optional 
            `True` to plot colorbar.
        returnMesh: bool, optional
            if True method returns a merged mesh. 
        cropMaxDepth : bool, optional
            If True, region below fine mesh depth (fmd) will be cropped.
        clipCorners : bool, optional
            If 'True', triangles from bottom corners will be cropped (only if the whole mesh is not shown).
        pvshow : bool, optional
            Set to False to not call `Plotter.show()` (useful for subplots).
        kwargs : -
            Keyword arguments to be passed to Mesh.show() class.
        """
        
        try:
            import pyvista as pv
        except Exception:
            print('ERROR: pyvista is needed to show pseudo 3D meshes. Use pip install pyvista')
            return
        
        if meshList is None:
            meshList = [deepcopy(p.mesh) for p in self.projs]
        
        if ax is None:
            ax = pv.Plotter()
            
        kwargs['pvshow'] = False # don't invoke show after each mesh added
        
        elecList = self.split3DGrid()  # split the electrodes to lines in 3D space   
                
        if returnMesh:
            meshOutList = []
            
        for proj, elecdf, mesh in zip(self.projs, elecList, meshList):
            if mesh is None:
                print('Mesh undefined for this project!')
                continue

            if cropMesh:
                node_x = mesh.node[:,0]
                node_z = mesh.node[:,2]
                xmin = np.min(node_x)
                xmax = np.max(node_x)
                zmin = np.min(node_z)
                zmax = np.max(node_z)
                (xsurf, zsurf) = mesh.extractSurface()
                if cropMaxDepth and proj.fmd is not None:
                    xfmd, zfmd = xsurf[::-1], zsurf[::-1] - proj.fmd
                    verts = np.c_[np.r_[xmin, xmin, xsurf, xmax, xmax, xfmd, xmin],
                                  np.r_[zmin, zmax, zsurf, zmax, zmin, zfmd, zmin]]
                else:
                    verts = np.c_[np.r_[xmin, xmin, xsurf, xmax, xmax, xmin],
                                  np.r_[zmin, zmax, zsurf, zmax, zmin, zmin]]
                mesh = mesh.crop(verts)
            
            if clipCorners and proj.param['num_xz_poly'] != 0: # not clipping the corners of a mesh outside of the survey area!
                elec_x = mesh.elec[:,0]
                elec_z = mesh.elec[:,2]
                elec_xmin = np.min(elec_x)
                elec_xmax = np.max(elec_x)
                if cropMaxDepth and proj.fmd is not None:
                    zminl = elec_z[elec_x.argmin()] - proj.fmd
                    zminr = elec_z[elec_x.argmax()] - proj.fmd
                elif cropMaxDepth is False and proj.fmd is not None:
                    zminl = zminr = np.min(elec_z) - proj.fmd
                else:
                    zminl = zminr = zmin
                zmaxl = elec_z[elec_x.argmin()]
                zmaxr = elec_z[elec_x.argmax()]
                ll = np.abs(zmaxl - zminl)
                lr = np.abs(zmaxr - zminr)
                lx = np.abs(elec_xmin - elec_xmax)
                if ll >= lx/4:
                    ll = lx/4
                if lr >= lx/4:
                    lr = lx/4
    
                # surf bound to elecs
                elec_surf = np.c_[xsurf, zsurf][(np.abs(xsurf - elec_xmin)).argmin():
                                                (np.abs(xsurf - elec_xmax)).argmin(),:]
                idxl = (np.abs(elec_surf[:,0] - ll)).argmin()
                idxr = (np.abs(elec_surf[:,0] - np.abs(elec_xmax - lr))).argmin() + 1       
                xtrapbot = elec_surf[idxl:idxr,0]
                if cropMaxDepth and proj.fmd is not None:
                    ztrapbot = elec_surf[idxl:idxr,1] - proj.fmd
                else:
                    ztrapbot = np.ones_like(xtrapbot) * zminl
                    
                proj.trapeziod = np.c_[np.r_[elec_xmin, xsurf, elec_xmax, xtrapbot[::-1], elec_xmin],
                                  np.r_[zmaxl, zsurf, zmaxr, ztrapbot[::-1], zmaxl]]
                mesh = mesh.crop(proj.trapeziod)
            
            else:
                proj.trapeziod = None # make sure trapeziod mask is clear

            meshMoved = mt.moveMesh2D(meshObject=mesh, elecLocal=proj.elec, elecGrid=elecdf)
            
            if cropMesh:
                limits = np.c_[meshMoved.elec[:,0][~meshMoved.iremote], meshMoved.elec[:,1][~meshMoved.iremote]]
            else:
                limits = np.c_[meshMoved.node[:,0], meshMoved.node[:,1]]

            xlim = findminmax(limits[:,0],1)
            ylim = findminmax(limits[:,1],1)
            zlim = proj.zlim if cropMesh or cropMaxDepth else None
            meshMoved.ndims = 3 # overwrite dimension to use show3D() method
            
            # best to let pseudo 3d algorithm choose extents? 
            kwargs['xlim'] = xlim 
            kwargs['ylim'] = ylim 
            kwargs['zlim'] = zlim 

            meshMoved.show(ax=ax, color_map=color_map, color_bar=color_bar,
                           darkMode=self.darkMode, **kwargs)
            if returnMesh:
                meshOutList.append(meshMoved)
        
        if pvshow:
            ax.show() # call plotter.show()
        
        if returnMesh:
            self.pseudo3DMeshResultList = meshOutList
            meshMerged = mt.mergeMeshes(meshOutList)
            meshMerged.ndims = 3
            self.pseudo3DMeshResult = meshMerged
    
    
    def _setPseudo3DParam(self, targetProjParams):
        """Set params from master Project to a target Project.
            IMPORTANT: some params (e.g., mesh, reg0, etc.) are excluded
        
        Parameters
        ----------
        targetProjParams : dict
            Target Project instance params
        
        Returns
        -------
        targetProjParams : dict
            Target Project instance params are set from self (master Project)
        """
        keys = ['num_xz_poly', 'num_xy_poly','a_wgt', 'b_wgt', 'lineTitle', 'job_type',  
                'flux_type', 'singular_type', 'res_matrix', 'scale', 'patch_x', 'patch_z',
                'inverse_type', 'target_decrease', 'qual_ratio', 'data_type', 'zmin', 'zmax',
                'reg_mode', 'tolerance', 'max_iter', 'error_mod', 'alpha_aniso',
                'alpha_s', 'min_error', 'rho_min', 'rho_max', 'mesh_type']
        
        for key in keys:
            if key in self.param.keys():
                targetProjParams[key] = self.param[key]
        if 'node_elec' in self.param:
            targetProjParams['node_elec'] = [self.param['node_elec'][0].tolist(),
                                    self.param['node_elec'][1].tolist()]
        if 'xz_poly_table' in self.param:
            targetProjParams['xz_poly_table'] = self.param['xz_poly_table'].tolist()
        if 'xy_poly_table' in self.param:
            targetProjParams['xy_poly_table'] = self.param['xy_poly_table'].tolist()
        
        return targetProjParams
    
    def _checkMemory(self,parallel=False, ncores=None):
        """
        Check the memory allocation of the system that could be used for inversion. 
        Adjusts the number of parallel threads to match amount of RAM available. 

        Parameters
        ----------
        parallel : bool, optional
            If True, ResIPy can run in parallel in time-lapse or batch modes. 
            The default is False.
        ncores : int, optional
            Number of parallel threads allowed. The default is None.

        Raises
        ------
        Exception
            If memory estimate required for inversion exceeds the total system memory 
            then an error is thrown! 

        Returns
        -------
        ncores: int
            Number of parallel threads adjusted for system memory 

        """
        
        if ncores is None: 
            ncores = sysinfo['physicalCpuCount']
            
        # get memory estimates 
        memInv = self._estimateMemory(dump=donothing)
        memTot = sysinfo['totalMemory'] #use the total memory (its more conservative)
        
        # check if static inversion possible 
        if memInv > memTot: 
            raise Exception('The amount of memory required for inversion exceeds that of known system resources (please make mesh coarser)')
            return 0 
         
        # if not parallel we can exit here 
        if not parallel:
            return 1 
        
        maxInversions = 1 
        adjustTrigger = False 
        
        if len(self.surveys) > 1:
            nsurveys = len(self.surveys) - 1
            # number of surveys that can possibly be inverted in parallel is 
            # n - 1 becuase of the baseline which is inverted independently 
            maxInversions = min([nsurveys,ncores]) # maximum number of inversions at one time 
                            
            # provide the user with some safety net 
            while (maxInversions*memInv) > memTot: 
                adjustTrigger = True 
                ncores -= 1 
                maxInversions = ncores # maximum number of inversions at one time 
                if ncores == 1:
                    break 
                    
        # flag to the user they're using too many cores 
        if adjustTrigger: 
            warnings.warn('Adjusting parallel threads (to %i) to match RAM available'%ncores)
                    
        # check if static inversion possible 
        if (memInv*maxInversions) > memTot: 
            raise Exception('The amount of memory required for time-lapse inversion exceeds that of known system resources')
                    
        return ncores 
            
    
    def invertPseudo3D(self, invLog=None, runParallel=False, **kwargs):
        """Run pseudo3D inversions.
        
        Parameters
        ----------
        invLog : function, optional
            Passes project inversion outputs.
        runParallel : bool
            if True, inversions will run in parallel based on number of CPU cores.
        kwargs : -
            Keyword arguments to be passed to invert().
        """
        # check memory, retrieve optimum number of cores 
        ncores = self._checkMemory(runParallel,None)
        
        # also, as each proj will be treated independently, we dont want the user 
        # calling the parallel functionalities of Project.invert in the pseudo3D case 
        if 'parallel' in kwargs.keys():
            kwargs['parallel'] = False 
        if 'ncores' in kwargs.keys():
            kwargs['ncores'] = None  
        
        # kill management
        self.procs = []
        self.proc = ProcsManagement(self)
        self._updatePseudo3DSurvey() # make sure we have set all attributes
        self.meshResults = [] # clean meshResults list
        for proj in self.projs: # preparing inversion params
            proj.param = self._setPseudo3DParam(proj.param)
        
        if runParallel is False: # non-parallel inversion
            for proj in self.projs:
                self.projectPseudo3D = proj # get functions for UI
                proj.invert(**kwargs)
                self.procs.append(proj.proc)
                if self.proc.killFlag is True:
                    break
                if invLog is not None:
                    invLog(proj.dirname, proj.typ)
            
        else: # parallel inversion
            if invLog is None:
                def invLog(x):
                    print(x, end='')
            
            # create R2.exe path
            exeName = self.typ + '.exe'
            exePath = os.path.join(self.apiPath, 'exe', exeName)
    
            # create .in and protocol.dat files
            wds = []
            for proj in self.projs:
                wds.append(proj.dirname)
                invLog('Writing .in file and protocol.dat for {} survey... '.format(
                    proj.surveys[0].name))
                proj.write2in() # R2.in
                proj.write2protocol() # protocol.dat
                invLog('done\n')
            wds2 = wds.copy()
    
            # create workers directory
            if ncores > sysinfo['cpuCount']:
                raise ValueError('Number of cores larger than available')
    
            if OS == 'Windows':
                cmd = [exePath]
            elif OS == 'Darwin':
                # winetxt = 'wine'
                # if getMacOSVersion():
                #     winetxt = 'wine64'
                winetxt = whichWineMac()
                winePath = []
                wine_path = Popen(['which', winetxt], stdout=PIPE, shell=False, universal_newlines=True)
                for stdout_line in iter(wine_path.stdout.readline, ''):
                    winePath.append(stdout_line)
                if winePath != []:
                    cmd = ['%s' % (winePath[0].strip('\n')), exePath]
                else:
                    cmd = [wPath, exePath]
            else:
                #if platform.machine() == 'aarch64':
                #    cmd = [exePath.replace('.exe', '_aarch64')]
                if platform.machine() in ['armv7l', 'aarch64']:
                    cmd = [exePath.replace('.exe', '_armv7l')]
                else:
                    cmd = ['wine',exePath]
    
            if OS == 'Windows':
                startupinfo = subprocess.STARTUPINFO()
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    
            # run them all in parallel as child processes
            invLog('----------- PARALLEL INVERSION BEGINS ----------\n')
            def dumpOutput(out):
                for line in iter(out.readline, ''):
                    invLog(line.rstrip() + '\n')
                out.close()
    
            # create essential attribute
            self.irunParallel2 = True
    
            # run in // (http://code.activestate.com/recipes/577376-simple-way-to-execute-multiple-process-in-parallel/)
            # In an infinite loop, will run an number of process (according to the number of cores)
            # the loop will check when they finish and start new ones.
            def done(p):
                return p.poll() is not None

            c = 0
            invLog('\r{:.0f}/{:.0f} inversions completed'.format(c, len(wds2)))
            while self.irunParallel2:
                while wds and len(self.procs) < ncores:
                    wd = wds.pop()
                    if OS == 'Windows':
                        p = Popen(cmd, cwd=wd, stdout=PIPE, shell=False, universal_newlines=True, startupinfo=startupinfo)
                    else:
                        p = Popen(cmd, cwd=wd, stdout=PIPE, shell=False, universal_newlines=True)
                    self.procs.append(p)
    
                for p in self.procs:
                    if done(p):
                        self.procs.remove(p)
                        c = c+1
                        # TODO get RMS and iteration number here ?
                        invLog('\r{:.0f}/{:.0f} inversions completed'.format(c, len(wds2)))
    
                if not self.procs and not wds:
                    invLog('\n')
                    break
                else:
                    time.sleep(0.05)

        self.invLog = '' # clearing the inversion log for saving
        if self.proc.killFlag is False: # make sure we haven't killed the processes
            for proj, survey in zip(self.projs, self.surveys):
                if runParallel:
                    proj.getInvError()
                    proj.getResults()
                self.meshResults.append(deepcopy(proj.meshResults[0]))
                survey.df = proj.surveys[0].df.copy() # to populate inversion error outputs
                survey.dfInvErrOutputOrigin = survey.df.copy()
                survey.isequenceInvErrOutputReset = survey.isequence.copy()
                if proj.invLog == '': # in case of parallel inversion
                    with open(os.path.join(proj.dirname, proj.typ + '.out'),'r') as f:
                        proj.invLog = f.read()
                self.invLog += '###>>> Dataset: ' + survey.name + proj.invLog + '\n'

        print('----------- END OF INVERSION IN // ----------')

    
    def showPseudo3DResults(self, cropMesh=False, **kwargs):
        """Show 2D Inversions in 3D view
        
        Parameters
        ----------
        cropMesh : bool, optional
            If True, 2D mesh will be bound to electrodes and zlim.
        kwargs : -
            Keyword arguments to be passed to showPseudo3DMesh().
        """
        meshResults = [deepcopy(p.meshResults[0]) for p in self.projs]
        self.showPseudo3DMesh(color_bar=True, cropMesh=cropMesh, meshList=meshResults, **kwargs)
        
    
    @classmethod
    def _createProjects4Pseudo3D(cls, dirname, invtyp='R2'):
        """Create a Project instance for future use of meshing and inversion.
        
        Parameters
        ----------
        dirname : str
            directory where 2D line will be dealt with (meshing, inversion, etc.)
        invtyp : str
            'R2' - inverting 2D resistivity
            'cR2' - inverting 2D induced polarization
        
        Returns
        -------
        ProjInstance : Instance of Project class
        """
        ProjInstance = cls(dirname=dirname, typ=invtyp)
        ProjInstance.dirname = dirname
        shutil.rmtree(os.path.join(dirname, 'invdir')) # we don't want this invdir anymore
        return ProjInstance


    def showPseudo(self, index=0, vmin=None, vmax=None, ax=None, **kwargs):
        """Plot pseudo-section with dots.

        Parameters
        ----------
        index : int, optional
            Index of the survey to be used for the pseudo-section (in case of
            timelapse or batch).
        vmin : float, optional
            Minimum value for colorscale.
        vmax : float, optional
            Maximum value for colorscale.
        ax : matplotlib.Axes, optional
            If specified, axis along which to plot the graph.
        **kwargs : optional
            Passed to `Survey.showPseudo()`.
        """
        self.surveys[index].showPseudo(vmin=vmin, vmax=vmax, ax=ax, **kwargs)


    def showPseudoIP(self, index=0, vmin=None, vmax=None, ax=None, **kwargs):
        """Plot pseudo-section with dots for IP data.

        Parameters
        ----------
        index : int, optional
            Index of the survey to be used for the pseudo-section (in case of
            timelapse or batch).
        vmin : float, optional
            Minimum value for colorscale.
        vmax : float, optional
            Maximum value for colorscale.
        ax : matplotlib.Axes, optional
            If specified, axis along which to plot the graph.
        **kwargs : optional
            Passed to `Survey.showPseudoIP()`.
        """
        self.surveys[index].showPseudoIP(vmin=vmin, vmax=vmax, ax=ax, **kwargs)


    def matchSurveys(self):
        """Will trim all surveys to get them ready for difference inversion
        so that each PAIRS of (background, surveyX) have the same number of
        quadrupoles. We do not take all quadrupoles in common among all
        surveys as this is not needed and if there is a small survey, it would
        reduce all other larger surveys.
        """
        print('Matching quadrupoles between pairs of (background, surveyX) for difference inversion...', end='')
        t0 = time.time()
            
        df0 = self.surveys[0].df.reset_index(drop=True)
        df0['tlindex0'] = df0.index.astype(int)
        df0['irecip0'] = df0.irecip.copy()
        ie0 = np.ones(df0.shape[0], dtype=bool)
        indexes = [(ie0, ie0)]  # array of tuple
        for survey in self.surveys[1:]:
            df = survey.df.reset_index(drop=True)
            df['tlindex'] = df.index.astype(int)
            dfm = pd.merge(df0[['a', 'b', 'm', 'n', 'tlindex0','irecip0']],
                           df[['a', 'b', 'm', 'n', 'tlindex','irecip']],
                           how='inner', on=['a', 'b', 'm', 'n'])
            # don't keep reciprocals
            dfm = dfm[dfm[['irecip0', 'irecip']].ge(0).all(axis=1)]
            ie0 = np.zeros(df0.shape[0], dtype=bool)
            ie0[dfm['tlindex0'].values] = True
            ie = np.zeros(df.shape[0], dtype=bool)
            ie[dfm['tlindex'].values] = True
            # dont keep reciprocals too 
#            for i in dfm['tlindex0'].values:
#                if df0.irecip0[i] < 0: 
#                    ie0[i] = False 
#            for i in dfm['tlindex'].values:
#                if df.irecip[i] < 0: 
#                    ie[i] = False  
            indexes.append((ie0, ie))

        print('done in {:.3}s'.format(time.time()-t0))

        return indexes


    def showError(self, index=0, ax=None):
        """Plot the reciprocal errors.

        Parameters
        ----------
        index : int, optional
            Index of the survey to plot. If `index == -1` then all combined
            data of all survey will be plotted together. Default is to plot
            the first survey (`index==0`).
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.

        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if index == -2: # show them all from bigSurvey
            self.bigSurvey.showError(ax=ax)
        else:
            self.surveys[index].showError(ax=ax)


    def showErrorDist(self, index=0, ax=None):
        """Calculate and plots reciprocal error probablity histogram.
        Good data will have a bell shape (normal) distribution where most datapoints have near
        zero reciprocal error.

        Parameters
        ----------
        index : int, optional
            Index of the survey to plot. Default is first survey `index == 0`.
            If `index == -2` then the error distribution of the combined data
            will be plotted.
        ax : Matplotlib.Axes
            If specified, the graph will be plotted against it.
        """
        if index == -2: # show them all from bigSurvey
            self.bigSurvey.showErrorDist(ax=ax)
        else:
            self.surveys[index].showErrorDist(ax=ax)


    def filterManual(self, index=-1, ax=None, **kwargs):
        """Interactive manually filters the data visually. The manually selected
        points index are stored in `Survey.iselect` or `Survey.eselect``if it is
        an electrodes. Use `Survey.filterData()` to filter them out for a single
        survey. Or `R2._filterSimilarQuads()` to filter quadrupoles amongs all
        `R2.surveys`.
        """
        if self.typ[-1] == 't':
            kwargs['flag3d'] = True 
        if index == -1:
            for s in self.surveys:
                s.iselect = np.zeros(s.df.shape[0], dtype=bool)
                s.eselect = np.zeros(len(s.elec), dtype=bool)
            self.surveys[0].filterManual(ax=ax, darkMode=self.darkMode, **kwargs)
        else:
            self.surveys[index].filterManual(ax=ax, darkMode=self.darkMode, **kwargs)


    def filterDummy(self, index=-1):
        """Remove measurements where abs(a-b) != abs(m-n) (likely to be dummy
        measurements added for speed).

        Parameters
        ----------
        index : int, optional
            Index of the survey to process. If `index == -1` (default) then the
            processing is applied on all survey independantly.
        """
        if index == -1:
            for s in self.surveys:
                s.filterDummy()
        else:
            self.surveys[index].filterDummy()
            
    def fixLegendItems(self,ax):
        """Display legend items with survey names in the case of fitting 
        individual error models to multiple data sets. 

        Parameters
        ----------
        ax : matplotlib axes 
            Axes handle with reciprocal error plots in it. 

        """
        newlegend = [] # list to store new legend entries 
        legends = ax.get_legend().get_texts() # return text objects of the legend entries 
        c = 0 # legend count number 
        i = 0 # survey iteration number 
        maxleg = len(self.surveys)*3 # maximum number of legends that should be possible 
        for a in legends: 
            if c%3 == 0 and c!=0:
                i+=1 
            name = self.surveys[i].name 
            newlegend.append('{:s} {:s}'.format(name, a.get_text()))
            c+=1 # add one every time loop is run 
            if c == maxleg:
                break # TODO: maybe delete surplus legend items in future? 
            
        ax.get_legend().remove()
        ax.legend(newlegend, fontsize=10) # replace legend entries with survey names appended 
        ax.set_title(ax.get_title().split('\n')[0])


    def fitErrorLin(self, index=-1, ax=None):
        """Fit a linear relationship to the resistivity data.

        Parameters
        ----------
        index : int, optional
            Index of the survey to fit. If `index == -1` (default) then the fit
            is done on all surveys independantly.
            If `ìndex == -2` then the fit is done on the combined surveys.
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.

        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if index == -2: # apply to combined survey
            self.bigSurvey.fitErrorLin(ax=ax)
            for s in self.surveys:
                s.df['resError'] = self.bigSurvey.errorModel(s.df)
        elif index == -1: # apply to each
            for s in self.surveys:
                s.fitErrorLin(ax=ax)
            # redo the legend
            if ax is not None and len(self.surveys) > 1:
                self.fixLegendItems(ax)
        else:
            self.surveys[index].fitErrorLin(ax=ax)
        self.pinfo['Error model'] = True 


    def fitErrorPwl(self, index=-1, ax=None):
        """Fit an power law to the resistivity data.

        Parameters
        ----------
        index : int, optional
            Index of the survey to fit. If `index == -1` (default) then the fit
            is done on all surveys independantly.
            If `ìndex == -2` then the fit is done on the combined surveys.
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.

        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if index == -2: # apply to combined data of bigSurvey
            self.bigSurvey.fitErrorPwl(ax=ax)
            for s in self.surveys:
                s.df['resError'] = self.bigSurvey.errorModel(s.df)
        elif index == -1: # apply to each
            for s in self.surveys:
                s.fitErrorPwl(ax=ax)
            # redo the legend
            if ax is not None and len(self.surveys) > 1:
                self.fixLegendItems(ax)
        else:
            self.surveys[index].fitErrorPwl(ax=ax)
        self.pinfo['Error model'] = True 

    def fitErrorPwlEnv(self, index=-1, ax=None):
        """Fit an power law to the resistivity data with an envelope fit.

        Parameters
        ----------
        index : int, optional
            Index of the survey to fit. If `index == -1` (default) then the fit
            is done on all surveys independantly.
            If `ìndex == -2` then the fit is done on the combined surveys.
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.

        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if index == -2: # apply to combined data of bigSurvey
            self.bigSurvey.fitErrorPwlEnv(ax=ax)
            for s in self.surveys:
                s.df['resError'] = self.bigSurvey.errorModel(s.df)
        elif index == -1: # apply to each
            for s in self.surveys:
                s.fitErrorPwlEnv(ax=ax)
            # redo the legend
            if ax is not None and len(self.surveys) > 1:
                self.fixLegendItems(ax)
        else:
            self.surveys[index].fitErrorPwlEnv(ax=ax)
        self.pinfo['Error model'] = True 
        
    def fitErrorLME(self, index=-1, ax=None, rpath=None, iplot=True): # :pragma: no cover
        """Fit a linear mixed effect (LME) model by having the electrodes as
        as grouping variables.

        Parameters
        ----------
        Index of the survey to fit. If `index == -1` (default) then the fit
            is done on all surveys independantly.
            If `ìndex == -2` then the fit is done on the combined surveys.
        ax : matplotlib.Axes, optional
            If specified, the graph will be plotted against this axis,
            otherwise a new figure will be created.
        rpath : str, optional
            Path of the directory with R (for Windows only).
        iplot : bool, optional
            If `True` plot it.
        """
        if index == -2: # apply to combined data of bigSurvey
            print('ERROR : LME survey can not be fitted on combined data.')
        elif index == -1: # apply to each
            for s in self.surveys:
                s.fitErrorLME(ax=ax, rpath=rpath, iplot=iplot)
            # redo the legend
            if ax is not None and len(self.surveys) > 1:
                self.fixLegendItems(ax)
        else:
            self.surveys[index].fitErrorLME(ax=ax, rpath=rpath, iplot=iplot)
        self.pinfo['Error model'] = True 


    def showErrorIP(self, index=0, ax=None):
        """Plot the reciprocal phase discrepancies against the reciprocal mean
        transfer resistance.

        Parameters
        ----------
        index : int, optional
            Index of the survey to show. Default is the first survey
            `index == 0`. If `ìndex == -2` then the combined data from all
            surveys are shown.
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.

        Returns
        -------
        fig : matplotlib.Figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if index == -2: # show with combined data of bigSurvey
            self.bigSurvey.showErrorIP(ax=ax)
        else:
            self.surveys[index].showErrorIP(ax=ax)


    def fitErrorPwlIP(self, index=-1, ax=None):
        """Plot the reciprocal phase errors with a power-law fit.

        Parameters
        ----------
        index : int, optional
            Index of the survey to fit. If `index == -1` (default) then the fit
            is done on all surveys independantly.
            If `ìndex == -2` then the fit is done on the combined surveys.
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.

        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if index == -2: # apply to combined data of bigSurvey
            self.bigSurvey.fitErrorPwlIP(ax=ax)
            for s in self.surveys:
                s.df['phaseError'] = self.bigSurvey.phaseErrorModel(s.df)
        elif index == -1: # apply to each
            for s in self.surveys:
                s.fitErrorPwlIP(ax=ax)
            # redo the legend
            if ax is not None and len(self.surveys) > 1:
                self.fixLegendItems(ax)
        else:
            self.surveys[index].fitErrorPwlIP(ax=ax)
        self.pinfo['Error model'] = True 


    def fitErrorParabolaIP(self, index=-1, ax=None):
        """Plot the reciprocal phase errors with a parabola fit.

        Parameters
        ----------
        index : int, optional
            Index of the survey to fit. If `index == -1` (default) then the fit
            is done on all surveys independantly.
            If `ìndex == -2` then the fit is done on the combined surveys.
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.

        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        if index == -2: # apply to combined data of bigSurvey
            self.bigSurvey.fitErrorParabolaIP(ax=ax)
            for s in self.surveys:
                s.df['phaseError'] = self.bigSurvey.phaseErrorModel(s.df)
        elif index == -1: # apply to each
            for s in self.surveys:
                s.fitErrorParabolaIP(ax=ax)
            # redo the legend
            if ax is not None and len(self.surveys) > 1:
                self.fixLegendItems(ax)
        else:
            self.surveys[index].fitErrorParabolaIP(ax=ax)
        self.pinfo['Error model'] = True 



    def showHeatmap(self, index=0, ax=None):
        """Plot a phase heatmap (x = M, y = A and value = -phi) based on:
        Orozco, A. F., K. H. Williams, and A. Kemna (2013),
        Time-lapse spectral induced polarization imaging of stimulated uranium bioremediation,
        Near Surf. Geophys., 11(5), 531–544, doi:10.3997/1873-0604.2013020)

        Parameters
        ----------
        index : int, optional
            Index of the survey to fit. Default is the first survey
            `index == 0`.
        ax : matplotlib axis, optional
            If specified, graph will be plotted on the given axis.

        Returns
        -------
        fig : matplotlib figure, optional
            If ax is not specified, the function will return a figure object.
        """
        self.surveys[index].showHeatmap(ax=ax)

    
    def checkTxSign(self):
        """Checking and correcting the polarity of the transfer resistances (flat 2D surveys only !)."""
        for s in self.surveys:
            # if np.all(s.df['resist'].values > 0): # TODO why?! 
            s.checkTxSign()


    def _filterSimilarQuad(self, quads):
        """Filter out similar quadrupole based on iselect (Survey.filterManual)
        from the specified survey.
        
        Parameters
        ----------
        quads : array
            2D array with an ABMN quadrupole per row.
        """
        totalRemoved = 0
        for s in self.surveys:
            i2keep = np.ones(s.df.shape[0], dtype=bool)
            for quad in quads:
                ie = (s.df[['a','b','m','n']].values == quad).all(1)
                i2keep = i2keep & ~ie
            s.filterData(i2keep)
            totalRemoved += np.sum(~i2keep)
        return totalRemoved


    def filterRangeIP(self, index=-1, phimin=None, phimax=None):
        """Filter IP data according to a specified range.

        Parameters
        ----------
        index : int, optional
            Index of the survey to fit. If `index == -1` (default)
            then the fit is done on all surveys independantly.
        phimin : float, optional
            Minimium phase angle [mrad].
        phimax : float, optional
            Maximum phase angle [mrad].
        """
        if index == -1: # apply to each
            for s in self.surveys:
                s.filterRangeIP(phimin, phimax)
        else:
            self.surveys[index].filterRangeIP(phimin, phimax)


    def filterRecipIP(self, index=0):
        """Remove reciprocal for IP data ONLY. Additional arguments to be
        passed to :func: `~resipy.Survey.filterRecipIP`.
        """
        if index == -2: # apply to combined data of bigSurvey
            self.bigSurvey.filterRecipIP()
            for s in self.surveys:
                s.filterRecipIP()
        elif index == -1: # apply to each
            for s in self.surveys:
                s.filterRecipIP()
        else:
            self.surveys[index].filterRecipIP()


    def filterNested(self, index=-1):
        """Removes nested measurements:
        Where M or N are in between A and B.

        Parameters
        ----------
        index : int, optional
            Index of the survey to fit. If `index == -1` (default) then the fit
            is done on all surveys independantly.
            If `ìndex == -2` then the fit is done on the combined surveys.
        """
        if index == -1: # apply to each
            for s in self.surveys:
                s.filterNested()
        else:
            self.surveys[index].filterNested()


    def addFilteredIP(self):
        """Add filtered IP to the dataframes.
        """
        for s in self.surveys:
            s.addFilteredIP()


    def filterDCA(self, index=-1, dump=None):
        """Execute DCA filtering. Decay Curve Analysis (DCA) based on.
        Flores Orozco, A., Gallistl, J., Bücker, M., & Williams, K. H. (2017).,
        Decay curve analysis for data error quantification in time-domain induced polarization imaging.,
        Geophysics, 83(2), 1–48. https://doi.org/10.1190/geo2016-0714.1

        Parameters
        ----------
        index : int, optional
            Index of the survey to use for processing. Default `index == -1`
            will apply the processing to all surveys.
        dump : function, optional
            Function onto pass the progress.
        """
        if index == -1:
            for s in self.surveys:
                s.filterDCA(dump=dump)
        else:
            self.surveys[index].filterDCA(dump=dump)


    def filterElec(self, elec=[], index=-1):
        """Filter out measurements associated with specific electrodes.

        Parameters
        ----------
        elec : list
            List of electrode number to be removed.
        index : int, optional
            Index of the survey on which to apply the processing. If the
            processing is to be applied to all surveys then specifiy
            `index=-1` (default).
        """
        numRemoved = 0
        if index == -1: # apply to all surveys
            for s in self.surveys:
                numRemoved += s.filterElec(elec)
        else:
            numRemoved = self.surveys[index].filterElec(elec)
        return numRemoved
                

    def filterRecip(self, percent=20, index=-1):
        """Filter on reciprocal errors.

        Parameters
        ----------
        index : int, optional
            Index of the survey on which to apply the processing. If the
            processing is to be applied to all surveys then specifiy
            `index=-1` (default).
        percent : float, optional
            Percentage of reciprocal error above witch a measurement will be
            discarded. 20% by default.
        """
        numRemoved = 0
        if index == -1: # apply to all surveys
            for s in self.surveys:
                numRemoved += s.filterRecip(percent)
        else:
            numRemoved = self.surveys[index].filterRecip(percent)
        return numRemoved
    
    
    def filterStack(self, percent=2, index=-1):
        """Filter on stacking (dev) errors.

        Parameters
        ----------
        index : int, optional
            Index of the survey on which to apply the processing. If the
            processing is to be applied to all surveys then specifiy
            `index=-1` (default).
        percent : float, optional
            Percentage of stacking error above witch a measurement will be
            discarded. 2% by default.
        """
        numRemoved = 0
        if index == -1: # apply to all surveys
            for s in self.surveys:
                numRemoved += s.filterStack(percent)
        else:
            numRemoved = self.surveys[index].filterStack(percent)
        return numRemoved


    def filterUnpaired(self, index=-1):
        """Remove quadrupoles that don't have reciprocals. This might
        remove dummy measurements added for sequence optimization.

        Parameters
        ----------
        index : int, optional
            Index of the survey on which to apply the processing. If the
            processing is to be applied to all surveys then specifiy
            `index=-1` (default).
        """
        numRemoved = 0
        if index == -1: # apply to all surveys
            for s in self.surveys:
                numRemoved += s.filterUnpaired()
        else:
            numRemoved = self.surveys[index].filterUnpaired()
        return numRemoved



    def filterNegative(self):
        """Remove negative apparent resistivty values
        """
        for s in self.surveys:
            s.filterNegative()
            
    
    def filterAppResist(self, index=-1, vmin=None, vmax=None):
        """Filter measurements by apparent resistivity for surface surveys 
        
        Parameters
        ----------
        vmin : float, optional
            Minimum value.
        vmax : float, optional
            Maximum value.
        index : int, optional
            Index of the survey on which to apply the processing. If the
            processing is to be applied to all surveys then specifiy
            `index=-1` (default).
        """
        numRemoved = 0
        if index == -1: # apply to all surveys
            for s in self.surveys:
                numRemoved += s.filterAppResist(vmin=vmin, vmax=vmax)
        else:
            numRemoved = self.surveys[index].filterAppResist(vmin=vmin, vmax=vmax)
        return numRemoved


    def filterTransferRes(self, index=-1, vmin=None, vmax=None):
        """Filter measurements by transfer resistance. 
        
        Parameters
        ----------
        vmin : float, optional
            Minimum value.
        vmax : float, optional
            Maximum value.
        index : int, optional
            Index of the survey on which to apply the processing. If the
            processing is to be applied to all surveys then specifiy
            `index=-1` (default).
        """
        numRemoved = 0
        if index == -1: # apply to all surveys
            for s in self.surveys:
                numRemoved += s.filterTransferRes(vmin=vmin, vmax=vmax)
        else:
            numRemoved = self.surveys[index].filterTransferRes(vmin=vmin, vmax=vmax)
        return numRemoved
    
    
    def filterContRes(self, index=-1, vmin=None, vmax=None):
        """Filter measurements by contact resistance if available. 
        
        Parameters
        ----------
        vmin : float, optional
            Minimum value in unit of cR column.
        vmax : float, optional
            Maximum value in unit of cR column.
        index : int, optional
            Index of the survey on which to apply the processing. If the
            processing is to be applied to all surveys then specifiy
            `index=-1` (default).
        """
        numRemoved = 0
        if index == -1: # apply to all surveys
            for s in self.surveys:
                numRemoved += s.filterContRes(vmin=vmin, vmax=vmax)
        else:
            numRemoved = self.surveys[index].filterContRes(vmin=vmin, vmax=vmax)
        return numRemoved
    
    def computeReciprocal(self,alg='Bisection Search',forceSign=False):
        """
        Compute Reciprocals and store them in self.surveys[0:n].df. 

        Parameters
        ----------
        alg : str, optional
            Algorithm used to compute reciprocals. Choose between 
            'Bisection Search', 'Pandas Merge' or 'Array Expansion'.
            The default is 'Bisection Search', other string are casted to
            'Array Expansion'.
        forceSign: bool, optional
            Force reciprocal and forward measurements to have the same 
            polarity regarding the calculation of the reciprocal errors. 
            Default is False. 
            
        Notes
        -----
        If the relevant cython code cannot be found then the function falls 
        back to using the pandas merge approach. 
        """
        for s in self.surveys:
            s.computeReciprocal(alg,forceSign)
            
    def computeFineMeshDepth(self):
        """Compute the Fine Mesh Depth (FMD) based on electrode
        positions and the larger dipole spacing. Express as a positive number,
        it represents the relative vertical distance to extend the fine mesh
        region below the surface.
        """
        dfelec = self.elec.copy()
        dfelec = dfelec[~self.elec['remote']] # discard remote electrode for this
        elec = dfelec[['x','y','z']].values
        if (self.typ == 'R2') | (self.typ == 'cR2'): # 2D survey:
            if (len(self.surveys) > 0):# & (self.iForward == False):
                lookupDict = dict(zip(dfelec['label'], np.arange(dfelec.shape[0])))
                array = self.surveys[0].df[['a','b','m','n']].replace(lookupDict).values.copy().astype(int) # strings don't have max/min
                maxDist = np.max(np.abs(elec[array[:,0]-np.min(array[:,0]),0] - elec[array[:,2]-np.min(array[:,2]),0])) # max dipole separation
                self.fmd = (1/3)*maxDist
            else: # if it's a forward model for instance
                self.fmd = (1/3)*(np.max(elec[:,0]) - np.min(elec[:,0]))

        else: # for 3D survey
            dist = mt.findDist(elec[:,0], elec[:,1], elec[:,2])
            self.fmd = (1/3)*np.max(dist)

        if self.elec['buried'].sum() > 0:
            # catch where buried electrodes are present as the fmd needs adjusting in this case 
            if self.elec['buried'].sum() == self.elec.shape[0]:
                # find max distance 
                dist = mt.findDist(elec[:,0], elec[:,1], elec[:,2])
                self.fmd = (1/3)*np.max(dist)
            else: # surface given by max z elec
                self.fmd = np.abs(np.max(elec[:,2])  - np.min(elec[:,2])) + (0.5*self.fmd)
        # print('Fine Mesh Depth (relative to the surface): {:.2f} m'.format(self.fmd))


    def createMesh(self, typ='default', buried=None, surface=None, cl_factor=2,
                   cl=-1, dump=None, res0=100, show_output=False, fmd=None,
                   remote=None, refine=0, **kwargs):
        """Create a mesh.
        
        Parameters
        ----------
        typ : str, optional
            Type of mesh.
            For 2D:
                - 'quad': quadrilateral (fast to build)
                - 'trian': triangular (fast to run) - default
                - 'circle': close circular mesh
            For 3D: 
                - 'tetra': tetrahedral mesh for half-space - default
                - 'cylinder': column build with tetrahedral
                - 'prism': column build with prism
                - 'tank': closed geometry with tetrahedra
        buried : numpy.array, optional
            Boolean array of electrodes that are buried. Should be the same
            length as `R2.elec`
        surface : numpy.array, optional
            Array with two or three columns x, y (optional) and elevation for
            additional surface points.
        cl_factor : float, optional
            Characteristic length factor. Only used for triangular mesh to allow
            mesh to be refined close the electrodes and then expand.
        cl : float, optional
            Characteristic length that define the mesh size around the
            electrodes.
        dump : function, optional
            Function to which pass the output during mesh generation. `print()`
             is the default.
        res0 : float, optional
            Starting resistivity for mesh elements.
        show_output : bool, optional
            If `True`, the output of gmsh will be shown on screen.
        fmd : float, optional
            Depth of fine region specifies as a positive number relative to the mesh surface.
        remote : bool, optional
            Boolean array of electrodes that are remote (ie not real). Should be the same
            length as `Project.elec`.
        refine : int, optional
            Number times the mesh will be refined. Refinement split the triangles
            or the tetrahedra but keep the same number of parameter for the inversion.
            This helps having a more accurate forward response and a faster inversion
            (as the number of elements does not increase). Only available for
            triangles or tetrahedral mesh.
        kwargs : -, optional
            Keyword arguments to be passed to mesh generation schemes
            Specific for 'tank mesh':
                - origin : list of 3 floats
                    Origin in X,Y,Z of one of the tank corner.
                - dimension : list of 3 floats
                    Dimension from the corner on how to extend the tank.
            Specific for 'cylinder mesh':
                - zlim : list of 2 int
                For the bottom and top of the column along the Z axis.
        """     
        if dump is None:
            if show_output:
                def dump(x):
                    print(x, end='')
            else:
                def dump(x):
                    pass
                
        if typ == 'default':
            if self.typ == 'R2' or self.typ == 'cR2': # it's a 2D mesh
                typ = 'trian'
            else:
                typ = 'tetra'
                
        self.meshParams = {'typ':typ, 'buried':buried, 'surface':surface,
                           'cl_factor':cl_factor, 'cl':cl, 'dump':dump,
                           'res0': res0, 'show_output':show_output,
                           'refine':refine,'fmd':fmd}
        for key in kwargs.keys():
            self.meshParams[key] = kwargs[key]
            
        meshtypename = 'default'
        if kwargs is not None:
            self.meshParams.update(kwargs)
        
        # flag for if mesh gets refined during mesh creation 
        refined = False 

        # define electrode types
        if buried is not None:
            if len(buried) == self.elec.shape[0]:
                self.elec['buried'] = buried
            else:
                print('length of argument "buried" ({:s}) does not match length'
                      ' of self.elec ({:d})'.format(len(buried), self.elec.shape[0]))
                
        elec_x = self.elec['x'].values
        elec_y = self.elec['y'].values
        elec_z = self.elec['z'].values
        elec_type = np.repeat('electrode',len(elec_x))
        elec_type[self.elec['buried'].values] = 'buried'
        elec_type[self.elec['remote'].values] = 'remote'
        elecLabels = self.elec['label'].values.astype(str)
        
        
        # assign possible topography (surface)
        if surface is not None:
            if surface.shape[1] == 2:
                self.topo = pd.DataFrame(surface, columns=['x','z'])
            else:
                self.topo = pd.DataFrame(surface, columns=['x','y','z'])
                
            # error checking on surface 
            for column in self.topo.columns:
                if any(np.isnan(self.topo[column].values)):
                    raise Exception('There are nan values in surface topography %s column!'%column)
            
            if self.coordLocal: 
                # need to convert surface to local coordinates 
                surface = surface.copy() 
                if surface.shape[1] == 2: 
                    utmy = np.zeros(surface.shape[0])
                else: 
                    utmy = surface[:,1]
                utmx = surface[:,0]
                x0 = self.coordParam['x0']
                y0 = self.coordParam['y0']
                a = self.coordParam['a']
                localx, localy = rotGridData(utmx, utmy, x0, y0, a)
                surface[:,0] = localx 
                if surface.shape[1] == 3:
                    surface[:,1] = localy 
                    # reassign topo 
                    self.topo = pd.DataFrame(surface, columns=['x','y','z'])
                else:
                    self.topo = pd.DataFrame(surface, columns=['x','z'])

        
        # estimate depth of fine mesh
        if fmd is None:
            self.computeFineMeshDepth()
        else:
            self.fmd = fmd
        
        if typ == 'quad':
            print('Creating quadrilateral mesh...', end='')
            if 'fmd' not in kwargs.keys():
                kwargs['fmd'] = self.fmd 
            meshtypename = 'Quadrilateral'
            surface_x = self.topo['x'].values if surface is not None else None
            surface_z = self.topo['z'].values if surface is not None else None
            mesh,meshx,meshy,topo,e_nodes = mt.quadMesh(elec_x,elec_z,list(elec_type),
                                                         surface_x=surface_x, surface_z=surface_z,
                                                         **kwargs)   #generate quad mesh
            self.param['mesh_type'] = 6
            e_nodes = np.array(mesh.eNodes) + 1 # +1 because of indexing staring at 0 in python
            self.param['node_elec'] = [elecLabels, e_nodes.astype(int)]

            if 'regions' in self.param: # allow to create a new mesh then rerun inversion
                del self.param['regions']
            if 'num_regions' in self.param:
                del self.param['num_regions']
        else:
            geom_input = {}
            
            if surface is not None:
                geom_input['surface'] = [self.topo['x'].values,
                                         self.topo['z'].values]
                
            if 'geom_input' in kwargs:
                geom_input.update(kwargs['geom_input'])
                kwargs.pop('geom_input')
            
            if self.elec['buried'].values.all() and surface is None:
                # all electrodes buried and no surface given
                self.wholespace = True
            else:
                self.wholespace = False 

            elec_type = elec_type.tolist()

            def setMeshProc(a): # little function to be passed to meshTools fct
            # to get the variable outputed by Popen (allow mesh killing in UI)
                self.mproc = a
                
            with cd(self.dirname):#change to working directory so that mesh files written in working directory
                self.mproc = None
                if typ == 'trian':
                    print('Creating triangular mesh...', end='')
                    meshtypename = 'Triangular'
                    mesh = mt.triMesh(elec_x,elec_z,elec_type,geom_input,
                                 path=os.path.join(self.apiPath, 'exe'),
                                 cl_factor=cl_factor,
                                 cl=cl, dump=dump, show_output=show_output,
                                 fmd=self.fmd, whole_space=self.wholespace,
                                 handle=setMeshProc, **kwargs)
                elif typ == 'circle':
                    print('Creating circular mesh...NOT IMPLEMENTED YET', end='')
                    meshtypename = 'Triangular'
                    mesh = mt.circularMesh(np.c_[elec_x, elec_y, elec_z],
                                             path=os.path.join(self.apiPath, 'exe'),
                                             cl=cl, dump=dump, show_output=show_output,
                                             handle=setMeshProc, **kwargs)
                    self.param['num_xy_poly'] = 0
                    
                elif typ == 'tetra':
                    print('Creating tetrahedral mesh...', end='')    
                    meshtypename = 'Tetrahedral'
                    if cl == -1:
                        dist = cdist(self.elec[~self.elec['remote']][['x','y']].values) # half the minimal electrode distance
                        cl = np.min(dist[dist != 0])/4 
                    mesh = mt.tetraMesh(elec_x, elec_y, elec_z,elec_type,
                                 path=os.path.join(self.apiPath, 'exe'),
                                 surface_refinement=surface,
                                 cl_factor=cl_factor,
                                 cl=cl, dump=dump, show_output=show_output,
                                 fmd=self.fmd, whole_space=self.wholespace,
                                 handle=setMeshProc, **kwargs)
                elif typ == 'prism':
                    print('Creating prism mesh...', end='')
                    meshtypename = 'Prism'
                    mesh = mt.prismMesh(elec_x, elec_y, elec_z,
                                         path=os.path.join(self.apiPath, 'exe'),
                                         cl=cl, dump=dump, show_output=show_output,
                                         handle=setMeshProc, **kwargs)
                    self.param['num_xz_poly'] = 0
                    # self.iadvanced = False # as current implimentation of advanced mesh doesnt work with prisms 
                elif typ == 'cylinder':
                    print('Creating cylinder mesh...', end='')
                    meshtypename = 'Tetrahredal (cylindrical)'
                    mesh = mt.cylinderMesh(elec_x, elec_y, elec_z,
                                             path=os.path.join(self.apiPath, 'exe'),
                                             cl=cl, cl_factor=cl_factor, dump=dump, 
                                             show_output=show_output,
                                             handle=setMeshProc, **kwargs)
                    self.param['num_xz_poly'] = 0
                    # self.iadvanced = False # ditto 
                elif typ == 'tank':
                    print('Creating tank mesh...', end='')
                    meshtypename = 'Tetrahedral (tank)'
                    mesh = mt.tankMesh(np.c_[elec_x, elec_y, elec_z],
                                             path=os.path.join(self.apiPath, 'exe'),
                                             cl=cl, dump=dump, show_output=show_output,
                                             handle=setMeshProc, **kwargs)
                    self.param['num_xz_poly'] = 0
            self.mproc = None # mesh successfully done so let's put this back to None in case it isn't already

            # mesh refinement
            if (typ != 'prism') and (typ != 'quad'):
                for l in range(refine):
                    print('refining...', end='')
                    mesh = mesh.refine()
                refined = True 
                
            self.param['mesh_type'] = 3
            e_nodes = np.array(mesh.eNodes) + 1 # +1 because of indexing staring at 0 in python
            self.param['node_elec'] = [elecLabels, e_nodes.astype(int)]

        self.mesh = mesh
        self.param['mesh'] = mesh
        self.param['num_regions'] = 0
        self.param['res0File'] = 'res0.dat'

        numel = self.mesh.numel
        self.mesh.addAttribute(np.ones(numel)*res0, 'res0') # default starting resisivity [Ohm.m]
        self.mesh.addAttribute(np.ones(numel)*0, 'phase0') # default starting phase [mrad]
        self.mesh.addAttribute(np.ones(numel, dtype=int), 'zones')
        self.mesh.addAttribute(np.zeros(numel, dtype=float), 'iter')
        if not refined: # dont want to re-allocate parameter number it has already been assigned  
            ## is bugged in 2D, needs fixing! ## 
            self.mesh.addAttribute(np.arange(numel)+1,'param') # param = 0 if fixed
        elif self.typ=='R2':
            self.mesh.addAttribute(np.arange(numel)+1,'param')
        
        # check ram 
        memInv = self._estimateMemory(dump=dump)
        self.param['reqMemory'] = getSysStat()[2] - memInv # if negative then we need more RAM
        
        # deal with remote electrodes in mesh object 
        self.mesh.iremote = self.elec['remote'].values
        
        # define zlim
        if (typ == 'tank') | (typ == 'cylinder') | (typ == 'prism'):
            zmin = np.min(self.mesh.node[:,2])
            zmax = np.max(self.mesh.node[:,2])
            extent = zmax - zmin
            self.zlim = [zmin - 0.1*extent, zmax + 0.1*extent]
        else:
            if surface is not None:
                zlimTop = np.max([np.max(elec_z), np.max(surface[:,-1])])
            else:
                zlimTop = np.max(elec_z)
            if all(self.elec['buried']) and surface is None: # whole mesh
                zlimBot = np.min(elec_z)
            elif any(self.elec['buried']):
                if surface is None:
                    zlimBot = np.min(elec_z[~self.elec['buried']]) - self.fmd 
                else:
                    zlimBot = np.min(self.topo['z'].values) - self.fmd 
            else:
                zlimBot = np.min(elec_z[~self.elec['remote']]) - self.fmd 
                
            self.zlim = [zlimBot, zlimTop]
        self._computePolyTable()
        print('done ({:d} elements)'.format(self.mesh.df.shape[0]))
        self.pinfo['Number of Elements'] = self.mesh.numel 
        self.pinfo['Number of Nodes']=self.mesh.numnp 
        self.pinfo['Mesh Type'] = meshtypename
        self.fwdErrModel = False # reset model error as mesh has changed 
        
        
    def _computePolyTable(self):
        # define num_xz_poly or num_xy_poly
        elec = self.elec[~self.elec['remote']][['x','y','z']].values
        elec_x, elec_y, elec_z = elec[:,0], elec[:,1], elec[:,2]
        self.param['num_xz_poly'] = 0 
        self.param['num_xy_poly'] = 0 
        
        if (self.typ == 'R2') | (self.typ == 'cR2'):
            self.param['num_xz_poly'] = 5
            if all(self.elec['buried']): # we don't know if there is a surface
            # or not so we trust the zlim computation done in createMesh()
                zmax = self.zlim[1]
                zmin = self.zlim[0]
            else:
                zmax = np.max(elec_z)
                zmin = np.min(elec_z) - self.fmd 
            xmin, xmax = np.min(elec_x), np.max(elec_x)
            xz_poly_table = np.array([
            [xmin, zmax],
            [xmax, zmax],
            [xmax, zmin],
            [xmin, zmin],
            [xmin, zmax]])
            self.param['xz_poly_table'] = xz_poly_table
        elif all(np.array(elec_y)==elec_y[0]): #catch 2d line case 
            self.param['num_xy_poly'] = 5
            xmin, xmax = np.min(elec_x), np.max(elec_x)
            ymin, ymax = np.min(elec_y)-(self.fmd/5), np.max(elec_y)+(self.fmd/5)
            zmin, zmax = np.min(elec_z)-self.fmd, np.max(elec_z)
            xy_poly_table = np.array([
            [xmin, ymax],
            [xmax, ymax],
            [xmax, ymin],
            [xmin, ymin],
            [xmin, ymax]])
            self.param['zmin'] = zmin
            self.param['zmax'] = zmax
            self.param['xy_poly_table'] = xy_poly_table
        else:
            self.param['num_xy_poly'] = 5
            xmin, xmax = np.min(elec_x), np.max(elec_x)
            ymin, ymax = np.min(elec_y), np.max(elec_y)
            zmin, zmax = np.min(elec_z)-self.fmd, np.max(elec_z)
            xy_poly_table = np.array([
            [xmin, ymax],
            [xmax, ymax],
            [xmax, ymin],
            [xmin, ymin],
            [xmin, ymax]])
            self.param['zmin'] = zmin
            self.param['zmax'] = zmax
            self.param['xy_poly_table'] = xy_poly_table

    
    def _defineZlim(self):
        """Computes zlim"""
        if self.fmd is None:
            self.computeFineMeshDepth()
        zlimMax = self.elec['z'].max()
        zlimMin = self.elec['z'].min() - self.fmd
        self.zlim = [zlimMin, zlimMax]
        

    def importMesh(self, file_path, mesh_type=None, node_pos=None, elec=None,
                   order_nodes=True, res0=100):
        """Import mesh from .vtk / .msh / .dat, rather than having ResIPy
        create one for you.

        Parameters
        ----------
        file_path : str
            File path mapping to the mesh file
        mesh_type : str
            Not used anymore. 
        node_pos : array like, optional
            Array of ints referencing the electrode nodes. If left as none no electrodes
            will be added to the mesh class. Consider using mesh.moveElecNodes()
            to add nodes to mesh using their xyz coordinates.
        elec : array, optional
            N*3 numpy array of electrode x,y,z coordinates. Electrode node positions
            will be computed by finding the nearest nodes to the relevant coordinates.
        res0 : float, optional
            Starting resistivity for mesh elements.
        """
        if (self.typ == 'R3t') or (self.typ == 'cR3t'):
            flag_3D = True
        else:
            flag_3D = False
        self.mesh = mt.readMesh(file_path, node_pos=node_pos, 
                                          order_nodes=order_nodes)

        self.meshParams = {'typ':'imported'}

        # recover region based on resistivity
        if file_path[-4:] == '.vtk':
            if 'Resistivity(ohm.m)' in self.mesh.df.columns:
                res = self.mesh.df['Resistivity(ohm.m)']
                ures = np.unique(res)
                for i, reg in enumerate(ures):
                    ie = res == reg
                    self.mesh.df.loc[ie, 'region'] = i+1
            
        if elec is not None:
            self.mesh.moveElecNodes(elec[:,0], elec[:,1], elec[:,2])

        #add the electrodes to the R2 class
        if elec is not None or node_pos is not None: # then electrode positions should be known
            self.setElec(self.mesh.elec)
        else:
            try:
                elec = self.elec[['x','y','z']].values
                self.mesh.moveElecNodes(elec[:,0],elec[:,1],elec[:,2])
            except AttributeError:
                warnings.warn("No electrode nodes associated with mesh! Electrode positions are unknown!")

        #R2 class mesh handling
        e_nodes = np.array(self.mesh.eNodes) + 1 # +1 because of indexing staring at 0 in python
        self.param['mesh'] = self.mesh
        # if mesh_type == 'quad':
        #     self.param['mesh_type'] = 6
        #     colx = self.mesh.quadMeshNp() # convert nodes into column indexes
        #     self.param['node_elec'] = np.c_[1+np.arange(len(e_nodes)), np.array(colx), np.ones((len(e_nodes,1)))].astype(int)
            #will only work for assuming electrodes are a surface array
        meshtypename = 'default'
        if self.mesh.type2VertsNo() == 4:
            if flag_3D:
                self.param['mesh_type'] = 4 # tetra mesh
                meshtypename = 'Tetrahedral'
            else:
                self.param['mesh_type'] = 6 # general quad mesh
                meshtypename = 'Quadrilateral'
        else:
            if flag_3D:
                self.param['mesh_type'] = 6 # prism mesh 
                meshtypename = 'Prism'
            else:
                self.param['mesh_type'] = 3 # triangular mesh
                meshtypename = 'Triangular'
        self.param['node_elec'] = [self.elec['label'].values, e_nodes.astype(int)]

        # checking
        if len(np.unique(e_nodes)) < len(e_nodes):
            warnings.warn('Some electrodes are positionned on the same nodes : e_nodes=' + str(e_nodes))
        
        # make regions continuous
        regions = self.mesh.df['region']
        uregions = np.unique(regions)
        iregions = np.arange(len(uregions)) + 1
        dico = dict(zip(uregions, iregions))
        self.mesh.df['region'] = [dico[a] for a in regions]
        
        self.param['num_regions'] = 0
        self.param['res0File'] = 'res0.dat'
        numel = self.mesh.numel
        if 'res0' not in self.mesh.df.columns:
            self.mesh.df['res0'] = np.ones(numel)*res0 # default starting resisivity [Ohm.m]
        if 'phase0' not in self.mesh.df.columns:
            self.mesh.df['phase0'] = np.ones(numel)*0 # default starting phase [mrad]
        if 'zones' not in self.mesh.df.columns:
            self.mesh.df['zones'] = np.ones(numel, dtype=int)
        else:
            self.mesh.df['zones'] = self.mesh.df['zones'].astype(int) # for some reason read as float
        if 'iter' not in self.mesh.df.columns:
            self.mesh.df['iter'] = np.zeros(numel, dtype=float)
        self.mesh.iremote = self.elec['remote'].values

        if self.typ == 'R3t' or self.typ == 'cR3t':
            self.mesh.datAdv(os.path.join(self.dirname, 'mesh3d.dat'), iadvanced=self.iadvanced)
        else:
            self.mesh.dat(os.path.join(self.dirname, 'mesh.dat'))

        # define zlim
        if self.zlim is None:
            self._defineZlim()
        self._computePolyTable()
        
        #flag up mesh in pinfo 
        self.pinfo['Number of Elements'] = self.mesh.numel 
        self.pinfo['Number of Nodes']=self.mesh.numnp 
        self.pinfo['Mesh Type'] = meshtypename
        
        _ = self._estimateMemory(dump=donothing)
        
        

    def showMesh(self, ax=None, **kwargs):
        """Display the mesh and handles best default values. 
        """
        if self.mesh is None:
            raise Exception('Mesh undefined')
        else:
            elec = self.elec[['x','y']].values
            elec = elec[~self.elec['remote'].values,:]
            if 'zlim' not in kwargs.keys():
                kwargs['zlim'] = self.zlim
             
            # best to use x/y limits provided by project class becuase it knows if electrodes are remote 
            if 'xlim' not in kwargs.keys():
                # scale up electrode coordinates by 20% to give some padding around where the electrodes are 
                kwargs['xlim'] = findminmax(elec[:,0]) 
                # can use surface topography too if that has been provided 
                if len(self.topo) > 4:
                    kwargs['xlim'] = findminmax(self.topo.x,1)
            if 'ylim' not in kwargs.keys() and self.typ[-1] == 't':
                kwargs['ylim'] = findminmax(elec[:,1])
                if len(self.topo) > 4:
                    kwargs['ylim'] = findminmax(self.topo.y,1)
                
            if 'color_map' not in kwargs.keys(): # pick a color map based on display type
                if self.typ[-1] == 't':
                    kwargs['color_map'] = 'Greys'
                else:
                    kwargs['color_map'] = 'gray'
                    
            if 'attr' not in kwargs.keys():
                kwargs['attr'] = 'region' # this will show regions by default
                
            if 'color_bar' not in kwargs.keys():
                if np.unique(np.array(self.mesh.df['region'])).shape[0] > 1:
                    kwargs['color_bar'] = True # show colorbar for multiple regions
                    if kwargs['attr'] == 'region' and kwargs['color_map'] == 'Greys':
                        kwargs['color_map'] = 'Spectral'
                else:
                    kwargs['color_bar'] = False
            
            if 'typ' in self.meshParams.keys():
                for a in ['tank', 'cylinder', 'prism']:
                    if self.meshParams['typ'] == a:
                        kwargs['clipping'] = False 
        
            self.mesh.show(ax=ax, darkMode=self.darkMode, **kwargs)
            
    def refineMesh(self):
        self.mesh = self.mesh.refine()
        
    def res0fromPseudo(self):
        """
        Generate initial starting model from pseudo section.

        Returns
        -------
        res0: array like 
            Initial resistivity for each cell in the mesh. 

        """
        if len(self.surveys) == 0:
            warnings.warn('Cannot set res0 from pseudo section without first loading in data')
            return 
        if self.mesh is None:
            warnings.warn('Cannot set res0 from pseudo section without first creating a mesh')
            return 
        px, py, pz = self.surveys[0]._computePseudoDepth(flag3d=True)
        df = self.surveys[0].df
        app = df['K'].values*df['resist'].values
        
        mx, my, mz = self.mesh.df.X, self.mesh.df.Y, self.mesh.df.Z
        # xnew,ynew,znew,xknown, yknown, zknown, iknown
        res0 = nearest3d(mx, my, mz, px, py, pz, app)
        res0[res0 < 0] = 1.0 # cap minimum resistivity value 
        # self.mesh.df.loc[:,'res0'] = res0 #  perhaps dont assign res0 here? 
        return res0 
        

    def write2in(self, param={}, err=None):
        """Create configuration file for inversion. Write mesh.dat and res0.dat.

        Parameters
        ----------
        param : dict
            Dictionnary of parameters and values for the inversion settings.
        err: bool
            Used to overide self.err. 
        """
        typ = self.typ
        if err is None: 
            err = self.err 
        if (err is True) and ('a_wgt' not in self.param):
            self.param['a_wgt'] = 0
            self.param['b_wgt'] = 0
        elif (typ == 'R2') or (typ == 'R3t'): # DC case
            if 'a_wgt' not in self.param:
                self.param['a_wgt'] = 0.001
            if 'b_wgt' not in self.param:
                self.param['b_wgt'] = 0.02
        elif (typ == 'cR2') | (typ == 'cR3t'): # IP case
            if 'a_wgt' not in self.param:
                self.param['a_wgt'] = 0.02 # variance for magnitude (no more offset)
            if 'b_wgt' not in self.param:
                self.param['b_wgt'] = 2 # mrad
        
        #catch infinite z limits 
        if 'zmin' in self.param and self.param['zmin'] == np.inf:
            self.param['zmin'] = np.min(self.mesh.node[:,2]) - 10 
        if 'zmax' in self.param and self.param['zmax'] == np.inf:
            self.param['zmax'] = np.max(self.mesh.node[:,2]) + 10 

        # all those parameters are default but the user can change them and call
        # write2in again
        for p in param:
            self.param[p] = param[p]

        if self.iTimeLapse == True:
            # set background survey parameters first
            refdir = os.path.join(self.dirname, 'ref')
            if os.path.exists(refdir) == False:
                os.mkdir(refdir)
            param = self.param.copy()
            if err:
                param['a_wgt'] = 0
                param['b_wgt'] = 0
            else: # default DC case as timelapse not supported for IP yet
                if 'a_wgt' not in param:#this allows previously assigned values to be
                    param['a_wgt'] = 0.001 # written to the reference.in config file
                if 'b_wgt' not in param:
                    param['b_wgt'] = 0.02
            param['reg_mode'] = 0 # set by default in ui.py too
            param['res0File'] = 'res0.dat'
            if (self.typ == 'R2') or (self.typ == 'cR2'):
                param['num_xz_poly'] = 0
            else:
                param['num_xy_poly'] = 0
                param['inverse_type'] = 0 # normal regularisation
                param['zmin'] = np.min(self.mesh.node[:,2]) - 10 # we want to keep the whole mesh for background regularisation
                param['zmax'] = np.max(self.mesh.node[:,2]) + 10
            
            if 'baseline_target_decrease' in param.keys():
                param['target_decrease'] = param['baseline_target_decrease']
                
            self.configFile = write2in(param, refdir, typ=typ) # background survey
            
            # now prepare the actual timelapse settings
            self.param['num_regions'] = 0
            if 'reg_mode' not in self.param.keys():
                self.param['reg_mode'] = 2
            self.param['res0File'] = 'Start_res.dat'
            write2in(self.param, self.dirname, typ=typ) # actual time-lapse
        else:
            self.configFile = write2in(self.param, self.dirname, typ=typ)

        # writing mesh.dat
        ifixed = np.array(self.mesh.df['param']) == 0
        if np.sum(ifixed) > 0: # fixed element need to be at the end
            self.mesh.orderElem()
            self.mesh.writeAttr('res0', os.path.join(self.dirname, 'Start_res.dat')) # needs rewriting if the elements are reordered. 
            
        if typ == 'R3t' or typ == 'cR3t':
            self.mesh.datAdv(os.path.join(self.dirname, 'mesh3d.dat'), iadvanced=self.iadvanced)
        else:
            self.mesh.dat(os.path.join(self.dirname, 'mesh.dat'))
        
        # write the res0.dat needed for starting resistivity
        if self.iForward is True: # we will invert results from forward
            # inversion so we need to start from a homogeneous model
            print('All non fixed parameters reset to 100 Ohm.m and 0 mrad, '
                  'as the survey to be inverted is from a forward model.')
            ifixed = self.mesh.df['param'] == 0
            res0 = np.array(self.mesh.df['res0'])
            phase0 = np.array(self.mesh.df['phase0'])
            res0f = res0.copy()
            phase0f = phase0.copy()
            res0f[~ifixed] = 100
            phase0f[~ifixed] = 0
            self.mesh.df['res0'] = list(res0f)
            self.mesh.df['phase0'] = list(phase0f)

        if (self.typ == 'cR2') or (self.typ == 'cR3t'):
            r = np.array(self.mesh.df['res0'])
            phase = np.array(self.mesh.df['phase0'])
            centroids = self.mesh.elmCentre.copy()
            centroids2 = centroids[:,[0,2]] if self.typ[-1] != 't' else centroids
            x = np.c_[centroids2,
                      r,
                      phase, # mrad
                      np.log10(r),
                      np.log10(np.cos(-phase/1000)/np.log10(r)), #log10(real conductivity)
                      np.log10(np.sin(-phase/1000)/np.log10(r))] #log10(imaginary conductivity)
            np.savetxt(os.path.join(self.dirname, 'res0.dat'), x)
        else:
            self.mesh.writeAttr('res0', os.path.join(self.dirname, 'res0.dat'))
        
        if self.iForward: # restore initial res0 and phase0 so that user can 
        # rerun the forward model with a different sequence for instance
            self.mesh.df['res0'] = list(res0)
            self.mesh.df['phase0'] = list(phase0)
            

    def write2protocol(self, err=None, errTot=False, fm0=None, **kwargs):
        """Write a protocol.dat file for the inversion code.

        Parameters
        ----------
        err : bool, optional
            If `True` error columns will be written in protocol.dat provided
            an error model has been fitted or errors have been imported.
        errTot : bool, optional
            If `True`, it will compute the modelling error due to the mesh and
            add it to the error from an error model.
        fm0 : numpy.array of float, optional
            Only for 3D time-lapse with reg_mode == 2 (difference inversion).
            Response of the inverted reference survey according to sequence of
            the reference survey as transfer resistance (Ohm).
        **kwargs : optional
            To be passed to `Survey.write2protocol()`.
        """
        if self.typ[0] == 'c':
            ipBool = True
        else:
            ipBool = False

        if self.typ == 'R2' or self.typ == 'cR2':
            threed = False
        else:
            threed = True
            
        if err is None:
            err = self.err
            
        errTyp = self.errTyp # either 'global' (default) or 'survey'

        # important changing sign of resistivity and quadrupoles so to work
        # with complex resistivity
        if self.typ[0] == 'c':
            for s in self.surveys:
                ie = s.df['resist'].values < 0
                m = s.df['m'].values.copy()
                n = s.df['n'].values.copy()
                s.df.loc[ie, 'm'] = n[ie]
                s.df.loc[ie, 'n'] = m[ie]
                # let's change the sign as cR2 will take the log of it anyway
                # and we are dealing with a magnitude here, not a resistivity
                s.df.loc[ie, 'resist'] = s.df.loc[ie, 'resist'].values*-1
                s.df.loc[ie, 'recipMean'] = s.df.loc[ie, 'recipMean'].values*-1
        
        # check transfer resistance sign
        if 'checkTxSign' in self.param.keys() and self.param['checkTxSign'] is True:
            self.checkTxSign()
            
        # for time-lapse inversion ------------------------------
        if self.iTimeLapse is True:
            if 'reg_mode' not in self.param.keys():
                self.param['reg_mode'] = 2 # by default it's timelapse (difference)
            if self.param['reg_mode'] == 2: # it's a difference inversion
                indexes = self.matchSurveys()
                # if (self.typ == 'R3t') | (self.typ == 'cR3t'):
                #     if fm0 is not None:
                #         fm0 = fm0[indexes[0] # we crop it so it has the same
                #         # shape as all quad in common --> done when writing it below
            else:
                indexes = [(None, None)]*len(self.surveys)
            # a bit simplistic but assign error to all based on Transfer resistance
            # let's assume it's False all the time for now
            content = ''
            df0 = self.surveys[0].df[['a','b','m','n','resist','recipMean','resError']]
            df0 = df0.rename(columns={'resist':'resist0', 'recipMean':'recipMean0','resError':'resError0'})
            for i, s in enumerate(self.surveys):
                if 'resist0' in s.df.columns:
                    s.df = s.df.drop('resist0', axis=1)
                if 'recipMean0' in s.df.columns:
                    s.df = s.df.drop('recipMean0', axis=1)
                if self.param['reg_mode'] == 2 and i > 0:
                    # modelling error already accounted for in the baseline model! 
                    errTot = False 
                    
                s.df = pd.merge(s.df, df0, on=['a','b','m','n'], how='left')
                # resError and phaseError should already have been populated
                # handle the case when SOME survey were fitted but not all
                # then we use the bigSurvey default fit to fullfill them
                if err is True and errTyp == 'global':
                    if np.sum(np.isnan(s.df['resError'])) != 0: # there is some nan
                        print('Survey {:s} has no fitted error model, default to combined fit.'.format(s.name))
                        if self.bigSurvey.errorModel is None:
                            self.bigSurvey.fitErrorPwl() # default fit
                        s.df['resError'] = self.bigSurvey.errorModel(s.df)
                    if self.typ[0] == 'c' and np.sum(np.isnan(s.df['phaseError'])) != 0: # there is some nan
                        print('Survey {:s} has no fitted IP error model, default to combined fit.'.format(s.name))
                        if self.bigSurvey.errorModel is None:
                            self.bigSurvey.fitErrorPwlIP()
                        s.df['phaseError'] = self.bigSurvey.phaseErrorModel(s.df)
                        
                # if not it means that the 'resError' columns has already
                # been populated when the files has been imported

                res0Bool = False if self.param['reg_mode'] == 1 else True
                # use the measurements that are matched with both the reference and timelapse inversion 
                fm00 = fm0[indexes[i][0]] if fm0 is not None else fm0 
                protocol = s.write2protocol('', err=err, errTot=errTot, res0=res0Bool,
                                            ip=False, # no IP timelapse possible for now
                                            isubset=indexes[i][1], threed=threed,
                                            fm0=fm00)
                
                # need to do a guassian propogration of errors if the case we have 
                # errors and difference inversion 
                if err and self.param['reg_mode'] == 2 and errTyp!='global':
                    resError0 = df0['resError0'][indexes[i][0]].values 
                    resError = protocol['resError'].values 
                    protocol.loc[:,'resError'] = np.sqrt(resError**2 + resError0**2)
                
                if i == 0:
                    refdir = os.path.join(self.dirname, 'ref')
                    if os.path.exists(refdir) == False:
                        os.mkdir(refdir)
                    if 'mesh.dat' in os.listdir(self.dirname):
                        shutil.copy(os.path.join(self.dirname, 'mesh.dat'),
                                    os.path.join(refdir, 'mesh.dat'))
                    if 'mesh3d.dat' in os.listdir(self.dirname):
                        shutil.copy(os.path.join(self.dirname, 'mesh3d.dat'),
                                    os.path.join(refdir, 'mesh3d.dat'))
                    s.write2protocol(os.path.join(refdir, 'protocol.dat'), err=err, 
                                     errTot=errTot, threed=threed) # no subset for background, just use all
                else:
                    content = content + str(protocol.shape[0]) + '\n'
                    content = content + protocol.to_csv(sep='\t', header=False, 
                                                        index=False, lineterminator='\n')

            with open(os.path.join(self.dirname, 'protocol.dat'), 'w') as f:
                f.write(content)

        # for batch inversion -------------------
        elif self.iBatch is True:
            content = ''
            for i, s in enumerate(self.surveys):
                # resError and phaseError should already have been populated
                # handle the case when SOME survey were fitted but not all
                # then we use the bigSurvey default fit to fullfill them
                if err is True and errTyp == 'global':
                    if np.sum(np.isnan(s.df['resError'])) != 0: # there is some nan
                        print('Survey {:s} has no fitted error model, default to combined fit.'.format(s.name))
                        if self.bigSurvey.errorModel is None:
                            self.bigSurvey.fitErrorPwl() # default fit
                        s.df['resError'] = self.bigSurvey.errorModel(s.df)
                    if self.typ[0] == 'c' and np.sum(np.isnan(s.df['phaseError'])) != 0: # there is some nan
                        print('Survey {:s} has no fitted IP error model, default to combined fit.'.format(s.name))
                        if self.bigSurvey.phaseErrorModel is None:
                            self.bigSurvey.fitErrorPwlIP()
                        s.df['phaseError'] = self.bigSurvey.phaseErrorModel(s.df)
                    # if not it means that the 'resError' columns has already
                    # been populated when the files has been imported
                df = s.write2protocol(outputname='', err=err, ip=ipBool, errTot=errTot, threed=threed)
                content = content + str(len(df)) + '\n'
                content = content + df.to_csv(sep='\t', header=False, index=False)
            with open(os.path.join(self.dirname, 'protocol.dat'), 'w') as f:
                f.write(content)

        # for normal inversion (one survey) --------------------------
        else:
            self.surveys[0].write2protocol(os.path.join(self.dirname, 'protocol.dat'),
                        err=err, ip=ipBool, errTot=errTot, threed=threed)


    def runR2(self, dirname='', dump=None):
        """Run the executable in charge of the inversion.

        Parameters
        ----------
        dirname : str, optional
            Path of the directory where to run the inversion code.
        dump : function, optional
            Function to print the output of the invrsion code while running.
        """
        if dump is None:
            def dump(x):
                print(x, end='')
                
        # run R2.exe
        exeName = self.typ + '.exe'
        if dirname == '':
            dirname = self.dirname

        # get R2.exe path
        with cd(dirname):
            exePath = os.path.join(self.apiPath, 'exe', exeName)
    
            if OS == 'Windows':
                cmd = [exePath]
            elif OS == 'Darwin':
                # winetxt = 'wine'
                # if getMacOSVersion():
                #     winetxt = 'wine64'
                winetxt = whichWineMac()
                winePath = []
                wine_path = Popen(['which', winetxt], stdout=PIPE, stderr=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
                for stdout_line in iter(wine_path.stdout.readline, ''):
                    winePath.append(stdout_line)
                if winePath != []:
                    cmd = ['%s' % (winePath[0].strip('\n')), exePath]
                else:
                    cmd = [wPath, exePath]
            
            else:  # linux here
                # check if running on raspberrypi
                #if platform.machine() == 'aarch64':
                #    cmd = [exePath.replace('.exe', '_aarch64')]
                if platform.machine() in ['armv7l', 'aarch64']:
                    cmd = [exePath.replace('.exe', '_armv7l')]
                else:
                    cmd = ['wine',exePath]
    
            if OS == 'Windows':
                startupinfo = subprocess.STARTUPINFO()
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    
            def execute(cmd):
                if OS == 'Windows':
                    self.proc = subprocess.Popen(cmd, stdout=PIPE, stderr=PIPE, shell=False, universal_newlines=True, startupinfo=startupinfo)
                else:
                    self.proc = subprocess.Popen(cmd, stdout=PIPE, stderr=PIPE, shell=False, universal_newlines=True)
                for stdout_line in iter(self.proc.stdout.readline, ""):
                    yield stdout_line
                self.proc.stdout.close()
                return_code = self.proc.wait()
                if return_code:
                    print('error on return_code')
            for text in execute(cmd):
                dump(text)


    def runParallel(self, dirname=None, dump=None, iMoveElec=False,
                    ncores=None, rmDirTree=True):
        """Run several instances of R2 in parallel according to the number of
        cores available.

        Parameters
        ----------
        dirname : str, optional
            Path of the working directory.
        dump : function, optional
            Function to be passed to `R2.runR2()` for printing output during
            inversion.
        iMoveElec : bool, optional
            If `True` will move electrodes according to their position in each
            `Survey` object.
        ncores : int, optional
            Number or cores to use. If None, the maximum number of cores
            available will be used.
        rmDirTree: bool, optional
            Remove excess directories and files created during parallel.
            Default is True.
        """
        if dirname is None:
            dirname = self.dirname
            
        if dump is None:
            def dump(x):
                print(x, end='')

        if self.iTimeLapse is True and self.iBatch is False:
            surveys = self.surveys[1:]
        else:
            surveys = self.surveys

        # create R2.exe path
        exeName = self.typ + '.exe'
        exePath = os.path.join(self.apiPath, 'exe', exeName)


        # split the protocol.dat
        # in pandas >= 1.4.0 header=None with first row with one column (nb of meas)
        # causes ParseError. to fix it we first read the number of rows from line 2
        with open(os.path.join(self.dirname, 'protocol.dat'), 'r') as f:
            f.readline()  # first wow, we don't care
            nline = len(f.readline().split('\t'))
        dfall = pd.read_csv(os.path.join(self.dirname, 'protocol.dat'),
                            sep='\t', header=None, names=np.arange(nline))
        
        # the line where the last column is NaN is a line where a new dataset start
        idf = list(np.where(np.isnan(dfall[dfall.columns[-1]].values))[0])
        idf.append(len(dfall))
        dfs = [dfall.loc[idf[i]:idf[i+1]-1,:] for i in range(len(idf)-1)]

        # writing all protocol.dat
        files = []
        for s, df in zip(surveys, dfs):
            outputname = os.path.join(dirname, 'protocol_' + s.name + '.dat')
            files.append(outputname)
            df.to_csv(outputname, sep='\t', header=False, index=False)
            # header with line count already included

        # if iMoveElec is True, writing different R2.in
        if iMoveElec is True:
            dump('Electrodes position will be updated for each survey\n')
            for s in self.surveys:
                # print(s.name, '...', end='')
                elec = s.elec[['x','y','z']].values
                e_nodes = self.mesh.moveElecNodes(elec[:,0], elec[:,1], elec[:,2])
                self.param['node_elec'][1] = e_nodes + 1 # WE MUST ADD ONE due indexing differences between python and fortran
                if int(self.mesh.cell_type[0])==8 or int(self.mesh.cell_type[0])==9:#elements are quads
                    colx = self.mesh.quadMeshNp() # so find x column indexes instead. Wont support change in electrode elevation
                    self.param['node_elec'][1] = colx
                self.param['inverse_type'] = 1 # regularise against a background model
                #self.param['reg_mode'] = 1
                write2in(self.param, self.dirname, self.typ)
                r2file = os.path.join(self.dirname, self.typ + '.in')
                shutil.move(r2file, r2file.replace('.in', '_' + s.name + '.in'))
                # print('done')

        # create workers directory
        ncoresAvailable = sysinfo['cpuCount']
        if ncores is None: # and self.ncores is None:
            ncores = sysinfo['physicalCpuCount']
        else:
            if ncores > ncoresAvailable:
                raise ValueError('Number of cores larger than available')
        dump('Using %i logical processors'%ncores)


        def prepare(wd, fname):
            # copying usefull files from the main directory
            toMove = ['mesh.dat', 'mesh3d.dat','R2.in','cR2.in',
                      'R3t.in', 'cR3t.in', 'res0.dat','resistivity.dat',
                      'Start_res.dat']
            for f in toMove:
                file = os.path.join(dirname, f)
                if os.path.exists(file):
                    shutil.copy(file, os.path.join(wd, f))

            # copy the protocol.dat
            shutil.copy(fname, os.path.join(wd, 'protocol.dat'))
            name = os.path.basename(fname).replace('.dat', '').replace('protocol_','')
            if iMoveElec is True:
                r2inFile = os.path.join(os.path.dirname(fname),
                                        self.typ + '_' + name + '.in')
                shutil.copy(r2inFile, os.path.join(wd, self.typ + '.in'))

        if OS == 'Windows':
            cmd = ['powershell', '-Command', '&', '"'+exePath+'"']
        elif OS == 'Darwin':
            # winetxt = 'wine'
            # if getMacOSVersion():
            #     winetxt = 'wine64'
            winetxt = whichWineMac()
            winePath = []
            wine_path = Popen(['which', winetxt], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
            for stdout_line in iter(wine_path.stdout.readline, ''):
                winePath.append(stdout_line)
            if winePath != []:
                cmd = ['%s' % (winePath[0].strip('\n')), exePath]
            else:
                cmd = [wPath, exePath]
        else:
            # check if running on raspberrypi
            if platform.machine == 'aarch64':
                cmd = [exePath.replace('.exe', '_aarch64')]
            if platform.machine == 'arm':
                cmd = [exePath.replace('.exe', '_arm')]
            else:
                cmd = ['wine',exePath]

        if OS == 'Windows':
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW

        def retrieve(wd, fname):
            # moving inversion results back
            name = os.path.basename(fname).replace('.dat', '').replace('protocol_','')
            originalDir = self.dirname
            toMove = ['f001_res.dat', 'f001_res.vtk', 'f001_err.dat',
                      'f001_sen.dat', 'f001_diffres.dat',
                      'f001.dat', 'f001.sen', 'f001.err', 'f001.vtk'] # all 3D stuff
            for f in toMove:
                if os.path.exists(os.path.join(wd, f)):
                    shutil.move(os.path.join(wd, f),
                                os.path.join(originalDir, f.replace('f001', name)))
            shutil.move(os.path.join(wd, self.typ + '.out'),
                        os.path.join(originalDir, self.typ + '_' + name + '.out'))
            shutil.move(os.path.join(wd, 'electrodes.dat'),
                        os.path.join(originalDir, 'electrodes_' + name + '.dat'))
            shutil.move(os.path.join(wd, 'electrodes.vtk'),
                        os.path.join(originalDir, 'electrodes_' + name + '.vtk'))

        # create all the working directories
        wds = []
        for i, f in enumerate(files):
            wd = os.path.join(self.dirname, str(i+1))
            if os.path.exists(wd):
                shutil.rmtree(wd)
            os.mkdir(wd)
            prepare(wd, f)
            wds.append(wd)
        wds2 = wds.copy()

        # run them all in parallel as child processes
        def dumpOutput(out):
            for line in iter(out.readline, ''):
                dump(line.rstrip() + '\n')
            out.close()

        # create essential attribute
        self.irunParallel2 = True
        self.procs = []
        ts = []

        # kill management
        self.proc = ProcsManagement(self)

        # run in // (http://code.activestate.com/recipes/577376-simple-way-to-execute-multiple-process-in-parallel/)
        # In an infinite loop, will run an number of process (according to the number of cores)
        # the loop will check when they finish and start new ones.
        def done(p):
            return p.poll() is not None

        c = 0
        dump('\r{:.0f}/{:.0f} inversions completed'.format(c, len(wds2)))
        while self.irunParallel2:
            while wds and len(self.procs) < ncores:
                wd = wds.pop()
                # NOTE: we need to PIPE stdout and sterr to avoid being
                # flooded by R2 output. 'All ok' is only printed in stderr
                if OS == 'Windows':
                    p = Popen(cmd, cwd=wd, shell=False, universal_newlines=True, startupinfo=startupinfo, stdout=PIPE, stderr=PIPE)
                else:
                    p = Popen(cmd, cwd=wd, shell=False, universal_newlines=True, stdout=PIPE, stderr=PIPE)
                self.procs.append(p)
                # t = Thread(target=dumpOutput, args=(p.stdout,))
                # t.daemon = True # thread dies with the program
                # t.start()
                # ts.append(t)

            for p in self.procs:
                if done(p):
                    #print('------done!!', p)
                    self.procs.remove(p)
                    c = c+1
                    dump('\r{:.0f}/{:.0f} inversions completed'.format(c, len(wds2)))

            if not self.procs and not wds:
                dump('\n')
                break
            else:
                time.sleep(0.05)

        for wd, f in zip(wds2, files):
            try:
                retrieve(wd, f)
            except Exception as e:
                print('Error retrieving for ', wd, ':', e)
                pass


        # get the files as it was a sequential inversion
        toRename = ['_res.dat', '_res.vtk', '_err.dat', '_sen.dat', '_diffres.dat']
        r2outText = ''
        for i, s in enumerate(surveys):
            for ext in toRename:
                originalFile = os.path.join(dirname,  s.name + ext)
                newFile = os.path.join(dirname, 'f' + str(i+1).zfill(3) + ext)
                if os.path.exists(originalFile):
                    shutil.move(originalFile, newFile)
            r2outFile = os.path.join(dirname, self.typ + '_' + s.name + '.out')
            with open(r2outFile, 'r') as f:
                r2outText = r2outText + f.read()
            os.remove(r2outFile)
        with open(os.path.join(dirname, self.typ + '.out'), 'w') as f:
            f.write(r2outText)
        
        # remove electrodes files if iMoveElec is False
        if iMoveElec is False:
            for i, s in enumerate(surveys[:-1]):
                os.remove(os.path.join(dirname, 'electrodes_' + s.name + '.vtk'))
                os.remove(os.path.join(dirname, 'electrodes_' + s.name + '.dat'))
            shutil.move(os.path.join(dirname, 'electrodes_' + surveys[-1].name + '.vtk'),
                        os.path.join(dirname, 'electrodes.vtk'))
            shutil.move(os.path.join(dirname, 'electrodes_' + surveys[-1].name + '.dat'),
                        os.path.join(dirname, 'electrodes.dat'))      

        # delete the dirs and the files
        if rmDirTree:
            [shutil.rmtree(d) for d in wds2]
            [os.remove(f) for f in files]

        print('----------- END OF INVERSION IN // ----------')



    def invert(self, param={}, iplot=False, dump=None, err=None, modErr=False,
               parallel=False, iMoveElec=False, ncores=None,
               rmDirTree=True, modelDOI=False, errResCol=None, errIPCol=None, 
               pseudoStart = True):
        """Invert the data, first generate R2.in file, then run
        inversion using appropriate wrapper, then return results.  

        Parameters
        ----------
        param : dict, optional
            Dictionary of parameters for inversion. Will be passed to
            `R2.write2in()`.
        iplot : bool, optional
            If `True`, will plot the results of the inversion using
            `R2.showResults()`.
        dump : function, optinal
            Function to print the output of the inversion. To be passed to
            `R2.runR2()`.
        err : bool, optional
            If 'True' reciprocal error model will be used to weight the 
            inversion. Can be used to overide self.err. Default is None (in 
            which case self.err sets the inclusion of error weighting). 
        modErr : bool, optional
            If `True`, the model error will be compute and added before
            inversion.
        parallel : bool, optional
            If `True`, batch and time-lapse survey will be inverted in //. No
            output will be display during inversion.
        iMoveElec : bool, optional
            If `True`, then different electrode location will be used for
            the different surveys. Electrodes location are specified in the
            `Survey` object. Only for parallel inversion for now.
        ncores : int, optional
            If `parallel==True` then ncores is the number of cores to use (by
            default all the cores available are used).
        rmDirTree : bool, optional
            Remove excess directories and files created during parallel inversion
        modelDOI : bool, optional
            If `True`, the Depth of Investigation will be model by reinverting
            the data on with an initial res0 different of an order of magnitude.
            Note that this option is only available for *single* survey.
        errResCol : str, optional
            If no reciprocal data is present but an error column is supplied,
            this error column can be used in the inversion. Note that the column
            must contain a value for each quadrupole.
        errIPCol : str, optioanl
            If no reciprocal data, a user supplied column can be read as IP
            error.
        pseudoStart : bool, optional 
            Flag to start the inversion using using resisitivities from the pseudo section.
            Default is True, in which case the starting resistivities will be estimated
            from the psuedo section values. 
        """
        if dump is None:
            def dump(x):
                print(x, end='')
                
        # check inversion is doable first with the amount of total amount of RAM available to the machine 
        ncores = self._checkMemory(parallel,ncores)
                
        # clean meshResults list
        self.meshResults = []
        
        # create mesh if not already done
        if 'mesh' not in self.param:
            self.createMesh()
            
        # clean previous iterations
        for f in os.listdir(self.dirname):
            if f[:3] == 'f00':
                os.remove(os.path.join(self.dirname, f))
                
        # deal with potential user supply error (res)
        if errResCol is not None:
            for survey in self.surveys:
                if errResCol in survey.df.columns:
                    if survey.df[errResCol].notnull().sum() == survey.df.shape[0]:
                        survey.df['resError'] = survey.df[errResCol]
                        err = True  # we assume that if the user supplied a column it's to use it as error
                    else:
                        dump('ERROR: ' + errResCol + ' does not contain a value for all quadrupoles.')
                else:
                    dump('ERROR: ' + errResCol + ' is not found in the dataframe.')

        # deal with potential user supply error (ip)
        if errIPCol is not None:
            for survey in self.surveys:
                if errIPCol in survey.df.columns:
                    if survey.df[errIPCol].notnull().sum() == survey.df.shape[0]:
                        survey.df['phaseError'] = survey.df[errIPCol]
                        err = True  # we assume that if the user supplied a column it's to use it as error
                    else:
                        dump('ERROR: ' + errIPCol + ' does not contain a value for all quadrupoles.')
                else:
                    dump('ERROR: ' + errIPCol + ' is not found in the dataframe.')
            
        # run Oldenburg and Li DOI estimation
        if modelDOI is True:
            sensScaled = self.modelDOI(dump=dump)
            
        if pseudoStart: 
            dump('Looking up starting resistivities from pseudo section')
            self.mesh.df.loc[:,'res0'] = self.res0fromPseudo() # res0 assigned

        # compute modelling error if selected
        if modErr is True and self.fwdErrModel is False: #check no error model exists
            # dump('Computing error model... ')
            self.computeModelError(dump=dump)
            dump('done\n')
            errTot = True
        elif modErr is True and self.fwdErrModel:
            # avoid computing error model again if it has already been run.
            errTot = True
        else:
            errTot = False
            
        # catch case where errTot requested but no reciprocal error model is present 
        if errTot: 
            s = self.surveys[0]
            if np.sum(np.isnan(s.df['resError'])) == len(s.df): 
                if 'a_wgt' in self.param.keys(): 
                    if self.param['a_wgt'] > 0: 
                        self.estimateError(
                            a_wgt = self.param['a_wgt'], 
                            b_wgt = self.param['b_wgt'])

        # write configuration file
        dump('Writing .in file and protocol.dat... ')
        self.write2in(param=param, err=err) # R2.in
        self.write2protocol(err=err, errTot=errTot) # protocol.dat
        dump('done\n')

        # runs inversion
        if self.iTimeLapse == True and self.referenceMdl==False:
            dump('------------ INVERTING REFERENCE SURVEY ---------------\n')
            refdir = os.path.join(self.dirname, 'ref')
            shutil.move(os.path.join(self.dirname,'res0.dat'),
                        os.path.join(refdir, 'res0.dat'))
            self.write2in(param=param,err=err)
            self.runR2(refdir, dump=dump) # this line actually runs R2
            shutil.copy(os.path.join(refdir, 'f001_res.dat'),
                        os.path.join(self.dirname, 'Start_res.dat'))
            res0 = np.loadtxt(os.path.join(refdir, 'f001_res.dat'))[:,3] # reference resistivity  
            
            if ((self.typ == 'R3t') | (self.typ == 'cR3t')) & (self.param['reg_mode'] == 2):
                dump('----------------- Computing d-d0+f(m0) ---------------\n')
                # as per v3.2 of R3t we need to compute MANUALLY d-d0+f(m0)
                # this is done automatically in R2 and cR2
                # unfortunatly R3t does not output the transfer resistances of the final model 
                # (only apparent resistivities) , therefore forward modelling the baseline 
                # resistivities seems like the best option 
                self.sequence = self.surveys[0].df[['a','b','m','n']].values
                surveysBackup = self.surveys.copy()
                # res0Backup = self.mesh.df['res0'].values.copy() # might want to reset res0 
                self.mesh.df.loc[:,'res0'] = res0 # set res 0 to baseline inversion 
                iForwardBackup = self.iForward
                self.forward()
                self.iForward = iForwardBackup
                fm0 = self.surveys[0].df['resist'].values.copy()
                self.sequence = None
                self.surveys = surveysBackup
                self.write2protocol(err=err, errTot=False, fm0=fm0) # rewrite them with d-d0+f(m0)
        elif self.iTimeLapse == True and self.referenceMdl==True:
            print('Note: Skipping reference inversion, as reference model has already been assigned')

        dump('\n--------------------- MAIN INVERSION ------------------\n')
        if parallel is True and (self.iTimeLapse is True or self.iBatch is True):
            self.runParallel(dump=dump, iMoveElec=iMoveElec, ncores=ncores, rmDirTree=rmDirTree)
        else:
            self.runR2(dump=dump)
            
        # extract inversion errors
        try: # this is in the case getInvError() is called after the file .err is
            # created by R2 but before it is populated (when killing the run)
            self.getInvError()
            self.getResults()
            if modelDOI is True:
                for m in self.meshResults:
                    m.df['doiSens'] = sensScaled
            
            # read final R2.out
            with open(os.path.join(self.dirname, self.typ + '.out'),'r') as f:
                self.invLog += f.read() + '\n'
                
        except Exception as e:
            print('Could not retrieve files maybe inversion failed')
            print('Error: ', e)
            return

        if iplot is True:
            if self.iForward:
                self.showResults(index=1)
            else:
                self.showResults()
                
                


    def modelDOI(self, dump=None):
        """Will rerun the inversion with a background constrain (alpha_s) with
        the normal background and then a background 10 times more resistive.
        From the two different inversion a senstivity limit will be computed.
        """
        if dump is None:
            def dump(x):
                print(x, end='')
                
        # we store the value of the iForward flag is this flag is True
        # res0 is automatically ovewritten when doing the inversion
        iForward = self.iForward  # store the value
        self.iForward = False  # set to false so write2in won't overwrite us
        print('All non fixed parameters reset to 100 Ohm.m and 0 mrad, '
              'as the survey to be inverted is from a forward model.')
                
        # backup for normal inversion (0 : original, 1 : normal background, 2: background *10)
        ifixed = self.mesh.df['param'] == 0
        res0 = np.array(self.mesh.df['res0'])
        res0f = res0.copy()
        res0f[~ifixed] = 100
        self.mesh.df['res0'] = list(res0f)
        param0 = self.param.copy()
        self.param['reg_mode'] = 1 # we need constrain to background
        typ0 = self.typ
        if self.typ[0] == 'c':
            self.typ = self.typ[1:]
        iTimeLapse0 = self.iTimeLapse
        self.iTimeLapse = False
        surveys0 = self.surveys.copy()
        self.surveys = [surveys0[0]] # just use first survey
        self.write2in()
        self.write2protocol()
        
        # build the cropping polygon
        if self.param['num_xz_poly'] != 0:
            path = mpath.Path(self.param['xz_poly_table'])
            iselect = path.contains_points(np.c_[self.mesh.elmCentre[:,0], self.mesh.elmCentre[:,2]])
        else:
            iselect = np.ones(len(self.mesh.elmCentre[:,0]), dtype=bool)
            
        # clean function
        def cleandir():
            dirname = self.dirname
            os.remove(os.path.join(dirname, 'res0.dat'))
            for f in os.listdir(dirname):
                if f[:3] == 'f00':
                    os.remove(os.path.join(dirname, f))
        
        # run first background constrained inversion
        dump('===== modelDOI: Running background constrained inversion with initial resistivity =====\n')
        res1 = res0
        self.mesh.df['res0b'] = res1
        self.mesh.writeAttr('res0b', os.path.join(self.dirname,'res0.dat'))
        self.runR2(dump=dump) # re-run inversion
        self.getResults()
        mesh1 = self.meshResults[0]
        cleandir()
        
        # run second background constrained inversion
        dump('===== modelDOI: Running background constrained inversion with initial resistivity * 10 =====\n')
        res2 = res0 * 10
        self.mesh.df['res0b'] = list(res2)
        self.mesh.writeAttr('res0b', os.path.join(self.dirname,'res0.dat'))
        self.runR2(dump=dump) # re-run inversion
        self.getResults()
        mesh2 = self.meshResults[0]
        cleandir()
        os.remove(os.path.join(self.dirname, 'R2.in'))
        
        # sensitivity = difference between final inversion / difference init values
        res_names = np.array(['Resistivity','Resistivity(Ohm-m)','Resistivity(ohm.m)'])
        res_name = res_names[np.in1d(res_names, list(self.meshResults[0].df.keys()))][0]
        invValues1 = np.array(mesh1.df[res_name])
        invValues2 = np.array(mesh2.df[res_name])
        sens = (invValues1 - invValues2)/(res1[iselect]-res2[iselect])
        sensScaled = np.abs(sens)
        self.doiComputed = True
        
        # restore
        self.meshResults = []
        self.param = param0
        self.typ = typ0
        self.surveys = surveys0
        self.iTimeLapse = iTimeLapse0
        self.mesh.df['res0'] = list(res0)
        self.iForward = iForward # restore value
        # .in and protocol will be written again in R2.invert()
            
        return sensScaled
        
    def _clipContour(self, ax, cax, cropMaxDepth=False, clipCorners=False):
        """Clip contours using mesh bound and surface if available.
        
        Parameters
        ----------
        ax : matplotlib.Axes
            Axis.
        cax : matplotlib.Artist
            Artist (Collection) returned.
        cropMaxDepth : bool, optional
            If 'True', area below fmd will be cropped out.
        clipCorners : bool, optional
            If 'True', triangles from bottom corners will be cropped (only if the whole mesh is not shown).
        """
        
        def patcher(verts):
            # cliping using a patch (https://stackoverflow.com/questions/25688573/matplotlib-set-clip-path-requires-patch-to-be-plotted)
            poly_codes = [mpath.Path.MOVETO] + (len(verts) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
            path = mpath.Path(verts, poly_codes)
            patch = mpatches.PathPatch(path, facecolor='none', edgecolor='none')
            ax.add_patch(patch) # need to add so it knows the transform
            cax.set_clip_path(patch)
        
        # mask outer region
        node_x = self.mesh.node[:,0]
        node_z = self.mesh.node[:,2]
        xmin = np.min(node_x)
        xmax = np.max(node_x)
        zmin = np.min(node_z)
        zmax = np.max(node_z)
        
        (xsurf, zsurf) = self.mesh.extractSurface()
        if cropMaxDepth and self.fmd is not None:
            xfmd, zfmd = xsurf[::-1], zsurf[::-1] - self.fmd
            verts = np.c_[np.r_[xmin, xmin, xsurf, xmax, xmax, xfmd, xmin],
                          np.r_[zmin, zmax, zsurf, zmax, zmin, zfmd, zmin]]
        else:
            verts = np.c_[np.r_[xmin, xmin, xsurf, xmax, xmax, xmin],
                          np.r_[zmin, zmax, zsurf, zmax, zmin, zmin]]
        patcher(verts)

        if clipCorners and self.param['num_xz_poly'] != 0: # not clipping the corners of a mesh outside of the survey area!
            elec_x = self.mesh.elec[:,0]
            elec_z = self.mesh.elec[:,2]
            elec_xmin = np.min(elec_x)
            elec_xmax = np.max(elec_x)
            if cropMaxDepth and self.fmd is not None:
                zminl = elec_z[elec_x.argmin()] - self.fmd
                zminr = elec_z[elec_x.argmax()] - self.fmd
            elif cropMaxDepth is False and self.fmd is not None:
                zminl = zminr = np.min(elec_z) - self.fmd
            else:
                zminl = zminr = zmin
            zmaxl = elec_z[elec_x.argmin()]
            zmaxr = elec_z[elec_x.argmax()]
            ll = np.abs(zmaxl - zminl)
            lr = np.abs(zmaxr - zminr)
            lx = np.abs(elec_xmin - elec_xmax)
            if ll >= lx/4:
                ll = lx/4
            if lr >= lx/4:
                lr = lx/4

            # surf bound to elecs
            elec_surf = np.c_[xsurf, zsurf][(np.abs(xsurf - elec_xmin)).argmin():
                                            (np.abs(xsurf - elec_xmax)).argmin(),:]
            idxl = (np.abs(elec_surf[:,0] - ll)).argmin()
            idxr = (np.abs(elec_surf[:,0] - np.abs(elec_xmax - lr))).argmin() + 1       
            xtrapbot = elec_surf[idxl:idxr,0]
            if cropMaxDepth and self.fmd is not None:
                ztrapbot = elec_surf[idxl:idxr,1] - self.fmd
            else:
                ztrapbot = np.ones_like(xtrapbot) * zminl
                
            self.trapeziod = np.c_[np.r_[elec_xmin, xsurf, elec_xmax, xtrapbot[::-1], elec_xmin],
                                   np.r_[zmaxl, zsurf, zmaxr, ztrapbot[::-1], zmaxl]]
            patcher(self.trapeziod)
        
        else:
            self.trapeziod = None # make sure trapeziod mask is clear
     

    def showResults(self, index=0, ax=None, edge_color='none', attr='',
                    sens=True, color_map='viridis', zlim=None, clabel=None,
                    doi=False, doiSens=False, contour=False, cropMaxDepth=True,
                    clipContour=True, clipCorners=False, use_pyvista=True, background_color=(0.8,0.8,0.8),
                    pvslices=([],[],[]), pvspline=None, pvthreshold=None, pvgrid=True,
                    pvcontour=[], pvdelaunay3d=False, volume=None, **kwargs):
        """Show the inverteds section.

        Parameters
        ----------
        index : int, optional
            Index of the inverted section (mainly in the case of time-lapse
            inversion). If index == -1, then all 2D survey will be plotted
            on a 3D grid.
        ax : matplotlib axis, optional
            If specified, the inverted graph will be plotted agains `ax`.
        edge_color : str, optional
            Color of the edges of the mesh.
        attr : str, optional
            Name of the attribute to be plotted.
        sens : bool, optional
            If `True` and if sensitivity is available, it will be plotted as
            a white transparent shade on top of the inverted section.
        color_map : str, optional
            Name of the colormap to be used.
        clabel : str, optional
            Label of the colorbar (by default the label is the value of `attr`).
        doi : bool, optional
            If True, it will draw a dotted red line corresponding to 0.02 from the
            Oldenburg and Li method. Note that `R2.modeDOI()` needs to be run
            for that.
        doiSens : bool, optional
            If True, it will draw a dashed line corresponding to 0.001 of the maximum
            of the log10 sensitivity.
        contour : bool, optional
            If True, contours will be plotted.
        cropMaxDepth : bool, optional
            If True, the mesh will be clipped with at a depth following the surface.
            If False, the mesh will be clipped at the maximum depth available.
            This doesn't have any effect if clipContour is False.
        clipContour : bool, optional
            If True, the contour of the area of interest will be clipped (default).
        clipCorners : bool, optional
            If 'True', triangles from bottom corners will be cropped (only if the whole mesh is not shown).
        use_pyvista : bool, optional
            (3D only) Use visual toolkit backend for displaying 3D mesh, note that pyvista
            must be installed for this to work. 
        background_color : tuple, optional 
            (3D only) Background color assigned to pyvista plotter object when created. Not yet
            supported for matplotlib axis handles. 
        pvslices : tuple of list of float, optional
            (3D only) Determine the X, Y, Z slices. e.g.: ([3], [], [-3, -4]) will add
            a slice normal to X in 3 and two slices normal to Z in -3 and -4.
        pvspline : 'elec' or numpy array, optional
            (3D only) If 'elec' mesh will be sliced along the electrodes path otherwise 
            an array of X, Y, Z of points on a path to slice the mesh along that path is needed.
        pvthreshold : list of two floats, optional
            (3D only) Keep values between pvthreshold[0] and pvthreshold[1].
        pvgrid : bool, optional
            (3D only) Show grid or not.
        pvcontour: list of float, optional
            (3D only) Values of the isosurface to be plotted.
        pvdelaunay3d : bool, optional
            If `True` a "Delaunay 3D" triangulation filter will be applied on the mesh.
        volume : float, optional
            If not 'None' then volume of selected region of the mesh will be printed onto the pyvista plot.
        """
        if len(self.meshResults) == 0:
            self.getResults()
        if (attr == '') & (self.typ[0] != 'c'):
            attr = 'Resistivity(log10)'
        if (attr == '') & (self.typ[0] == 'c'):
            attr = 'Sigma_real(log10)'
        keys = list(self.meshResults[index].df.keys())
        if attr not in keys:
            attr = keys[3]
            print('Attribute not found, revert to {:s}'.format(attr))
            
        #check for clipping 
        if 'typ' in self.meshParams.keys() and 'clipping' not in kwargs:
            for a in ['tank', 'cylinder', 'prism']:
                if self.meshParams['typ'] == a:
                    kwargs['clipping'] = False 
                    
        elec = self.elec[['x','y']].values
        elec = elec[~self.elec['remote'].values,:]
         
        # best to use x/y limits provided by project class becuase it knows if electrodes are remote 
        if 'xlim' not in kwargs.keys():
            # scale up electrode coordinates by 20% to give some padding around where the electrodes are 
            kwargs['xlim'] = findminmax(elec[:,0]) 
            # can use surface topography too if that has been provided 
            if len(self.topo) > 4:
                kwargs['xlim'] = findminmax(self.topo.x,1)
        if 'ylim' not in kwargs.keys() and self.typ[-1] == 't':
            kwargs['ylim'] = findminmax(elec[:,1])
            if len(self.topo) > 4:
                kwargs['ylim'] = findminmax(self.topo.y,1)
            
        if len(self.meshResults) > 0:
            mesh = self.meshResults[index]
            if self.typ[-1] == '2' and index != -1: # 2D case
                if self.pseudo3DSurvey is not None and self.projs != []: # we have pseudo 3D survey
                    self.param['num_xz_poly'] = self.projs[index].param['num_xz_poly']
                    self.mesh = mesh # should update this based on current mesh to get right limits
                    self.zlim = self.projs[index].zlim
                    self.fmd = self.projs[index].fmd
                if zlim is None:
                    zlim = self.zlim
                mesh.show(ax=ax, edge_color=edge_color, darkMode=self.darkMode,
                            attr=attr, sens=sens, color_map=color_map,
                            zlim=zlim, clabel=clabel, contour=contour, **kwargs)
                if doi is True: # DOI based on Oldenburg and Li
                    if self.doiComputed is True: 
                        z = np.array(mesh.df['doiSens'])
                        levels = [0.2]
                        linestyle = ':'
                    else:
                        raise ValueError('Rerun the inversion with `modelDOI=True` first or use `doiSens`.')
                if doiSens is True: # DOI based on log10(sensitivity)
                    if 'Sensitivity(log10)' in mesh.df.keys():
                        z = np.array(mesh.df['Sensitivity(log10)'])
                        levels=[np.log10(0.001*(10**np.nanmax(z)))]
                        linestyle = '--'
                    else:
                        doiSens = False
                if (clipContour) & (self.topo.shape[0] == 0) & (all(self.elec['buried'])):
                    # it's a whole space mesh, clipContour is not needed for that but contour can be drawn
                    cropMaxDepth = False
                if doi is True or doiSens is True:
                    xc, yc = mesh.elmCentre[:,0], mesh.elmCentre[:,2]
                    triang = tri.Triangulation(xc, yc)
                    cont = mesh.ax.tricontour(triang, z, levels=levels, colors='k', linestyles=linestyle)
                    if clipContour:
                        self._clipContour(mesh.ax, cont, clipCorners=clipCorners)
                if clipContour:
                    self._clipContour(mesh.ax, mesh.cax, cropMaxDepth=cropMaxDepth, clipCorners=clipCorners)
            elif self.typ[-1] == '2' and index == -1: # 3D grid of 2D surveys (pseudo 3D)
                self.showPseudo3DResults(ax=ax, edge_color=edge_color,
                    attr=attr, color_map=color_map, clabel=clabel, returnMesh=True,
                    use_pyvista=use_pyvista, background_color=background_color,
                    pvslices=pvslices, pvspline=pvspline, pvthreshold=pvthreshold, pvgrid=pvgrid,
                    pvcontour=pvcontour, cropMaxDepth=cropMaxDepth, clipCorners=clipCorners, **kwargs)
            else: # 3D case
                if zlim is None:
                    zlim = self.zlim
                    # zlim = [np.min(mesh.node[:,2]), np.max(mesh.node[:,2])]
                if cropMaxDepth and self.fmd is not None and zlim is None:
                    zlim[0] = self.elec[self.elec['remote']]['z'].min() - self.fmd # TODO not sure about that
                mesh.show(ax=ax, edge_color=edge_color,
                        attr=attr, color_map=color_map, clabel=clabel,
                        zlim=zlim, use_pyvista=use_pyvista, background_color=background_color,
                        pvslices=pvslices, pvspline=pvspline, pvthreshold=pvthreshold, pvgrid=pvgrid,
                        pvcontour=pvcontour, pvdelaunay3d=pvdelaunay3d, darkMode=self.darkMode, volume=volume, **kwargs)
                
        else:
            raise ValueError('len(R2.meshResults) == 0, no inversion results parsed.')



    def getResults(self, dirname=None):
        """Collect inverted results after running the inversion and adding
        them to `R2.meshResults` list.
        
        Parameters
        ----------
        dirname : str, optional
            If specified, dirname will be used as the working directory (this
            is needed for R2.loadResults()). Default is self.dirname.
        """
        if dirname is None:
            dirname = self.dirname
        self.getRMS(dirname)
        self.computeRMSP(dirname)
        idone = 0
        ifailed = 0
        self.meshResults = [] # make sure we empty the list first
        if self.iTimeLapse: # grab reference mesh 
            fname = os.path.join(dirname, 'ref', 'f001_res.vtk')
            mesh0 = mt.vtk_import(fname, order_nodes=False)
            mesh0.mesh_title = self.surveys[0].name
            elec = self.surveys[0].elec.copy()
            mesh0.setElec(elec['x'].values, elec['y'].values, elec['z'].values)
            mesh0.iremote = elec['remote'].values
            self.meshResults.append(mesh0)
            idone += 1
            
        if self.iForward is True:
            initMesh = mt.vtk_import(os.path.join(dirname, 'fwd','forward_model.vtk'), order_nodes=False)
            elec_x = self.elec['x'].values
            elec_y = self.elec['y'].values
            elec_z = self.elec['z'].values
            initMesh.setElec(elec_x, elec_y, elec_z)
            initMesh.mesh_title = 'Forward Model'
            self.meshResults.append(initMesh)

        for i in range(len(self.surveys)):
            if self.iTimeLapse is True:
                j = i + 1
            else:
                j = i
            fname = os.path.join(dirname, 'f' + str(i+1).zfill(3) + '_res.vtk')
            if os.path.exists(fname):
                try:
                    mesh = mt.vtk_import(fname, order_nodes=False)
                    mesh.mesh_title = self.surveys[j].name
                    elec = self.surveys[j].elec.copy()
                    mesh.setElec(elec['x'].values, elec['y'].values, elec['z'].values)
                    mesh.iremote = elec['remote'].values
                    self.meshResults.append(mesh) # this will be very memory intensive to put all meshes into a list for long time lapse surveys
                    #TODO : Rethink storage of timelapse results 
                    idone += 1
                except Exception as e:
                    ifailed += 1
                    # print(e)
                    # if inversion fails in time-lapse it's that the initial
                    # model is good enough to explain the data (a_wgt/b_wgt
                    # error too low) so we can replace it by the initial model
                    if self.iTimeLapse:
                        self.meshResults.append(mesh0) 
                        self.meshResults[-1].mesh_title = self.surveys[j].name
                print('\r{:d}/{:d} results parsed ({:d} ok; {:d} failed)'.format(
                    j+1, len(self.surveys), idone, ifailed), end='')
            else:
                pass
                #break
                
        self.pinfo['Number of Surveys'] = len(self.meshResults)
        print('')

        # compute conductivity in mS/m
        res_names = np.array(['Resistivity','Resistivity(Ohm-m)','Resistivity(ohm.m)', 'Magnitude(ohm.m)'])
        for mesh in self.meshResults:
            res_name = res_names[np.in1d(res_names, list(mesh.df.keys()))][0]
            mesh.df['Conductivity(mS/m)'] = 1000/np.array(mesh.df[res_name])
        if self.typ[0] == 'c' and self.surveys[0].kFactor != 1: # if kFactor is 1 then probably phase is provided and we shouldn't estimate chargeability
            for mesh in self.meshResults:
                mesh.df['Chargeability(mV/V)'] = np.array(mesh.df['Phase(mrad)'])/-self.surveys[0].kFactor
        
        # compute difference in percent for time lapse surveys 
        if (self.iTimeLapse is True):# and (self.param['reg_mode'] == 1):
            # even with reg_mode == 2 when the inversion converged by overshooting
            # it won't output 'difference(percent)' attribute, so let's compute for all TL
            try:
                self.postProcTl() 
            except Exception as e:
                print('failed to compute difference: ', e)
                pass
            
            
    def getRMS(self,dirname=None):
        if dirname is None:
            dirname = self.dirname 
        fpath = os.path.join(dirname, self.typ+'.out')
        if not os.path.exists(fpath):
            return 
        fh = open(fpath,'r')
        lines = fh.readlines() 
        rms = 0 
        iterations = 0 
        for line in lines: 
            if 'Final RMS Misfit' in line: 
                rms = float(line.split()[-1])
            if 'Iteration' in line:
                iterations += 1 
        if rms == 0:
            self.pinfo['Convergence of inverse solution'] = False 
        else:
            self.pinfo['Convergence of inverse solution'] = True 
        self.pinfo['Number of iterations'] = iterations 
        self.pinfo['RMS of inverse solution'] = rms 
        return rms 
    
    
    def computeRMSP(self, dirname=None):
        """
        Compute RMS as a percentage value (estimate), this is to mimic the behaviour of
        something like ResInv. 
        
        Parameters
        ----------
        dirname : str, optional
            Directory of working inversion. The default is the project working
            directory. 

        Returns
        -------
        rmsp: float 
            RMS of all surveys as a percentage 

        """
        if dirname is None:
            dirname = self.dirname 
        
        # need to find all instances of err.dat files 
        errfiles = []
        if 'ref' in os.listdir(dirname):
            fpath = os.path.join(dirname,'ref','f001_err.dat')
            if os.path.exists(fpath):
                errfiles.append(fpath)
        for f in os.listdir(dirname):
            if f.endswith('err.dat'):
                fpath = os.path.join(dirname, f)
                errfiles.append(fpath)
        
        iobs = 5 # observation column 
        ical = 6 # calculated measurement column 
        if '3' in self.typ: 
            iobs += 4 
            ical += 4 
        
        errdf = np.genfromtxt(errfiles[0], skip_header=1)
        _ = errfiles.pop(0)
        for f in errfiles:
            _df = np.genfromtxt(f, skip_header=1)
            errdf = np.vstack([errdf, _df])
        
        robs = errdf[:, iobs]
        rcal = errdf[:, ical]
        N = len(errdf)
        
        rmsp_sq = np.sum((100*((robs-rcal)/robs))**2)/N
        rmsp = np.sqrt(rmsp_sq)
        
        self.pinfo['RMS error as percentage estimate'] = float(rmsp)  
        return rmsp 
        
        
    def computeVol(self, attr='Resistivity(ohm.m)', vmin=None , vmax=None, index=0):
        """Given a 3D dataset, calculates volume based on mesh type.
        Note: The majority of this function's code is also seen in the meshTools' cellArea function as well. But that calculates all the mesh not the mesh results. 
            
        Parameters
        ----------
        attr: str, optional 
            Attribute displayed in mesh results where apropriate cell vertices are
            assigned a location. These locations are utilized for volume calculations
        vmin: float, optional
            User assigned minimum threshold value of attribute
        vmax: float, optional
            User asigned maximum threshold value of attribute
        index: int, optional
            meshResults index. In case of multi-file inversion/modeling (e.g., timelapse)
            
        Returns
        -------
        volume: float
            The volume of selected mesh region in cubic meters
        """
        mesh = self.meshResults[index].copy()
        con_mat = mesh.connection
        elm_vol = np.zeros(mesh.numel)
        node_x = mesh.node[:,0]
        node_y = mesh.node[:,1]
        node_z = mesh.node[:,2]
        # now find volume according to cell type
        if int(mesh.cell_type[0]) == 11: # elements are voxels
            for i in range(mesh.numel):
                p = mesh.node[con_mat[i]].T
                #compute volume (which is a bit of an approximation)
                dx = abs(max(p[0]) - min(p[0]))
                dy = abs(max(p[1]) - min(p[1]))
                dz = abs(max(p[2]) - min(p[2]))
                elm_vol[i] = dx*dy*dz
        elif int(mesh.cell_type[0]) == 10:# elements are tetrahedra 
            for i in range(mesh.numel):
                p = mesh.node[con_mat[i]].T    
                P = p[:,0] # point 1 
                Q = p[:,1] # point 2 
                R = p[:,2] # point 3
                pq = Q - P # p to q vector 
                pr = R - P # p to r vector 
                v = np.cross(pq, pr)
                S = p[:,3] # point 4
                #calculate height
                PS= (S[0]-P[0]), (S[1]-P[1]), (S[2]-P[2])
                #Calculate Volume
                elm_vol[i]= (np.dot(v,PS))/6
                 
        elif int(mesh.cell_type[0]) == 13: # elements are 3d wedges
            for i in range(mesh.numel):
                n1=(node_x[con_mat[i][0]], node_y[con_mat[i][0]], node_z[con_mat[i][0]])
                n2=(node_x[con_mat[i][1]], node_y[con_mat[i][1]], node_z[con_mat[i][1]])
                n3=(node_x[con_mat[i][2]], node_y[con_mat[i][2]], node_z[con_mat[i][2]])
                n4=(node_x[con_mat[i][3]], node_y[con_mat[i][3]], node_z[con_mat[i][3]])
                n5=(node_x[con_mat[i][4]], node_y[con_mat[i][4]], node_z[con_mat[i][4]])
                n6=(node_x[con_mat[i][5]], node_y[con_mat[i][5]], node_z[con_mat[i][5]])
                #compute wedge volume by computing face area first
                base=(((n1[0]-n2[0])**2) + ((n1[1]-n2[1])**2))**0.5
                mid_pt=((n1[0]+n2[0])/2,(n1[1]+n2[1])/2)
                height=(((mid_pt[0]-n3[0])**2) + ((mid_pt[1]-n3[1])**2))**0.5
                area = 0.5*base*height
                p = np.array((n1,n2,n3,n4,n5,n6)).T
                dz = abs(max(p[2]) - min(p[2]))
                elm_vol[i] = area * dz

        mesh.df['Volume'] = elm_vol
        if vmin is None:
            vmin = np.min(mesh.df[attr].values)
        if vmax is None:
            vmax = np.max(mesh.df[attr].values)
        sel_attr = mesh.df[(mesh.df[attr] >= vmin) & (mesh.df[attr] <= vmax)]
        volume = np.sum(sel_attr['Volume'].values)
        
        return volume
    
    def loadResults(self, invdir):
        """Given working directory, will attempt to load the results of an
        already run inversion.
        
        Parameters
        ----------
        invdir : str
            Path to the inversion directory.
        
        Note
        ----
        This does not load the data files neither the mesh nor the settings so
        you can't run another inversion with that. It's just for display.
        """
        # get all files 
        files = os.listdir(invdir)
        
        # get typ
        self.typ = [f.split('.')[0] for f in files if f[-3:] == '.in'][0]

        # detect if time-lapse and assume reg_mode == 0
        if 'ref' in files:
            self.iTimeLapse = True
            self.param['reg_mode'] = 0
            
        # load surveys
        self.surveys = []
        if (self.typ == 'cR2') or (self.typ == 'cR3t'):
            ftype = 'ProtocolIP'
        else:
            ftype = 'ProtocolDC'
        if self.iTimeLapse:
            # split the protocol.dat
            # in pandas >= 1.4.0 header=None with first row with one column (nb of meas)
            # causes ParseError. to fix it we first read the number of rows from line 2
            with open(os.path.join(invdir, 'protocol.dat'), 'r') as f:
                f.readline()  # first wow, we don't care
                nline = len(f.readline().split('\t'))
            dfall = pd.read_csv(os.path.join(invdir, 'protocol.dat'),
                                sep='\t', header=None, names=np.arange(nline))
            idf = list(np.where(np.isnan(dfall[dfall.columns[-1]].values))[0])
            idf.append(len(dfall))
            dfs = [dfall.loc[idf[i]:idf[i+1]-1,:] for i in range(len(idf)-1)]
    
            # writing all protocol.dat
            files = []
            for i, df in enumerate(dfs):
                outputname = os.path.join(self.dirname, 'survey{:03d}.dat'.format(i))
                files.append(outputname)
                if self.typ[-1] == 't':
                    df2 = df[1:].astype({0:int, 1:int, 2:int,
                                         3:int, 4:int, 5:int,
                                         6:int, 7:int, 8:int})
                    df2 = df2.drop(10, axis=1) # discard res0 as parser doesn't support it
                else:
                    df2 = df[1:].astype({0:int, 1:int, 2:int,
                                         3:int, 4:int})
                    df2 = df2.drop(6, axis=1) # discard res0

                with open(outputname, 'w') as f:
                    f.write('{:d}\n'.format(int(df.values[0, 0])))
                    df2.to_csv(f, sep='\t', header=False, index=False)
                # header with line count already included
            
            fnames = [os.path.join(invdir, 'ref', 'protocol.dat')] + files
            self.createTimeLapseSurvey(fnames, ftype=ftype)
            self.param['reg_mode'] = 0 # assumed
        else:
            self.createSurvey(os.path.join(invdir, 'protocol.dat'), ftype=ftype)

        # get electrodes
        elec = pd.read_csv(os.path.join(invdir, 'electrodes.dat'), delim_whitespace=True, header=None)
        self.setElec(elec.values) # assuming none are buried
 
        # load mesh
        self.importMesh(os.path.join(invdir, 'mesh.msh'))
        # self.mesh = mt.vtk_import(os.path.join(invdir, 'f001_res.vtk'))
        
        # get results
        self.getResults(invdir)
    
    
    
    def getR2out(self):
        """Reat the .out file and parse its content.
        
        Returns
        -------
        Dataframe with the dataset name, and the RMS decrease for each iteration.
        """
        fname = os.path.join(self.dirname, self.typ + '.out')
        lines = []
        if self.iTimeLapse:    
            with open(os.path.join(self.dirname, 'ref', self.typ + '.out'), 'r') as f:
                lines += f.readlines()
        with open(fname, 'r') as f:
            lines += f.readlines()
        name = ''
        idataset = 0
        iiter = 0
        resRMS = np.nan
        phaseRMS = np.nan
        read = np.nan
        rejected = np.nan
        irow = 0
        df = pd.DataFrame(columns=['name', 'dataset', 'iteration', 'resRMS',
                                   'phaseRMS', 'read', 'rejected', 'success'])
        for x in lines:
            success = 'N/A'
            line = x.split()
            if len(line) > 1:
                if line[0] == 'Iteration':
                    iiter += 1
                elif (line[0] == 'Measurements') & (line[1] == 'read:'):
                    read = int(line[2])
                    rejected = int(line[5])
                elif line[0] == 'Final':
                    resRMS = float(line[3])
                    df.loc[irow, :] = [name, idataset, iiter, resRMS, phaseRMS,
                                       read, rejected, success]
                    irow += 1
                elif line[0] == 'FATAL:':
                    resRMS = np.nan
                elif line[0] == 'Processing':
                    iiter = 0
                    idataset += 1
                    if idataset <= len(self.surveys):
                        name = self.surveys[idataset - 1].name
                    else:
                        name = 'dataset{:03.0f}'.format(idataset)
                elif line[0] == 'End' and self.typ != 'R3t':
                    idataset -= 1
        df = df.apply(pd.to_numeric, errors='ignore').reset_index(drop=True)
        return df


    def showRMS(self, index=0, ax=None):
        """Show the RMS decrease for each iteration.
        
        Parameters
        ----------
        index : int, optional
            Index of the dataset for which to plot the RMS.
        ax : matplotlib axis, optional
            If provided, the graph will be plotted against it.
        """
        df = self.getR2out()
        idatasets = np.unique(df['dataset'])
        if ax is None:
            fig, ax = plt.subplots()
        ax.set_title(self.surveys[index].name)
        offset = 0
        for i in idatasets:
            ie = df['dataset'] == i
            ax.plot(offset + df[ie]['iteration'], df[ie]['resRMS'], '.-')
            offset += np.sum(ie)
        ax.set_xlabel('Iterations')
        ax.set_ylabel('RMS misfit')
        ax.set_xticks([])
        

    # def showSection(self, fname='', ax=None, ilog10=True, isen=False, figsize=(8,3)):
    #     """Show inverted section based on the `_res.dat``file instead of the
    #     `.vtk`.

    #     Parameters
    #     ----------
    #     fname : str, optional
    #         Name of the inverted `.dat` file produced by the inversion.
    #     ax : matplotlib axis, optional
    #         If specified, the graph will be plotted along `ax`.
    #     ilog10 : bool, optional
    #         If `True`, the log10 of the resistivity will be used.
    #     isen : bool, optional
    #         If `True`, sensitivity will be displayed as white transparent
    #         shade on top of the inverted section.
    #     figsize : tuple, optional
    #         Size of the figure.
    #     """
    #     print('showSection called (to be discarded in the futur)')
    #     if fname == '':
    #         fname = os.path.join(self.dirname, 'f001.dat')
    #     res = pd.read_csv(fname, delimiter=' *', header=None, engine='python').values
    #     lenx = len(np.unique(res[:,0]))
    #     leny = len(np.unique(res[:,1]))
    #     x = res[:,0].reshape((leny, lenx), order='F')
    #     y = res[:,1].reshape((leny, lenx), order='F')
    #     z = res[:,2].reshape((leny, lenx), order='F')
    #     if isen:
    #         sen = pd.read_csv(fname.replace('res','sen'), delimiter=' *', header=None, engine='python').values
    #         lenx = len(np.unique(sen[:,0]))
    #         leny = len(np.unique(sen[:,1]))
    #         zs = sen[:,2].reshape((leny, lenx), order='F')
    #         zs = np.log10(zs)
    #         zs -= np.min(zs)
    #         alpha = zs/np.max(zs)
    #         print(np.max(alpha), np.min(alpha))
    #     if ilog10:
    #         z = np.log10(z)
    #     if ax is None:
    #         fig, ax = plt.subplots(figsize=figsize)
    #     else:
    #         fig = ax.get_figure()
    #     cax = ax.pcolormesh(x, y, z)
    #     ax.plot(self.elec[:,0], self.elec[:,2], 'ko')
    #     cbar = fig.colorbar(cax, ax=ax)
    #     if ilog10:
    #         cbar.set_label(r'$\log_{10}(\rho) [\Omega.m]$')
    #     else:
    #         cbar.set_label(r'$\rho [\Omega.m]$')
    #     ax.set_ylabel('Depth [m]')
    #     ax.set_xlabel('Distance [m]')


    def addRegion(self, xz, res0=100, phase0=1, blocky=False, fixed=False,
                  ax=None, iplot=False):
        """Add region according to a polyline defined by `xz` and assign it
        the starting resistivity `res0`.

        Parameters
        ----------
        xz : array
            Array with two columns for the x and y coordinates.
        res0 : float, optional
            Resistivity values of the defined area.
        phase0 : float, optional
            Read only if you choose the cR2 option. Phase value of the defined
            area in mrad
        blocky : bool, optional
            If `True` the boundary of the region will be blocky if inversion
            is block inversion.
        fixed : bool, optional
            If `True`, the inversion will keep the starting resistivity of this
            region.
        ax : matplotlib.axes.Axes
            If not `None`, the region will be plotted against this axes.
        iplot : bool, optional
            If `True` , the updated mesh with the region will be plotted.
        """

        centroids = self.mesh.elmCentre[:,[0,2]]
        path = mpath.Path(np.array(xz))
        idx = path.contains_points(centroids)

        region = np.array(self.mesh.df['region'])
        regid = np.max(region) + 1
        region[idx] = regid
        self.mesh.df['region'] = region
        resist0 = np.array(self.mesh.df['res0'])
        resist0[idx] = res0
        self.mesh.df['res0'] = resist0
        phase = np.array(self.mesh.df['phase0'])
        phase[idx] = phase0
        self.mesh.df['phase0'] = phase

        # define zone
        if blocky is True:
            zones = np.array(self.mesh.df['zones'])
            zones[idx] = regid
            self.mesh.df['zones'] = zones

        # define fixed area
        if fixed is True:
            paramFixed = np.array(self.mesh.df['param'])
            paramFixed[idx] = 0
            self.mesh.df['param'] = list(paramFixed)
            print('fixed {:d} elements'.format(np.sum(paramFixed == 0)))

        if iplot is True:
            self.showMesh()


    def resetRegions(self):
        """Just reset all regions already drawn. Shouldn't be needed as
        the `self.runR2()` automatically use a homogenous model when starting
        for inversion. The only purpose of this is to use an inhomogeous
        starting model to invert data from forward modelling.
        """
        self.mesh.df['region'] = np.ones(self.mesh.numel)
        self.mesh.df['res0'] = np.ones(self.mesh.numel)*100 # set back as default


    def computeAttribute(self, formula, name, dump=None):
        """Compute a new attribute for each meshResults.
        
        NOTE: this function present a security risk as it allows execution
        of python code inputed by the user. It should be used with caution.
        
        Parameters
        ----------
        formula : str
            Formula as a string. All attribute already in the mesh object are
            available from the x dictionary and can be sed as x['nameOfAttribute'].
            e.g. : "1/x['Resistivity']" will compute the inverse of the 'Resistivity'
            attribute if this attribute exists.
            Operators available are addition (+), soustraction (-), division (/)
            multiplication (*) and exponent (**). Parenthesis are also taken into account.
            numpy function such as np.sqrt() or np.cos() are also available.
            e.g. : "(np.sqrt(x['Resistivity']) + 100)*1000"
            Note: use a mixture of double and single quote such as double for
            the entire string and single for indexing the dictionnary x['keyName'].
        name : str
            Name of the new attribute computed.
        dump : function, optional
            Function to write stdout.
        """
        if dump is None:
            def dump(x):
                print(x, end='')
                
        for i, m in enumerate(self.meshResults):
            cols = m.df.columns
            vals = [m.df[c].values for c in cols]
            x = dict(zip(cols, vals))
            # DANGER ZONE =================================
            try:
                m.df[name] = eval(formula)
                dump('{:s} computation successful on meshResults[{:d}]\n'.format(name, i))
            except Exception as e:
                dump('{:s} computation failed on meshResults[{:d}]: {:s}\n'.format(name, i, str(e)))
                pass
            # DANGER ZONE =================================
            


    def createModel(self, ax=None, dump=None, typ='poly', addAction=None):
        """Interactive model creation for forward modelling.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to which the graph will be plotted.
        dump : function, optional
            Function that outputs messages from the interactive model creation.
        typ : str
            Type of selection either `poly` for polyline or `rect` for
            rectangle.
        addAction : function
            Function to be called once the selection is finished (design for
            GUI purpose).

        Returns
        -------
        fig : matplotlib.figure
            If `ax` is `None`, will return a figure.
        """
        if dump is None:
            def dump(x):
                print(x, end='')
                
        if self.mesh is None:
            print('will create a mesh before')
            self.createMesh()
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        def callback(idx): # :pragram: no cover
            region = np.array(self.mesh.df['region']).copy()
            regid = np.max(region) + 1
            print('nb elements selected:', np.sum(idx), 'in region', regid)
            region[idx] = regid
            self.mesh.df['region'] = region
            self.mesh.draw(attr='region')
            if addAction is not None:
                addAction()
        self.mesh.atribute_title = 'Regions'
        if self.zlim is None:
            self._defineZlim()
        self.mesh.show(attr='region', ax=ax, zlim=self.zlim, darkMode=self.darkMode)
        # we need to assign a selector to self otherwise it's not used
        self.selector = SelectPoints(ax, self.mesh.elmCentre[:,[0,2]],
                                     typ=typ, callback=callback)
        if ax is None:
            return fig



    def designModel(self, ax=None, dump=print, typ='poly', addAction=None, fmd=None):
        """Interactive model design for forward modelling (triangular only).
        As opposite to R2.createModel(). R2.designModel() allows to draw mesh
        region **before** meshing. This allows to have straight boundaries for
        triangular mesh.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to which the graph will be plotted.
        dump : function, optional
            Function that outputs messages from the interactive model creation.
        typ : str
            Type of selection either `poly` for polyline or `rect` for
            rectangle.
        addAction : function
            Function to be called once the selection is finished (design for
            GUI purpose).
        fmd : float, optional
            Depth of of interest specifies as a relative positive number from the surface.

        Returns
        -------
        fig : matplotlib.figure
            If `ax` is `None`, will return a figure.
        """
        if fmd is None:
            self.computeFineMeshDepth()
        else:
            self.fmd = fmd

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        
        dump('Fine mesh depth in deisgned model is %f'%float(self.fmd))
        self.geom_input = {}
        elecColor = 'ko' if self.darkMode is False else 'wo'
        ax.plot(self.elec['x'], self.elec['z'], elecColor, label='electrode')
        ax.set_ylim([np.min(self.elec['z']) - self.fmd, np.max(self.elec['z'])])
        ax.set_xlim(np.min(self.elec['x']), np.max(self.elec['x']))
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('Elevation [m]')
        def callback(): # :pragma: no cover
            vert = np.array(self.selector.vertices)
            self.geom_input['polygon' + str(len(self.geom_input)+1)] = [vert[:-1,0].tolist(), vert[:-1,1].tolist()]
            ax.plot(vert[:,0], vert[:,1], '.-')
            if addAction is not None:
                addAction()
        # we need to assign a selector to self otherwise it's not used
        self.selector = SelectPoints(ax, typ=typ, callback=callback)
#        surveyLength = np.max(self.elec[:,0]) - np.min(self.elec[:,0])
        self.selector.xmin = np.min(self.elec['x'])# - 10 * surveyLength
        self.selector.xmax = np.max(self.elec['x'])# + 10 * surveyLength
        if ax is None:
            return fig


    def createModelMesh(self, **kwargs):
        """Create a triangular mesh given the designed geometry by
        R2.designModel().

        Parameters
        ----------
        All parameters to be passed are similar to `R2.createMesh()`.
        """
        geom_input = self.geom_input
        self.createMesh(typ='trian', geom_input=geom_input, **kwargs)


    def setStartingRes(self, regionValues={}, zoneValues={}, fixedValues={}, ipValues={}):
        """Assign starting resitivity values.

        Parameters
        ----------
        regionValues : dict, optionals where it left off? I'm assuming it could be forced by importing the baseline data se
            Dictionnary with key being the region number and the value being
            the resistivity in [Ohm.m].
        zoneValues : dict, optional
            Dictionnary with key being the region number and the zone number.
            There would be no smoothing between the zones if 'block inversion'
            is selected (`inversion_type` = 4).
        fixedValues : dict, optional
            Dictionnary with key being the region number and a boolean value if
            we want to fix the resistivity of the zone to the starting one.
            Note that it only works for triangular mesh for now.
        ipValues : dict, optional
            Dictionnary with key being the region number and the values beeing
            the phase [mrad].

        Note
        ----
        Region 0 is the background region. It has zone=1, and fixed=False
        """
        regions = np.array(self.mesh.df['region'])
        res0 = np.array(self.mesh.df['res0']).copy()
        for key in regionValues.keys():
            idx = regions == int(key)
            res0[idx] = regionValues[key]
            self.param['res0File'] = 'res0.dat'
        self.mesh.df['res0'] = res0
        # print('regionValues:',regionValues)

        zones = np.array(self.mesh.df['zones']).copy()
        for key in zoneValues.keys():
            idx = regions == int(key)
            zones[idx] = zoneValues[key]
        self.mesh.df['zones'] = zones

        fixed = self.mesh.df['param'].values.copy()
        for key in fixedValues.keys():
            idx = regions == int(key)
            if fixedValues[key] == True:
                fixed[idx] = 0
        self.mesh.df['param'] = fixed

        phase0 = np.array(self.mesh.df['phase0']).copy()
        for key in ipValues.keys():
            idx = regions == int(key)
            phase0[idx] = ipValues[key]
        self.mesh.df['phase0'] = phase0


    def setRefModel(self, res0):
        """Set the reference model according to a previous inversion, avoids
        the need to invert reference model again for timelapse workflows.
        In contrast to `R2.setStartingRes()` which assign resistivity to group
        of elements, this method requires a vector of the same length as the 
        number of elements. This enables, notably to manually perform consecutive
        background constrained inversion.

        Parameters
        ----------
        res0: array like
            Array of resistivity values, ideally from a previous inversion. The
            length of this array should be the same as the number of elements.
        """
        try:
            self.mesh.addAttribute(res0,'res0')
        except AttributeError:
            print('Cant set reference model without first assigning/creating a mesh')
            return
        self.param['reg_mode'] = 1 # ensure inversion is background regularised
        if self.typ[-1] =='t':
            self.param['inverse_type']=1
        self.param['res0File'] = 'Start_res.dat'
        self.param['num_regions'] = 0
        self.mesh.writeAttr('res0',os.path.join(self.dirname,'Start_res.dat'))
        self.referenceMdl = True
        print('Reference model successfully assigned')


    def _seqIdxFromLabel(self):
        lines = [int(a.split(' ')[0]) for a in self.elec['label'].values]
        uline = np.unique(lines)
        nline = len(uline)
        seqIdx = []
        for i in range(nline):
            lidx = np.argwhere(lines==uline[i])
            idx = [0]*len(lidx)
            for j in range(len(lidx)):
                idx[j]=lidx[j][0]
            seqIdx.append(np.array(idx))
        return seqIdx
    
    
    def _genSeqIdx(self):
        if '3' in self.typ: 
            if not self.hasElecString():#then find it automatically 
                seqIdx = self.detectStrings()
            else:
                seqIdx = self._seqIdxFromLabel()#use electrode strings 
        else:
            seqIdx = [np.arange(len(self.elec))]
            
        return seqIdx
    
    
    def createSequence(self, params=[('dpdp',1,8,1,8)], seqIdx=None, dump=print):
        """Creates a forward modelling sequence, see examples below for usage.

        Parameters
        ----------
        params : list of tuple, optional
            Each tuple is the form (<array_name>, param1, param2, ...)
            Types of sequences available are : 	
            - 'dipole-dipole' (or 'dpdp')
            	- 'wenner' (or 'w') 
            	- 'wenner-schlumberger' (or 'schlum', 'ws')
            	- 'multigradient' (or 'mg')
            - 'cross' (or 'xbh', 'xli')
            	- 'custom' 
            if 'custom' is chosen, param1 should be a string of file path to a .csv
            file containing a custom sequence with 4 columns (a, b, m, n) containing 
            forward model sequence.
        seqIdx: list of array like, optional
            Each entry in list contains electrode indices (not label and string)
            for a given electrode string which is to be sequenced. The advantage
            of a list means that sequences can be of different lengths. 
        
        Examples
        --------
        >>> k = Project()
        >>> k.setElec(np.c_[np.linspace(0,5.75, 24), np.zeros((24, 2))])
        >>> k.createMesh(typ='trian')
        >>> k.createSequence([('dpdp', 1, 8, 1, 8), ('wenner', 1, 4), ('ws', 1, 8, 1, 8)]) # dipole-dipole sequence
        >>> k.createSequence([('custom', '<path to sequence file>/sequence.csv')]) # importing a custom sequence
        >>> seqIdx = [[0,1,2,3],[4,5,6,7],[8,9,10,11,12]]
        """
        availconfig = [
            'dipole-dipole', 'dpdp', 
            'wenner', 'w',
            'wenner-schlumberger', 'schlum', 'ws', 
            'multigradient', 'mg', 
            'cross','xbh', 'xli', 'equat-dp', 
            'custom', 'custSeq']
        
        msg = 'Sorry, given config parameter {:s} is not recognised!'
        msg += 'Available configs are = \n'
        for config in availconfig:
            msg += '\t%s\n'%config 
        for p in params: 
            if p[0] not in availconfig: 
                warnings.warn(msg.format(p[0])) 

        if seqIdx is None: #(not been set, so use electrode strings)
            seqIdx = self._genSeqIdx()
                
        self.sequenceGenerator = Generator(self.elec, seqIdx)
        
        sequence = self.sequenceGenerator.generate(params, dump=dump)
            
        self.sequence = sequence
        print('{:d} quadrupoles generated.'.format(self.sequence.shape[0]))
        return seqIdx
        
    def _createSequence(self, params=[('dpdp1', 1, 8)], seqIdx=None,
                       *kwargs):
        """Creates a forward modelling sequence, see examples below for usage.

        Parameters
        ----------
        params : list of tuple, optional
            Each tuple is the form (<array_name>, param1, param2, ...)
            Types of sequences available are : 'dpdp1','dpdp2','wenner_alpha',
            'wenner_beta', 'wenner_gamma', 'schlum1', 'schlum2', 'multigrad',
            'custSeq'.
            if 'custSeq' is chosen, param1 should be a string of file path to a .csv
            file containing a custom sequence with 4 columns (a, b, m, n) containing 
            forward model sequence.
        seqIdx: list of array like, optional
            Each entry in list contains electrode indices (not label and string)
            for a given electrode string which is to be sequenced. The advantage
            of a list means that sequences can be of different lengths. 
        
        Examples
        --------
        >>> k = Project()
        >>> k.setElec(np.c_[np.linspace(0,5.75, 24), np.zeros((24, 2))])
        >>> k.createMesh(typ='trian')
        >>> k.createSequence([('dpdp1', 1, 8), ('wenner_alpha', 1), ('wenner_alpha', 2)]) # dipole-dipole sequence
        >>> k.createSequence([('custSeq', '<path to sequence file>/sequence.csv')]) # importing a custom sequence
        >>> seqIdx = [[0,1,2,3],[4,5,6,7],[8,9,10,11,12]]
        """
        
        fdico = {'dpdp1': dpdp1, # dict referencing the seqeunce gen functions 
              'dpdp2': dpdp2,
              'wenner': wenner,
              'wenner_alpha': wenner_alpha,
              'wenner_beta': wenner_beta,
              'wenner_gamma': wenner_gamma,
              'schlum1': schlum1,
              'schlum2': schlum2,
              'multigrad': multigrad,
              'custSeq': addCustSeq}

        self.custSeq = False # reset the flag
        if (self.typ == 'R3t') | (self.typ == 'cR3t'): #its 3D
            for p in params:
                if p[0] == 'custSeq' and p[1] != '':
                    self.custSeq = True # we can't have custom sequence mixed with auto generated sequences in 3D - too complicated!
                    print('Custom sequence detected. Skipping auto sequence generation...')
                    sequence = addCustSeq(p[1]).astype(str)
                    if len(sequence[0][0].split(' ')) == 1: # no string, so let's add them
                        for i in range(sequence.shape[0]):
                            for j in range(sequence.shape[1]):
                                sequence[i,j] = '1 '+ str(sequence[i,j])
                    break # only one custom sequence allowed? 
                
            if self.custSeq is False:                    
                #determine sequence index if not already given 
                if seqIdx is None: #(not been set, so use electrode strings)
                    if not self.hasElecString():#then find it automatically 
                        seqIdx = self.detectStrings(*kwargs)
                        #raise ValueError('Electrode strings have not been set')
                    else:
                        seqIdx = self._seqIdxFromLabel()#use electrode strings 
    
                elif type(seqIdx) != list: # check we have a list 
                    raise TypeError('Expected list type argument for seqIdx')
                
                #sequentially go through and create sequences for each electrode string
                qs = []
                nseq = len(seqIdx) # number of sequences (ie num of lines strings)
                for i in range(nseq):
                    nelec = len(seqIdx[i])
                    slabels = self.elec['label'].values[seqIdx[i]] # string sequence labels
                    for p in params:
                        if p[0] != 'custSeq': # ignore custSeq as already dealt with
                            pok = [int(p[i]) for i in np.arange(1, len(p))] # make sure all are int
                            str_seq = fdico[p[0]](nelec, *pok).values.astype(int) # return np array 
                            qs.append(slabels[str_seq-1])
                            
                sequence = np.vstack(qs)
                if not self.hasElecString():
                    # if it's 3D, we need the line number for the sequence 
                    # (so all electrode on line 1)
                    for i in range(sequence.shape[0]):
                        for j in range(sequence.shape[1]):
                            sequence[i,j] = '1 '+ sequence[i,j]
   
                    
        else: # it's 2D 
            nelec = self.elec.shape[0] 
            qs = []
            for p in params:
                if p[0] == 'custSeq':
                    try:
                        qs.append(addCustSeq(p[1]))
                    except Exception as e:
                        print('error when importing custom sequence:', e)
                else:
                    pok = [int(p[i]) for i in np.arange(1, len(p))] # make sure all are int
                    qs.append(fdico[p[0]](nelec, *pok).values.astype(int))
            sequence = np.vstack(qs)
            
            # detecing quadrupoles using out of bound electrodes
            iabove = (sequence > self.elec.shape[0]).any(1)
            sequence = sequence[~iabove,:]
            
        self.sequence = sequence
        print('{:d} quadrupoles generated.'.format(self.sequence.shape[0]))
        return seqIdx

    def saveSequence(self, fname='sequence.csv', ftype='generic', integer = True, 
                     reciprocals = True, split = False, multichannel = True,  
                     condition = True, maxcha = 8, dump=print):
        """Save sequence as .csv file. Ex

        Parameters
        ----------
        fname : str, optional
            Path where to save the sequence.
        ftype: str, optional 
            Flag type of command file. Default is 'generic'. 
        integer : bool, optional
            Flag to convert sequence into integers before export. The default is True.
            String numbers will be appended to the start of the electrode indexes if 
            present. 
        reciprocals : bool, optional
            Flag to add reciprocal measurements. The default is True.
        split: bool, optional
            Flag to split reciprocal measurements into seperate files, default 
            is False. 
        multichannel : bool, optional
            Flag to convert measurements to a multichannel. The default is True.
        condition : bool, optional
            Flag to condition measurements for avoiding IP effects. The default is True.
        maxcha: int, optional 
            Maximum number of active channels of the resistivity instrument (normally
            8 for modern instruments). 
        """

        if self.sequence is not None:
            self.sequenceGenerator.exportSequence(
                fname, ftype, integer, split,
                reciprocals, multichannel, 
                condition, maxcha)
            

    def importElec(self, fname=''):
        """Import electrodes positions.

        Parameters
        ----------
        fname : str
            Path of the CSV  (or file containing the electrodes positions. It 
            should contains 3 columns maximum with the X, Y, Z positions 
            of the electrodes.
        """
        if fname.endswith('.geom'): # read in geometry file (special case) 
            elec = geomParser(fname)
            self.setElec(elec)
            return 
        
        with open(fname, 'r') as f:
            try:
                float(f.readline().split(',')[0])
                header = None
            except Exception:
                header = 'infer'
        df = pd.read_csv(fname, header=header)
        if header is None:
            elec = df.values
        else:
            elec = df
        self.setElec(elec)
                

    def importSequence(self, fname=''):
        """Import sequence for forward modelling.

        Parameters
        ----------
        fname : str
            Path of the CSV file to be imported. The file must have 4 columns with headers (a, b, m, n) containing 4 electrodes numbers.
        """
        seq = pd.read_csv(fname, header=0)
        keys = seq.keys()
        
        #check for headers 
        trigger = False 
        for a in ['a','b','m','n']:
            if a not in keys:
                print('Column "%s" not in sequence file'%a)
                trigger =True 
        if trigger:
            raise Exception('Missing headers in sequence file!')
                
        if seq.shape[1] != 4:
            raise ValueError('The file should be a CSV file with headers (a, b, m, n) with exactly 4 columns with electrode numbers.')
        else:
            self.sequence = seq
            
        #do check for electrode line numbers in the case of 3D 
        if self.typ[-1] == 't':
            if not isinstance(self.sequence['a'][0],str): # then add line numbers 
                # print('adding line numbers')
                surrogate = pd.DataFrame()
                a = ['1']*len(seq)
                b = ['1']*len(seq)
                m = ['1']*len(seq)
                n = ['1']*len(seq)
                for i in range(len(seq)):
                    a[i] = '1 '+ str(seq['a'][i])
                    b[i] = '1 '+ str(seq['b'][i])
                    m[i] = '1 '+ str(seq['m'][i])
                    n[i] = '1 '+ str(seq['n'][i])
                surrogate['a'] = a
                surrogate['b'] = b
                surrogate['m'] = m
                surrogate['n'] = n
                self.sequence = surrogate 


    def saveErrorData(self, fname):
        """Save quadruople, resistance, phase and their respective reciprocal
        errors as .csv file.

        Parameters
        ----------
        fname : str
            Path where to save the file.
        """
        cols = np.array(['a','b','m','n','resist','recipMean','recipError','resError',
                         'phase','reci_IP_err','phaseError'])
        if self.iTimeLapse is True:
            df = self.bigSurvey.df
        else:
            df = self.surveys[0].df
        ie = [c in df.columns for c in cols]
        dff = df[cols[ie]]
        dff = dff.rename(columns = {'resist':'Resistance [ohm]', 'recipError':'Resistance_err [ohm]',
                                    'resError':'Fit Resistance_err [ohm]','phase':'Phase [mRad]',
                                    'reci_IP_err':'Phase_err [mRad]','phaseError':'Fit Phase_err [mRad]'})
        dff.to_csv(fname, index=False)


    # def saveFilteredData(self, fname, savetyp='Res2DInv (*.dat)'):
    #     """Save filtered data in formats to be used outside ResIPy (e.g. Res2DInv).

    #     Parameters
    #     ----------
    #     fname : str
    #         Path where to save the file.
    #     savetyp : str, optional
    #         Saving format. To be determined in GUI. Default: Res2DInv (.dat)
    #     """
    #     elec = self.elec[['x','y','z']].values
    #     spacing = elec[1,0] - elec[0,0] # TODO (gb) not sure if this is needed
    #     for s, i in zip(self.surveys, range(len(self.surveys))):
    #         df = s.df.query('irecip >=0').copy() # not saving reciprocal data
    #         # if spacing == None:
    #         #     spacing = elec[1,0]-elec[0,0] # for batch surveys the spacing can differ and not follow user input
    #         # else:
    #         #     spacing = spacing
    #         # df[['a','b','m','n']] *= spacing
    #         lookupDict = dict(zip(self.elec['label'], self.elec['x'].values))
    #         data = df[['a','b','m','n']].replace(lookupDict).values
    #         for i, a in enumerate(['a','b','m','n']):
    #             df[a] = data[:,i] # this way should avoid warnings 
    #         # df[['a','b','m','n']] = data
    #         if savetyp == 'Res2DInv (*.dat)':
    #             param = {'num_meas': df.shape[0],
    #                      'lineTitle': self.param['lineTitle'],
    #                      'spacing': spacing}
    #             write2Res2DInv(param, fname, df, elec, self.typ)
    #         elif savetyp == 'Comma Separated Values (*.csv)':
    #             write2csv(fname, df, elec, self.typ)
    #         elif savetyp == 'E4D survey file (*.srv)':
    #             writeSrv(fname, df, elec)

    #         fname = fname[:-4] + str(i) + fname[-4:] # to iterate file numbers in case of timelapse survey


    def forward(self, noise=0.0, noiseIP=0.0, iplot=False, dump=None):
        """Operates forward modelling.

        Parameters
        ----------
        noise : float, optional 0% <= noise <= 100%
            Noise level in percent from a Gaussian distribution that should be
            applied on the forward apparent resistivities obtained.
        noiseIP : float, optional
            Absolute noise level in mrad from a Gaussian distribution that should be applied
            on the forward phase values obtained.
        iplot : bool, optional
            If `True` will plot the pseudo section after the forward modelling.
        dump : function, optional
            Function to print information messages when running the forward model.
        """
        if dump is None:
            def dump(x):
                print(x, end='')
                
        fwdDir = os.path.join(self.dirname, 'fwd')
        if os.path.exists(fwdDir):
            shutil.rmtree(fwdDir)
        os.mkdir(fwdDir)

        # no need to order the mesh in forward as zone and param are not read

        if self.typ[0] == 'c':
            r = np.array(self.mesh.df['res0'])
            phase = np.array(self.mesh.df['phase0'])
            centroids = self.mesh.elmCentre.copy()
            centroids2 = centroids[:,[0,2]] if self.typ[-1] != 't' else centroids
            x = np.c_[centroids2,
                      r,
                      phase, # mrad
                      np.log10(r),
                      np.log10(np.cos(-phase/1000)/np.log10(r)), #log10(real conductivity)
                      np.log10(np.sin(-phase/1000)/np.log10(r))] #log10(imaginary conductivity)
            np.savetxt(os.path.join(fwdDir, 'resistivity.dat'), x)
        else:
            self.mesh.writeAttr('res0', os.path.join(fwdDir,'resistivity.dat'))

        # write mesh.dat (no ordering of elements needed in forward mode)
        if (self.typ == 'R2') | (self.typ == 'cR2'):
            self.mesh.dat(os.path.join(fwdDir, 'mesh.dat'))
        else:
            self.mesh.datAdv(os.path.join(fwdDir, 'mesh3d.dat'), iadvanced=self.iadvanced)
            
        # write the forward .in file
        if self.custSeq is True: # only in case of 3D custom sequence
            # replacing labels with <1 elec_num> for simplicity - ignoring the initial self.elec['label'] values
            ones = np.ones_like(self.elec['x'].values.astype(str))
            newLabels = np.char.add(np.char.replace(ones, '1', '1 '), np.arange(1, len(ones)+1).astype(str))
            self.elec['label'] = newLabels
            self.param['node_elec'][0] = newLabels
            
        dump('Writing .in file and mesh.dat... ')
        fparam = self.param.copy()
        fparam['job_type'] = 0
        fparam['num_regions'] = 0
        fparam['res0File'] = 'resistivity.dat' # just starting resistivity
        
        write2in(fparam, fwdDir, typ=self.typ)
        dump('done\n')

        # write the protocol.dat (that contains the sequence)
        if self.sequence is None:
            dump('Creating sequence... ')
            self.createSequence()
            dump('done\n')
        elif isinstance(self.sequence, pd.core.frame.DataFrame):
            if len(self.sequence.columns) != 4: 
                raise Exception('Sequence should have 4 columns')
        elif isinstance(self.sequence,np.ndarray): # check sequence is a matrix of string etc
            if self.sequence.shape[1] != 4: 
                raise Exception('Sequence should have 4 columns')
            if isinstance(self.sequence[0,0], np.float32) or isinstance(self.sequence[0,0], np.float64):
                self.sequence = np.asarray(self.sequence,dtype=int)
            if isinstance(self.sequence[0,0], np.int32) or isinstance(self.sequence[0,0], np.int64):
                self.sequence = np.char.mod('%d',self.sequence) # reassign to string 
        else: 
            raise Exception('Sequence should be a numpy.ndarray or pandas.DataFrame')
            
        dump('Writing protocol.dat... ')
        seq = self.sequence

        ## let's check if IP that we have a positive geometric factor
        ## EH? Why do this - Jimmy
        ## commented out for now as not a generalised solution and did not work 
        ## for 3D problems. 
        # if self.typ[0] == 'c': # NOTE this doesn't work for borehole
        #     lookupDict = dict(zip(self.elec['label'], np.arange(self.elec.shape[0])))
        #     seqdf = pd.DataFrame(seq).rename(columns = {0:'a', 1:'b', 2:'m', 3:'n'})
        #     array = seqdf[['a','b','m','n']].astype(str).replace(lookupDict).values.astype(int)
        #     elec = self.elec[['x','y','z']].values.copy()
            
        #     aposx = elec[:,0][array[:,0]]
        #     aposy = elec[:,1][array[:,0]]
        #     aposz = elec[:,2][array[:,0]]
            
        #     bposx = elec[:,0][array[:,1]]
        #     bposy = elec[:,1][array[:,1]]
        #     bposz = elec[:,2][array[:,1]]
            
        #     mposx = elec[:,0][array[:,2]]
        #     mposy = elec[:,1][array[:,2]]
        #     mposz = elec[:,2][array[:,2]]
            
        #     nposx = elec[:,0][array[:,3]]
        #     nposy = elec[:,1][array[:,3]]
        #     nposz = elec[:,2][array[:,3]]
            
        #     AM = np.sqrt((aposx-mposx)**2 + (aposy-mposy)**2 + (aposz-mposz)**2)
        #     BM = np.sqrt((bposx-mposx)**2 + (bposy-mposy)**2 + (bposz-mposz)**2)
        #     AN = np.sqrt((aposx-nposx)**2 + (aposy-nposy)**2 + (aposz-nposz)**2)
        #     BN = np.sqrt((bposx-nposx)**2 + (bposy-nposy)**2 + (bposz-nposz)**2)
            
        #     K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
            
        #     ie = K < 0
        #     seq2 = seq.copy()
        #     seq[ie,2] = seq2[ie,3] # swap if K is < 0
        #     seq[ie,3] = seq2[ie,2]
            
        protocol = pd.DataFrame(np.c_[1+np.arange(seq.shape[0]),seq],
                                columns=['index','a','b','m','n'])
        
        # if it's 3D, we add the line number (all electrode on line 1)
        if ((self.typ == 'R3t') or (self.typ == 'cR3t')):
            if len(protocol['a'].values[0].split()) == 1: # we don't have string number
                for c in ['a','b','m','n']: 
                    protocol[c] = '1 ' + protocol[c]
        
        outputname = os.path.join(fwdDir, 'protocol.dat')
        with open(outputname, 'w') as f:
            f.write(str(len(protocol)) + '\n')
        with open(outputname, 'a') as f:
            protocol.to_csv(f, sep='\t', header=False, index=False,lineterminator='\n')
        dump('done\n')

        # fun the inversion
        dump('Running forward model... ')
        self.runR2(fwdDir, dump=dump) # this will copy the R2.exe inside as well
        self.iForward = True

        # create a protocol.dat file (overwrite the method)
        def addnoise(x, level=0.05):
            return x + np.random.randn(1)*x*level

        def addnoiseIP(x, level=2):
            return x + np.random.randn(1)*level

        addnoise = np.vectorize(addnoise)
        addnoiseIP = np.vectorize(addnoiseIP)
        self.noise = noise # percentage noise e.g. 5 -> 5% noise
        self.noiseIP = noiseIP #absolute noise in mrad, following convention of cR2
        
        #fmd = self.fmd#.copy()
        elec = self.elec.copy()
        if self.typ[-1]=='t' and not self.hasElecString():
            #need to add elec strings to labels if in 3D
            for i in range(elec.shape[0]):
                elec.loc[i,'label'] = '1 ' + elec['label'][i]
            
        self.surveys = [] # need to flush it (so no timeLapse forward)
        if self.typ[0] == 'c':
            self.createSurvey(os.path.join(fwdDir, self.typ + '_forward.dat'), 
                              ftype='forwardProtocolIP',
                              compRecip=False) # dont compute reciprocals as that will be done after adding noise (see lines below)
        else:
            self.createSurvey(os.path.join(fwdDir, self.typ + '_forward.dat'), 
                              ftype='forwardProtocolDC',
                              compRecip=False)
            
        # NOTE the 'ip' columns here is in PHASE not in chargeability
        self.surveys[0].kFactor = 1 # kFactor by default is = 1 now, though wouldn't hurt to have this here!
        self.surveys[0].df['resist'] = addnoise(self.surveys[0].df['resist'].values, self.noise/100)
        self.surveys[0].df['ip'] = addnoiseIP(self.surveys[0].df['ip'].values, self.noiseIP)
        self.surveys[0].computeReciprocal() # to recreate the other columns
        self.setElec(elec) # using R2.createSurvey() overwrite self.elec so we need to set it back
        self.surveys[0].computeK()  # need to recompute K with the new electrode given by setElec()
        self.surveys[0].df['app'] = self.surveys[0].df['K']*self.surveys[0].df['resist']  # and recompute app
        # self.fmd = fmd      

        # recompute doi (don't actually otherwise zlim is jumping)
        # self.computeFineMeshDepth()
        # self.zlim[0] = np.min(elec['z']) - self.fmd
        if iplot is True:
            self.showPseudo()
            
        # save results with noise 
        if noise > 0 or noiseIP > 0:
            outputname = os.path.join(fwdDir, self.typ + '_forward_w_noise.dat')
            if self.typ[0] == 'c':
                self.surveys[0].write2protocol(outputname,ip=True)
            else: 
                self.surveys[0].write2protocol(outputname)
                
        _ = self._estimateMemory(dump=donothing)
        
        dump('Forward modelling done.')
        
    def saveForwardModelResult(self,fname):
        """
        Save the result of a forward model run to a specific file name/ 

        Parameters
        ----------
        fname : str
            path to file. 
        """
        if self.iForward: 
            _ = self.surveys[0].write2protocol(fname)
        else:
            print('No forward model has been run!')

    def createModelErrorMesh(self, **kwargs):
        """Create an homogeneous mesh to compute modelling error.

        Same arguments as `Project.createMesh()`.
        """
        if self.wholespace:
            raise Exception('Model Error mesh function meant for halfspace problems only')
    
        cases = ['tetra','trian','quad'] # use cases where modelling error is a half space 
        if self.meshParams['typ'] not in cases: 
            raise Exception('Modelling error mesh not avialable for this mesh type yet')
            
        # backup
        elecZ = self.elec['z'].values.copy()
        mesh = self.mesh.copy() if self.mesh is not None else None
        zlim = self.zlim.copy()
        param = self.param.copy()
        
        kwargs['model_err'] = True # force model_err argument to be true 
        
        self.createMesh(**kwargs)
        self.modErrMesh = self.mesh.copy()
        self.modErrMeshNE = self.param['node_elec'].copy()
        self.param['num_regions'] = 0
        
        # restore
        self.mesh = mesh
        self.elec['z'] = elecZ
        self.zlim = zlim
        self.param = param
        

    def estimateError(self, a_wgt=0.001, b_wgt=0.02):
        """Estimate reciprocal error data for data with no reciprocals for each
        survey, using the same routine present in R2. This allows for the 
        additional inclusion of modelling errors. It could be used when the user
        want to assign invidual errors based on a_wgt/b_wgt. This action is irreversable.

        Parameters
        ----------
        a_wgt: float, optional
            a_wgt documented in the R2 documentation
        b_wgt: float, optional
            b_wgt documented in the R2 documentation
        """
        for s in self.surveys:
            s.estimateError(a_wgt=a_wgt, b_wgt=b_wgt)
            

    def addFlatError(self, percent=2.5):# TODO (gb) why would we want that? - see below (jb)
        """Add a flat percentage error to resistivity data (for each survey in
        the class). This action is irreversable. Use case is for simulating 
        modelling errors in unconventional surveys. 

        resError = res*(percent/100) + resError

        Parameters
        ----------
        percent : float
            Error in percent.
        """
        for s in self.surveys:
            s.addPerError(percent)
            

    def computeModelError(self, rmTree=False, dump=None, refRes=100, refPhs=-10):
        """Compute modelling error associated with the mesh for a half space
        problem. This is computed on a flat mesh with the same meshing parameters
        as the project mesh. In the case of non-conventional surveys a different
        modelling error scheme will be used. 

        Parameters
        ----------
        rmTree : bool
            Remove the working directory used for the error modelling. Default
            is True.
        dump : bool 
            Function to direct the output of the modelling error process, by default
            the outputs are printed to the console. 
        refRes: float, int 
            Reference resistivity value (normally 100)
        refPhs: float, int 
            Reference phase value (normally -10)
            
        """
        if dump is None:
            def dump(x):
                print(x, end='')
        
        if self.mesh is None:
            raise Exception('Cannot compute modelling error without a mesh, make a mesh first!')
                
        dump('Doing error modelling...\n')
        
        cases = ['tetra','trian','quad'] # use cases where modelling error is a half space 
        halfspace = True 
        if 'typ' in self.meshParams.keys():
            if self.meshParams['typ'] not in cases and not self.wholespace: 
                halfspace = False 
                
        node_elec = None # we need this as the node_elec with topo and without might be different
        if not halfspace:
            dump('Refining the mesh for the case of a generic problem...\n\n')
            mesh = self.mesh.refine() # refine the mesh (maybe even need 2x as much? kinda slow though)  
        elif any(self.elec['z'].values != 0): # so we have topography in this case and need a new mesh 
            dump('Creating mesh without ANY topography...\n\n')
            meshParams = self.meshParams.copy()
            meshParams['model_err'] = True 
            # if 'geom_input' in meshParams: # dont use geometry from here because it'll likley be incompatible on the flat mesh
            #     meshParams['geom_input'] = {}
                
            # NB: Don't use the same dump from the mesh creation tab in the UI 
            if 'show_output' in meshParams: 
                if meshParams['show_output']:
                    meshParams['dump'] = dump 
                else: 
                    meshParams['dump'] = None 
                    
            self.createModelErrorMesh(**meshParams)
            node_elec = self.modErrMeshNE
            mesh = self.modErrMesh # create flat mesh
        else:
            mesh = self.mesh # use same mesh

        # create working directory
        dump('Making temporary error modelling directory...\n')
        fwdDir = os.path.join(self.dirname, 'err')
        if os.path.exists(fwdDir):
            shutil.rmtree(fwdDir)
        os.mkdir(fwdDir)

        # write the resistivity.dat and fparam
        fparam = self.param.copy()
        fparam['job_type'] = 0
        
        if '2' in self.typ:
            ndim = 2
            mname = 'mesh.dat'
            meshPath = os.path.join(fwdDir, mname)
            mesh.dat(meshPath)
            appResCol = 6
            trCol = 5 
        else:
            ndim = 3
            mname = 'mesh3d.dat'
            meshPath = os.path.join(fwdDir, mname)
            mesh.datAdv(meshPath, iadvanced=self.iadvanced) # use advanced mesh format if 3D 
            appResCol = 10
            trCol = 9
        
        #check to see if an extra column needs adding in the case of IP surveys 
        ncol = 1 
        phaseCol = 0 
        ipflag = False 
        if 'c' in self.typ: 
            ncol = 2
            appResCol = trCol + 2 
            phaseCol = trCol + 1 
            ipflag = True 
    
        #make starting resistivity file 
        
        resFile = np.zeros((mesh.numel, ndim+ncol),dtype=float) # centroid x, y, z, res0, (phase)
        resFile[:,ndim] = refRes
        if ipflag:  
            resFile[:,ndim+1] = refPhs 
        np.savetxt(os.path.join(fwdDir, 'resistivity.dat'), resFile,
                   fmt='%.3f')

        if node_elec is not None: # then we need to overwrite it
            fparam['node_elec'] = node_elec
        fparam['num_regions'] = 0
        fparam['res0File'] = 'resistivity.dat'
        
        write2in(fparam, fwdDir, typ=self.typ)

        # write the protocol.dat based on measured sequence
        seq = self.surveys[0].df[['a','b','m','n']].values
        if len(self.surveys) > 0: # multiple survey here
            seq = []
            for s in self.surveys:
                seq.append(s.df[['a','b','m','n']].values)
            seq = np.vstack(seq).astype(str)
            seq = np.unique(seq, axis=0)
        protocol = pd.DataFrame(np.c_[1+np.arange(seq.shape[0]),seq],
                                columns=['index','a','b','m','n'])
        # seems to be a problem here if the first electrode is not one. #### TODO #### 
        if (self.typ == 'R3t') | (self.typ == 'cR3t'): # it's a 3D survey
            if len(protocol['a'].values[0].split()) == 1: # we don't have string number
                for c in ['a','b','m','n']: 
                    protocol[c] = '1 ' + protocol[c]

        outputname = os.path.join(fwdDir, 'protocol.dat')
        with open(outputname, 'w') as f:
            f.write(str(len(protocol)) + '\n')
        with open(outputname, 'a') as f:
            protocol.to_csv(f, sep='\t', header=False, index=False)
        
        # run the inversion
        dump('Running R* code (forward mode)...\n')
        self.runR2(fwdDir,dump) # run R*.exe in inside working directory 
        
        if halfspace: 
            # assume that apparent resistivities should equal 100 in the case of a halfspace 
            ap = np.genfromtxt(os.path.join(fwdDir, self.typ + '_forward.dat'), skip_header=1)[:,appResCol]
            modErrRel = np.abs(refRes-ap)/100 
        else: # need to run code again with reference mesh 
            # grab reference transfer resistances made on refined mesh 
            tr0 = np.genfromtxt(os.path.join(fwdDir, self.typ + '_forward.dat'), skip_header=1)[:,trCol]
            
            # keep files for refined mesh 
            shutil.move(os.path.join(fwdDir,mname),
                        os.path.join(fwdDir,mname.replace('.dat','_refined.dat')))
            shutil.move(os.path.join(fwdDir,self.typ + '_forward.dat'),
                        os.path.join(fwdDir,self.typ + '_forward_refined.dat'))
            shutil.move(os.path.join(fwdDir,'resistivity.dat'),
                        os.path.join(fwdDir,'resistivity_refined.dat'))

            # write new mesh that is used for the actual inversion
            if ndim == 2: 
                self.mesh.dat(meshPath) # (use class mesh not the refined one)
            else: # use advanced mesh format if 3D   
                self.mesh.datAdv(meshPath, iadvanced=self.iadvanced)   
            # write a new starting resistivity file
            resFile = np.zeros((self.mesh.numel, ndim+ncol)) # centroid x, y, z, res0, (phase)
            resFile[:,ndim] = refRes
            np.savetxt(os.path.join(fwdDir, 'resistivity.dat'), 
                       resFile,
                       fmt='%.3f')
            
            dump('Running R* code for reference (unrefined) mesh ...\n')
            self.runR2(fwdDir,dump)
            # grab transfer resistances for the reference mesh 
            tr1 = np.genfromtxt(os.path.join(fwdDir, self.typ + '_forward.dat'), skip_header=1)[:,trCol]
            
            k = refRes/tr0 # compute geometric factor using refined mesh transfer resistances 
            ap0 = tr0*k # apparent resistivities for refined mesh 
            ap1 = tr1*k # apparent resistivities for unrefined mesh (use same k as before)
            modErrRel = np.abs(ap0-ap1)/100 # compute relative modelling error 
            
        # append modelling error to survey dataframes 
        # NOTE: In order to compute absolute modelling error multiply it by survey transfer resistances  
        # bug fixed on 18/10/2024
        dump('Adding to modelling error to data structures...\n')
        dferr = pd.DataFrame(seq, columns=['a','b','m','n'])
        dferr['modErrRel'] = modErrRel 
        
        if ipflag :
            simPhase = np.genfromtxt(os.path.join(fwdDir, self.typ + '_forward.dat'), skip_header=1)[:,phaseCol]
            # todo: how to compute phase modelling errors???
            dferr['modErrPhaseRel'] = np.abs(refPhs-simPhase)/100
            
        
        for s in self.surveys:
            if 'modErrRel' in s.df.columns:
                s.df.drop('modErrRel', axis=1)
            if 'modErr' in s.df.columns:
                s.df.drop('modErr', axis=1)
            if 'modErrPhaseRel' in s.df.columns:
                s.df.drop('modErrPhaseRel', axis=1)
            if 'modErrPhase' in s.df.columns:
                s.df.drop('modErrPhase', axis=1)
            s.df = pd.merge(s.df, dferr, on=['a','b','m','n'], how='inner')
            s.df['modErr'] = (s.df['modErrRel']*s.df['resist']).abs() 
            if not ipflag: # continue in the case of resistivity surveys 
                continue 
            # otherwise compute phase related errors ... 
        
        if rmTree:# eventually delete the directory to spare space
            dump('Removing temporary error directory...\n')
            shutil.rmtree(fwdDir)

        self.fwdErrModel = True # class now has a forward error model.
        dump('Modelling error done!\n')

    def showIter(self, index=-2, ax=None, modelDOI=False, cropMaxDepth=False):
        """Dispay temporary inverted section after each iteration.

        Parameters
        ----------
        index : int, optional
            Iteration number to show.
        ax : matplotib axis, optional
            If specified, the graph will be plotted along `ax`.
        modelDOI : bool, optional
            As modelDOI() is always computed using R2 (not cR2), this tells the
            method to look for an R2 looking iteration file.
        cropMaxDepth : bool, optional
            if True, below max depth will be cropped
        """
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        files = os.listdir(self.dirname)
        fs = []
        for f in files:
            if (f[-8:] == '_res.dat') & ((len(f) == 16) | (len(f) == 12)):
                fs.append(f)

        fs = sorted(fs)
        if len(fs) > 1: # the last file is always open and not filled with data
            x = pd.read_csv(os.path.join(self.dirname, fs[index]), header=None, delim_whitespace=True).values
            if self.typ[-1]=='t' and self.elec['buried'].sum()==0: # if 3D get the iteration mesh - and no buried electrodes [HOT FIX] #TODO: extract 3D xbh top surface
                try: # let it pass if surface extraction failed for showing iteration plot
                    if self.surfaceIdx is None: # no surface index
                        self.iterMesh = mt.readMesh(os.path.join(self.dirname, fs[index]).replace('.dat','.vtk')) 
                        _, self.surfaceIdx = self.iterMesh.extractSurface(return_idx=True) # get surface index 
                    x = x[self.surfaceIdx,:]
                except:
                    pass
            if x.shape[0] > 0:
                try: # sometimes i have this part to fail, not sure why
                    triang = tri.Triangulation(x[:,0], x[:,1])
                except ValueError:
                    warnings.warn('Failed to plot iteration mesh')
                    return # if the triangulation fails then fall out of function 
                
                if self.typ[0] == 'c' and modelDOI is False:
                    z = x[:,4]
                else:
                    z = x[:,3] # modelDOI is always computed with R2 not cR2

                cax = ax.tricontourf(triang, z, extend='both')
                if (self.typ == 'R2') or (self.typ == 'cR2'):
                    print('topot', self.topo)
                    print('elec', self.elec['buried'].sum())
                    if (self.topo.shape[0] != 0) | (self.elec['buried'].sum() < self.elec.shape[0]):
                        # just checking we don't have a whole space mesh
                        self._clipContour(ax, cax, cropMaxDepth)
                fig.colorbar(cax, ax=ax, label=r'$\log_{10}\rho$ [$\Omega$.m]')
                elec = self.elec[~self.elec['remote']][['x','y','z']].values
                if self.typ[-1] == 't': # adjust for 3D 
                    yelec = elec[:,1]
                    ylabel = 'Distance [m]'
                    ylim = [np.min(elec[:,1]), np.max(elec[:,1])]
                else:
                    yelec = elec[:,2]
                    ylabel = 'Elevation [m]'
                    ylim = self.zlim
                    
                elecColor = 'ko' if self.darkMode is False else 'wo'
                ax.plot(elec[:,0], yelec, elecColor, markersize=4)
                ax.set_aspect('equal')
                ax.set_xlabel('Distance [m]')
                ax.set_ylabel(ylabel)
                ax.set_xlim([np.min(elec[:,0]), np.max(elec[:,0])])
                ax.set_ylim(ylim)



    def saveInvPlots(self, outputdir=None, **kwargs):
        """Save all plots to output (or working directory). Parameters
        are passed to the `showResults()` method.

        Parameters
        ----------
        outputdir : str, optional
            Path of the output directory. Default is the working directory.
        """
        if outputdir is None:
            outputdir = self.dirname
        if len(self.meshResults) == 0:
            self.getResults()
        for i in range(len(self.meshResults)):
            fig, ax = plt.subplots()
            self.showResults(index=i, ax=ax, **kwargs)
            fname = self.meshResults[i].mesh_title
            fig.savefig(os.path.join(outputdir, fname + '.png'))
        
    
    def getInvError(self):
        a = 1 if self.iTimeLapse else 0
        try:
            if self.typ == 'cR2' or self.typ == 'R2':
                dfs = []
                if self.iTimeLapse:
                    fname = os.path.join(self.dirname, 'ref/f001_err.dat')
                    if os.path.exists(fname):                    
                        df = pd.read_csv(fname, delim_whitespace=True)
                        dfs.append(df)
                for i in range(len(self.surveys)-a):
                    fname = os.path.join(self.dirname, 'f{:03.0f}_err.dat'.format(i+1))
                    if os.path.exists(fname):
                        df = pd.read_csv(fname, delim_whitespace=True)
                        dfs.append(df)
            elif self.typ == 'R3t':
                dfs = []
                if self.iTimeLapse:
                    fname = os.path.join(self.dirname, 'ref/f001_err.dat')
                    if os.path.exists(fname):
                        err = np.genfromtxt(fname, skip_header=1)
                        df = pd.DataFrame(err[:,[0,1,2,3,4,5,6,7,8]],
                                          columns=['sa','P+','sb','P-','sm','C+','sn','C-', 'Normalised_Error'])
                        dfs.append(df)
                for i in range(len(self.surveys)-a):
                    fname = os.path.join(self.dirname, 'f{:03.0f}_err.dat'.format(i+1))
                    if os.path.exists(fname):
                        err = np.genfromtxt(fname, skip_header=1)
                        df = pd.DataFrame(err[:,[0,1,2,3,4,5,6,7,8]],
                                          columns=['sa','P+','sb','P-','sm','C+','sn','C-', 'Normalised_Error'])
                        dfs.append(df);
            else: # TODO cR3t header needs to be standardized
                dfs = []
                if self.iTimeLapse:
                    fname = os.path.join(self.dirname, 'ref/f001_err.dat')
                    if os.path.exists(fname):
                        err = np.genfromtxt(fname, skip_header=1)
                        df = pd.DataFrame(err[:,[0,1,2,3,4,5,6,7,8,11,12]],
                                          columns=['sa','P+','sb','P-','sm','C+','sn','C-', 'Normalised_Error', 'Observed_Phase', 'Calculated_Phase'])
                        dfs.append(df)
                for i in range(len(self.surveys)-a):
                    fname = os.path.join(self.dirname, 'f{:03.0f}_err.dat'.format(i+1))
                    if os.path.exists(fname):
                        err = np.genfromtxt(fname, skip_header=1)
                        df = pd.DataFrame(err[:,[0,1,2,3,4,5,6,7,8,11,12]],
                                          columns=['sa','P+','sb','P-','sm','C+','sn','C-', 'Normalised_Error', 'Observed_Phase', 'Calculated_Phase'])
                        dfs.append(df)
                        
            dfs2 = []
            check = True # check if first value of survey frames are in line number and electrode fmt 
            if len(self.surveys[0].df['a'][0].split()) == 1:
                check = False 
            for df in dfs:
                df = df.astype({'P+':int, 'P-':int, 'C+':int, 'C-':int})
                if check and (self.typ == 'R3t' or self.typ == 'cR3t'):
                    df = df.astype({'sa':int, 'sb':int, 'sm':int, 'sn':int})
                    df['P+'] = df['sa'].astype(str) + ' ' + df['P+'].astype(str)
                    df['P-'] = df['sb'].astype(str) + ' ' + df['P-'].astype(str)
                    df['C+'] = df['sm'].astype(str) + ' ' + df['C+'].astype(str)
                    df['C-'] = df['sn'].astype(str) + ' ' + df['C-'].astype(str)
                else:
                    df = df.astype({'P+':str, 'P-':str, 'C+':str, 'C-':str})
                dfs2.append(df)
            dfs = dfs2
        except Exception as e:
            return # this code is error prone (mainly to empty dataframe error)
        # merge the columns to each survey dataframe
        if  np.sum([df.shape[0] > 0 for df in dfs]) != len(self.surveys):
            print('error in reading error files (do not exists or empty')
            return # this check the number of dfs AND the fact that they are not empty
        
        normerr = np.array([])
        for s, df in zip(self.surveys, dfs):
            normerr = np.vstack(df['Normalised_Error'].values)
            df = df.rename(columns=dict(zip(['P+','P-','C+','C-', 'Normalised_Error'], ['a','b','m','n', 'resInvError'])))
            cols = ['a','b','m','n','resInvError']
            
            if (self.typ == 'cR2') | (self.typ == 'cR3t'):
                df['phaseInvMisfit'] = np.abs(df['Observed_Phase'] - df['Calculated_Phase'])
                cols += ['phaseInvMisfit']
            if 'resInvError' in s.df.columns:
                s.df = s.df.drop('resInvError', axis=1)
            if 'phaseInvMisfit' in s.df.columns:
                s.df = s.df.drop('phaseInvMisfit', axis=1)
            s.df = pd.merge(s.df, df[cols], on=['a','b','m','n'], how='left')
            s.dfInvErrOutputOrigin = s.df.copy() # for being able to reset post processing filters
            s.isequenceInvErrOutputReset = s.isequence.copy()
        self.pinfo['Median normalised error'] = float(np.median(np.abs(normerr)))

                    

    def showPseudoInvError(self, index=0, ax=None, vmin=None, vmax=None, elec=True):
        """Plot pseudo section of errors from file `f001_err.dat`.

        Parameters
        ----------
        index : int, optional
            Index of the survey (if time-lapse or batch). Default `index == 0`.
        ax : matplotlib axis
            If specified, the graph will be plotted against `ax`.
        vmin : float, optional
            Min value.
        vmax : float, optional
            Max value.
        elec : bool, optional
            If `True`, the electrodes are displayed and can be used for filtering.
        """
        self.surveys[index].filterManual(attr='resInvError', vmin=vmin, vmax=vmax,
                    ax=ax, log=False, label='Normalised Error', elec=elec, darkMode=self.darkMode)



    def showPseudoInvErrorIP(self, index=0, ax=None, vmin=None, vmax=None):
        """Display normalized phase error.

        Parameters
        ----------
        index : int, optional
            Index of the survey (if time-lapse or batch). Default `index == 0`.
        ax : matplotlib axis
            If specified, the graph will be plotted against `ax`.
        vmin : float, optional
            Min value.
        vmax : float, optional
            Max value.
        """
        self.surveys[index].filterManual(attr='phaseInvMisfit', vmin=vmin, vmax=vmax,
                    ax=ax, log=False, label='Phase misfit [mrad]', darkMode=self.darkMode)
        

    def showInvError(self, index=0, ax=None):
        """Display inversion error by measurment numbers.
        
        Parameters
        ----------
        index : int, optional
            Index of survey (if time-lapse or batch). Default `index == 0`.
        ax : matplotlib axis
            If provided, the graph will be plotted against this axis.
        """
        self.surveys[index].showInvError(ax=ax)


    def filterInvError(self, index=-1, vmin=None, vmax=None):
        """Remove measurements where inversion error is outside of the defined range.

        Parameters
        ----------
        vmin : float, optional
            minimum value of normalized error below which data will be discarded.
        vmax : float, optional
            maximum value of normalized error above which data will be discarded.
        index : int, optional
            Index of the survey to process. If `index == -1` (default) then the
            processing is applied on all survey independantly.
        """
        if index == -1:
            for s in self.surveys:
                s.filterInvError(vmin=vmin, vmax=vmax)
        else:
            self.surveys[index].filterInvError(vmin=vmin, vmax=vmax) 
        
        
    def saveMeshVtk(self, outputname=None):
        """Save mesh as .vtk to be viewed in paraview.

        Parameters
        ----------
        outputname : str, optional
            Output path of the .vtk produced. By default the mesh is saved in
            the working directory `self.dirname` as `mesh.vtk`.
        """
        warnings.warn('The function is deprecated, use saveMesh() instead.',
                      DeprecationWarning)
        if outputname is None:
            outputname = os.path.join(self.dirname, 'mesh.vtk')
        self.mesh.vtk(outputname)



    def saveMesh(self, outputname=None):
        """Save mesh as .vtk to be viewed in paraview.

        Parameters
        ----------
        outputname : str, optional
            Output path with extension. Available mesh format are:
                - .vtk (Paraview)
                - .node (Tetgen)
                - .dat (R* codes)
            If not provided the mesh is saved in the working directory
            as mesh.vtk.
        """
        if outputname is None:
            self.mesh.vtk(os.path.join(self.dirname, 'mesh.vtk'))
        else:
            if outputname.lower()[-4:] == '.vtk':
                self.mesh.vtk(outputname)
            elif outputname.lower()[-5:] == '.node':
                self.mesh.exportTetgenMesh(prefix=outputname.replace('.node',''))
            elif outputname.lower()[-4:] == '.dat':
                self.mesh.dat(outputname)
            else:
                raise ValueError('mesh export format not recognized. Try either .vtk, .node or .dat.')
                
                
    def exportMesh(self, outputname=None, voxel=False, _forceLocal=False):
        """Export mesh as a different file format, with coordinate conversion
        if Project.coordLocal set to True. 

        Parameters
        ----------
        outputname : str, optional
            Output path with extension. Available mesh format are:
                - .vtk (Paraview)
                - .node (Tetgen)
                - .dat (R* codes)
            If not provided the mesh is saved in the working directory
            as mesh.vtk.
        voxel: bool, optional 
            If True, mesh will be converted to voxel format. 
        _forceLocal: bool, optional
            Forces outputs in the local grid format is self.coordLocal == True. 
        """
        if outputname is None: 
            outputname = 'mesh.vtk'
            
        if isinstance(outputname,str):
            # check for extension 
            if len(outputname.split('.'))<2: 
                outputname = outputname + '.vtk'
                
        # force local grid output in the UI
        if _forceLocal and self.coordLocal:
            self.mesh.exportMesh(outputname,False,self.coordParam, voxel)
            return 
            
        self.mesh.exportMesh(outputname, None, self.coordLocal,self.coordParam, voxel)
        
    def exportElec(self, outputname=None, _forceLocal=False):
        """Export electrodes as a different file format, with coordinate conversion
        if Project.coordLocal set to True. 

        Parameters
        ----------
        outputname : str, optional
            Output path with extension. Available mesh format are:
                - .vtk (Paraview)
                - .vtp (Paraview)
                - .dat (R* codes)
                - .csv (general)
            If not provided the electrode coordinates are saved in the working directory
            as electrodes.vtk.
        _forceLocal: bool, optional
            Forces outputs in the local grid format is self.coordLocal == True. 
        """
        if self.elec is None: # cant run function with electrodes object being populated 
            warnings.warn('No electrode dataframe has been set, cannot export electrodes!')
            return # so return 
        
        if outputname is None: #' set default file name for output 
            outputname = 'electrodes.vtk'
            
        if not isinstance(outputname, str):
            raise Exception('Output name must be a string!')
        
        if len(outputname.split('.')) == 1: 
            outputname += '.vtk'

        eleccopy = self.elec.copy() 
        
        # force local grid output in the UI
        if not _forceLocal and self.coordLocal:
            localx = eleccopy.x.values
            localy = eleccopy.y.values 
            eleccopy['x'], eleccopy['y'] = invRotGridData(
                localx, localy, 
                self.coordParam['x0'],
                self.coordParam['y0'],
                self.coordParam['a']
                )
            
        if outputname.endswith('.csv'):
            eleccopy.to_csv(outputname,index=False)
        elif outputname.endswith('.dat'):
            eleccopy.to_csv(outputname,index=False,sep='\t')
        elif outputname.endswith('.vtk'):
            mt.points2vtk(eleccopy.x, eleccopy.y, eleccopy.z, outputname, 'Electrodes')
        elif outputname.endswith('.vtp'): # todo, add this in meshTools 
            mt.points2vtp(eleccopy.x, eleccopy.y, eleccopy.z, outputname, 'Electrodes')
        else:
            raise Exception('Unrecognised output extension type')
            

    def _toParaview(self, fname,  paraview_loc=None): # pragma: no cover
        """Open file in paraview (might not work if paraview is not in the PATH,
        in this case, pass parview location as `paraview_loc`).

        Parameters
        ----------
        fname : str
            Path of the .vtk file to be opened.
        paraview_loc: str, optional
            **Windows ONLY** maps to the executable paraview.exe. The program
            will attempt to find the location of the paraview install if not given.
        """

        if OS == "Windows":
            if paraview_loc is None:
                found,cmd_line = self.mesh.findParaview()
                if not found:
                    print("Cannot find paraview location")
                    return
                cmd_line = '"' + cmd_line + '" ' + fname
            elif isinstance(paraview_loc,str):
                cmd_line = '"' + paraview_loc + '" ' + fname
            else:
                print("Cannot find where paraview is installed")
                return
        else:
            cmd_line = 'paraview ' + fname

        try:#try and launch paraview
            #Popen([cmd_line, os.path.join(self.dirname, fname)])
            os.popen(cmd_line)
        except PermissionError:
            print("Your operating system has blocked launching Paraview")
            #windows doesnt like calling paraview from python for some reason
            #will need to look into this further.


    def showMeshInParaview(self, paraview_loc=None): # pragma: no cover
        """Show the mesh in paraview (mostly useful for 3D surveys.

        Parameters
        ----------
        paraview_loc: str, optional
            **Windows ONLY** maps to the excuatable paraview.exe. The program
            will attempt to find the location of the paraview install if not given.
        """
        print('Saving mesh as vtk...', end='')
        self.saveMeshVtk() # save in default dirname
        print('done.\n Launching paraview.')
        self._toParaview(os.path.join(self.dirname, 'mesh.vtk'),
                         paraview_loc=paraview_loc)


    def showInParaview(self, index=0, paraview_loc=None): # pragma: no cover
        """Open paraview to display the .vtk file.

        Parameters
        ----------
        index: int, optional
            Timestep to be shown in paraview (for an individual survey this 1).
        paraview_loc: str, optional
            **Windows ONLY** maps to the excuatable paraview.exe. The program
            will attempt to find the location of the paraview install if not given.
        """
        fname = 'f{:03d}_res.vtk'.format(index+1)
        self._toParaview(os.path.join(self.dirname, fname), paraview_loc=paraview_loc)


    def showSlice(self, index=0, ax=None, attr=None, axis='z', vmin=None, vmax=None):
        """Show slice of 3D mesh interactively.
        
        Parameters
        ----------
        index : int, optional
            Index of the survey. Default is first survey (index == 0).
        ax : matplotlib.Axes, optional
            Axis on which to plot the graph.
        attr : str, optional
            Attribute to plot. Default is 'Resistivity(ohm.m)'.
        axis : str, optional
            Either 'x', 'y', or 'z' (default).
        vmin : float, optional
            Minimum value for colorbar.
        vmax : float, optional
            Maximum value for colorbar.
        """
        if attr is None:
            attr = list(self.meshResults[index].df.keys())[0]
        self.meshResults[index].showSlice(
                attr=attr, axis=axis, ax=ax, vmin=vmin, vmax=vmax)

        
    ## Sorting electrode numbers ##
    # def shuntIndexes(self):
    #     """Shunt electrode indexes to start at 1.
    #     """
    #     debug=True
    #     if len(self.surveys)>1:
    #         debug=False
    #     for i in range(len(self.surveys)):
    #         self.surveys[i].shuntIndexes(debug=debug)

    def normElecIdx(self): # pragma: no cover
        """Normalise electrode indexes to start at 1 in consective and ascending order.
        """
        debug = True
        if len(self.surveys)>1:
            debug=False
        for i in range(len(self.surveys)):
            self.surveys[i].normElecIdx(debug=debug)

    ## make 3d coordinates for a 2d line in to 2d ##
    def elec2distance(self, yDominant=False, iMoveElec=False):
        """
        Convert 3D electrode XY coordinates into just X coordinates. Use for
        2D lines only!
        If self.elec has been set then each survey will use the electrodes set
        in the R2 master class. If not then the R2 master class will take on the
        elec values set for the first survey in a sequence.

        Parameters
        ----------
        yDominant: bool, optional
            If electrodes are prodimently spaced in the y direction then set yDominant
            to True.
        iMoveElec: bool, optional
            If moving electrodes are present then set to True, so that the same
            electrode positions are not given to each survey.
        """
        if self.typ == 'R3t' or self.typ == 'cR3t':
            raise ValueError("Cannot compress 3D survey coordinates to 2D for a 3D survey type.")

        if not iMoveElec:#check if elec has been assigned already
            try: # if electrodes are set in the R2 class then use these for each survey
                for i in range(len(self.surveys)):
                    self.surveys[i].elec = self.elec
            except AttributeError: #if not already set then assume the electrodes are set for each survey
                pass

        if yDominant: # swap x and y around in raw coordinates
            for i in range(len(self.surveys)):
                elec = self.surveys[i].elec[['x','y','z']].values.copy()
                x = elec[:,0]
                y = elec[:,1]
                self.surveys[i].elec['x'] = y
                self.surveys[i].elec['y'] = x

        for i in range(len(self.surveys)):
            self.surveys[i].elec2distance() # go through each survey and compute electrode
        self.elec = None
        self.setElec(self.surveys[0].elec)

    @staticmethod
    def topo2distance(x,y,z): # pragma: no cover
        """Convert 3d xy data in pure x lateral distance.
        Use for 2D data only!
        """
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        
        idx = np.argsort(x) # sort by x axis
        x_sorted = x[idx]
        y_sorted = y[idx]
        z_sorted = z[idx]
        
        x_abs = np.zeros_like(x, dtype=float)
        #the first entry should be x = 0
        for i in range(len(x)-1):
            delta_x = x_sorted[i] - x_sorted[i+1]
            delta_y = y_sorted[i] - y_sorted[i+1]
            sq_dist = delta_x**2 + delta_y**2
            x_abs[i+1] =  x_abs[i] + np.sqrt(sq_dist)
        
        # return values in the order in which they came
        new_topo = np.zeros((len(x),3), dtype=float)        
        for i in range(len(z)):
            put_back = idx[i] # index to where to the value back again
            new_topo[put_back,0] = x_abs[i]
            new_topo[put_back,2] =  z_sorted[i]
    
        return new_topo

# WIP
#    def timelapseErrorModel(self, ax=None):
#        """Fit an power law to time-lapse datasets.
#
#        Parameters
#        ----------
#        ax : matplotlib axis, optional
#            If specified, graph will be plotted on the given axis.
#
#        Returns
#        -------
#        fig : matplotlib figure, optional
#            If ax is not specified, the function will return a figure object.
#        """
#        if ax is None:
#            fig, ax = plt.subplots()
#        numbins = 20
#
#        if 'recipMean' not in self.df.columns:
#            self.computeReciprocal()
#        dfg = self.df[self.df['irecip'] > 0]
#        binsize = int(len(dfg['recipMean'])/numbins)
#        error_input = np.abs(dfg[['recipMean', 'recipError']]).sort_values(by='recipMean').reset_index(drop=True) # Sorting data based on R_avg
#        bins = np.zeros((numbins,2))
#        for i in range(numbins): # bining
#            ns=i*binsize
#            ne=ns+binsize-1
#            bins[i,0] = error_input['recipMean'].iloc[ns:ne].mean()
#            bins[i,1] = error_input['recipError'].iloc[ns:ne].mean()
##        print(bins)
##        print(np.sum(np.isnan(bins)))
##        print(np.sum(np.isinf(bins)))
##        coefs= np.linalg.lstsq(np.vstack([np.ones(len(bins[:,0])), np.log(bins[:,0])]).T, np.log(bins[:,1]), rcond=None)[0] # calculating fitting coefficients (a,m)
#        coefs = np.polyfit(np.log(bins[:,0]), np.log(bins[:,1]), 1)[::-1] #order is of coefs is opposite to lstqd
#        R_error_predict = np.exp(coefs[0])*(bins[:,0]**coefs[1]) # error prediction based of power law model
#        ax.plot(np.abs(dfg['recipMean']),np.abs(dfg['recipError']), '+', label = "Raw")
#        ax.plot(bins[:,0],bins[:,1],'o',label="Bin Means")
#        ax.plot(bins[:,0],R_error_predict,'r', label="Power Law Fit")
#        ax.set_xscale('log')
#        ax.set_yscale('log')
#        # lines above are work around to https://github.com/matplotlib/matplotlib/issues/5541/
#        ax.set_ylabel(r'$R_{error} [\Omega]$')
#        ax.set_xlabel(r'$R_{avg} [\Omega]$')
#        ax.legend(loc='best', frameon=True)
#        R2= self.R_sqr(np.log(bins[:,1]),np.log(R_error_predict))
#        a1 = np.exp(coefs[0])
#        a2 = coefs[1]
##        a3 = np.exp(coefs[0])
##        a4 = coefs[1]
#        print('Error model is R_err = {:.2f} R_avg^{:.3f} (R^2 = {:.4f})'.format(a1,a2,R2))
#        if a1 > 0.001:
#            ax.set_title('Multi bin power-law resistance error plot\n' + r'$R_{{error}}$ = {:.3f}$R_{{avg}}^{{{:.3f}}}$ (R$^2$ = {:.3f})'.format(a1,a2,R2))
#        else:
#            ax.set_title('Multi bin power-law resistance error plot\n' + r'$R_{{error}}$ = {:.2e}$R_{{avg}}^{{{:.3e}}}$ (R$^2$ = {:.3f})'.format(a1,a2,R2))
#        self.df['resError'] = a1*(np.abs(self.df['recipMean'])**a2)
#        def errorModel(df):
#            x = df['recipMean'].values
#            return a1*(np.abs(x)**a2)
#        self.errorModel = errorModel
##        self.errorModel = lambda x : a1*(np.abs(x)**a2)
#        if ax is None:
#            return fig



    # def computeCond(self): # automatically done in getResults()
    #     """Compute conductivities from resistivities for the ERT mesh
    #     """
    #     if self.typ=='R3t' or self.typ=='cR3t':
    #         res_name = 'Resistivity'
    #     else:
    #         res_name = 'Resistivity(Ohm-m)'
    #     for i in range(len(self.meshResults)):
    #         self.meshResults[i].computeReciprocal(res_name,'Conductivity(S/m)')


    def postProcTl(self):
        """
        Post processing for time-lapse surveys. 
        
        Compute the absolute and the relative difference in resistivity
        between inverted surveys.
        """
        if not self.iTimeLapse:
            raise Exception("Difference calculation only available for time-lapse surveys.")
        if len(self.meshResults) == 0:
            self.getResults()

        # create an index for the values inside of the zone of interest
        # needed as the reference survey is not cropped by default
        inside = np.ones(self.meshResults[0].numel, dtype=bool)
        if (self.typ == 'R3t') or (self.typ == 'cR3t'):
            pname = 'num_xy_poly'
            tableName = 'xy_poly_table'
        else:
            pname = 'num_xz_poly'
            tableName = 'xz_poly_table'
        if self.param[pname] > 0:
            meshx = np.array(self.meshResults[0].elmCentre[:,0])
            meshy = np.array(self.meshResults[0].elmCentre[:,1])
            meshz = np.array(self.meshResults[0].elmCentre[:,2])
            # poly = (self.param['xz_poly_table'][:,0],
                    # self.param['xz_poly_table'][:,1])
            path = mpath.Path(self.param[tableName])

            if self.typ[-2]=='3':
                # inside1 = iip.isinpolygon(meshx, meshy, poly)
                inside1 = path.contains_points(np.c_[meshx, meshy])
                inside2 = (meshz > self.param['zmin']) & (meshz < self.param['zmax'])
                inside = inside1 & inside2
            else:
                # inside = iip.isinpolygon(meshx, meshz, poly)
                inside = path.contains_points(np.c_[meshx, meshz])
                
        # compute absolute and relative difference in resistivity
        res_names = np.array(['Resistivity','Resistivity(Ohm-m)','Resistivity(ohm.m)'])
        res_name = res_names[np.in1d(res_names, list(self.meshResults[0].df.keys()))][0]
        res0 = np.array(self.meshResults[0].df[res_name])[inside]
        
        # nb: the first survey will have 0 % change
        self.meshResults[0].addAttribute(np.zeros(self.meshResults[0].numel), 'difference(percent)') 
        for i in range(1, len(self.meshResults)):
            if 'difference(percent)' not in self.meshResults[i].df.columns:
                try:
                    res = np.array(self.meshResults[i].df[res_name])
                    self.meshResults[i].addAttribute(res - res0, 'diff(Resistivity)')
                    self.meshResults[i].addAttribute((res-res0)/res0*100, 'difference(percent)')
                except Exception as e:
                    print('error in computing difference:', e)
                    pass
                
        # crop the mesh0 object
        if self.iTimeLapse:
            self.meshResults[0] = self.meshResults[0].crop(self.param['xz_poly_table'])
        

    def saveVtks(self, dirname=None):
        """Save vtk files of inversion results to a specified directory.

        Parameters
        ----------
        dirname: str
            Directory in which results will be saved. Default is the working directory.
        """
        if dirname is None:
            dirname = self.dirname
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        amtContent = ''
        if len(self.meshResults) == 0:
            self.getResults()
        count=0
        
        # check param for truncation variables 
        for key in ['num_xz_poly','num_xy_poly']:
            if key not in self.param.keys():
                self.param[key] = 0 
        for key in ['zmin','zmax']:
            if key not in self.param.keys():
                self.param[key] = None 
            
        # loop through and export meshes as vtk files 
        for mesh, s in zip(self.meshResults, self.surveys):
            count+=1
            meshcopy = mesh.copy()
            if self.iTimeLapse:
                fname = 'time_step{:0>4}.vtk'.format(count)
            else:
                fname = mesh.mesh_title + '.vtk'
            file_path = os.path.join(dirname, fname) 
            if self.trapeziod is not None and self.pseudo3DMeshResult is None:
                meshcopy = meshcopy.crop(self.trapeziod)
            elif '3' in self.typ and self.param['num_xy_poly'] > 2: 
                meshcopy = meshcopy.crop(self.param['xy_poly_table'])
                meshcopy.elec = None 
                if self.param['zmin'] is None or self.param['zmax'] is None:
                    pass # cant truncate mesh with out z limits 
                else: 
                    meshcopy = meshcopy.truncateMesh(zlim=[self.param['zmin'],
                                                           self.param['zmax']])
            elif self.pseudo3DMeshResult is not None and self.projs[count-1].trapeziod is not None:
                meshcopy = meshcopy.crop(self.projs[count-1].trapeziod)
                
            meshcopy.vtk(file_path, title=mesh.mesh_title) # save to vtk 
            amtContent += "\tannotations.append('%s')\n"%mesh.mesh_title
            if self.pseudo3DMeshResultList is not None:
                file_path = os.path.join(dirname, mesh.mesh_title + '_3D.vtk')
                self.pseudo3DMeshResultList[count-1].vtk(file_path, title=mesh.mesh_title)
        if self.pseudo3DMeshResult is not None: 
            self.pseudo3DMeshResult.mesh_title = 'Pseudo_3D_result'
            file_path = os.path.join(dirname, self.pseudo3DMeshResult.mesh_title + '.vtk')
            self.pseudo3DMeshResult.vtk(file_path, title='Pseudo_3D_result')
            amtContent += "\tannotations.append('%s')\n"%self.pseudo3DMeshResult.mesh_title

        fh = open(os.path.join(dirname,'amt_track.py'),'w')
        fh.write('from paraview.simple import *\n')
        fh.write('def start_cue(self):\n')
        fh.write('\tglobal annotations\n')
        fh.write('\tglobal maxIndex\n')
        fh.write("\ttextSource = paraview.simple.FindSource('Text1')\n")
        fh.write("\tif textSource is None:\n")
        fh.write("\t\tText()\n")
        fh.write("\tannotations = []\n")
        fh.write(amtContent)
        fh.write('\tmaxIndex=len(annotations)\n')
        fh.write('def tick(self):\n')
        fh.write('\tglobal annotations\n')
        fh.write('\tglobal maxIndex\n')
        fh.write('\tindex = int( self.GetClockTime() )\n')
        fh.write('\tif index >= maxIndex :\n')
        fh.write('\t\tindex = maxIndex - 1\n')
        fh.write("\ttextSource = paraview.simple.FindSource('Text1')\n")
        fh.write('\ttextSource.Text = annotations[index]\n')
        fh.write('def end_cue(self): pass\n')
        fh.close()

    def exportMeshResults(self, dirname=None, ftype='.vtu', voxel=False, 
                          _forceLocal=False):
        """Save mesh files of inversion results to a specified directory. If 
        results are in a local grid, they will be converted back into thier 
        respective utm or nationa grid coordinates. 

        Parameters
        ----------
        dirname: str
            Directory in which results will be saved. Default is the working directory.
        ftype: str
            Type of file extension. 
        voxel: bool, optional
            Force export to be in a voxel format. 
        _forceLocal: bool, optional
            Force output to be in a local grid if self.coordLocal is True. Meant 
            for use with the UI. 
        """
            
        if dirname is None:
            dirname = self.dirname
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        amtContent = ''
        if len(self.meshResults) == 0:
            self.getResults()
        count=0
        ext = '.vtk'
        if 'vtk' in ftype.lower():
            ext = '.vtk'
        elif 'vts' in ftype.lower():
            ext = '.vts'
        elif 'vtu' in ftype.lower():
            ext = '.vtu'
        elif 'dat' in ftype.lower(): 
            ext = '.dat'
        elif 'node' in ftype.lower(): 
            ext = '.node'
        elif 'xyz' in ftype.lower(): 
            ext = '.xyz'
        elif 'csv' in ftype.lower(): 
            ext = '.csv'
        # add file types as appropiate! 
        
        # force local grid output in the UI
        coordLocalFlipped = False 
        if _forceLocal and self.coordLocal:
            self.coordLocal = False 
            coordLocalFlipped = True 
        
        
        # check param for truncation variables 
        for key in ['num_xz_poly','num_xy_poly']:
            if key not in self.param.keys():
                self.param[key] = 0 
        for key in ['zmin','zmax']:
            if key not in self.param.keys():
                self.param[key] = None 
                
        # grab surface of mesh in 3d surveys to better inform voxel mesh creation        
        if '3' in self.typ and voxel:
            surfaceMesh = self.mesh.extractSurface() 
        else:
            surfaceMesh = None 
    
        # loop through and export meshes as vtu (or other) files 
        for mesh, s in zip(self.meshResults, self.surveys):
            count+=1
            meshcopy = mesh.copy()
            
            if self.iTimeLapse and voxel:
                fname = 'time_step{:0>4}'.format(count) + '_voxel' + ext 
            elif self.iTimeLapse:
                fname = 'time_step{:0>4}'.format(count) + ext 
            elif voxel:
                fname = mesh.mesh_title + '_voxel' + ext 
            else:
                fname = mesh.mesh_title + ext 
            file_path = os.path.join(dirname, fname) 
            if self.trapeziod is not None and self.pseudo3DMeshResult is None:
                meshcopy = meshcopy.crop(self.trapeziod)
            elif '3' in self.typ and self.param['num_xy_poly'] > 2: 
                meshcopy = meshcopy.crop(self.param['xy_poly_table'])
                meshcopy.elec = None 
                if self.param['zmin'] is None or self.param['zmax'] is None:
                    pass # cant truncate mesh with out z limits 
                else: 
                    meshcopy = meshcopy.truncateMesh(zlim=[self.param['zmin'],
                                                           self.param['zmax']])
            elif self.pseudo3DMeshResult is not None and self.projs[count-1].trapeziod is not None:
                meshcopy = meshcopy.crop(self.projs[count-1].trapeziod)
                
            # convert to utm if needed (not for psuedo3d, that'll be handled differently)
            if self.coordLocal and self.pseudo3DMeshResult is None: 
                coordLocal = True 
            else:
                coordLocal = False 
            
            # add electrode types 
            meshcopy.elec = s.elec[['x','y','z']].values 
            elec_type =['surface']*len(s.elec)
            s.elec.reset_index(inplace=True, drop=True)
            for i in range(len(s.elec)):
                if s.elec['buried'][i]:
                    elec_type[i] = 'buried'
                elif s.elec['remote'][i]:
                    elec_type[i] = 'remote'
            meshcopy.elec_type = elec_type  
            
            meshcopy.surfaceMesh = surfaceMesh
            meshcopy.exportMesh(file_path, ftype.replace('.',''),
                                coordLocal, self.coordParam, 
                                voxel, self.meshParams) # save to vtk 
            amtContent += "\tannotations.append('%s')\n"%mesh.mesh_title
            if self.pseudo3DMeshResultList is not None:
                file_path = os.path.join(dirname, mesh.mesh_title + '_3D' + ext)
                self.pseudo3DMeshResultList[count-1].exportMesh(file_path, 
                                                              ftype.replace('.',''),
                                                              self.coordLocal,
                                                              self.coordParam,
                                                              False)
                
        if self.pseudo3DMeshResult is not None: 
            self.pseudo3DMeshResult.mesh_title = 'Pseudo_3D_result'
            self.pseudo3DMeshResult.elec = self.pseudo3DSurvey.elec[['x','y','z']].values 
            elec_type =['surface']*len(self.pseudo3DSurvey.elec)
            for i in range(len(self.pseudo3DSurvey.elec)):
                if self.pseudo3DSurvey.elec['buried'][i]:
                    elec_type[i] = 'buried'
                elif self.pseudo3DSurvey.elec['remote'][i]:
                    elec_type[i] = 'remote'
            self.pseudo3DMeshResult.elec_type = elec_type  
            file_path = os.path.join(dirname, self.pseudo3DMeshResult.mesh_title + ext)
            self.pseudo3DMeshResult.exportMesh(file_path, ftype.replace('.',''), 
                                             self.coordLocal, self.coordParam,
                                             voxel)
            amtContent += "\tannotations.append('%s')\n"%self.pseudo3DMeshResult.mesh_title

        if 'v' in ext: 
            # write out paraview animation track 
            fh = open(os.path.join(dirname,'amt_track.py'),'w')
            fh.write('from paraview.simple import *\n')
            fh.write('def start_cue(self):\n')
            fh.write('\tglobal annotations\n')
            fh.write('\tglobal maxIndex\n')
            fh.write("\ttextSource = paraview.simple.FindSource('Text1')\n")
            fh.write("\tif textSource is None:\n")
            fh.write("\t\tText()\n")
            fh.write("\tannotations = []\n")
            fh.write(amtContent)
            fh.write('\tmaxIndex=len(annotations)\n')
            fh.write('def tick(self):\n')
            fh.write('\tglobal annotations\n')
            fh.write('\tglobal maxIndex\n')
            fh.write('\tindex = int( self.GetClockTime() )\n')
            fh.write('\tif index >= maxIndex :\n')
            fh.write('\t\tindex = maxIndex - 1\n')
            fh.write("\ttextSource = paraview.simple.FindSource('Text1')\n")
            fh.write('\ttextSource.Text = annotations[index]\n')
            fh.write('def end_cue(self): pass\n')
            fh.close()
            
        # force local grid output in the UI
        if _forceLocal and coordLocalFlipped:
            self.coordLocal = True 


    def saveCwd(self, outputdir):
        """Save all ouputs (_res.dat, .vtk, ...) from the working directory
        generated during inversion to the designated directory.

        Parameters
        ----------
        outputdir : str
            Path to the directory to save the files.
        """
        wd = os.path.join(outputdir, 'wd')
        if os.path.exists(wd):
            shutil.rmtree(wd)
        shutil.copytree(self.dirname, wd)
        
    def saveData(self, outputdir):
        """
        Saves all data generated by ResIPy in the current working directory 
        to a specified folder. Serves same function as Project.saveCwd()
        but retained for backwards compability. 

        Parameters
        ----------
        outputdir : str
            Path to the directory to save the files.

        """
        self.saveCwd(outputdir)
        
    
    
    
    
    # def saveFilteredData(self, fname, savetyp='Res2DInv (*.dat)'):
    #     """Save filtered data in formats to be used outside ResIPy (e.g. Res2DInv).

    #     Parameters
    #     ----------
    #     fname : str
    #         Path where to save the file.
    #     savetyp : str, optional
    #         Saving format. To be determined in GUI. Default: Res2DInv (.dat)
    #     """
    #     elec = self.elec[['x','y','z']].values
    #     spacing = elec[1,0] - elec[0,0] # TODO (gb) not sure if this is needed
    #     for s, i in zip(self.surveys, range(len(self.surveys))):
    #         df = s.df.query('irecip >=0').copy() # not saving reciprocal data
    #         # if spacing == None:
    #         #     spacing = elec[1,0]-elec[0,0] # for batch surveys the spacing can differ and not follow user input
    #         # else:
    #         #     spacing = spacing
    #         # df[['a','b','m','n']] *= spacing
    #         lookupDict = dict(zip(self.elec['label'], self.elec['x'].values))
    #         data = df[['a','b','m','n']].replace(lookupDict).values
    #         for i, a in enumerate(['a','b','m','n']):
    #             df[a] = data[:,i] # this way should avoid warnings 
    #         # df[['a','b','m','n']] = data
    #         if savetyp == 'Res2DInv (*.dat)':
    #             param = {'num_meas': df.shape[0],
    #                      'lineTitle': self.param['lineTitle'],
    #                      'spacing': spacing}
    #             write2Res2DInv(param, fname, df, elec, self.typ)
    #         elif savetyp == 'Comma Separated Values (*.csv)':
    #             write2csv(fname, df, elec, self.typ)
    #         elif savetyp == 'E4D survey file (*.srv)':
    #             writeSrv(fname, df, elec)

    #         fname = fname[:-4] + str(i) + fname[-4:] # to iterate file numbers in case of timelapse survey

    
    
    def exportData(self, fname, savetyp='Res2DInv (*.dat)', err=False, recip=False):
        """Save filtered data in formats to be used outside ResIPy (e.g. Res2DInv).

        Parameters
        ----------
        fname : str
            Path where to save the file.
        savetyp : str, optional
            Saving format. To be determined in GUI. Default: Res2DInv (.dat).
            Also possible: Protocol (*.dat).
        err : bool, optional
            Flag to include errors. The default is False.
        recip : bool, optional
            Flag to include reciprocals. The default is False.
        """
        flag3d = False 
        if '3' in self.typ:
            flag3d = True
        
        ipFlag = False
        if 'c' in self.typ:
            ipFlag = True
        
        elec = self.elec[['x','y','z']].values
        for s, i in zip(self.surveys, range(len(self.surveys))):
            if recip:
                isubset = [True]*len(s.df)
                df = s.df
            else:
                isubset = None
                df = s.df.query('irecip >=0').copy() # not saving reciprocal data
            
            if savetyp == 'Protocol (*.dat)':
                s.write2protocol(outputname=fname, err=err, isubset=isubset, threed=flag3d, ip=ipFlag)
            else:
                lookupDict = dict(zip(self.elec['label'], np.arange(self.elec.shape[0])))
                data = df[['a','b','m','n']].replace(lookupDict).values
                for i, a in enumerate(['a','b','m','n']):
                    df[a] = data[:,i] 
                df[['a','b','m','n']] = df[['a','b','m','n']] + 1
                if savetyp == 'Res2DInv (*.dat)':
                    if 'lineTitle' not in self.param:
                        self.param['lineTitle'] = 'MyBeautifulLine'
                    spacing = elec[1,0] - elec[0,0]
                    param = {'num_meas': df.shape[0],
                             'lineTitle': self.param['lineTitle'],
                             'spacing': spacing}
                    write2Res2DInv(param, fname, df, elec, self.typ)
                elif savetyp == 'Comma Separated Values (*.csv)':
                    write2csv(fname, df, elec, self.typ)
                elif savetyp == 'E4D survey file (*.srv)':
                    writeSrv(fname, df, elec)

            fname = fname[:-4] + str(i) + fname[-4:] # to iterate file numbers in case of timelapse survey
    
    
    
    
    
    # def exportData(self, outputname=None, ftype='protocol', err=False, 
    #                recip=False):
    #     """
    #     Export preconditioned data used by ResIPy into another format (
    #     different to project.saveData).

    #     Parameters
    #     ----------
    #     outputname : str, optional
    #         Outputname. The default is None. If not given then the function falls
    #         back on the survey name. If set then a number will be appended to the 
    #         file name in the case of time-lapse (or batch) surveys. 
    #     ftype : str, optional
    #         Export File type, choose from either protocol, srv, csv, ResInv. 
    #         The default is 'protocol'.
    #     err : bool, optional
    #         Flag to include errors. The default is False.
    #     recip : bool, optional
    #         Flag to include reciprocals. The default is False.
    #     """
        
    #     # check the ftype 
    #     ext = '.dat'
    #     if ftype == 'protocol':
    #         ext == '.dat'
    #     elif ftype == 'ResInv':
    #         ext = '.dat'
    #     elif ftype == 'srv':
    #         ext = '.srv'
    #     elif ftype == 'csv':
    #         ext == '.csv'
    #     else:
    #         raise Exception('unknown file type')
            
    #     flag3d = False 
    #     if '3' in self.typ:
    #         flag3d = True 
            
    #     for i,s in enumerate(self.surveys): 
    #         if recip:
    #             isubset = [True]*len(s.df)
    #         else:
    #             isubset = None 
                
    #         protocol = s.write2protocol(err=err, isubset=isubset, threed=flag3d)
            
    #         # decide on file outputname 
    #         if outputname is None: 
    #             fout = s.name 
    #             if fout == '':
    #                 fout = 'SurveyData'
    #                 if len(self.surveys) > 1:
    #                     fout = 'SurveyData_{:0>3d}'.format(i)
    #         else:
    #             fout = outputname 
    #             if fout.endswith(ext):
    #                 fout = fout.replace(ext,'')
    #             if len(self.surveys) > 1: 
    #                 fout += '_{:0>3d}'.format(i)
                    
    #         if not fout.endswith(ext):
    #             fout += ext 
                
    #         # write to file 
    #         if ftype == 'protocol':
    #             protocol.to_csv(fout, index=False, sep='\t', lineterminator='\n')
    #         elif ftype == 'ResInv':
    #             param = {'lineTitle':s.name}
    #             write2Res2DInv(param, fout, protocol, 
    #                            self.elec[['x','y','z']].values)
    #         elif ftype =='srv':
    #             writeSrv(fout, protocol, self.elec[['x','y','z']].values)
    #         elif ftype =='csv':
    #             protocol.to_csv(fout, index=False, sep=',', lineterminator='\n')
                
    #     return 
                

    def showParam(self):
        """Print parameters in `R2.param` dictionary.
        """
        [print(key) for i,key in enumerate(self.param)]


    def filterZeroMeasSurveys(self):
        """Filter out badly behaved surveys, where after all other QC no measurements
        are actually left."""
        count=0
        survey_len = [len(self.surveys[i].df) for i in range(len(self.surveys))]
        while min(survey_len)==0:
            survey_len = [len(self.surveys[i].df) for i in range(len(self.surveys))]
            if min(survey_len)==0:
                bad_idx = np.argmin(np.array(survey_len))
                del(self.surveys[bad_idx])
                count += 1
        print("%i surveys removed as they had no measurements!"%count)
        
    def _estimateMemory(self, dump=print, debug=False):
        """More accurate calculation of the amount of memory required
        to run a forward model or inversion. 
        
        Nb: currently only experimental 

        Parameters
        ----------
        dump : function, optional
            stdout direction, ie where to print outputs 
        debug : TYPE, optional
            If true all the variable values are printed to console. The 
            default is False.

        Returns
        -------
        memInv : float
            Memory needed for inverse problem in gigabytes 

        """
        #NB: using variable names from Andy's codes 
        if self.mesh is None:
            dump('A mesh is required before a memory usage estimate can be made')
            return 0
        if len(self.surveys) == 0:
            dump('A survey needs to imported before a memory usage estimate can be made')
            return 0
        else: # number of measurements computation 
            nmeas = []
            for s in self.surveys:
                df = s.df 
                if 'irecip' in df.columns: 
                    ie = df['irecip'].values >= 0 # count the number of measurements actually put to file 
                    nmeas.append(sum(ie))
                else:
                    nmeas.append(len(df))
            num_ind_meas = np.mean(nmeas) # number of measurements 
            mnum_ind_meas = max(nmeas) # maximum number of measurements (more conservative estimate)
                    
        #nsize A estimation - describes number of connected nodes 
        if self.mesh.NsizeA is None:
            self.mesh.computeNconnec()
        nsizeA = self.mesh.NsizeA 
        
        #other mesh parameters 
        numnp = self.mesh.numnp
        numel = self.mesh.numel
        npere = self.mesh.type2VertsNo()
        if 'param' in self.mesh.df.keys():
            num_param = len(np.unique(self.mesh.df['param'].values))
        else:
            num_param = self.mesh.numel
            
        nfaces = self.mesh.type2FaceNo()
        
        #electrodes 
        num_electrodes = len(self.elec)
        
        # special case of pseudo 3D 
        if len(self.projs) > 0: 
            _num_param = [0]*len(self.projs)
            _num_elec = [0]*len(self.projs)
            for i,proj in enumerate(self.projs):
                _num_param[i] = proj.mesh.numel 
                _num_elec[i] = len(proj.elec)
            # convservative estimate, take the maximums! 
            num_param = max(_num_param)
            num_electrodes = max(_num_elec)
        
        # this refers to the roughness matrix in 2D problems
        if 'inverse_type' in self.param and self.param['inverse_type'] == 2:
            numRterm = 13 
        else:
            numRterm = 5 
        
        if '3' in self.typ: # do 3D calculation 
            memDP=numnp*(8+num_electrodes)+nsizeA+numel+mnum_ind_meas * 2
            memR=0
            memI=numnp*2+(npere+2)*numel+numnp+1+nsizeA+num_electrodes*3+mnum_ind_meas*12
            memL=numel*2
            
            memFwd = (memL + memI*4 + memR*4 + memDP*8)/1.0e9
              
            # now estimate for inverse case 
            memDP=memDP+num_param*9+mnum_ind_meas*(num_param+6) 
            memR=memR+(num_param*nfaces)
            memI=memI+num_param*nfaces    
                
        else: # do 2D calculation
            memDP=(numnp)*(5+num_electrodes)+nsizeA+numel+mnum_ind_meas*3+mnum_ind_meas       
            memR=0
            memI=nsizeA+numel*6+numnp*4+mnum_ind_meas*8
            memL=numel+numnp
            
            memFwd = (memL + memI*4 + memR*4 + memDP*8)/1.0e9

            # now estimate for inverse case 
            memDP=memDP+numel+num_param*10+mnum_ind_meas*(num_param+7) 
            memR=memR+num_param*numRterm
            memI=memI+num_param*numRterm
            memL=memL+mnum_ind_meas
        
        memInv = (memL + memI*4 + memR*4 + memDP*8)/1.0e9
        
        dump('ResIPy Estimated RAM usage = %f Gb'%memInv)
        
        self.pinfo['Estimated RAM for forward solution (Gb)'] = float(memFwd)
        self.pinfo['Estimated RAM for inverse solution (Gb)'] = float(memInv) 
            
        
        avialMemory = psutil.virtual_memory().available*9.31e-10
        if memInv >= avialMemory:
            dump('*** It is likely that more RAM is required for inversion! ***\n'
                 '*** Make a coarser mesh ***')
        
        if debug: #print everything out 
            print('numnp = %i'%numnp)
            print('numel = %i'%numel)
            print('nsizA = %i'%nsizeA)
            print('num_param = %i'%num_param)
            print('num_electrodes = %i'%num_electrodes)
            print('num_ind_meas = %i'%num_ind_meas)
            print('npere = %i'%npere)
            print('nfaces = %i'%nfaces)
            print('memDP = %i'%memDP)
            print('memR = %i'%memR)
            print('memI = %i'%memI)
            print('memL = %i'%memL)
            
        return memInv 
    
    def _estimateMemoryJac(self, dump=print):
        """Estimates the memory needed by inversion code to formulate 
        a jacobian matrix

        Parameters
        ----------
        dump : function, optional
            stdout direction, ie where to print outputs 

        Returns
        -------
        Gb : float
            Memory needed for jacobian formulation in gigabytes 

        """
        #NB: using variable names from Andy's codes 
        if self.mesh is None:
            print('A mesh is required before a memory usage estimate can be made')
            return 0
        if len(self.surveys) == 0:
            print('A survey needs to imported before a memory usage estimate can be made')
            return 0
        else: # number of measurements computation 
            nmeas = []
            for s in self.surveys:
                df = s.df 
                ie = df['irecip'].values >= 0 # count the number of measurements actually put to file 
                nmeas.append(sum(ie))
            num_ind_meas=np.mean(nmeas)
        
        numel = self.mesh.numel
        
        Gb=(numel*num_ind_meas*8)/1.0e9
        dump('ResIPy Estimated RAM usage = %f Gb'%Gb)  
        
        avialMemory = psutil.virtual_memory().available*9.31e-10
        if Gb >= avialMemory:
            dump('*** It is likely that more RAM is required for inversion! ***\n'
                 '*** Make a coarser mesh ***')
        return Gb     
    
    @staticmethod
    def setNcores(ncores):
        """Set the number of cores to use for big calculations, for now 
        this value only to mesh generation and calculations done on 3D meshes.

        Parameters
        ----------
        ncores : int
            Number of cores to use 
        Raises
        ------
        EnvironmentError
            Raised if ncores is more than that avialable

        """
        if ncores>sysinfo['cpuCount']:
            raise EnvironmentError('More cores requested than detected/avialable')
        if type(ncores) != int:
            raise TypeError('Expected ncores as an int, but got %s'%str(type(ncores)))
        mt.ncores = ncores
        
# for backward compatibility, retain the main class called R2
class R2(Project): # :pragma: no cover
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
#%% Resolution matrix calculation
# compute covariance matrix on Nvidia GPU / multi core processor 
def __readSize(fname): # :pragma: no cover
    fh = open(fname)
    line1 = fh.readline().split()
    fh.close()
    return [int(k) for k in line1]

def __readJacob(fname): # :pragma: no cover
    array = []
    fh = open(fname)
    header = fh.readline()
    size = [int(k) for k in header.split()]
    lines = fh.readlines()
    fh.close()
    for line in lines:
        array += [float(x) for x in line.split()]
    return tuple(size), array

def __readRm(fname, jsize, rsize): # :pragma: no cover
    Rmap = np.genfromtxt(fname.replace('R','Rindex'),
                         skip_header=1,dtype=int)-1
    
    Rvals = np.genfromtxt(fname,
                          skip_header=1)
    
    Rn = np.zeros((jsize[1],jsize[1]))
    for i in range(rsize[0]):
        for j in range(rsize[1]):
            k = Rmap[i,j]
            if k!=-1:
                Rn[i,k] = Rvals[i,j]
    
    return Rn 

def __getAlpha(fname): # :pragma: no cover
    #return final reported alpha value, file should be the .out file from andy's code
    fh = open(fname,'r')
    lines = fh.readlines()
    fh.close()
    idx = []
    c = 0
    for line in lines:
        if line.find('Alpha:') != -1:
            idx.append(c)
        c+=1
    if len(idx) == 0:
        raise ValueError("Can't find alpha line")
        
    fnl = lines[max(idx)].split()
    alpha = float(fnl[1])
    return alpha


def cudaRm(invdir): # :pragma: no cover
    """Compute Resolution and Covariance matrix for 2D problems using nVIDIA GPU. 

    Parameters
    ----------
    invdir : string 
        Inversion directory used by R2.

    Returns
    -------
    covar : nd array 
        Values along the diagonal of the coviarance matrix.
    remat : nd array 
        Values along the diagonal of the Resolution matrix..

    """
    import cupy as cp # import cupy (needs cuda enabled pc)
    
    mempool = cp.get_default_memory_pool()
    pinned_mempool = cp.get_default_pinned_memory_pool()
    # read in jacobian
    jsize, jdata = __readJacob(os.path.join(invdir,'f001_J.dat'))
    Jn = np.array(jdata,dtype=np.float32).reshape(jsize) # numpy equivalent 
    J = cp.array(Jn,dtype=np.float32)
    
    # read in data Weighting matrix
    protocol = np.genfromtxt(os.path.join(invdir,'f001_err.dat'),
                             skip_header=1)
    
    Wd = cp.array(np.diag(protocol[:,8]),dtype=np.float32)
    
    # read in model roughness matrix 
    rsize = __readSize(os.path.join(invdir,'f001_R.dat'))
    
    Rn = __readRm(os.path.join(invdir,'f001_R.dat'), jsize, rsize)
    R = cp.array(Rn,dtype=np.float32)
    #construct A and b on GPU 
    files = os.listdir(invdir)
    for f in files:
        if f.endswith('.out'):
            alpha = __getAlpha(os.path.join(invdir,f))
            break
    S = cp.matmul(cp.matmul(J.T,Wd.T), cp.matmul(Wd,J))
    A = S + alpha*R #Form A (Menke et al, 2015)
    
    #get rid of parameters we dont need anymore to free up memory 
    J = None
    Wd = None
    R = None 
    mempool.free_all_blocks()
    
    Cm = cp.linalg.inv(A) # solve inverse of A to get covariance matrix 
    ResM = Cm*S
    A = None
    mempool.free_all_blocks()
    
    # retrieve outputs as numpy arrays 
    covar = np.diagonal(Cm.get())
    remat = np.diagonal(ResM.get())
    
    #finally clear memory once again 
    Cm = None
    S = None
    ResM = None
    mempool.free_all_blocks()
    pinned_mempool.free_all_blocks()
    
    return covar, remat

def parallelRm(invdir): # :pragma: no cover
    """Compute Resolution and Covariance matrix for 2D problems using multicore CPU. 
    Behaves the same as cudaRm but uses numpy / openBlas. 

    Parameters
    ----------
    invdir : string 
        Inversion directory used by R2.

    Returns
    -------
    covar : nd array 
        Values along the diagonal of the coviarance matrix.
    remat : nd array 
        Values along the diagonal of the Resolution matrix..

    """
    # read in jacobian
    jsize, jdata = __readJacob(os.path.join(invdir,'f001_J.dat'))
    Jn = np.array(jdata,dtype=np.float32).reshape(jsize) # numpy equivalent 
    
    # read in data Weighting matrix
    protocol = np.genfromtxt(os.path.join(invdir,'f001_err.dat'),
                             skip_header=1)
    
    Wd = np.diag(protocol[:,8])
    
    # read in model roughness matrix 
    rsize = __readSize(os.path.join(invdir,'f001_R.dat'))
    
    Rn = __readRm(os.path.join(invdir,'f001_R.dat'), jsize, rsize)
    #construct A and b on GPU 
    files = os.listdir(invdir)
    for f in files:
        if f.endswith('.out'):
            alpha = __getAlpha(os.path.join(invdir,f))
            break
    S = np.matmul(np.matmul(Jn.T,Wd.T), np.matmul(Wd,Jn))
    A = S + alpha*Rn #Form A (Menke et al, 2015)
    
    #get rid of parameters we dont need anymore to free up memory 
    Jn = None
    Wd = None
    Rn = None 
    
    Cm = np.linalg.inv(A) # solve inverse of A to get covariance matrix (should use multiple cores)
    ResM = Cm*S
    A = None
    
    # retrieve outputs as numpy arrays 
    covar = np.diagonal(Cm)
    remat = np.diagonal(ResM)
    
    #finally clear memory once again 
    Cm = None
    S = None
    ResM = None
    
    return covar, remat

#%% deprecated funcions

    def pseudoIP(self, index=0, vmin=None, vmax=None, ax=None, **kwargs): # pragma: no cover
        warnings.warn('The function is deprecated, use showPseudoIP() instead.',
                      DeprecationWarning)
        self.showPseudoIP(index=index, vmin=vmin, vmax=vmax, ax=ax, **kwargs)

    def plotError(self, index=0, ax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use showError() instead.',
                      DeprecationWarning)
        self.showError(index=index, ax=ax)

    def errorDist(self, index=0, ax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use showErrorDist() instead.',
                      DeprecationWarning)
        self.showErrorDist(index=index, ax=ax)

    def removeDummy(self, index=-1): # pragma: no cover
        warnings.warn('This function is deprecated, use filterDummy() instead.',
                      DeprecationWarning)
        self.filterDummy(index=index)

    def linfit(self, index=-1, ax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use fitErrorLin() instead.',
                      DeprecationWarning)
        self.fitErrorLin(index=index, ax=ax)


    def pwlfit(self, index=-1, ax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use fitErrorPwl() instead.',
                      DeprecationWarning)
        self.fitErrorPwl(index=index, ax=ax)

    def lmefit(self, index=-1, ax=None, rpath=None, iplot=True): # pragma: no cover
        warnings.warn('This function is deprecated, use fitErrorLME() instead.',
                      DeprecationWarning)
        self.fitErrorLME(index=index, ax=ax, rpath=rpath, iplot=iplot)

    def phaseplotError(self, index=0, ax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use showErrorIP() instead.',
                      DeprecationWarning)
        self.showErrorIP(index=index, ax=ax)

    def plotIPFit(self, index=-1, ax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use fitErrorPwlIP() instead.',
                      DeprecationWarning)
        self.fitErrorPwlIP(index=index, ax=ax)

    def plotIPFitParabola(self, index=-1, ax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use fitErrorParabolaIP() instead.',
                      DeprecationWarning)
        self.fitErrorParabolaIP(index=index, ax=ax)

    def heatmap(self, index=0, ax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use showHeatmap() instead.',
                      DeprecationWarning)
        self.showHeatmap(index=index, ax=ax)

    def removenested(self, index=-1): # pragma: no cover
        warnings.warn('This function is deprecated, use filterNested() instead.',
                      DeprecationWarning)
        self.filterNested(index=index)

    def dca(self, index=-1, dump=None): # pragma: no cover
        warnings.warn('This function is deprecated, use filterDCA() instead.',
                      DeprecationWarning)
        self.filterDCA(index=index, dump=dump)

    def removerecip(self, index=0): # pragma: no cover
        warnings.warn('This function is deprecated, use filterRecip() instead.',
                      DeprecationWarning)
        self.filterRecip(index=index)

    def iprangefilt(self, phimin, phimax, index=-1): # pragma: no cover
        warnings.warn('This function is deprecated, use filterRangeIP() instead.',
                      DeprecationWarning)
        self.filterRangeIP(phimin, phimax, index=index)
        
    def removeUnpaired(self, index=-1): # pragma: no cover
        warnings.warn('This function is deprecated, use filterUnpaired() instead.',
                      DeprecationWarning)
        n = self.filterUnpaired(index=index)
        return n
    
    def removeneg(self): # pragma: no cover
        warnings.warn('This function is deprecated, use filterNegative() instead.',
                      DeprecationWarning)
        self.filterNegative()

    def assignRes0(self, regionValues={}, zoneValues={}, fixedValues={}, ipValues={}): # pragma: no cover
        warnings.warn('This function is deprecated, use setStartingRes() instead.',
                      DeprecationWarning)
        self.setStartingRes(regionValues=regionValues, zoneValues=zoneValues, fixedValues=fixedValues, ipValues=ipValues)


    def assignRefModel(self, res0): # pragma: no cover
        warnings.warn('This function is deprecated, use setRefModel() instead.',
                      DeprecationWarning)
        self.setRefModel(res0=res0)


    def createModellingMesh(self, typ='default', buried=None, surface=None, cl_factor=2,
                   cl=-1, dump=print, res0=100, show_output=True, doi=None, **kwargs): # pragma: no cover
        warnings.warn('This function is deprecated, use createModelErrorMesh() instead.',
                      DeprecationWarning)
        self.createModelErrorMesh(typ=typ, buried=buried, surface=surface, cl_factor=cl_factor,
                                  cl=cl, dump=dump, res0=res0, show_output=show_output, doi=doi, **kwargs)

    def estError(self, a_wgt=0.001, b_wgt=0.02): # pragma: no cover
        warnings.warn('This function is deprecated, use estimateError() instead.',
                      DeprecationWarning)
        self.estimateError(a_wgt=a_wgt, b_wgt=b_wgt)

    def pseudoError(self, ax=None, vmin=None, vmax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use showPseudInvError() instead.',
                      DeprecationWarning)
        self.showPseudoInvError(ax=ax, vmin=vmin, vmax=vmax)


    def pseudoErrorIP(self, ax=None, vmin=None, vmax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use showErrorIP() instead.',
                      DeprecationWarning)
        self.showPseudoErrorIP(ax=ax, vmin=vmin, vmax=vmax)


    def showInversionErrors(self, ax=None): # pragma: no cover
        warnings.warn('This function is deprecated, use showInvError() instead.',
                      DeprecationWarning)
        self.showInvError(ax=ax)


    # def compCond(self): # pragma: no cover
    #     warnings.warn('This function is deprecated, use computeCond() instead.',
    #                   DeprecationWarning)
    #     self.computeCond()


    def pseudo(self, index=0, vmin=None, vmax=None, ax=None, **kwargs):  # pragma: no cover
        warnings.warn('The function is deprecated, use showPseudo() instead.',
                      DeprecationWarning)
        self.showPseudo(index=index, vmin=vmin, vmax=vmax, ax=ax, **kwargs)
