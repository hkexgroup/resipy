# -*- coding: utf-8 -*-
"""
Main R2 class, wraps the other ResIPy modules (API) in to an object orientated approach
@author: Guillaume, Sina, Jimmy and Paul
"""
ResIPy_version = '1.1.5' # ResIPy version (semantic versionning in use) 

#import relevant modules 
import os, sys, shutil, platform, warnings, time # python standard libs
from subprocess import PIPE, call, Popen
import subprocess
import numpy as np # import default 3rd party libaries (can be downloaded from conda repositry, incl with winpython)
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from multiprocessing import Pool, Process, Queue
from threading import Thread

OS = platform.system()
sys.path.append(os.path.relpath('..'))

#import ResIPy resipy packages 
from resipy.Survey import Survey
from resipy.r2in import write2in
import resipy.meshTools as mt
import resipy.isinpolygon as iip
from resipy.template import parallelScript, startAnmt, endAnmt
from resipy.protocol import (dpdp1, dpdp2, wenner_alpha, wenner_beta, wenner,
                          wenner_gamma, schlum1, schlum2, multigrad)
from resipy.SelectPoints import SelectPoints

apiPath = os.path.abspath(os.path.join(os.path.abspath(__file__), '../'))
print('API path = ', apiPath)
print('ResIPy version = ',str(ResIPy_version))

        
def workerInversion(path, dump, exePath, qin, iMoveElec=False):
    ''' Take protocol.dat from the queue (`qin`) and run inversion in its own
    working directory and then put the inverted files back to their origin.
    This function is called by `R2.runParallel()` and shouldn't be used on its
    own.
    
    Parameters
    ----------
    path : str
        Path of the working directory.
    dump : callable
        The output of the running inversion will be passed as text to it. By
        default the `print()` function is used.
    exePath: str
        Abslolute path of the executable to be called to run the inversion.
    qin : queue
        Queue containing filename of protocol.dat file to be inverted for each
        survey.
    iMoveElec : boolean, optional
        If `True`, each protocol.dat file will have it's own .in file that
        can provide different electrode positions.
    '''
    os.chdir(path)
    typ = os.path.basename(exePath).replace('.exe','')
    
    for fname in iter(qin.get, 'stop'):
        # copy the protocol.dat
        shutil.copy(fname, os.path.join(path, 'protocol.dat'))
        name = os.path.basename(fname).replace('.dat', '').replace('protocol_','')
        if iMoveElec is True:
            exeName = os.path.basename(exePath).replace('.exe','')
            r2inFile = os.path.join(os.path.dirname(fname),
                                    exeName + '_' + name + '.in')
            shutil.copy(r2inFile, os.path.join(path, exeName + '.in'))
            
        # run inversion
        if OS == 'Windows':
            cmd = [exePath]
        elif OS == 'Darwin':
            winePath = []
            wine_path = Popen(['which', 'wine'], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
            for stdout_line in iter(wine_path.stdout.readline, ''):
                winePath.append(stdout_line)
            if winePath != []:
                cmd = ['%s' % (winePath[0].strip('\n')), exePath]
            else:
                cmd = ['/usr/local/bin/wine', exePath]
        else:
            cmd = ['wine',exePath]
            
        if OS == 'Windows':
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        
        def execute(cmd):
            if OS == 'Windows':
                proc = Popen(cmd, stdout=PIPE, shell=False, universal_newlines=True, startupinfo=startupinfo)
            else:
                proc = Popen(cmd, stdout=PIPE, shell=False, universal_newlines=True)                
            for stdout_line in iter(proc.stdout.readline, ""):
                yield stdout_line
            proc.stdout.close()
            return_code = proc.wait()
            if return_code:
                print('error on return_code')

        for text in execute(cmd):
                dump(text.rstrip())
    
        # moving inversion results back
        originalDir = os.path.dirname(fname)
        toMove = ['f001_res.dat', 'f001_res.vtk', 'f001_err.dat',
                  'f001_sen.dat', 'f001_diffres.dat',
                  'f001.dat', 'f001.sen', 'f001.err', 'f001.vtk'] # all 3D stuff
        for f in toMove:
            if os.path.exists(os.path.join(path, f)):
                shutil.move(os.path.join(path, f),
                            os.path.join(originalDir, f.replace('f001', name)))
        shutil.move(os.path.join(path, typ + '.out'),
                    os.path.join(originalDir, typ + '_' + name + '.out'))
        shutil.move(os.path.join(path, 'electrodes.dat'),
                    os.path.join(originalDir, 'electrodes_' + name + '.dat'))
        shutil.move(os.path.join(path, 'electrodes.vtk'),
                    os.path.join(originalDir, 'electrodes_' + name + '.vtk'))


# small useful function for reading and writing mesh.dat
def readMeshDat(fname):
    """ Read mesh.dat or mesh3d.dat and returns elements and nodes.
    """
    with open(fname, 'r') as f:
        x = f.readline().split()
    numel = int(x[0])
    elems = np.genfromtxt(fname, skip_header=1, max_rows=numel)
    if fname[-6:] == '3d.dat': # it's a 3D mesh
        skip_footer = 1
    else:
        skip_footer = 0
    nodes = np.genfromtxt(fname, skip_header=numel+1, skip_footer=skip_footer)
    return elems, nodes


def writeMeshDat(fname, elems, nodes, extraHeader='', footer='1'):
    """ Write mesh.dat/mesh3d.dat provided elements and nodes at least.
    """
    numel = len(elems)
    nnodes = len(nodes)
    threed = nodes.shape[1] == 4 # it's a 3D mesh
    if threed is True:
        extraHeader = '\t1\t0\t4'
    with open(fname, 'w') as f:
        f.write('{:.0f} {:.0f}{}\n'.format(numel, nnodes, extraHeader))
    with open(fname, 'ab') as f:
        np.savetxt(f, elems, fmt='%.0f')
        if threed is True:
            np.savetxt(f, nodes, fmt='%.0f %f %f %f')
        else:
            np.savetxt(f, nodes, fmt='%.0f %f %f')
    if threed is True: # for 3D only
        with open(fname, 'a') as f:
            f.write(footer)
                

# distance matrix function for 2D (numpy based from https://stackoverflow.com/questions/22720864/efficiently-calculating-a-euclidean-distance-matrix-using-numpy)
def cdist(a):
    z = np.array([complex(x[0], x[1]) for x in a])
    return np.abs(z[...,np.newaxis]-z)
            
            
class R2(object): # R2 master class instanciated by the GUI
    """ Master class to handle all processing around the inversion codes.
    
    Parameters
    ----------
    dirname : str, optional
        Path of the working directory. Can also be set using `R2.setwd()`.
    typ : str, optional
        Either `R2` or `R3t` for 3D. Complex equivalents are `cR2` and `cR3t`.
        Automatically infered when creating the survey.
    """ 
    def __init__(self, dirname='', typ='R2'):
        self.apiPath = os.path.dirname(os.path.abspath(__file__)) # directory of the code
        if dirname == '':
            dirname = os.path.join(self.apiPath, 'invdir')
        else:
            dirname = os.path.abspath(dirname)
            
        print('Working directory is:', dirname)
        self.setwd(dirname) # working directory (for the datas)
        self.elec = None # will be assigned when creating a survey
        self.surveys = [] # list of survey object
        self.surveysInfo = [] # info about surveys (date)
        self.mesh = None # mesh object (one per R2 instance)
        self.param = {} # dict configuration variables for inversion
        self.configFile = ''
        self.typ = typ # or cR2 or R3, cR3
        self.err = False # if we want error in protocol.dat or not
        self.iBorehole = False # to tell the software to not plot pseudoSection
        self.iTimeLapse = False # to enable timelapse inversion
        self.iBatch = False # to enable batch inversion
        self.meshResults = [] # contains vtk mesh object of inverted section
        self.sequence = None # quadrupoles sequence if forward model
        self.resist0 = None # initial resistivity
        self.iForward = False # if True, it will use the output of the forward
        # to run an inversion (and so need to reset the regions before this)
        self.doi = None # depth of investigation below the surface [in survey units]
        self.proc = None # where the process to run R2/cR2 will be
        self.zlim = None # zlim to plot the mesh by default (from max(elec, topo) to min(doi, elec))
        
    def setwd(self, dirname):
        """ Set the working directory.
        
        Parameters
        ----------
        dirname : str
            Path of the working directory.
        """
        # first check if directory exists
        if os.path.exists(dirname): # ok it exists, let's clear it
            print('clearing the dirname')
            # get rid of some stuff
            files = os.listdir(dirname)
            if 'ref' in files: # only for timelapse survey
                try:
                    shutil.rmtree(os.path.join(dirname, 'ref'))
                except PermissionError:
                    warnings.warn("OS reports reference inversion directory already in use, try changing the working directory")
            if 'err' in files: # only for error modelling
                try:
                    shutil.rmtree(os.path.join(dirname, 'err')) 
                except PermissionError:
                    warnings.warn("OS reports forward modelling directory already in use, try changing the working directory")
            for f in files:
                if (f[:9] == 'protocol_') or (f[:11] == 'electrodes_'):
                    os.remove(os.path.join(dirname , f))
            files2remove = ['R2.in','cR2.in','R3t.in', 'cR3t.in',
                            'R2.out','cR2.out','R3t.out','cR3t.out',
                            'mesh.dat','r100.dat','res0.dat','Start_res.dat',
                            'protocol.dat']
            for f in files2remove:
                if f in files:
                    try:
                        os.remove(os.path.join(dirname, f))
                    except Exception as e:
                        raise Exception("Error setting up working directories:"+ str(e) + "\n...try changing the working inversion directory")
            def isNumber(x):
                try:
                    float(x)
                    return True
                except:
                    return False
            for f in files:
                if (f[0] == 'f') & (isNumber(f[1:3]) is True):
                    os.remove(os.path.join(dirname, f))
        else:
            print('creating the dirname')
            os.mkdir(dirname)
        self.dirname = os.path.abspath(dirname)
    
    
    def setElec(self, elec, elecList=None):
        """ Set electrodes.
        
        Parameters
        ----------
        elec : numpy array
            Array of NxM dimensions. N = number of electodes, M = 2 for x,y or
            M = 3 if x,y,z coordinates are supplied.
        elecList : list, optional
            If not None then elec is ignored in favour of elecList. This option 
            is to be used in the advanced use case where electrodes move which 
            each survey. Each entry of the list is a numpy array the same format
            of 'elec', which is then assigned to each survey class. 
        """
        if elecList is None:
            ok = False
            if self.elec is not None: # then check the shape
                if elec.shape[0] == self.elec.shape[0]:
                    ok = True
                elif self.iForward: # in case of forward modelling, changing the number of electrodes is allowed
                    ok = True
                else:
                    print('ERROR : elec, does not match shape from Survey;')
            else:
                self.elec = elec # first assignement of electrodes
            if ok:
                if elec.shape[1] == 2:
                    self.elec[:,[0,2]] = elec
                    for s in self.surveys:
                        s.elec[:,[0,2]] = elec
                else:
                    self.elec = elec
                    for s in self.surveys:
                        s.elec = elec
        else:
            #some error checking 
            try:
                num_surveys = len(self.surveys) # number of surveys
                if len(elecList) != num_surveys:
                    raise ValueError("The number of electrode matrices must match the number of surveys")
            except AttributeError:
                raise AttributeError("No Survey attribute assocaited with R2 class, make sure you create a survey first")
            
            initElec = elecList[0]
            self.elec = np.zeros((len(initElec),3))
            if elecList[0].shape[1] == 2:
                self.elec[:,[0,2]] = initElec # set R2 class electrodes to initial electrode positions
                for i in range(num_surveys):
                    self.surveys[i].elec[:,[0,2]] = elecList[i] # if 2D set electrode x and elevation coordinates only
            else:
                self.elec = initElec
                for i in range(num_surveys):
                    self.surveys[i].elec = elecList[i] # set survey electrodes to each electrode coordinate
    
        
    def setBorehole(self, val=False):
        """ Set all surveys in borehole type if `True` is passed.
        """
        self.iBorehole = val
        for s in self.surveys:
            s.iBorehole = val
        
        
    def setTitle(self,linetitle):
        """Set the title of the survey name when inverting data. Input is a string.
        """
        if isinstance(linetitle,str):
            self.param['lineTitle']=linetitle
        else:
            print("Cannot set Survey title as input is not a string")
    
    
    def createSurvey(self, fname='', ftype='Syscal', info={}, spacing=None,
                     parser=None):
        """ Read electrodes and quadrupoles data and return 
        a survey object.
        
        Parameters
        ----------
        fname : str
            Filename to be parsed.
        ftype : str, optional
            Type of file to be parsed. Either 'Syscal' or 'Protocol'.
        info : dict, optional
            Dictionnary of info about the survey.
        spacing : float, optional
            Electrode spacing to be passed to the parser function.
        parser : function, optional
            A parser function to be passed to `Survey` constructor.
        """    
        self.surveys.append(Survey(fname, ftype, spacing=spacing, parser=parser))
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
            self.elec = self.surveys[0].elec
            
            # attribute method of Survey object to R2
#            self.pseudoIP = self.surveys[0].pseudoIP
#            self.pseudo = self.surveys[0].pseudo
            self.plotError = self.surveys[0].plotError
            self.errorDist = self.surveys[0].errorDist
            self.linfit = self.surveys[0].linfit
            self.lmefit = self.surveys[0].lmefit
            self.pwlfit = self.surveys[0].pwlfit
            self.phaseplotError = self.surveys[0].phaseplotError
            self.plotIPFit = self.surveys[0].plotIPFit
            self.plotIPFitParabola = self.surveys[0].plotIPFitParabola
            self.heatmap = self.surveys[0].heatmap
            self.iprangefilt = self.surveys[0].iprangefilt
            self.removerecip = self.surveys[0].removerecip
            self.removenested = self.surveys[0].removenested
            self.addFilteredIP = self.surveys[0].addFilteredIP
        
        
    def createBatchSurvey(self, dirname, ftype='Syscal', info={}, spacing=None,
                          parser=None, isurveys=[], dump=print):
        """ Read multiples files from a folders (sorted by alphabetical order).
        
        Parameters
        ----------
        dirname : str
            Directory with files to be parsed.
        ftype : str, optional
            Type of file to be parsed. Either 'Syscal' or 'Protocol'.
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
        """  
        self.createTimeLapseSurvey(dirname=dirname, ftype=ftype, info=info,
                                   spacing=spacing, isurveys=isurveys, 
                                   parser=parser, dump=dump)
        self.iTimeLapse = False
        self.iBatch = True
        self.setBorehole(self.iBorehole)


    def createTimeLapseSurvey(self, dirname, ftype='Syscal', info={},
                              spacing=None, parser=None, isurveys=[],
                              dump=print):
        """ Read electrodes and quadrupoles data and return 
        a survey object.
        
        Parameters
        ----------
        dirname : str
            Directory with files to be parsed.
        ftype : str, optional
            Type of file to be parsed. Either 'Syscal' or 'Protocol'.
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
        """    
        self.iTimeLapse = True
        self.iTimeLapseReciprocal = [] # true if survey has reciprocal
        files = np.sort(os.listdir(dirname))

        for f in files:
            self.createSurvey(os.path.join(dirname, f), ftype=ftype, parser=parser, spacing=spacing)
            haveReciprocal = all(self.surveys[-1].df['irecip'].values == 0)
            self.iTimeLapseReciprocal.append(haveReciprocal)
            dump(f + ' imported')
            print('---------', f, 'imported')
            # all surveys are imported whatever their length, they will be matched
            # later if reg_mode == 2 (difference inversion)
        self.iTimeLapseReciprocal = np.array(self.iTimeLapseReciprocal)
        self.elec = self.surveys[0].elec
        self.setBorehole(self.iBorehole)
        
        # create bigSurvey
        print('creating bigSurvey')
        self.bigSurvey = Survey(os.path.join(dirname, files[0]), ftype=ftype, spacing=spacing)
        # then override the df
        if len(isurveys) == 0: # assume all surveys would be use for error modelling
            isurveys = np.ones(len(self.surveys), dtype=bool)
        isurveys = np.where(isurveys)[0] # convert to indices
        df = self.bigSurvey.df.copy()
        c = 0
        for i in isurveys:
            df2 = self.surveys[i].df
            ipos = df2['irecip'].values > 0
            ineg = df2['irecip'].values < 0
            df2.loc[ipos, 'irecip'] = df2[ipos]['irecip'] + c
            df2.loc[ineg, 'irecip'] = df2[ineg]['irecip'] - c
            df = df.append(df2)
            c = c + df2.shape[0]
        self.bigSurvey.df = df.copy() # override it
        self.bigSurvey.dfOrigin = df.copy()
        self.bigSurvey.ndata = df.shape[0]
#        self.pseudo = self.surveys[0].pseudo # just display first pseudo section

        self.plotError = self.bigSurvey.plotError
        self.errorDist = self.bigSurvey.errorDist
        self.linfit = self.bigSurvey.linfit
        self.lmefit = self.bigSurvey.lmefit
        self.pwlfit = self.bigSurvey.pwlfit
        self.phaseplotError = self.bigSurvey.phaseplotError
        self.plotIPFit = self.bigSurvey.plotIPFit
        self.plotIPFitParabola = self.bigSurvey.plotIPFitParabola
        
        print("%i survey files imported"%len(self.surveys))
        
        
    def pseudo(self, index=0, vmin=None, vmax=None, ax=None, **kwargs):
        """ Plot pseudo-section with dots.
        
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
            Passed to `Survey.pseudo()`.
        """
        self.surveys[index].pseudo(vmin=vmin, vmax=vmax, ax=ax, **kwargs)
        
    
    def pseudoIP(self, index=0, vmin=None, vmax=None, ax=None, **kwargs):
        """ Plot pseudo-section with dots for IP data.
        
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
            Passed to `Survey.pseudo()`.
        """
        self.surveys[index].pseudoIP(vmin=vmin, vmax=vmax, ax=ax, **kwargs)
        
        
    def matchSurveys(self):
        """ Will trim all surveys to get them ready for difference inversion
        where all datasets must have the same number of quadrupoles.
        """
        print('Matching quadrupoles between surveys for difference inversion ...', end='')
        t0 = time.time()
        dfs = [s.df for s in self.surveys]
        
        # sort all dataframe (should already be the case)
        dfs2 = []
        for df in dfs:
            dfs2.append(df)#.sort_values(by=['a','b','m','n']).reset_index(drop=True))
        
        # concatenate columns of string
        def cols2str(cols):
            cols = cols.astype(str)
            x = cols[:,0]
            for i in range(1, cols.shape[1]):
                x = np.core.defchararray.add(x, cols[:,i])
            return x
        
        # get measurements common to all surveys
        df0 = dfs2[0]
        x0 = cols2str(df0[['a','b','m','n']].values.astype(int))
        icommon = np.ones(len(x0), dtype=bool)
        for df in dfs2[1:]:
            x = cols2str(df[['a','b','m','n']].values.astype(int))
            ie = np.in1d(x0, x)
            icommon = icommon & ie
        print(np.sum(icommon), 'in common...', end='')
        
        # create boolean index to match those measurements
        indexes = []
        xcommon = x0[icommon]
        for df in dfs2:
            x = cols2str(df[['a','b','m','n']].values)
            indexes.append(np.in1d(x, xcommon))
        
        print('done in {:.5}s'.format(time.time()-t0)) 
    
        return indexes               
                    
    
    def filterElec(self, elec=[]):
        """ Filter out specific electrodes given in all surveys.
        
        Parameters
        ----------
        elec : list
            List of electrode number to be removed.
        
        """
        for e in elec:
            for i, s in enumerate(self.surveys):
                i2keep = (s.df[['a','b','m','n']].values != e).all(1)
                s.filterData(i2keep)
                print(np.sum(~i2keep), '/', len(i2keep), 'quadrupoles removed in survey', i+1)
    
    
    def filterRecip(self, percent=20):
        """ Filter on reciprocal errors.
        
        Parameters
        ----------
        percent : float, optional
            Percentage of reciprocal error above witch a measurement will be
            discarded. 20% by default.
        """
        numRemoved = 0
        for s in self.surveys:
            numRemoved += s.filterRecip(percent)
        return numRemoved
    
    
    def removeUnpaired(self):
        """ Remove quadrupoles that don't have reciprocals. This might
        remove dummy measurements added for sequence optimization.
        """
        numRemoved = 0
        for s in self.surveys:
            numRemoved += s.removeUnpaired()
        return numRemoved
            
        
    def computeDOI(self):
        """ Compute the Depth Of Investigation (DOI).
        """
        elec = self.elec.copy()
        if all(self.elec[:,1] == 0): # 2D survey:
            if len(self.surveys) > 0:
                array = self.surveys[0].df[['a','b','m','n']].values.copy().astype(int)
                maxDist = np.max(np.abs(elec[array[:,0]-1,0] - elec[array[:,2]-1,0])) # max dipole separation
                self.doi = np.min(elec[:,2])-2/3*maxDist
            else:
                self.doi = np.min(elec[:,2])-2/3*(np.max(elec[:,0]) - np.min(elec[:,0]))

            # set num_xy_poly
            self.param['num_xy_poly'] = 5
            ymax = np.max(self.elec[:,2])
            ymin = self.doi
            xmin, xmax = np.min(self.elec[:,0]), np.max(self.elec[:,0])
            xy_poly_table = np.array([
            [xmin, ymax],
            [xmax, ymax],
            [xmax, ymin],
            [xmin, ymin],
            [xmin, ymax]])
            self.param['xy_poly_table'] = xy_poly_table
        
        else: # for 3D survey
            dist = np.zeros((len(elec), len(elec)))
            for i, el1 in enumerate(elec):
                dist[:,i] = np.sqrt(np.sum((el1[None,:] - elec)**2, axis=1))
            self.doi = np.min(elec[:,2])-2/3*np.max(dist)
            
            # set num_xy_poly
            self.param['num_xy_poly'] = 5
            xmin, xmax = np.min(self.elec[:,0]), np.max(self.elec[:,0])
            ymin, ymax = np.min(self.elec[:,1]), np.max(self.elec[:,1])
            zmin, zmax = self.doi, np.max(self.elec[:,2])
            if (self.typ == 'R2') | (self.typ == 'cR2'): # 2D
                xy_poly_table = np.array([
                [xmin, zmax],
                [xmax, zmax],
                [xmax, zmin],
                [xmin, zmin],
                [xmin, zmax]])
            else:
                xy_poly_table = np.array([
                [xmin, ymax],
                [xmax, ymax],
                [xmax, ymin],
                [xmin, ymin],
                [xmin, ymax]])
                self.param['zmin'] = zmin
                self.param['zmax'] = zmax
            self.param['xy_poly_table'] = xy_poly_table 
        print('computed DOI : {:.2f}'.format(self.doi))
        
    
    def createMesh(self, typ='default', buried=None, surface=None, cl_factor=2,
                   cl=-1, dump=print, res0=100, show_output=True, **kwargs):
        """ Create a mesh.
        
        Parameters
        ----------
        typ : str, optional
            Type of mesh. Eithter 'quad' or 'trian'. If no topography, 'quad'
            mesh will be chosen.
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
        kwargs: -
            Keyword arguments to be passed to mesh generation schemes 
        """
        self.computeDOI()
        
        if typ == 'default':
            if all(self.elec[:,1] == 0): # it's a 2D mesh
                typ = 'quad'
                print('Using a quadrilateral mesh.')
            else:
                typ = 'tetra'
                print('Using a tetrahedral mesh.')
        print(typ)
        
        if typ == 'quad':
            elec = self.elec.copy()
            elec_x = self.elec[:,0]
            elec_z = self.elec[:,2]
            #add buried electrodes? 
            elec_type = np.repeat('electrode',len(elec_x))
            if (buried is not None
                    and elec.shape[0] == len(buried)
                    and np.sum(buried) != 0):
                elec_type[buried]='buried'
                
            elec_type = elec_type.tolist()
            surface_x = surface[:,0] if surface is not None else None
            surface_z = surface[:,1] if surface is not None else None
            mesh,meshx,meshy,topo,e_nodes = mt.quad_mesh(elec_x,elec_z,elec_type,
                                                         surface_x=surface_x, surface_z=surface_z,
                                                         **kwargs)   #generate quad mesh   
            #update parameters accordingly 
            self.param['meshx'] = meshx
            self.param['meshy'] = meshy
            self.param['topo'] = topo
            self.param['mesh_type'] = 4
            self.param['node_elec'] = np.c_[1+np.arange(len(e_nodes[0])), e_nodes[0], e_nodes[1]].astype(int)
            
            if 'regions' in self.param: # allow to create a new mesh then rerun inversion
                del self.param['regions']
            if 'num_regions' in self.param:
                del self.param['num_regions']
        elif typ == 'trian' or typ == 'tetra':
            elec = self.elec.copy()
            geom_input = {}
            elec_x = self.elec[:,0]
            elec_y = self.elec[:,1]
            elec_z = self.elec[:,2]
            elec_type = np.repeat('electrode',len(elec_x))
            if (buried is not None
                    and elec.shape[0] == len(buried)
                    and np.sum(buried) != 0):
                elec_type[buried]='buried'
            
            if surface is not None:
                if surface.shape[1] == 2:
                    geom_input['surface'] = [surface[:,0], surface[:,1]]
                else:
                    geom_input['surface'] = [surface[:,0], surface[:,2]]
            
            whole_space = False
            if buried is not None:
                if np.sum(buried) == len(buried) and surface is None:
                    # all electrodes buried and no surface given
                    whole_space = True

            elec_type = elec_type.tolist()

            #print('elec_type', elec_type)
            ui_dir = os.getcwd()#current working directory (usually the one the ui is running in)
            os.chdir(self.dirname)#change to working directory so that mesh files written in working directory 
#            try:
            if typ == 'trian':
                mesh = mt.tri_mesh(elec_x,elec_z,elec_type,geom_input,
                             path=os.path.join(self.apiPath, 'exe'),
                             cl_factor=cl_factor,
                             cl=cl, dump=dump, show_output=show_output,
                             doi=self.doi-np.max(elec_z), whole_space=whole_space,
                             **kwargs)
            if typ == 'tetra': # TODO add buried
                if cl == -1:
                    dist = cdist(self.elec[:,:2])/2 # half the minimal electrode distance
                    cl = np.min(dist[dist != 0])
                elec_type = None # for now
                mesh = mt.tetra_mesh(elec_x, elec_y, elec_z,elec_type,
                             path=os.path.join(self.apiPath, 'exe'),
                             surface_refinement=surface,
                             cl_factor=cl_factor,
                             cl=cl, dump=dump, show_output=show_output,
                             doi=self.doi-np.max(elec_z), whole_space=whole_space,
                             **kwargs)
#            except Exception as e:
#                print("Mesh generation failed :", e)   
#                pass
            os.chdir(ui_dir)#change back to original directory
            
            self.param['mesh_type'] = 3
#            self.param['num_regions'] = len(mesh.regions)
#            regs = np.array(np.array(mesh.regions))[:,1:]
#            if self.typ == 'R2':
#                regions = np.c_[regs, np.ones(regs.shape[0])*50]
#            if self.typ == 'cR2':
#                regions = np.c_[regs, np.ones(regs.shape[0])*50, np.ones(regs.shape[0])*-0.1]
#            self.param['regions'] = regions
#            self.param['num_xy_poly'] = 5
            # define xy_poly_table (still need to do it here because param['meshx'] is undefined if triangular mesh)
#            doi = np.abs(self.elec[0,0]-self.elec[-1,0])/2
#            xmax = np.max(self.elec[:,0])
#            xmin = np.min(self.elec[:,0])
#            zmax = np.max(self.elec[:,1])
#            zmin = np.min(self.elec[:,1])-self.doi
#            xy_poly_table = np.array([
#            [xmin, zmax],
#            [xmax, zmax],
#            [xmax, zmin],
#            [xmin, zmin],
#            [xmin, zmax]])
#            self.param['xy_poly_table'] = xy_poly_table
#            e_nodes = np.arange(len(self.elec))+1
            e_nodes = mesh.e_nodes + 1 # +1 because of indexing staring at 0 in python
#            if buried is not None:
#                if np.sum(buried) > 0:
#                    enodes = np.zeros(len(e_nodes), dtype=int)
#                    nburied = np.sum(buried)
#                    enodes[~buried] = e_nodes[:-nburied]
#                    enodes[buried] = e_nodes[-nburied:]
#                    e_nodes = enodes
            self.param['node_elec'] = np.c_[1+np.arange(len(e_nodes)), e_nodes].astype(int)
        
        
        self.mesh = mesh
        self.param['mesh'] = mesh
        
        self.param['num_regions'] = 0
        self.param['res0File'] = 'res0.dat'
        
        numel = self.mesh.num_elms
        self.mesh.add_attribute(np.ones(numel)*res0, 'res0') # default starting resisivity [Ohm.m]
        self.mesh.add_attribute(np.ones(numel)*0, 'phase0') # default starting phase [mrad]
        self.mesh.add_attribute(np.ones(numel, dtype=int), 'zones')
        self.mesh.add_attribute(np.zeros(numel, dtype=bool), 'fixed')
        self.mesh.add_attribute(np.zeros(numel, dtype=float), 'iter')
        
        name = 'mesh.dat'
        if self.typ == 'R3t' or self.typ == 'cR3t':
            name = 'mesh3d.dat'
        file_path = os.path.join(self.dirname, name)
        self.mesh.write_dat(file_path)
        
        self.regid = 1 # 1 is the background (no 0 region)
        self.regions = np.ones(len(self.mesh.elm_centre[0]))
        self.resist0 = np.ones(len(self.regions))*100
        
        # define zlim
        if surface is not None:
            zlimMax = np.max([np.max(elec[:,2]), np.max(surface[:,1])])
        else:
            zlimMax = np.max(elec[:,2])
        zlimMin = np.min([np.min(elec[:,2]), self.doi])
        self.zlim = [zlimMin, zlimMax]

        
    def importMesh(self,file_path,mesh_type='tetra',node_pos=None,elec=None,
                   flag_3D=False, res0=100):
        """
        Import mesh from .vtk / .msh / .dat, rather than having <ResIPy> create
        one for you.
        
        Parameters
        ------------
        file_path: str
            File path mapping to the mesh file
        mesh_type: str
            Type of mesh, 'quad', 'trian', 'tetra'
        node_pos: array like, optional
            Array of ints referencing the electrode nodes. If left as none no electrodes 
            will be added to the mesh class. Consider using mesh.move_elec_nodes()
            to add nodes to mesh using their xyz coordinates.
        elec: numpy array, optional
            N*3 numpy array of electrode x,y,z coordinates. Electrode node positions
            will be computed by finding the nearest nodes to the relevant coordinates.
        flag_3D: bool, optional
            Make this true for 3D meshes if importing .msh type. 
        res0 : float, optional 
            Starting resistivity for mesh elements. 
        Returns
        -----------
        mesh: class 
            Added to R2 class
        """
        if (self.typ == 'R3t') or (self.typ == 'cR3t'):
            flag_3D = True
        else:
            flag_3D = False
        self.mesh = mt.custom_mesh_import(file_path, node_pos=node_pos, flag_3D=flag_3D)
        if elec is not None:
            self.mesh.move_elec_nodes(elec[:,0],elec[:,1],elec[:,2])
        
        #add the electrodes to the R2 class
        if elec is not None or node_pos is not None: # then electrode positions should be known
            self.elec = np.array((self.mesh.elec_x, self.mesh.elec_y, self.mesh.elec_z)).T
        else:
            try:
                elec = self.elec
                self.mesh.move_elec_nodes(elec[:,0],elec[:,1],elec[:,2])
            except AttributeError:
                warnings.warn("No electrode nodes associated with mesh! Electrode positions are unknown!")
          
        #R2 class mesh handling 
        e_nodes = np.array(self.mesh.e_nodes) + 1 # +1 because of indexing staring at 0 in python
        self.param['mesh'] = self.mesh
        if mesh_type == 'quad':
            self.param['mesh_type'] = 4
            colx = self.mesh.quadMeshNp() # convert nodes into column indexes 
            self.param['node_elec'] = np.c_[1+np.arange(len(e_nodes)), np.array(colx), np.ones((len(e_nodes,1)))].astype(int)
            #will only work for assuming electrodes are a surface array 
        else:
            self.param['mesh_type'] = 3
            self.param['node_elec'] = np.c_[1+np.arange(len(e_nodes)), e_nodes].astype(int)
        
        
        self.param['num_regions'] = 0
        self.param['res0File'] = 'res0.dat'
        numel = self.mesh.num_elms
        self.mesh.add_attribute(np.ones(numel)*res0, 'res0') # default starting resisivity [Ohm.m]
        self.mesh.add_attribute(np.ones(numel)*0, 'phase0') # default starting phase [mrad]
        self.mesh.add_attribute(np.ones(numel, dtype=int), 'zones')
        self.mesh.add_attribute(np.zeros(numel, dtype=bool), 'fixed')
        self.mesh.add_attribute(np.zeros(numel, dtype=float), 'iter')
        
        name = 'mesh.dat'
        if self.typ == 'R3t' or self.typ == 'cR3t':
            name = 'mesh3d.dat'
        file_path = os.path.join(self.dirname, name)
        self.mesh.write_dat(file_path)
        
        self.regid = 1 # 1 is the background (no 0 region)
        self.regions = np.ones(len(self.mesh.elm_centre[0]))
        self.resist0 = np.ones(len(self.regions))*100
        
        # define zlim
        if self.doi == None:
            self.computeDOI()
        zlimMax = np.max(self.elec[:,2])
        zlimMin = np.min([np.min(self.elec[:,2]), self.doi])
        self.zlim = [zlimMin, zlimMax]
        
        
    def showMesh(self, ax=None):
        """ Display the mesh.
        """
        if self.mesh is None:
            raise Exception('Mesh undefined')
        else:
#            if self.typ[-2] == '3':
#                self.mesh.show(ax=ax, color_bar=False) # show the whole 3D mesh
                # not just the ROI -> maybe we just want to show ROI actually ... TODO
#            else:
            self.mesh.show(ax=ax, color_bar=False, zlim=self.zlim)
    
    
    def write2in(self, param={}):
        """ Create configuration file for inversion.
        
        Parameters
        ----------
        param : dict
            Dictionnary of parameters and values for the inversion settings.
        """
        typ = self.typ
        if (self.err is True) and ('a_wgt' not in self.param):
            self.param['a_wgt'] = 0
            self.param['b_wgt'] = 0
        elif typ[0] != 'c': # DC case
            if 'a_wgt' not in self.param:
                self.param['a_wgt'] = 0.01
            if 'b_wgt' not in self.param:
                self.param['b_wgt'] = 0.02
        elif typ == 'cR2': # TODO what about cR3 ?
            if 'a_wgt' not in self.param:
                self.param['a_wgt'] = 0.02 # variance for magnitude (no more offset)
            if 'b_wgt' not in self.param:
                self.param['b_wgt'] = 2 # mrad
            
        
        if self.param['mesh_type'] == 4:
            self.param['zones'] = self.mesh.attr_cache['zones']
            #TODO reshape it to the form of the mesh
                
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
            if self.err:
                param['a_wgt'] = 0
                param['b_wgt'] = 0
            else:
                if 'a_wgt' not in param:#this allows previously assigned values to be 
                    param['a_wgt'] = 0.01 # written to the reference.in config file
                if 'b_wgt' not in param:
                    param['b_wgt'] = 0.02
            param['num_xy_poly'] = 0
            param['reg_mode'] = 0 # set by default in ui.py too
            param['res0File'] = 'res0.dat'
            if self.typ[-2] == '3':
                print("Writing background inversion config for 3D inversion!")
                param['inverse_type'] = 0 # normal regulurisation
                param['zmin'] = min(self.mesh.node_z)-10 # we want to keep the whole mesh for background regularisation
                param['zmax'] = max(self.mesh.node_z)+10
            self.configFile = write2in(param, refdir, typ=typ) # background survey
            
            # now prepare the actual timelapse settings
            self.param['num_regions'] = 0
            if 'reg_mode' not in self.param.keys():
                self.param['reg_mode'] = 2
            self.param['res0File'] = 'Start_res.dat'
            if self.typ == 'R3t' or self.typ == 'cR3t':
                # for R3t/cR3t, there is no regularization mode but just
                # an inversion_type variable that we need to overwrite based
                # on the reg_mode
                if self.param['reg_mode'] == 2: # notice regularistion mode needs to change 
                    self.param['inverse_type'] = 2 # difference inversion
                else:
                    self.param['inverse_type'] = 1 # constrained to background
            write2in(self.param, self.dirname, typ=typ) # actual time-lapse
        else:
            self.configFile = write2in(self.param, self.dirname, typ=typ)
        
        # write the res0.dat needed for starting resistivity
        if self.iForward is True: # we will invert results from forward
            # inversion so we need to start from a homogeneous model
            print('Setting a homogeneous background model as the survey to \
                  be inverted is from a forward model already.')
            res0 = np.ones(self.mesh.num_elms)*100 # default starting resistivity [Ohm.m]
            self.mesh.add_attribute(res0, 'r100')
            phase2 = np.ones(self.mesh.num_elms)*0
            self.mesh.add_attribute(phase2, 'phase2')
            self.mesh.attr_cache['fixed'] = np.zeros(self.mesh.num_elms, dtype=bool)
            
            if self.typ[0] == 'c' : # we're dealing with IP here !
                r = np.array(self.mesh.attr_cache['r100'])
                phase = np.array(self.mesh.attr_cache['phase2'])
                centroids = np.array(self.mesh.elm_centre).T
                centroids2 = centroids[:,[0,2]] if self.typ[-1] != 't' else centroids
                x = np.c_[centroids2,
                          r,
                          phase, # mrad
                          np.log10(r),
                          np.log10(1/r),
                          np.log10(-10**np.log10(1/r)*phase/1000)]
                np.savetxt(os.path.join(self.dirname, 'res0.dat'), x)
            else:
                self.mesh.write_attr('r100', 'res0.dat', self.dirname)

            
        else: # if we invert field data, we allow the user to define prior
            # knowledge of the resistivity structure            
            if self.typ[0] == 'c' : # we're dealing with IP here !
                r = np.array(self.mesh.attr_cache['res0'])
                phase = np.array(self.mesh.attr_cache['phase0'])
                centroids = np.array(self.mesh.elm_centre).T
                centroids2 = centroids[:,[0,2]] if self.typ[-1] != 't' else centroids
                x = np.c_[centroids2,
                          r,
                          phase, # mrad
                          np.log10(r),
                          np.log10(1/r),
                          np.log10(-10**np.log10(1/r)*phase/1000)]
                np.savetxt(os.path.join(self.dirname, 'res0.dat'), x)
            else:
                self.mesh.write_attr('res0', 'res0.dat', self.dirname)

                
        # rewriting mesh.dat
        paramFixed = 1+ np.arange(self.mesh.num_elms)
        ifixed = np.array(self.mesh.attr_cache['fixed'])
        paramFixed[ifixed] = 0
        
        name = 'mesh.dat'
        if self.typ == 'R3t' or self.typ == 'cR3t':
            name = 'mesh3d.dat'
        self.mesh.write_dat(os.path.join(self.dirname, name),
                            zone = self.mesh.attr_cache['zones'],
                            param = paramFixed)

        # NOTE if fixed elements, they need to be at the end !
        
        # TODO not sure the sorting fixed element issue if for 3D as well
       
        meshFile = os.path.join(self.dirname, name)
        elems, nodes = readMeshDat(meshFile)
        ifixed = elems[:,-2] == 0
        elems2 = np.r_[elems[~ifixed,:], elems[ifixed,:]]
        elems2[:,0] = 1 + np.arange(elems2.shape[0])
        ifixed2 = elems2[:,-2] == 0
        elems2[~ifixed2,-2] = 1 + np.arange(np.sum(~ifixed2))
        writeMeshDat(meshFile, elems2, nodes)
        
        res0File = os.path.join(self.dirname, 'res0.dat')
        resistivityFile = os.path.join(self.dirname, 'resistivity.dat')
        fnames = [res0File, resistivityFile]
        for f in fnames:
            if os.path.exists(f):
                x = np.genfromtxt(f)
                x2 = np.r_[x[~ifixed,:], x[ifixed,:]]
                np.savetxt(f, x2)
                

    def write2protocol(self, err=None, errTot=False, **kwargs):
        """ Write a protocol.dat file for the inversion code.
        
        Parameters
        ----------
        err : bool, optional
            If `True` error columns will be written in protocol.dat provided
            an error model has been fitted or error have been imported.
        errTot : bool, optional
            If `True`, it will compute the modelling error due to the mesh and
            add it to the error from an error model.
        **kwargs : optional
            To be passed to `Survey.write2protocol()`.
        """
        if self.typ[0] == 'c':
            ipBool = True
        else:
            ipBool = False

        if err is None:
            err = self.err
        
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

        # for time-lapse inversion ------------------------------
        if self.iTimeLapse is True:
            if 'reg_mode' not in self.param.keys():
                self.param['reg_mode'] = 2 # by default it's timelapse (difference)
            if self.param['reg_mode'] == 2: # it's a difference inversion
                indexes = self.matchSurveys()
            else:
                indexes = [None]*len(self.surveys)
            # a bit simplistic but assign error to all based on Transfer resistance
            # let's assume it's False all the time for now
            content = ''
            df0 = self.surveys[0].df[['a','b','m','n','resist','recipMean']]
            df0 = df0.rename(columns={'resist':'resist0', 'recipMean':'recipMean0'})
            for i, s in enumerate(self.surveys[1:]):
                if 'resist0' in s.df.columns:
                    s.df = s.df.drop('resist0', axis=1)
                if 'recipMean0' in s.df.columns:
                    s.df = s.df.drop('recipMean0', axis=1)
                s.df = pd.merge(s.df, df0, on=['a','b','m','n'], how='left')
                if err is True:
                    s.df['resError'] = self.bigSurvey.errorModel(s.df)
                res0Bool = False if self.param['reg_mode'] == 1 else True
                protocol = s.write2protocol('', err=err, res0=res0Bool,
                                            isubset=indexes[i+1])
                content = content + str(protocol.shape[0]) + '\n'
                content = content + protocol.to_csv(sep='\t', header=False, index=False)
                
                if i == 0:
                    refdir = os.path.join(self.dirname, 'ref')
                    if os.path.exists(refdir) == False:
                        os.mkdir(refdir)
                    if 'mesh.dat' in os.listdir(self.dirname):
                        shutil.copy(os.path.join(self.dirname, 'mesh.dat'),
                                os.path.join(self.dirname, 'ref', 'mesh.dat'))
                    if 'mesh3d.dat' in os.listdir(self.dirname):
                        shutil.copy(os.path.join(self.dirname, 'mesh3d.dat'),
                                os.path.join(self.dirname, 'ref', 'mesh3d.dat'))
                    s.write2protocol(os.path.join(refdir, 'protocol.dat'),err=err) # no subset for background, just use all
            with open(os.path.join(self.dirname, 'protocol.dat'), 'w') as f:
                f.write(content)
        
        # for batch inversion -------------------
        elif self.iBatch is True:
            content = ''
            for i, s in enumerate(self.surveys):
                if err is True:
                    s.df['resError'] = self.bigSurvey.errorModel(s.df)
                df = s.write2protocol(outputname='', err=err, ip=ipBool, errTot=errTot)
                content = content + str(len(df)) + '\n'
                content = content + df.to_csv(sep='\t', header=False, index=False)
            with open(os.path.join(self.dirname, 'protocol.dat'), 'w') as f:
                f.write(content)
        
        # for normal inversion (one survey) --------------------------
        else:
            self.surveys[0].write2protocol(os.path.join(self.dirname, 'protocol.dat'),
                        err=err, ip=ipBool, errTot=errTot)
        
        
    def runR2(self, dirname='', dump=print):
        """ Run the executable in charge of the inversion.
        
        Parameters
        ----------
        dirname : str, optional
            Path of the directory where to run the inversion code.
        dump : function, optional
            Function to print the output of the invrsion code while running.
        """
        # run R2.exe
        exeName = self.typ + '.exe'
        cwd = os.getcwd()
        if dirname == '':
            dirname = self.dirname
        os.chdir(dirname)
        
        # get R2.exe path
        exePath = os.path.join(self.apiPath, 'exe', exeName)
            
        if OS == 'Windows':
            cmd = [exePath]
        elif OS == 'Darwin':
            winePath = []
            wine_path = Popen(['which', 'wine'], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
            for stdout_line in iter(wine_path.stdout.readline, ''):
                winePath.append(stdout_line)
            if winePath != []:
                cmd = ['%s' % (winePath[0].strip('\n')), exePath]
            else:
                cmd = ['/usr/local/bin/wine', exePath]
        else:
            cmd = ['wine',exePath]
        
        if OS == 'Windows':
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
        
        def execute(cmd):
            if OS == 'Windows':
                self.proc = subprocess.Popen(cmd, stdout=PIPE, shell=False, universal_newlines=True, startupinfo=startupinfo)
            else:
                self.proc = subprocess.Popen(cmd, stdout=PIPE, shell=False, universal_newlines=True)                
            for stdout_line in iter(self.proc.stdout.readline, ""):
                yield stdout_line
            self.proc.stdout.close()
            return_code = self.proc.wait()
            if return_code:
                print('error on return_code')        
        for text in execute(cmd):
                dump(text.rstrip())
        os.chdir(cwd)
    
    def runDistributed(self, dirname=None, dump=print, iMoveElec=False, ncores=None):
        """ run R2 in // according to the number of cores available but in a 
        non-concurrent way (!= runParallel) -> this works on Windows
        
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
        """
        if dirname is None:
            dirname = self.dirname # use default working directory
        
        if self.iTimeLapse is True and self.iBatch is False:
            surveys = self.surveys[1:] # don't take the first survey
        else:
            surveys = self.surveys
            
        # create R2/R3t/cR2/cR3.exe path
        exeName = self.typ + '.exe'
        exePath = os.path.join(self.apiPath, 'exe', exeName)
        
        # split the protocol.dat
        dfall = pd.read_csv(os.path.join(self.dirname, 'protocol.dat'),
                            sep='\t', header=None, engine='python').reset_index()
        idf = list(np.where(np.isnan(dfall[dfall.columns[-1]].values))[0])
        idf.append(len(dfall))
        dfs = [dfall.loc[idf[i]:idf[i+1]-1,:] for i in range(len(idf)-1)]
        
        # trick if moving electrodes for each survey (expand e_nodes)
        node_elec = []
        c = 0
        for s, df in zip(surveys, dfs):
            elec = s.elec
            if int(self.mesh.cell_type[0])==8 or int(self.mesh.cell_type[0])==9:#elements are quads
                colx = self.mesh.quadMeshNp() # so find x column indexes instead. Wont support change in electrode elevation
                node_elec.append(colx)
            else:
                e_nodes = self.mesh.move_elec_nodes(elec[:,0], elec[:,1], elec[:,2])
                node_elec.append(e_nodes+1) #add one to be consistent with fortran indexing
            if (self.typ == 'R3t') or (self.typ == 'cR3t'):
                df.iloc[1:, [2,4,6,8]] = df.iloc[1:, [2,4,6,8]] + c
            else:
                df.iloc[1:, 1:5] = df.iloc[1:, 1:5] + c                
            c += len(elec)
        node_elec = np.hstack(node_elec)
        node_elec = np.c_[np.arange(len(node_elec))+1, node_elec,
                          np.ones(len(node_elec))].astype(int)
        if self.param['node_elec'].shape[1] == 3:
            self.param['node_elec'] = node_elec
        else:
            self.param['node_elec'] = node_elec[:,:-1]
        
        # write a new R2.in
        write2in(self.param, self.dirname, self.typ)
        
        # prepare the worker directory
        if ncores is None:
            ncores = mt.systemCheck()['core_count']
        perWorker = int(len(dfs)/ncores)
        if len(dfs) % ncores > 0:
            perWorker = perWorker + 1 # not sure about that
        print(perWorker, 'surveys distributed per worker.')
        dirs = []
        c = 0
        assignDir = []
        for i in range(ncores):
            cend = c + perWorker
            if cend > len(dfs):
                cend = len(dfs)
            assignDir.append(cend-c)
            # creating the working directory
            wd = os.path.join(dirname, str(i+1))
            dirs.append(wd)
            if os.path.exists(wd):
                shutil.rmtree(wd)
            os.mkdir(wd)
            
            # copying usefull files from the main directory
            toMove = ['mesh.dat', 'mesh3d.dat','R2.in','cR2.in',
                      'R3t.in', 'cR3t.in', 'res0.dat','resistivity.dat', 
                      'Start_res.dat']
            for f in toMove:
                fname = os.path.join(dirname, f)
                if os.path.exists(fname):
                    shutil.copy(fname, os.path.join(wd, f))
            
            # reformed protocol.dat
            df = pd.concat(dfs[c:cend])
            print(c, '->', cend, 'in', wd)
            outputname = os.path.join(wd, 'protocol.dat')
            df.to_csv(outputname, sep='\t', header=False, index=False)
            # header with line count already included
            
            c = cend
            if cend == len(dfs):
                break # stop the loop here
        
        # run them all in parallel as child processes
        def dumpOutput(out):
            for line in iter(out.readline, ''):
                dump(line.rstrip())
            out.close()
            
        if OS == 'Windows':
            cmd = [exePath]
        elif OS == 'Darwin':
            winePath = []
            wine_path = Popen(['which', 'wine'], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
            for stdout_line in iter(wine_path.stdout.readline, ''):
                winePath.append(stdout_line)
            if winePath != []:
                cmd = ['%s' % (winePath[0].strip('\n')), exePath]
            else:
                cmd = ['/usr/local/bin/wine', exePath]
        else:
            cmd = ['wine', exePath]
            
        if OS == 'Windows':
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            
        procs = []
        ts = []
        for wd in dirs:
            if OS == 'Windows':
                p = Popen(cmd, stdout=PIPE, cwd=wd, shell=False, universal_newlines=True, startupinfo=startupinfo)
            else:
                p = Popen(cmd, stdout=PIPE, cwd=wd, shell=False, universal_newlines=True)                
            procs.append(p)
            t = Thread(target=dumpOutput, args=(p.stdout,))
            t.daemon = True # thread dies with the program
            t.start()
            ts.append(t)

        class ProcsManagement(object): # little class to handle the kill
            def __init__(self, procs):
                self.procs = procs
            def kill(self):
                for p in self.procs:
                    p.terminate()
                    
        self.proc = ProcsManagement(procs)
        
        for p, t in zip(procs, ts): # block until all processes finished
            p.wait()
                
        # get the files as if it was a sequential inversion
        if self.typ=='R3t' or self.typ=='cR3t':
            toRename = ['.dat', '.vtk', '.err', '.sen', '_diffres.dat']
        else:
            toRename = ['_res.dat', '_res.vtk', '_err.dat', '_sen.dat', '_diffres.dat']
        r2outText = ''
        c = 0
        for i, d in enumerate(dirs):
            for j in range(assignDir[i]):
                c += 1
                fname = 'f' + str(j+1).zfill(3)
                for ext in toRename:
                    originalFile = os.path.join(d,  fname + ext)
                    newFile = os.path.join(dirname, 'f' + str(c).zfill(3) + ext)
                    if os.path.exists(originalFile):
                        shutil.move(originalFile, newFile)
                r2outFile = os.path.join(d, self.typ + '.out')
                with open(r2outFile, 'r') as f:
                    r2outText = r2outText + f.read()
        with open(os.path.join(dirname, self.typ + '.out'), 'w') as f:
            f.write(r2outText)
        
        # delete the dirs and the files
        [shutil.rmtree(d) for d in dirs]
        
        print('----------- END OF DISTRIBUTED INVERSION ----------')
        
        
    
    def runParallel(self, dirname=None, dump=print, iMoveElec=False, 
                    ncores=None, rmDirTree=True):
        """ Run R2 in // according to the number of cores available.
        
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
        
        if self.iTimeLapse is True and self.iBatch is False:
            surveys = self.surveys[1:]
        else:
            surveys = self.surveys
            
        # create R2.exe path
        exeName = self.typ + '.exe'
        exePath = os.path.join(self.apiPath, 'exe', exeName)
        
        # split the protocol.dat
        dfall = pd.read_csv(os.path.join(self.dirname, 'protocol.dat'),
                            sep='\t', header=None, engine='python').reset_index()
        
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
            print('Electrodes position will be updated for each survey')
            for s in self.surveys:
                print(s.name, '...', end='')
                elec = s.elec
                e_nodes = self.mesh.move_elec_nodes(elec[:,0], elec[:,1], elec[:,2])
                self.param['node_elec'][:,1] = e_nodes + 1 # WE MUST ADD ONE due indexing differences between python and fortran
                if int(self.mesh.cell_type[0])==8 or int(self.mesh.cell_type[0])==9:#elements are quads
                    colx = self.mesh.quadMeshNp() # so find x column indexes instead. Wont support change in electrode elevation
                    self.param['node_elec'][:,1] = colx
                self.param['inverse_type'] = 1 # regularise against a background model 
                #self.param['reg_mode'] = 1
                write2in(self.param, self.dirname, self.typ)
                r2file = os.path.join(self.dirname, self.typ + '.in')
                shutil.move(r2file, r2file.replace('.in', '_' + s.name + '.in'))
                print('done')
        queueIn = Queue() # queue
        
        # create workers directory
        ncoresAvailable = ncores = mt.systemCheck()['core_count']
        if ncores is None:
            ncores = ncoresAvailable
        else:
            if ncores > ncoresAvailable:
                raise ValueError('Number of cores larger than available')
        
        procs = []
        dirs = []
        for i in range(ncores):
            # creating the working directory
            wd = os.path.join(dirname, str(i+1))
            dirs.append(wd)
            if os.path.exists(wd):
                shutil.rmtree(wd)
            os.mkdir(wd)
            
            # copying usefull files from the main directory
            toMove = ['mesh.dat', 'mesh3d.dat','R2.in','cR2.in',
                      'R3t.in', 'cR3t.in', 'res0.dat','resistivity.dat', 
                      'Start_res.dat']
            for f in toMove:
                fname = os.path.join(dirname, f)
                if os.path.exists(fname):
                    shutil.copy(fname, os.path.join(wd, f))
            
            # creating the process
            procs.append(Process(target=workerInversion,
                                 args=(wd, dump, exePath, queueIn, iMoveElec)))
            procs[-1].start()
            
        # feed the queue
        for f in files:
            queueIn.put(f) # this will trigger the inversion
        
        # when finished stop the processes
        for i in range(ncores):
            queueIn.put('stop')
        
        for p in procs: # this blocks until all processes have finished
            p.join()
        
        class ProcsManagement(object): # little class to handle the kill
            def __init__(self, procs):
                self.procs = procs
            def kill(self):
                for p in self.procs:
                    p.terminate()
                    
        self.proc = ProcsManagement(procs)
        
        # get the files as it was a sequential inversion
        if self.typ=='R3t' or self.typ=='cR3t':
            toRename = ['.dat', '.vtk', '.err', '.sen', '_diffres.dat']
        else:
            toRename = ['_res.dat', '_res.vtk', '_err.dat', '_sen.dat', '_diffres.dat']
        r2outText = ''
        for i, s in enumerate(surveys):
            for ext in toRename:
                originalFile = os.path.join(dirname,  s.name + ext)
                newFile = os.path.join(dirname, 'f' + str(i+1).zfill(3) + ext)
                if os.path.exists(originalFile):
                    shutil.move(originalFile, newFile)
            r2outFile = os.path.join(dirname, self.typ + '_' + s.name + '.out')
            print(r2outFile)
            with open(r2outFile, 'r') as f:
                r2outText = r2outText + f.read()
            os.remove(r2outFile)
        with open(os.path.join(dirname, self.typ + '.out'), 'w') as f:
            f.write(r2outText)
        
        # delete the dirs and the files
        if rmDirTree:
            [shutil.rmtree(d) for d in dirs]
            [os.remove(f) for f in files]
        
        print('----------- END OF INVERSION IN // ----------')
    
    
    def runParallel2(self, dirname=None, dump=print, iMoveElec=False, 
                    ncores=None, rmDirTree=True):
        """ Run R2 in // according to the number of cores available.
        
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
        
        if self.iTimeLapse is True and self.iBatch is False:
            surveys = self.surveys[1:]
        else:
            surveys = self.surveys
            
        # create R2.exe path
        exeName = self.typ + '.exe'
        exePath = os.path.join(self.apiPath, 'exe', exeName)

        
        # split the protocol.dat
        dfall = pd.read_csv(os.path.join(self.dirname, 'protocol.dat'),
                            sep='\t', header=None, engine='python').reset_index()
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
            print('Electrodes position will be updated for each survey')
            for s in self.surveys:
                print(s.name, '...', end='')
                elec = s.elec
                e_nodes = self.mesh.move_elec_nodes(elec[:,0], elec[:,1], elec[:,2])
                self.param['node_elec'][:,1] = e_nodes + 1 # WE MUST ADD ONE due indexing differences between python and fortran
                if int(self.mesh.cell_type[0])==8 or int(self.mesh.cell_type[0])==9:#elements are quads
                    colx = self.mesh.quadMeshNp() # so find x column indexes instead. Wont support change in electrode elevation
                    self.param['node_elec'][:,1] = colx
                self.param['inverse_type'] = 1 # regularise against a background model 
                #self.param['reg_mode'] = 1
                write2in(self.param, self.dirname, self.typ)
                r2file = os.path.join(self.dirname, self.typ + '.in')
                shutil.move(r2file, r2file.replace('.in', '_' + s.name + '.in'))
                print('done')

        # create workers directory
        ncoresAvailable = ncores = mt.systemCheck()['core_count']
        if ncores is None:
            ncores = ncoresAvailable
        else:
            if ncores > ncoresAvailable:
                raise ValueError('Number of cores larger than available')
        
        
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
            cmd = [exePath]
        elif OS == 'Darwin':
            winePath = []
            wine_path = Popen(['which', 'wine'], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
            for stdout_line in iter(wine_path.stdout.readline, ''):
                winePath.append(stdout_line)
            if winePath != []:
                cmd = ['%s' % (winePath[0].strip('\n')), exePath]
            else:
                cmd = ['/usr/local/bin/wine', exePath]
        else:
            cmd = ['wine',exePath]

        if OS == 'Windows':
            startupinfo = subprocess.STARTUPINFO()
            startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
#        
#        def execute(cmd):
#            if OS == 'Windows':
#                proc = Popen(cmd, stdout=PIPE, shell=False, universal_newlines=True, startupinfo=startupinfo)
#            else:
#                proc = Popen(cmd, stdout=PIPE, shell=False, universal_newlines=True)                
#            for stdout_line in iter(proc.stdout.readline, ""):
#                yield stdout_line
#            proc.stdout.close()
#            return_code = proc.wait()
#            if return_code:
#                print('error on return_code')
#
#        for text in execute(cmd):
#                dump(text.rstrip())
#    

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
                dump(line.rstrip())
            out.close()
        
        # create essential attribute
        self.irunParallel2 = True
        self.procs = []
        
        # kill management
        class ProcsManagement(object): # little class to handle the kill
            def __init__(self, r2object):
                self.r2 = r2object
            def kill(self):
                print('killing ...')
                self.r2.irunParallel2 = False # this will end the infinite loop
                print('kk')
                procs = self.r2.procs # and kill the running processes
                for p in procs:
                    p.terminate()
                print('all done')
                    
        self.proc = ProcsManagement(self)
        
        # run in // (http://code.activestate.com/recipes/577376-simple-way-to-execute-multiple-process-in-parallel/)
        # In an infinite loop, will run an number of process (according to the number of cores)
        # the loop will check when they finish and start new ones.
        def done(p):
            return p.poll() is not None
#        def success(p):
#            return p.returncode == 0 # this doesn't work so well in compiled windows binaries
#        def fail():
#            sys.exit(1)
                
#        ts = []
        c = 0
        dump('{:.0f}/{:.0f} inversions completed'.format(c, len(wds2)))
        while self.irunParallel2:
            while wds and len(self.procs) < ncores:
                wd = wds.pop()
#                print('task', wd)
                if OS == 'Windows':
                    p = Popen(cmd, cwd=wd, stdout=PIPE, shell=False, universal_newlines=True, startupinfo=startupinfo)
                else:
                    p = Popen(cmd, cwd=wd, stdout=PIPE, shell=False, universal_newlines=True) 
                self.procs.append(p)
#                t = Thread(target=dumpOutput, args=(p.stdout,))
#                t.daemon = True # thread dies with the program
#                t.start()
#                ts.append(t)
    
            for p in self.procs:
                if done(p):
                    self.procs.remove(p)
                    c = c+1
                    dump('{:.0f}/{:.0f} inversions completed'.format(c, len(wds2)))
    
            if not self.procs and not wds:
                dump('')
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
        if self.typ=='R3t' or self.typ=='cR3t':
            toRename = ['.dat', '.vtk', '.err', '.sen', '_diffres.dat']
        else:
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
        
        # delete the dirs and the files
        if rmDirTree:
            [shutil.rmtree(d) for d in wds2]
            [os.remove(f) for f in files]
        
        print('----------- END OF INVERSION IN // ----------')
    
    
    def runParallelWindows(self, dirname=None, dump=print, iMoveElec=False, 
                    ncores=None, rmDirTree=False):
        """ Run R2 in // according to the number of cores available.
        
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
            Remove excess directories and files created during parallel inversion
        """
        if dirname is None:
            dirname = self.dirname
        
        if self.iTimeLapse is True and self.iBatch is False:
            surveys = self.surveys[1:] # skips out first survey as this should be inverted seperately as a baseline
        else:
            surveys = self.surveys
            
        # create R2.exe path
        exeName = self.typ + '.exe'
        exePath = os.path.join(self.apiPath, 'exe', exeName)
        
        # split the protocol.dat
        dfall = pd.read_csv(os.path.join(self.dirname, 'protocol.dat'),
                            sep='\t', header=None, engine='python').reset_index()
        
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
            print('Electrodes position will be updated for each survey')
            for s in self.surveys:
                print(s.name, '...', end='')
                elec = s.elec
                e_nodes = self.mesh.move_elec_nodes(elec[:,0], elec[:,1], elec[:,2])
                self.param['node_elec'][:,1] = e_nodes + 1 # WE MUST ADD ONE due indexing differences between python and fortran
                if int(self.mesh.cell_type[0])==8 or int(self.mesh.cell_type[0])==9:#elements are quads
                    colx = self.mesh.quadMeshNp() # so find x column indexes instead. Wont support change in electrode elevation
                    self.param['node_elec'][:,1] = colx
                self.param['inverse_type'] = 1 # regularise against a background model 
                #self.param['reg_mode'] = 1
                write2in(self.param, self.dirname, self.typ)
                r2file = os.path.join(self.dirname, self.typ + '.in')
                shutil.move(r2file, r2file.replace('.in', '_' + s.name + '.in')) # each file has different node positions
                print('done')   
                
        print('Creating working directories for each inversion...',end='')
        workerDirs=['']*len(surveys) # number of surveys, preallocate list 
        toMove = ['mesh.dat', 'mesh3d.dat','R2.in','cR2.in',
                  'R3t.in', 'cR3t.in', 'res0.dat','resistivity.dat', 
                  'Start_res.dat']
        
        for i,s in enumerate(surveys):
            workerDirs[i] = os.path.join(self.dirname,'%i'%(i+1)) # add one so starting directory ==1 
            if not os.path.exists(workerDirs[i]): # make inversion directory 
                os.mkdir(workerDirs[i])
            for f in toMove: # copy mesh file, starting resistivity etc... 
                fname = os.path.join(self.dirname, f)
                if os.path.exists(fname):
                    shutil.copy(fname, os.path.join(workerDirs[i], f))
            #copy over protocol file 
            fname = os.path.join(self.dirname, 'protocol_' + s.name + '.dat')
            shutil.copy(fname,os.path.join(workerDirs[i], 'protocol.dat'))
            #copy .in file if moving electrodes occur
            if iMoveElec:
                fname = os.path.join(self.dirname, self.typ + '_' + s.name + '.in')
                shutil.copy(fname,os.path.join(workerDirs[i], self.typ +'.in'))
        print('Done')
                        
        # create workers directory
        ncoresAvailable = ncores = mt.systemCheck()['core_count']
        if ncores is None:
            ncores = ncoresAvailable
        else:
            if ncores > ncoresAvailable:
                raise ValueError('Number of cores larger than available')
        
        #we need to add some formating strings to worker directories in the python script
        drive = self.dirname[0]#get the working drive for the inversion directory 
        formattedDirs = str(workerDirs).replace('\\\\','\\').replace(',',',\n').replace("'"+drive+":","r'"+drive+":")
        
        #write a parallised python scrip using an imported template 
        parallel = parallelScript.format(formattedDirs,
                                         "r'"+exePath+"'",ncores)#format template
        
        #according to docs found here: 
        #https://docs.python.org/2/library/multiprocessing.html#multiprocessing-programming
        #we need to spawn a whole new python interpreter in order to safely run 
        #multiprocessing within windows as the os doesnt support forking. 
        #A workaround I found is here:
        #https://stackoverflow.com/questions/45110287/workaround-for-using-name-main-in-python-multiprocessing
        
        #write python parallel script to file 
        fh = open('parallelScript.py','w')
        fh.write(parallel)
        fh.close()
        
        #now to run the actual inversion                 
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW   
        python_interpreter = sys.executable
        filename = 'parallel.log'
        with open(filename, 'w') as writer, open(filename, 'r', 1) as reader:
            p = subprocess.Popen([python_interpreter, 'parallelScript.py'], 
                                 stdout=writer, 
                                 universal_newlines=True, 
                                 startupinfo=startupinfo)
            while p.poll() is None:
                sys.stdout.write(reader.read())# Update what is printed to screen
                time.sleep(0.05)
            sys.stdout.write(reader.read())#read whats left after the process has finished 
            
            return_code = p.wait()
            if return_code:
                print('error on return_code')    
         
        class ProcsManagement(object): # little class to handle the kill
            def __init__(self, procs):
                self.procs = procs
            def kill(self):
                for p in self.procs:
                    p.terminate()          
        self.proc = ProcsManagement([p])
          
        # delete the dirs and the files
        if rmDirTree:
            [shutil.rmtree(d) for d in workerDirs]
            #[os.remove(f) for f in files]
        
        print('------------- END OF PARALLISED INVERSION -------------')
        
        
    def invert(self, param={}, iplot=False, dump=print, modErr=False,
               parallel=False, iMoveElec=False, ncores=None, forceParallel=False):
        """ Invert the data, first generate R2.in file, then run
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
            default all the cores available are used).)
        """
        # clean meshResults list
        self.meshResults = []
        
        # create mesh if not already done
        if 'mesh' not in self.param:
            dump('Create Rectangular mesh...')
            self.createMesh()
            dump('done\n')
            
        # compute modelling error if selected
        if modErr is True:
            dump('Computing error model ...')
            self.computeModelError()
            dump('done\n')
            errTot = True
        else:
            errTot = False
        
        # write configuration file
        dump('Writing .in file and protocol.dat ...')
        self.write2in(param=param) # R2.in
        self.write2protocol(errTot=errTot) # protocol.dat
        #check to make sure the number of electrodes in the protocal matches the
        #the number of electrodes.
        df = self.surveys[0].df
        check = np.array((df['a'],df['b'],df['m'],df['n']))
        if len(self.elec) < np.max(check): # Make sure there are not more electrodes locations in the schedule file than in R2 class
            raise Exception("The number of electrodes given to ResIPy (%i) does not match the number of electrodes parsed in the scheduling file (%i)."%(len(self.elec),np.max(check)))
        dump('done\n')
        
        # runs inversion
        if self.iTimeLapse == True:
            dump('------------ INVERTING REFERENCE SURVEY ---------------\n')
            refdir = os.path.join(self.dirname, 'ref')
            shutil.move(os.path.join(self.dirname,'res0.dat'),
                        os.path.join(refdir, 'res0.dat'))
            self.write2in(param=param)
            self.runR2(refdir, dump=dump) # this line actaully runs R2
            if self.typ=='R3t' or self.typ=='cR3t':
                shutil.copy(os.path.join(refdir, 'f001.dat'),
                            os.path.join(self.dirname, 'Start_res.dat'))
            else:
                shutil.copy(os.path.join(refdir, 'f001_res.dat'),
                            os.path.join(self.dirname, 'Start_res.dat'))
  
        dump('--------------------- MAIN INVERSION ------------------\n')
        if parallel is True and (self.iTimeLapse is True or self.iBatch is True):
#            if platform.system() == "Windows": # different method needed on windows due to lack of forking
#                if forceParallel:
#                    self.runParallelWindows(dump=dump, iMoveElec=iMoveElec,ncores=ncores)
#                else:
#                    self.runDistributed(dump=dump,iMoveElec=iMoveElec,ncores=ncores)
#            else:
            self.runParallel2(dump=dump, iMoveElec=iMoveElec, ncores=ncores)
        else:
            self.runR2(dump=dump)
        
        if iplot is True:
            self.showResults()
            
            
    def showResults(self, index=0, ax=None, edge_color='none', attr='',
                    sens=True, color_map='viridis', zlim=None, clabel=None,
                    **kwargs):
        """ Show the inverteds section.
        
        Parameters
        ----------
        index : int, optional
            Index of the inverted section (mainly in the case of time-lapse
            inversion)
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
        """
        if len(self.meshResults) == 0:
            self.getResults()
        if (attr == '') & (self.typ[0] != 'c'):
            attr = 'Resistivity(log10)'
        if (attr == '') & (self.typ[0] == 'c'):
            attr = 'Sigma_real(log10)'
        keys = list(self.meshResults[index].attr_cache.keys())
        if attr not in keys:
            print('Attribute not found, revert to default')
            attr = keys[0]
        if len(self.meshResults) > 0:
            if zlim is None:
                zlim = self.zlim
            if self.typ[-1] == '2': # 2D case
                self.meshResults[index].show(ax=ax, edge_color=edge_color,
                                attr=attr, sens=sens, color_map=color_map,
                                zlim=zlim, clabel=clabel, **kwargs)
            else: # 3D case
                self.meshResults[index].show(ax=ax,
                            attr=attr, color_map=color_map, clabel=clabel,
                            **kwargs)
        else:
            print('Unexpected Error')

    
    def getResults(self):
        """ Collect inverted results after running the inversion and adding
        them to `R2.meshResults` list.
        """
        self.meshResults = [] # make sure we empty the list first
        if self.iTimeLapse == True:
            if self.typ[-2] == '3':
                fresults = os.path.join(self.dirname, 'ref', 'f001.vtk')
            else:
                fresults = os.path.join(self.dirname, 'ref', 'f001_res.vtk')
            print('reading ref', fresults)
            mesh = mt.vtk_import(fresults)
            mesh.mesh_title = self.surveys[0].name
            mesh.elec_x = self.elec[:,0]
            mesh.elec_y = self.elec[:,1]
            mesh.elec_z = self.elec[:,2]
            mesh.surface = self.mesh.surface
            self.meshResults.append(mesh)
        if self.iForward is True:
            initMesh = mt.vtk_import(os.path.join(self.dirname, 'fwd','forward_model.vtk'))
            initMesh.elec_x = self.elec[:,0]
            initMesh.elec_y = self.elec[:,1]
            initMesh.elec_z = self.elec[:,2]
            initMesh.surface = self.mesh.surface
            self.meshResults.append(initMesh)
            
        for i in range(len(self.surveys)):
            if self.iTimeLapse is True:
                j = i + 1
            else:
                j = i
            if self.typ[-2] == '3':
                fresults = os.path.join(self.dirname, 'f' + str(i+1).zfill(3) + '.vtk')
            else:
                fresults = os.path.join(self.dirname, 'f' + str(i+1).zfill(3) + '_res.vtk')
            if os.path.exists(fresults):
                print('reading ', fresults, '...', end='')
                try:
                    mesh = mt.vtk_import(fresults)
                    mesh.mesh_title = self.surveys[j].name
                    mesh.elec_x = self.surveys[j].elec[:,0]
                    mesh.elec_y = self.surveys[j].elec[:,1]
                    mesh.elec_z = self.surveys[j].elec[:,2]
                    mesh.surface = self.mesh.surface
                    self.meshResults.append(mesh)
                    print('done')
                except Exception as e:
                    print('failed', e)
            else:
                break
        
        # compute conductivity in mS/m
        for mesh in self.meshResults:
            if 'Resistivity(Ohm-m)' in mesh.attr_cache.keys():
                mesh.attr_cache['Conductivity(mS/m)'] = 1000/np.array(mesh.attr_cache['Resistivity(Ohm-m)'])
        
        # compute difference in percent in case of reg_mode == 1
        if (self.iTimeLapse is True) and (self.param['reg_mode'] == 1):
            try: 
                self.compDiff()
            except:
                pass
#            resRef = np.array(self.meshResults[0].attr_cache['Resistivity(Ohm-m)']) 
#NOTE: this will not work as the reference array will be bigger than the timesteps if the mesh is cropped
#            for mesh in self.meshResults[1:]:
#                res = np.array(mesh.attr_cache['Resistivity(Ohm-m)'])
#                mesh.attr_cache['difference(percent)'] = (res-resRef)/resRef*100
            
            
    def showSection(self, fname='', ax=None, ilog10=True, isen=False, figsize=(8,3)):
        """ Show inverted section based on the `_res.dat``file instead of the
        `.vtk`.
        
        Parameters
        ----------
        fname : str, optional
            Name of the inverted `.dat` file produced by the inversion.
        ax : matplotlib axis, optional
            If specified, the graph will be plotted along `ax`.
        ilog10 : bool, optional
            If `True`, the log10 of the resistivity will be used.
        isen : bool, optional
            If `True`, sensitivity will be displayed as white transparent
            shade on top of the inverted section.
        figsize : tuple, optional
            Size of the figure.
        """
        print('showSection called (to be discarded in the futur)')
        if fname == '':
            fname = os.path.join(self.dirname, 'f001.dat')
        res = pd.read_csv(fname, delimiter=' *', header=None, engine='python').values
        lenx = len(np.unique(res[:,0]))
        leny = len(np.unique(res[:,1]))
        x = res[:,0].reshape((leny, lenx), order='F')
        y = res[:,1].reshape((leny, lenx), order='F')
        z = res[:,2].reshape((leny, lenx), order='F')
        if isen:
            sen = pd.read_csv(fname.replace('res','sen'), delimiter=' *', header=None, engine='python').values
            lenx = len(np.unique(sen[:,0]))
            leny = len(np.unique(sen[:,1]))
            zs = sen[:,2].reshape((leny, lenx), order='F')
            zs = np.log10(zs)
            zs -= np.min(zs)
            alpha = zs/np.max(zs)
            print(np.max(alpha), np.min(alpha))
        if ilog10:
            z = np.log10(z)
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()
        cax = ax.pcolormesh(x, y, z)
        ax.plot(self.elec[:,0], self.elec[:,2], 'ko')
        cbar = fig.colorbar(cax, ax=ax)
        if ilog10:
            cbar.set_label(r'$\log_{10}(\rho) [\Omega.m]$')
        else:
            cbar.set_label(r'$\rho [\Omega.m]$')
        ax.set_ylabel('Depth [m]')
        ax.set_xlabel('Distance [m]')

    
    def addRegion(self, xy, res0=100, phase0=1, blocky=False, fixed=False, ax=None):
        """ Add region according to a polyline defined by `xy` and assign it
        the starting resistivity `res0`.
        
        Parameters
        ----------
        xy : array
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
        """
        if ax is None:
            fig, ax = plt.subplots()
        self.mesh.show(ax=ax)
        selector = SelectPoints(ax, np.array(self.mesh.elm_centre).T[:,[0,2]],
                                typ='poly') # LIMITED FOR 2D case
        selector.setVertices(xy)
        selector.getPointsInside()
        idx = selector.iselect
        self.regid = self.regid + 1
        self.regions[idx] = self.regid
        self.mesh.cell_attributes = list(self.regions) # overwriting regions
        self.resist0[idx] = res0
        self.mesh.attr_cache['res0'] = self.resist0 # hard way to do it
        phase = self.mesh.attr_cache['phase0'].copy()
        phase[idx] = phase0
        self.mesh.attr_cache['phase0'] = phase
        
        # define zone
        if blocky is True:
            zones = self.mesh.attr_cache['zones'].copy()
            zones[idx] = self.regid
            self.mesh.attr_cache['zones'] = zones
        
        # define fixed area
        if fixed is True:
            print('fixing')
            paramFixed = self.mesh.attr_cache['fixed'].copy()
            paramFixed[idx] = True
            self.mesh.attr_cache['fixed'] = paramFixed
            print('sum = ', np.sum(paramFixed == True))
        
        
    def resetRegions(self):
        """ Just reset all regions already draw. Shouldn't be needed as 
        the `self.runR2()` automatically use a homogenous model as starting
        for inversion. The only purpose of this is to use an inhomogeous
        starting model to invert data from forward modelling.
        """
        self.regid = 1
        self.regions.fill(1)
        self.mesh.attr_cache['res0'] = np.ones(len(self.regions))*100 # set back as default
        
        
    def createModel(self, ax=None, dump=print, typ='poly', addAction=None):
        """ Interactive model creation for forward modelling.
        
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
        if self.mesh is None:
            print('will create a mesh before')
            self.createMesh()
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        def callback(idx):
            print('nb elements selected:', np.sum(idx))
            self.regid = self.regid + 1
            self.regions[idx] = self.regid            
            self.mesh.cell_attributes = list(self.regions) # overwritin regions            
            self.mesh.draw()
            if addAction is not None:
                addAction()
        self.mesh.atribute_title = 'Material'
        self.mesh.show(ax=ax, zlim=self.zlim)
        # we need to assign a selector to self otherwise it's not used
        self.selector = SelectPoints(ax, np.array(self.mesh.elm_centre).T[:,[0,2]],
                                     typ=typ, callback=callback)
        if ax is None:
            return fig
            
    
    def assignRes0(self, regionValues={}, zoneValues={}, fixedValues={}, ipValues={}):
        """ Assign starting resitivity values.
        
        Parameters
        ----------
        regionValues : dict
            Dictionnary with key being the region number and the value being
            the resistivity in [Ohm.m].
        zoneValues : dict
            Dictionnary with key being the region number and the zone number.
            There would be no smoothing between the zones if 'block inversion'
            is selected (`inversion_type` = 4).
        fixedValues : dict
            Dictionnary with key being the region number and a boolean value if
            we want to fix the resistivity of the zone to the starting one.
            Note that it only works for triangular mesh for now.
        ipValues : dict
            Dictionnary with key being the region number and the values beeing
            the phase [mrad].
        
        Note
        ----
        Region 0 is the background region. It has zone=1, and fixed=False
        """
        res0 = np.array(self.mesh.attr_cache['res0']).copy()
        for key in regionValues.keys():
            idx = self.regions == key
            res0[idx] = regionValues[key]
        self.mesh.attr_cache['res0'] = res0
        
        zones = np.array(self.mesh.attr_cache['zones']).copy()
        for key in zoneValues.keys():
            idx = self.regions == key
            zones[idx] = zoneValues[key]
        self.mesh.attr_cache['zones'] = zones
        
        fixed = np.array(self.mesh.attr_cache['fixed']).copy()
        for key in fixedValues.keys():
            idx = self.regions == key
            fixed[idx] = fixedValues[key]
        self.mesh.attr_cache['fixed'] = fixed
        
        phase0 = np.array(self.mesh.attr_cache['phase0']).copy()
        for key in ipValues.keys():
            idx = self.regions == key
            phase0[idx] = ipValues[key]
        self.mesh.attr_cache['phase0'] = phase0
        
        print('assignRes0-------------', np.sum(fixed), np.sum(phase0))
        

    def createSequence(self, params=[('dpdp1', 1, 8)]):
        """ Create a dipole-dipole sequence.
        
        Parameters
        ----------
        params : list of tuple, optional
            Each tuple is the form (<array_name>, param1, param2, ...)
            Types of sequences available are : 'dpdp1','dpdp2','wenner_alpha',
            'wenner_beta', 'wenner_gamma', 'schlum1', 'schlum2', 'multigrad'.
        
        Examples
        --------
        >>> k = R2()
        >>> k.setElec(np.c_[np.linspace(0,5.75, 24), np.zeros((24, 2))])
        >>> k.createMesh(typ='trian')
        >>> k.createSequence([('dpdp1', 1, 8), ('wenner_alpha', 1), ('wenner_alpha', 2)])
        """
        qs = []
        nelec = len(self.elec)
        fdico = {'dpdp1': dpdp1,
              'dpdp2': dpdp2,
              'wenner': wenner,
              'wenner_alpha': wenner_alpha,
              'wenner_beta': wenner_beta,
              'wenner_gamma': wenner_gamma,
              'schlum1': schlum1,
              'schlum2': schlum2,
              'multigrad': multigrad}
        
        for p in params:
            pok = [int(p[i]) for i in np.arange(1, len(p))] # make sure all are int
            qs.append(fdico[p[0]](nelec, *pok).values.astype(int))
        self.sequence = np.vstack(qs)
    
    
    def importElec(self, fname=''):
        """ Import electrodes positions.
        
        Parameters
        ----------
        fname : str
            Path of the CSV file containing the electrodes positions. It should contains 3 columns maximum with the X, Y, Z positions of the electrodes.
        """
        elec = pd.read_csv(fname, header=None).values
        if elec.shape[1] > 3:
            raise ValueError('The file should have no more than 3 columns')
        else:
            self.setElec(elec)            
    
    
    def importSequence(self, fname=''):
        """ Import sequence for forward modelling.
        
        Parameters
        ----------
        fname : str
            Path of the CSV file to be imported. The file shouldn't have any headers just 4 columns with the 4 electrodes numbers.
        """
        seq = pd.read_csv(fname, header=None)
        if seq.shape[1] != 4:
            raise ValueError('The file should be a CSV file wihtout headers with exactly 4 columns with electrode numbers.')
        else:
            self.sequence = seq
    
    
    def forward(self, noise=0.0, noiseIP=0.0, iplot=False, dump=print):
        """ Operates forward modelling.
        
        Parameters
        ----------
        noise : float, optional 0 <= noise <= 1
            Noise level from a Gaussian distribution that should be applied
            on the forward apparent resistivities obtained. 
        noiseIP : float, optional
            Absolute noise level in mrad from a Gaussian distribution that should be applied
            on the forward phase values obtained. 
        iplot : bool, optional
            If `True` will plot the pseudo section after the forward modelling.
        dump : function, optional
            Function to print information messages when running the forward model.
        """
        fwdDir = os.path.join(self.dirname, 'fwd')
        if os.path.exists(fwdDir):
            shutil.rmtree(fwdDir)
        os.mkdir(fwdDir)
             
        if self.typ[0] == 'c':
            r = np.array(self.mesh.attr_cache['res0'])
            phase = np.array(self.mesh.attr_cache['phase0'])
            centroids = np.array(self.mesh.elm_centre).T
            centroids2 = centroids[:,[0,2]] if self.typ[-1] != 't' else centroids
            x = np.c_[centroids2,
                      r,
                      phase, # mrad
                      np.log10(r),
                      np.log10(1/r),
                      np.log10(-10**np.log10(1/r)*phase/1000)]
            np.savetxt(os.path.join(fwdDir, 'resistivity.dat'), x)
        else:
            self.mesh.write_attr('res0', 'resistivity.dat', fwdDir)
        
        if os.path.exists(os.path.join(self.dirname, 'mesh.dat')) is True:
            shutil.copy(os.path.join(self.dirname, 'mesh.dat'),
                        os.path.join(fwdDir, 'mesh.dat'))
        if os.path.exists(os.path.join(self.dirname, 'mesh3d.dat')) is True:
            shutil.copy(os.path.join(self.dirname, 'mesh3d.dat'),
                        os.path.join(fwdDir, 'mesh3d.dat'))
        
        # write the forward .in file
        dump('Writing .in file...')
        fparam = self.param.copy()
        fparam['job_type'] = 0
        fparam['num_regions'] = 0
        fparam['res0File'] = 'resistivity.dat' # just starting resistivity
        write2in(fparam, fwdDir, typ=self.typ)
        dump('done\n')
        
        # write the protocol.dat (that contains the sequence)
        if self.sequence is None:
            dump('Creating sequence ...')
            self.createSequence()
            dump('done\n')
        dump('Writing protocol.dat ...')
        seq = self.sequence
        
        # let's check if IP that we have a positive geometric factor
        if self.typ[0] == 'c': # NOTE this does'nt work for borehole
            elecpos = self.elec[:,0].copy() # and works only for 2D
            array = seq.copy()
            apos = elecpos[array[:,0]-1]
            bpos = elecpos[array[:,1]-1]
            mpos = elecpos[array[:,2]-1]
            npos = elecpos[array[:,3]-1]
            AM = np.abs(apos-mpos)
            BM = np.abs(bpos-mpos)
            AN = np.abs(apos-npos)
            BN = np.abs(bpos-npos)
            K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
            ie = K < 0
            seq2 = seq.copy()
            seq[ie,2] = seq2[ie,3] # swap if K is < 0
            seq[ie,3] = seq2[ie,2]
            
        protocol = pd.DataFrame(np.c_[1+np.arange(seq.shape[0]),seq])            
        outputname = os.path.join(fwdDir, 'protocol.dat')
        with open(outputname, 'w') as f:
            f.write(str(len(protocol)) + '\n')
        with open(outputname, 'a') as f:
            protocol.to_csv(f, sep='\t', header=False, index=False)
        dump('done\n')
        
        # fun the inversion
        dump('Running forward model')
        self.runR2(fwdDir, dump=dump) # this will copy the R2.exe inside as well
        self.iForward = True
        
        # create a protocol.dat file (overwrite the method)
        def addnoise(x, level=0.05):
            return x + np.random.randn(1)*x*level
        
        def addnoiseIP(x, level=2):
            return x + np.random.randn(1)*level
        
        addnoise = np.vectorize(addnoise)
        addnoiseIP = np.vectorize(addnoiseIP)
        self.noise = noise #proportional noise, e.g. 0.05 = 5% noise
        self.noiseIP = noiseIP #absolute noise in mrad, following convention of cR2
        
        elec = self.elec.copy()
        self.surveys = [] # need to flush it (so no timeLapse forward)
        if self.typ[0] == 'c':
            self.createSurvey(os.path.join(fwdDir, self.typ + '_forward.dat'), ftype='ProtocolIP')
        else:
            self.createSurvey(os.path.join(fwdDir, self.typ + '_forward.dat'), ftype='forwardProtocolDC')
        # NOTE the 'ip' columns here is in PHASE not in chargeability
        self.surveys[0].kFactor = 1 # kFactor by default is = 1 now, though wouldn't hurt to have this here!
        self.surveys[0].df['resist'] = addnoise(self.surveys[0].df['resist'].values, self.noise)
        self.surveys[0].df['ip'] = addnoiseIP(self.surveys[0].df['ip'].values, self.noiseIP)
        self.setElec(elec) # using R2.createSurvey() overwrite self.elec so we need to set it back
        
        self.pseudo()
        dump('Forward modelling done.')

        
        
    def computeModelError(self):
        """ Compute modelling error due to the mesh.
        """
        if self.mesh is None:
            raise ValueError('You fist need to generate a mesh to compute the modelling error.')
            return
        fwdDir = os.path.join(self.dirname, 'err')
        if os.path.exists(fwdDir):
            shutil.rmtree(fwdDir)
        os.mkdir(fwdDir)
        
        # write the resistivity.dat and fparam
        fparam = self.param.copy()
        fparam['job_type'] = 0
        centroids = np.array(self.mesh.elm_centre).T
        if self.param['mesh_type'] == 4:
            fparam['num_regions'] = 1
            maxElem = centroids.shape[0]
            fparam['regions'] = np.array([[1, maxElem, 100]])
        else:
            if '2' in self.typ:
                n = 2
            else:
                n = 3
            resFile = np.zeros((centroids.shape[0],n+1)) # centroix x, y, z, res0
            resFile[:,-1] = 100
            np.savetxt(os.path.join(fwdDir, 'resistivity.dat'), resFile,
                       fmt='%.3f')
            if os.path.exists(os.path.join(self.dirname, 'mesh.dat')) is True:
                shutil.copy(os.path.join(self.dirname, 'mesh.dat'),
                            os.path.join(fwdDir, 'mesh.dat'))
            if os.path.exists(os.path.join(self.dirname, 'mesh3d.dat')) is True:
                shutil.copy(os.path.join(self.dirname, 'mesh3d.dat'),
                            os.path.join(fwdDir, 'mesh3d.dat'))
            fparam['num_regions'] = 0
            fparam['res0File'] = 'resistivity.dat'
        write2in(fparam, fwdDir, typ=self.typ)
        
        # write the protocol.dat based on measured sequence
        seq = self.surveys[0].df[['a','b','m','n']].values
        protocol = pd.DataFrame(np.c_[1+np.arange(seq.shape[0]),seq])
        if all(self.elec[:,1] == 0) is False: # it's a 3D survey
            protocol.insert(1, 'sa', 1)
            protocol.insert(3, 'sb', 1)
            protocol.insert(5, 'sm', 1)
            protocol.insert(7, 'sn', 1)
        outputname = os.path.join(fwdDir, 'protocol.dat')
        with open(outputname, 'w') as f:
            f.write(str(len(protocol)) + '\n')
        with open(outputname, 'a') as f:
            protocol.to_csv(f, sep='\t', header=False, index=False)
    
        # fun the inversion
        self.runR2(fwdDir) # this will copy the R2.exe inside as well
        
        # get error model
        x = np.genfromtxt(os.path.join(fwdDir, self.typ + '_forward.dat'), skip_header=1)
        modErr = np.abs(100-x[:,-1])/100
        for s in self.surveys:
            s.df['modErr'] = modErr
        
        # eventually delete the directory to spare space
        shutil.rmtree(fwdDir)
        
    
    def showIter(self, index=-2, ax=None):
        """ Dispay temporary inverted section after each iteration.
        
        Parameters
        ----------
        index : int, optional
            Iteration number to show.
        ax : matplotib axis, optional
            If specified, the graph will be plotted along `ax`.
        """
        if ax is None:
            fig, ax = plt.subplots()
            iplot = True
        else:
            fig = ax.figure
            iplot = False
        files = os.listdir(self.dirname)
        fs = []
        for f in files:
            if (f[-8:] == '_res.dat') & ((len(f) == 16) or (len(f) == 12)):
                fs.append(f)
        fs = sorted(fs)
        if len(fs) > 1: # the last file is always open and not filled with data
            if self.param['mesh_type'] == 10:
                self.showSection(os.path.join(self.dirname, fs[index]), ax=ax)
                # TODO change that to full meshTools
                
            else:
                x = np.genfromtxt(os.path.join(self.dirname, fs[index]))
                if x.shape[0] > 0:
                    triang = tri.Triangulation(x[:,0],x[:,1])
                    cax = ax.tricontourf(triang, x[:,3], extend='both')
                    # TODO might want to crop surface here as well
                    fig.colorbar(cax, ax=ax, label=r'$\rho$ [$\Omega$.m]')
                    ax.plot(self.elec[:,0], self.elec[:,2], 'ko')
                    ax.set_aspect('equal')
                    ax.set_xlabel('Distance [m]')
                    ax.set_ylabel('Elevation [m]')
                    if iplot is True:
                        fig.show()

    
    def saveInvPlots(self, outputdir=None, **kwargs):
        """ Save all plots to output (or working directory). Parameters
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
        else:
            for i in range(len(self.meshResults)):
                kwargs2 = kwargs.copy()
                fig, ax = plt.subplots()
                if 'ylim' not in kwargs2:
                    ylim = [self.doi, np.max(self.elec[:,2])]
                    kwargs2 = dict(kwargs2, ylim=ylim)
                if 'color_map' not in kwargs2:
                    kwargs2 = dict(kwargs2, color_map='viridis')
                if 'attr' in kwargs2:
                    if kwargs2['attr'] not in list(self.meshResults[i].attr_cache.keys()):
                        kwargs2['attr'] = 'Resistivity(log10)' 
                self.meshResults[i].show(ax=ax, **kwargs2)
                fname = self.surveys[i].name
                fig.savefig(os.path.join(outputdir, fname + '.png'))
    
    
    def getInvError(self):
        """ Collect inversion error from _err.dat or .err file after inversion.
        
        Returns
        -------
        array : numpy.array
            Contains the quadrupoles.
        errors : numpy.array
            Vector of normalized error.
        """
#        if self.typ == 'R2': # old format
#            err = np.genfromtxt(os.path.join(self.dirname, 'f001_err.dat'), skip_header=1)        
#            array = err[:,[-2,-1,-4,-3]].astype(int)
#            errors = err[:,0]
        if self.typ == 'cR2' or self.typ == 'R2':
            df = pd.read_csv(os.path.join(self.dirname, 'f001_err.dat'), delim_whitespace=True)
            array = np.array([df['C+'],df['C-'],df['P+'],df['P-']],dtype=int).T
            errors = np.array(df['Normalised_Error'])
        elif self.typ == 'R3t' or self.typ == 'cR3t':
            err = np.genfromtxt(os.path.join(self.dirname, 'f001.err'), skip_header=1)        
            array = err[:,[-3,-1,-7,-5]].astype(int)
            errors = err[:,0]
            
        return array, errors
    
            
    def pseudoError(self, ax=None):
        """ Plot pseudo section of errors from file `f001_err.dat`.
        
        Parameters
        ----------
        ax : matplotlib axis
            If specified, the graph will be plotted against `ax`.
        """
        array, errors = self.getInvError()
            
        spacing = np.diff(self.elec[[0,1],0])
        pseudo(array, errors, spacing, ax=ax, label='Normalized Errors', log=False, geom=False, contour=False)
    
    
    def pseudoErrorIP(self, ax=None):
        """ Display normalized phase error.
        """
        if self.typ == 'cR2':
            df = pd.read_csv(os.path.join(self.dirname, 'f001_err.dat'), delim_whitespace=True)   
            array = np.array([df['C+'],df['C-'],df['P+'],df['P-']],dtype=int)
            errors = np.array(df['Calculated_Phase']-df['Observed_Phase'])
        spacing = np.diff(self.elec[[0,1],0])
        pseudo(array.T, errors, spacing, ax=ax, label='Normalized Errors', log=False, geom=False, contour=False)
    
        
    def showInversionErrors(self, ax=None):
        """ Display inversion error by measurment numbers.
        """
#        if self.typ == 'R2':
#            file_path = os.path.join(self.dirname, 'f001_err.dat')
#            err = np.genfromtxt(file_path,skip_header=1)
#            errors = err[:,0]
#        if self.typ == 'cR2':
#            file_path = os.path.join(self.dirname, 'f001_err.dat')
#            err = np.genfromtxt(file_path,skip_header=1)
#            errors = err[:,4]
#        if self.typ == 'R3t':
#            file_path = os.path.join(self.dirname, 'f001.err')
#            err = np.genfromtxt(file_path,skip_header=1)
#            errors = err[:,0]
#        if self.typ == 'cR3t':
#            file_path = os.path.join(self.dirname, 'f001.err')
#            err = np.genfromtxt(file_path,skip_header=1)
#            errors = err[:,4]
        _, errors = self.getInvError()
        measurement_no = np.arange(1,len(errors)+1)
        #make figure
        if ax is None: 
            fig, ax = plt.subplots() 
        ax.scatter(measurement_no,errors)
        ax.set_ylabel("Normalised Error")
        ax.set_xlabel("Measurement Number")
        #add diagnositic lines
        y_pos_limit = (3,3)
        y_neg_limit = (-3,-3)
        baseline = (0,0)
        ax.plot((1,measurement_no[-1]),y_pos_limit,'r--')
        ax.plot((1,measurement_no[-1]),y_neg_limit,'r--')
        ax.plot((1,measurement_no[-1]),baseline,'k--')


    def showInParaview(self, index=0, paraview_loc=None):
        """ Open paraview to display the .vtk file.
        
        Parameters
        ----------
        index: int, optional
            Timestep to be shown in paraview (for an individual survey this 1).
        paraview_loc: str, optional
            **Windows ONLY** maps to the excuatable paraview.exe. The program 
            will attempt to find the location of the paraview install if not given. 
        """
        if self.typ[-1] == '2':
            fname = 'f{:03d}_res.vtk'.format(index+1)
        else:
            fname = 'f{:03d}.vtk'.format(index+1)  
        if OS == "Windows":
            if paraview_loc is None:
                found,cmd_line = self.mesh.findParaview()
                if not found:
                    print("Cannot find paraview location")
                    return
                cmd_line = '"' + cmd_line + '" ' + os.path.join(self.dirname, fname)
            elif isinstance(paraview_loc,str):
                cmd_line = '"' + paraview_loc + '" ' + os.path.join(self.dirname, fname)
            else:
                print("Cannot find where paraview is installed")
                return
        else:
            cmd_line = 'paraview ' + os.path.join(self.dirname, fname)
            
        try:#try and launch paraview
            #Popen([cmd_line, os.path.join(self.dirname, fname)])
            os.popen(cmd_line)
        except PermissionError:
            print("Your operating system has blocked launching Paraview")
            #windows doesnt like calling paraview from python for some reason
            #will need to look into this further. 


    def showSlice(self, index=0, ax=None, attr=None, axis='z'): 
        """ Show slice of 3D mesh interactively.
        """
        if attr is None:
            attr = list(self.meshResults[index].attr_cache.keys())[0]
        self.meshResults[index].showSlice(
                attr=attr, axis=axis)
        
    ## Sorting electrode numbers ## 
    def shuntIndexes(self):
        """
        Shunt electrode indexes to start at 1. 
        """
        debug=True
        if len(self.surveys)>1:
            debug=False
        for i in range(len(self.surveys)):
            self.surveys[i].shuntIndexes(debug=debug)
        
    def normElecIdx(self):
        """
        Normalise electrode indexes to start at 1 in consective and ascending order. 
        """
        debug = True
        if len(self.surveys)>1:
            debug=False
        for i in range(len(self.surveys)):
            self.surveys[i].normElecIdx(debug=debug)
    
    ## make 3d coordinates for a 2d line in to 2d ##     
    def elecXY2elecX(self,yDominant=False,iMoveElec=False):
        """
        Convert 3D electrode XY coordinates into just X coordinates. Use for 
        2D lines only! 
        If self.elec has been set then each survey will use the electrodes set 
        in the R2 master class. If not then the R2 master class will take on the
        elec values set for the first survey in a sequence. 
        
        Parameters
        -----------
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
                elec = self.surveys[i].elec.copy()
                x = elec[:,0]
                y = elec[:,1]
                self.surveys[i].elec[:,0]=y
                self.surveys[i].elec[:,1]=x
        
        for i in range(len(self.surveys)):
            self.surveys[i].elec2distance() # go through each survey and compute electrode
        self.elec = self.surveys[0].elec
        
    def compCond(self):
        """Compute conductivities from resistivities for the ERT mesh
        """
        if self.typ=='R3t' or self.typ=='cR3t':
            res_name = 'Resistivity'
        else:
            res_name = 'Resistivity(Ohm-m)'
        for i in range(len(self.meshResults)):
            self.meshResults[i].reciprocal(res_name,'Conductivity(S/m)')
            
            
    def compDiff(self):
        """Compute difference in meshResult parameters 
        """
        if not self.iTimeLapse:
            raise Exception("Difference calculation only available for time lapse surveys")
        if len(self.meshResults) == 0:
            self.getResults()
        
        crop=False
        if len(self.param['xy_poly_table'])>0:
            meshx = np.array(self.meshResults[0].elm_centre[0])
            meshy = np.array(self.meshResults[0].elm_centre[1])
            meshz = np.array(self.meshResults[0].elm_centre[2])
            crop=True
            if self.typ[-2]=='3':
                inside1 = iip.isinpolygon(meshx,meshy,(self.param['xy_poly_table'][:,0],self.param['xy_poly_table'][:,1]))
                inside2 = (meshz > self.param['zmin']) & (meshz < self.param['zmax'])
                inside = (inside1==True) & (inside2==True)
            else:
                inside = iip.isinpolygon(meshx,meshz,(self.param['xy_poly_table'][:,0],self.param['xy_poly_table'][:,1]))
                        
        num_attr = len(self.meshResults[0].attr_cache)
        num_elm = self.meshResults[0].num_elms
        baselines = np.zeros((num_attr,num_elm))
        for i, key in enumerate(self.meshResults[0].attr_cache):
            baselines[i,:] = self.meshResults[0].attr_cache[key]
        change = np.zeros_like(baselines)
        new_keys = []
        baseline_keys = []
        for j, key in enumerate(self.meshResults[0].attr_cache):
            new_keys.append('Difference('+key+')')
            baseline_keys.append(key)
        for j, key in enumerate(new_keys):
            self.meshResults[0].add_attribute(change[j,:],key)
        
        #filter baseline to just the measurements left over after cropping the mesh
        if crop:
            baselines = baselines[:,inside]
    
        problem = 0
        for i in range(1,len(self.meshResults)): 
            step = self.meshResults[i]
            new_keys = []
            count = 0
            change = np.zeros_like(baselines)
            for j, key in enumerate(baseline_keys):
                try:
                    change[count,:] = (np.array(step.attr_cache[key])-baselines[count,:])/baselines[count,:] * 100
                except KeyError:
                    problem+=1
                new_keys.append('Difference('+key+')')
                count += 1
            count = 0
            for j, key in enumerate(new_keys):
                self.meshResults[i].add_attribute(change[count,:],key)
                count += 1
        if problem>0:
            print('Had a problem computing differences for %i attributes'%problem)
                

    def saveVtks(self,dirname,prefix='ResIPyoutput'):
        """Save vtk files of inversion results to a specified directory. Format
        for file names will be 'prefix'xxx.vtk. A python script will also be saved
        to the relevant directory 

        Parameters
        ------------
        dirname: str
            Direcotry in which results will be saved
        prefix: str, optional
            Characters appended to the front of each file name, ie by default
            files will be named "ResIPyoutput"+"xxx.vtk", where x is the survey
            number. For timelapse surveys "...001.vtk" will be the baseline 
            survey.
        """   
        amtContent = startAnmt 
        if len(self.meshResults) == 0:
            self.getResults()
        count=0
        for mesh, s in zip(self.meshResults, self.surveys):
            count+=1
            file_path = os.path.join(dirname, prefix + '{:03d}.vtk'.format(count))
            mesh.write_vtk(file_path,title=mesh.mesh_title)
            amtContent += "\tannotations.append('%s')\n"%mesh.mesh_title
        amtContent += endAnmt 
        fh = open(os.path.join(dirname,'amt_track.py'),'w')
        fh.write(amtContent)
        fh.close()
        
    def showParam(self):
        """ Print parameters in `R2.param` dictionary.
        """
        [print(key) for i,key in enumerate(self.param)]
    
    
    def filterZeroMeasSurveys(self):
        """ Filter out badly behaved surveys, where after all other QC no measurements 
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


def pseudo(array, resist, spacing, label='', ax=None, contour=False, log=True, geom=True):
    array = np.sort(array, axis=1) # for better presentation, especially Wenner arrays
    nelec = np.max(array)
    elecpos = np.arange(0, spacing*nelec, spacing)
    resist = resist
    
    if geom: # compute and applied geometric factor
        apos = elecpos[array[:,0]-1]
        bpos = elecpos[array[:,1]-1]
        mpos = elecpos[array[:,2]-1]
        npos = elecpos[array[:,3]-1]
        AM = np.abs(apos-mpos)
        BM = np.abs(bpos-mpos)
        AN = np.abs(apos-npos)
        BN = np.abs(bpos-npos)
        K = 2*np.pi/((1/AM)-(1/BM)-(1/AN)+(1/BN)) # geometric factor
        resist = resist*K

#    array = np.sort(array, axis=1)
        
    if log:
        resist = np.sign(resist)*np.log10(np.abs(resist))
    if label == '':
        if log:
            label = r'$\log_{10}(\rho_a)$ [$\Omega.m$]'
        else:
            label = r'$\rho_a$ [$\Omega.m$]'
    
    cmiddle = np.min([elecpos[array[:,0]-1], elecpos[array[:,1]-1]], axis=0) \
        + np.abs(elecpos[array[:,0]-1]-elecpos[array[:,1]-1])/2
    pmiddle = np.min([elecpos[array[:,2]-1], elecpos[array[:,3]-1]], axis=0) \
        + np.abs(elecpos[array[:,2]-1]-elecpos[array[:,3]-1])/2
    xpos = np.min([cmiddle, pmiddle], axis=0) + np.abs(cmiddle-pmiddle)/2
    ypos = - np.sqrt(2)/2*np.abs(cmiddle-pmiddle)
    
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    cax = ax.scatter(xpos, ypos, c=resist, s=70)#, norm=mpl.colors.LogNorm())
    cbar = fig.colorbar(cax, ax=ax)
    cbar.set_label(label)
    ax.set_title('Pseudo Section')
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Pseudo depth [m]')

