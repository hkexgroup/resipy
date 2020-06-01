# -*- coding: utf-8 -*-
"""
Main R2 class, wraps the other ResIPy modules (API) in to an object orientated approach
@author: Guillaume, Sina, Jimmy and Paul
"""
ResIPy_version = '2.2.0' # ResIPy version (semantic versionning in use)

#import relevant modules
import os, sys, shutil, platform, warnings, time # python standard libs
from subprocess import PIPE, call, Popen

# used to download the binaries
import requests
import hashlib

import subprocess
import numpy as np # import default 3rd party libaries (can be downloaded from conda repositry, incl with winpython)
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
#import resipy.interpolation as interp # for cropSurface()
import matplotlib.patches as mpatches
import matplotlib.path as mpath

OS = platform.system()
sys.path.append(os.path.relpath('..'))

#import ResIPy resipy packages
from resipy.Survey import Survey
from resipy.r2in import write2in
import resipy.meshTools as mt
from resipy.meshTools import cropSurface
import resipy.geomTools as iip
from resipy.template import parallelScript, startAnmt, endAnmt
from resipy.protocol import (dpdp1, dpdp2, wenner_alpha, wenner_beta, wenner,
                          wenner_gamma, schlum1, schlum2, multigrad)
from resipy.SelectPoints import SelectPoints
from resipy.saveData import (write2Res2DInv, write2csv, writeSrv)

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
              'b483d7001b57e148453b33412d614dd95e71db21',
              '3fda8d15377974b10cafa5b32398d6e2d2e40345',
              # '577c02cf87bcd2d64cccff14919d607e79ff761a',
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
            response = requests.get("https://gitlab.com/hkex/pyr2/-/raw/master/src/resipy/exe/" + exe)
            with open(fname, 'wb') as f:
                f.write(response.content)
            print('done')
        else:
            print('{:s} found and up to date.'.format(exe))
                
checkExe(os.path.join(apiPath, 'exe'))
            
#%% wine check
def wineCheck():
    #check operating system
    OpSys=platform.system()
    #detect wine
    if OpSys == 'Linux':
        p = Popen("wine --version", stdout=PIPE, shell=True)
        is_wine = str(p.stdout.readline())
        if is_wine.find("wine") == -1:
            print('wine could not be found on your system. resipy needs wine to run the inversion. You can install wine by running `sudo apt-get install wine-stable`.')
        else:
            pass

    elif OpSys == 'Darwin':
        try:
            winePath = []
            wine_path = Popen(['which', 'wine'], stdout=PIPE, shell=False, universal_newlines=True)#.communicate()[0]
            for stdout_line in iter(wine_path.stdout.readline, ''):
                winePath.append(stdout_line)
            if winePath != []:
                is_wine = Popen(['%s' % (winePath[0].strip('\n')), '--version'], stdout=PIPE, shell = False, universal_newlines=True)
            else:
                is_wine = Popen(['/usr/local/bin/wine','--version'], stdout=PIPE, shell = False, universal_newlines=True)

        except:
            print('wine could not be found on your system. resipy needs wine to run the inversion. You can install wine by running `brew install wine`.')

wineCheck()


#%% useful functions

# small useful function for reading and writing mesh.dat
# def readMeshDat(fname):
#     """Read mesh.dat or mesh3d.dat and returns elements, nodes, idirichlet.
#     """
#     with open(fname, 'r') as f:
#         x = f.readline().split()
#     numel = int(x[0])
#     nnodes = int(x[1])
#     idirichlet = int(x[2])
#     elems = np.genfromtxt(fname, skip_header=1, max_rows=numel)
#     if fname[-6:] == '3d.dat': # it's a 3D mesh
#         idirichlet = np.genfromtxt(fname, skip_header=numel+nnodes+1, skip_footer=0)[0]
#         skip_footer = 1
#     else:
#         skip_footer = 0
#     nodes = np.genfromtxt(fname, skip_header=numel+1, skip_footer=skip_footer)
#     return elems, nodes, idirichlet


# def writeMeshDat(fname, elems, nodes, extraHeader='', footer='1', idirichlet=1):
#     """Write mesh.dat/mesh3d.dat provided elements and nodes at least.
#     """
#     numel = len(elems)
#     nnodes = len(nodes)
#     threed = nodes.shape[1] == 4 # it's a 3D mesh
#     if threed is True:
#         extraHeader = '\t1\t0\t4'
#     with open(fname, 'w') as f:
#         f.write('{:.0f} {:.0f} {:.0f}{}\n'.format(numel, nnodes, idirichlet, extraHeader))
#     with open(fname, 'ab') as f:
#         np.savetxt(f, elems, fmt='%.0f')
#         if threed is True:
#             np.savetxt(f, nodes, fmt='%.0f %f %f %f')
#         else:
#             np.savetxt(f, nodes, fmt='%.0f %f %f')
#     if threed is True: # for 3D only
#         with open(fname, 'a') as f:
#             f.write(footer)

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)
    
    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)
    
    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

# distance matrix function for 2D (numpy based from https://stackoverflow.com/questions/22720864/efficiently-calculating-a-euclidean-distance-matrix-using-numpy)
def cdist(a):
    z = np.array([complex(x[0], x[1]) for x in a])
    return np.abs(z[...,np.newaxis]-z)



#%% main R2 class
class R2(object): # R2 master class instanciated by the GUI
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
        self.fmd = None # depth of investigation below the surface [in survey units]
        self.proc = None # where the process to run R2/cR2 will be
        self.zlim = None # zlim to plot the mesh by default (from max(elec, topo) to min(doi, elec))
        self.geom_input = {} # dictionnary used to create the mesh
        
        # attributes needed for independant error model for timelapse/batch inversion
        self.referenceMdl = False # is there a starting reference model already?
        self.fwdErrModel = False # is there is a impoforward modelling error already (due to the mesh)?
        self.errTyp = 'global'# type of error model to be used in batch and timelapse surveys


    def setBorehole(self, val=False):
        """Set all surveys in borehole type if `True` is passed.
        """
        self.iBorehole = val
        for s in self.surveys:
            s.iBorehole = val


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
                elec['buried'] = elec['buried'].astype(bool)
        if 'y' not in elec.columns:
            elec['y'] = 0 # all on same line by default
        if 'z' not in elec.columns:
            elec['z'] = 0 # all flat topo by default
        if 'remote' not in elec.columns:
            elec['remote'] = False # all non remote by default
        if 'buried' not in elec.columns:
            elec['buried'] = False # all surface elec by default
        if 'label' not in elec.columns:
            elec['label'] = (1 + np.arange(elec.shape[0])).astype(str) # all elec ordered and start at 1
        elec = elec.astype({'x':float, 'y':float, 'z':float, 'buried':bool, 'remote':bool, 'label':str})
        
        return elec
    

    def setElec(self, elec, elecList=None):
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
        """
        if elecList is None: # same electrode set shared by all surveys (most common case)
            ok = False
            elec = self._num2elec(elec)
            if self.elec is not None: # electrode already inferred when parsing data
                if self.iForward: # in case of forward modelling, changing the number of electrodes is allowed
                    ok = True
                else: # check intersection of labels
                    if len(self.surveys) > 0:
                        s1 = np.unique(elec['label'].values)
                        s2 = np.unique(self.surveys[0].df[['a','b','m','n']].values.flatten())
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
                # identification remote electrode
                remote_flags = [-9999999, -999999, -99999,-9999,-999,
                            9999999, 999999, 99999, 9999, 999] # values asssociated with remote electrodes
                iremote = np.in1d(elec['x'].values, remote_flags)
                iremote = np.isinf(elec[['x','y','z']].values).any(1) | iremote
                elec.loc[:, 'remote'] = iremote
                if np.sum(iremote) > 0:
                    print('Detected {:d} remote electrode.'.format(np.sum(iremote)))
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
            self.elec = self._num2elec(np.zeros((len(initElec),3)))
            for i, survey in enumerate(self.surveys):
                survey.elec = self._num2elec(elecList[i])
        
        if len(self.surveys) > 0:
            self.computeFineMeshDepth()


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

        # create save directory
        name = os.path.basename(fname)
        savedir = os.path.join(self.dirname, name)
        if os.path.exists(savedir):
            shutil.rmtree(savedir)
        os.mkdir(savedir)
        
        # add files to it
        self.mesh.write_vtk(os.path.join(savedir, 'mesh.vtk'))
        self.elec.to_csv(os.path.join(savedir, 'elec.csv'))
        for i, survey in enumerate(self.surveys):
            f = os.path.join(savedir, 'survey{:d}'.format(i))
            survey.df.to_csv(f + '-df.csv', index=False)
            survey.elec.to_csv(f + '-elec.csv')
            self.meshResults[i].write_vtk(f + '.vtk')

        # TODO add flags (borehole, timelapse)
        # TODO add all param
        
        # zip the directory, move it and clean
        with ZipFile(fname + '.resipy', 'w') as fz:
            for file in os.listdir(savedir):
                fz.write(os.path.join(savedir, file), file)
        shutil.rmtree(savedir)
        
        
        
    def loadProject(self, fname):
        """Load data from project file.
        """
        from zipfile import ZipFile, ZipInfo
        
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
        fs = [f for f in os.listdir(savedir) if f[:6] == 'survey']
        meshResults = []
        for i in range(len(fs)//2):
            f = os.path.join(savedir, 'survey{:d}'.format(i))
            df = pd.read_csv(f + '-df.csv')
            elec = pd.read_csv(f + '-elec.csv').values
            meshResults.append(mt.vtk_import(f + '.vtk'))
                

    def createSurvey(self, fname='', ftype='Syscal', info={}, spacing=None, 
                     parser=None, debug=True):
        """Read electrodes and quadrupoles data and return 
        a survey object.

        Parameters
        ----------
        fname : str
            Filename to be parsed.
        ftype : str, optional
            Type of file to be parsed. Either 'Syscal','ProtocolDC','Res2Dinv',
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
        """
        self.surveys.append(Survey(fname, ftype, spacing=spacing, parser=parser, debug=debug))
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
            
            
    def addData(self, **kwargs):
        """Adds data to the survey - used usually to add reciprocal datasets
        
        Parameters
        ----------
        **kwargs: Keyword arguments to be passed to Survey.addData()
        """
        
        self.surveys[0].addData(**kwargs)


    def createBatchSurvey(self, dirname, ftype='Syscal', info={}, spacing=None,
                          parser=None, isurveys=[], dump=None, debug=False):
        """Read multiples files from a folders (sorted by alphabetical order).

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
            self.createSurvey(f, ftype=ftype, parser=parser, spacing=spacing, debug=debug)
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
        self.bigSurvey = Survey(files[0], ftype=ftype, spacing=spacing, debug=False)
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
            df = df.append(df2, sort=True) # sort to silence the future warning if concatenation axis is not aligned
            c = c + df2.shape[0]
        self.bigSurvey.df = df.copy() # override it
        self.bigSurvey.dfOrigin = df.copy()
        self.bigSurvey.ndata = df.shape[0]



    def create3DSurvey(self, fname, lineSpacing=1, zigzag=False, ftype='Syscal',
                       name=None, parser=None):
        """Create a 3D survey based on 2D regularly spaced surveys.
        
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
            Type of the survey to choose which parser to use.
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

        surveys = []
        for fname in fnames:
            surveys.append(Survey(fname, ftype=ftype, parser=parser))
        survey0 = surveys[0]
        
        # check this is a regular grid
        nelec = survey0.elec.shape[0]
        for s in surveys:
            if s.elec.shape[0] != nelec:
                raise ValueError('Survey {:s} has {:d} electrodes while the first survey has {:d}.'
                                 'All surveys should have the same number of electrodes.'.format(s.name, s.elec.shape[0], nelec))
        # build global electrodes and merged dataframe
        elec = []
        dfs = []
        for i, s in enumerate(surveys):
            e = s.elec.copy()
            e.loc[:, 'y'] = i*lineSpacing
            prefix = '{:d} '.format(i+1)
            e.loc[:, 'label'] = prefix + e['label']
            elec.append(e)
            df = s.df.copy()
            df.loc[:,['a','b','m','n']] = prefix + df[['a','b','m','n']]
            dfs.append(df)
        elec = pd.concat(elec, axis=0, sort=False).reset_index(drop=True)
        dfm = pd.concat(dfs, axis=0, sort=False).reset_index(drop=True)
        
        survey0.elec = elec
        survey0.df = dfm
        survey0.name = '3Dfrom2Dlines' if name is None else name
        self.surveys= [survey0]
        self.elec = None
        self.setElec(elec)
        self.setBorehole(self.iBorehole)
        


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
        where all datasets must have the same number of quadrupoles.
        """
        print('Matching quadrupoles between surveys for difference inversion...', end='')
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
        x0 = cols2str(df0[['a','b','m','n']].values)
        icommon = np.ones(len(x0), dtype=bool)
        for df in dfs2[1:]:
            x = cols2str(df[['a','b','m','n']].values)
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
        if index == -1:
            for s in self.surveys:
                s.iselect = np.zeros(s.df.shape[0], dtype=bool)
                s.eselect = np.zeros(len(s.elec), dtype=bool)
            self.surveys[0].filterManual(ax=ax, **kwargs)
        else:
            self.surveys[index].filterManual(ax=ax, **kwargs)


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
            newlegend = []
            if ax is not None and len(self.surveys) > 1:
                legends = ax.get_legend().get_texts()
                for a in legends:
                    newlegend.append('{:s} {:s}'.format(s.name, a.get_text()))
                ax.get_legend().remove()
                ax.legend(newlegend, fontsize=10)
                ax.set_title(ax.get_title().split('\n')[0])
        else:
            self.surveys[index].fitErrorLin(ax=ax)


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
            newlegend = []
            if ax is not None and len(self.surveys) > 1:
                legends = ax.get_legend().get_texts()
                for a in legends:
                    newlegend.append('{:s} {:s}'.format(s.name, a.get_text()))
                ax.get_legend().remove()
                ax.legend(newlegend, fontsize=10)
                ax.set_title(ax.get_title().split('\n')[0])
        else:
            self.surveys[index].fitErrorPwl(ax=ax)


    def fitErrorLME(self, index=-1, ax=None, rpath=None, iplot=True):
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
            newlegend = []
            if ax is not None and len(self.surveys) > 1:
                legends = ax.get_legend().get_texts()
                for a in legends:
                    newlegend.append('{:s} {:s}'.format(s.name, a.get_text()))
                ax.get_legend().remove()
                ax.legend(newlegend, fontsize=10)
        else:
            self.surveys[index].fitErrorLME(ax=ax, rpath=rpath, iplot=iplot)


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
            newlegend = []
            if ax is not None and len(self.surveys) > 1:
                legends = ax.get_legend().get_texts()
                for a in legends:
                    newlegend.append('{:s} {:s}'.format(s.name, a.get_text()))
                ax.get_legend().remove()
                ax.legend(newlegend, fontsize=10)
                ax.set_title(ax.get_title().split('\n')[0])
        else:
            self.surveys[index].fitErrorPwlIP(ax=ax)


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
            newlegend = []
            if ax is not None and len(self.surveys) > 1:
                legends = ax.get_legend().get_texts()
                for a in legends:
                    newlegend.append('{:s} {:s}'.format(s.name, a.get_text()))
                ax.get_legend().remove()
                ax.legend(newlegend, fontsize=10)
                ax.set_title(ax.get_title().split('\n')[0])
        else:
            self.surveys[index].fitErrorParabolaIP(ax=ax)



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
        -----------
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
        -----------
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
            if len(self.surveys) > 0:
                lookupDict = dict(zip(dfelec['label'], np.arange(dfelec.shape[0])))
                array = self.surveys[0].df[['a','b','m','n']].replace(lookupDict).values
                maxDist = np.max(np.abs(elec[array[:,0]-np.min(array[:,0]),0] - elec[array[:,2]-np.min(array[:,2]),0])) # max dipole separation
                self.fmd = (1/3)*maxDist
            else: # if it's a forward model for instance
                self.fmd = (1/3)*(np.max(elec[:,0]) - np.min(elec[:,0]))

        else: # for 3D survey
            dist = np.zeros((len(elec), len(elec)))
            for i, el1 in enumerate(elec):
                dist[:,i] = np.sqrt(np.sum((el1[None,:] - elec)**2, axis=1))
            self.fmd = (1/3)*np.max(dist)

        if self.elec['buried'].sum() > 0:
            # catch where buried electrodes are present as the fmd needs adjusting in this case 
            if self.elec['buried'].sum() == self.elec.shape[0]:
                # if all buried, we assume surface at 0 m
                self.fmd = np.abs(0 - np.min(elec[:,2])) + 1
            else: # surface given by max z elec
                self.fmd = np.abs(np.max(elec[:,2])  - np.min(elec[:,2])) + (0.5*self.fmd)
        print('Fine Mesh Depth (relative to the surface): {:.2f} m'.format(self.fmd))



    def createMesh(self, typ='default', buried=None, surface=None, cl_factor=2,
                   cl=-1, dump=None, res0=100, show_output=False, fmd=None,
                   remote=None, refine=0, **kwargs):
        """Create a mesh.

        Parameters
        ----------
        typ : str, optional
            Type of mesh. Either 'quad' or 'trian' in the case of 2d surveys.
            By default, 'trian' is chosen for 2D and 'tetra' is used for 
            3D surveys, but 'prism' can be used for column type experiments. 
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
            length as `R2.elec`.
        refine : int, optional
            Number times the mesh will be refined. Refinement split the triangles
            or the tetrahedra but keep the same number of parameter for the inversion.
            This helps having a more accurate forward response and a faster inversion
            (as the number of elements does not increase). Only available for
            triangles or tetrahedral mesh.
        kwargs : -
            Keyword arguments to be passed to mesh generation schemes
        """        
        if dump is None:
            if show_output:
                def dump(x):
                    print(x, end='')
            else:
                def dump(x):
                    pass
        self.meshParams = {'typ':typ, 'buried':buried, 'surface':surface,
                           'cl_factor':cl_factor, 'cl':cl, 'dump':dump,
                           'res0': res0, 'show_output':show_output}#, 'dfm':doi}
        if kwargs is not None:
            self.meshParams.update(kwargs)

        if typ == 'default':
            if self.typ == 'R2' or self.typ == 'cR2': # it's a 2D mesh
                typ = 'trian'
            else:
                typ = 'tetra'

        #check if remote electrodes present?
        if (self.elec['remote'].sum() > 0) & (typ == 'quad'):
            dump('remote electrode is not supported in quadrilateral mesh for now, please use triangular mesh instead.')
            return

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
        elecLabels = self.elec['label'].values
        
        # estimate depth of fine mesh
        if fmd is None:
            self.computeFineMeshDepth()
        else:
            self.fmd = fmd

        if typ == 'quad':
            print('Creating quadrilateral mesh...', end='')
            surface_x = surface[:,0] if surface is not None else None
            surface_z = surface[:,2] if surface is not None else None
            mesh,meshx,meshy,topo,e_nodes = mt.quad_mesh(elec_x,elec_z,list(elec_type),
                                                         surface_x=surface_x, surface_z=surface_z,
                                                         **kwargs)   #generate quad mesh
            self.param['mesh_type'] = 6
            e_nodes = np.array(mesh.e_nodes) + 1 # +1 because of indexing staring at 0 in python
            self.param['node_elec'] = [elecLabels, e_nodes.astype(int)]

            if 'regions' in self.param: # allow to create a new mesh then rerun inversion
                del self.param['regions']
            if 'num_regions' in self.param:
                del self.param['num_regions']
        elif typ == 'trian' or typ == 'tetra' or typ=='prism':
            geom_input = {}

            if surface is not None:
                if surface.shape[1] == 2:
                    geom_input['surface'] = [surface[:,0], surface[:,1]]
                else:
                    geom_input['surface'] = [surface[:,0], surface[:,2]]

            if 'geom_input' in kwargs:
                geom_input.update(kwargs['geom_input'])
                kwargs.pop('geom_input')

            whole_space = False
            if buried is not None:
                if np.sum(buried) == len(buried) and surface is None:
                    # all electrodes buried and no surface given
                    whole_space = True

            elec_type = elec_type.tolist()

            with cd(self.dirname):#change to working directory so that mesh files written in working directory
                if typ == 'trian':
                    print('Creating triangular mesh...', end='')
                    mesh = mt.tri_mesh(elec_x,elec_z,elec_type,geom_input,
                                 path=os.path.join(self.apiPath, 'exe'),
                                 cl_factor=cl_factor,
                                 cl=cl, dump=dump, show_output=show_output,
                                 fmd=self.fmd, whole_space=whole_space,
                                 **kwargs)
                if typ == 'tetra': # TODO add buried
                    print('Creating tetrahedral mesh...', end='')    
                    if cl == -1:
                        dist = cdist(self.elec[~self.elec['remote']][['x','y']].values)/2 # half the minimal electrode distance
                        cl = np.min(dist[dist != 0])
                    mesh = mt.tetra_mesh(elec_x, elec_y, elec_z,elec_type,
                                 path=os.path.join(self.apiPath, 'exe'),
                                 surface_refinement=surface,
                                 cl_factor=cl_factor,
                                 cl=cl, dump=dump, show_output=show_output,
                                 fmd=self.fmd, whole_space=whole_space,
                                 **kwargs)
                if typ=='prism':
                    print('Creating prism mesh...', end='')
                    mesh = mt.prism_mesh(elec_x, elec_y, elec_z,
                                         path=os.path.join(self.apiPath, 'exe'),
                                         cl=cl, dump=dump, show_output=show_output,
                                         **kwargs)
                    self.param['num_xz_poly'] = 0
                    

            # mesh refinement
            if (typ == 'trian') | (typ == 'tetra'):
                for l in range(refine):
                    print('Refining...', end='')
                    mesh.refine()
            
            self.param['mesh_type'] = 3
            e_nodes = np.array(mesh.e_nodes) + 1 # +1 because of indexing staring at 0 in python
            self.param['node_elec'] = [elecLabels, e_nodes.astype(int)]

        self.mesh = mesh
        self.param['mesh'] = mesh
        self.param['num_regions'] = 0
        self.param['res0File'] = 'res0.dat'
        numel = self.mesh.num_elms
        self.mesh.add_attribute(np.ones(numel)*res0, 'res0') # default starting resisivity [Ohm.m]
        self.mesh.add_attribute(np.ones(numel)*0, 'phase0') # default starting phase [mrad]
        self.mesh.add_attribute(np.ones(numel, dtype=int), 'zones')
        self.mesh.add_attribute(np.zeros(numel, dtype=float), 'iter')
        self.mesh.add_attribute(np.arange(numel)+1,'param') # param = 0 if fixed

        # define zlim
        if surface is not None:
            zlimTop = np.max([np.max(elec_z), np.max(surface[:,-1])])
        else:
            if all(self.elec['buried']): # if all buried we assume surface at 0 m
                zlimTop = 0
            else:
                zlimTop = np.max(elec_z)
        zlimBot = np.min(elec_z)-self.fmd # if fmd is correct (defined as positive number
        # from surface then it is as well the zlimMin)
        self.zlim = [zlimBot, zlimTop]
        
        self._computePolyTable()
        print('done')
        
        
    def _computePolyTable(self):
        # define num_xz_poly or num_xy_poly
        elec = self.elec[~self.elec['remote']][['x','y','z']].values
        elec_x, elec_y, elec_z = elec[:,0], elec[:,1], elec[:,2]
        if (self.typ == 'R2') | (self.typ == 'cR2'):
            self.param['num_xz_poly'] = 5
            if self.elec['buried'].sum() == self.elec.shape[0]:
                zmax = 0
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
        else:
            self.param['num_xy_poly'] = 5
            xmin, xmax = np.min(elec_x), np.max(elec_x)
            ymin, ymax = np.min(elec_y), np.max(elec_y)
            zmin, zmax = np.min(elec_z)-self.fmd, np.max(elec_z)
            xz_poly_table = np.array([
            [xmin, ymax],
            [xmax, ymax],
            [xmax, ymin],
            [xmin, ymin],
            [xmin, ymax]])
            self.param['zmin'] = zmin
            self.param['zmax'] = zmax
            self.param['xy_poly_table'] = xz_poly_table

        

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
        self.mesh = mt.custom_mesh_import(file_path, node_pos=node_pos, 
                                          order_nodes=order_nodes)
        if elec is not None:
            self.mesh.moveElecNodes(elec[:,0], elec[:,1], elec[:,2])

        #add the electrodes to the R2 class
        if elec is not None or node_pos is not None: # then electrode positions should be known
            self.setElec(np.array((self.mesh.elec_x, self.mesh.elec_y, self.mesh.elec_z)).T)
        else:
            try:
                elec = self.elec[['x','y','z']].values
                self.mesh.moveElecNodes(elec[:,0],elec[:,1],elec[:,2])
            except AttributeError:
                warnings.warn("No electrode nodes associated with mesh! Electrode positions are unknown!")

        #R2 class mesh handling
        e_nodes = np.array(self.mesh.e_nodes) + 1 # +1 because of indexing staring at 0 in python
        self.param['mesh'] = self.mesh
        # if mesh_type == 'quad':
        #     self.param['mesh_type'] = 6
        #     colx = self.mesh.quadMeshNp() # convert nodes into column indexes
        #     self.param['node_elec'] = np.c_[1+np.arange(len(e_nodes)), np.array(colx), np.ones((len(e_nodes,1)))].astype(int)
            #will only work for assuming electrodes are a surface array
        if self.mesh.type2VertsNo() == 4:
            if flag_3D:
                self.param['mesh_type'] = 4 # tetra mesh
            else:
                self.param['mesh_type'] = 6 # general quad mesh
        else:
            if flag_3D:
                self.param['mesh_type'] = 6 # prism mesh 
            else:
                self.param['mesh_type'] = 3 # triangular mesh
        self.param['node_elec'] = [self.elec['label'].values, e_nodes.astype(int)]

        # checking
        if len(np.unique(e_nodes)) < len(e_nodes):
            raise ValueError('Some electrodes are positionned on the same nodes !')
        
        # make regions continuous
        regions = self.mesh.attr_cache['region']
        uregions = np.unique(regions)
        iregions = np.arange(len(uregions)) + 1
        dico = dict(zip(uregions, iregions))
        self.mesh.attr_cache['region'] = [dico[a] for a in regions]
        
        self.param['num_regions'] = 0
        self.param['res0File'] = 'res0.dat'
        numel = self.mesh.num_elms
        self.mesh.add_attribute(np.ones(numel)*res0, 'res0') # default starting resisivity [Ohm.m]
        self.mesh.add_attribute(np.ones(numel)*0, 'phase0') # default starting phase [mrad]
        self.mesh.add_attribute(np.ones(numel, dtype=int), 'zones')
        self.mesh.add_attribute(np.zeros(numel, dtype=float), 'iter')

        name = 'mesh.dat'
        if self.typ == 'R3t' or self.typ == 'cR3t':
            name = 'mesh3d.dat'
        file_path = os.path.join(self.dirname, name)
        self.mesh.write_dat(file_path)

        # define zlim
        if self.fmd is None:
            self.computeFineMeshDepth()
        zlimMax = self.elec['z'].max()
        zlimMin = self.elec['z'].min() - self.fmd
        self.zlim = [zlimMin, zlimMax]
        
        self._computePolyTable()



    def showMesh(self, ax=None, **kwargs):
        """Display the mesh.
        """
        if self.mesh is None:
            raise Exception('Mesh undefined')
        else:
            if 'zlim' not in kwargs.keys():
                kwargs['zlim'] = self.zlim
            if 'color_map' not in kwargs.keys():
                kwargs['color_map'] = 'Spectral'
            if 'attr' not in kwargs.keys():
                kwargs['attr'] = 'region' # this will print regions
            if 'color_bar' not in kwargs.keys():
                if np.unique(np.array(self.mesh.attr_cache['region'])).shape[0] > 1:
                    kwargs['color_bar'] = True # show colorbar for multiple regions
                else:
                    kwargs['color_bar'] = False
        
            self.mesh.show(ax=ax, **kwargs)


    def write2in(self, param={}):
        """Create configuration file for inversion. Write mesh.dat and res0.dat.

        Parameters
        ----------
        param : dict
            Dictionnary of parameters and values for the inversion settings.
        """
        typ = self.typ
        if (self.err is True) and ('a_wgt' not in self.param):
            self.param['a_wgt'] = 0
            self.param['b_wgt'] = 0
        elif (typ == 'R2') or (typ == 'R3t'): # DC case
            if 'a_wgt' not in self.param:
                self.param['a_wgt'] = 0.01
            if 'b_wgt' not in self.param:
                self.param['b_wgt'] = 0.02
        elif (typ == 'cR2') | (typ == 'cR3t'): # IP case
            if 'a_wgt' not in self.param:
                self.param['a_wgt'] = 0.02 # variance for magnitude (no more offset)
            if 'b_wgt' not in self.param:
                self.param['b_wgt'] = 2 # mrad

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
            else: # default DC case as timelapse not supported for IP yet
                if 'a_wgt' not in param:#this allows previously assigned values to be
                    param['a_wgt'] = 0.01 # written to the reference.in config file
                if 'b_wgt' not in param:
                    param['b_wgt'] = 0.02
            param['reg_mode'] = 0 # set by default in ui.py too
            param['res0File'] = 'res0.dat'
            if (self.typ == 'R2') or (self.typ == 'cR2'):
                param['num_xz_poly'] = 0
            else:
                param['num_xy_poly'] = 0
                param['inverse_type'] = 0 # normal regularisation
                param['zmin'] = np.min(self.mesh.node_z) - 10 # we want to keep the whole mesh for background regularisation
                param['zmax'] = np.max(self.mesh.node_z) + 10
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
        ifixed = np.array(self.mesh.attr_cache['param']) == 0
        if np.sum(ifixed) > 0: # fixed element need to be at the end
            self.mesh.orderElem()
        name = 'mesh.dat'
        if (typ == 'R3t') | (typ == 'cR3t'):
            name = 'mesh3d.dat'
        self.mesh.write_dat(os.path.join(self.dirname, name))
        
        # write the res0.dat needed for starting resistivity
        if self.iForward is True: # we will invert results from forward
            # inversion so we need to start from a homogeneous model
            print('All non fixed parameters reset to 100 Ohm.m and 0 mrad, '
                  'as the survey to be inverted is from a forward model.')
            ifixed = self.mesh.attr_cache['param'] == 0
            res0 = np.array(self.mesh.attr_cache['res0'])
            phase0 = np.array(self.mesh.attr_cache['phase0'])
            res0f = res0.copy()
            phase0f = phase0.copy()
            res0f[~ifixed] = 100
            phase0f[~ifixed] = 0
            self.mesh.attr_cache['res0'] = list(res0f)
            self.mesh.attr_cache['phase0'] = list(phase0f)

        if (self.typ == 'cR2') or (self.typ == 'cR3t'):
            r = np.array(self.mesh.attr_cache['res0'])
            phase = np.array(self.mesh.attr_cache['phase0'])
            centroids = np.array(self.mesh.elm_centre).T
            centroids2 = centroids[:,[0,2]] if self.typ[-1] != 't' else centroids
            x = np.c_[centroids2,
                      r,
                      phase, # mrad
                      np.log10(r),
                      np.log10(np.cos(-phase/1000)/np.log10(r)), #log10(real conductivity)
                      np.log10(np.sin(-phase/1000)/np.log10(r))] #log10(imaginary conductivity)
            np.savetxt(os.path.join(self.dirname, 'res0.dat'), x)
        else:
            self.mesh.write_attr('res0', os.path.join(self.dirname, 'res0.dat'))
        
        if self.iForward: # restore initial res0 and phase0 so that user can 
        # rerun the forward model with a different sequence for instance
            self.mesh.attr_cache['res0'] = list(res0)
            self.mesh.attr_cache['phase0'] = list(phase0)
            

    def write2protocol(self, err=None, errTot=False, **kwargs):
        """Write a protocol.dat file for the inversion code.

        Parameters
        ----------
        err : bool, optional
            If `True` error columns will be written in protocol.dat provided
            an error model has been fitted or errors have been imported.
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
            else:
                indexes = [None]*len(self.surveys)
            # a bit simplistic but assign error to all based on Transfer resistance
            # let's assume it's False all the time for now
            content = ''
            df0 = self.surveys[0].df[['a','b','m','n','resist','recipMean']]
            df0 = df0.rename(columns={'resist':'resist0', 'recipMean':'recipMean0'})
            for i, s in enumerate(self.surveys):
                if 'resist0' in s.df.columns:
                    s.df = s.df.drop('resist0', axis=1)
                if 'recipMean0' in s.df.columns:
                    s.df = s.df.drop('recipMean0', axis=1)
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
                protocol = s.write2protocol('', err=err, errTot=errTot, res0=res0Bool,
                                            ip=False, # no IP timelapse possible for now
                                            isubset=indexes[i], threed=threed)
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
                    s.write2protocol(os.path.join(refdir, 'protocol.dat'), err=err, threed=threed) # no subset for background, just use all
                else:
                    content = content + str(protocol.shape[0]) + '\n'
                    content = content + protocol.to_csv(sep='\t', header=False, index=False)

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

        # kill management
        class ProcsManagement(object): # little class to handle the kill
            def __init__(self, r2object):
                self.r2 = r2object
            def kill(self):
                print('killing...')
                self.r2.irunParallel2 = False # this will end the infinite loop
                procs = self.r2.procs # and kill the running processes
                for p in procs:
                    p.terminate()
                print('all done!')

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
        dump('\r{:.0f}/{:.0f} inversions completed'.format(c, len(wds2)))
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
                    # TODO get RMS and iteration number here ?
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
        # TODO should now be consistent
        # if self.typ=='R3t' or self.typ=='cR3t':
            # toRename = ['.dat', '.vtk', '.err', '.sen', '_diffres.dat']
        # else:
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



    def invert(self, param={}, iplot=False, dump=None, modErr=False,
               parallel=False, iMoveElec=False, ncores=None,
               rmDirTree=True, modelDOI=False):
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
        """
        if dump is None:
            def dump(x):
                print(x, end='')
                
        # clean meshResults list
        self.meshResults = []
        
        # create mesh if not already done
        if 'mesh' not in self.param:
            self.createMesh()
            
        # run Oldenburg and Li DOI estimation
        if modelDOI is True:
            sensScaled = self.modelDOI(dump=dump)

        # compute modelling error if selected
        if modErr is True and self.fwdErrModel is False: #check no error model exists
            dump('Computing error model... ')
            self.computeModelError()
            dump('done!\n')
            errTot = True
        elif modErr is True and self.fwdErrModel:
            # avoid computing error model again if it has already been run.
            errTot = True
        else:
            errTot = False

        # write configuration file
        dump('Writing .in file and protocol.dat... ')
        self.write2in(param=param) # R2.in
        self.write2protocol(errTot=errTot) # protocol.dat
        dump('done!\n')

        # runs inversion
        if self.iTimeLapse == True and self.referenceMdl==False:
            dump('------------ INVERTING REFERENCE SURVEY ---------------\n')
            refdir = os.path.join(self.dirname, 'ref')
            shutil.move(os.path.join(self.dirname,'res0.dat'),
                        os.path.join(refdir, 'res0.dat'))
            self.write2in(param=param)
            self.runR2(refdir, dump=dump) # this line actually runs R2
            shutil.copy(os.path.join(refdir, 'f001_res.dat'),
                        os.path.join(self.dirname, 'Start_res.dat'))
        elif self.iTimeLapse == True and self.referenceMdl==True:
            print('Note: Skipping reference inversion, as reference model has already been assigned')

        dump('--------------------- MAIN INVERSION ------------------\n')
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
                    m.attr_cache['doiSens'] = sensScaled
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
                
        # backup for normal inversion (0 : original, 1 : normal background, 2: background *10)
        res0 = np.array(self.mesh.attr_cache['res0'])
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
            iselect = path.contains_points(np.c_[self.mesh.elm_centre[0], self.mesh.elm_centre[2]])
        else:
            iselect = np.ones(len(self.mesh.elm_centre[0]), dtype=bool)
            
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
        self.mesh.attr_cache['res0b'] = list(res1)
        self.mesh.write_attr('res0b', os.path.join(self.dirname,'res0.dat'))
        self.runR2(dump=dump) # re-run inversion
        self.getResults()
        mesh1 = self.meshResults[0]
        cleandir()
        
        # run second background constrained inversion
        dump('===== modelDOI: Running background constrained inversion with initial resistivity * 10 =====\n')
        res2 = res0 * 10
        self.mesh.attr_cache['res0b'] = list(res2)
        self.mesh.write_attr('res0b', os.path.join(self.dirname,'res0.dat'))
        self.runR2(dump=dump) # re-run inversion
        self.getResults()
        mesh2 = self.meshResults[0]
        cleandir()
        os.remove(os.path.join(self.dirname, 'R2.in'))
        
        # sensitivity = difference between final inversion / difference init values
        res_names = np.array(['Resistivity','Resistivity(Ohm-m)','Resistivity(ohm.m)'])
        res_name = res_names[np.in1d(res_names, list(self.meshResults[0].attr_cache.keys()))][0]
        invValues1 = np.array(mesh1.attr_cache[res_name])
        invValues2 = np.array(mesh2.attr_cache[res_name])
        sens = (invValues1 - invValues2)/(res1[iselect]-res2[iselect])
        sensScaled = np.abs(sens)
#        mesh0.attr_cache['doiSens'] = sensScaled # add attribute to original mesh
        self.doiComputed = True
        
        # restore
        self.meshResults = []
        self.param = param0
        self.typ = typ0
        self.surveys = surveys0
        self.iTimeLapse = iTimeLapse0
        # .in and protocol will be written again in R2.invert()
        
        return sensScaled
        
    
    
    def _clipContour(self, ax, collections, cropMaxDepth=False):
        """Clip contours using mesh bound and surface if available.
        
        Parameters
        ----------
        ax : matplotlib.Axes
            Axis.
        collections : matplotlib.collections
            Matplotlib collection.
        cropMaxDepth : bool, optional
            If 'True', area below fmd will be cropped out.
        """
        # mask outer region
        xmin = np.min(self.mesh.node_x)
        xmax = np.max(self.mesh.node_x)
        zmin = np.min(self.mesh.node_z)
        zmax = np.max(self.mesh.node_z)
        if self.mesh.surface is not None:
            xsurf, zsurf = self.mesh.surface[:,0], self.mesh.surface[:,1]
            if cropMaxDepth and self.fmd is not None:
                xfmd, zfmd = self.mesh.surface[:,0][::-1], self.mesh.surface[:,1][::-1] - self.fmd
                verts = np.c_[np.r_[xmin, xmin, xsurf, xmax, xmax, xfmd, xmin],
                              np.r_[zmin, zmax, zsurf, zmax, zmin, zfmd, zmin]]
            else:
                verts = np.c_[np.r_[xmin, xmin, xsurf, xmax, xmax, xmin],
                              np.r_[zmin, zmax, zsurf, zmax, zmin, zmin]]
                
        else:
            verts = np.c_[np.r_[xmin, xmin, xmax, xmax, xmin],
                          np.r_[zmin, zmax, zmax, zmin, zmin]]                
        # cliping using a patch (https://stackoverflow.com/questions/25688573/matplotlib-set-clip-path-requires-patch-to-be-plotted)
        poly_codes = [mpath.Path.MOVETO] + (len(verts) - 2) * [mpath.Path.LINETO] + [mpath.Path.CLOSEPOLY]
        path = mpath.Path(verts, poly_codes)
        patch = mpatches.PathPatch(path, facecolor='none', edgecolor='none')
        ax.add_patch(patch) # need to add so it knows the transform
        for col in collections:
            col.set_clip_path(patch)
            
        

    def showResults(self, index=0, ax=None, edge_color='none', attr='',
                    sens=True, color_map='viridis', zlim=None, clabel=None,
                    doi=False, doiSens=False, contour=False, cropMaxDepth=True,
                    **kwargs):
        """Show the inverteds section.

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
        doi : bool, optional
            If True, it will draw a dotted red line corresponding to 0.02 from the
            Oldenburg and Li method. Note that `R2.modeDOI()` needs to be run
            for that.
        doiSens : bool, optional
            If True, it will draw a dashed line corresponding to 0.001 of the maximum
            of the log10 sensitivity.
        contour : bool, optional
            If True, contours will be plotted.
        """
        if len(self.meshResults) == 0:
            self.getResults()
        if (attr == '') & (self.typ[0] != 'c'):
            attr = 'Resistivity(log10)'
        if (attr == '') & (self.typ[0] == 'c'):
            attr = 'Sigma_real(log10)'
        keys = list(self.meshResults[index].attr_cache.keys())
        if attr not in keys:
            attr = keys[3]
            print('Attribute not found, revert to {:s}'.format(attr))
        if len(self.meshResults) > 0:
            mesh = self.meshResults[index]
            if self.typ[-1] == '2': # 2D case
                if zlim is None:
                    zlim = self.zlim
                mesh.show(ax=ax, edge_color=edge_color,
                            attr=attr, sens=sens, color_map=color_map,
                            zlim=zlim, clabel=clabel, contour=contour, **kwargs)
                if doi is True: # DOI based on Oldenburg and Li
                    if self.doiComputed is True: 
                        z = np.array(mesh.attr_cache['doiSens'])
                        levels = [0.2]
                        linestyle = ':'
                    else:
                        raise ValueError('Rerun the inversion with `modelDOI=True` first or use `doiSens`.')
                if doiSens is True: # DOI based on log10(sensitivity)
                    if 'Sensitivity(log10)' in mesh.attr_cache.keys():
                        z = np.array(mesh.attr_cache['Sensitivity(log10)'])
                        levels=[np.log10(0.001*(10**np.nanmax(z)))]
                        linestyle = '--'
                    else:
                        doiSens = False
                if doi is True or doiSens is True:
                    xc, yc = np.array(mesh.elm_centre[0]), np.array(mesh.elm_centre[2])
                    triang = tri.Triangulation(xc, yc)
                    cont = mesh.ax.tricontour(triang, z, levels=levels, colors='k', linestyles=linestyle)
                    self._clipContour(mesh.ax, cont.collections)
                colls = mesh.cax.collections if contour == True else [mesh.cax]
                self._clipContour(mesh.ax, colls, cropMaxDepth=cropMaxDepth)
            else: # 3D case
                if zlim is None:
                    zlim = [np.min(mesh.node_z), np.max(mesh.node_z)]
                if cropMaxDepth and self.fmd is not None:
                    zlim[0] = self.elec['z'].min() - self.fmd
                mesh.show(ax=ax, edge_color=edge_color,
                        attr=attr, color_map=color_map, clabel=clabel,
                        zlim=zlim, **kwargs)
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
        idone = 0
        ifailed = 0
        self.meshResults = [] # make sure we empty the list first
        if self.iTimeLapse == True:
            fname = os.path.join(dirname, 'ref', 'f001_res.vtk')
            mesh0 = mt.vtk_import(fname, order_nodes=False)
            mesh0.mesh_title = self.surveys[0].name
            elec = self.surveys[0].elec.copy()
            ie = ~elec['remote'].values
            mesh0.elec_x = elec[ie]['x'].values
            mesh0.elec_y = elec[ie]['y'].values
            mesh0.elec_z = elec[ie]['z'].values
            mesh0.surface = self.mesh.surface
            self.meshResults.append(mesh0)
            idone += 1
        if self.iForward is True:
            initMesh = mt.vtk_import(os.path.join(dirname, 'fwd','forward_model.vtk'), order_nodes=False)
            initMesh.elec_x = self.elec['x'].values
            initMesh.elec_y = self.elec['y'].values
            initMesh.elec_z = self.elec['z'].values
            initMesh.surface = self.mesh.surface
            initMesh.mesh_title = 'Initial Model'
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
                    ie = ~elec['remote'].values
                    mesh.elec_x = elec[ie]['x'].values
                    mesh.elec_y = elec[ie]['y'].values
                    mesh.elec_z = elec[ie]['z'].values
                    mesh.surface = self.mesh.surface
                    self.meshResults.append(mesh)
                    idone += 1
                except Exception:
                    ifailed += 1
                    # if inversion fails in time-lapse it's that the initial
                    # model is good enough to explain the data (a_wgt/b_wgt
                    # error too low) so we can replace it by the initial model
                    if self.iTimeLapse:
                        self.meshResults.append(mesh0) # TODO not sure
                print('\r{:d}/{:d} results parsed ({:d} ok; {:d} failed)'.format(
                    j+1, len(self.surveys), idone, ifailed), end='')
            else:
                pass
                #break
        print('')

        # compute conductivity in mS/m
        res_names = np.array(['Resistivity','Resistivity(Ohm-m)','Resistivity(ohm.m)', 'Magnitude(ohm.m)'])
        for mesh in self.meshResults:
            res_name = res_names[np.in1d(res_names, list(mesh.attr_cache.keys()))][0]
            mesh.attr_cache['Conductivity(mS/m)'] = 1000/np.array(mesh.attr_cache[res_name])
        if self.typ[0] == 'c' and self.surveys[0].kFactor != 1: # if kFactor is 1 then probably phase is provided and we shouldn't estimate chargeability
            for mesh in self.meshResults:
                mesh.attr_cache['Chargeability(mV/V)'] = np.array(mesh.attr_cache['Phase(mrad)'])/-self.surveys[0].kFactor
        # compute difference in percent in case of reg_mode == 1
        if (self.iTimeLapse is True) and (self.param['reg_mode'] == 1):
            try:
                self.computeDiff()
            except Exception as e:
                print('failed to compute difference: ', e)
                pass
    
    
    
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
            dfall = pd.read_csv(os.path.join(invdir, 'protocol.dat'),
                                sep='\t', header=None, engine='python').reset_index()
            idf = list(np.where(np.isnan(dfall[dfall.columns[-1]].values))[0])
            idf.append(len(dfall))
            dfs = [dfall.loc[idf[i]:idf[i+1]-1,:] for i in range(len(idf)-1)]
    
            # writing all protocol.dat
            files = []
            for i, df in enumerate(dfs):
                outputname = os.path.join(self.dirname, 'protocol_{:03d}.dat'.format(i))
                files.append(outputname)
                if self.typ[-1] == 't':
                    df2 = df[1:].astype({'level_0':int, 'level_1':int, 'level_2':int,
                                         'level_3':int, 'level_4':int, 'level_5':int,
                                         'level_6':int, 'level_7':int, 'level_8':int})
                    df2 = df2.drop('level_10', axis=1) # discard res0 as parser doesn't support it
                else:
                    df2 = df[1:].astype({'level_0':int, 'level_1':int, 'level_2':int,
                                         'level_3':int, 'level_4':int})
                    df2 = df2.drop('level_6', axis=1) # discard res0

                with open(outputname, 'w') as f:
                    f.write('{:d}\n'.format(df['level_0'].values[0]))
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
        with open(fname, 'r') as f:
            lines = f.readlines()
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
                        name = self.surveys[idataset-1].name
                    else:
                        name = 'dataset{:03.0f}'.format(idataset)
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
#        selector = SelectPoints(ax, np.array(self.mesh.elm_centre).T[:,[0,2]],
#                                typ='poly', iplot=iplot) # LIMITED FOR 2D case
#        selector.setVertices(xz)
#        selector.getPointsInside()
#        idx = selector.iselect

        centroids = np.array(self.mesh.elm_centre).T[:,[0,2]]
        path = mpath.Path(np.array(xz))
        idx = path.contains_points(centroids)

        region = np.array(self.mesh.attr_cache['region'])
        regid = np.max(region) + 1
        region[idx] = regid
        self.mesh.attr_cache['region'] = region
        resist0 = np.array(self.mesh.attr_cache['res0'])
        resist0[idx] = res0
        self.mesh.attr_cache['res0'] = resist0
        phase = np.array(self.mesh.attr_cache['phase0'])
        phase[idx] = phase0
        self.mesh.attr_cache['phase0'] = phase

        # define zone
        if blocky is True:
            zones = np.array(self.mesh.attr_cache['zones'])
            zones[idx] = regid
            self.mesh.attr_cache['zones'] = zones

        # define fixed area
        if fixed is True:
            paramFixed = np.array(self.mesh.attr_cache['param'])
            paramFixed[idx] = 0
            self.mesh.attr_cache['param'] = list(paramFixed)
            print('fixed {:d} elements'.format(np.sum(paramFixed == 0)))

        if iplot is True:
            self.showMesh()


    def resetRegions(self):
        """Just reset all regions already drawn. Shouldn't be needed as
        the `self.runR2()` automatically use a homogenous model when starting
        for inversion. The only purpose of this is to use an inhomogeous
        starting model to invert data from forward modelling.
        """
        self.mesh.attr_cache['region'] = np.ones(self.mesh.num_elms)
        self.mesh.attr_cache['res0'] = np.ones(self.mesh.num_elms)*100 # set back as default


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

        def callback(idx):
            region = np.array(self.mesh.attr_cache['region']).copy()
            regid = np.max(region) + 1
            print('nb elements selected:', np.sum(idx), 'in region', regid)
            region[idx] = regid
            self.mesh.attr_cache['region'] = region
            self.mesh.draw(attr='region')
            if addAction is not None:
                addAction()
        self.mesh.atribute_title = 'Regions'
        self.mesh.show(attr='region', ax=ax, zlim=self.zlim)
        # we need to assign a selector to self otherwise it's not used
        self.selector = SelectPoints(ax, np.array(self.mesh.elm_centre).T[:,[0,2]],
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

        self.geom_input = {}
        ax.plot(self.elec['x'], self.elec['z'], 'ko', label='electrode')
        ax.set_ylim([np.min(self.elec['z']) - self.fmd, np.max(self.elec['z'])])
        ax.set_xlim(np.min(self.elec['x']), np.max(self.elec['x']))
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('Elevation [m]')
        def callback():
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
        regionValues : dict, optional
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
        regions = np.array(self.mesh.attr_cache['region'])
        res0 = np.array(self.mesh.attr_cache['res0']).copy()
        for key in regionValues.keys():
            idx = regions == key
            res0[idx] = regionValues[key]
        self.mesh.attr_cache['res0'] = res0
        print('regionValues:',regionValues)

        zones = np.array(self.mesh.attr_cache['zones']).copy()
        for key in zoneValues.keys():
            idx = regions == key
            zones[idx] = zoneValues[key]
        self.mesh.attr_cache['zones'] = zones

        # fixed = np.array(self.mesh.attr_cache['param']).copy()
        fixed = np.arange(self.mesh.num_elms)+1
        for key in fixedValues.keys():
            idx = regions == key
            if fixedValues[key] == True:
                fixed[idx] = 0
        self.mesh.attr_cache['param'] = fixed

        phase0 = np.array(self.mesh.attr_cache['phase0']).copy()
        for key in ipValues.keys():
            idx = regions == key
            phase0[idx] = ipValues[key]
        self.mesh.attr_cache['phase0'] = phase0


    def setRefModel(self, res0):
        """Set the reference model according to a previous inversion, avoids
        the need to invert reference model again for timelapse workflows.
        In contrast to `R2.setStartingRes()` which assign resistivity to group
        of elements, this method requires a vector of the same length as the 
        number of elements. This enables, notably to manually perform consecutive
        background constrained inversion.

        Parameters
        -------------
        res0: array like
            Array of resistivity values, ideally from a previous inversion. The
            length of this array should be the same as the number of elements.
        """
        try:
            self.mesh.add_attribute(res0,'res0')
        except AttributeError:
            print('Cant set reference model without first assigning/creating a mesh')
            return
        self.param['reg_mode'] = 1 # ensure inversion is background regularised
        if self.typ[-1] =='t':
            self.param['inverse_type']=1
        self.param['res0File'] = 'Start_res.dat'
        self.param['num_regions'] = 0
        self.mesh.write_attr('res0',os.path.join(self.dirname,'Start_res.dat'))
        self.referenceMdl = True
        print('Reference model successfully assigned')


    def createSequence(self, params=[('dpdp1', 1, 8)]):
        """Create a dipole-dipole sequence.

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
        nelec = self.elec.shape[0]
        def addCustSeq(fname):
            seq = pd.read_csv(fname, header=0)
            if seq.shape[1] != 4:
                raise ValueError('The file should be a CSV file wihtout headers with exactly 4 columns with electrode numbers.')
            else:
                return seq.values
        fdico = {'dpdp1': dpdp1,
              'dpdp2': dpdp2,
              'wenner': wenner,
              'wenner_alpha': wenner_alpha,
              'wenner_beta': wenner_beta,
              'wenner_gamma': wenner_gamma,
              'schlum1': schlum1,
              'schlum2': schlum2,
              'multigrad': multigrad,
              'custSeq': addCustSeq}

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


    def saveSequence(self, fname=''):
        """Save sequence as .csv file.

        Parameters
        ----------
        fname : str, optional
            Path where to save the sequence.
        """
        if self.sequence is not None:
            df = pd.DataFrame(self.sequence, columns=['a','b','m','n'])
            df.to_csv(fname, index=False)
            

    def importElec(self, fname=''):
        """Import electrodes positions.

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
            elec = df.values
        else:
            elec = df
        self.setElec(elec)
                

    def importSequence(self, fname=''):
        """Import sequence for forward modelling.

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


    def saveFilteredData(self, fname, savetyp='Res2DInv (*.dat)'):
        """Save filtered data in formats to be used outside ResIPy (e.g. Res2DInv).

        Parameters
        ----------
        fname : str
            Path where to save the file.
        savetyp : str, optional
            Saving format. To be determined in GUI.
            Default: Res2DInv (*.dat)
        """
        elec = self.elec[['x','y','z']].values
        spacing = elec[1,0] - elec[0,0] # TODO (gb) not sure if this is needed
        for s, i in zip(self.surveys, range(len(self.surveys))):
            df = s.df.query('irecip >=0') # not saving reciprocal data
            # if spacing == None:
            #     spacing = elec[1,0]-elec[0,0] # for batch surveys the spacing can differ and not follow user input
            # else:
            #     spacing = spacing
            # df[['a','b','m','n']] *= spacing
            lookupDict = dict(zip(self.elec['label'], self.elec['x'].values))
            data = df[['a','b','m','n']].replace(lookupDict).values
            df.loc[:,['a','b','m','n']] = data
            if savetyp == 'Res2DInv (*.dat)':
                param = {'num_meas': df.shape[0],
                         'lineTitle': self.param['lineTitle'],
                         'spacing': spacing}
                write2Res2DInv(param, fname, df, elec, self.typ)
            elif savetyp == 'Comma Separated Values (*.csv)':
                write2csv(fname, df, elec, self.typ)
            elif savetyp == 'E4D survey file (*.srv)':
                writeSrv(fname, df, elec)

            fname = fname[:-4] + str(i) + fname[-4:] # to iterate file numbers in case of timelapse survey


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
            r = np.array(self.mesh.attr_cache['res0'])
            phase = np.array(self.mesh.attr_cache['phase0'])
            centroids = np.array(self.mesh.elm_centre).T
            centroids2 = centroids[:,[0,2]] if self.typ[-1] != 't' else centroids
            x = np.c_[centroids2,
                      r,
                      phase, # mrad
                      np.log10(r),
                      np.log10(np.cos(-phase/1000)/np.log10(r)), #log10(real conductivity)
                      np.log10(np.sin(-phase/1000)/np.log10(r))] #log10(imaginary conductivity)
            np.savetxt(os.path.join(fwdDir, 'resistivity.dat'), x)
        else:
            self.mesh.write_attr('res0', os.path.join(fwdDir,'resistivity.dat'))

        # write mesh.dat (no ordering of elements needed in forward mode)
        if (self.typ == 'R2') | (self.typ == 'cR2'):
            self.mesh.write_dat(os.path.join(fwdDir, 'mesh.dat'))
        else:
            self.mesh.write_dat(os.path.join(fwdDir, 'mesh3d.dat'))

        # write the forward .in file
        dump('Writing .in file and mesh.dat... ')
        fparam = self.param.copy()
        fparam['job_type'] = 0
        fparam['num_regions'] = 0
        fparam['res0File'] = 'resistivity.dat' # just starting resistivity
        
        write2in(fparam, fwdDir, typ=self.typ)
        dump('done!\n')

        # write the protocol.dat (that contains the sequence)
        if self.sequence is None:
            dump('Creating sequence... ')
            self.createSequence()
            dump('done!\n')
        dump('Writing protocol.dat... ')
        seq = self.sequence

        # let's check if IP that we have a positive geometric factor
        if self.typ[0] == 'c': # NOTE this doesn't work for borehole
            elecpos = self.elec['x'].values # and works only for 2D
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
        # if it's 3D, we add the line number (all electrode on line 1)
        if self.typ[-2] == '3':
            protocol.insert(1, 'sa', 1)
            protocol.insert(3, 'sb', 1)
            protocol.insert(5, 'sm', 1)
            protocol.insert(7, 'sn', 1)  
            
        outputname = os.path.join(fwdDir, 'protocol.dat')
        with open(outputname, 'w') as f:
            f.write(str(len(protocol)) + '\n')
        with open(outputname, 'a') as f:
            protocol.to_csv(f, sep='\t', header=False, index=False)
        dump('done!\n')

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

        elec = self.elec[['x','y','z']].values
        self.surveys = [] # need to flush it (so no timeLapse forward)
        if self.typ[0] == 'c':
            self.createSurvey(os.path.join(fwdDir, self.typ + '_forward.dat'), ftype='forwardProtocolIP')
        else:
            self.createSurvey(os.path.join(fwdDir, self.typ + '_forward.dat'), ftype='forwardProtocolDC')
        # NOTE the 'ip' columns here is in PHASE not in chargeability
        self.surveys[0].kFactor = 1 # kFactor by default is = 1 now, though wouldn't hurt to have this here!
        self.surveys[0].df['resist'] = addnoise(self.surveys[0].df['resist'].values, self.noise/100)
        self.surveys[0].df['ip'] = addnoiseIP(self.surveys[0].df['ip'].values, self.noiseIP)
        self.surveys[0].computeReciprocal() # to recreate the other columns
        self.setElec(elec) # using R2.createSurvey() overwrite self.elec so we need to set it back

        # recompute doi
        self.computeFineMeshDepth()
        self.zlim[0] = np.min(elec[:,2]) - self.fmd

        if iplot is True:
            self.showPseudo()
        dump('Forward modelling done.')


    def createModelErrorMesh(self, **kwargs):
        """Create an homogeneous mesh to compute modelling error.

        Same arguments as `R2.createMesh()`.
        """
        # backup
        elecZ = self.elec['z'].values.copy()
        mesh = self.mesh.copy() if self.mesh is not None else None
        zlim = self.zlim.copy()
        param = self.param.copy()
        
        # create FLAT homogeneous mesh
        self.elec['z'] = 0
        self.createMesh(**kwargs)
        self.modErrMesh = self.mesh.copy()
        self.modErrMeshNE = self.param['node_elec'].copy()
        self.param['num_regions'] = 0
        
        # restore
        self.mesh = mesh
        self.elec['z'] = elecZ
        self.zlim = zlim
        self.param = param
        

    def estimateError(self, a_wgt=0.01, b_wgt=0.02):
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
            

    def addFlatError(self, percent=2.5):# TODO (gb) why would we want that?
        """Add a flat percentage error to resistivity data (for each survey in
        the class). This action is irreversable.

        resError = res*(percent/100) + resError

        Parameters
        ----------
        percent : float
            Error in percent.
        """
        for s in self.surveys:
            s.addPerError(percent)
            

    def computeModelError(self, rmTree=True):
        """Compute modelling error associated with the mesh.
        This is computed on a flat triangular or tetrahedral mesh.

        Parameters
        ----------
        rmTree : bool
            Remove the working directory used for the error modelling. Default
            is True.
        """
        node_elec = None # we need this as the node_elec with topo and without might be different
        if all(self.elec['z'].values == 0) is False: # so we have topography
            print('New mesh created with flat topo...', end='')
            if 'interp_method' in self.meshParams:
                del self.meshParams['interp_method']
            self.createModelErrorMesh(**self.meshParams)
            node_elec = self.modErrMeshNE
            mesh = self.modErrMesh # create flat mesh
        else:
            mesh = self.mesh # use same mesh

        # create working directory
        fwdDir = os.path.join(self.dirname, 'err')
        if os.path.exists(fwdDir):
            shutil.rmtree(fwdDir)
        os.mkdir(fwdDir)

        # write the resistivity.dat and fparam
        fparam = self.param.copy()
        fparam['job_type'] = 0
        centroids = np.array(mesh.elm_centre).T
        if self.param['mesh_type'] == 6:
            fparam['num_regions'] = 1
            maxElem = centroids.shape[0]
            fparam['regions'] = np.array([[1, maxElem, 100]])
        else:
            if (self.typ == 'R2') | (self.typ == 'cR2'):
                n = 2
                name = 'mesh.dat'
            else:
                n = 3
                name = 'mesh3d.dat'
            resFile = np.zeros((centroids.shape[0],n+1)) # centroix x, y, z, res0
            resFile[:,-1] = 100
            np.savetxt(os.path.join(fwdDir, 'resistivity.dat'), resFile,
                       fmt='%.3f')
            file_path = os.path.join(fwdDir, name)
            mesh.write_dat(file_path)
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
        if (self.typ == 'R3t') | (self.typ == 'cR3t'): # it's a 3D survey
            if len(protocol['a'].values[0].split()) == 1: # we don't have string number
                for c in ['a','b','m','n']: 
                    protocol.loc[:, c] = '1 ' + protocol[c]
            # protocol.insert(1, 'sa', 1)
            # protocol.insert(3, 'sb', 1)
            # protocol.insert(5, 'sm', 1)
            # protocol.insert(7, 'sn', 1)
        outputname = os.path.join(fwdDir, 'protocol.dat')
        with open(outputname, 'w') as f:
            f.write(str(len(protocol)) + '\n')
        with open(outputname, 'a') as f:
            protocol.to_csv(f, sep='\t', header=False, index=False)

        # run the inversion
        self.runR2(fwdDir) # this will copy the R2.exe inside as well

        # get error model
        # if (self.typ == 'R3t') | (self.typ == 'cR3t'):
        #     try:
            # x = np.genfromtxt(os.path.join(fwdDir, self.typ + '_forward.dat'), skip_header=0)
        #     except:#try just reading in the last 2 columns instead
        #         fh = open(os.path.join(fwdDir, self.typ + '.fwd'))
        #         no_meas = len(protocol)
        #         trans_res = [0]*no_meas
        #         app_res = [0]*no_meas
        #         for i in range(no_meas):
        #             line = fh.readline().split()
        #             trans_res[i] = float(line[-2])
        #             app_res[i] = float(line[-1])
        #         x = np.array((trans_res,app_res)).T
        #         fh.close()

        # else:
        x = np.genfromtxt(os.path.join(fwdDir, self.typ + '_forward.dat'), skip_header=1)
        modErr = np.abs(100-x[:,-1])/100
        dferr = pd.DataFrame(seq, columns=['a','b','m','n'])
        dferr['modErr'] = modErr
        for s in self.surveys:
            if 'modErr' in s.df:
                s.df.drop('modErr', axis=1)
            s.df = pd.merge(s.df, dferr, on=['a','b','m','n'], how='inner')

        if rmTree:# eventually delete the directory to spare space
            shutil.rmtree(fwdDir)

        self.fwdErrModel = True # class now has a forward error model.



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
            x = pd.read_csv(os.path.join(self.dirname, fs[index]), delim_whitespace=True).values
            if x.shape[0] > 0:
                triang = tri.Triangulation(x[:,0], x[:,1])
                if self.typ[0] == 'c' and modelDOI is False:
                    z = x[:,4]
                else:
                    z = x[:,3] # modelDOI is always computed with R2 not cR2
#                cax = ax.tricontourf(triang, z, extend='both')
                
#                if self.mesh.surface is not None:
#                    xf, yf = self.mesh.surface[:,0], self.mesh.surface[:,1]
#                    xc, yc, zc = x[:,0], x[:,1], z
#                    zf = interp.nearest(xf, yf, xc, yc, zc) # interpolate before overiding xc and yc
#                    xc = np.r_[xc, xf]
#                    yc = np.r_[yc, yf]
#                    zc = np.r_[zc, zf]
#                    triang = tri.Triangulation(xc, yc) # build grid based on centroids
#                    try:
#                        triang.set_mask(~cropSurface(triang, self.mesh.surface[:,0], self.mesh.surface[:,1]))
#                    except Exception as e:
#                        print('Error in R2.showIter() for contouring: ', e)
#                else:
#                    zc = z.copy()
                cax = ax.tricontourf(triang, z, extend='both')
                self._clipContour(ax, cax.collections, cropMaxDepth)
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
                    
                ax.plot(elec[:,0], yelec, 'ko', markersize=4)
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
                        dfs.append(df)
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
            for df in dfs:
                df = df.astype({'P+':int, 'P-':int, 'C+':int, 'C-':int})
                if (self.typ == 'R3t') | (self.typ == 'cR3t'):
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
        for s, df in zip(self.surveys, dfs):
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
                    ax=ax, log=False, label='Normalised Error', elec=elec)



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
                    ax=ax, log=False, label='Phase misfit')
        

    def showInvError(self, index=0, ax=None):
        """Display inversion error by measurment numbers.
        
        Parameters
        ----------
        index : int, optional
            Index of survey (if time-lapse or batch). Default `index == 0`.
        ax : matplotlib axis
            If provided, the graph will be plotted against this axis.
        """
        errors = self.surveys[index].df['resInvError'].values
        errors = errors[~np.isnan(errors)]
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


    def saveMeshVtk(self, outputname=None):
        """Save mesh as .vtk to be viewed in paraview.

        Parameters
        ----------
        outputname : str, optional
            Output path of the .vtk produced. By default the mesh is saved in
            the working directory `self.dirname` as `mesh.vtk`.
        """
        if outputname is None:
            outputname = os.path.join(self.dirname, 'mesh.vtk')
        self.mesh.write_vtk(outputname)


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
        if self.typ[-1] == '2':
            fname = 'f{:03d}_res.vtk'.format(index+1)
        else:
            fname = 'f{:03d}.vtk'.format(index+1)
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
            attr = list(self.meshResults[index].attr_cache.keys())[0]
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
                elec = self.surveys[i].elec[['x','y','z']].values.copy()
                x = elec[:,0]
                y = elec[:,1]
                self.surveys[i].elec.loc[:,'x'] = y
                self.surveys[i].elec.loc[:,'y'] = x

        for i in range(len(self.surveys)):
            self.surveys[i].elec2distance() # go through each survey and compute electrode
        self.elec = None
        self.setElec(self.surveys[0].elec)

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


    def computeDiff(self):
        """Compute the absolute and the relative difference in resistivity
        between inverted surveys.
        """
        if not self.iTimeLapse:
            raise Exception("Difference calculation only available for time-lapse surveys.")
        if len(self.meshResults) == 0:
            self.getResults()

        # create an index for the values inside of the zone of interest
        # needed as the reference survey is not cropped by default
        inside = np.ones(self.meshResults[0].num_elms, dtype=bool)
        if self.param['num_xz_poly'] > 0:
            meshx = np.array(self.meshResults[0].elm_centre[0])
            meshy = np.array(self.meshResults[0].elm_centre[1])
            meshz = np.array(self.meshResults[0].elm_centre[2])
            # poly = (self.param['xz_poly_table'][:,0],
                    # self.param['xz_poly_table'][:,1])
            path = mpath.Path(self.param['xz_poly_table'])

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
        res_name = res_names[np.in1d(res_names, list(self.meshResults[0].attr_cache.keys()))][0]
        res0 = np.array(self.meshResults[0].attr_cache[res_name])[inside]
        for i in range(1, len(self.meshResults)):
            try:
                res = np.array(self.meshResults[i].attr_cache[res_name])
                self.meshResults[i].add_attribute(res - res0, 'diff(Resistivity)')
                self.meshResults[i].add_attribute((res-res0)/res0*100, 'difference(percent)')
            except Exception as e:
                print('error in computing difference:', e)
                pass
        
        # num_attr = len(self.meshResults[0].attr_cache)
        # num_elm = self.meshResults[0].num_elms
        # baselines = np.zeros((num_attr,num_elm))
        # for i, key in enumerate(self.meshResults[0].attr_cache):
        #     baselines[i,:] = self.meshResults[0].attr_cache[key]
        # change = np.zeros_like(baselines)
        # new_keys = []
        # baseline_keys = []
        # for j, key in enumerate(self.meshResults[0].attr_cache):
        #     new_keys.append('Difference('+key+')')
        #     baseline_keys.append(key)
        # for j, key in enumerate(new_keys):
        #     self.meshResults[0].add_attribute(change[j,:],key)

        # #filter baseline to just the measurements left over after cropping the mesh
        # if crop:
        #     baselines = baselines[:,inside]

        # problem = 0
        # for i in range(1,len(self.meshResults)):
        #     step = self.meshResults[i]
        #     new_keys = []
        #     count = 0
        #     change = np.zeros_like(baselines)
        #     for j, key in enumerate(baseline_keys):
        #         try:
        #             change[count,:] = (np.array(step.attr_cache[key])-baselines[count,:])/baselines[count,:] * 100
        #         except KeyError:
        #             problem+=1
        #         new_keys.append('Difference('+key+')')
        #         count += 1
        #     count = 0
        #     for j, key in enumerate(new_keys):
        #         self.meshResults[i].add_attribute(change[count,:],key)
        #         count += 1
        # if problem>0:
        #     print('Had a problem computing differences for %i attributes'%problem)



    def saveVtks(self, dirname=None):
        """Save vtk files of inversion results to a specified directory.

        Parameters
        ------------
        dirname: str
            Directory in which results will be saved. Default is the working directory.
        """
        if dirname is None:
            dirname = self.dirname
        amtContent = startAnmt
        if len(self.meshResults) == 0:
            self.getResults()
        count=0
        for mesh, s in zip(self.meshResults, self.surveys):
            count+=1
            file_path = os.path.join(dirname, mesh.mesh_title + '.vtk')
            mesh.write_vtk(file_path, title=mesh.mesh_title)
            amtContent += "\tannotations.append('%s')\n"%mesh.mesh_title
        amtContent += endAnmt
        fh = open(os.path.join(dirname,'amt_track.py'),'w')
        fh.write(amtContent)
        fh.close()


    def saveData(self, outputdir):
        """Save all data (_res.dat, .vtk, ...) from the working directory
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

    def estError(self, a_wgt=0.01, b_wgt=0.02): # pragma: no cover
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

    def compDiff(self): # pragma: no cover
        warnings.warn('This function is deprecated, use computeDiff() instead.',
                      DeprecationWarning)
        self.computeDiff()


    def pseudo(self, index=0, vmin=None, vmax=None, ax=None, **kwargs):  # pragma: no cover
        warnings.warn('The function is deprecated, use showPseudo() instead.',
                      DeprecationWarning)
        self.showPseudo(index=index, vmin=vmin, vmax=vmax, ax=ax, **kwargs)
