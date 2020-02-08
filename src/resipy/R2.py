# -*- coding: utf-8 -*-
"""
Main R2 class, wraps the other ResIPy modules (API) in to an object orientated approach
@author: Guillaume, Sina, Jimmy and Paul
"""
ResIPy_version = '2.0.2' # ResIPy version (semantic versionning in use)

#import relevant modules
import os, sys, shutil, platform, warnings, time # python standard libs
from subprocess import PIPE, call, Popen
import subprocess
import numpy as np # import default 3rd party libaries (can be downloaded from conda repositry, incl with winpython)
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import resipy.interpolation as interp # for cropSurface()
import matplotlib.patches as mpatches
import matplotlib.path as mpath

OS = platform.system()
sys.path.append(os.path.relpath('..'))

#import ResIPy resipy packages
from resipy.Survey import Survey
from resipy.r2in import write2in
import resipy.meshTools as mt
from resipy.meshTools import cropSurface
import resipy.isinpolygon as iip
from resipy.template import parallelScript, startAnmt, endAnmt
from resipy.protocol import (dpdp1, dpdp2, wenner_alpha, wenner_beta, wenner,
                          wenner_gamma, schlum1, schlum2, multigrad)
from resipy.SelectPoints import SelectPoints
from resipy.saveData import (write2Res2DInv, write2csv)

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
def readMeshDat(fname):
    """Read mesh.dat or mesh3d.dat and returns elements and nodes.
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
    """Write mesh.dat/mesh3d.dat provided elements and nodes at least.
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


def pseudo(array, resist, spacing, label='', ax=None, contour=False, log=True,
           geom=True, vmin=None, vmax=None):
    print('=======Hey I am usefull you know !')
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
    ypos = np.sqrt(2)/2*np.abs(cmiddle-pmiddle)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    ax.invert_yaxis() # to remove negative sign in y axis
    cax = ax.scatter(xpos, ypos, c=resist, s=70, vmin=vmin, vmax=vmax)#, norm=mpl.colors.LogNorm())
    cbar = fig.colorbar(cax, ax=ax)
    cbar.set_label(label)
    ax.set_title('Pseudo Section')
    ax.set_xlabel('Distance [m]')
    ax.set_ylabel('Pseudo depth [m]')




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
        self.geom_input = {} # dictionnary used to create the mesh
        self.iremote = None # populate on electrode import, True if electrode is remote
        self.iburied = None # populate on electrode import, True if electrode is buried
        
        # attributes needed for independant error model for timelapse/batch inversion
        self.referenceMdl = False # is there a starting reference model already?
        self.fwdErrMdl = False # is there is a forward modelling error already (due to the mesh)?
        self.errTyp = 'global'# type of error model to be used in batch and timelapse surveys


    def setBorehole(self, val=False):
        """Set all surveys in borehole type if `True` is passed.
        """
        self.iBorehole = val
        for s in self.surveys:
            s.iBorehole = val



    def setElec(self, elec, elecList=None):
        """Set electrodes. Automatically identified remote electrode.

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
        
        # identified remote electrode
        remote_flags = [-9999999, -999999, -99999,-9999,-999,
                    9999999, 999999, 99999, 9999, 999] # values asssociated with remote electrodes
        self.iremote = np.in1d(self.elec[:,0], remote_flags)
        if np.sum(self.iremote) > 0:
            print('Detected {:d} remote electrode.'.format(np.sum(self.iremote)))
            for s in self.surveys:
                s.iremote = self.iremote
    
        if len(self.surveys) > 0:
            self.computeDOI()



    def setwd(self, dirname):
        """Set the working directory.

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



    def setTitle(self,linetitle):
        """Set the title of the survey name when inverting data. Input is a string.
        """
        if isinstance(linetitle, str):
            self.param['lineTitle'] = linetitle
        else:
            print("Cannot set Survey title as input is not a string")



    def createSurvey(self, fname='', ftype='Syscal', info={}, spacing=None, parser=None):
        """Read electrodes and quadrupoles data and return 
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
            self.elec = None
            self.setElec(self.surveys[0].elec)


    def createBatchSurvey(self, dirname, ftype='Syscal', info={}, spacing=None,
                          parser=None, isurveys=[], dump=print):
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
        """
        self.iTimeLapse = True
        self.iTimeLapseReciprocal = [] # true if survey has reciprocal
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


        for f in files:
            self.createSurvey(f, ftype=ftype, parser=parser, spacing=spacing)
            haveReciprocal = all(self.surveys[-1].df['irecip'].values == 0)
            self.iTimeLapseReciprocal.append(haveReciprocal)
            dump(f + ' imported')
            print('---------', f, 'imported')
            # all surveys are imported whatever their length, they will be matched
            # later if reg_mode == 2 (difference inversion)
        self.iTimeLapseReciprocal = np.array(self.iTimeLapseReciprocal)
        self.elec = None
        self.setElec(self.surveys[0].elec)
        self.setBorehole(self.iBorehole)

        # create bigSurvey (useful if we want to fit a single error model
        # based on the combined data of all the surveys)
        print('creating bigSurvey')
        self.bigSurvey = Survey(files[0], ftype=ftype, spacing=spacing)
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
            df = df.append(df2, sort=True) # sort to silence the future warning if concatenation axis is not aligned
            c = c + df2.shape[0]
        self.bigSurvey.df = df.copy() # override it
        self.bigSurvey.dfOrigin = df.copy()
        self.bigSurvey.ndata = df.shape[0]

        print("{:d} survey files imported".format(len(self.surveys)))


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
            e[:,1] = i*lineSpacing
            elec.append(e)
            df = s.df.copy()
            df.loc[:,['a','b','m','n']] = df.loc[:,['a','b','m','n']] + i*nelec
            dfs.append(df)
        elec = np.vstack(elec)
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
                s.df['phaseError'] = self.bigSurvey.filterRecipIP(s.df['ip'])
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


    def filterDCA(self, index=-1, dump=print):
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


    def filterElec(self, elec=[]):
        """Filter out specific electrodes given in all surveys.

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


    def filterRecip(self, percent=20,index=-1):
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


    def computeDOI(self):
        """Compute the Depth Of Investigation (DOI) based on electrode
        positions and the larger dipole spacing.
        """
        elec = self.elec.copy()[~self.iremote,:]
        if self.typ == 'R2' or self.typ == 'cR2': # 2D survey:
            if len(self.surveys) > 0:
                array = self.surveys[0].df[['a','b','m','n']].values.copy().astype(int)
                maxDist = np.max(np.abs(elec[array[:,0]-np.min(array[:,0]),0] - elec[array[:,2]-np.min(array[:,2]),0])) # max dipole separation
                self.doi = np.min(elec[:,2])-2/3*maxDist
            else:
                self.doi = np.min(elec[:,2])-2/3*(np.max(elec[:,0]) - np.min(elec[:,0]))

            # set num_xy_poly
            self.param['num_xy_poly'] = 5
            ymax = np.max(elec[:,2])
            ymin = self.doi
            xmin, xmax = np.min(elec[:,0]), np.max(elec[:,0])
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
            xmin, xmax = np.min(elec[:,0]), np.max(elec[:,0])
            ymin, ymax = np.min(elec[:,1]), np.max(elec[:,1])
            zmin, zmax = self.doi, np.max(elec[:,2])
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
                   cl=-1, dump=print, res0=100, show_output=True, doi=None,
                   remote=None, **kwargs):
        """Create a mesh.

        Parameters
        ----------
        typ : str, optional
            Type of mesh. Eithter 'quad' or 'trian' in the case of 2d surveys.
            If no topography, 'quad' mesh will be chosen; 'tetra' is used for 
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
        remote : bool, optional
            Boolean array of electrodes that are remote (ie not real). Should be the same
            length as `R2.elec`.
        kwargs : -
            Keyword arguments to be passed to mesh generation schemes
        """
        self.meshParams = {'typ':typ, 'buried':buried, 'surface':surface,
                           'cl_factor':cl_factor, 'cl':cl, 'dump':dump,
                           'res0': res0, 'show_output':show_output, 'doi':doi}
        if kwargs is not None:
            self.meshParams.update(kwargs)

        if doi is None:# compute depth of investigation if it is not given
            self.computeDOI()
        else:
            self.doi = doi

        if typ == 'default':
            if self.typ == 'R2' or self.typ == 'cR2': # it's a 2D mesh
                typ = 'quad'
                print('Using a quadrilateral mesh.')
            else:
                typ = 'tetra'
                print('Using a tetrahedral mesh.')

        #check if remote electrodes present?
        if remote is None: # automatic detection
            remote = self.iremote
            if np.sum(remote) > 0:
                print('remote electrode detected')
                if typ == 'quad':
                    print('remote electrode is not supported in quadrilateral mesh for now, please use triangular mesh instead.')

        if typ == 'quad':
            elec = self.elec.copy()
            elec_x = self.elec[:,0]
            elec_z = self.elec[:,2]
            #add buried electrodes?
            elec_type = np.repeat('electrode',len(elec_x))
            if buried is None and self.iburied is not None:
                buried = self.iburied
            if (buried is not None
                    and elec.shape[0] == len(buried)
                    and np.sum(buried) != 0):
                elec_type[buried]='buried'
            elec_type = elec_type.tolist()
            surface_x = surface[:,0] if surface is not None else None
            surface_z = surface[:,2] if surface is not None else None
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
        elif typ == 'trian' or typ == 'tetra' or typ=='prism':
            elec = self.elec.copy()
            geom_input = {}
            elec_x = self.elec[:,0]
            elec_y = self.elec[:,1]
            elec_z = self.elec[:,2]
            elec_type = np.repeat('electrode',len(elec_x))
            if buried is None and self.iburied is not None:
                buried = self.iburied
            if (buried is not None
                    and elec.shape[0] == len(buried)
                    and np.sum(buried) != 0):
                elec_type[buried] = 'buried'
            if remote is not None:
                elec_type[remote] = 'remote'

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

            ui_dir = os.getcwd()#current working directory (usually the one the ui is running in)
            os.chdir(self.dirname)#change to working directory so that mesh files written in working directory
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
                mesh = mt.tetra_mesh(elec_x, elec_y, elec_z,elec_type,
                             path=os.path.join(self.apiPath, 'exe'),
                             surface_refinement=surface,
                             cl_factor=cl_factor,
                             cl=cl, dump=dump, show_output=show_output,
                             doi=self.doi-np.max(elec_z), whole_space=whole_space,
                             **kwargs)
            if typ=='prism':
                mesh = mt.prism_mesh(elec_x, elec_y, elec_z,
                                     path=os.path.join(self.apiPath, 'exe'),
                                     cl=cl, dump=dump, show_output=show_output,
                                     **kwargs)
                self.param['num_xy_poly'] = 0
                
            os.chdir(ui_dir)#change back to original directory

            self.param['mesh_type'] = 3
            e_nodes = np.array(mesh.e_nodes) + 1 # +1 because of indexing staring at 0 in python
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

        try:
            self.regions = np.array(self.mesh.cell_attributes)
        except Exception as e:
            print('Error in self.regions=', e)
            self.regions = np.ones(len(self.mesh.elm_centre[0]))
        self.regid = len(np.unique(self.regions)) # no region 0
        self.resist0 = np.ones(len(self.regions))*res0

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

        # checking
        if len(np.unique(e_nodes)) < len(e_nodes):
            raise ValueError('Some electrodes are positionned on the same nodes !')

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
        """Display the mesh.
        """
        if self.mesh is None:
            raise Exception('Mesh undefined')
        else:
            self.mesh.show(ax=ax, color_bar=True, zlim=self.zlim)


    def write2in(self, param={}):
        """Create configuration file for inversion.

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
            #paramFixed = self.mesh.parameteriseMesh()
        self.mesh.write_dat(os.path.join(self.dirname, name),
                            zone = self.mesh.attr_cache['zones'],
                            param = paramFixed)

        # NOTE if fixed elements, they need to be at the end !

        # TODO not sure the sorting fixed element issue if for 3D as well

        if self.mesh.ndims == 2:
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
                    print('+++++++++', s.name)
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
                            self.bigSurvey.fitErrroPwl() # default fit
                        s.df['resError'] = self.bigSurvey.errorModel(s.df)
                    if self.typ[0] == 'c' and np.sum(np.isnan(s.df['phaseError'])) != 0: # there is some nan
                        print('Survey {:s} has no fitted IP error model, default to combined fit.'.format(s.name))
                        if self.bigSurvey.errorModel is None:
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


    def runR2(self, dirname='', dump=print):
        """Run the executable in charge of the inversion.

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



    def runParallel(self, dirname=None, dump=print, iMoveElec=False,
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
                    # TODO get RMS and iteration number here ?
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



    def invert(self, param={}, iplot=False, dump=print, modErr=False,
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
        # clean meshResults list
        self.meshResults = []
            
        # create mesh if not already done
        if 'mesh' not in self.param:
            dump('Create Rectangular mesh...')
            self.createMesh()
            dump('done!\n')
            
        # run Oldenburg and Li DOI estimation
        if modelDOI is True:
#            self.param['num_xy_poly'] = 0 # we need full mesh to match the doiSens
            sensScaled = self.modelDOI(dump=dump)

        # compute modelling error if selected
        if modErr is True and self.fwdErrMdl is False: #check no error model exists
            dump('Computing error model ...')
            self.computeModelError()
            dump('done!\n')
            errTot = True
        elif modErr is True and self.fwdErrMdl:
            # aviod computing error model again if it has already been run.
            errTot = True
        else:
            errTot = False

        # write configuration file
        dump('Writing .in file and protocol.dat...', end='\n')
        self.write2in(param=param) # R2.in
        self.write2protocol(errTot=errTot) # protocol.dat
        #check to make sure the number of electrodes in the protocal matches the
        #the number of electrodes.
        df = self.surveys[0].df
        check = np.array((df['a'],df['b'],df['m'],df['n']))
        if len(self.elec) < np.max(check): # Make sure there are not more electrodes locations in the schedule file than in R2 class
            raise Exception("The number of electrodes given to ResIPy (%i) does not match the number of electrodes parsed in the scheduling file (%i)."%(len(self.elec),np.max(check)))
        dump('done!\n')

        # runs inversion
        if self.iTimeLapse == True and self.referenceMdl==False:
            dump('------------ INVERTING REFERENCE SURVEY ---------------\n')
            refdir = os.path.join(self.dirname, 'ref')
            shutil.move(os.path.join(self.dirname,'res0.dat'),
                        os.path.join(refdir, 'res0.dat'))
            self.write2in(param=param)
            self.runR2(refdir, dump=dump) # this line actually runs R2
            if self.typ=='R3t' or self.typ=='cR3t':
                shutil.copy(os.path.join(refdir, 'f001.dat'),
                            os.path.join(self.dirname, 'Start_res.dat'))
            else:
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
        except:
            print('Could not retrieve files maybe inversion failed')
            return

        if iplot is True:
            self.showResults()


    def modelDOI(self, dump=print):
        """Will rerun the inversion with a background constrain (alpha_s) with
        the normal background and then a background 10 times more resistive.
        From the two different inversion a senstivity limit will be computed.
        """
        # backup normal inversion (0 : original, 1 : normal background, 2: background *10)
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
        if self.param['num_xy_poly'] != 0:
            path = mpath.Path(self.param['xy_poly_table'])
            iselect = path.contains_points(np.c_[self.mesh.elm_centre[0], self.mesh.elm_centre[2]])
            print(np.sum(iselect), len(iselect))
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
        self.mesh.write_attr('res0b', 'res0.dat', self.dirname)
        self.runR2(dump=dump) # re-run inversion
        self.getResults()
        mesh1 = self.meshResults[0]
        cleandir()
        
        # run second background constrained inversion
        dump('===== modelDOI: Running background constrained inversion with initial resistivity * 10 =====\n')
        res2 = res0 * 10
        self.mesh.attr_cache['res0b'] = list(res2)
        self.mesh.write_attr('res0b', 'res0.dat', self.dirname)
        self.runR2(dump=dump) # re-run inversion
        self.getResults()
        mesh2 = self.meshResults[0]
        cleandir()
        os.remove(os.path.join(self.dirname, 'R2.in'))
        
        # sensitivity = difference between final inversion / difference init values
        invValues1 = np.array(mesh1.attr_cache['Resistivity(Ohm-m)'])
        invValues2 = np.array(mesh2.attr_cache['Resistivity(Ohm-m)'])
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
        
    
    
    def _clipContour(self, ax, cont):
        """Clip contours using mesh bound and surface if available.
        
        Parameters
        ----------
        ax : matplotlib.Axes
            Axis.
        cont : matplotlib.collections
            Collection of contours.
        """
        # mask outer region
        xmin = np.min(self.mesh.node_x)
        xmax = np.max(self.mesh.node_x)
        zmin = np.min(self.mesh.node_z)
        zmax = np.max(self.mesh.node_z)
        if self.mesh.surface is not None:
            xsurf, zsurf = self.mesh.surface[:,0], self.mesh.surface[:,1]
            verts = np.c_[np.r_[xmin, xmin, xsurf, xmax, xmax, xmin],
                          np.r_[zmin, zmax, zsurf, zmax, zmin, zmin]]
        else:
            verts = np.c_[np.r_[xmin, xmin, xmax, xmax, xmin],
                          np.r_[zmin, zmax, zmax, zmin, zmin]]                
        # cliping using a patch (https://stackoverflow.com/questions/25688573/matplotlib-set-clip-path-requires-patch-to-be-plotted)
        path = mpath.Path(verts)
        patch = mpatches.PathPatch(path, facecolor='none', edgecolor='none')
        ax.add_patch(patch) # need to add so it knows the transform
        for col in cont.collections:
            col.set_clip_path(patch)
                    
        

    def showResults(self, index=0, ax=None, edge_color='none', attr='',
                    sens=True, color_map='viridis', zlim=None, clabel=None,
                    doi=False, doiSens=False, **kwargs):
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
                mesh = self.meshResults[index]
                if doi is True:
                    if self.doiComputed is True: # DOI based on Oldenburg and Li
                        z = np.array(mesh.attr_cache['doiSens'])
                        levels = [0.2]
                        linestyle = ':'
                    else:
                        raise ValueError('Rerun the inversion with `modelDOI=True` first or use `doiSens`.')
                if doiSens is True: # DOI based on log10(sensitivity)
                    z = np.array(mesh.attr_cache['Sensitivity(log10)'])
                    levels=[np.log10(0.001*(10**np.nanmax(z)))]
                    linestyle = '--'
                
                if doi is True or doiSens is True:
                    # plotting of the sensitivity contour (need to cropSurface as well)
                    xc, yc = np.array(mesh.elm_centre[0]), np.array(mesh.elm_centre[2])
#                    if self.mesh.surface is not None:
#                        zc = z
#                        xf, yf = self.mesh.surface[:,0], self.mesh.surface[:,1]
#                        zf = interp.nearest(xf, yf, xc, yc, zc) # interpolate before overiding xc and yc
#                        xc = np.r_[xc, xf]
#                        yc = np.r_[yc, yf]
#                        zc = np.r_[zc, zf]
#                        triang = tri.Triangulation(xc, yc) # build grid based on centroid
#                        try:
#                            triang.set_mask(~cropSurface(triang, self.mesh.surface[:,0], self.mesh.surface[:,1]))
#                        except Exception as e:
#                            print('Error in cropSurface for contouring: ', e)
#                    else:
#                        triang = tri.Triangulation(xc, yc)
#                        zc = z
                    triang = tri.Triangulation(xc, yc)
                    cont = mesh.ax.tricontour(triang, z, levels=levels, colors='k', linestyles=linestyle)
                    self._clipContour(mesh.ax, cont)
            else: # 3D case
                self.meshResults[index].show(ax=ax,
                            attr=attr, color_map=color_map, clabel=clabel,
                            **kwargs)
        else:
            raise ValueError('len(R2.meshResults) == 0, no inversion results parsed.')



    def getResults(self):
        """Collect inverted results after running the inversion and adding
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
            elec = self.elec.copy()
            iremote = self.iremote
            mesh.elec_x = elec[~iremote,0]
            mesh.elec_y = elec[~iremote,1]
            mesh.elec_z = elec[~iremote,2]
            mesh.surface = self.mesh.surface
            self.meshResults.append(mesh)
        if self.iForward is True:
            initMesh = mt.vtk_import(os.path.join(self.dirname, 'fwd','forward_model.vtk'))
            initMesh.elec_x = self.elec[:,0]
            initMesh.elec_y = self.elec[:,1]
            initMesh.elec_z = self.elec[:,2]
            initMesh.surface = self.mesh.surface
            initMesh.mesh_title = 'Initial Model'
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
                    elec = self.surveys[j].elec.copy()
                    elec = self.elec.copy()
                    iremote = self.iremote
                    mesh.elec_x = elec[~iremote,0]
                    mesh.elec_y = elec[~iremote,1]
                    mesh.elec_z = elec[~iremote,2]
                    mesh.surface = self.mesh.surface
                    self.meshResults.append(mesh)
                    print('done')
                except Exception as e:
                    print('failed', e)
            else:
                pass
                #break

        # compute conductivity in mS/m
        for mesh in self.meshResults:
            if 'Resistivity(Ohm-m)' in mesh.attr_cache.keys():
                mesh.attr_cache['Conductivity(mS/m)'] = 1000/np.array(mesh.attr_cache['Resistivity(Ohm-m)'])

        # compute difference in percent in case of reg_mode == 1
        if (self.iTimeLapse is True) and (self.param['reg_mode'] == 1):
            try:
                self.computeDiff()
            except:
                pass
#            resRef = np.array(self.meshResults[0].attr_cache['Resistivity(Ohm-m)'])
#NOTE: this will not work as the reference array will be bigger than the timesteps if the mesh is cropped
#            for mesh in self.meshResults[1:]:
#                res = np.array(mesh.attr_cache['Resistivity(Ohm-m)'])
#                mesh.attr_cache['difference(percent)'] = (res-resRef)/resRef*100

    
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
        irow = 0
        df = pd.DataFrame(columns=['name', 'dataset', 'iteration', 'resRMS', 'phaseRMS', 'success'])
        for x in lines:
            success = 'N/A'
            line = x.split()
            if len(line) > 1:
                if line[0] == 'Iteration':
                    iiter += 1
                elif (line[0] == 'Measurements') & (line[1] == 'read:'):
                    c = float(line[2])
                    d = float(line[5])
                elif line[0] == 'Final':
                    resRMS = float(line[3])
                    df.loc[irow, :] = [name, idataset, iiter, resRMS, phaseRMS, success]
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
        ax.set_xticks([],[])
        

    def showSection(self, fname='', ax=None, ilog10=True, isen=False, figsize=(8,3)):
        """Show inverted section based on the `_res.dat``file instead of the
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


    def addRegion(self, xy, res0=100, phase0=1, blocky=False, fixed=False,
                  ax=None, iplot=False):
        """Add region according to a polyline defined by `xy` and assign it
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
        iplot : bool, optional
            If `True` , the updated mesh with the region will be plotted.
        """
#        selector = SelectPoints(ax, np.array(self.mesh.elm_centre).T[:,[0,2]],
#                                typ='poly', iplot=iplot) # LIMITED FOR 2D case
#        selector.setVertices(xy)
#        selector.getPointsInside()
#        idx = selector.iselect

        centroids = np.array(self.mesh.elm_centre).T[:,[0,2]]
        path = mpath.Path(np.array(xy))
        idx = path.contains_points(centroids)

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

        if iplot is True:
            self.showMesh()


    def resetRegions(self):
        """Just reset all regions already draw. Shouldn't be needed as
        the `self.runR2()` automatically use a homogenous model as starting
        for inversion. The only purpose of this is to use an inhomogeous
        starting model to invert data from forward modelling.
        """
        self.regid = 1
        self.regions.fill(1)
        self.mesh.attr_cache['res0'] = np.ones(len(self.regions))*100 # set back as default


    def createModel(self, ax=None, dump=print, typ='poly', addAction=None):
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
        if self.mesh is None:
            print('will create a mesh before')
            self.createMesh()
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        def callback(idx):
            self.regid = self.regid + 1
            print('nb elements selected:', np.sum(idx), 'in region', self.regid)
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



    def designModel(self, ax=None, dump=print, typ='poly', addAction=None):
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

        Returns
        -------
        fig : matplotlib.figure
            If `ax` is `None`, will return a figure.
        """
        if self.doi is None:
            self.computeDOI()

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        self.geom_input = {}
        ax.plot(self.elec[:,0], self.elec[:,2], 'ko', label='electrode')
        ax.set_ylim([self.doi, np.max(self.elec[:,1])])
        ax.set_xlim(np.min(self.elec[:,0]), np.max(self.elec[:,0]))
        def callback():
            vert = np.array(self.selector.vertices)
            self.geom_input['polygon' + str(len(self.geom_input)+1)] = [vert[:-1,0].tolist(), vert[:-1,1].tolist()]
            ax.plot(vert[:,0], vert[:,1], '.-')
            if addAction is not None:
                addAction()
        # we need to assign a selector to self otherwise it's not used
        self.selector = SelectPoints(ax, typ=typ, callback=callback)
#        surveyLength = np.max(self.elec[:,0]) - np.min(self.elec[:,0])
        self.selector.xmin = np.min(self.elec[:,0])# - 10 * surveyLength
        self.selector.xmax = np.max(self.elec[:,0])# + 10 * surveyLength
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



    def setRefModel(self, res0):
        """Set the reference model according to a previous inversion, avoids
        the need to invert reference model again for timelapse workflows.

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
        self.mesh.write_attr('res0',file_name='Start_res.dat',file_path=self.dirname)
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
        nelec = len(self.elec)
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
            if p[0] is 'custSeq':
                try:
                    qs.append(addCustSeq(p[1]))
                except Exception as e:
                    print('error when importing custom sequence:', e)
            else:
                pok = [int(p[i]) for i in np.arange(1, len(p))] # make sure all are int
                qs.append(fdico[p[0]](nelec, *pok).values.astype(int))
        self.sequence = np.vstack(qs)


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
        if df.shape[1] > 4:
            raise ValueError('The file should have no more than 4 columns')
        else:
            if header is not None:
                elec = df[['x','y','z']].values
            else:
                elec = df.values
            self.setElec(elec)
            if 'buried' in df.columns:
                self.iburied = df['buried'].values.astype(bool)
            else:
                self.iburied = None
                
                

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


    def saveFilteredData(self, fname, elec, savetyp='Res2DInv (*.dat)', spacing=None):
        """Save filtered data in formats to be used outside ResIPy (e.g. Res2DInv).

        Parameters
        ----------
        fname : str
            Path where to save the file.
        elec : Array
            Array containing topography information.
        savetyp : str, optional
            Saving format. To be determined in GUI.
            Default: Res2DInv (*.dat)
        """
        for s, i in zip(self.surveys, range(len(self.surveys))):
            df = s.df.query('irecip >=0') # not saving reciprocal data
            if spacing == None:
                spacing = elec[1,0]-elec[0,0] # for batch surveys the spacing can differ and not follow user input
            else:
                spacing = spacing
            df[['a','b','m','n']] *= spacing
            if savetyp == 'Res2DInv (*.dat)':
                param = {'num_meas':len(df),
                         'lineTitle':self.param['lineTitle'],
                         'spacing':spacing}
                write2Res2DInv(param, fname, df, elec, self.typ)
            elif savetyp == 'Comma Separated Values (*.csv)':
                write2csv(fname, df, elec, self.typ)

            fname = fname[:-4]+str(i)+fname[-4:] # to iterate file numbers in case of timelapse survey


    def forward(self, noise=0.0, noiseIP=0.0, iplot=False, dump=print):
        """Operates forward modelling.

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
        dump('Writing .in file...', end='\n')
        fparam = self.param.copy()
        fparam['job_type'] = 0
        fparam['num_regions'] = 0
        fparam['res0File'] = 'resistivity.dat' # just starting resistivity

        write2in(fparam, fwdDir, typ=self.typ)
        dump('done!\n')

        # write the protocol.dat (that contains the sequence)
        if self.sequence is None:
            dump('Creating sequence ...')
            self.createSequence()
            dump('done!\n')
        dump('Writing protocol.dat ...')
        seq = self.sequence

        # let's check if IP that we have a positive geometric factor
        if self.typ[0] == 'c': # NOTE this doesn't work for borehole
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
        elif self.typ[-2] == '3':
            self.createSurvey(os.path.join(fwdDir, self.typ + '.fwd'), ftype='forwardProtocolDC')
        else:
            self.createSurvey(os.path.join(fwdDir, self.typ + '_forward.dat'), ftype='forwardProtocolDC')
        # NOTE the 'ip' columns here is in PHASE not in chargeability
        self.surveys[0].kFactor = 1 # kFactor by default is = 1 now, though wouldn't hurt to have this here!
        self.surveys[0].df['resist'] = addnoise(self.surveys[0].df['resist'].values, self.noise)
        self.surveys[0].df['ip'] = addnoiseIP(self.surveys[0].df['ip'].values, self.noiseIP)
        self.setElec(elec) # using R2.createSurvey() overwrite self.elec so we need to set it back

        # recompute doi
        self.computeDOI()
        self.zlim[0] = self.doi

        if iplot is True:
            self.showPseudo()
        dump('Forward modelling done.')



    def createModelErrorMesh(self, typ='default', buried=None, surface=None, cl_factor=2,
                   cl=-1, dump=print, res0=100, show_output=True, doi=None, **kwargs):
        """Create an homogeneous mesh to compute modelling error.

        Same arguments as `R2.createMesh()`.
        """
        fix_me = False
        try:
            old_attr_cache = self.mesh.attr_cache.copy()
            fix_me = True
        except AttributeError:
            pass
        if typ == 'quad':
            elec = self.elec.copy()
            elec_x = self.elec[:,0]
            elec_z = np.zeros(len(elec_x))# THIS IS KEY HERE we need a flat surface
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

            if 'regions' in self.param: # allow to create a new mesh then rerun inversion
                del self.param['regions']
            if 'num_regions' in self.param:
                del self.param['num_regions']
        elif typ == 'trian' or typ == 'tetra':
            elec = self.elec.copy()
            geom_input = {}
            elec_x = self.elec[:,0]
            elec_y = self.elec[:,1]
            elec_z = np.zeros(len(elec_y)) # THIS IS KEY HERE we need a flat surface
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
                             interp_method = None, # this aviods doing a lengthy interpolation
                             cl_factor=cl_factor,
                             cl=cl, dump=dump, show_output=show_output,
                             doi=self.doi-np.max(elec_z), whole_space=whole_space,
                             **kwargs)
            os.chdir(ui_dir)#change back to original directory
            e_nodes = mesh.e_nodes + 1 # +1 because of indexing staring at 0 in python
            self.modErrMeshNE = np.c_[1+np.arange(len(e_nodes)), e_nodes].astype(int)

        self.modErrMesh = mesh

        self.param['num_regions'] = 0

        numel = self.modErrMesh.num_elms
        self.modErrMesh.add_attribute(np.ones(numel)*res0, 'res0') # default starting resisivity [Ohm.m]
        self.modErrMesh.add_attribute(np.ones(numel)*0, 'phase0') # default starting phase [mrad]
        self.modErrMesh.add_attribute(np.ones(numel, dtype=int), 'zones')
        self.modErrMesh.add_attribute(np.zeros(numel, dtype=bool), 'fixed')
        self.modErrMesh.add_attribute(np.zeros(numel, dtype=float), 'iter')
        if fix_me:
            self.mesh.attr_cache = old_attr_cache


    def estimateError(self, a_wgt=0.01, b_wgt=0.02):
        """Estimate reciprocal error data for data with no recipricols for each
        survey, using the same routine present in R2. This allows for the additional inclusion
        of modelling errors. This action is irreversable.

        Parameters
        ------------
        a_wgt: float, optional
            a_wgt documented in the R2 documentation
        b_wgt: float, optional
            b_wgt documented in the R2 documentation
        """
        for s in self.surveys:
            s.estimateError(a_wgt=a_wgt,b_wgt=b_wgt)

    def addFlatError(self,pnct=2.5):
        """Add a flat percentage error to resistivity data (for each survey in
        the class). This action is irreversable.

        resError = res*(pnct/100) + resError

        Parameters
        --------------
        pnct: float
            Error in percent
        """
        for s in self.surveys:
            s.addPerError(pnct)

    def computeModelError(self,rm_tree=True):
        """Compute modelling error assocaited with the mesh.
        This is computed on a flat tetrahedral mesh.

        Parameters
        ------------
        rm_tree: bool
            Remove the working directory used for the error modelling. Default
            is True.
        """
        try:#bug fix for overwriting attr_cache in mesh object
            attr_cache = self.mesh.attr_cache.copy()
            #for some reason the mesh.attr_cache is dynamically linked to the modelling mesh, and i cant figure out why
        except AttributeError:
            print("No mesh already in place")

        node_elec = None # we need this as the node_elec with topo and without might be different
        if all(self.elec[:,2] == 0) is False: # so we have topography
            print('A new mesh will be created as the surface is not flat.')
            try:
                _ = self.meshParams
            except AttributeError:
                self.meshParams = {}
            if 'interp_method' in self.meshParams:
                del self.meshParams['interp_method']
            self.createModelErrorMesh(**self.meshParams)
            node_elec = self.modErrMeshNE
            mesh = self.modErrMesh
            fix_me = True
        else:
            mesh = self.mesh
            fix_me = False

        fwdDir = os.path.join(self.dirname, 'err')
        if os.path.exists(fwdDir):
            shutil.rmtree(fwdDir)
        os.mkdir(fwdDir)

        # write the resistivity.dat and fparam
        fparam = self.param.copy()
        fparam['job_type'] = 0
        centroids = np.array(mesh.elm_centre).T
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
            name = 'mesh.dat'
            if self.typ == 'R3t' or self.typ == 'cR3t':
                name = 'mesh3d.dat'
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
            seq = np.vstack(seq)
            seq = np.unique(seq, axis=0)
        protocol = pd.DataFrame(np.c_[1+np.arange(seq.shape[0]),seq])
        if self.typ == 'R3t' or self.typ == 'cR3t': # it's a 3D survey
            protocol.insert(1, 'sa', 1)
            protocol.insert(3, 'sb', 1)
            protocol.insert(5, 'sm', 1)
            protocol.insert(7, 'sn', 1)
        outputname = os.path.join(fwdDir, 'protocol.dat')
        with open(outputname, 'w') as f:
            f.write(str(len(protocol)) + '\n')
        with open(outputname, 'a') as f:
            protocol.to_csv(f, sep='\t', header=False, index=False)

        # run the inversion
        self.runR2(fwdDir) # this will copy the R2.exe inside as well

        # get error model
        if self.typ[-2] == '3':
            try:
                x = np.genfromtxt(os.path.join(fwdDir, self.typ + '.fwd'), skip_header=0)
            except:#try just reading in the last 2 columns instead
                fh = open(os.path.join(fwdDir, self.typ + '.fwd'))
                no_meas = len(protocol)
                trans_res = [0]*no_meas
                app_res = [0]*no_meas
                for i in range(no_meas):
                    line = fh.readline().split()
                    trans_res[i] = float(line[-2])
                    app_res[i] = float(line[-1])
                x = np.array((trans_res,app_res)).T
                fh.close()

        else:
            x = np.genfromtxt(os.path.join(fwdDir, self.typ + '_forward.dat'), skip_header=1)
        modErr = np.abs(100-x[:,-1])/100
        dferr = pd.DataFrame(np.c_[seq, modErr], columns=['a','b','m','n','modErr'])
        for s in self.surveys:
            if 'modErr' in s.df:
                s.df.drop('modErr', axis=1)
            s.df = pd.merge(s.df, dferr, on=['a','b','m','n'], how='inner')

        if rm_tree:# eventually delete the directory to spare space
            shutil.rmtree(fwdDir)

        if fix_me: #apply fix to sort attr_cache inside mesh object
            self.mesh.attr_cache = attr_cache.copy()

        self.fwdErrMdl = True # class now has a forward error model.


    def showIter(self, index=-2, ax=None, modelDOI=False):
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
#            if self.param['mesh_type'] == 10:
#                self.showSection(os.path.join(self.dirname, fs[index]), ax=ax)
#                # TODO change that to full meshTools?
#
#            else:
#            x = np.genfromtxt(os.path.join(self.dirname, fs[index])) # too sensitive to empty columns of cR2 output
            x = pd.read_csv(os.path.join(self.dirname, fs[index]), delim_whitespace=True).values
            if x.shape[0] > 0:
                triang = tri.Triangulation(x[:,0],x[:,1])
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
                self._clipContour(ax, cax)
                fig.colorbar(cax, ax=ax, label=r'$\rho$ [$\Omega$.m]')
                ax.plot(self.elec[:,0], self.elec[:,2], 'ko', markersize=4)
                ax.set_aspect('equal')
                ax.set_xlabel('Distance [m]')
                ax.set_ylabel('Elevation [m]')
                ax.set_xlim([np.min(self.elec[~self.iremote,:][:,0]), np.max(self.elec[~self.iremote,:][:,0])])
                ax.set_ylim(self.zlim)
                if iplot is True:
                    fig.show()


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
#                kwargs2 = kwargs.copy()
#                fig, ax = plt.subplots()
#                if 'ylim' not in kwargs2:
#                    ylim = [self.doi, np.max(self.elec[:,2])]
#                    kwargs2 = dict(kwargs2, ylim=ylim)
#                if 'color_map' not in kwargs2:
#                    kwargs2 = dict(kwargs2, color_map='viridis')
#                if 'attr' in kwargs2:
#                    if kwargs2['attr'] not in list(self.meshResults[i].attr_cache.keys()):
#                        kwargs2['attr'] = 'Resistivity(log10)'
#                self.meshResults[i].show(ax=ax, **kwargs2)
            fig, ax = plt.subplots()
            self.showResults(index=i, ax=ax, **kwargs)
            fname = self.meshResults[i].mesh_title
            fig.savefig(os.path.join(outputdir, fname + '.png'))


#    def getInvError(self, index=0):
#        """Collect inversion error from _err.dat or .err file after inversion.
#
#        Parameters
#        ----------
#        index : int, optional
#            Index of the survey (if Timelapse or batch). Default is 0.
#            
#        Returns
#        -------
#        array : numpy.array
#            Contains the quadrupoles.
#        errors : numpy.array
#            Vector of normalized error.
#        """
#        if self.typ == 'cR2' or self.typ == 'R2':
#            df = pd.read_csv(os.path.join(self.dirname, 'f{:03.0f}_err.dat'.format(index+1)), delim_whitespace=True)
#            array = np.array([df['C+'],df['C-'],df['P+'],df['P-']],dtype=int).T
#            errors = np.array(df['Normalised_Error'])
#        elif self.typ == 'R3t' or self.typ == 'cR3t':
#            err = np.genfromtxt(os.path.join(self.dirname, 'f{:03.0f}.err'.format(index+1)), skip_header=1)
#            array = err[:,[-3,-1,-7,-5]].astype(int)
#            errors = err[:,0]
#
#        return array, errors
    
    
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
            elif self.typ == 'R3t' or self.typ == 'cR3t':
                dfs = []
                if self.iTimeLapse:
                    fname = os.path.join(self.dirname, 'ref/f001.err')
                    if os.path.exists(fname):
                        err = np.genfromtxt(fname, skip_header=1)
                        df = pd.DataFrame(err[:,[-3, -1, -7, -5, 0]],
                                          columns=['P+','P-','C+','C-', 'Normalised_Error'])
                        dfs.append(df)
                for i in range(len(self.surveys)-a):
                    fname = os.path.join(self.dirname, 'f{:03.0f}.err'.format(i+1))
                    if os.path.exists(fname):
                        err = np.genfromtxt(fname, skip_header=1)
                        df = pd.DataFrame(err[:,[-3, -1, -7, -5, 0]],
                                          columns=['P+','P-','C+','C-', 'Normalised_Error'])
                        dfs.append(df)
            #TODO not implemented for cR3t and phase misfit
        except:
            return # this code is error prone (mainly to empty dataframe error)
        
        # merge the columns to each survey dataframe
        if  np.sum([df.shape[0] > 0 for df in dfs]) != len(self.surveys):
            print('error in reading error files (do not exists or empty')
            return # this check the number of dfs AND the fact that they are not empty
        for s, df in zip(self.surveys, dfs):
            if self.typ == 'cR2': #TODO figure out why Andy's code produce different f001_err.dat files
                df = df.rename(columns=dict(zip(['C+','C-','P+','P-', 'Normalised_Error'], ['a','b','m','n', 'resInvError']))) #there is something wrong here. R2 and cR2 produce different f001_err.dat! 'P+','P-','C+','C-' are different!!
            elif self.typ == 'R2':
                df = df.rename(columns=dict(zip(['P+','P-','C+','C-', 'Normalised_Error'], ['a','b','m','n', 'resInvError'])))
            else: # for 3D ones
                df = df.rename(columns=dict(zip(['C+','C-','P+','P-', 'Normalised_Error'], ['a','b','m','n', 'resInvError'])))
            cols = ['a','b','m','n','resInvError']
            if self.typ == 'cR2':
                df['phaseInvMisfit'] = np.abs(df['Observed_Phase'] - df['Calculated_Phase'])
                cols += ['phaseInvMisfit']
            if 'resInvError' in s.df.columns:
                s.df = s.df.drop('resInvError', axis=1)
            if 'phaseInvMisfit' in s.df.columns:
                s.df = s.df.drop('phaseInvMisfit', axis=1)
            s.df = pd.merge(s.df, df[cols], on=['a','b','m','n'], how='left')
        # TODO assign the errors to normal and reciprocal ? in case we use recipMean only ? 
        # This error has nothing to do with reciprocity!

                    

    def showPseudoInvError(self, index=0, ax=None, vmin=None, vmax=None):
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
        """
        self.surveys[index].filterManual(attr='resInvError', vmin=vmin, vmax=vmax,
                    ax=ax, geom=False, log=False)


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
                    ax=ax, geom=False, log=False)
        

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


    def _toParaview(self, fname,  paraview_loc=None):
        """Open file in paraview.

        Parameters
        ----------
        fname : str
            Path of the .vtk file to be opened.
        paraview_loc: str, optional
            **Windows ONLY** maps to the excuatable paraview.exe. The program
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


    def showMeshInParaview(self, paraview_loc=None):
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


    def showInParaview(self, index=0, paraview_loc=None):
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


    def showSlice(self, index=0, ax=None, attr=None, axis='z'):
        """Show slice of 3D mesh interactively.
        """
        if attr is None:
            attr = list(self.meshResults[index].attr_cache.keys())[0]
        self.meshResults[index].showSlice(
                attr=attr, axis=axis)

    ## Sorting electrode numbers ##
    def shuntIndexes(self):
        """Shunt electrode indexes to start at 1.
        """
        debug=True
        if len(self.surveys)>1:
            debug=False
        for i in range(len(self.surveys)):
            self.surveys[i].shuntIndexes(debug=debug)

    def normElecIdx(self):
        """Normalise electrode indexes to start at 1 in consective and ascending order.
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
                self.surveys[i].elec[:,0] = y
                self.surveys[i].elec[:,1] = x

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



    def computeCond(self):
        """Compute conductivities from resistivities for the ERT mesh
        """
        if self.typ=='R3t' or self.typ=='cR3t':
            res_name = 'Resistivity'
        else:
            res_name = 'Resistivity(Ohm-m)'
        for i in range(len(self.meshResults)):
            self.meshResults[i].computeReciprocal(res_name,'Conductivity(S/m)')


    def computeDiff(self):
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

    def pseudoIP(self, index=0, vmin=None, vmax=None, ax=None, **kwargs):
        warnings.warn('The function is deprecated, use showPseudoIP() instead.',
                      DeprecationWarning)
        self.showPseudoIP(index=index, vmin=vmin, vmax=vmax, ax=ax, **kwargs)

    def plotError(self, index=0, ax=None):
        warnings.warn('This function is deprecated, use showError() instead.',
                      DeprecationWarning)
        self.showError(index=index, ax=ax)

    def errorDist(self, index=0, ax=None):
        warnings.warn('This function is deprecated, use showErrorDist() instead.',
                      DeprecationWarning)
        self.showErrorDist(index=index, ax=ax)

    def removeDummy(self, index=-1):
        warnings.warn('This function is deprecated, use filterDummy() instead.',
                      DeprecationWarning)
        self.filterDummy(index=index)

    def linfit(self, index=-1, ax=None):
        warnings.warn('This function is deprecated, use fitErrorLin() instead.',
                      DeprecationWarning)
        self.fitErrorLin(index=index, ax=ax)


    def pwlfit(self, index=-1, ax=None):
        warnings.warn('This function is deprecated, use fitErrorPwl() instead.',
                      DeprecationWarning)
        self.fitErrorPwl(index=index, ax=ax)

    def lmefit(self, index=-1, ax=None, rpath=None, iplot=True):
        warnings.warn('This function is deprecated, use fitErrorLME() instead.',
                      DeprecationWarning)
        self.fitErrorLME(index=index, ax=ax, rpath=rpath, iplot=iplot)

    def phaseplotError(self, index=0, ax=None):
        warnings.warn('This function is deprecated, use showErrorIP() instead.',
                      DeprecationWarning)
        self.showErrorIP(index=index, ax=ax)

    def plotIPFit(self, index=-1, ax=None):
        warnings.warn('This function is deprecated, use fitErrorPwlIP() instead.',
                      DeprecationWarning)
        self.fitErrorPwlIP(index=index, ax=ax)

    def plotIPFitParabola(self, index=-1, ax=None):
        warnings.warn('This function is deprecated, use fitErrorParabolaIP() instead.',
                      DeprecationWarning)
        self.fitErrorParabolaIP(index=index, ax=ax)

    def heatmap(self, index=0, ax=None):
        warnings.warn('This function is deprecated, use showHeatmap() instead.',
                      DeprecationWarning)
        self.showHeatmap(index=index, ax=ax)

    def removenested(self, index=-1):
        warnings.warn('This function is deprecated, use filterNested() instead.',
                      DeprecationWarning)
        self.filterNested(index=index)

    def dca(self, index=-1, dump=print):
        warnings.warn('This function is deprecated, use filterDCA() instead.',
                      DeprecationWarning)
        self.filterDCA(index=index, dump=dump)

    def removerecip(self, index=0):
        warnings.warn('This function is deprecated, use filterRecip() instead.',
                      DeprecationWarning)
        self.filterRecip(index=index)

    def iprangefilt(self, phimin, phimax, index=-1):
        warnings.warn('This function is deprecated, use filterRangeIP() instead.',
                      DeprecationWarning)
        self.filterRangeIP(phimin, phimax, index=index)
        
    def removeUnpaired(self, index=-1):
        warnings.warn('This function is deprecated, use filterUnpaired() instead.',
                      DeprecationWarning)
        n = self.filterUnpaired(index=index)
        return n
    
    def removeneg(self):
        warnings.warn('This function is deprecated, use filterNegative() instead.',
                      DeprecationWarning)
        self.filterNegative()

    def assignRes0(self, regionValues={}, zoneValues={}, fixedValues={}, ipValues={}):
        warnings.warn('This function is deprecated, use setStartingRes() instead.',
                      DeprecationWarning)
        self.setStartingRes(regionValues=regionValues, zoneValues=zoneValues, fixedValues=fixedValues, ipValues=ipValues)


    def assignRefModel(self, res0):
        warnings.warn('This function is deprecated, use setRefModel() instead.',
                      DeprecationWarning)
        self.setRefModel(res0=res0)


    def createModellingMesh(self, typ='default', buried=None, surface=None, cl_factor=2,
                   cl=-1, dump=print, res0=100, show_output=True, doi=None, **kwargs):
        warnings.warn('This function is deprecated, use createModelErrorMesh() instead.',
                      DeprecationWarning)
        self.createModelErrorMesh(typ=typ, buried=buried, surface=surface, cl_factor=cl_factor,
                                  cl=cl, dump=dump, res0=res0, show_output=show_output, doi=doi, **kwargs)

    def estError(self, a_wgt=0.01, b_wgt=0.02):
        warnings.warn('This function is deprecated, use estimateError() instead.',
                      DeprecationWarning)
        self.estimateError(a_wgt=a_wgt, b_wgt=b_wgt)

    def pseudoError(self, ax=None, vmin=None, vmax=None):
        warnings.warn('This function is deprecated, use showPseudInvError() instead.',
                      DeprecationWarning)
        self.showPseudoInvError(ax=ax, vmin=vmin, vmax=vmax)


    def pseudoErrorIP(self, ax=None, vmin=None, vmax=None):
        warnings.warn('This function is deprecated, use showErrorIP() instead.',
                      DeprecationWarning)
        self.showPseudoErrorIP(ax=ax, vmin=vmin, vmax=vmax)


    def showInversionErrors(self, ax=None):
        warnings.warn('This function is deprecated, use showInvError() instead.',
                      DeprecationWarning)
        self.showInvError(ax=ax)


    def compCond(self):
        warnings.warn('This function is deprecated, use computeCond() instead.',
                      DeprecationWarning)
        self.computeCond()

    def compDiff(self):
        warnings.warn('This function is deprecated, use computeDiff() instead.',
                      DeprecationWarning)
        self.computeDiff()


    def pseudo(self, index=0, vmin=None, vmax=None, ax=None, **kwargs):
        warnings.warn('The function is deprecated, use showPseudo() instead.',
                      DeprecationWarning)
        self.showPseudo(index=index, vmin=vmin, vmax=vmax, ax=ax, **kwargs)
