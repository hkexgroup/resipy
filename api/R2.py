# -*- coding: utf-8 -*-
"""
Created on Wed May 30 16:48:54 2018

@author: jkl
"""

import numpy as np
import pandas as pd
import os
import sys
import shutil
import platform
import matplotlib.pyplot as plt
from subprocess import PIPE, call, Popen
OS = platform.system()

sys.path.append(os.path.relpath('..'))

from api.Survey import Survey
from api.r2in import write2in
import api.meshTools as mt
from api.meshTools import Mesh_obj, tri_mesh, dat_import
from api.sequenceHelper import ddskip
from api.SelectPoints import SelectPoints
from api.post_processing import import_R2_err_dat, disp_R2_errors

apiPath = os.path.abspath(os.path.join(os.path.abspath(__file__), '../'))
print('API path = ', apiPath)

class R2(object): # R2 master class instanciated by the GUI
    """ Master class to handle all processing around the inversion codes.
    """
    def __init__(self, dirname='', typ='R2'):
        """ Create an R2 object.
        
        Parameters
        ----------
        dirnaname : str, optional
            Path of the working directory. Can also be set using `R2.setwd()`.
        """
        self.apiPath = os.path.dirname(os.path.abspath(__file__)) # directory of the code
        if dirname == '':
            dirname = os.path.join(self.apiPath, 'invdir')
        print('Working directory is:', dirname)
        self.setwd(dirname) # working directory (for the datas)
        self.surveys = [] # list of survey object
        self.surveysInfo = [] # info about surveys (date)
        self.mesh = None # mesh object (one per R2 instance)
        self.param = {} # dict configuration variables for inversion
        self.configFile = ''
        self.typ = typ # or cR2 or R3, cR3
        self.errTyp = 'none' # type of error to add for DC
        self.errTypIP = 'none' # type of error to add for IP phase
        self.iBorehole = False # to tell the software to not plot pseudoSection
        self.iTimeLapse = False # to enable timelapse inversion
        self.iBatch = False # to enable batch inversion
        self.meshResults = [] # contains vtk mesh object of inverted section
        self.sequence = None # quadrupoles sequence if forward model
        self.resist0 = None # initial resistivity
        self.iForward = False # if True, it will use the output of the forward
        # to run an inversion (and so need to reset the regions before this)
        
    def setwd(self, dirname):
        """ Set the working directory.
        
        Parameters
        ----------
        dirname : str
            Path of the working directory.
        """
        # get rid of some stuff
        files = os.listdir(dirname)
        if 'ref' in files: # only for timelapse survey
            shutil.rmtree(os.path.join(dirname, 'ref'))
        if 'err' in files: # only for error modelling
            shutil.rmtree(os.path.join(dirname, 'err'))
        if 'R2.exe' in files:
            os.remove(os.path.join(dirname, 'R2.exe'))
        if 'cR2.exe' in files:
            os.remove(os.path.join(dirname, 'cR2.exe'))
        if 'mesh.dat' in files:
            os.remove(os.path.join(dirname, 'mesh.dat'))
        if 'Start_res.dat' in files:
            os.remove(os.path.join(dirname, 'Start_res.dat'))
        def isNumber(x):
            try:
                float(x)
                return True
            except:
                return False
        for f in files:
            if (f[0] == 'f') & (isNumber(f[1:3]) is True):
                os.remove(os.path.join(dirname, f))
        self.dirname = dirname
    
    
    def setBorehole(self, val=False):
        """ Set all surveys in borehole type if `True` is passed.
        """
        self.iBorehole = val
        for s in self.surveys:
            s.iBorehole = val
    
    
    def createSurvey(self, fname='', ftype='Syscal', info={}, spacing=None, parser=None):
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
        
        # define electrode position according to first survey
        if len(self.surveys) == 1:
            self.elec = self.surveys[0].elec
            
            # attribute method of Survey object to R2
            self.pseudoIP = self.surveys[0].pseudoIP
            self.pseudo = self.surveys[0].pseudo
            self.plotError = self.surveys[0].plotError
            self.linfit = self.surveys[0].linfit
            self.lmefit = self.surveys[0].lmefit
            self.pwlfit = self.surveys[0].pwlfit
            self.phaseplotError = self.surveys[0].phaseplotError
            self.plotIPFit = self.surveys[0].plotIPFit
            self.heatmap = self.surveys[0].heatmap
            self.iprangefilt = self.surveys[0].iprangefilt
            self.removerecip = self.surveys[0].removerecip
            self.removenested = self.surveys[0].removenested
        
    def createBatchSurvey(self, dirname, ftype='Syscal', info={}, spacing=None, parser=None, isurveys=[]):
        """ Read multiples files from a folders (sorted by alphabetical order).
        
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
        isurveys : list
            List of surveys index that will be used for error modelling and so
            reciprocal measurements. By default all surveys are used.
        """  
        self.createTimeLapseSurvey(dirname=dirname, ftype=ftype, info=info, spacing=spacing, isurveys=isurveys, parser=parser)
        self.iTimeLapse = False
        self.iBatch = True

    def createTimeLapseSurvey(self, dirname, ftype='Syscal', info={}, spacing=None, parser=None, isurveys=[]):
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
        isurveys : list
            List of surveys index that will be used for error modelling and so
            reciprocal measurements. By default all surveys are used.
        """    
        self.iTimeLapse = True
        self.iTimeLapseReciprocal = [] # true if survey has reciprocal
        files = np.sort(os.listdir(dirname))
        for f in files:
            self.createSurvey(os.path.join(dirname, f), ftype=ftype, parser=parser)
            haveReciprocal = all(self.surveys[-1].df['irecip'].values == 0)
            self.iTimeLapseReciprocal.append(haveReciprocal)
            print('---------', f, 'imported')
            if len(self.surveys) == 1:
                ltime = len(self.surveys[0].df)
            if len(self.surveys) > 1:
                if len(self.surveys[-1].df) != ltime:
                    print('ERROR:', f, 'survey doesn\'t have the same length')
                    return
        self.iTimeLapseReciprocal = np.array(self.iTimeLapseReciprocal)
        self.elec = self.surveys[0].elec
        
        
        # create bigSurvey
        print('creating bigSurvey')
        self.bigSurvey = Survey(os.path.join(dirname, files[0]), ftype=ftype)
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
        self.pseudo = self.surveys[0].pseudo # just display first pseudo section
            
        self.plotError = self.bigSurvey.plotError
        self.linfit = self.bigSurvey.linfit
        self.lmefit = self.bigSurvey.lmefit
        self.pwlfit = self.bigSurvey.pwlfit
        self.phaseplotError = self.bigSurvey.phaseplotError
        self.plotIPFit = self.bigSurvey.plotIPFit
    
#    def pseudo(self, **kwargs):
#        """ Create a pseudo section.
#        
#        Parameters
#        ----------
#        **kwargs :
#            To be passed to `R2.pseudoCallback()`.
#        """
#        if self.iBorehole == True:
#            print('NOT PLOTTING PSEUDO FOR BOREHOLE FOR NOW')
#        else:
#            self.pseudoCallback(**kwargs)
    
    def createMesh(self, typ='default', buried=None, surface=None, cl_factor=2,
                   cl=-1, dump=print, **kwargs):
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
            Array with two columns x and y for additional surface points.
        cl_factor : float, optional
            Characteristic length factor. Only used for triangular mesh to allow
            mesh to be refined close the electrodes and then expand.
        cl : float, optional
            Characteristic length that define the mesh size around the
            electrodes.
        dump : function, optional
            Function to which pass the output during mesh generation. `print()`
             is the default.
        """
        if typ == 'default':
            if self.elec[:,2].sum() == 0:
                typ = 'quad'
                print('Using a quadrilateral mesh')
            else:
                typ = 'trian'
                print('Using a triangular mesh')
        if typ == 'quad':
#            mesh = QuadMesh(elec, nnode=4)
            elec_x = self.elec[:,0]
            elec_y = self.elec[:,1]
            mesh,meshx,meshy,topo,e_nodes = mt.quad_mesh(elec_x,elec_y)#,**kwargs)
#            mesh = QuadMesh()
#            meshx, meshy, topo, e_nodes = mesh.createMesh(elec=self.elec, **kwargs)            
            self.param['meshx'] = meshx
            self.param['meshy'] = meshy
            self.param['topo'] = topo
            self.param['mesh_type'] = 4
            self.param['node_elec'] = np.c_[1+np.arange(len(e_nodes)), e_nodes, np.ones(len(e_nodes))].astype(int)
            if 'regions' in self.param: # allow to create a new mesh then rerun inversion
                del self.param['regions']
            if 'num_regions' in self.param:
                del self.param['num_regions']
        if typ == 'trian':
            elec = self.elec.copy()
            geom_input = {}
            if (buried is not None
                    and elec.shape[0] == len(buried)
                    and np.sum(buried) != 0):
                iburied = buried == True
                geom_input['electrode'] = [self.elec[~iburied, 0],
                                           self.elec[~iburied, 1]]
                geom_input['buried1'] = [self.elec[iburied, 0],
                                       self.elec[iburied, 1]]
            else:
                geom_input['electrode'] = [self.elec[:,0],
                                           self.elec[:,1]]
            if surface is not None:
                geom_input['surface'] = [surface[:,0], surface[:,1]]
            mesh = tri_mesh(geom_input,
                             path=os.path.join(self.apiPath, 'exe'),
                             cl_factor=cl_factor,
                             cl=cl, dump=dump, show_output=True)
            self.param['mesh_type'] = 3
#            self.param['num_regions'] = len(mesh.regions)
#            regs = np.array(np.array(mesh.regions))[:,1:]
#            if self.typ == 'R2':
#                regions = np.c_[regs, np.ones(regs.shape[0])*50]
#            if self.typ == 'cR2':
#                regions = np.c_[regs, np.ones(regs.shape[0])*50, np.ones(regs.shape[0])*-0.1]
#            self.param['regions'] = regions
            self.param['num_xy_poly'] = 5
            # define xy_poly_table (still need to do it here because param['meshx'] is undefined if triangular mesh)
            doi = np.abs(self.elec[0,0]-self.elec[-1,0])/2
            ymax = np.max(self.elec[:,1])
            ymin = np.min(self.elec[:,1])-doi
            xy_poly_table = np.array([
            [self.elec[0,0], ymax],
            [self.elec[-1,0], ymax],
            [self.elec[-1,0], ymin],
            [self.elec[0,0], ymin],
            [self.elec[0,0], ymax]])
            self.param['xy_poly_table'] = xy_poly_table
            e_nodes = np.arange(len(self.elec))+1
            self.param['node_elec'] = np.c_[1+np.arange(len(e_nodes)), e_nodes, np.ones(len(e_nodes))].astype(int)
        self.mesh = mesh
        self.param['mesh'] = mesh
        
        self.param['num_regions'] = 0
        self.param['res0File'] = 'res0.dat'
        
        numel = self.mesh.num_elms
        self.mesh.add_attribute(np.ones(numel)*100, 'res0') # default starting resisivity [Ohm.m]
        self.mesh.add_attribute(np.ones(numel, dtype=int), 'zones')
        self.mesh.add_attribute(np.zeros(numel, dtype=bool), 'fixed')
        self.mesh.add_attribute(np.zeros(numel, dtype=float), 'iter')
        #write mesh to working directory - edit by jamyd91 
        file_path = os.path.join(self.dirname,'mesh.dat')
        self.mesh.write_dat(file_path)
        
        
        self.regid = 0
        self.regions = np.zeros(len(self.mesh.elm_centre[0]))
        self.resist0 = np.ones(len(self.regions))*100
        
        
    def showMesh(self, ax=None):
        """ Display the mesh.
        """
        if self.mesh is None:
            raise Exception('Mesh undefined')
        else:
#            xlim = (np.min(self.elec[:,0]-20, np.max(self.elec[:,0])))
#            ylim = (0, 110) # TODO
#            self.mesh.show(xlim=xlim, ylim=ylim) # add ax argument
            self.mesh.show(ax=ax, color_bar=False)
    
    def write2in(self, param={}, typ=''):
        """ Create configuration file for inversion.
        
        Parameters
        ----------
        param : dict
            Dictionnary of parameters and values for the inversion settings.
        typ : str, optional
            Type of inversion. By default given by `R2.typ`.
        """
        if typ == '':
            typ = self.typ
        if all(self.surveys[0].df['irecip'].values == 0):
            if 'a_wgt' not in self.param:
                self.param['a_wgt'] = 0.01
            if 'b_wft' not in self.param:
                self.param['b_wgt'] = 0.02
        if typ == 'cR2':
            if self.errTypIP != 'none': # we have individual errors
                if 'b_wgt' not in self.param:
                    self.param['b_wgt'] = 0
                if 'c_wgt' not in self.param:
                    self.param['c_wgt'] = 0
                if 'a_wgt' not in self.param:
                    self.param['a_wgt'] = 0.01 # not sure of that (Gui)
            else:
                if 'c_wgt' not in self.param:
                    self.param['c_wgt'] = 1 # better if set by user !!
                if 'd_wgt' not in self.param:
                    self.param['d_wgt'] = 2
        
        if self.param['mesh_type'] == 4:
            self.param['zones'] = self.mesh.attr_cache['zones']
            #TODO reshape it to the form of the mesh
            
        # all those parameters are default but the user can change them and call
        # write2in again
        for p in param:
            self.param[p] = param[p]
        
        # overwite the cropping of the mesh, too many issue (though we need it
        # to display it in paraview)
        self.param['num_xy_poly'] = 0
        
        if self.iTimeLapse == True:
            refdir = os.path.join(self.dirname, 'ref')
            if os.path.exists(refdir) == False:
                os.mkdir(refdir)
            param = self.param.copy()
            param['a_wgt'] = 0.01
            param['b_wgt'] = 0.02
            param['num_xy_poly'] = 0
            param['reg_mode'] = 0 # set by default in ui.py too
            self.configFile = write2in(param, refdir, typ=typ)
            param = self.param.copy()
            param['num_regions'] = 0
            param['reg_mode'] = 2
            param['res0File'] = 'Start_res.dat'
            write2in(param, self.dirname, typ=typ)
        else:
            self.configFile = write2in(self.param, self.dirname, typ=typ)
        
        # write the res0.dat needed for starting resistivity
        if self.iForward is True: # we will invert results from forward
            # inversion so we need to start from a homogeneous model
            print('Setting a homogeneous background model as the survey to \
                  be inverted is from a forward model already.')
            res0 = np.ones(self.mesh.num_elms)*100 # default starting resistivity [Ohm.m]
            self.mesh.add_attribute(res0, 'r100')
            self.mesh.write_attr('r100', 'res0.dat', self.dirname)
        else: # if we invert field data, we allow the user to define prior
            # knowledge of the resistivity structure
            self.mesh.write_attr('res0', 'res0.dat', self.dirname)
        
        # rewriting mesh.dat
        paramFixed = 1+ np.arange(self.mesh.num_elms)
        ifixed = np.array(self.mesh.attr_cache['fixed'])
        paramFixed[ifixed] = 0
        
#        print('SUM OF PARAMFIXED = ', np.sum(paramFixed == 0))
        self.mesh.write_dat(os.path.join(self.dirname, 'mesh.dat'),
                            zone = self.mesh.attr_cache['zones'],
                            param = paramFixed)

        # NOTE if fixed elements, they need to be at the end !
        def readMeshDat(fname):
            with open(fname, 'r') as f:
                x = f.readline().split()
            numel = int(x[0])
            elems = np.genfromtxt(fname, skip_header=1, max_rows=numel)
            nodes = np.genfromtxt(fname, skip_header=numel+1)
            return elems, nodes
        
        def writeMeshDat(fname, elems, nodes):
            numel = len(elems)
            nnodes = len(nodes)
            with open(fname, 'w') as f:
                f.write('{:.0f} {:.0f}\n'.format(numel, nnodes))
            with open(fname, 'a') as f:
                np.savetxt(f, elems, fmt='%.0f')
                np.savetxt(f, nodes, fmt='%.0f %f %f')
        
        meshFile = os.path.join(self.dirname, 'mesh.dat')
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
        
        # UPDATE MESH
        self.mesh = dat_import(meshFile)
        
        # We NEED parameter = elemement number if changed AND
        # all fixed element (with parameter = 0) at the end of the element
        # matrix. BOTH conditions needs to be filled.
    

    def write2protocol(self, errTyp='', errTypIP='', errTot=False, **kwargs):
        """ Write a protocol.dat file for the inversion code.
        
        Parameters
        ----------
        errTyp : str
            Type of the DC error. Either 'pwl', 'lin', 'obs'.
        errTypIP : str
            Type of the IP error. Either 'pwl'.
        """
        if self.typ == 'R2':
            ipBool = False
        elif self.typ == 'cR2':
            ipBool = True
        else:
            print('NOT IMPLEMENTED YET')

        if errTyp == '':
            errTyp = self.errTyp
#        if ipBool == True:
        if errTypIP == '':
            errTypIP = self.errTypIP
        
        # important changing sign of resistivity and quadrupoles so to work
        # with complex resistivity
        if self.typ[0] == 'c':
            for s in self.surveys:
                ie = s.df['resist'].values < 0
                m = s.df['m'].values.copy()
                n = s.df['n'].values.copy()
                s.df.loc[ie, 'm'] = n[ie]
                s.df.loc[ie, 'n'] = m[ie]
                s.df.loc[ie, 'resist'] = s.df.loc[ie, 'resist'].values*-1
                s.df.loc[ie, 'recipMean'] = s.df.loc[ie, 'recipMean'].values*-1

        if self.iTimeLapse is True:
            # a bit simplistic but assign error to all based on Transfer resistance
#            allHaveReciprocal = all(self.iTimeLapseReciprocal == True)
            # let's assume it's False all the time for now
            content = ''
            for i, s in enumerate(self.surveys[1:]):
                content = content + str(len(s.df)) + '\n'
                s.df['resist0'] = self.surveys[0].df['resist']
                if errTyp != 'none': # there is an error model
                    s.df['error'] = self.bigSurvey.errorModel(s.df['resist'].values)
                    s.df['index'] = np.arange(1, len(s.df)+1)
                    content = content + s.df[['index','a','b','m','n','resist', 'resist0','error']].to_csv(sep='\t', header=False, index=False)
                else:
                    s.df['index'] = np.arange(1, len(s.df)+1)
                    content = content + s.df[['index','a','b','m','n','resist', 'resist0']].to_csv(sep='\t', header=False, index=False)
                if i == 0:
                    refdir = os.path.join(self.dirname, 'ref')
                    if os.path.exists(refdir) == False:
                        os.mkdir(refdir)
                    if 'mesh.dat' in os.listdir(self.dirname):
                        shutil.copy(os.path.join(self.dirname, 'mesh.dat'),
                                os.path.join(self.dirname, 'ref', 'mesh.dat'))
                    with open(os.path.join(refdir, 'protocol.dat'), 'w') as f:
                        f.write(content) # write the protocol for the reference file
            with open(os.path.join(self.dirname, 'protocol.dat'), 'w') as f:
                f.write(content)
                
        elif self.iBatch is True:
            print('print in iBatch format')
            content = ''
            for i, s in enumerate(self.surveys):
                if self.errTyp != 'none':
                    s.df[self.errTyp+'Error'] = self.bigSurvey.errorModel(s.df['resist'].values)
                df = s.write2protocol(outputname='',
                    errTyp=errTyp, ip=ipBool, errTypIP=errTypIP, errTot=errTot)
                content = content + str(len(df)) + '\n'
                content = content + df.to_csv(sep='\t', header=False, index=False)
#            content = ''
#            for i, s in enumerate(self.surveys[0:]):
#                #TODO only take mean of reciprocal no ?
#                content = content + str(len(s.df)) + '\n'
#                if errTyp != 'none': # there is an error model
#                    s.df['error'] = self.bigSurvey.errorModel(s.df['resist'].values)
#                    s.df['index'] = np.arange(1, len(s.df)+1)
#                    content = content + s.df[['index','a','b','m','n','resist','error']].to_csv(sep='\t', header=False, index=False)
#                else:
#                    s.df['index'] = np.arange(1, len(s.df)+1)
#                    content = content + s.df[['index','a','b','m','n','resist']].to_csv(sep='\t', header=False, index=False)
            with open(os.path.join(self.dirname, 'protocol.dat'), 'w') as f:
                f.write(content)
        
        else:
            self.surveys[0].write2protocol(os.path.join(self.dirname, 'protocol.dat'),
                    errTyp=errTyp, ip=ipBool, errTypIP=errTypIP, errTot=errTot)
        
        
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
        targetName = os.path.join(dirname, exeName)
        
        # copy R2.exe
        if ~os.path.exists(targetName):
            shutil.copy(os.path.join(self.apiPath, 'exe', exeName), targetName)  
        
            
        if OS == 'Windows':
            cmd = [exeName]
        else:
            cmd = ['wine',exeName]
            
#        p = Popen(cmd, stdout=PIPE, shell=False)
#        while p.poll() is None:
#            line = p.stdout.readline().rstrip()
#            dump(line.decode('utf-8'))
        
        def execute(cmd):
            popen = Popen(cmd, stdout=PIPE, shell=False, universal_newlines=True)
            for stdout_line in iter(popen.stdout.readline, ""):
                yield stdout_line
            popen.stdout.close()
            return_code = popen.wait()
            if return_code:
                print('error on return_code')
        
        for text in execute(cmd):
                dump(text.rstrip())

        os.chdir(cwd)
        
        
    def invert(self, param={}, iplot=False, dump=print, modErr=False):
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
        dump('done\n')
        
        # runs inversion
        if self.iTimeLapse == True:
            dump('---------- Inverting background/reference model ---------\n')
            refdir = os.path.join(self.dirname, 'ref')
            shutil.move(os.path.join(self.dirname,'res0.dat'),
                        os.path.join(refdir, 'res0.dat'))
            self.runR2(refdir, dump=dump)
            shutil.copy(os.path.join(refdir, 'f001_res.dat'),
                    os.path.join(self.dirname, 'Start_res.dat'))
        dump('-------- Main inversion ---------------\n')
        self.runR2(dump=dump)
        
        if iplot is True:
            self.showResults()

    
    def showResults(self, index=0, ax=None, edge_color='none', attr='', sens=True, color_map='viridis', **kwargs):
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
        """
        if (attr == '') & (self.typ == 'R2'):
            attr = 'Resistivity(log10)'
        if (attr == '') & (self.typ == 'cR2'):
            attr = 'Sigma_real(log10)'
        if len(self.meshResults) == 0:
            self.getResults()
        if len(self.meshResults) > 0:
            self.meshResults[index].show(ax=ax, edge_color=edge_color, attr=attr, sens=sens, color_map=color_map, **kwargs)
        else:
            print('Unexpected Error')

    
    def getResults(self):
        """ Collect inverted results after running the inversion and adding
        them to `R2.meshResults` list.
        """
        if self.typ == 'R2':
            if self.iTimeLapse == True:
                fresults = os.path.join(self.dirname, 'ref', 'f001_res.vtk')
                print('reading ref', fresults)
                mesh = mt.vtk_import(fresults)
                mesh.elec_x = self.elec[:,0]
                mesh.elec_y = self.elec[:,1]
                self.meshResults.append(mesh)
            for i in range(100):
                fresults = os.path.join(self.dirname, 'f' + str(i+1).zfill(3) + '_res.vtk')
                if os.path.exists(fresults):
                    print('reading ', fresults)
                    mesh = mt.vtk_import(fresults)
                    mesh.elec_x = self.elec[:,0]
                    mesh.elec_y = self.elec[:,1]
                    self.meshResults.append(mesh)
                else:
                    break
        if self.typ == 'cR2':
            fresults = os.path.join(self.dirname, 'f001.vtk')
            print('reading ref', fresults)
            mesh = mt.vtk_import(fresults)
            mesh.elec_x = self.elec[:,0]
            mesh.elec_y = self.elec[:,1]
            self.meshResults.append(mesh)

            
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
        print('showSection called')
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
        #            xs = sen[:,0].reshape((leny, lenx), order='F')
        #            ys = sen[:,1].reshape((leny, lenx), order='F')
            zs = sen[:,2].reshape((leny, lenx), order='F')
            zs = np.log10(zs)
            zs -= np.min(zs)
            alpha = zs/np.max(zs)
        #            alpha[alpha < 0] = 0
            print(np.max(alpha), np.min(alpha))
        if ilog10:
            z = np.log10(z)
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()
        cax = ax.pcolormesh(x, y, z)
        ax.plot(self.elec[:,0], self.elec[:,1], 'ko')
    #    fig.canvas.draw() # need to draw the figure to have the cax.get_facecolors()
    #    print(cax.get_facecolors().shape)
    #    print(alpha.flatten().shape)
    #    for a in cax.get_facecolors():
    #        a[3] = 0
        #for a, b in zip(cax.get_facecolors(), alpha.flatten()):
        #    a[3] = 0.5
        #    print(a)
    #    fig.canvas.draw()
        cbar = fig.colorbar(cax, ax=ax)
        if ilog10:
            cbar.set_label(r'$\log_{10}(\rho) [\Omega.m]$')
        else:
            cbar.set_label(r'$\rho [\Omega.m]$')
        ax.set_ylabel('Depth [m]')
        ax.set_xlabel('Distance [m]')
#        fig.tight_layout()
    #    fig.show()
#        return fig
    
    def addRegion(self, xy, res0, blocky=False, fixed=False, ax=None):
        """ Add region according to a polyline defined by `xy` and assign it
        the starting resistivity `res0`.
        
        Parameters
        ----------
        xy : array
            Array with two columns for the x and y coordinates.
        res0 : float
            Resistivity values of the defined area.
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
        selector = SelectPoints(ax, np.array(self.mesh.elm_centre).T,
                                typ='poly')
        selector.setVertices(xy)
        selector.getPointsInside()
        idx = selector.iselect
        self.regid = self.regid + 1
        self.regions[idx] = self.regid
        self.mesh.cell_attributes = list(self.regions) # overwriting regions
        self.resist0[idx] = res0
        self.mesh.attr_cache['res0'] = self.resist0 # hard way to do it
        
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
        self.regid = 0
        self.regions.fill(0)
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
        # regions definition
#        self.mesh.add_attribute()
        def callback(idx):
            print('nb elements selected:', np.sum(idx))
#                res = input('Input the resistivity [Ohm.m] of the section:')
#                print(res)
            self.regid = self.regid + 1
            self.regions[idx] = self.regid
            
            self.mesh.cell_attributes = list(self.regions) # overwritin regions
#            res0list = np.unique(self.regions) # TODO ask for resistvity in the table
#            self.mesh.assign_zone_attribute(self.regions,res0list,'res0')
            
            self.mesh.draw()
            if addAction is not None:
                addAction()
        self.mesh.show(ax=ax)
        # we need to assign a selector to self otherwise it's not used
        self.selector = SelectPoints(ax, np.array(self.mesh.elm_centre).T,
                                     typ=typ, callback=callback)
        if ax is None:
            return fig
            
        
    def assignRes0(self, regionValues={}, zoneValues={}, fixedValues={}):
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
        
#        for key in regionValues.keys():
#            if key in self.regions:
#                res0list.append(regionValues[key])
#            else:
#                res0list.append(100) # default resistivity
#
#        zoneList = []
#        for key in zoneValues.keys():
#            if key in self.regions:
#                zoneList.append(zoneValues[key])
#            else:
#                zoneList.append(1) # default zone
#        
#        fixedList = []
#        for key in fixedValues.keys():
#            if key in self.regions:
#                fixedList.append(fixedValues[key])
#            else:
#                fixedList.append(False) # default fixed
#
#        self.mesh.assign_zone_attribute(self.regions, res0list, 'res0')
#        self.mesh.assign_zone_attribute(self.regions, zoneList, 'zones')
#        #if self.param['mesh_type'] == 3:
#        
#        print('fixed = in assignRes0', fixedList, zoneList)
#        self.mesh.assign_zone_attribute(self.regions, fixedList, 'fixed')
        

    def createSequence(self, skipDepths=[(0, 10)]):
        """ Create a dipole-dipole sequence.
        
        Parameters
        ----------
        skipDepths : list of tuple, optional
            Each tuple in the list is of the form `(skip, depths)`. The `skip` is the number of electrode between the A B and M N electrode. The `depths` is the number of quadrupole which will have the same current electrode (same A B). The higher this number, the deeper the investigation.
        """
        qs = []
        nelec = len(self.elec)
        for sd in skipDepths:
            elec, quad = ddskip(nelec, spacing=1, skip=int(sd[0]), depth=int(sd[1]))
            qs.append(quad)
        qs = np.vstack(qs)
        self.sequence = qs
    
    
    def importElec(self, fname=''):
        """ Import electrodes positions.
        
        Parameters
        ----------
        fname : str
            Path of the CSV file containing the electrodes positions. It should contains 3 columns maximum with the X, Y, Z positions of the electrodes.
        """
        elec = pd.read_csv(fname)
        if elec.shape[1] > 3:
            raise ValueError('The file should have no more than 3 columsn')
        else:
            self.elec = elec
            
    
    def importSequence(self, fname=''):
        """ Import sequence for forward modelling.
        
        Parameters
        ----------
        fname : str
            Path of the CSV file to be imported. The file shouldn't have any headers just 4 columns with the 4 electrodes numbers.
        """
        seq = pd.read_csv(fname)
        if seq.shape[1] != 4:
            raise ValueError('The file should be a CSV file wihtout headers with exactly 4 columns with electrode numbers.')
        else:
            self.sequence = seq
    
    
    def forward(self, noise=0.00, iplot=False, dump=print):
        """ Operates forward modelling.
        
        Parameters
        ----------
        noise : float, optional 0 <= noise <= 1
            Noise level from a Gaussian distribution that should be applied
            on the forward apparent resistivities obtained. 
        iplot : bool, optional
            If `True` will plot the pseudo section after the forward modelling.
        dump : function, optional
            Function to print information messages when running the forward model.
        """
        fwdDir = os.path.join(self.dirname, 'fwd')
        if os.path.exists(fwdDir):
            shutil.rmtree(fwdDir)
        os.mkdir(fwdDir)
        
        # write the resistivity.dat
#        centroids = np.array(self.mesh.elm_centre).T
#        resFile = np.zeros((centroids.shape[0],3)) # centroix x, y, z, res0
##        resFile[:,:centroids.shape[1]] = centroids
#        resFile[:,-1] = self.resist0
#        np.savetxt(os.path.join(fwdDir, 'resistivity.dat'), resFile,
#                   fmt='%.3f')
        
        self.mesh.write_attr('res0', 'resistivity.dat', fwdDir)
        
        if os.path.exists(os.path.join(self.dirname, 'mesh.dat')) is True:
            shutil.copy(os.path.join(self.dirname, 'mesh.dat'),
                        os.path.join(fwdDir, 'mesh.dat'))
        
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
            return np.random.normal(x,level,1)
        addnoise = np.vectorize(addnoise)
        self.noise = noise
        
        elec = self.elec.copy()
        self.surveys = [] # need to flush it (so no timeLapse forward)
        self.createSurvey(os.path.join(fwdDir, 'R2_forward.dat'), ftype='Protocol')
        self.surveys[0].df['resist'] = addnoise(self.surveys[0].df['resist'], noise)
        self.elec = elec
        
        self.write2protocol()
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
            shutil.copy(os.path.join(self.dirname, 'mesh.dat'),
                        os.path.join(fwdDir, 'mesh.dat'))
            fparam['num_regions'] = 0
            fparam['res0File'] = 'resistivity.dat'
        write2in(fparam, fwdDir, typ=self.typ)
        
        # write the protocol.dat based on measured sequence
        seq = self.surveys[0].df[['a','b','m','n']].values
        protocol = pd.DataFrame(np.c_[1+np.arange(seq.shape[0]),seq])
        outputname = os.path.join(fwdDir, 'protocol.dat')
        with open(outputname, 'w') as f:
            f.write(str(len(protocol)) + '\n')
        with open(outputname, 'a') as f:
            protocol.to_csv(f, sep='\t', header=False, index=False)
    
        # fun the inversion
        self.runR2(fwdDir) # this will copy the R2.exe inside as well
        
        # get error model
        x = np.genfromtxt(os.path.join(fwdDir, 'R2_forward.dat'), skip_header=1)
        modErr = np.abs(100-x[:,-1])/100
        self.surveys[0].df['modErr'] = modErr
        
        # eventually delete the directory to space space
        
        
    
    def showIter(self, ax=None, index=-1):
        """ Dispay temporary inverted section after each iteration.
        
        Parameters
        ----------
        ax : matplotib axis, optional
            If specified, the graph will be plotted along `ax`.
        index : int, optional
            Iteration number to show.
        """
        files = os.listdir(self.dirname)
        fs = []
        for f in files:
            if (f[-8:] == '_res.dat') & (len(f) == 16):
                fs.append(f)
        fs = sorted(fs)
        print(fs)
        if len(fs) > 0:
            if self.param['mesh_type'] == 10:
                self.showSection(os.path.join(self.dirname, fs[index]), ax=ax)
                # TODO change that to full meshTools
                
            else:
                x = np.genfromtxt(os.path.join(self.dirname, fs[index]))
#                iterNumber = fs[-1].split('_')[0].split('.')[1]
#                attrName = '$log_{10}(\rho)$ [Ohm.m] (iter {:.0f})'.format(iterNumber) # not sure it is log10
#                print('iterNumber = ', iterNumber, 'name=', attrName)
#                self.mesh.add_attr_dict({'iter':x[:,-2]})
                self.mesh.attr_cache['iter'] = x[:,-2]
#                self.mesh.draw(attr=attrName)
                self.mesh.show(ax=ax, attr='iter', edge_color='none', color_map='viridis')
                
               
    def pseudoError(self, ax=None):
        """ Plot pseudo section of errors from file `f001_err.dat`.
        
        Parameters
        ----------
        ax : matplotlib axis
            If specified, the graph will be plotted against `ax`.
        """
        if self.typ == 'R2':
            err = np.genfromtxt(os.path.join(self.dirname, 'f001_err.dat'), skip_header=1)        
            array = err[:,[-2,-1,-4,-3]].astype(int)
            errors = err[:,0]
        if self.typ == 'cR2':
            err = pd.read_fwf(os.path.join(self.dirname, 'f001.err')).values       
            array = err[:,[0,1,2,3]].astype(int)
            errors = err[:,4]
        
        spacing = np.diff(self.elec[[0,1],0])
        pseudo(array, errors, spacing, ax=ax, label='Normalized Errors', log=False, geom=False, contour=False)
    
    def showInversionErrors(self, ax=None):
        """ Display inversion error by measurment numbers.
        """
        # TODO make it work with cR2
        file_path = os.path.join(self.dirname, 'f001_err.dat')
        error_info = import_R2_err_dat(file_path)
        disp_R2_errors(error_info, ax=ax)



def pseudo(array, resist, spacing, label='', ax=None, contour=False, log=True, geom=True):
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

        
#%% test code
#os.chdir('/media/jkl/data/phd/tmp/r2gui/')
#k = R2('/media/jkl/data/phd/tmp/r2gui/api/invdir')
#k.typ = 'cR2'
#k = R2
#k.createSurvey('api/test/syscalFile.csv', ftype='Syscal')
#k.createMesh(typ='trian')
#k.write2in()
#k.computeModelError()
#k.pwlfit()
#k.write2protocol(errTyp='pwl', errTot=True)
#k.invert(modErr=True)
#k.createSurvey('api/test/rifleday8.csv', ftype='Syscal')
#k.pwlfit()
#k.invert(iplot=True)
#k.showIter()
#k.showResults()
#k.surveys[0].dca()
#k.pseudo(contour=True)
#k.linfit(iplot=True)
#k.pwlfit()
#k.errTyp='obs'
#k.lmefit(iplot=True)
#k.createMesh(typ='quad', elemx=8)
#k.createMesh(typ='trian')
#k.mesh.show()
#fig, ax = plt.subplots()
#fig.suptitle('kkk')
#k.mesh.show(ax=ax)
#k.write2in()
#k.plotIPFit()
#k.errTyp = 'pwl'
#k.errTypIP = 'pwl'
#k.invert(iplot=True)
#k.showIter()
#k.showResults(edge_color='k')
#k.pseudoError()
#k.showSection()
#fig, ax = plt.subplots()
#fig.suptitle('hkk')
#k.showResults()
#k.showResults(edge_color='none', sens=True)
#k.showResults(attr=attr[0])
#fig, ax = plt.subplots()
#fig.suptitle('kkk')
#k.showResults(ax=ax)
#print(os.path.dirname(os.path.realpath(__file__)))


#fresults = os.path.join('./test/f001_res.vtk')
#if os.path.isfile(fresults):
#    print('kk')
#    mesh_dict=mt.vtk_import(fresults)#makes a dictionary of a mesh 
#    mesh = Mesh_obj.mesh_dict2obj(mesh_dict)# this is a mesh_obj class instance 
#    mesh.show()
#

#%% test for DC with topo
#k = R2()
#k.createSurvey('api/test/syscalFileTopo.csv')
#topo = np.genfromtxt('api/test/elecTopo.csv', delimiter=',')
#k.elec = topo
#k.createMesh(typ='trian',cl=0.1, cl_factor=5)
#k.showMesh()
#k.invert()
#
#k.showIter(index=0)

#%% check order by building centroid


    
#%% test for IP
#os.chdir('/media/jkl/data/phd/tmp/r2gui/')
#k = R2()
#k.typ = 'cR2'
#k.createSurvey('api/test/rifleday8.csv', ftype='Syscal')
#k.createSurvey('api/test/syscalFileIP.csv')
#k.write2protocol()
#k.invert()
#k.pseudoError()

#%% test for timelapse inversion
#os.chdir('/media/jkl/data/phd/tmp/r2gui/')
#k = R2()
#k.createTimeLapseSurvey(os.path.join(k.dirname, '../test/testTimelapse'))
#k.createBatchSurvey(os.path.join(k.dirname, '../test/testTimelapse'))
#k.linfit()
#k.pwlfit()
#k.errTyp = 'pwl'
#k.param['a_wgt'] = 0
#k.param['b_wgt'] = 0
#k.createMesh(typ='trian', cl_factor=5, cl=0.05)
#k.showMesh()
#k.write2in()
#k.write2protocol()
#k.invert(iplot=False)
#k.showResults(index=0)
#k.showResults(index=1)
#k.showSection(os.path.join(k.dirname, 'f001_res.vtk'))
#k.showSection(os.path.join(k.dirname, 'f002_res.vtk'))

#%% test mesh
#os.chdir('/media/jkl/data/phd/tmp/r2gui/')
#k = R2()
#k.createSurvey('api/test/syscalFile.csv', ftype='Syscal')
#k.elec[:,1] = np.tile(np.arange(0,-12,-1),2)
#k.elec[:,0] = np.repeat([0,8], 12)
#k.elec[1:11,:2] = np.c_[np.ones(10)*2, np.linspace(0,-4,10)]
#buried = np.ones(k.elec.shape[0], dtype=bool)
#buried[[0,12]] = False
#surface = np.array([[-2,0],[10,0]])
#k.createMesh(typ='trian', buried=buried, cl=0.5, cl_factor=5) # works well
#k.showMesh()
#k.invert()


#%% forward modelling
#def readMeshDat(fname):
#    with open(fname, 'r') as f:
#        x = f.readline().split()
#    numel = int(x[0])
#    numnode = int(x[1])
#    elems = np.genfromtxt(fname, skip_header=1, max_rows=numel)
#    nodes = np.genfromtxt(fname, skip_header=numel+1)
#    return elems, nodes

#os.chdir('/media/jkl/data/phd/tmp/r2gui/')
#plt.close('all')
#k = R2()
#k.createSurvey('api/test/syscalFile.csv')
#k.elec = np.c_[np.linspace(0,5.75, 24), np.zeros((24, 2))]
#k.createMesh(typ='trian')
#
## full API function
#k.addRegion(np.array([[1,0],[2,0],[2,-0.5],[1,-0.5],[1,0]]), 10)
#k.addRegion(np.array([[3,-0.5],[3.5,-0.5],[3.5,-1],[3,-1],[3,-0.5]]), 20, blocky=True, fixed=True)
#k.addRegion(np.array([[4,0],[5,0],[5,-0.5],[4,-0.5],[4,0]]), 30, blocky=True, fixed=False)

## full GUI function
#k.createModel() # manually define 3 regions
#k.assignRes0({1:500, 2:20, 3:30}, {1:1, 2:2, 3:1}, {1:False, 2:False, 3:True})

#k.forward(iplot=True, noise=0.0)
##k.resetRegions() # don't need to do this anymore you need to reset regions as they are used for starting model
#k.invert(iplot=False)
#k.showResults(attr='Resistivity(Ohm-m)', sens=False)


#%% test Sina
#os.chdir('/media/jkl/data/phd/tmp/r2gui/')
#k = R2('/media/jkl/data/phd/tmp/r2gui/api/invdir/')
#k.createSurvey('/home/jkl/Downloads/D10_N.csv')
#k.surveys[0].addData('/home/jkl/Downloads/D10_R.csv')
#k.invert()


#%% test Paul River
#k = R2()
#k.createSurvey('api/test/primeFile.dat', ftype='BGS Prime')
#elec = np.genfromtxt('api/test/primePos.csv', delimiter=',')
#k.elec[:,0] = elec[:,0]
#k.elec[:,1] = elec[:,1]
#k.elec = elec
#k.surveys[0].manualFiltering()
#k.createMesh(typ='quad', elemx=4)
#k.showMesh()


#%% check mesh
#k.createMesh('trian')
#k.param['num_xy_poly'] = 0
#k.invert()

##%% test mesh issue
#def readMeshDat(fname):
#    with open(fname, 'r') as f:
#        x = f.readline().split()
#    numel = int(x[0])
#    elems = np.genfromtxt(fname, skip_header=1, max_rows=numel).astype(int)
#    nodes = np.genfromtxt(fname, skip_header=numel+1)
#    return elems, nodes
#
#def computeCentroid(elems, nodes):
#    elx = nodes[elems,1]
#    ely = nodes[elems,2]
#    cx = np.sum(elx, axis=1)/3
#    cy = np.sum(ely, axis=1)/3
#    return np.c_[cx, cy]
#    
#elems, nodes = readMeshDat('api/invdir/mesh.dat')
#centroids = computeCentroid(elems[:,1:4]-1, nodes)
#np.savetxt('api/invdir/centroids.dat', centroids)
#df = pd.read_fwf('api/invdir/f001.001_res.dat', header=None)
#df2 = pd.read_fwf('api/invdir/f001_res.dat', header=None)
#x = df.values[:,:2] - centroids
#print(np.sum(x, axis=0))
#
#x2 = df.values[:,:2] - df2.values[:,:2]
#np.sum(x2, axis=0) # ok same format for iteration and non iteration


#%% mesh investigation part 2
#k = R2()
#k.createSurvey('api/test/syscalFile.csv')
#k.createMesh('quad')
#k.showMesh()
#k.param['num_xy_poly'] = 0
#k.invert()
#k.showResults(contour=True)
#k.showIter(index=0)
#mesh2 = mt.dat_import('api/invdir/mesh.dat')
#mesh2.write_dat('api/invdir/mesh2.dat')
#
#elems, nodes = readMeshDat('api/invdir/mesh.dat')
#elems2, nodes2 = readMeshDat('api/invdir/mesh2.dat')
