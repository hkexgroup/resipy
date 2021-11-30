#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import sys
sys.setrecursionlimit(10000) # to prevent the display issue if so many measurements are rejected at once.
import time
print(
'''
============================================================
d8888b. d88888b .d8888. d888888b d8888b. db    db
88  `8D 88'     88'  YP   `88'   88  `8D `8b  d8'
88oobY' 88ooooo `8bo.      88    88oodD'  `8bd8' 
88`8b   88~~~~~   `Y8b.    88    88~~~      88   
88 `88. 88.     db   8D   .88.   88         88   
88   YD Y88888P `8888Y' Y888888P 88         YP   
============================================================
''')
from PyQt5.QtWidgets import (QMainWindow, QSplashScreen, QApplication, QPushButton, QWidget,
    QTabWidget, QVBoxLayout, QGridLayout, QLabel, QLineEdit, QMessageBox, QSplitter,
    QFileDialog, QCheckBox, QComboBox, QTextEdit, QSlider, QHBoxLayout, QFrame, 
    QTableWidget, QFormLayout, QTableWidgetItem, QHeaderView, QProgressBar, QDialog,
    QStackedLayout, QRadioButton, QGroupBox, QTextBrowser, QMenu)#, QAction, qApp, QButtonGroup, QListWidget, QShortcut)
from PyQt5.QtGui import QIcon, QPixmap, QIntValidator, QDoubleValidator, QColor, QPalette#, QKeySequence
from PyQt5.QtCore import QThread, pyqtSignal, QTimer, QUrl, QObject#, QProcess#, QSize
from PyQt5.QtCore import Qt
from functools import partial
QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True) # for high dpi display
QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

import threading
import traceback

# COMMENT below is using UI automatic testing
# from resipy.Project import ResIPy_version
# import matplotlib
# matplotlib.use('Qt5Agg')
# from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
# from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
# from matplotlib.figure import Figure
# import numpy as np
# import pandas as pd
# from datetime import datetime
# from matplotlib import rcParams
# rcParams.update({'font.size': 11}) # CHANGE HERE for graph font size
# import platform
# OS = platform.system()
# from subprocess import PIPE, Popen
# from urllib import request as urlRequest
# import webbrowser
#try:
#    import pyvista as pv
#    pvfound = True
#except:
#    pvfound = False
#    print('WARNING: pyvista not found, 3D plotting capabilities will be limited.')
# from resipy.Project import R2
# from resipy.r2help import r2help



# debug options
DEBUG = True # set to false to not display message in the console
def pdebug(*args, **kwargs):
    if DEBUG:
        print('DEBUG:', *args, **kwargs)
    else:
        pass

#%% General crash ERROR

def errorMessage(etype, value, tb):
    print('ERROR begin:')
    traceback.print_exception(etype, value, tb)
    print('ERROR end.')
    errorMsg = traceback.format_exception(etype, value, tb,limit=None, chain=True)
    finalError ='Build: %s\n' % ResIPy_version # adding ResIPy version and system info to the error
    for key in sysinfo: # adding system info to the error
        finalError += str(key) + ': ' + str(sysinfo[key]) + '\n'
    for errs in errorMsg:
        finalError += errs
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText("<b>Critical error:</b>")
    msg.setInformativeText('''Please see the detailed error below.<br>You can report the errors at:<p><a href='https://gitlab.com/hkex/resipy/issues'>https://gitlab.com/hkex/resipy/issues</a></p><br>''')
    msg.setWindowTitle("Error!")
    msg.setDetailedText('%s' % (finalError))
    msg.setStandardButtons(QMessageBox.Retry)
    msg.exec_()

def catchErrors():
    sys.excepthook = errorMessage
    threading.Thread.__init__
    
    
#%%get relative path of images
def resource_path(relative_path):
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)


class customThread(QThread):
    ''' class needed to run out of main thread computation and then emit
    signal to the main GUI thread to display message box.
    '''
    signal = pyqtSignal('PyQt_PyObject')

    def __init__(self, func):
        QThread.__init__(self)
        self.func = func

    def run(self):
        output = self.func()
        self.signal.emit(output) # inform the main thread of the output


# worker class to run the inversion in another thread
# https://realpython.com/python-pyqt-qthread/
class Worker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal('PyQt_PyObject')
    
    def __init__(self, project, kwargs, pseudo=False):
        QThread.__init__(self)
        self.project = project
        self.kwargs = kwargs
        self.pseudo = pseudo

    def run(self):
        """Long-running task."""
        def emitLog(x):
            self.progress.emit(x)
        self.kwargs['dump'] = emitLog
        if self.pseudo == False:
            self.project.invert(**self.kwargs)
        else:
            self.project.invertPseudo3D(**self.kwargs)
        self.finished.emit()
        
# small code to see where are all the directories
frozen = 'not'
if getattr(sys, 'frozen', False):
        # we are running in a bundle
        frozen = 'ever so'
        bundle_dir = sys._MEIPASS
else:
        # we are running in a normal Python environment
        bundle_dir = os.path.dirname(os.path.abspath(__file__))
pdebug( 'we are',frozen,'frozen')
pdebug( 'bundle dir is', bundle_dir )
pdebug( 'sys.argv[0] is', sys.argv[0] )
pdebug( 'sys.executable is', sys.executable )
pdebug( 'os.getcwd is', os.getcwd() )


#%% MatplotlibWidget class
class MatplotlibWidget(QWidget):
    def __init__(self, parent=None, figure=None, navi=False, itight=False,
                 threed=False, aspect='auto'):
        super(MatplotlibWidget, self).__init__(parent) # we can pass a figure but we can replot on it when
        # pushing on a button (I didn't find a way to do it) while with the axes, you can still clear it and
        # plot again on them
        self.callback = None
        clearIt = False
        if figure is None:
            clearIt = True
            figure = Figure(tight_layout=itight) # tight_layout will be called internally for each draw
            self.canvas = FigureCanvasQTAgg(figure)
            if threed is True:
                axes = figure.add_subplot(111, projection='3d')
            else:
                axes = figure.add_subplot(111)
        else:
            axes = figure.get_axes()
        self.figure = figure
        if resipySettings.param['dark'] == 'True':
            self.figure.set_facecolor((0.2,0.2,0.2)) # just making figures's face color grey not the whole plot 
        self.axis = axes
        self.xlim = None
        self.ylim = None
        self.xlim0 = None
        self.ylim0 = None
        self.aspect = aspect
        self.layoutVertical = QVBoxLayout(self)
        self.layoutVertical.addWidget(self.canvas)#, stretch = 1, alignment=Qt.AlignCenter)
        if clearIt is True:
            self.clear()

        if navi is True:
            self.navi_toolbar = NavigationToolbar(self.canvas, self)
            self.navi_toolbar.setMaximumHeight(30)
            self.navi_toolbar.actions()[0].triggered.connect(self.getHome)
            self.layoutVertical.addWidget(self.navi_toolbar)        


    def setMinMax(self, vmin=None, vmax=None):
        coll = self.axis.collections[0]
        if vmin is None:
            vmin = np.nanmin(coll.get_array())
        else:
            vmin = float(vmin)
        if vmax is None:
            vmax = np.nanmax(coll.get_array())
        else:
            vmax = float(vmax)
        coll.set_clim(vmin, vmax)
        self.canvas.draw()


    def plot(self, callback, aspect = None, threed=False):
        ''' call a callback plot function and give it the ax to plot to
        '''
        self.figure.clear() # need to clear the figure with the colorbar as well
        if threed is False:
            ax = self.figure.add_subplot(111)
        else:
            ax = self.figure.add_subplot(111, projection='3d')            
        self.callback = callback
        callback(ax=ax)
        if aspect == None:
            aspect = self.aspect
        ax.set_aspect(aspect)
        self.canvas.draw()

    def setCallback(self, callback):
        self.callback = callback
        self.xlim0 = None
        self.ylim0 = None
        self.xlim = None
        self.ylim = None
    
    ''' to keep the zoom level we register event on xylim change and store it
    as attribute of the MatplotlibWidget. At first drawn of the figure by
    replot() we also set up xlim0 and ylim0 to be triggered by the 'home'
    button. If another replot is triggered the plot will programmatically
    restore the previous zoom level using getHome()
    '''

    def on_xlims_change(self, ax):
        self.xlim = ax.get_xlim()
        
    def on_ylims_change(self, ax):
        self.ylim = ax.get_ylim()

    def setHome(self, ax):
        xlim = ax.get_xlim()
        if (xlim[0] != 0) | (xlim[1] != 1): # it means callback actually plot smth
            self.xlim0 = ax.get_xlim()
            self.ylim0 = ax.get_ylim()

    def getHome(self):
        self.axis.set_xlim(self.xlim0)
        self.axis.set_ylim(self.ylim0)
        self.canvas.draw()
        
    def restoreZoom(self):
        self.figure.axes[0].set_xlim(self.xlim)
        self.figure.axes[0].set_ylim(self.ylim)

    def replot(self, threed=False, aspect=None, **kwargs):
        self.figure.clear()
        if threed is False:
            ax = self.figure.add_subplot(111)
        else:
            ax = self.figure.add_subplot(111, projection='3d')
        self.axis = ax
        self.callback(ax=ax, **kwargs)
        if self.xlim0 is None:
            self.setHome(ax)
        ax.callbacks.connect('xlim_changed', self.on_xlims_change)
        ax.callbacks.connect('ylim_changed', self.on_ylims_change)
        if aspect == None:
            aspect = self.aspect
        ax.set_aspect(aspect)
        if self.xlim0 is not None:
            self.restoreZoom()
        self.canvas.draw()

    def clear(self):
        self.axis.clear()
        self.figure.clear()
        self.canvas.draw()


#%% Main class
class App(QMainWindow):
    def __init__(self, parent=None):
        super().__init__()
        # do the checks for wine and updates in seperate thread
        tupdate = customThread(self.updateChecker)
        tupdate.signal.connect(self.updateCheckerShow)
        tupdate.start()
        
        twine = customThread(self.checkWine)
        twine.signal.connect(self.checkWineShow)
        twine.start()
        
        self.setWindowTitle('ResIPy v' + ResIPy_version)
        self.setGeometry(100,100,1100,600)
        self.newwd = os.path.join(bundle_dir, 'resipy')

        # UI attributes (sync with self.project attributes when instantiated)
        self.project = None
        self.loadedProjectFlag = False # to prevent removing loaded files from working directory
        self.typ = 'R2'
        self.parser = None
        self.fname = None
        self.iBatch = False
        self.iTimeLapse = False
        self.iBorehole = False
        self.iForward = False
        self.inputPhaseFlag = False
        self.iCropping = True # by default crop the mesh
        self.num_xz_poly = None # to store the values
        self.tempElec = None # place holder to compare the new electrode agains
        self.clip = None # to store the clipped mesh for 3D forward
        self.apiLog = '' # store API calls
        self.writeLog('# -------- ResIPy log file -------\n')
        
        if frozen == 'not':
            self.datadir = os.path.join(bundle_dir, './examples')
        else:
            self.datadir = os.path.join(bundle_dir, 'resipy', 'examples')
        pdebug('self.datadir=', self.datadir)
        self.plotAspect = 'equal'
        self.iDesign = False # boolean to know if the mesh has been designed before meshing

        self.tableWidget = QWidget()
        self.layout = QVBoxLayout()
        self.tabs = QTabWidget()
        
        
        #%% message bar below the UI for info and error message
        self.timer = QTimer() # a better method to get rid of expired messages in status bar
        if resipySettings.param['dark'] == 'True':
            initErrMsgColor = 'white'
        else:
            initErrMsgColor = 'black'
        self.errorLabel = QLabel('<i style="color:%s">Error messages will be displayed here</i>' % initErrMsgColor)
        self.errorLabel.setWordWrap(True)
        QApplication.processEvents()



        #%% tab 1 importing data
        self.tabImporting = QTabWidget()
        self.tabs.addTab(self.tabImporting, 'Importing')

        self.tabImportingData = QWidget()
        self.tabImporting.addTab(self.tabImportingData, 'Data')
        self.tabImportingDataLayout = QVBoxLayout()
        
        # restart a new project
        # self.restartBtn = QPushButton('Restart')
        # self.restartBtn.setAutoDefault(True)
        # self.restartBtn.clicked.connect(self.restartFunc)
        # self.restartBtn.setToolTip('Press to reset all tabs and start a new survey.')
        
        # saving project
        def saveProjectBtnFunc():
            fname, _ = QFileDialog.getSaveFileName(self.tabImportingData,'Save Project',
                                                   self.datadir, '.resipy')
            if fname != '':
                self.loadingWidget('Saving project...')
                self.project.saveProject(fname)
                self.loadingWidget(exitflag=True)
                self.infoDump('Project successfully saved.')
        # self.saveProjectBtn = QPushButton('Save Project')
        # self.saveProjectBtn.setToolTip('Save files and results to load them again in ResIPy later.')
        # self.saveProjectBtn.clicked.connect(saveProjectBtnFunc)
        
        # loading project
        def loadBatchProject(): # for batch/time-lapse/pseudo 3D
            self.fnamesCombo.clear()
            self.psContourCheck.setEnabled(True)
            self.mergElecCheck.setEnabled(True)
            for s in self.project.surveys: # we already loaded the results
                self.fnamesCombo.addItem(s.name)
                self.errFitPlotIndexList.append(0)
                self.iperrFitPlotIndexList.append(0)
            self.errorCombosShow(True)
            errorCombosFill(self.prepFnamesComboboxes)
            if self.pseudo3DCheck.isChecked():
                self.fnamesCombo.hide()
                self.fnamesComboLabel.hide()
            else:
                self.fnamesCombo.show()
                self.fnamesComboLabel.show()

        
        def loadProjectBtnFunc():
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingData,'Open File', self.datadir, '*.resipy')
            if fname != '':
                try:
                    self.loadingWidget('Loading Project, please wait...', False)
                    k = Project(dirname=self.newwd)
                    k.loadProject(fname)
                    self.loadedProjectFlag = True
                    
                    # set flags that call restartFunc()
                    if k.iTimeLapse:
                        self.timeLapseCheck.setChecked(True)
                    if k.iBatch:
                        self.batchCheck.setChecked(True)
                    if k.pseudo3DSurvey is not None:
                        self.pseudo3DCheck.setChecked(True)
                    if (k.typ == 'R2') | (k.typ == 'cR2'):
                        self.m2DRadio.setChecked(True)
                    else:
                        self.m3DRadio.setChecked(True)
                    if k.iForward:
                        self.fwdRadio.setChecked(True)
                    else:
                        self.invRadio.setChecked(True)
                    
                    # pre-processing default flags
                    self.recipErrDataIndex = -1
                    self.errFitApplyToAll = True
                    self.errFitDataIndex = -1
                    self.iperrFitDataIndex = -1
                    self.errFitType.setCurrentIndex(0)
                    self.iperrFitType.setCurrentIndex(0)
                    self.dcaProgress.setValue(0)
                    self.errFitPlotIndexList = []
                    self.iperrFitPlotIndexList = []
                    self.phaseFiltDataIndex = -1
                    
                    self.project = k # assign the project object
                    self.project.darkMode = eval(resipySettings.param['dark']) # set matplotlib theme
                    if self.pseudo3DCheck.isChecked():
                        self.elecTable.initTable(self.project.pseudo3DSurvey.elec) # load electrodes back in
                    else:
                        self.elecTable.initTable(self.project.elec) # load electrodes back in
                    self.tempElec = self.elecTable.getTable()
                    if self.project.iBorehole:
                        self.boreholeCheck.setChecked(True)
                    
                    if not self.project.topo.empty:
                        if self.project.topo.shape[0] > 0:
                            self.topoTable.initTable(self.project.topo.values)
                    else:
                        self.topoTable.initTable(np.array([['',''],['','']]))
                    
                    # display pseudo-section and filtering graphs
                    self.settingUI()
                    
                    # load batch/time-lapse/pseudo 3D results 
                    if self.batchCheck.isChecked() or self.timeLapseCheck.isChecked():
                        loadBatchProject()
                        
                    # load IP
                    if (self.project.typ == 'cR2') | (self.project.typ == 'cR3t'):
                        self.ipCheck.setChecked(True)
                        self.ipCheck.setEnabled(True)
                    else:
                        self.ipCheck.setChecked(False)
                        self.ipCheck.setEnabled(False)
                        
                    # display mesh
                    if self.project.mesh is not None:
                        if (self.project.typ == 'R2') | (self.project.typ == 'cR2'):
                            replotMesh()
                        else:
                            if pvfound:
                                self.mesh3Dplotter.clear() # clear all actors 
                                if len(np.unique(self.project.mesh.df['region'])) == 1:
                                    color_map = 'Greys'
                                    color_bar = False
                                else:
                                    color_map = 'Spectral'
                                    color_bar = True
                                self.project.showMesh(ax=self.mesh3Dplotter, color_map=color_map, color_bar=color_bar)
                            else:
                                self.mwMesh3D.plot(self.project.showMesh, threed=True)
                            self.meshOutputStack.setCurrentIndex(2)
                    
                    # load mesh regions back in
                    if self.project.mesh is not None:
                        self.regionTable.nrow = 0
                        self.regionTable.setRowCount(0)
                        for i in range(len(np.unique(self.project.mesh.df['region']))):
                            self.regionTable.addRow()
                            ie = self.project.mesh.df['region'] == i+1
                            row = self.project.mesh.df[ie]
                            self.regionTable.setItem(i, 0, QTableWidgetItem(str(row['res0'].values[0])))
                            self.regionTable.setItem(i, 1, QTableWidgetItem(str(row['phase0'].values[0])))
                            self.regionTable.setItem(i, 2, QTableWidgetItem(str(row['zones'].values[0])))
                            if row['param'].values[0] == 0:
                                self.regionTable.cellWidget(i,3).findChildren(QCheckBox)[0].setChecked(True)
                    
                    # restore changed inversion settings
    #                 self.modelDOICheck
    #                 self.checkTxSign
                    if 'flux_type' in self.project.param:
                        if self.project.param['flux_type'] == 3:
                            self.flux_type.setCurrentIndex(0)
                        else:
                            self.flux_type.setCurrentIndex(1)
                    if 'singular_type' in self.project.param:
                        if self.project.param['singular_type'] == 1:
                            self.singular_type.setChecked(True)
                        else:
                            self.singular_type.setChecked(False)
                    if 'res_matrix' in self.project.param:
                        self.res_matrix.setCurrentIndex(self.project.param['res_matrix'])
                    if 'scale' in self.project.param:
                        self.scale.setText(str(self.project.param['scale']))
                    if 'patch_x' in self.project.param:
                        self.patch_x.setText(str(self.project.param['patch_x']))
                    if 'patch_z' in self.project.param:
                        self.patch_z.setText(str(self.project.param['patch_z']))
                    if 'inv_type' in self.project.param:
                        self.inv_type.setCurrentIndex(self.project.param['inv_type'])
                    if 'target_decrease' in self.project.param:
                        self.target_decrease.setText(str(self.project.param['target_decrease']))
                    if 'data_type' in self.project.param:
                        self.data_type.setCurrentIndex(self.project.param['data_type'])
                    if 'reg_mode' in self.project.param:
                        self.reg_mode.setCurrentIndex(self.project.param['reg_mode'])
                    if 'tolerance' in self.project.param:
                        self.tolerance.setText(str(self.project.param['tolerance']))
                    if 'max_iter' in self.project.param:
                        self.max_iterations.setText(str(self.project.param['max_iter']))
                    if 'error_mod' in self.project.param:
                        self.error_mod.setCurrentIndex(self.project.param['error_mod'])
                    if 'alpha_aniso' in self.project.param:
                        self.alpha_aniso.setText(str(self.project.param['alpha_aniso']))
                    if 'alpha_s' in self.project.param:
                        self.alpha_s.setText(str(self.project.param['alpha_s']))
                    if 'min_error' in self.project.param:
                        self.min_error.setText(str(self.project.param['min_error']))
                    if 'a_wgt' in self.project.param:
                        self.a_wgt.setText(str(self.project.param['a_wgt']))
                    if 'b_wgt' in self.project.param:
                        self.b_wgt.setText(str(self.project.param['b_wgt']))
                    if 'rho_min' in self.project.param:
                        self.rho_min.setText(str(self.project.param['rho_min']))
                    if 'rho_max' in self.project.param:
                        self.rho_max.setText(str(self.project.param['rho_max']))
                    
                    # display inversion results
                    if len(self.project.meshResults) > 0:
                        if self.fwdRadio.isChecked() and len(self.project.meshResults) == 1: # a fwd-only project is saved
                            pass
                        else:
                            self.displayInvertedResults()
                            prepareInvError()
                            self.logText.setText(self.project.invLog)
                    self.loadingWidget(exitflag=True)
                    self.loadedProjectFlag = False # now if user toggles stuff, working directory should be cleaned.
                    # activate tabs
                    self.activateTabs(True)
                    if self.project.iForward:
                        self.tabs.setTabEnabled(4, True)
                        self.tabs.setTabEnabled(5, True)
                    self.tabs.setTabEnabled(6, True)
                    
                    # write log and let the world knows!
                    self.writeLog('k.loadProject("{:s}")'.format(fname))
                    self.infoDump('Project successfully loaded.')
                except Exception as e:
                    self.loadingWidget(exitflag=True)
                    self.errorDump('Error in loading project: ' + str(e))

        
        # instead, let's make a hamburger menu in the top right corner
        self.optionMenu = QMenu()
        self.optionMenu.addAction('Load Project', loadProjectBtnFunc)
        self.optionMenu.addAction('Save Project', saveProjectBtnFunc)
        self.optionMenu.addAction('Restart Project', self.restartFunc)
        self.optionMenu.addAction('Save API log', self.saveLog)
        themeMode = 'Light theme' if resipySettings.param['dark'] == 'True' else 'Dark theme'
        self.optionMenu.addAction(themeMode, self.darkModeFunc)
        self.optionMenu.addAction('Restart ResIPy', self.restartGUI)
        self.hamBtn = QPushButton('Options')
        self.tabs.setCornerWidget(self.hamBtn, Qt.TopRightCorner)        
        self.hamBtn.setMenu(self.optionMenu)
        
        
        def dimSurvey():
            if self.m2DRadio.isChecked():
                self.typ = self.typ.replace('3t','2')
                if self.project is not None:
                    self.project.typ = self.project.typ.replace('3t','2')
                    self.writeLog('k.typ = {:s}'.format(self.project.typ))
                # importing tab
                self.elecTable.setColumnHidden(2, True)
                self.topoTable.setColumnHidden(2, True)
                self.regular3DCheck.setChecked(False)
                self.pseudo3DCheck.setVisible(True)
                if pvfound:
                    self.pseudo3DCheck.setEnabled(True)
                else:
                    self.pseudo3DCheck.setEnabled(False)
                if self.project is not None:
                    if self.project.elec is not None:
                        self.elecTable.initTable(self.project.elec)
                self.elecLineEdit.setEnabled(False)
                self.elecLineSpacingEdit.setEnabled(False)
                self.elecLineEdit.setText('1')
                self.boreholeCheck.setChecked(False)
                self.boreholeCheck.setEnabled(True)
                self.regular3DCheck.setVisible(False)
                self.pParams['threed'] = False
                self.pParamsIP['threed'] = False
                self.pseudoLayout.setCurrentIndex(0)
                # self.ipCheck.setEnabled(True)
                
                #Pre-processing tab
                self.recipErrorBottomTabs.setTabEnabled(0, True)
                self.recipErrorBottomTabs.setCurrentIndex(0)
                self.recipErrorSaveBtn.setVisible(True)
                self.tabPreProcessing.setCurrentIndex(0)
                self.tabPreProcessing.setTabEnabled(0, True)
                self.mwManualFiltering.show()
                try:
                    if not self.project.surveys[0].df.empty:
                        self.tabs.setTabEnabled(1, True)
                except:
                    pass
                
                # mesh tab
                self.meshQuadGroup.setVisible(True)
                self.meshTrianGroup.setVisible(True)
                self.mesh3DGroup.setVisible(False)
                self.meshTetraGroup.setVisible(False)
                self.meshTankGroup.setVisible(False)
                self.meshCylinderGroup.setVisible(False)
                self.meshCustom2dGroup.setVisible(True)
                self.meshCustom3dGroup.setVisible(False)
                self.instructionLabel.setVisible(True)
                self.meshAspectBtn.setVisible(True)
                self.resetMeshBtn.setVisible(True)
                self.instructionLabel3D.setVisible(False)
                self.select3DRegionBtn.setVisible(False)
                self.add3DRegionBtn.setVisible(False)
                self.fin3DRegionBtn.setVisible(False)
                
                # inversion settings
                show3DOptions(False)
                if self.project is not None:
                    showIpOptions(self.typ[0] == 'c')

                # inversion tab
                self.contourCheck.setVisible(True)
                self.clipCornersCheck.setVisible(True)
                self.doiSensCheck.setVisible(True)
                self.sensWidget.setVisible(True)
                show3DInvOptions(False)
                
                # post-processing
                self.errorGraphs.setTabEnabled(0, True)
                self.errorGraphs.setCurrentIndex(0)
            else:
                self.typ = self.typ.replace('2','3t')
                if self.project is not None:
                    self.project.typ = self.project.typ.replace('2', '3t')
                    self.writeLog('k.typ = {:s}'.format(self.project.typ))

                # importing tab
                self.elecLineEdit.setEnabled(True)
                self.elecLineSpacingEdit.setEnabled(True)
                self.elecTable.setColumnHidden(2, False)
                self.topoTable.setColumnHidden(2, False)
                if self.project is not None:
                    if self.project.elec is not None:
                        self.elecTable.initTable(self.project.elec)
                self.invRadio.setChecked(True)
                # self.boreholeCheck.setChecked(True) # to disable pseudo-section
                self.boreholeCheck.setEnabled(False)
                self.regular3DCheck.setVisible(True)
                self.pseudo3DCheck.setChecked(False)
                self.pseudo3DCheck.setVisible(False)
                self.pParams['threed'] = True
                self.pParamsIP['threed'] = True
                if pvfound:
                    self.pseudoLayout.setCurrentIndex(1)
                # self.ipCheck.setEnabled(False) # TODO disabling cR3t for now
                
                #Pre-processing tab
                self.recipErrorBottomTabs.setTabEnabled(0, False)
                self.recipErrorSaveBtn.setVisible(False)
                self.mwManualFiltering.hide() # manual filtering pseudo section is not optimized for 3D
                
                try:
                    if all(self.project.surveys[0].df['irecip'].values == 0) is False:
                        self.recipErrorBottomTabs.setCurrentIndex(1)
                except:
                    pass

                # mesh tab
                self.meshQuadGroup.setVisible(False)
                self.meshTrianGroup.setVisible(False)
                self.mesh3DGroup.setVisible(True)
                self.mesh3DBtn.disconnect()
                self.mesh3DBtn.clicked.connect(meshTetraFunc)
                self.mesh3DCombo.disconnect()
                self.mesh3DCombo.setCurrentIndex(0)
                self.mesh3DCombo.currentIndexChanged.connect(mesh3DComboFunc)
                self.meshTetraGroup.setVisible(True)
                self.meshTankGroup.setVisible(False)
                self.meshCylinderGroup.setVisible(False)
                self.meshCustom3dGroup.setVisible(False)
                self.meshCustom2dGroup.setVisible(False)
                self.instructionLabel.setVisible(False)
                self.meshAspectBtn.setVisible(False)
                self.resetMeshBtn.setVisible(False)
                self.instructionLabel3D.setVisible(True)
                self.select3DRegionBtn.setVisible(True)
                self.add3DRegionBtn.setVisible(True)
                self.fin3DRegionBtn.setVisible(True)

                # inversion settings
                show3DOptions(True)
                if self.project is not None:
                    showIpOptions(self.typ[0] == 'c')

                # inversion tab
                self.contourCheck.setVisible(False)
                self.clipCornersCheck.setVisible(False)
                self.doiSensCheck.setVisible(False)
                self.sensWidget.setVisible(False)
                show3DInvOptions(True)
                
                # post-processing
                self.errorGraphs.setTabEnabled(0, False)
                self.errorGraphs.setCurrentIndex(1)
            try: # to force update the pseudo sections
                if self.project is not None:
                    self.plotPseudo()
                    if self.project.typ[0] == 'c':
                        self.plotPseudoIP()
            except:
                pass

        self.m2DRadio = QRadioButton('2D')
        self.m2DRadio.setChecked(True)
        self.m2DRadio.toggled.connect(dimSurvey)
        self.m3DRadio = QRadioButton('3D')
        self.m3DRadio.setChecked(False)
        self.m3DRadio.toggled.connect(dimSurvey)
        self.dimLayout = QHBoxLayout()
        self.dimLayout.addWidget(self.m2DRadio)
        self.dimLayout.addWidget(self.m3DRadio)
        self.dimLayout.setContentsMargins(0,0,0,0)
        self.dimGroup = QGroupBox()
        self.dimGroup.setLayout(self.dimLayout)
        self.dimGroup.setFlat(True)
        self.dimGroup.setContentsMargins(0,0,0,0)
        self.dimGroup.setStyleSheet('QGroupBox{border: 0px;'
                                'border-style:inset;}')

        # meta data (title and date of survey)
        self.titleLabel = QLabel('Title')
        self.titleEdit = QLineEdit()
        self.titleEdit.setText('My beautiful survey')
        self.titleEdit.setToolTip('This title will be used in the ".in" file.')

        self.dateLabel = QLabel('Date')
        self.dateEdit = QLineEdit()
        self.dateEdit.setText(datetime.now().strftime('%Y-%m-%d')) # get today date
        self.dateEdit.setToolTip('This date will be used in the ".in" file.')

        def timeLapseCheckFunc(state):
            if state == Qt.Checked:
                self.iTimeLapse = True
                if self.project is not None: # if there is already an R2 object
                    self.restartFunc()
                    self.reg_mode.setCurrentIndex(2)
                self.importDataBtn.setText('Import multiple datasets')
                self.importDataBtn.clicked.disconnect()
                self.importDataBtn.clicked.connect(getdir)
                self.batchCheck.setEnabled(False)
                self.pseudo3DCheck.setEnabled(False)
                self.pseudo3DCheck.setChecked(False)
            else:
                self.iTimeLapse = False
                if self.project is not None:
                    self.restartFunc()
                    self.reg_mode.setCurrentIndex(0)
                self.importDataBtn.setText('Import Data')
                self.importDataBtn.clicked.disconnect()
                self.importDataBtn.clicked.connect(importDataBtnFunc)
                self.batchCheck.setEnabled(True)
                self.pseudo3DCheck.setEnabled(True)

        self.timeLapseCheck = QCheckBox('Time-lapse Survey')
        self.timeLapseCheck.stateChanged.connect(timeLapseCheckFunc)
        self.timeLapseCheck.setToolTip('Check to import time-lapse datasets and enable time-lapse inversion.')

        def batchCheckFunc(state):
            if state == Qt.Checked:
                self.iBatch = True
                if self.project is not None:
                    self.restartFunc()
                self.importDataBtn.setText('Import multiple datasets')
                self.importDataBtn.clicked.disconnect()
                self.importDataBtn.clicked.connect(getdir)
                self.timeLapseCheck.setEnabled(False)
            else:
                self.iBatch = False
                if self.project is not None:
                    self.restartFunc()
                self.importDataBtn.setText('Import Data')
                self.importDataBtn.clicked.disconnect()
                self.importDataBtn.clicked.connect(importDataBtnFunc)
                self.timeLapseCheck.setEnabled(True)
        self.batchCheck = QCheckBox('Batch Inversion')
        self.batchCheck.stateChanged.connect(batchCheckFunc)
        self.batchCheck.setToolTip('Check if you want to invert multiple surveys with the same settings and same electrodes.')
        
        def boreholeCheckFunc(state):
            if state == Qt.Checked:
                self.iBorehole = True
                self.errorGraphs.setTabEnabled(0, False)
                if self.project is not None:
                    self.project.setBorehole(True)
                    self.writeLog('k.setBorehole(True)')
            else:
                self.iBorehole = False
                self.errorGraphs.setTabEnabled(0, True)
                if self.project is not None:
                    self.project.setBorehole(False)
                    self.writeLog('k.setBorehole(False)')
            try:
                if self.fname is not None:
                    self.plotPseudo()
                    self.plotPseudoIP()
            except:
                pass
        self.boreholeCheck = QCheckBox('Unconventional Survey')
        self.boreholeCheck.stateChanged.connect(boreholeCheckFunc)
        self.boreholeCheck.setToolTip('Check if you have an unconventional survey (e.g. boreholes).\nThis will just change the pseudo-section.')

        def pseudo3DFunc(state):
            if state == Qt.Checked:
                self.lineSpacing.setVisible(True)
                self.lineSpacingLabel.setVisible(True)
                self.create3DBtn.setVisible(True)
                self.create3DBtn.clicked.disconnect()
                self.create3DBtn.clicked.connect(createPseudo3DFunc)
                self.importDataBtn.setVisible(False)
                self.batchCheck.setChecked(True)
                self.batchCheck.setEnabled(False)
                self.mergElecCheck.setEnabled(False)
                self.timeLapseCheck.setEnabled(False)
                self.fwdRadio.setEnabled(False)
                self.pseudoLayout.setCurrentIndex(1)
                self.elecLineEdit.setEnabled(True)
                self.elecLineSpacingEdit.setEnabled(True)
                self.elecTable.setColumnHidden(2, False)
                self.topoTable.hide() # we can't add extra topo points in here, we don't know which mesh they belong to
                self.regionTable.hide() # we can't add regions at this level - future: set/design/create individual 2D meshes then create a pseudo 3D mesh
                self.topoBtn.hide()
                self.topoAddRowBtn.hide()
                self.elecLabelTextP2 = '<br><b>IMPORTANT:</b> Labels <b>Must</b> be defined in <b><font color="red">"Line Number [space] Electrode Number"</font></b> format.'
                self.designModelBtn.setEnabled(False)
                self.importCustomMeshBtn2.setEnabled(False)
                self.resetMeshBtn.hide()
                self.instructionLabel.hide()
                self.meshAspectBtn.hide()
                self.meshSubLayout.setStretch(0,100)
                self.meshSubLayout.setStretch(1,0)
                self.meshOutputStack.setCurrentIndex(2)
            
            else:
                self.lineSpacing.setVisible(False)
                self.lineSpacingLabel.setVisible(False)
                self.create3DBtn.setVisible(False)
                self.importDataBtn.setVisible(True)
                self.batchCheck.setChecked(False)
                self.batchCheck.setEnabled(True)
                self.mergElecCheck.setEnabled(True)
                self.timeLapseCheck.setEnabled(True)
                self.fwdRadio.setEnabled(True)
                self.pseudoLayout.setCurrentIndex(0)
                self.elecLineEdit.setEnabled(False)
                self.elecLineSpacingEdit.setEnabled(False)
                self.elecTable.setColumnHidden(2, True)
                self.topoTable.show()
                self.regionTable.show()
                self.topoBtn.show()
                self.topoAddRowBtn.show()
                self.elecLabelTextP2 = ''
                self.designModelBtn.setEnabled(True)
                self.importCustomMeshBtn2.setEnabled(True)
                self.resetMeshBtn.show()
                self.instructionLabel.show()
                self.meshAspectBtn.show()
                self.meshSubLayout.setStretch(0,70)
                self.meshSubLayout.setStretch(1,30)
                self.meshOutputStack.setCurrentIndex(1)
            self.elecLabel.setText(self.elecLabelTextP1 + self.elecLabelTextP2)
                
        self.pseudo3DCheck = QCheckBox('Pseudo 3D inversion from 2D lines')
        self.pseudo3DCheck.setToolTip('Arrange 2D inversion lines in a 3D grid. ResIPy will assume batch 2D surveys are imported.\n'
                                      'Electrodes can be manually set in "Electrodes (XYZ/Topo)" tab.\n'
                                      '***beta*** some ResIPy features may not work properly.\n\n'
                                      'IMPORTANT: electrode labels MUST be imported and have "<line number> <electrode number>" format.')
        self.pseudo3DCheck.stateChanged.connect(pseudo3DFunc)
        if pvfound is False: # We can't (don't want to) have this feature without pyvista!!
            self.pseudo3DCheck.setEnabled(False)
        
        def regular3DFunc(state):
            if state == Qt.Checked:
                self.lineSpacing.setVisible(True)
                self.lineSpacingLabel.setVisible(True)
                self.importDataBtn.setVisible(False)
                self.create3DBtn.setVisible(True)
                self.create3DBtn.clicked.disconnect()
                self.create3DBtn.clicked.connect(create3DFunc) 
            else:
                self.lineSpacing.setVisible(False)
                self.lineSpacingLabel.setVisible(False)
                self.create3DBtn.setVisible(False)
                self.importDataBtn.setVisible(True)
                
        self.regular3DCheck = QCheckBox('3D survey from 2D lines')
        self.regular3DCheck.setToolTip('3D survey from 2D lines - electrodes can be manually set in "Electrodes (XYZ/Topo)" tab')
        self.regular3DCheck.stateChanged.connect(regular3DFunc)
        self.regular3DCheck.setVisible(False)
        
        # select inverse or forward model
        def fwdRadioFunc():
            self.iForward = True
            self.ftypeCombo.setEnabled(False)
            self.invRadio.setChecked(False)
            self.tabs.setTabEnabled(1,False)
            self.tabs.setTabEnabled(3, True)
            self.importDataBtn.setEnabled(False)
            self.timeLapseCheck.setEnabled(False)
            self.batchCheck.setEnabled(False)
            self.pseudo3DCheck.setEnabled(False)
            self.pseudo3DCheck.setChecked(False)
            self.tabImporting.setTabEnabled(2,False) # no custom parser needed
            if self.loadedProjectFlag is False: # loading a project doesn't need a restart
                self.restartFunc() # let's first from previous inversion
            self.nbElecEdit.setEnabled(True)
            self.tabImporting.setTabEnabled(1, True) # here because restartFunc() set it to False
            self.ipCheck.setEnabled(True)
            self.psContourCheck.setEnabled(False)
            self.mergElecCheck.setEnabled(False)
            self.activateTabs(True)
        self.fwdRadio = QRadioButton('Forward')
        self.fwdRadio.setChecked(False)
        self.fwdRadio.toggled.connect(fwdRadioFunc)
        self.fwdRadio.setToolTip('To create a model, a sequence and see what output you can obtain.')

        def invRadioFunc():
            self.iForward = False
            self.ftypeCombo.setEnabled(True)
            self.fwdRadio.setChecked(False)
            self.spacingEdit.setReadOnly(False)
            self.tabs.setTabEnabled(1,True)
            self.tabs.setTabEnabled(3, False)
            self.tabImporting.setTabEnabled(1, False)
            self.importDataBtn.setEnabled(True)
            self.nbElecEdit.setEnabled(False)
            self.timeLapseCheck.setEnabled(True)
            self.ipCheck.setEnabled(False)
            self.tabImporting.setTabEnabled(2,True)
            self.batchCheck.setEnabled(True)
            self.pseudo3DCheck.setEnabled(True)
            self.timeLapseCheck.setChecked(False) # not checked by default
            self.batchCheck.setChecked(False) # not checked by default
            self.activateTabs(False)
        self.invRadio = QRadioButton('Inverse')
        self.invRadio.setChecked(True)
        self.invRadio.toggled.connect(invRadioFunc)
        self.invRadio.setToolTip('To invert data that is already collected.')
        
        self.dimInvLayout = QHBoxLayout()
        self.dimInvLayout.addWidget(self.fwdRadio)
        self.dimInvLayout.addWidget(self.invRadio)
        self.dimInvLayout.setContentsMargins(0,0,0,0)
        self.dimInvGroup = QGroupBox()
        self.dimInvGroup.setLayout(self.dimInvLayout)
        self.dimInvGroup.setFlat(True)
        self.dimInvGroup.setContentsMargins(0,0,0,0)
        self.dimInvGroup.setStyleSheet('QGroupBox{border: 0px;'
                                'border-style:inset;}')


        # ask for working directory, and survey file to input
        def getwd():
            fdir = QFileDialog.getExistingDirectory(self.tabImportingData, 'Choose Working Directory')
            if fdir != '':
                self.newwd = fdir
                if self.project is not None:
                    self.project.setwd(fdir)
                    self.writeLog('k.setwd("{:s}")'.format(fdir))
                print('Working directory = ', fdir)
                self.wdBtn.setText(os.path.basename(self.newwd))
        self.wdBtn = QPushButton('Working directory:' + os.path.basename(self.newwd))
        self.wdBtn.setAutoDefault(True)
        self.wdBtn.clicked.connect(getwd)
        self.wdBtn.setToolTip('Select the working directory, containing your data\nThe working directory will automatically have all the necessary files for the inversion (e.g. R2.in, R2.exe, protocol.dat, f001_res.vtk, etc.)')

        self.ftype = 'ProtocolDC' # by default
        self.fformat = 'DAT (Tab delimited) (*.dat)' # default

        def ftypeComboFunc(index):
            if index == 0:
                self.ftype = 'ProtocolDC'
                self.fformat = 'DAT (Tab delimited) (*.dat *.DAT)'
            elif index == 1:
                self.ftype = 'Syscal'
                self.fformat = 'Comma Separated Values (*.csv *.CSV)'            
            elif index == 2:
                self.ftype = 'ProtocolIP'
                self.fformat = 'DAT (Tab delimited) (*.dat *.DAT)'
            elif index == 3:
                self.ftype = 'ResInv'
                self.fformat = 'DAT (*.dat *.DAT)'
            elif index == 4:
                self.ftype = 'BGS Prime'
                self.fformat = 'DAT (*.dat *.DAT *.tab *.TAB)'
            elif index == 5:
                self.ftype = 'Sting'
                self.fformat = 'Sting (*.stg *.STG)'
            elif index == 6:
                self.ftype = 'ABEM-Lund'
                self.fformat = 'OHM (*.OHM *.ohm)'
            elif index == 7:
                self.ftype = 'Lippmann'
                self.fformat = 'TX0 (*.tx0 *.TX0);;Text (*.txt *.TXT)'
            elif index == 8:
                self.ftype = 'ARES'
                self.fformat = 'ARES (*.2dm *.2DM)'
            elif index == 9:
                self.ftype = 'BERT'
                self.fformat = 'BERT (*.dat *.DAT *.ohm *.OHM)'
            elif index == 10:
                self.ftype = 'E4D'
                self.fformat = 'srv (*.srv *.SRV)'
            elif index == 11:
                self.ftype = 'DAS-1'
                self.fformat = 'Data (*.dat *.data *.DAT *.DATA *.txt)'
            elif index == 12:
                self.ftype = 'Custom'
                self.tabImporting.setCurrentIndex(2) # switch to the custom parser
            else:
                self.ftype = '' # let to be guessed
        self.ftypeComboLabel = QLabel('File format:')
        self.ftypeComboLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.ftypeCombo = QComboBox()
        self.ftypeCombo.addItem('Protocol DC')
        self.ftypeCombo.addItem('Syscal')
        self.ftypeCombo.addItem('Protocol IP')
        self.ftypeCombo.addItem('ResInv (2D/3D)')
        self.ftypeCombo.addItem('BGS Prime')
        self.ftypeCombo.addItem('Sting')
        self.ftypeCombo.addItem('ABEM-Lund')
        self.ftypeCombo.addItem('Lippmann')
        self.ftypeCombo.addItem('ARES (beta)')
        self.ftypeCombo.addItem('BERT')
        self.ftypeCombo.addItem('E4D')
        self.ftypeCombo.addItem('DAS-1')
        self.ftypeCombo.addItem('Custom')
        self.ftypeCombo.activated.connect(ftypeComboFunc)
        self.ftypeCombo.setFixedWidth(150)
        self.ftypeCombo.setToolTip('Select data format.')

        self.spacingEdit = QLineEdit()
        self.spacingEdit.setValidator(QDoubleValidator())
        self.spacingEdit.setText('-1.0') # -1 let it search for the spacing
        self.spacingEdit.setFixedWidth(80)
        self.spacingEdit.setToolTip('Electrode spacing.')

        def getdir():
            fnames, _ = QFileDialog.getOpenFileNames(self.tabImportingData, 'Select file(s)', self.datadir, self.fformat)
            # fdir = QFileDialog.getExistingDirectory(self.tabImportingData, 'Choose the directory containing the data', directory=self.datadir)
            
            if fnames != []:
                fdir = os.path.dirname(fnames[0])
                self.restartFunc()
                self.datadir = os.path.dirname(fdir)
                try:
                    self.loadingWidget('Loading data, please wait...', False)
                    if self.project.iBatch is False:
                        if len(fnames) < 2:
                            self.errorDump('at least two files needed for timelapse.')
                            return
                        self.project.createTimeLapseSurvey(fnames, ftype=self.ftype, dump=self.infoDump)
                        self.writeLog('k.createTimeLapseSurvey({:s}, ftype="{:s}")'.format(str(fnames), self.ftype))
                        self.ipCheck.setEnabled(False) # TODO enable IP for timelapse
                        self.infoDump('Time-lapse survey created.')
                    else:
                        self.project.createBatchSurvey(fnames, ftype=self.ftype, dump=self.infoDump)
                        self.writeLog('k.createBatchSurvey({:s}, ftype="{:s}")'.format(str(fnames), self.ftype))
                        self.ipCheck.setEnabled(True)
                        self.infoDump('Batch survey created.')
                    self.fnamesCombo.clear()
                    self.psContourCheck.setEnabled(True)
                    self.mergElecCheck.setEnabled(True)
                    
                    for s in self.project.surveys:
                        self.fnamesCombo.addItem(s.name)
                        self.errFitPlotIndexList.append(0)
                        self.iperrFitPlotIndexList.append(0)
                    self.errorCombosShow(True)
                    errorCombosFill(self.prepFnamesComboboxes)
                    self.fnamesCombo.show()
                    self.fnamesComboLabel.show()
                    self.importDataBtn.setText(os.path.basename(fdir) + ' (Press to change)')
                    self.elecTable.initTable(self.project.elec)
                    self.tabImporting.setTabEnabled(1,True)
                    self.invNowBtn.setEnabled(True)
                    self.nbElecEdit.setText(str(self.project.elec.shape[0]))
                    if all(self.project.surveys[0].df['irecip'].values == 0):
                        self.recipOrNoRecipShow(recipPresence=False)
                    else:
                        self.recipOrNoRecipShow(recipPresence=True)
                        self.tabPreProcessing.setTabEnabled(2, True)
                        self.filterAttrCombo.addItem('Reciprocal Error')
                        self.plotError()
                        self.errHist()
                    self.plotManualFiltering()
                    self.plotPseudo()
                    self.loadingWidget(exitflag=True)
                    self.activateTabs(True)
                    if 'dev' in self.project.surveys[0].df.columns:
                        self.filterAttrCombo.addItem('Stacking Error (Dev.)')
                    if 'ip' in self.project.surveys[0].df.columns and self.iTimeLapse is False:
                        if np.sum(self.project.surveys[0].df['ip'].values) > 0 or np.sum(self.project.surveys[0].df['ip'].values) < 0: # np.sum(self.project.surveys[0].df['ip'].values) !=0 will result in error if all the IP values are set to NaN
                            self.ipCheck.setChecked(True)
                        if self.ftype == 'Syscal':
                            self.dcaButton.setEnabled(True)
                            self.dcaProgress.setEnabled(True)
                except Exception as e:
                    self.loadingWidget(exitflag=True)
                    print('Error in getdir(): ', e)
                    self.errorDump('File format is not recognized (one or all files!)')
        
        self.loadingDialog = QDialog()
        self.loadingDialogTxtWidget = QLabel()
        self.loadingDialogTxtWidget.setStyleSheet("font-size: 16px")
        loadingDialogLayout = QHBoxLayout()
        loadingDialogLayout.addWidget(self.loadingDialogTxtWidget)
        self.loadingDialog.setLayout(loadingDialogLayout)
        self.loadingDialog.setWindowFlags(Qt.FramelessWindowHint)
        

        def importDataBtnFunc():
            pdebug('importDataBtnFunc: ftype=', self.ftype)
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingData,'Open File', self.datadir, self.fformat)
            if fname != '':
                self.restartFunc()
                self.datadir = os.path.dirname(fname)
                self.importFile(fname)
                
        self.importDataBtn = QPushButton('Import Data')
        self.importDataBtn.setAutoDefault(True)
        self.importDataBtn.clicked.connect(importDataBtnFunc)
        self.importDataBtn.setToolTip('Select input files (time-lapse: select all the datasets, sorted based on your criteria in the system (e.g., Date modified, Name)).')

        def importDataRecipBtnFunc(): # import reciprocal file
            fnameRecip, _ = QFileDialog.getOpenFileName(self.tabImportingData,'Open File', self.datadir, self.fformat)
            if fnameRecip != '':
                try:
                    self.loadingWidget('Loading data, please wait...', False)
                    self.importDataRecipBtn.setText(os.path.basename(fnameRecip))
                    self.project.addData(fname=fnameRecip, ftype=self.ftype, parser=self.parser)
                    self.writeLog('k.addData(fname={:s}, ftype="{:s}", parser={:s})'.format(
                        fnameRecip, self.ftype, str(self.parser)))
                    if all(self.project.surveys[0].df['irecip'].values == 0) is False:
                        self.recipOrNoRecipShow(recipPresence = True)
                        self.filterAttrCombo.addItem('Reciprocal Error')
                        self.tabPreProcessing.setTabEnabled(2, True) # no point in doing error processing if there is no reciprocal
                        self.plotError()
                        if self.ipCheck.checkState() == Qt.Checked:
                            self.tabPreProcessing.setTabEnabled(3, True)
                            self.recipFiltBtn.setEnabled(True)
                            phaseplotError()
                            heatRaw()
                            heatFilter()    
                        self.errHist()
                        if self.m3DRadio.isChecked():
                            self.recipErrorBottomTabs.setCurrentIndex(1)
                    self.plotPseudo()
                    if self.project.typ[0] == 'c':
                        self.plotPseudoIP()
                    self.plotManualFiltering()
                    self.loadingWidget(exitflag=True)
                    self.infoDump(fnameRecip + ' imported successfully')
                except Exception as e:
                    self.loadingWidget(exitflag=True)
                    pdebug('importFile: ERROR:', e)
                    self.errorDump('Importation failed. File is not being recognized. \
                              Make sure you have selected the right file type.')
                    pass
                    
                
        self.importDataRecipBtn = QPushButton('If you have a reciprocal dataset upload it here')
        self.importDataRecipBtn.setAutoDefault(True)
        self.importDataRecipBtn.clicked.connect(importDataRecipBtnFunc)
        self.importDataRecipBtn.hide()
        self.importDataRecipBtn.setToolTip('Import file with reciprocal measurements (not mandatory).')

        self.lineSpacing = QLineEdit('1')
        self.lineSpacing.setValidator(QDoubleValidator())
        self.lineSpacing.setMaximumWidth(100)
        self.lineSpacing.setFixedWidth(120)
        self.lineSpacing.setVisible(False)
        self.lineSpacingLabel = QLabel('Line spacing [m]:')
        self.lineSpacingLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.lineSpacingLabel.setVisible(False)
        self.lineSpacingLabel.setFixedWidth(120)
        
        def createPseudo3DFunc():
            fnames, _ = QFileDialog.getOpenFileNames(self.tabImportingData, 'Select file(s)', self.datadir, self.fformat)            
            if fnames != []:
                fdir = os.path.dirname(fnames[0])
                self.restartFunc()
                self.datadir = os.path.dirname(fdir)
                val = float(self.lineSpacing.text())
                if len(fnames) <= 1:
                    self.errorDump('More than one dataset is needed!')
                    return
                try:
                    self.loadingWidget('Loading data, please wait...', False)
                    self.project.createPseudo3DSurvey(fnames, lineSpacing=val, ftype=self.ftype, parser=self.parser)
                    self.writeLog('k.createPseudo3DSurvey({:s}, lineSpacing={:s}, ftype="{:s}", parser={:s})'.format(
                        str(fnames), str(val), self.ftype, str(self.parser)))
                    self.infoDump('Pseudo 3D survey from 2D lines created.')
                    self.ipCheck.setEnabled(True)
                    self.create3DBtn.setText(os.path.basename(fdir) + ' (Press to change)')
                    if 'magErr' in self.project.pseudo3DSurvey.df.columns:
                        self.a_wgt.setText('0.0')
                        self.b_wgt.setText('0.0')
                    self.importDataRecipBtn.hide()
                    self.psContourCheck.setEnabled(True)
                    self.mergElecCheck.setEnabled(False)
                    for s in self.project.surveys:
                        self.errFitPlotIndexList.append(0)
                        self.iperrFitPlotIndexList.append(0)
                    self.errorCombosShow(True)
                    errorCombosFill(self.prepFnamesComboboxes)
                    self.elecTable.initTable(self.project.elec)
                    self.tabImporting.setTabEnabled(1,True)
                    if all(self.project.pseudo3DSurvey.df['irecip'].values == 0): # enabling some filtering capacities
                        self.recipOrNoRecipShow(recipPresence=False)
                    else:
                        self.recipOrNoRecipShow(recipPresence=True)
                        self.tabPreProcessing.setTabEnabled(2, True)
                        self.filterAttrCombo.addItem('Reciprocal Error')
                        self.plotError()
                        self.errHist()
                    self.plotManualFiltering()
                    if 'dev' in self.project.pseudo3DSurvey.df.columns:
                        self.filterAttrCombo.addItem('Stacking Error (Dev.)')
                    if 'ip' in self.project.pseudo3DSurvey.df.columns:
                        if np.sum(self.project.pseudo3DSurvey.df['ip'].values) > 0 or np.sum(self.project.surveys[0].df['ip'].values) < 0: # np.sum(self.project.surveys[0].df['ip'].values) !=0 will result in error if all the IP values are set to NaN
                            self.ipCheck.setChecked(True)
                            for proj in self.project.projs:
                                proj.typ = 'cR2'
                        if self.ftype == 'Syscal':
                            self.dcaButton.setEnabled(True)
                            self.dcaProgress.setEnabled(True)               
                    self.nbElecEdit.setText(str(self.project.elec.shape[0]))
                    self.elecDxEdit.setText('{:.2f}'.format(np.diff(self.project.elec[~self.project.elec['remote']]['x'].values[:2])[0]))
                    self.loadingWidget(exitflag=True)
                    self.invNowBtn.setEnabled(True)
                    self.activateTabs(True)
                    self.updateElec()
                    self.plotPseudo()
                    if np.isnan(self.project.elec[['x','y','z']].values).any():
                        self.topoInterpBtnFunc()
                except Exception as e:
                    self.loadingWidget(exitflag=True)
                    print('Error in createPseudo3DFunc(): ', e)
                    self.errorDump('File format is not recognized or not all files have similar format!')
            
        
        def create3DFunc():
            fnames, _ = QFileDialog.getOpenFileNames(self.tabImportingData, 'Select file(s)', self.datadir, self.fformat)            
            if fnames != []:
                fdir = os.path.dirname(fnames[0])
                self.restartFunc()
                self.datadir = os.path.dirname(fdir)
                val = float(self.lineSpacing.text())
                
                try:
                    self.loadingWidget('Loading data, please wait...', False)
                    self.project.create3DSurvey(fnames, lineSpacing=val, ftype=self.ftype, parser=self.parser)
                    self.writeLog('k.create3DSurvey({:s}, lineSpacing={:s}, ftype="{:s}", parser={:s})'.format(
                        str(fnames), str(val), self.ftype, str(self.parser)))
                    self.infoDump('3D survey from regular 2D lines created.')
                    self.ipCheck.setEnabled(True)
                    self.psContourCheck.setEnabled(True)
                    self.mergElecCheck.setEnabled(True)
                    self.create3DBtn.setText(os.path.basename(fdir) + ' (Press to change)')
                    # self.calcAspectRatio()
                    if 'magErr' in self.project.surveys[0].df.columns:
                        self.a_wgt.setText('0.0')
                        self.b_wgt.setText('0.0')
                    self.recipOrNoRecipShow(recipPresence = True)
                    self.importDataRecipBtn.hide()
                    self.tabPreProcessing.setTabEnabled(2, True)
                    self.plotError()
                    self.errHist()
                    self.plotManualFiltering()
                    self.elecTable.initTable(self.project.elec)
                    self.tabImporting.setTabEnabled(1,True)
                    if all(self.project.surveys[0].df['irecip'].values == 0): # enabling some filtering capacities
                        self.recipOrNoRecipShow(recipPresence=False)
                    else:
                        self.recipOrNoRecipShow(recipPresence=True)
                        self.tabPreProcessing.setTabEnabled(2, True)
                        self.filterAttrCombo.addItem('Reciprocal Error')
                        self.plotError()
                        self.errHist()
                    if 'dev' in self.project.surveys[0].df.columns:
                        self.filterAttrCombo.addItem('Stacking Error (Dev.)')
                    if 'ip' in self.project.surveys[0].df.columns:
                        if np.sum(self.project.surveys[0].df['ip'].values) > 0 or np.sum(self.project.surveys[0].df['ip'].values) < 0: # np.sum(self.project.surveys[0].df['ip'].values) !=0 will result in error if all the IP values are set to NaN
                            self.ipCheck.setChecked(True)
                        if self.ftype == 'Syscal':
                            self.dcaButton.setEnabled(True)
                            self.dcaProgress.setEnabled(True)               
                    self.loadingWidget(exitflag=True)
                    self.plotPseudo()
                    self.invNowBtn.setEnabled(True)
                    self.activateTabs(True)
                    self.nbElecEdit.setText(str(self.project.elec.shape[0]))
                    self.elecDxEdit.setText('{:.2f}'.format(np.diff(self.project.elec[~self.project.elec['remote']]['x'].values[:2])[0]))
                    self.fnamesCombo.hide()
                    self.fnamesComboLabel.hide()
                    if np.isnan(self.project.elec[['x','y','z']].values).any():
                        self.topoInterpBtnFunc()
                except Exception as e:
                    self.loadingWidget(exitflag=True)
                    print('Error in create3DFunc(): ', e)
                    self.errorDump('File format is not recognized or not all files have similar format!')
        
        self.create3DBtn = QPushButton('Select 2D lines')
        self.create3DBtn.clicked.connect(create3DFunc)
        self.create3DBtn.setVisible(False)


        def invNowBtnFunc():
            self.tabs.setCurrentIndex(5) # jump to inversion tab
            self.invertBtn.animateClick() # invert
        self.invNowBtn = QPushButton('Invert')
        self.invNowBtn.setStyleSheet('background-color: green; color:black')
        self.invNowBtn.setFixedWidth(150)
        self.invNowBtn.setAutoDefault(True)
        self.invNowBtn.clicked.connect(invNowBtnFunc)
        self.invNowBtn.setEnabled(False)
        self.invNowBtn.setToolTip('Invert with default settings. This will redirect you to the inversion tab.')

        def ipCheckFunc(state):
            if state  == Qt.Checked:
                if self.project.typ == 'R2' or self.project.typ == 'R3t':
                    self.project.typ = 'c' + self.project.typ
                    self.typ = 'c' + self.typ
                else:
                    self.typ = self.project.typ
                showIpOptions(True)
                [p.setVisible(True) for p in [self.pvminIPLabel, self.pvminIP, self.pvmaxIPLabel, self.pvmaxIP]]
#                self.timeLapseCheck.setEnabled(False)
                if self.project.iForward == True:
                    self.mwFwdPseudoIP.setVisible(True)
                    self.noiseLabelIP.show()
                    self.noiseEditIP.show()
                else:
                    self.mwPseudoIP.setVisible(True)
                    self.pseudoFrameIP.setVisible(True)
                    self.plotPseudoIP()
                    self.tabPreProcessing.setTabEnabled(1, True)
                    if all(self.project.surveys[0].df['irecip'].values == 0) is False:
                        phaseplotError()
                        self.tabPreProcessing.setTabEnabled(3, True) # no reciprocity = no IP error model
                        self.recipFiltBtn.setEnabled(True)
                    heatRaw()
    #                self.project.surveys[0].filterDataIP_plot = self.project.surveys[0].filterDataIP_plotOrig
                    self.project.surveys[0].filterDataIP = self.project.surveys[0].df
                    heatFilter()
                self.regionTable.setColumnHidden(1, False)

            else:
                self.project.typ = self.project.typ[1:]
                self.typ = self.typ[1:]
                showIpOptions(False)
                [p.setVisible(False) for p in [self.pvminIPLabel, self.pvminIP, self.pvmaxIPLabel, self.pvmaxIP]]
#                self.timeLapseCheck.setEnabled(True)
                self.mwPseudoIP.setVisible(False)
                self.pseudoFrameIP.setVisible(False)
                self.tabPreProcessing.setTabEnabled(1, False)
                self.tabPreProcessing.setTabEnabled(3, False)
                self.regionTable.setColumnHidden(1, True)
                if self.project.iForward == True:
                    self.mwFwdPseudo.setVisible(True)
                    self.mwFwdPseudoIP.setVisible(False)
                    self.noiseLabelIP.hide()
                    self.noiseEditIP.hide()
            pdebug('ipCheckFunc: mode =', self.project.typ)

        self.ipCheck = QCheckBox('Induced Polarization')
        self.ipCheck.stateChanged.connect(ipCheckFunc)
        self.ipCheck.setEnabled(False)
        self.ipCheck.setToolTip('Check if you have IP data or want IP forward modeling')

        self.fnamesComboLabel = QLabel('Choose a dataset to plot:')
        self.fnamesComboLabel.hide()
        
        def fnamesComboFunc(index):
            self.pParams['index'] = index
            self.pParamsIP['index'] = index
            self.plotPseudo()
            if self.project.typ[0] == 'c':
                self.plotPseudoIP()
                self.dcaProgress.setValue(0)

        self.fnamesCombo = QComboBox()
        self.fnamesCombo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.fnamesCombo.setMinimumWidth(150)
        self.fnamesCombo.activated.connect(fnamesComboFunc)
        self.fnamesCombo.hide()
        
        def psContourFunc(state):
            if state  == Qt.Checked:
                self.pParams['contour'] = True
                self.pParamsIP['contour'] = True
            else:
                self.pParams['contour'] = False
                self.pParamsIP['contour'] = False
            self.plotPseudo()
            if self.project.typ[0] == 'c':
                self.plotPseudoIP()
                
        def mergeElecs(state):
            if self.project is not None and self.project.elec is not None:
                if state  == Qt.Checked:
                    self.elecBackUp = self.project.elec.copy() # backing up electrodes and surveys if checkbox is unchecked
                    self.surveysBackUp = self.project.surveys.copy()
                    if self.project.mergeElec():
                        self.infoDump('Merging close electrodes successful!')
                    else:
                        self.errorDump('Merging close electrodes unsuccessful!')
                else:
                    self.project.elec = self.elecBackUp.copy() # backing up electrodes and surveys if checkbox is unchecked
                    self.project.surveys = self.surveysBackUp.copy()
                
                self.elecTable.initTable(self.project.elec)
                self.tempElec = None
                self.nbElecEdit.setText(str(self.project.elec.shape[0]))
                self.elecDxEdit.setText('{:.2f}'.format(np.diff(self.project.elec[~self.project.elec['remote']]['x'].values[:2])[0]))
                self.updateElec()
        
        # display options for pseudo-sections
        self.mergElecCheck = QCheckBox('Merge close electrodes')
        self.mergElecCheck.stateChanged.connect(mergeElecs)
        self.mergElecCheck.setEnabled(False)
        self.mergElecCheck.setToolTip('Merge electrodes that are very close (0.01 of average electrode spacing)')
        
        self.psContourCheck = QCheckBox('Contour')
        self.psContourCheck.stateChanged.connect(psContourFunc)
        self.psContourCheck.setEnabled(False)
        self.psContourCheck.setToolTip('Check/uncheck to contour pseudo section plots (in 3D mode "Delaunay 3D" is used).')
        
        self.pvminLabel = QLabel('ρ<sub>min</sub>')
        self.pvmin = QLineEdit()
        self.pvmin.setValidator(QDoubleValidator())
        
        self.pvmaxLabel = QLabel('ρ<sub>max</sub>')
        self.pvmax = QLineEdit()
        self.pvmax.setValidator(QDoubleValidator())
        
        self.pvminIPLabel = QLabel('ΙP<sub>min</sub>')
        self.pvminIPLabel.setVisible(False)
        self.pvminIP = QLineEdit()
        self.pvminIP.setValidator(QDoubleValidator())
        self.pvminIP.setVisible(False)
        
        self.pvmaxIPLabel = QLabel('IP<sub>max</sub>')
        self.pvmaxIPLabel.setVisible(False)
        self.pvmaxIP = QLineEdit()
        self.pvmaxIP.setValidator(QDoubleValidator())
        self.pvmaxIP.setVisible(False)
        
        self.pParams = {'index':0, 'vmin':None, 'vmax':None, 'threed':False}
        self.pParamsIP = {'index':0, 'vmin':None, 'vmax':None, 'threed':False}
        def prescaleBtnFunc():
            if self.project is not None:
                self.pParams['vmin'] = float(self.pvmin.text()) if self.pvmin.text() != '' else None
                self.pParams['vmax'] = float(self.pvmax.text()) if self.pvmax.text() != '' else None
                self.pParamsIP['vmin'] = float(self.pvminIP.text()) if self.pvminIP.text() != '' else None
                self.pParamsIP['vmax'] = float(self.pvmaxIP.text()) if self.pvmaxIP.text() != '' else None    
                self.plotPseudo()
                if self.project.typ[0] == 'c':
                    self.plotPseudoIP()
            QApplication.processEvents()
        self.prescaleBtn = QPushButton('Apply')
        self.prescaleBtn.setAutoDefault(True)
        self.prescaleBtn.clicked.connect(prescaleBtnFunc)
    
        # tab focus events functions
        self.currentImportTab = 0
        self.currentTab = 0
        def logImportTab(index):
            if index != 1 and self.currentImportTab == 1: # elec tab looses focus
                self.updateElec()
            self.currentImportTab = index
        def logTab(index):
            if index != 0 and self.currentImportTab == 1:
                # elec tab still visible but we moved out to another
                # higher level tab
                self.tabImporting.setCurrentIndex(0)
            self.currentTab = index
        self.tabImporting.currentChanged.connect(logImportTab)
        self.tabs.currentChanged.connect(logTab)

        if pvfound:
            self.pseudoFrame = QFrame()
            vlayout = QVBoxLayout()
            self.pseudoPlotter = QtInteractor(self.pseudoFrame)
            vlayout.addWidget(self.pseudoPlotter.interactor)
            self.pseudoFrame.setLayout(vlayout)
        self.mwPseudo = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        
        if pvfound:
            self.pseudoFrameIP = QFrame()
            vlayout = QVBoxLayout()
            self.pseudoPlotterIP = QtInteractor(self.pseudoFrameIP)
            vlayout.addWidget(self.pseudoPlotterIP.interactor)
            self.pseudoFrameIP.setLayout(vlayout)
            self.pseudoFrameIP.setVisible(False)
        self.mwPseudoIP = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        self.mwPseudoIP.setVisible(False)

        # layout
        self.hbox1 = QHBoxLayout()
        self.hbox1.addWidget(self.dimGroup)
        self.hbox1.addWidget(self.titleLabel)
        self.hbox1.addWidget(self.titleEdit)
        self.hbox1.addWidget(self.dateLabel)
        self.hbox1.addWidget(self.dateEdit)
        # self.hbox1.addWidget(self.saveProjectBtn)
        # self.hbox1.addWidget(self.loadProjectBtn)
        # self.hbox1.addWidget(self.restartBtn)

        self.hbox2 = QHBoxLayout()
        self.hbox2.addWidget(self.dimInvGroup)
        self.hbox2.addWidget(self.timeLapseCheck)
        self.hbox2.addWidget(self.batchCheck)
        self.hbox2.addWidget(self.boreholeCheck)
        self.hbox2.addWidget(self.pseudo3DCheck)
        self.hbox2.addWidget(self.regular3DCheck)

        self.hbox4 = QHBoxLayout()
        self.hbox4.addWidget(self.wdBtn)
        self.hbox4.addWidget(self.ftypeComboLabel)
        self.hbox4.addWidget(self.ftypeCombo)
#        self.hbox4.addWidget(self.spacingEdit)
        self.hbox4.addWidget(self.importDataBtn)
        self.hbox4.addWidget(self.importDataRecipBtn)
        self.hbox4.addWidget(self.lineSpacingLabel)
        self.hbox4.addWidget(self.lineSpacing)
        self.hbox4.addWidget(self.create3DBtn)
        self.hbox4.addWidget(self.invNowBtn)
        
        self.hbox5 = QHBoxLayout()
        self.hbox5.setAlignment(Qt.AlignRight)
        self.hbox5.addWidget(self.ipCheck, Qt.AlignLeft)
        self.hbox5.addWidget(self.mergElecCheck)
        self.hbox5.addWidget(self.psContourCheck)
        self.hbox5.addWidget(self.pvminLabel)
        self.hbox5.addWidget(self.pvmin)
        self.hbox5.addWidget(self.pvmaxLabel)
        self.hbox5.addWidget(self.pvmax)
        self.hbox5.addWidget(self.pvminIPLabel)
        self.hbox5.addWidget(self.pvminIP)
        self.hbox5.addWidget(self.pvmaxIPLabel)
        self.hbox5.addWidget(self.pvmaxIP)
        self.hbox5.addWidget(self.prescaleBtn)
        self.hbox5.addWidget(self.fnamesComboLabel)
        self.hbox5.addWidget(self.fnamesCombo)

        self.metaLayout = QVBoxLayout()
        self.metaLayout.addLayout(self.hbox1)
        self.metaLayout.addLayout(self.hbox2)
        self.metaLayout.addLayout(self.hbox4)
        self.metaLayout.addLayout(self.hbox5)
        self.tabImportingDataLayout.addLayout(self.metaLayout, 40)

        self.pseudoLayout = QStackedLayout()
        self.mwLayoutWidget = QWidget()
        self.mwLayout = QHBoxLayout()
        self.mwLayout.addWidget(self.mwPseudo, 50)
        self.mwLayout.addWidget(self.mwPseudoIP, 50)
        self.mwLayoutWidget.setLayout(self.mwLayout)
        self.pseudoLayout.addWidget(self.mwLayoutWidget)
        if pvfound:
            self.pvLayoutWidget = QWidget()
            self.pvLayout = QHBoxLayout()
            self.pvLayout.addWidget(self.pseudoFrame, 50)
            self.pvLayout.addWidget(self.pseudoFrameIP, 50)
            self.pvLayoutWidget.setLayout(self.pvLayout)
            self.pseudoLayout.addWidget(self.pvLayoutWidget)
        self.tabImportingDataLayout.addLayout(self.pseudoLayout, 60)
        self.tabImportingData.setLayout(self.tabImportingDataLayout)


        #%% sub tab with elec and topo informations
        self.tabImportingTopo = QWidget()
        self.tabImporting.addTab(self.tabImportingTopo, 'Electrodes (XYZ/Topo)')

        # electrode table
        class ElecTable(QTableWidget):
            def __init__(self, nrow=2, selfInit=False, parent=None):
                """ if selfInit is true, it will automatically add rows if tt
                is bigger than the actual rows
                """
                ncol = 5
                super(ElecTable, self).__init__(nrow, ncol)
                self.nrow = nrow
                self.ncol = ncol
                self.parent = parent
                self.headers = headers=['Label','X','Y','Z','Buried']
                self.selfInit = selfInit # if True, table is auto-expanded
                self.initTable(np.array([['1','2'],['',''],['',''],['','']]))
                self.horizontalHeader().sortIndicatorChanged.connect(self.setAllBuried)
                self.setHorizontalHeaderLabels(headers)
                self.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
                self.useNarray = False
                self.xyz = self.getTable()[['x','y','z']].values

            def addRow(self):
                self.nrow += 1
                self.setRowCount(self.nrow)
    
            def setBuried(self, vals=None):
                if vals is None:
                    vals = np.zeros(self.nrow, dtype=bool)
                for i in range(len(vals)):
                    self.checkBoxWidget = QWidget()
                    checkBoxLayout = QHBoxLayout()
                    checkBoxLayout.setContentsMargins(5,5,5,5)
                    checkBoxLayout.setAlignment(Qt.AlignCenter)
                    self.buriedCheck = QCheckBox()
                    self.buriedCheck.setChecked(bool(vals[i]))
                    checkBoxLayout.addWidget(self.buriedCheck)
                    self.checkBoxWidget.setLayout(checkBoxLayout)
                    self.setCellWidget(i, 4, self.checkBoxWidget) # column 4 == Buried

            def setAllBuried(self, colIndex):
                if colIndex == 4:
                    for i in range(self.nrow):
                        self.buriedCheck = self.cellWidget(i, 4).findChildren(QCheckBox)[0]
                        if self.buriedCheck.isChecked() is True:
                            self.buriedCheck.setChecked(False)
                        else:
                            self.buriedCheck.setChecked(True)
                            

            def keyPressEvent(self, e):
                if (e.modifiers() == Qt.ControlModifier) & (e.key() == Qt.Key_V):
                    cell = self.selectedIndexes()[0]
                    c0, r0 = cell.column(), cell.row()
                    self.paste(c0, r0)
                elif e.modifiers() != Qt.ControlModifier:
                    cell = self.selectedIndexes()[0]
                    c0, r0 = cell.column(), cell.row()
                    self.editItem(self.item(r0,c0))

            def paste(self, c0, r0):
                text = QApplication.clipboard().text() # get clipboard
                tt = []
                for row in text.split('\n'):
                    trow = row.split()
                    if len(trow) > 0:
                        tt.append(trow)
                tt = np.array(tt)
                if self.isColumnHidden(2) == True: # 'y' is hidden
                    tt = np.c_[tt[:,0], np.zeros(tt.shape[0]), tt[:,1:]]
                if self.isColumnHidden(0) == True: # label hidden
                    tt = np.c_[(1 + np.arange(tt.shape[0])).astype(str), tt]
                pdebug('elecTable.paste():', tt)
                if np.sum(tt.shape) > 0:
                    if self.selfInit is True:
                        self.initTable(tt)
                    else:
                        self.setTable(tt, c0, r0)

            def initTable(self, tt):
                pdebug('elecTable.initTable():', tt)
                self.clear() # this clears out the headers as well
                self.nrow = tt.shape[0]
                print(self.nrow)
                if self.nrow > 10000: # set a hard cap on the number of rows displayed to avoid segmentation fault 
                    self.nrow = 10000
                    self.useNarray = True # use the numpy data array to index the data 
                    pdebug('The maximum number of qt table entries have been exceeded, shown table is capped at 10000')
                self.setRowCount(self.nrow)
                self.setHorizontalHeaderLabels(self.headers)
                self.setTable(tt[0:self.nrow])
                self.setBuried()
            
            def setTable(self, tt, c0=0, r0=0):
                if type(tt) == type(pd.DataFrame()):
                    tt = tt[['label','x','y','z']].values
                pdebug('elecTable.setTable():', self.nrow, self.ncol, tt.shape)
                nanFlag = False # to inform the user that their topography is messed up!
                for i in range(c0, min([self.ncol, c0+tt.shape[1]])):
                    for j in range(r0, min([self.nrow, r0+tt.shape[0]])):
                        self.setItem(j,i,QTableWidgetItem(str(tt[j-r0, i-c0])))
                        if self.item(j,i).text() == 'nan': #tried np.isnan(tt[j-r0, i-c0]) == True but as cell is initially an empty string, encountered error
                            self.item(j,i).setBackground(QColor(255,0,0))
                            nanFlag = True
                if nanFlag == True:
                    self.parent.infoDump('Missing topography points found! Use "Interpolate missing topo"')


            def getTable(self):
                """Return dataframe.
                """
                cols = ['label','x','y','z','buried']
                df = pd.DataFrame(np.zeros((self.nrow, self.ncol)), columns=cols)
                for i in range(self.nrow):
                    for j in range(self.ncol):
                        col = cols[j]
                        if j == 0:
                            df.loc[i,col] = self.item(i,j).text()
                        elif j < 4:
                            if self.item(i,j) is None:
                                df.loc[i, col] = np.nan
                            else:
                                t = self.item(i,j).text()
                                df.loc[i, col] = float(t) if t != '' else np.nan
                        elif j == 4:
                            self.buriedCheck = self.cellWidget(i, 4).findChildren(QCheckBox)[0]
                            if self.buriedCheck.isChecked() is True:
                                df.loc[i, 'buried'] = True
                            else:
                                df.loc[i, 'buried'] = False
                return df.astype({'label':str, 'x':float, 'y':float, 'z':float, 'buried':bool}).reset_index(drop=True)
            

            def readTable(self, fname, nbElec=None):
                # identify if we have header (recommended) or not
                with open(fname, 'r') as f:
                    line = f.readline().split(',')[0]
                try:
                    float(line)
                    header = None
                    pdebug('elecTable.readTable: no header')
                except Exception:
                    header = 'infer'
                    pdebug('elecTable.readTable: header provided')
                df = pd.read_csv(fname, header=header)
                newcols = [str(a).lower() for a in df.columns]
                df = df.rename(columns = dict(zip(df.columns, newcols)))
                
                # let's ensure we always have 4 columns
                columns = ['x','y','z','buried']
                df2 = pd.DataFrame(np.zeros((df.shape[0], 4)), columns=columns)
                if header is not None:
                    for i, l in enumerate(columns):
                        if l in df.columns:
                            df2[l] = df[l]
                else: # provided for backward compatibility
                    if df.shape[1] == 2:# assume XZ
                        df2['x'] = df.values[:,0]
                        df2['z'] = df.values[:,1]
                    elif df.shape[1] == 3:
                        if 'y' in self.headers: # assume XYZ
                            df2['x'] = df.values[:,0]
                            df2['y'] = df.values[:,1]
                            df2['z'] = df.values[:,2]
                        else: # assume XZ buried
                            df2['x'] = df.values[:,0]
                            df2['z'] = df.values[:,1]
                            df2['buried'] = df.values[:,2]
                    elif df.shape[1] == 4: # XYZ buried
                        df2['x'] = df.values[:,0]
                        df2['y'] = df.values[:,1]
                        df2['z'] = df.values[:,2]
                        df2['buried'] = df.values[:,3]
                
                # 3 options for constructing label here... 
                if 'line' in df.columns and 'number' in df.columns:
                    #construct label column from line number and electrode number 
                    line = df.line.values # line number array 
                    elecid  = df.number.values # electrode number array 
                    label = ['1 1'] * len(df)
                    for i in range(len(df)):
                        label[i] = '%i %i'%(line[i],elecid[i])
                    df2.insert(0,'label',label)
                    pdebug('Constructed label column from individual line and electrode numbers')
                elif 'label' in df.columns: # use label column directly 
                    a = df['label'].values[0]
                    if isinstance(a, str) is False:
                        df['label'] = df['label'].astype(int).astype(str)
                    df2.insert(0, 'label', df['label'].values)
                else: # fallback if no label column provided/found  
                    df2.insert(0, 'label', (1 + np.arange(df.shape[0])).astype(str))
                    # if self.parent.r2.typ[-1] == 't': # 3D
                    #     df2_groups = df2.groupby('y') # finding number of lines based on y vals
                    #     df2_grouppedData = [df2_groups.get_group(x) for x in df2_groups.groups]
                    #     for i, group in enumerate(df2_grouppedData):
                    #         group['label'] = (np.arange(group.shape[0]) + 1).astype(str)
                    #         group.loc[:, 'label'] = str(i+1) + ' ' + group['label']
                    #     df2_temp = pd.concat(df2_grouppedData)
                    #     df2['label'] = df2_temp['label'].values
                    #     df2.loc[:, 'label'] = '1 ' + df2['label']

                
                pdebug('elecTable.readTable():\n', df2.head())
                self.xyz = df2[['x','y','z']].values # just the xyz numpy array, called when creating meshes 
                tt = df2.values
                if nbElec is not None:
                    if tt.shape[0] < nbElec:
                        self.parent.errorDump('The file must have exactly ' + \
                                  str(nbElec) + ' lines (same number as number of electrodes).')
                        return
                if 'Buried' in self.headers:
                    if 1 <= len(np.unique(tt[:,-1])) <= 2: # only 1 and 0
                        self.initTable(tt[:,:-1])
                        self.setBuried(tt[:,-1])
                    else:
                        self.parent.errorDump('the "buried" column should contains only 0 or 1')
                else:
                    self.initTable(tt) # topoTable doesn't have "buried" thus last column is not needed 


        self.elecTable = ElecTable(parent=self)
        self.elecTable.setColumnHidden(2, True)
        self.elecLabelTextP1 = '<i>Add electrode position. Use <code>Ctrl+V</code> to paste or import from CSV with headers matching: label,x,y,z,buried.' \
            'The last column (buried) is 1 if checked (= buried electrode) and 0 if not (=surface electrode).' \
                'You can also use the form below to generate ' \
                    'regular electrode spacing. <b>Click on the <font color="red">"Buried"</font> table header to check/unchek all.</b></i>'
        self.elecLabelTextP2 = ''
        self.elecLabel = QLabel(self.elecLabelTextP1)
        self.elecLabel.setWordWrap(True)

        def importElecBtnFunc():
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingTopo,'Open File', directory=self.datadir)
            if fname != '':
                if self.iForward is False:
                    nbElec = int(self.nbElecEdit.text())
                else:
                    nbElec = None
                self.elecTable.readTable(fname, nbElec=nbElec)
                self.writeLog('k.importElec("{:s}")'.format(fname))
        self.importElecBtn = QPushButton('Import from CSV files with headers: label (or line, number), x, y, z, buried')
        self.importElecBtn.setAutoDefault(True)
        self.importElecBtn.clicked.connect(importElecBtnFunc)
        
        self.nbElecEdit = QLineEdit()
        self.nbElecEdit.setValidator(QDoubleValidator())
        self.nbElecEdit.setEnabled(False)
        self.nbElecLabel = QLabel('Number of electrodes:')
        
        self.elecDxLabel = QLabel('X spacing:')
        self.elecDxEdit = QLineEdit('1.0')
        self.elecDxEdit.setValidator(QDoubleValidator())

        self.elecDzLabel = QLabel('Z spacing:')
        self.elecDzEdit = QLineEdit('0.0')
        self.elecDzEdit.setValidator(QDoubleValidator())
        
        self.elecLineLabel = QLabel('Nb Lines:')
        self.elecLineEdit = QLineEdit('1')
        self.elecLineEdit.setValidator(QIntValidator())
        self.elecLineEdit.setEnabled(False)
        
        self.elecLineSpacingLabel = QLabel('Line spacing:')
        self.elecLineSpacingEdit = QLineEdit('2')
        self.elecLineSpacingEdit.setValidator(QDoubleValidator())
        self.elecLineSpacingEdit.setEnabled(False)

        def elecGenButtonFunc():
            nb = int(self.nbElecEdit.text())
            dx = float(self.elecDxEdit.text())
            dz = float(self.elecDzEdit.text())
            nline = int(self.elecLineEdit.text())
            lineSpacing = float(self.elecLineSpacingEdit.text())
            
            # df = pd.DataFrame()
            # if self.project.iForward: # if we are in forward mode, we overwrite the
            # # electrode labels, otherwise, we keep them and just fill
            #     df['label'] = 1 + np.arange(nbElec) # label
            #     df['label'] = df['label'].astype(int).astype(str)
            # else:
            #     df['label'] = self.project.elec['label'].values
            # df['x'] = np.linspace(0.0, (nbElec-1)*dx, nbElec)
            # df['y'] = np.linspace(0.0, (nbElec-1)*dy, nbElec)
            # df['z'] = np.linspace(0.0, (nbElec-1)*dz, nbElec)
            
            try:
                self.project.generateElec(nb, dx, dz, nline, lineSpacing)
                self.writeLog('k.generateElec(nb={:s},dx={:s},dz={:s},nline={:s},lineSpacing={:s})'.format(
                str(nb), str(dx), str(dz), str(nline), str(lineSpacing)))
                self.elecTable.initTable(self.project.elec)
            except Exception as e:
                self.errorDump(e)
        self.elecGenButton = QPushButton('Generate')
        self.elecGenButton.setAutoDefault(True)
        self.elecGenButton.clicked.connect(elecGenButtonFunc)

        self.topoTable = ElecTable(selfInit=True, parent=self)
        self.topoTable.setColumnHidden(2, True) # Y
        self.topoTable.setColumnHidden(0, True) # elec label
        self.topoTable.setColumnHidden(4, True) # buried
        topoLabel = QLabel('<i>Add additional surface points. \
                           You can use <code>Ctrl+V</code> to paste directly \
                           into a cell.</i>')
                           
        def topoBtnFunc():
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingTopo,'Open File', directory=self.datadir)
            if fname != '':
                self.topoTable.readTable(fname)
        self.topoBtn = QPushButton('Import from CSV file with headers: x, y, z')
        self.topoBtn.setAutoDefault(True)
        self.topoBtn.clicked.connect(topoBtnFunc)
        
        def topoAddRowBtnFunc():
            self.topoTable.addRow()
        self.topoAddRowBtn = QPushButton('Add Row')
        self.topoAddRowBtn.clicked.connect(topoAddRowBtnFunc)
        
        
        self.topoInterpBtn = QPushButton('Interpolate missing topo')
        self.topoInterpBtn.clicked.connect(self.topoInterpBtnFunc)
        self.topoInterpBtn.setToolTip('The topography points will be used to find '
                                  'the Z coordinates of the electrodes using '
                                  'linear interpolation.')
        
        
        # layout
        self.topoLayout = QVBoxLayout()

        self.elecGenLayout = QHBoxLayout()
        self.elecGenLayout.addWidget(self.nbElecLabel)
        self.elecGenLayout.addWidget(self.nbElecEdit)
        self.elecGenLayout.addWidget(self.elecDxLabel)
        self.elecGenLayout.addWidget(self.elecDxEdit)
        self.elecGenLayout.addWidget(self.elecDzLabel)
        self.elecGenLayout.addWidget(self.elecDzEdit)
        self.elecGenLayout.addWidget(self.elecLineLabel)
        self.elecGenLayout.addWidget(self.elecLineEdit)
        self.elecGenLayout.addWidget(self.elecLineSpacingLabel)
        self.elecGenLayout.addWidget(self.elecLineSpacingEdit)
        self.elecGenLayout.addWidget(self.elecGenButton)
        self.topoLayout.addWidget(self.elecLabel)
        self.topoLayout.addLayout(self.elecGenLayout)
        self.topoLayout.addWidget(self.importElecBtn)
        self.topoLayout.addWidget(self.elecTable)
        
        self.topoLayout.addWidget(topoLabel)
        self.topoBtnLayout = QHBoxLayout()
        self.topoBtnLayout.addWidget(self.topoBtn, 70)
        self.topoBtnLayout.addWidget(self.topoAddRowBtn, 10)
        self.topoBtnLayout.addWidget(self.topoInterpBtn, 20)
        self.topoLayout.addLayout(self.topoBtnLayout)
        self.topoLayout.addWidget(self.topoTable)

        self.tabImportingTopo.setLayout(self.topoLayout)
        self.tabImporting.setTabEnabled(1, False)


        #%% sub tab for custom parser
        self.customParser = QWidget()
        self.tabImporting.addTab(self.customParser, 'Custom Parser')
        
        self.delimiter = ''
        def delimFunc(index):
            self.delimiterBox.hide()
            self.delimiterBox.clear()
            if index == 0:
                self.delimiter = ''
            if index == 1: 
                self.delimiter = ','
            elif index == 2: 
                self.delimiter = '\t'
            elif index == 3: 
                self.delimiter = '\s+'
            elif index == 4: 
                self.delimiter = ';'
            elif index == 5:
                self.delimiterBox.show()
                
        def customDelimFunc():
            self.delimiter = self.delimiterBox.text()

        self.delimCombo = QComboBox()
        self.delimCombo.setMinimumWidth(200)
        self.delimCombo.addItems(['Select delimiter','Comma','Tab','Space','Semicolon', 'Other'])
        self.delimCombo.activated.connect(delimFunc)
        
        self.delimiterBox = QLineEdit()
        self.delimiterBox.setFixedWidth(100)
        self.delimiterBox.editingFinished.connect(customDelimFunc)
        self.delimiterBox.hide()
        
        self.delimLabel = QLabel('Delimiter:')
#        delimiterEdit = QLineEdit('')
#        delimiterEdit.setToolTip(r'For tab delimited data use: \t')
        self.skipRowsLabel = QLabel('Number of header to skip:')
        self.skipRowsEdit = QLineEdit('0')
        self.skipRowsEdit.setValidator(QIntValidator())
        self.nrowsLabel = QLabel('Number of rows to read:')
        self.nrowsEdit = QLineEdit('')
        self.nrowsEdit.setValidator(QIntValidator())

        self.fnameManual = None
        def openFileBtnFunc(file):
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingTopo,'Open File', directory=self.datadir)
            if fname != '':
                self.fnameManual = fname
                self.openFileBtn.setText(fname + ' (Click to change)')
                parseBtnFunc()
                self.parseBtn.setEnabled(True)
        self.openFileBtn = QPushButton('Open File')
        self.openFileBtn.setAutoDefault(True)
        self.openFileBtn.clicked.connect(openFileBtnFunc)

        def parseBtnFunc():
            if self.fnameManual is None:
                self.errorDump('Select a file to parse first.')
                return
            try:
                delimiter = self.delimiter
                delimiter = None if delimiter == '' else delimiter
                skipRows = self.skipRowsEdit.text()
                skipRows = None if skipRows == '' else int(skipRows)
                nrows = self.nrowsEdit.text()
                nrows = None if nrows == '' else int(nrows)
                self.parserTable.readTable(self.fnameManual, delimiter=delimiter,
                                          skiprows=skipRows, nrows=nrows)
                fillBoxes(boxes) # last one is elecSpacingEdit
                self.infoDump('Parsing successful.')
            except ValueError as e:
                self.errorDump('Parsing error:' + str(e) + ' - Incorrect delimiter or skip headers')

        self.parseBtn = QPushButton('Reorder')
        self.parseBtn.setEnabled(False)
        self.parseBtn.setAutoDefault(True)
        self.parseBtn.clicked.connect(parseBtnFunc)

        # have qcombobox to be read for each columns
        self.aBoxLabel = QLabel('A (or C1):')
        self.bBoxLabel = QLabel('B (or C2):')
        self.mBoxLabel = QLabel('M (or P1):')
        self.nBoxLabel = QLabel('N (or P2):')
        self.vpBoxLabel = QLabel('Vp Potential Difference:')
        self.InBoxLabel = QLabel('In Current:')
        self.resistBoxLabel = QLabel('Transfer Resistance:')
#        ipStartBoxLabel = QLabel('IP start column') # we don't need these for now, since DCA only works with syscal files
#        ipEndBoxLabel = QLabel('IP end column')
        self.chargeabilityBoxLabel = QLabel('Chargeability:')
        self.phaseBoxLabel = QLabel('Phase shift:')
        self.devErrLabel = QLabel('Stacking Error (dev):')
        self.resErrLabel = QLabel('Resistance Error:')
        self.phaseLabel = QLabel('Phase/Chargeability Error:')
        
        
#        boxesLabels = [self.aBoxLabel, bBoxLabel, mBoxLabel, nBoxLabel, vpBoxLabel, InBoxLabel, resistBoxLabel, ipStartBoxLabel,
#                 ipEndBoxLabel, chargeabilityBoxLabel, phaseBoxLabel]#, elecSpacingLabel]
        boxesLabels = [self.aBoxLabel, self.bBoxLabel, self.mBoxLabel, 
                       self.nBoxLabel, self.vpBoxLabel, self.InBoxLabel,
                       self.resistBoxLabel, self.chargeabilityBoxLabel,
                       self.phaseBoxLabel, self.devErrLabel, self.resErrLabel, 
                       self.phaseLabel]

        self.aBox = QComboBox()
        self.bBox = QComboBox()
        self.mBox = QComboBox()
        self.nBox = QComboBox()
        self.vpBox = QComboBox()
        self.InBox = QComboBox()
        self.resistBox = QComboBox()
#        ipStartBox = QComboBox() # we don't need these for now, since DCA only works with syscal files
#        ipEndBox = QComboBox()
        self.chargeabilityBox = QComboBox()
        self.chargeabilityBox.setToolTip('input the column containing chargeability (mV/V) values')
        self.phaseBox = QComboBox()
        self.phaseBox.setToolTip('input the column containing phase shift (mRad) values')
        self.devErrBox = QComboBox()
        self.resErrBox = QComboBox()
        self.phaseErrBox = QComboBox()
#        elecSpacingEdit = QLineEdit('')
#        elecSpacingEdit.setEnabled(False)
#        elecSpacingEdit.setValidator(QDoubleValidator())
#        elecSpacingEdit.setFixedWidth(80)
#        elecSpacingEdit.setToolTip('Number to divide the selected columns to get electrode number.')

#        boxes = [aBox, bBox, mBox, nBox, vpBox, InBox, resistBox, ipStartBox,
#                 ipEndBox, chargeabilityBox, phaseBox]#, elecSpacingEdit]
        boxes = [self.aBox, self.bBox, self.mBox, self.nBox, self.vpBox,
                 self.InBox, self.resistBox, self.chargeabilityBox, 
                 self.phaseBox, self.devErrBox, self.resErrBox, self.phaseErrBox]#, elecSpacingEdit]

        def fillBoxes(bs):
            for b in bs:
                b.clear()
                choices = self.parserTable.headers
                for i, choice in enumerate(choices):
                    if i == 0:
                        b.addItem('Select if available')
                    b.addItem(str(choice))

        def getBoxes(bs):
            cols = []
            for b in bs:
                cols.append(b.currentIndex()-1) # -1 because of 1st item
            return np.array(cols)

        # add a qtableview or create a custom class
        class ParserTable(QTableWidget):
            def __init__(self, nrow=10, ncol=10):
                super(ParserTable, self).__init__(nrow, ncol)
                self.nrow = nrow
                self.ncol = ncol
                self.headers = []

            def readTable(self, fname='', delimiter=None, skiprows=None, nrows=None):
                if fname != '':
                    df = pd.read_csv(fname, delimiter=delimiter, skiprows=skiprows, nrows=nrows, encoding_errors='ignore')
                    df = df.reset_index() # in case all parse columns goes into the index (don't know why)
                    self.setRowCount(df.shape[0])
                    self.setColumnCount(df.shape[1])
                    self.headers = df.columns.values.astype(str) # make sure it's string
                    self.setHorizontalHeaderLabels(self.headers)
                    tt = df.values
                    self.setTable(tt)

            def setTable(self, tt):
                # paste clipboard to qtableView
                self.setRowCount(tt.shape[0])
                self.setColumnCount(tt.shape[1])
                self.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
                for i in range(tt.shape[1]):
                    for j in range(tt.shape[0]):
                        self.setItem(j,i,QTableWidgetItem(str(tt[j, i])))


        self.parserTable = ParserTable()

        def importBtnFunc():
            self.restartFunc()
            colIndex = []
            newHeaders = []
            vals = getBoxes([self.aBox, self.bBox, self.mBox, self.nBox])
            if (np.sum(vals > 0) == 4) & (len(np.unique(vals)) == 4):
                colIndex.append(vals)
                newHeaders.append(['a','b','m','n'])
            else:
                self.errorDump('Please select columns for electrodes.')
                return
            vals = getBoxes([self.resistBox])
            if vals[0] > 0:
                colIndex.append(vals)
                newHeaders.append(['resist'])
            else:
                vals = getBoxes([self.vpBox, self.InBox])
                if all(vals > 0) is True:
                    colIndex.append(vals)
                    newHeaders.append(['vp','i'])
                else:
                    self.errorDump('Please select columns for Vp, In or Resist.')
                    return
            vals = getBoxes([self.chargeabilityBox])
            if vals[0] > 0:
                colIndex.append(vals)
                newHeaders.append(['ip'])
                self.phiConvFactor.setText('1')
                self.phiConvFactor.setEnabled(True)
                self.phiConvFactorlabel.setEnabled(True)
                self.inputPhaseFlag = False
            else:
                vals = getBoxes([self.phaseBox])
                if vals[0] > 0:
                    colIndex.append(vals)
                    newHeaders.append(['ip'])
                    self.phiConvFactor.setText('')
                    self.phiConvFactor.setEnabled(False)
                    self.phiConvFactorlabel.setEnabled(False)
                    self.inputPhaseFlag = True
                else:
                    self.ipCheck.setChecked(False)
            vals = getBoxes([self.devErrBox])
            if vals[0] > 0:
                colIndex.append(vals)
                newHeaders.append(['dev'])
            vals = getBoxes([self.resErrBox])
            if vals[0] > 0:
                colIndex.append(vals)
                newHeaders.append(['magErr'])
            vals = getBoxes([self.phaseErrBox])
            if vals[0] > 0:
                colIndex.append(vals)
                newHeaders.append(['phiErr'])
                
            # currently not importing each IP columns (M1 -> M..) so no
            # decay curve analysis can be performed

            colIndex = np.hstack(colIndex)
            newHeaders = np.hstack(newHeaders)

            def parserFunc(fname):
                # retrieve usefull values
                try:
                    delimiter = self.delimiter
                    delimiter = None if delimiter == '' else delimiter
                    skipRows = self.skipRowsEdit.text()
                    skipRows = None if skipRows == '' else int(skipRows)
                    nrows = self.nrowsEdit.text()
                    nrows = None if nrows == '' else int(nrows)
                    espacing = None #if elecSpacingEdit.text() == '' else float(elecSpacingEdit.text())
    
                    # parse
                    print('delimiter=', delimiter)
                    df = pd.read_csv(fname, delimiter=delimiter, skiprows=skipRows, nrows=nrows)
                    df = df.reset_index() # solve issue all columns in index
                    oldHeaders = df.columns.values[colIndex]
                    df = df.rename(columns=dict(zip(oldHeaders, newHeaders)))
                    if 'resist' not in df.columns:
                        df['resist'] = df['vp']/df['i']
                    if 'ip' not in df.columns:
                        df['ip'] = 0
                    elif self.inputPhaseFlag == True:
                        df['ip'] *= -1 # if the input ip values are already phase, in custom parser only!
                        
                    ################ better method (adopted from syscalParser)
                    
                    array = df[['a','b','m','n']].values
            
                    # get unique electrode positions and create ordered labels for them
                    val = np.sort(np.unique(array.flatten()))
                    elecLabel = 1 + np.arange(len(val))
                    searchsoterdArr = np.searchsorted(val, array)
                    newval = elecLabel[searchsoterdArr] # magic ! https://stackoverflow.com/questions/47171356/replace-values-in-numpy-array-based-on-dictionary-and-avoid-overlap-between-new
                    df.loc[:,['a','b','m','n']] = newval # assign new label
                    zval = np.zeros_like(val)
                    yval = np.zeros_like(val) # 2D so Y values are all zeros
                    elec = np.c_[val, yval, zval]
                    remoteFlags_p1 = np.array([-9999999, -999999, -99999,-9999,-999])
                    remoteFlags_p2 = np.array([9999999, 999999, 99999, 9999, 999])
                    iremote_p1 = np.in1d(val, remoteFlags_p1)
                    elec[iremote_p1, 0] = -99999
                    iremote_p2 = np.in1d(val, remoteFlags_p2)
                    elec[iremote_p2, 0] = 99999
                    
                    ##### Old way
                        
                    # array = df[['a','b','m','n']].values.copy()
                    # arrayMin = np.min(np.unique(np.sort(array.flatten())))
                    # if arrayMin != 0:
                    #     array -= arrayMin
                    # if espacing is None:
                    #     espacing = np.unique(np.sort(array.flatten()))[1] - np.unique(np.sort(array.flatten()))[0]
                    # array = np.round(array/espacing+1).astype(int)
                    # df[['a','b','m','n']] = array
                    # imax = int(np.max(array))
                    # elec = np.zeros((imax,3))
                    # elec[:,0] = np.arange(0,imax)*espacing
                    
                    ################
                    
                    self.nbElecEdit.setText('%s' % (len(elec)))
    #                self.nbElecEdit.setEnabled(False)
                    self.elecDxEdit.setText('%s' % (espacing))
                    self.activateTabs(True)
                    return elec, df
                except:
                    self.activateTabs(False)
                    self.errorDump("Import Failed: 'nan' values must be removed before importation. Use the 'Number of rows to read or skip' to remove 'nan's.")
            self.parser = parserFunc

            # test the parser
            elec, df = parserFunc(self.fnameManual)
            pdebug('custom parserFunc: shapes = ', elec.shape, df.shape)

            if (self.project.iTimeLapse is False) & (self.project.iBatch is False):
                self.importFile(self.fnameManual)
            self.ftypeCombo.setCurrentIndex(12)
            self.tabImporting.setCurrentIndex(0)

        self.importBtn = QPushButton('Import Dataset')
        self.importBtn.setAutoDefault(True)
        self.importBtn.clicked.connect(importBtnFunc)


        # layout
        self.parserLayout = QVBoxLayout()
        self.parserOptions = QHBoxLayout()
        self.columnsAssign = QGridLayout()

        self.parserLayout.addWidget(self.openFileBtn)
        self.parserOptions.addWidget(self.delimLabel)
        self.parserOptions.addWidget(self.delimCombo)
        self.parserOptions.addWidget(self.delimiterBox)
        self.parserOptions.addWidget(self.skipRowsLabel)
        self.parserOptions.addWidget(self.skipRowsEdit)
        self.parserOptions.addWidget(self.nrowsLabel)
        self.parserOptions.addWidget(self.nrowsEdit)
        self.parserOptions.addWidget(self.parseBtn)
        self.parserLayout.addLayout(self.parserOptions)

        self.parserLayout.addWidget(self.parserTable)
        for i in range(len(boxes)):
            c = (i % 3)*2 # in 2*3 columns (with labels)
            r = int(i/3)
            self.columnsAssign.addWidget(boxesLabels[i], r, c, Qt.AlignRight)
            self.columnsAssign.addWidget(boxes[i], r, c+1)
        self.parserLayout.addLayout(self.columnsAssign)
        self.parserLayout.addWidget(self.importBtn)

        self.customParser.setLayout(self.parserLayout)


#%% ===================== pre-processing tab =====================
        self.tabPreProcessing = QTabWidget()
        self.tabs.addTab(self.tabPreProcessing, 'Pre-processing')
        self.tabs.setTabEnabled(1, False)

        
#%% reciprocal filtering tab (or just normal filtering if no reciprocal)
        self.recipErrorWidget  = QWidget()
        self.tabPreProcessing.addTab(self.recipErrorWidget , 'Reciprocal error analysis')
        self.tabPreProcessing.setTabEnabled(0, True)

        self.recipErrorLabel = QLabel('<b>Remove datapoints that have reciprocal error larger than what you prefer.</b><br>Either select (<i>click on the dots to select them</i>) the points on the pseudo section below or choose a percentage threshold or both!</br>')
        self.recipErrorLabel.setAlignment(Qt.AlignLeft)
        
        self.recipErrorfnamesComboLabel = QLabel('Select a dataset:')
        self.recipErrorfnamesComboLabel.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
        
        def recipErrorfnamesComboFunc(index):
            pdebug('recipErrorfnamesComboFunc()')
            if index == 0:
                self.recipErrApplyToAll = True
                self.recipErrDataIndex = -1
                self.plotManualFiltering(0)
                self.errHist(0)
            elif index > 0: # show/hide makes the index = -1
                self.recipErrApplyToAll = False
                self.plotManualFiltering(index-1)
                self.errHist(index-1)
                self.recipErrDataIndex = index-1    
        self.recipErrorfnamesCombo = QComboBox()
        self.recipErrorfnamesCombo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.recipErrorfnamesCombo.setMinimumWidth(150)
        self.recipErrorfnamesCombo.activated.connect(recipErrorfnamesComboFunc)

        def recipFilter(): # filter selected quad/elec or by reciprocal error
            pdebug('recipFilter()')    
            try:
                # compute index of the survey displayed (index == 0 if applyToEach)
                index = 0 if self.recipErrDataIndex < 0 else self.recipErrDataIndex
                numElecRemoved = np.sum(self.project.surveys[index].eselect)
                # don't need filterElec as the selected points are in iselect anyway
                numSelectRemoved = 0
                if self.project.iBatch or self.project.iTimeLapse:
                    if not self.recipErrApplyToAll:
                        numSelectRemoved += self.project.surveys[index].filterData(~self.project.surveys[index].iselect)
                        self.writeLog('# k.surveys[{:d}].filterData(i2keep)'
                                      ' # i2keep to be defined manually '.format(index))
                    else:
                        s = self.project.surveys[index] # index should be 0 as we select on the 1st survey only
                        quads = s.df[s.iselect][['a','b','m','n']].values
                        numSelectRemoved += self.project._filterSimilarQuad(quads) # this will remove the same quads to all surveys
                else:
                    numSelectRemoved += self.project.surveys[0].filterData(~self.project.surveys[0].iselect)
                if self.recipErrorInputLine.text() != '':
                    percent = float(self.recipErrorInputLine.text())
                    if self.filterAttrCombo.currentText() == 'Reciprocal Error':
                        numRecipRemoved = self.project.filterRecip(index=self.recipErrDataIndex, percent=percent)
                        self.writeLog('k.filterRecip(index={:d}, percent={:f})'.format(
                        self.recipErrDataIndex, percent))
                    elif self.filterAttrCombo.currentText() == 'Stacking Error (Dev.)':
                        numRecipRemoved = self.project.filterStack(index=self.recipErrDataIndex, percent=percent)
                        self.writeLog('k.filterStack(index={:d}, percent={:f})'.format(
                        self.recipErrDataIndex, percent))
                    if numElecRemoved != 0:
                        self.infoDump("%i measurements with greater than %3.1f%% %s, \
                                 %i selected electrodes and %i measurements removed!" % (numRecipRemoved,percent,self.filterAttrCombo.currentText(),numElecRemoved,numSelectRemoved))
                    else:
                        self.infoDump("%i measurements with greater than %3.1f%% %s and %i selected measurements removed!" % (numRecipRemoved,percent,
                                                                                                                              self.filterAttrCombo.currentText(),numSelectRemoved))
                else:
                    numRhoRangeRemoved = 0
                    rhoRangeText = ''
                    if self.rhoRangeMinInput.text() != '' or self.rhoRangeMaxInput.text() != '':
                        vmin = float(self.rhoRangeMinInput.text()) if self.rhoRangeMinInput.text() != '' else None
                        vmax = float(self.rhoRangeMaxInput.text()) if self.rhoRangeMaxInput.text() != '' else None
                        if self.filterAttrCombo.currentText() == 'Transfer Resistance':
                            numRhoRangeRemoved = self.project.filterTransferRes(index=self.recipErrDataIndex, vmin=vmin, vmax=vmax)
                            self.writeLog('k.filterTransferRes(index={:d}, vmin={}, vmax={})'.format(self.recipErrDataIndex, str(vmin), str(vmax)))
                        elif self.filterAttrCombo.currentText() == 'App. Resistivity':
                            numRhoRangeRemoved = self.project.filterAppResist(index=self.recipErrDataIndex, vmin=vmin, vmax=vmax)
                            self.writeLog('k.filterAppResist(index={:d}, vmin={}, vmax={})'.format(self.recipErrDataIndex, str(vmin), str(vmax)))
                        rhoRangeText = '%i measurements outside of the range and ' % numRhoRangeRemoved
                    if numElecRemoved != 0:
                        self.infoDump("%s%i selected electrodes and %i measurements removed!" % (rhoRangeText, numElecRemoved, numSelectRemoved))
                    else:
                        self.infoDump("%s%i selected measurements removed!" % (rhoRangeText, numSelectRemoved))
                if self.ipCheck.checkState() == Qt.Checked:
                    for s in self.project.surveys:
                        s.dfPhaseReset = s.df
                        s.filterDataIP = s.df
                    heatFilter()
                    self.iperrFitType.setCurrentIndex(0)
                    phaseplotError()
                self.errHist(self.recipErrDataIndex)
                self.plotManualFiltering(self.recipErrDataIndex)
                self.errFitType.setCurrentIndex(0)
                self.plotError()

            except ValueError as e:
                if self.ipCheck.checkState() != Qt.Checked:
                    self.errorDump(e)
                else:
                    self.errorDump('Index error! Reciprocal Filtering cannot be done after Phase Filtering.\n'
                              'Reset the filters and redo the filterings, first reciprocity then phase.')

        def resetRecipFilter():
            pdebug('resetRecipFilter()')
            numRestored = 0
            if self.recipErrApplyToAll:
                for s in self.project.surveys:
                    numRestored += len(s.dfReset) - len(s.df)
                    s.df = s.dfReset.copy()
            else:
                numRestored = len(self.project.surveys[self.recipErrDataIndex].dfReset) - len(self.project.surveys[self.recipErrDataIndex].df)
                self.project.surveys[self.recipErrDataIndex].df = self.project.surveys[self.recipErrDataIndex].dfReset.copy()
            if self.recipErrorInputLine.text() != '':
                self.errHist(self.recipErrDataIndex)
                self.recipErrorInputLine.setText('')
            if self.ipCheck.checkState() == Qt.Checked:
                if self.recipErrApplyToAll:
                    for s in self.project.surveys:
                        s.dfPhaseReset = s.dfReset.copy()
                        s.filterDataIP = s.dfReset.copy()
                else:
                    self.project.surveys[self.recipErrDataIndex].dfPhaseReset = self.project.surveys[self.recipErrDataIndex].dfReset.copy()
                    self.project.surveys[self.recipErrDataIndex].filterDataIP = self.project.surveys[self.recipErrDataIndex].dfReset.copy()
                heatFilter()
                self.iperrFitType.setCurrentIndex(0)
                phaseplotError()
                
            self.errHist(self.recipErrDataIndex)
            self.plotManualFiltering(self.recipErrDataIndex)
            self.errFitType.setCurrentIndex(0)
            self.plotError()
            self.recipErrorInputLine.setText('')
            self.rhoRangeMinInput.setText('')
            self.rhoRangeMaxInput.setText('')
            self.infoDump('%i measurements restored!' % numRestored)


        self.recipErrorInputLabel = QLabel('Threshold:')
        self.recipErrorInputLabel.hide()
        self.recipErrorInputLine = QLineEdit('')
        self.recipErrorInputLine.hide()
        self.recipErrorInputLine.setPlaceholderText('%')
        self.recipErrorInputLine.setFixedWidth(100)
        self.recipErrorInputLine.setValidator(QDoubleValidator())

        self.rhoRangeInputLabel = QLabel('Range:')
        
        self.rhoRangeMinInput = QLineEdit('')
        self.rhoRangeMinInput.setPlaceholderText('min')
        self.rhoRangeMinInput.setFixedWidth(80)
        self.rhoRangeMinInput.setValidator(QDoubleValidator())
        
        self.rhoRangeMaxInput = QLineEdit('')
        self.rhoRangeMaxInput.setPlaceholderText('max')
        self.rhoRangeMaxInput.setFixedWidth(80)
        self.rhoRangeMaxInput.setValidator(QDoubleValidator())
        
        
        def filterAttrComboFunc(index):
            self.plotManualFiltering(self.recipErrDataIndex)
            selection = self.filterAttrCombo.currentText()
            if selection in ['Transfer Resistance', 'App. Resistivity']:
                self.recipErrorInputLabel.hide()
                self.rhoRangeInputLabel.show()
                self.rhoRangeMinInput.show()
                self.rhoRangeMaxInput.show()
                self.recipErrorInputLine.setText('')
                self.rhoRangeMinInput.setText('')
                self.rhoRangeMaxInput.setText('')
                self.recipErrorInputLine.hide()
            else:
                self.recipErrorInputLabel.show()
                self.rhoRangeInputLabel.hide()
                self.rhoRangeMinInput.hide()
                self.rhoRangeMaxInput.hide()
                self.recipErrorInputLine.show()
                self.recipErrorInputLine.setText('')
                self.rhoRangeMinInput.setText('')
                self.rhoRangeMaxInput.setText('')

            
        self.filterAttrCombo = QComboBox()
        
        self.filterAttrComboItems = ['Transfer Resistance', 'App. Resistivity']
        self.filterAttrCombo.addItems(self.filterAttrComboItems)
        self.filterAttrCombo.activated.connect(filterAttrComboFunc)


        def recipErrorUnpairedFunc():
            pdebug('recipErrorUnpairedFunc()')
            index = -1 if self.recipErrDataIndex < 0 else self.recipErrDataIndex
            numRemoved = self.project.filterUnpaired(index=index)
            self.writeLog('k.filterUnpaired(index={:d})'.format(index))
            if self.ipCheck.checkState() == Qt.Checked:
                if self.recipErrApplyToAll:
                    for s in self.project.surveys:
                        s.dfPhaseReset = s.dfReset.copy()
                        s.filterDataIP = s.dfReset.copy()
                else:
                    self.project.surveys[self.recipErrDataIndex].dfPhaseReset = self.project.surveys[self.recipErrDataIndex].dfReset.copy()
                    self.project.surveys[self.recipErrDataIndex].filterDataIP = self.project.surveys[self.recipErrDataIndex].dfReset.copy()
                heatFilter()
                self.iperrFitType.setCurrentIndex(0)
                phaseplotError()
            self.errHist(self.recipErrDataIndex)
            self.plotManualFiltering(self.recipErrDataIndex)
            self.errFitType.setCurrentIndex(0)
            self.plotError()
            self.infoDump('%i unpaired quadrupoles removed!' % numRemoved)

        self.recipErrorUnpairedBtn = QPushButton('Remove Unpaired')
        self.recipErrorUnpairedBtn.setFixedWidth(150)
        self.recipErrorUnpairedBtn.setToolTip('Remove quadrupoles without reciprocals')
        self.recipErrorUnpairedBtn.clicked.connect(recipErrorUnpairedFunc)

        self.recipErrorPltBtn = QPushButton('Apply filters')
        self.recipErrorPltBtn.setToolTip('Removes measuremtns that have either greater reciprocal error than "Percent error threshold" or are manually selected or both!')
        self.recipErrorPltBtn.clicked.connect(recipFilter)
        self.recipErrorPltBtn.setFixedWidth(150)
        
        self.recipErrorResetBtn = QPushButton('Reset')
        self.recipErrorResetBtn.setStyleSheet("color: red")
        self.recipErrorResetBtn.setToolTip('This will restore all deleted measurements at this stage')
        self.recipErrorResetBtn.clicked.connect(resetRecipFilter)
        self.recipErrorResetBtn.setFixedWidth(150)
        
        def saveFilteredData():
            fname, savetyp = QFileDialog.getSaveFileName(self.tabImportingData,'Save Filtered Data', 
                                                         self.datadir, 'Res2DInv (*.dat);;Comma Separated Values (*.csv);;'
                                                         'E4D survey file (*.srv)') # can add Protocol (*.dat) and so on
            if fname != '':
                # elec = self.elecTable.getTable() # getting the topography info
                self.project.param['lineTitle'] = self.titleEdit.text()
                # if not (self.project.iBatch or self.project.iTimeLapse):
                    # spacing = float(self.elecDxEdit.text())
                # else:
                    # spacing = None
                # self.project.saveFilteredData(fname, elec, savetyp, spacing=spacing)
                self.project.saveFilteredData(fname, savetyp) # elec is already known to self.project because of self.updateElec()
                self.writeLog('k.saveFilteredData("{:s}", {:s})'.format(fname, savetyp))

        self.recipErrorSaveBtn = QPushButton('Save data')
        self.recipErrorSaveBtn.setStyleSheet("color: green")
        self.recipErrorSaveBtn.setToolTip('This will save the data in available formats (e.g., Res2DInv.dat)')
        self.recipErrorSaveBtn.clicked.connect(saveFilteredData)
        self.recipErrorSaveBtn.setFixedWidth(150)
       
        self.mwRecipError = MatplotlibWidget(navi=True, aspect='auto', itight=True)

        self.mwManualFiltering = MatplotlibWidget(navi=True, aspect='auto', itight=True)


        # layout
        self.recipErrorLayout = QVBoxLayout()

        self.recipErrorTopLayout = QVBoxLayout()
        self.recipErrorLayout.addLayout(self.recipErrorTopLayout, 0) # number is stretch factor

        self.recipErrorLabelLayout = QHBoxLayout()
        self.recipErrorLabelLayout.addWidget(self.recipErrorLabel, 1)
        self.recipErrorLabelLayout.addWidget(self.recipErrorfnamesComboLabel)
        self.recipErrorLabelLayout.addWidget(self.recipErrorfnamesCombo)
        self.recipErrorTopLayout.addLayout(self.recipErrorLabelLayout)
        
        self.recipErrorInputlayout = QHBoxLayout()
        self.recipErrorTopLayout.addLayout(self.recipErrorInputlayout)
        
        self.recipErrorInputLeftlayout = QHBoxLayout()
        self.recipErrorInputLeftlayout.setAlignment(Qt.AlignLeft)
        self.recipErrorInputLeftlayoutL = QHBoxLayout()
        self.recipErrorInputLeftlayoutL.setAlignment(Qt.AlignRight)
        self.recipErrorInputLeftlayoutL.addWidget(self.rhoRangeInputLabel)
        self.recipErrorInputLeftlayoutL.addWidget(self.recipErrorInputLabel)
        self.recipErrorInputLeftlayout.addLayout(self.recipErrorInputLeftlayoutL)
        
        self.recipErrorInputLineLayout = QHBoxLayout()
        self.recipErrorInputLineLayout.setAlignment(Qt.AlignLeft)
        self.recipErrorInputLineLayout.addWidget(self.rhoRangeMinInput)
        self.recipErrorInputLineLayout.addWidget(self.rhoRangeMaxInput)
        self.recipErrorInputLineLayout.addWidget(self.recipErrorInputLine)
        self.recipErrorInputLineLayout.addWidget(self.filterAttrCombo)
        self.recipErrorInputLeftlayout.addLayout(self.recipErrorInputLineLayout)

        self.recipErrorInputlayout.addLayout(self.recipErrorInputLeftlayout)

        self.recipErrorBtnLayout = QHBoxLayout()
        self.recipErrorBtnLayout.setAlignment(Qt.AlignRight)
        self.recipErrorBtnLayout.addWidget(self.recipErrorUnpairedBtn)
        self.recipErrorBtnLayout.addWidget(self.recipErrorPltBtn)
        self.recipErrorBtnLayout.addWidget(self.recipErrorResetBtn)
        self.recipErrorBtnLayout.addWidget(self.recipErrorSaveBtn)
        self.recipErrorInputlayout.addLayout(self.recipErrorBtnLayout, 1)
               
        #tab widgets for the graphs
        self.recipErrorBottomTabs = QTabWidget()

        self.recipErrorBottomLayout = QVBoxLayout()
        self.recipErrorPlotLayout = QVBoxLayout()
        self.recipErrorPlotLayout.addWidget(self.mwRecipError)

        self.recipErrorPseudoPlotLayout = QVBoxLayout()
        self.recipErrorPseudoPlotLayout.addWidget(self.mwManualFiltering)

        self.pseudoSectionPlotTab = QWidget()
        self.pseudoSectionPlotTab.setLayout(self.recipErrorPseudoPlotLayout)
        self.recipErrorBottomTabs.addTab(self.pseudoSectionPlotTab, 'Pseudo Section')

        self.errorHistogramPlotTab = QWidget()
        self.errorHistogramPlotTab.setLayout(self.recipErrorPlotLayout)
        self.recipErrorBottomTabs.addTab(self.errorHistogramPlotTab, 'Error Histogram')

        self.recipErrorBottomLayout.addWidget(self.recipErrorBottomTabs)
        self.recipErrorLayout.addLayout(self.recipErrorBottomLayout, 1) # '1' to keep the plot in largest strech

        self.recipErrorWidget .setLayout(self.recipErrorLayout)



#%% phase filtering tab
        self.ipfiltWidget = QWidget()
        self.tabPreProcessing.addTab(self.ipfiltWidget, 'Phase Filtering')
        self.tabPreProcessing.setTabEnabled(1, False)

        self.phasefiltlayout = QVBoxLayout()
        self.phaseLabelLayout = QHBoxLayout()
        
        self.phasefiltLabel = QLabel('<b>Filter the data based on the phase/IP measurements.</b><br>\
                                Below graphs show the status of filtered data versus raw input.')
        self.phasefiltLabel.setWordWrap(True)
        self.phasefiltLabel.setAlignment(Qt.AlignLeft)
        self.phaseLabelLayout.addWidget(self.phasefiltLabel, 1)
        
        self.phasefiltfnamesComboLabel = QLabel('Select a dataset:')
        self.phasefiltfnamesComboLabel.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
        self.phaseLabelLayout.addWidget(self.phasefiltfnamesComboLabel)
        
        def phasefiltfnamesComboFunc(index):
            if index == 0:
                self.phaseFiltDataIndex = -1
            elif index > 0:
                self.phaseFiltDataIndex = index-1
            if self.project.surveys != []:
                heatRaw()
                heatFilter()
        
        self.phasefiltfnamesCombo = QComboBox()
        self.phasefiltfnamesCombo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.phasefiltfnamesCombo.setMinimumWidth(150)
        self.phasefiltfnamesCombo.activated.connect(phasefiltfnamesComboFunc)
        self.phaseLabelLayout.addWidget(self.phasefiltfnamesCombo)
        self.phasefiltlayout.addLayout(self.phaseLabelLayout)

        def phirange():
            vmin = float(self.phivminEdit.text())
            vmax = float(self.phivmaxEdit.text())
            self.project.filterRangeIP(self.phaseFiltDataIndex, vmin, vmax)
            self.writeLog('k.filterRangeIP(index={:d}, vmin={:.2f}, vmax={:.2f})'.format(
                self.phaseFiltDataIndex, vmin, vmax))
            heatFilter()

        def removerecip():
            self.project.filterRecipIP(self.phaseFiltDataIndex)
            self.writeLog('k.filterRecipIP(index={:d})'.format(self.phaseFiltDataIndex))
            heatFilter()

        def removenested():
            self.project.filterNested(self.phaseFiltDataIndex)
            self.writeLog('k.filterNested(index={:d})'.format(self.phaseFiltDataIndex))
            heatFilter()

        def convFactK():
            if self.phaseFiltDataIndex == -1:
                for s in self.project.surveys:
                    s.kFactor = float(self.phiConvFactor.text())
            else:
                self.project.surveys[self.phaseFiltDataIndex].kFactor = float(self.phiConvFactor.text())
            heatFilter()
            heatRaw()

        self.phitoplayout = QHBoxLayout()
        self.phitoplayoutL = QHBoxLayout()
        self.phitoplayoutL.setAlignment(Qt.AlignLeft)
        self.phitoplayoutC = QHBoxLayout()
        self.phitoplayoutC.setAlignment(Qt.AlignRight)
        self.phitoplayoutR = QHBoxLayout()
        self.phitoplayoutR.setAlignment(Qt.AlignRight)
        self.phiConvFactorlabel = QLabel('Conversion factor k (φ = -kM):')
        self.phiConvFactorlabel.setToolTip('Assuming linear relationship.\nk = 1.2 is for IRIS Syscal devices\nThis equation is not used when importing phase data')
        self.phiConvFactor = QLineEdit()
        self.phiConvFactor.setFixedWidth(50)
        self.phiConvFactor.setValidator(QDoubleValidator())
        self.phiConvFactor.setText('1.2')
        self.phiConvFactor.setToolTip('Assuming linear relationship.\nk = 1.2 is for IRIS Syscal devices\nThis equation is not used when importing phase data')
        self.phiConvFactor.editingFinished.connect(convFactK)
        self.rangelabel = QLabel('Phase range filtering:')
        self.phivminlabel = QLabel('-φ<sub>min</sub>:')
        self.phivminEdit = QLineEdit()
        self.phivminEdit.setFixedWidth(50)
        self.phivminEdit.setValidator(QDoubleValidator())
        self.phivmaxlabel = QLabel('-φ<sub>max</sub>:')
        self.phivmaxEdit = QLineEdit()
        self.phivmaxEdit.setFixedWidth(50)
        self.phivmaxEdit.setValidator(QDoubleValidator())
        self.phivminEdit.setText('0')
        self.phivmaxEdit.setText('25')
        self.rangebutton = QPushButton('Apply')
        self.rangebutton.setFixedWidth(100)
        self.rangebutton.setAutoDefault(True)
        self.rangebutton.clicked.connect(phirange)

        self.recipFiltBtn = QPushButton('Remove reciprocals')
        self.recipFiltBtn.setFixedWidth(150)
        self.recipFiltBtn.setToolTip('Reciprocal measurements will not be considered for inversion in ResIPy.\nThis filter just visualize the removal')
        self.recipFiltBtn.setAutoDefault(True)
        self.recipFiltBtn.setEnabled(False)
        self.recipFiltBtn.clicked.connect(removerecip)

        self.nestedfilt = QPushButton('Remove nested')
        self.nestedfilt.setFixedWidth(150)
        self.nestedfilt.setToolTip('Measurments where M and/or N are inbetween A and B will be removed.\nNOTE: Wenner like arrays will also be affected')
        self.nestedfilt.setAutoDefault(True)
        self.nestedfilt.clicked.connect(removenested)

        self.phitoplayoutL.addWidget(self.phiConvFactorlabel)
        self.phitoplayoutL.addWidget(self.phiConvFactor)
        self.phitoplayoutC.addWidget(self.rangelabel)
        self.phitoplayoutR.addWidget(self.phivminlabel)
        self.phitoplayoutR.addWidget(self.phivminEdit)
        self.phitoplayoutR.addWidget(self.phivmaxlabel)
        self.phitoplayoutR.addWidget(self.phivmaxEdit)
        self.phitoplayoutR.addWidget(self.rangebutton)
        self.phitoplayoutR.addWidget(self.recipFiltBtn)
        self.phitoplayoutR.addWidget(self.nestedfilt)
        
        self.phitoplayout.addLayout(self.phitoplayoutL, 0)
        self.phitoplayout.addLayout(self.phitoplayoutC, 1)
        self.phitoplayout.addLayout(self.phitoplayoutR, 0)
        
        self.phasefiltlayout.addLayout(self.phitoplayout,0)

        def filt_reset():
            if self.phaseFiltDataIndex == -1:
                for s in self.project.surveys:
                    s.filterDataIP = s.dfPhaseReset.copy()
                    s.df = s.dfPhaseReset.copy()
                self.infoDump('Phase filters are now reset for all datasets!')
            else:
                self.project.surveys[self.phaseFiltDataIndex].filterDataIP = self.project.surveys[self.phaseFiltDataIndex].dfPhaseReset.copy()
                self.project.surveys[self.phaseFiltDataIndex].df = self.project.surveys[self.phaseFiltDataIndex].dfPhaseReset.copy()
                self.infoDump('Phase filters are now reset for selected dataset!')
            heatFilter()
            self.dcaProgress.setValue(0)
            

        def phiCbarRange():
            if self.phaseFiltDataIndex == -1:
                for s in self.project.surveys:
                    s.phiCbarmin = float(self.phiCbarminEdit.text())
                    s.phiCbarMax = float(self.phiCbarMaxEdit.text())
            else:
                self.project.surveys[self.phaseFiltDataIndex].phiCbarmin = float(self.phiCbarminEdit.text())
                self.project.surveys[self.phaseFiltDataIndex].phiCbarMax = float(self.phiCbarMaxEdit.text())
            heatFilter()
            heatRaw()

        def phiCbarDataRange():
            if self.phaseFiltDataIndex == -1:
                for s in self.project.surveys:
                    minDataIP = np.min(s.dfOrigin['ip'])
                    maxDataIP = np.max(s.dfOrigin['ip'])
                    if self.ftype == 'ProtocolIP':
                        s.phiCbarmin = -maxDataIP
                        s.phiCbarMax = -minDataIP
                    else:
                        s.phiCbarmin = minDataIP
                        s.phiCbarMax = maxDataIP
            else:
                minDataIP = np.min(self.project.surveys[self.phaseFiltDataIndex].dfOrigin['ip'])
                maxDataIP = np.max(self.project.surveys[self.phaseFiltDataIndex].dfOrigin['ip'])
                if self.ftype == 'ProtocolIP':
                    self.project.surveys[self.phaseFiltDataIndex].phiCbarmin = -maxDataIP
                    self.project.surveys[self.phaseFiltDataIndex].phiCbarMax = -minDataIP
                else:
                    self.project.surveys[self.phaseFiltDataIndex].phiCbarmin = minDataIP
                    self.project.surveys[self.phaseFiltDataIndex].phiCbarMax = maxDataIP
            heatFilter()
            heatRaw()

        self.resetlayout = QHBoxLayout()
        self.resetlayoutL = QHBoxLayout()
        self.resetlayoutL.setAlignment(Qt.AlignLeft)
        self.resetlayoutR = QHBoxLayout()
        self.resetlayoutR.setAlignment(Qt.AlignRight)
        
        self.phaseSavebtn = QPushButton('Save data')
        self.phaseSavebtn.setStyleSheet("color: green")
        self.phaseSavebtn.setToolTip('This will save the data in available formats (e.g. Res2DInv.dat)')
        self.phaseSavebtn.clicked.connect(saveFilteredData)
        self.phaseSavebtn.setFixedWidth(150)
        
        self.filtreset = QPushButton('Reset phase filters')
        self.filtreset.setStyleSheet("color: red")
        self.filtreset.setToolTip('Reset all the filtering.\nk factor is not affected')
        self.filtreset.setAutoDefault(True)
        self.filtreset.clicked.connect(filt_reset)
        self.filtreset.setFixedWidth(150)
        self.phiCbarminlabel = QLabel('Colorbar<sub>min</sub>: ')
        self.phiCbarminEdit = QLineEdit()
        self.phiCbarminEdit.setFixedWidth(50)
        self.phiCbarminEdit.setValidator(QDoubleValidator())
        self.phiCbarminEdit.setText('0')
        self.phiCbarMaxlabel = QLabel('Colorbar<sub>Max</sub>: ')
        self.phiCbarMaxEdit = QLineEdit()
        self.phiCbarMaxEdit.setFixedWidth(50)
        self.phiCbarMaxEdit.setValidator(QDoubleValidator())
        self.phiCbarMaxEdit.setText('25')
        self.phiCbarrangebutton = QPushButton('Apply')
        self.phiCbarrangebutton.setFixedWidth(100)
        self.phiCbarrangebutton.setToolTip('This is not a filtering step.')
        self.phiCbarrangebutton.setAutoDefault(True)
        self.phiCbarrangebutton.clicked.connect(phiCbarRange)
        self.phiCbarDatarangebutton = QPushButton('Raw data range')
        self.phiCbarDatarangebutton.setToolTip('This is not a filtering step.')
        self.phiCbarDatarangebutton.setAutoDefault(True)
        self.phiCbarDatarangebutton.clicked.connect(phiCbarDataRange)
        self.phiCbarDatarangebutton.setFixedWidth(150)
        self.resetlayoutL.addWidget(self.phiCbarminlabel)
        self.resetlayoutL.addWidget(self.phiCbarminEdit)
        self.resetlayoutL.addWidget(self.phiCbarMaxlabel)
        self.resetlayoutL.addWidget(self.phiCbarMaxEdit)
        self.resetlayoutL.addWidget(self.phiCbarrangebutton)
        self.resetlayoutL.addWidget(self.phiCbarDatarangebutton)
        self.resetlayoutR.addWidget(self.filtreset)
        self.resetlayoutR.addWidget(self.phaseSavebtn)
        
        self.resetlayout.addLayout(self.resetlayoutL, 0)
        self.resetlayout.addLayout(self.resetlayoutR, 1)
#        self.recipFiltBtn.clicked.connect("add function")


        self.ipfiltlayout = QHBoxLayout()

        def heatRaw():
            if self.phaseFiltDataIndex == -1:
                index = 0
            else:
                index = self.phaseFiltDataIndex
            self.project.surveys[index].filt_typ = 'Raw'
            raw_hmp.setCallback(self.project.showHeatmap)
            raw_hmp.replot(index=index)
            self.writeLog('k.showHeatmap()')

        def heatFilter():
            if self.phaseFiltDataIndex == -1:
                index = 0
            else:
                index = self.phaseFiltDataIndex
            self.project.surveys[index].filt_typ = 'Filtered'
            filt_hmp.setCallback(self.project.showHeatmap)
            filt_hmp.replot(index=index)
            self.writeLog('k.showHeatmap()')

        raw_hmp = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        filt_hmp = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        self.ipfiltlayout.addWidget(raw_hmp)
        self.ipfiltlayout.addWidget(filt_hmp)


        def dcaDump(val):
            self.dcaProgress.setValue(val)
            QApplication.processEvents()

        def dcaFiltering():
            try:
                self.dcaButton.setEnabled(False) # prevent multiple clicks!
                if len(self.project.surveys) > 1 and self.phaseFiltDataIndex == -1: # well, DCA is long!!
                    self.infoDump('DCA will run %s times. Please be patient!' % len(self.project.surveys))
                self.project.filterDCA(index=self.phaseFiltDataIndex, dump=dcaDump)
                self.writeLog('k.filterDCA(index={:d})'.format(self.phaseFiltDataIndex))
                heatFilter()
                self.dcaButton.setEnabled(True)
            except:
                self.errorDump('No decay curves found or incomplete set of decay curves! Export the data from "Prosys" with M1, M2, ... , M20 and TM1 tabs enabled.')
                self.dcaButton.setEnabled(True)

        self.dcaLayout = QHBoxLayout()
        self.dcaButton = QPushButton('DCA filtering')
        self.dcaButton.setToolTip('Decay Curve Analysis filtering.\nFor more see: Flores Orozco, et al. (2017), Decay curve analysis for data error quantification in\ntime-domain induced polarization imaging')
        self.dcaButton.setAutoDefault(True)
        self.dcaButton.clicked.connect(dcaFiltering)
        self.dcaProgress = QProgressBar()
        self.dcaButton.setEnabled(False)
        self.dcaProgress.setEnabled(False)
        self.dcaLayout.addWidget(self.dcaButton)
        self.dcaLayout.addWidget(self.dcaProgress)

        self.phasefiltlayout.addLayout(self.dcaLayout, 1)
        self.phasefiltlayout.addLayout(self.resetlayout, 2)
        self.phasefiltlayout.addLayout(self.ipfiltlayout, 3)


        # layout
        #TODO tidy up and put all layout here
        
        self.ipfiltWidget.setLayout(self.phasefiltlayout)


#%% resistance error modelling tab
        self.errorWidget = QWidget()
        self.tabPreProcessing.addTab(self.errorWidget, 'Resistance Error Model')
        self.tabPreProcessing.setTabEnabled(2, False)
        
        self.errFitLabel = QLabel('Select an error model from the drop-down menu. Once\
                             fitted, the model will generate an error for each quadrupoles\
                             (even the ones with no reciprocals). This error will\
                             be written in the <code>protocol.dat</code> file \
                             and used in the inversion if both <code>a_wgt</code> and\
                             <code>b_wgt</code> are both set to 0 (see \'Inversion settings\' tab).')
        self.errFitLabel.setWordWrap(True)
        self.errFitLabel.setToolTip('In case of batch/time-lapse inversion, <i>all</i> datesets must either have an error model \
                               or not have any error models (i.e., select separate error models for each individual dataset or "Apply to all"). \
                               ResIPy can handle batch data with mixture of different error models.')
        self.errFitLabel.setAlignment(Qt.AlignLeft)
        
        self.errFitfnamesComboLabel = QLabel('Select a dataset:')
        self.errFitfnamesComboLabel.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
        
        def errFitfnamesComboFunc(index):
            if index == 0: # fit on each, apply to each
                self.errFitDataIndex = -1
                self.plotError(0)
            elif index == 1: # fit on all combined, apply to each (bigSurvey)
                self.errFitDataIndex = -2
                self.plotError(-2)
            else:
                self.errFitDataIndex = index - 2
                self.plotError(index-2)
                self.errFitType.setCurrentIndex(self.errFitPlotIndexList[index-2])
                errFitTypeFunc(self.errFitPlotIndexList[index-2])
        self.errFitfnamesCombo = QComboBox()
        self.errFitfnamesCombo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.errFitfnamesCombo.setMinimumWidth(150)
        self.errFitfnamesCombo.activated.connect(errFitfnamesComboFunc)
        
        def errorModelSpecified():
            self.a_wgt.setText('0.0')
            a_wgtFunc()
            self.b_wgt.setText('0.0')
            b_wgtFunc()

        self.mwFitError = MatplotlibWidget(navi=True, aspect='auto', itight=True)

        def errFitTypeFunc(index):
            if len(self.project.surveys) == 0:
                return
            if index != 0:
                if self.errFitDataIndex == -1:
                    self.infoDump('Error model applied individually on all datasets')
                elif self.errFitDataIndex == -2:
                    self.infoDump('Error model fit on the combined datasets and then applied to all datasets.')
            if index == 0:
                if self.errFitDataIndex > 0:
                    self.plotError(self.errFitDataIndex)
                else:
                    self.plotError(0)
            elif index == 1:
                self.mwFitError.setCallback(self.project.fitErrorLin)
                self.mwFitError.replot(index=self.errFitDataIndex)
                self.writeLog('k.fitErrorLin(index={:d})'.format(self.errFitDataIndex))
                self.writeLog('k.err = True')
                self.project.err = True
            elif index == 2:
                self.mwFitError.setCallback(self.project.fitErrorPwl)
                self.mwFitError.replot(index=self.errFitDataIndex)
                self.writeLog('k.fitErrorPwl(index={:d})'.format(self.errFitDataIndex))
                self.writeLog('k.err = True')
                self.project.err = True
            elif index == 3:
                self.mwFitError.setCallback(self.project.fitErrorLME)
                self.mwFitError.replot(index=self.errFitDataIndex)
                self.writeLog('k.fitErrorLME(index={:d})'.format(self.errFitDataIndex))
                self.writeLog('k.err = True')
                self.project.err = True
                
            # record the type of fit for each survey
            if self.errFitDataIndex == -1: # same model for each
                self.errFitPlotIndexList = [index]*len(self.project.surveys)
            elif self.errFitDataIndex == -2: # same fit from bigSurvey apply on all
                pass
            elif self.errFitDataIndex >= 0:
                self.errFitPlotIndexList[self.errFitDataIndex] = index
            print('errFitPlotIndexList', self.errFitPlotIndexList)

            
            # if an error model is fitted we need to set a_wgt and b_wgt to 0
            if index == 0:
                self.a_wgt.setText('0.01')
                a_wgtFunc()
                self.b_wgt.setText('0.02')
                b_wgtFunc()
            else:
                self.a_wgt.setText('0.0')
                a_wgtFunc()
                self.b_wgt.setText('0.0')
                b_wgtFunc()
        self.errFitType = QComboBox()
        self.errFitType.addItem('Observed Errors')
        self.errFitType.addItem('Linear')
        self.errFitType.addItem('Power-law')
        if platform.system() == 'Linux':
            self.errFitType.addItem('Linear Mixed Effect (requires R and the lme4 package, dc surveys only for now)')
        self.errFitType.activated.connect(errFitTypeFunc)
        self.errFitType.setToolTip('Select an error model to use.')

        def saveErrBtnFunc():
            fname, _ = QFileDialog.getSaveFileName(self.tabImportingData,'Save error data file', self.datadir, 'Comma Separated Values (*.csv)')
            if fname != '':
                self.project.saveErrorData(fname)  
                self.writeLog('k.saveErrorData("{:s}")'.format(fname))
        self.saveErrBtn = QPushButton('Save Error Data')
        self.saveErrBtn.setStyleSheet("color: green")
        self.saveErrBtn.setFixedWidth(150)
        self.saveErrBtn.clicked.connect(saveErrBtnFunc)
        self.saveErrBtn.setToolTip('Save error data for DC and IP (if available) as .csv')
        
        # layout
        self.errorLayout = QVBoxLayout()
        self.errorLayout.setAlignment(Qt.AlignTop)
        
        self.errorTopLayout = QHBoxLayout()
        self.errorTopLayout.addWidget(self.errFitLabel, 1)
        self.errorTopLayout.addWidget(self.errFitfnamesComboLabel)
        self.errorTopLayout.addWidget(self.errFitfnamesCombo)
        self.errorLayout.addLayout(self.errorTopLayout)
        
        self.errFitLayout = QHBoxLayout()
        self.errFitLayout.addWidget(self.errFitType, 70)
        self.errFitLayout.addWidget(self.saveErrBtn, 30)
        self.errorLayout.addLayout(self.errFitLayout)
        
        self.errorPlotLayout = QVBoxLayout()
        self.errorPlotLayout.addWidget(self.mwFitError)
        self.errorLayout.addLayout(self.errorPlotLayout, 1)

        self.errorWidget.setLayout(self.errorLayout)
        

#%% IP error model tab
        self.ipWidget = QWidget()
        self.tabPreProcessing.addTab(self.ipWidget, 'Phase Error Model')
        self.tabPreProcessing.setTabEnabled(3, False)

        self.iperrFitLabel = QLabel('Select an error model from the drop-down menu. Once\
                     fitted, the model will generate an error for each quadrupoles\
                     (even the ones with no reciprocals). This error will\
                     be written in the <code>protocol.dat</code> file \
                     and used in the inversion if both <code>a_wgt</code> and\
                     <code>b_wgt</code> are both set to 0 (see \'Inversion settings\' tab).')
        self.iperrFitLabel.setWordWrap(True)
        self.iperrFitLabel.setToolTip('In case of batch/time-lapse inversion, <i>all</i> datesets must either have an error model \
                                 or not have any error models (i.e., select separate error models for each individual dataset or "Apply to all"). \
                                 ResIPy can handle batch data with mixture of different error models.')
        self.iperrFitLabel.setAlignment(Qt.AlignLeft)
        
        self.iperrFitfnamesComboLabel = QLabel('Select a dataset:')
        self.iperrFitfnamesComboLabel.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
        
        def iperrFitfnamesComboFunc(index):
            if index == 0:
                self.iperrFitDataIndex = -1 # fit each, apply each
                phaseplotError(0)
            elif index == 1:
                self.iperrFitDataIndex = -2 # fit combined, apply each
                phaseplotError(-2)
            elif index > 0: # show/hide makes the index = -1 so > 0 is necessary.
                self.iperrFitDataIndex = index-2
                phaseplotError(index-2)
                self.iperrFitType.setCurrentIndex(self.iperrFitPlotIndexList[index-2])
                iperrFitTypeFunc(self.iperrFitPlotIndexList[index-2])
        self.iperrFitfnamesCombo = QComboBox()
        self.iperrFitfnamesCombo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.iperrFitfnamesCombo.setMinimumWidth(150)
        self.iperrFitfnamesCombo.activated.connect(iperrFitfnamesComboFunc)
        
        self.mwIPFitError = MatplotlibWidget(navi=True, aspect='auto', itight=True)

        def phaseplotError(index=0):
            if len(self.project.surveys) == 0:
                return
            self.mwIPFitError.setCallback(self.project.showErrorIP)
            self.mwIPFitError.replot(index=index)
            self.project.err = False
            self.writeLog('k.showErrorIP(index={:d})'.format(index))
            self.writeLog('k.err = False')

        def iperrFitTypeFunc(index):
            if len(self.project.surveys) == 0:
                return
            if index != 0:
                if self.iperrFitDataIndex == -1:
                    self.infoDump('IP error model applied individually on all datasets')
                elif self.iperrFitDataIndex == -2:
                    self.infoDump('IP error model fit on the combined datasets and then applied to all datasets.')
            if index == 0:
                if self.iperrFitDataIndex > 0:
                    phaseplotError(self.iperrFitDataIndex)
                else:
                    phaseplotError(0)
            elif index == 1:
                self.mwIPFitError.setCallback(self.project.fitErrorPwlIP)
                self.mwIPFitError.replot(index=self.iperrFitDataIndex)
                self.project.err = True
                self.writeLog('k.fitErrorPwlIP(index={:d})'.format(self.iperrFitDataIndex))
                self.writeLog('k.err = True')
            elif index == 2:
                self.mwIPFitError.setCallback(self.project.fitErrorParabolaIP)
                self.mwIPFitError.replot(index=self.iperrFitDataIndex)
                self.project.err = True
                self.writeLog('k.fitErrorParabolaIP(index={:d})'.format(self.iperrFitDataIndex))
                self.writeLog('k.err = True')
                
            # record the type of fit for each survey
            if self.iperrFitDataIndex == -1: # same model for each
                self.iperrFitPlotIndexList = [index]*len(self.project.surveys)
            elif self.iperrFitDataIndex == -2: # same fit from bigSurvey apply on all
                pass
            elif self.iperrFitDataIndex >= 0:
                self.iperrFitPlotIndexList[self.iperrFitDataIndex] = index
            print('iperrFitPlotIndexList', self.iperrFitPlotIndexList)

                
            if index == 0:
                self.a_wgt.setText('0.02')
                a_wgtFunc()
                self.b_wgt.setText('2')
                b_wgtFunc()
            else:
                self.a_wgt.setText('0')
                a_wgtFunc()
                self.b_wgt.setText('0')
                b_wgtFunc()
        self.iperrFitType = QComboBox()
        self.iperrFitType.addItem('Observed discrepancies') 
        self.iperrFitType.addItem('Power law')
        self.iperrFitType.addItem('Parabola')
        self.iperrFitType.activated.connect(iperrFitTypeFunc)
        self.iperrFitType.setToolTip('Select an error model for IP.')
        
        self.saveIPErrBtn = QPushButton('Save Error Data')
        self.saveIPErrBtn.setStyleSheet("color: green")
        self.saveIPErrBtn.setFixedWidth(150)
        self.saveIPErrBtn.clicked.connect(saveErrBtnFunc)
        self.saveIPErrBtn.setToolTip('Save error data (DC and IP) as .csv')
        
        
        # layout
        self.ipLayout = QVBoxLayout()
        self.ipLayout.setAlignment(Qt.AlignTop)
        
        self.ipTopLayout = QHBoxLayout()
        self.ipTopLayout.addWidget(self.iperrFitLabel, 1)
        self.ipTopLayout.addWidget(self.iperrFitfnamesComboLabel)
        self.ipTopLayout.addWidget(self.iperrFitfnamesCombo)
        self.ipLayout.addLayout(self.ipTopLayout)

        self.errIPFitLayout = QHBoxLayout()
        self.errIPFitLayout.addWidget(self.iperrFitType, 70)
        self.errIPFitLayout.addWidget(self.saveIPErrBtn, 30)        
        self.ipLayout.addLayout(self.errIPFitLayout)
        
        self.ipErrPlotLayout = QVBoxLayout()
        self.ipErrPlotLayout.addWidget(self.mwIPFitError)
        self.ipLayout.addLayout(self.ipErrPlotLayout,1)

        self.ipWidget.setLayout(self.ipLayout)


#%% additional actions for pre-processing tab
        
        def errorCombosFill(comboboxes=[]): #filling pre-processing comboboxes with fnames
            for comboboxe in comboboxes:
                comboboxe.clear()
            
            for comboboxe in comboboxes:
                comboboxe.addItem('Apply to Each')
            
            for comboboxe in comboboxes[2:]: # only for error modelling
                comboboxe.addItem('Combined Fit')
            
            for s in self.project.surveys:
                for comboboxe in comboboxes:
                    comboboxe.addItem(s.name)
            
            for combo in comboboxes:
                combo.blockSignals(False)
                
#            self.phasefiltfnamesCombo.removeItem(0) #"Apply to all" won't work with phase filtering process.
                     
        self.errorCombosShow(False) #hiding all file selection comboboxes in pre-processing
        self.prepFnamesComboboxes = [self.recipErrorfnamesCombo, 
                                self.phasefiltfnamesCombo,
                                self.errFitfnamesCombo,
                                self.iperrFitfnamesCombo]



#%% ============================= mesh tab =======================
        self.tabMesh = QWidget()
        self.tabs.addTab(self.tabMesh, 'Mesh')
        self.tabs.setTabEnabled(2, False)
        
        self.fmdToolTip = 'Boundary depth (meters) below which the coarse (or background) mesh starts.\n'\
                     'Default = 2/3 max dipole length.'
        self.fmdLabel = QLabel('Fine/Coarse\nboundary depth')
        self.fmdLabel.setToolTip(self.fmdToolTip)
        self.fmdLabel.setAlignment(Qt.AlignCenter)
        self.fmdBox = QLineEdit()
        self.fmdBox.setAlignment(Qt.AlignCenter)
        self.fmdBox.setPlaceholderText('[m]')
        self.fmdBox.setToolTip(self.fmdToolTip)
        self.fmdBox.setValidator(QDoubleValidator())
        
        
        def replotMesh(aspect='equal'):
            if self.pseudo3DCheck.isChecked(): # for pseudo 3D inversion
                self.mesh3Dplotter.clear() # clear all actors 
                self.project.showPseudo3DMesh(ax=self.mesh3Dplotter, color_map='Greys', color_bar=False)
                self.meshOutputStack.setCurrentIndex(2)
            else: # 2D mesh
                def func(ax):
                    self.project.createModel(ax=ax, addAction=self.regionTable.addRow)
                self.mwMesh.plot(func, aspect=aspect)
                self.mwMesh.canvas.setFocusPolicy(Qt.ClickFocus) # allows the keypressevent to go to matplotlib
                self.mwMesh.canvas.setFocus() # set focus on the canvas
                self.meshOutputStack.setCurrentIndex(1)

        
        # design 2D features before meshing them, for instance for known
        # region of given resistivity. Note that the region must not overlap
        # and need to be distinct from each other (not touching each other)
        def designModel():
            self.iDesign = True
            # read electrodes locations
            elec = self.elecTable.getTable()
            if elec['x'].isna().sum() > 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            
            # plot the interactive model
            self.regionTable.reset()
            def func(ax):
                self.project.designModel(ax=ax, addAction=self.regionTable.addRow)
            self.mwMesh.plot(func, aspect = self.plotAspect)
            self.mwMesh.canvas.setFocusPolicy(Qt.ClickFocus) # allows the keypressevent to go to matplotlib
            self.mwMesh.canvas.setFocus() # set focus on the canvas
            # if self.iForward is False:
            #     self.regionTable.setColumnHidden(2, False) # show zone column
            #     self.regionTable.setColumnHidden(3, False) # show fixed column
            self.meshOutputStack.setCurrentIndex(1)

        self.designModelBtn = QPushButton('Design Model before meshing')
        self.designModelBtn.clicked.connect(designModel)



        # ------------ QUAD MESH -----------------
        # additional options for quadrilateral mesh
        self.nnodesLabel = QLabel('Number of elements\nbetween electrodes')
        self.nnodesSld = QSlider(Qt.Horizontal)
        self.nnodesSld.setMinimumWidth(50)
        self.nnodesSld.setMinimum(1)
        self.nnodesSld.setMaximum(10)
        self.nnodesSld.setValue(4)
        self.nnodesGrid = QGridLayout()
        # self.nnodesGrid.setContentsMargins(9,0,9,0)
        self.nnodesGrid.setSpacing(2)
        self.nnodes1Label = QLabel('1')
        self.nnodes10Label = QLabel('10')
        self.nnodes10Label.setAlignment(Qt.AlignRight)#| Qt.AlignVCenter)
        self.nnodesGrid.addWidget(self.nnodesSld, 0, 0, 1, 2)
        self.nnodesGrid.addWidget(self.nnodes1Label, 1,0,1,1)
        self.nnodesGrid.addWidget(self.nnodes10Label, 1,1,1,1)

        def meshQuadFunc():
            self.cropBelowFmd.setEnabled(True)
            self.regionTable.reset()
            elec = self.elecTable.getTable()
            if self.project.elec['x'].isna().sum() > 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            if self.topoTable.useNarray:
                surface = self.topoTable.xyz
            else:
                surface = self.topoTable.getTable()[['x','y','z']].values
            inan = ~np.isnan(surface[:,0])
            if np.sum(~inan) == surface.shape[0]:
                surface = None
            else:
                surface = surface[inan,:]
            
            nnodes = self.nnodesSld.value()
            try:
                fmd = np.abs(float(self.fmdBox.text())) if self.fmdBox.text() != '' else None
                if self.pseudo3DCheck.isChecked():
                    self.project.createMultiMesh(typ='quad', elemx=nnodes, surface=surface, fmd=fmd)
                    self.meshOutputStack.setCurrentIndex(2)
                    reqMemory = -1 if any([proj.param['reqMemory'] < 0 for proj in self.project.projs]) else 1
                else:
                    self.project.createMesh(typ='quad', elemx=nnodes, surface=surface, fmd=fmd)
                    self.writeLog('k.createMesh(typ="quad", elemx={:s}, surface={:s}, fmd={:s})'.format(
                        str(nnodes), str(surface), str(fmd)))
                    self.scale.setVisible(False)
                    self.scaleLabel.setVisible(False)
                    self.meshOutputStack.setCurrentIndex(1)
                    reqMemory = self.project.param['reqMemory']
                replotMesh()
                if reqMemory <= 0: # RAM requirement
                    self.errorDump('Make a coarser mesh!! It is likely that <b>more RAM is required</b> for inversion!')
                    self.ramRequiredLabel.show()
                else:
                    self.ramRequiredLabel.hide()
            except Exception as e:
                self.errorDump('Error creating the mesh: ' + str(e))
        self.meshQuad = QPushButton('Quadrilateral Mesh')
        self.meshQuad.setAutoDefault(True)
        self.meshQuad.clicked.connect(meshQuadFunc)
        self.meshQuad.setToolTip('Generate quadrilateral mesh.')


#%% --------------- TRIANGULAR MESH -------------------

        # additional options for triangular mesh
        self.clLabel = QLabel('Characteristic Length:')
        self.clLabel.setToolTip('Control the number and size of elements between electrodes.')
        self.clGrid = QGridLayout()
        self.clGrid.setContentsMargins(9,0,9,0)
        self.clGrid.setSpacing(2)
        self.clFineLabel = QLabel('Fine')
        self.clFineLabel.setStyleSheet('font:12px;')
        self.clCoarseLabel = QLabel('Coarse')
        self.clCoarseLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.clCoarseLabel.setStyleSheet('font:12px;')
        self.clSld = QSlider(Qt.Horizontal)
        self.clSld.setMinimum(1) # this depends on electrode spacing
        self.clSld.setMaximum(10)
        self.clSld.setValue(5)
        self.clGrid.addWidget(self.clSld, 0, 0, 1, 2)
        self.clGrid.addWidget(self.clFineLabel, 1,0,1,1)
        self.clGrid.addWidget(self.clCoarseLabel, 1,1,1,1)
        self.clFactorLabel = QLabel('Growth factor:')
        self.clFactorLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.clFactorLabel.setToolTip('Factor by which elements grow away from electrodes.')
        self.clFactorSld = QSlider(Qt.Horizontal)
        self.clFactorSld.setMaximumWidth(200)
        self.clFactorSld.setMinimum(1)
        self.clFactorSld.setMaximum(10)
        self.clFactorSld.setValue(4)
        self.refineTrianCheck = QCheckBox('Refine')
        self.refineTrianCheck.setToolTip('Refine the mesh for forward modelling'
                                    'without increasing the number of parameters.')

        def meshTrianFunc():
            if self.project.mproc is not None:
                print('killing')
                self.project.mproc.kill()
                self.project.mproc = None
                self.meshOutputStack.setCurrentIndex(1)
                self.meshTrianBtn.setText('Triangular Mesh')
                self.meshTrianBtn.setStyleSheet('background-color:orange; color:black')
                return
            else:
                self.meshTrianBtn.setText('Kill')
                self.meshTrianBtn.setStyleSheet('background-color:red; color:black')
            self.cropBelowFmd.setEnabled(True)
            elec = self.elecTable.getTable()
            if self.project.elec['x'].isna().sum() > 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            self.meshOutputStack.setCurrentIndex(0)
            QApplication.processEvents()
            self.meshLogText.clear()
            elecSpacing = np.sqrt(np.sum(np.diff(self.project.elec[~self.project.elec['remote']]['x'].values[:2])**2 +
                                         np.diff(self.project.elec[~self.project.elec['remote']]['z'].values[:2])**2))
            cl = float(self.clSld.value())/10*(elecSpacing-elecSpacing/8)
            cl_factor = self.clFactorSld.value()
            if self.topoTable.useNarray:
                surface = self.topoTable.xyz
            else:
                surface = self.topoTable.getTable()[['x','y','z']].values
            inan = ~np.isnan(surface[:,0])
            refine = 1 if self.refineTrianCheck.isChecked() else 0
            if np.sum(~inan) == surface.shape[0]:
                surface = None
            else:
                surface = surface[inan,:]
            fmd = np.abs(float(self.fmdBox.text())) if self.fmdBox.text() != '' else None
            pdebug('meshTrian(): fmd', fmd)
            pdebug('meshTrian(): elec:', self.project.elec)
            pdebug('meshTrian(): surface:', surface)
            self.writeLog('k.createMesh(typ="trian", geom_input={:s}, surface={:s}, cl={:f},'
                          ' cl_factor={:f}, refine={:s}, fmd={:s})'.format(
                              str(self.project.geom_input), str(surface), cl,
                              cl_factor, str(refine), str(fmd)))
            try:
                if self.pseudo3DCheck.isChecked():
                    self.loadingWidget('Creating the mesh. Please wait...')
                    self.project.createMultiMesh(typ='trian', surface=surface, cl=cl, cl_factor=cl_factor, 
                                                 show_output=True, dump=meshLogTextFunc, refine=refine, fmd=fmd)
                    self.loadingWidget(exitflag=True)
                    self.meshOutputStack.setCurrentIndex(2)
                    reqMemory = -1 if any([proj.param['reqMemory'] < 0 for proj in self.project.projs]) else 1
                else:
                    self.project.createModelMesh(surface=surface,
                                            cl=cl, cl_factor=cl_factor, show_output=True,
                                            dump=meshLogTextFunc, refine=refine, fmd=fmd)
                    # if self.iForward is False:
                    #     self.regionTable.setColumnHidden(2, False) # show zone column
                    #     self.regionTable.setColumnHidden(3, False) # show fixed column
                    self.scale.setVisible(True)
                    self.scaleLabel.setVisible(True)
                    self.meshOutputStack.setCurrentIndex(1)
                    reqMemory = self.project.param['reqMemory']
                if reqMemory <= 0: # RAM requirement
                    self.errorDump('Make a coarser mesh!! It is likely that <b>more RAM is required</b> for inversion!')
                    self.ramRequiredLabel.show()
                else:
                    self.ramRequiredLabel.hide()
            except Exception:
                pass # caused by killing the mesh process
            self.meshTrianBtn.setText('Triangular Mesh')
            self.meshTrianBtn.setStyleSheet('background-color:orange; color:black')
            replotMesh()
        self.meshTrianBtn = QPushButton('Triangular Mesh')
        self.meshTrianBtn.setAutoDefault(True)
        self.meshTrianBtn.setFocus()
        self.meshTrianBtn.clicked.connect(meshTrianFunc)
        self.meshTrianBtn.setToolTip('Generate triangular mesh.')
        self.meshTrianBtn.setStyleSheet('background-color:orange; color:black')



#%% ------------------------ TETRAHEDRAL MESH --------------------

        # additional options for tetrahedral mesh
        self.cl3ToolTip = 'Describes how big the nodes assocaited elements will be aroud the electrodes.\n' \
                     'Default: 1/2 the minimum electrode spacing (in meters).'
        self.cl3Label = QLabel('Characteristic Length:')
        self.cl3Label.setToolTip(self.cl3ToolTip)
        self.cl3Edit = QLineEdit()
        self.cl3Edit.setValidator(QDoubleValidator())
        self.cl3Edit.setText('')
        self.cl3Edit.setPlaceholderText('[m]')
        self.cl3Edit.setToolTip(self.cl3ToolTip)
        self.cl3FactorToolTip = 'Factor for the incremental size increase with depth in the mesh.\n' \
                           'Default: 8 (elements at the fine/coarse bondary depth are 8 times as big as those at the surface)'
        self.cl3FactorLabel = QLabel('Growth factor Top:')
        self.cl3FactorLabel.setToolTip(self.cl3FactorToolTip)
        self.cl3FactorEdit = QLineEdit()
        self.cl3FactorEdit.setToolTip(self.cl3FactorToolTip)
        self.cl3FactorEdit.setValidator(QDoubleValidator())
        self.cl3FactorEdit.setText('8')
        self.clnFactorToolTip = 'Factor applied to the characteristic length for the region below fine/coarse bondary depth to ' \
                           'compute a characteristic length for background region'
        self.clnFactorLabel = QLabel('Growth factor Bottom:')
        self.clnFactorLabel.setToolTip(self.clnFactorToolTip)
        self.clnFactorEdit = QLineEdit()
        self.clnFactorEdit.setToolTip(self.clnFactorToolTip)
        self.clnFactorEdit.setValidator(QDoubleValidator())
        self.clnFactorEdit.setText('100')
        self.refineTetraCheck = QCheckBox('Refine')
        self.refineTetraCheck.setToolTip('Refine the mesh for forward modelling'
                                    'without increasing the number of parameters.')
        
        def meshTetraFunc():
            self.regionTable.reset()
            if self.project.mproc is not None:
                print('killing')
                self.project.mproc.kill()
                self.project.mproc = None
                self.meshOutputStack.setCurrentIndex(1)
                self.mesh3DBtn.setText('Create Mesh')
                self.mesh3DBtn.setStyleSheet('background-color:orange; color:black')
                return
            else:
                self.mesh3DBtn.setText('Kill')
                self.mesh3DBtn.setStyleSheet('background-color:red; color:black')
            elec = self.elecTable.getTable()
            if self.topoTable.useNarray:
                topo = self.topoTable.xyz
            else:
                topo = self.topoTable.getTable()[['x','y','z']].values
            inan = ~np.isnan(topo[:,0])
            
            if self.project.elec['x'].isna().sum() > 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            elif all(elec['y'].values == 0) & all(topo[inan,1] == 0):
                self.errorDump('For 3D meshes, Y coordinates must be supplied for topo or elec at least.')
                return
            self.meshOutputStack.setCurrentIndex(0)
            QApplication.processEvents()
            self.meshLogText.clear()
            cl = -1 if self.cl3Edit.text() == '' else float(self.cl3Edit.text())
            cl_factor = float(self.cl3FactorEdit.text())
            cln_factor = float(self.clnFactorEdit.text()) if self.clnFactorEdit.text() != '' else 100
            refine = 1 if self.refineTetraCheck.isChecked() else 0
            if np.sum(~inan) == topo.shape[0]:
                topo = None
            else:
                topo = topo[inan,:]
            
            fmd = np.abs(float(self.fmdBox.text())) if self.fmdBox.text() != '' else None
            self.project.createMesh(typ='tetra', surface=topo, fmd=fmd,
                               cl=cl, cl_factor=cl_factor, dump=meshLogTextFunc,
                               cln_factor=cln_factor, refine=refine, show_output=True)
            self.writeLog('k.createMesh(typ="tetra", surface={:s}, fmd={:s}, cl={:.2f},'
                          ' cl_factor={:.2f}, cln_factor={:.2f}, refine={:d})'.format(
                              str(topo), str(fmd), cl, cl_factor, cln_factor, refine))
            if pvfound:
                self.mesh3Dplotter.clear() # clear all actors 
                self.project.showMesh(ax=self.mesh3Dplotter, color_map='Greys', color_bar=False)
            else:
                self.mwMesh3D.plot(self.project.showMesh, threed=True)
            self.writeLog('k.showMesh()')
            self.meshOutputStack.setCurrentIndex(2)
            if self.project.param['reqMemory'] <= 0: # RAM requirement
                self.errorDump('Make a coarser mesh!! It is likely that <b>more RAM is required</b> for inversion!')
                self.ramRequiredLabel.show()
            else:
                self.ramRequiredLabel.hide()
            self.mesh3DBtn.setText('Create Mesh')
            self.mesh3DBtn.setStyleSheet('background-color:orange; color:black')

        # self.meshTetraBtn = QPushButton('Tetrahedral Mesh')
        # self.meshTetraBtn.setAutoDefault(True)
        # self.meshTetraBtn.clicked.connect(meshTetraFunc)
        # self.meshTetraBtn.setToolTip('Generate tetrahedral mesh.')
        # self.meshTetraBtn.setStyleSheet('background-color:orange; color:black')


#%% ------------------------- TANK MESH -----------------------------
        # addition options for tank mesh
        self.clTankLabel = QLabel('Characteristic Length:')
        self.clTankEdit = QLineEdit()
        self.clTankEdit.setToolTip('Characterist Length for electrodes')
        self.originTankLabel = QLabel('Origin [X,Y,Z]:')
        self.originTankEdit = QLineEdit()
        self.originTankEdit.setToolTip('Origin of the mesh as X,Y,Z')
        self.dimensionTankLabel = QLabel('Dimensions [X,Y,Z]:')
        self.dimensionTankEdit = QLineEdit()
        self.dimensionTankEdit.setToolTip('Dimension of the mesh in X,Y,Z from origin')
        # self.refineTankLabel = QLabel('Refine:')
        self.refineTankCheck = QCheckBox('Refine')
        self.refineTankCheck.setToolTip('Refin mesh by splitting parameters into 6 elements')
        
        def tankDefaultText():
            try:
                elec = self.project.elec
                minvals = elec[['x','y','z']].min().values
                maxvals = elec[['x','y','z']].max().values
                ranvals = np.abs(maxvals - minvals)
                meanval = np.mean(ranvals)
                ranvals[ranvals == 0] = meanval
                origvals = minvals.copy()
                origvals[2] = minvals[2] - 0.1 * ranvals[2] # only changing z (or one side) is enough
                dimvals = ranvals.copy()
                dimvals[2] = ranvals[2] + 0.2 * ranvals[2] # dimension changes must be twice as origin displacement
                self.originTankEdit.setText(','.join(origvals.astype(str)))
                self.dimensionTankEdit.setText(','.join(dimvals.astype(str)))
            except:
                pass # probably the electrodes are not defined yet

        def meshTankFunc():
            # cosmetic update to say we are running
            self.regionTable.reset()
            if self.project.mproc is not None:
                print('killing')
                self.project.mproc.kill()
                self.project.mproc = None
                self.meshOutputStack.setCurrentIndex(1)
                self.mesh3DBtn.setText('Tank Mesh')
                self.mesh3DBtn.setStyleSheet('background-color:orange; color:black')
                return
            else:
                self.mesh3DBtn.setText('Kill')
                self.mesh3DBtn.setStyleSheet('background-color:red; color:black')

            # check electrodes are present
            elec = self.elecTable.getTable()
            if self.project.elec['x'].isna().sum() > 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            
            # retrieve parameters
            cl = -1 if self.clTankEdit.text() == '' else float(self.clTankEdit.text())
            if self.originTankEdit.text() != '':
                origin = [float(a) for a in self.originTankEdit.text().split(',')]
            else:
                origin = np.min(elec, axis=0)
            if self.dimensionTankEdit.text() != '':
                dimension = [float(a) for a in self.dimensionTankEdit.text().split(',')]
            else:
                dimension = np.max(elec, axis=0) - np.min(elec, axis=0)
            refine = 1 if self.refineTankCheck.isChecked() else 0

            # run the meshing
            self.meshOutputStack.setCurrentIndex(0)
            QApplication.processEvents()
            self.meshLogText.clear()
            self.project.createMesh(typ='tank',
                                    cl=cl,
                                    origin=origin,
                                    dimension=dimension,
                                    dump=meshLogTextFunc,
                                    refine=refine,
                                    show_output=True)
            self.writeLog('k.createMesh(typ="tank", cl={:.2e}, origin={:s}, dimension={:s},'
                          'refine={:d})'.format(cl, str(origin), str(dimension), refine))
            
            # display the mesh
            if pvfound:
                self.mesh3Dplotter.clear() # clear all actors 
                self.project.showMesh(ax=self.mesh3Dplotter, color_map='Greys', color_bar=False)
            else:
                self.mwMesh3D.plot(self.project.showMesh, threed=True)
            self.writeLog('k.showMesh()')
            self.meshOutputStack.setCurrentIndex(2)
            if self.project.param['reqMemory'] <= 0: # RAM requirement
                self.errorDump('Make a coarser mesh!! It is likely that <b>more RAM is required</b> for inversion!')
                self.ramRequiredLabel.show()
            else:
                self.ramRequiredLabel.hide()
            self.mesh3DBtn.setText('Tank Mesh')
            self.mesh3DBtn.setStyleSheet('background-color:orange; color:black')

        # self.meshTankBtn = QPushButton('Tank Mesh')
        # self.meshTankBtn.setAutoDefault(True)
        # self.meshTankBtn.clicked.connect(meshTankFunc)
        # self.meshTankBtn.setToolTip('Generate a tank mesh.')
        # self.meshTankBtn.setStyleSheet('background-color:orange; color:black')


#%% ------------------------- CYLINDER MESH -----------------------------
        # addition options for cylinder mesh
        self.clCylinderLabel = QLabel('Characteristic Length')
        self.clCylinderEdit = QLineEdit()
        self.clCylinderEdit.setToolTip('Characteristic Length for electrodes')
        # self.originCylinderEdit = QLineEdit()
        # self.originCylinderEdit.setToolTip('Origin of the mesh as X,Y,Z')
        # self.dimensionCylinderEdit = QLineEdit()
        # self.dimensionCylinderEdit.setToolTip('Dimension of the mesh in X,Y,Z from origin')
        # self.refineCylinderLabel = QLabel('Refine:')
        self.refineCylinderCheck = QCheckBox('Refine')
        self.refineCylinderCheck.setToolTip('Refine mesh by splitting parameters into 6 elements')
        self.zlimCylinderLabel = QLabel('zBottom, zTop [m]:')
        self.zlimCylinderEdit = QLineEdit()
        
        def cylinderEditDefaultText():
            try:
                elec = self.elecTable.getTable()
                zBottom = elec['z'].min()
                zTop = elec['z'].max()
                distance = 0.1*(np.abs(zTop - zBottom))
                if not np.isnan(zBottom) and not np.isnan(zTop):
                    self.zlimCylinderEdit.setText(str(zBottom-distance)+','+str(zTop+distance))
            except:
                pass # probably the electrodes are not defined yet

        def meshCylinderFunc():
            # cosmetic update to say we are running
            self.regionTable.reset()
            if self.project.mproc is not None:
                print('killing')
                self.project.mproc.kill()
                self.project.mproc = None
                self.meshOutputStack.setCurrentIndex(1)
                self.mesh3DBtn.setText('Cylinder Mesh')
                self.mesh3DBtn.setStyleSheet('background-color:orange; color:black')
                return
            else:
                self.mesh3DBtn.setText('Kill')
                self.mesh3DBtn.setStyleSheet('background-color:red; color:black')

            # check electrodes are present
            elec = self.elecTable.getTable()
            if self.project.elec['x'].isna().sum() > 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            
            # retrieve parameters
            cl = -1 if self.clCylinderEdit.text() == '' else float(self.clCylinderEdit.text())
            # if self.originCylinderEdit.text() != '':
            #     origin = [float(a) for a in self.originCylinderEdit.text().split(',')]
            # else:
            #     origin = np.min(elec, axis=0)
            # if self.dimensionCylinderEdit.text() != '':
            #     dimension = [float(a) for a in self.dimensionCylinderEdit.text().split(',')]
            # else:
            #     dimension = np.max(elec, axis=0) - np.min(elec, axis=0)
            # print('++++', origin, dimension)
            zlim = [float(a) for a in self.zlimCylinderEdit.text().split(',')]
            refine = 1 if self.refineCylinderCheck.isChecked() else 0

            # run the meshing
            self.meshOutputStack.setCurrentIndex(0)
            QApplication.processEvents()
            self.meshLogText.clear()
            self.project.createMesh(typ='cylinder',
                                    cl=cl,
                                    # origin=origin,
                                    # dimension=dimension,
                                    zlim=zlim,
                                    dump=meshLogTextFunc,
                                    refine=refine,
                                    show_output=True)
            self.writeLog('k.createMesh(typ="cylinder", cl={:.2e}, zlim={:s}'
                          ', refine={:d})'.format(cl, str(zlim), refine))
            
            # display the mesh
            if pvfound:
                self.mesh3Dplotter.clear() # clear all actors 
                self.project.showMesh(ax=self.mesh3Dplotter, color_map='Greys', color_bar=False)
            else:
                self.mwMesh3D.plot(self.project.showMesh, threed=True)
            self.writeLog('k.showMesh()')
            self.meshOutputStack.setCurrentIndex(2)
            if self.project.param['reqMemory'] <= 0: # RAM requirement
                self.errorDump('Make a coarser mesh!! It is likely that <b>more RAM is required</b> for inversion!')
                self.ramRequiredLabel.show()
            else:
                self.ramRequiredLabel.hide()
            self.mesh3DBtn.setText('Cylinder Mesh')
            self.mesh3DBtn.setStyleSheet('background-color:orange; color:black')

        # self.meshCylinderBtn = QPushButton('Cylinder Mesh')
        # self.meshCylinderBtn.setAutoDefault(True)
        # self.meshCylinderBtn.clicked.connect(meshCylinderFunc)
        # self.meshCylinderBtn.setToolTip('Generate a Cylinder mesh.')
        # self.meshCylinderBtn.setStyleSheet('background-color:orange; color:black')


        
#%% -------------------- ADDITIONAL 

        # general meshing button (label and connect change according to combo)
        self.mesh3DBtn = QPushButton('Create Mesh')
        self.mesh3DBtn.setAutoDefault(True)
        self.mesh3DBtn.clicked.connect(meshTetraFunc)
        self.mesh3DBtn.setToolTip('Generate tetrahedral mesh.')
        self.mesh3DBtn.setStyleSheet('background-color:orange; color:black')
        

        def mesh3DComboFunc(index):
            if index != 3:
                self.mesh3DBtn.setText('Create Mesh')
                self.mesh3DBtn.setToolTip('Generate tetrahedral mesh.')
            if index == 0:
                self.fmdGroup.setVisible(True)
                self.meshTetraGroup.setVisible(True)
                self.meshTankGroup.setVisible(False)
                self.meshCylinderGroup.setVisible(False)
                self.meshCustom3dGroup.setVisible(False)
                self.mesh3DBtn.disconnect()
                self.mesh3DBtn.clicked.connect(meshTetraFunc)
            elif index == 1:
                self.fmdGroup.setVisible(False)
                self.meshTetraGroup.setVisible(False)
                self.meshTankGroup.setVisible(True)
                self.meshCylinderGroup.setVisible(False)
                self.meshCustom3dGroup.setVisible(False)
                if self.originTankEdit.text() == '' and self.dimensionTankEdit.text() == '':
                    tankDefaultText()
                self.mesh3DBtn.disconnect()
                self.mesh3DBtn.clicked.connect(meshTankFunc)
                
            elif index == 2:
                self.fmdGroup.setVisible(False)
                self.meshTetraGroup.setVisible(False)
                self.meshTankGroup.setVisible(False)
                self.meshCylinderGroup.setVisible(True)
                self.meshCustom3dGroup.setVisible(False)
                if self.zlimCylinderEdit.text() == '':
                    cylinderEditDefaultText()
                self.mesh3DBtn.disconnect()
                self.mesh3DBtn.clicked.connect(meshCylinderFunc)
            elif index == 3:
                self.fmdGroup.setVisible(False)
                self.meshTetraGroup.setVisible(False)
                self.meshTankGroup.setVisible(False)
                self.meshCylinderGroup.setVisible(False)
                self.meshCustom3dGroup.setVisible(True)
                self.mesh3DBtn.setText('Import Mesh')
                self.mesh3DBtn.setToolTip('Import custom mesh.')
                self.mesh3DBtn.disconnect()
                self.mesh3DBtn.clicked.connect(importCustomMeshFunc)
            self.mesh3DBtn.setStyleSheet('background-color:orange; color:black')
        self.mesh3DCombo = QComboBox()
        self.mesh3DCombo.addItem('Half-space (tetra)')
        self.mesh3DCombo.addItem('Tank (box)')
        self.mesh3DCombo.addItem('Cylinder (column)')
        self.mesh3DCombo.addItem('Custom Mesh')
        self.mesh3DCombo.currentIndexChanged.connect(mesh3DComboFunc)

        def saveMeshBtnFunc():
            fname, _ = QFileDialog.getSaveFileName(self.tabMesh, 'Save mesh as .vtk, .dat or .node (tetgen) format', self.datadir)
            if fname != '':
                self.project.saveMesh(fname)
                self.writeLog('k.saveMesh("{:s}")'.format(fname))
                self.infoDump('Mesh saved to {:s}'.format(fname))
        self.saveMeshBtn = QPushButton('Save Mesh')
        self.saveMeshBtn.setToolTip('Save mesh as .vtk, .node (tetgen) or .dat')
        self.saveMeshBtn.clicked.connect(saveMeshBtnFunc)

        self.importCustomMeshLabel = QLabel('Import .msh or .vtk or .dat files.<br>'
                                            'The electrodes will be snapped to the closest node.')
        self.importCustomMeshLabel.setAlignment(Qt.AlignLeft)
        self.importCustomMeshLabel.setWordWrap(True)
        
        def importCustomMeshFunc():
            elec = self.elecTable.getTable()
            self.regionTable.reset()
            if self.project.elec['x'].isna().sum() > 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingData,'Open File', self.datadir)
            if fname != '':
                try:
                    self.project.importMesh(fname)
                    self.writeLog('k.importMesh("{:s}")'.format(fname))
                    print('mesh imported ... now displaying ... ')
                    if pvfound:
                        self.mesh3Dplotter.clear() # clear all actors 
                        if len(np.unique(self.project.mesh.df['region'])) == 1:
                            color_map = 'Greys'
                            color_bar = False
                        else:
                            color_map = 'Spectral'
                            color_bar = True
                        self.project.showMesh(ax=self.mesh3Dplotter, color_map=color_map, color_bar=color_bar)
                    else:
                        self.mwMesh3D.plot(self.project.showMesh, threed=True)
                    self.writeLog('k.showMesh()')
                    self.meshOutputStack.setCurrentIndex(2)
                    self.regionTable.nrow = 0
                    self.regionTable.setRowCount(0)
                    for i in range(len(np.unique(self.project.mesh.df['region']))):
                        self.regionTable.addRow()
                        ie = self.project.mesh.df['region'] == i+1
                        row = self.project.mesh.df[ie]
                        self.regionTable.setItem(i, 0, QTableWidgetItem(str(row['res0'].values[0])))
                        self.regionTable.setItem(i, 1, QTableWidgetItem(str(row['phase0'].values[0])))
                        self.regionTable.setItem(i, 2, QTableWidgetItem(str(row['zones'].values[0])))
                        if row['param'].values[0] == 0:
                            self.regionTable.cellWidget(i,3).findChildren(QCheckBox)[0].setChecked(True)
                except Exception as e:
                    self.errorDump('Error importing mesh' + str(e))

        self.importCustomMeshLabel2 = QLabel('Import .msh or .vtk or .dat files.')
        self.importCustomMeshLabel2.setAlignment(Qt.AlignCenter)
        self.importCustomMeshLabel2.setWordWrap(True)
        
        def importCustomMeshFunc2():
            print('using importCustomMeshFunc2')
            self.regionTable.reset()
            elec = self.elecTable.getTable()
            if self.project.elec['x'].isna().sum() > 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingData,'Open File', self.datadir)
            if fname != '':
                try:
                    self.project.importMesh(fname, mesh_type='trian')
                    if (self.project.typ == 'R3t') or (self.project.typ == 'cR3t'):
                        self.mwMesh.plot(self.project.showMesh, threed=True)
                        self.meshOutputStack.setCurrentIndex(2)
                        print('NO WAY THIS CAN HAPPEN!')
                    else:
                        replotMesh()
                        self.meshOutputStack.setCurrentIndex(1)
                    self.regionTable.nrow = 0
                    self.regionTable.setRowCount(0)
                    for i in range(len(np.unique(self.project.mesh.df['region']))):
                        self.regionTable.addRow()
                        ie = self.project.mesh.df['region'] == i+1
                        row = self.project.mesh.df[ie]
                        self.regionTable.setItem(i, 0, QTableWidgetItem(str(row['res0'].values[0])))
                        self.regionTable.setItem(i, 1, QTableWidgetItem(str(row['phase0'].values[0])))
                        self.regionTable.setItem(i, 2, QTableWidgetItem(str(row['zones'].values[0])))
                        if row['param'].values[0] == 0:
                            self.regionTable.cellWidget(i,3).findChildren(QCheckBox)[0].setChecked(True)
                except Exception as e:
                    self.errorDump('Error importing mesh' + str(e))
        self.importCustomMeshBtn2 = QPushButton('Import Custom Mesh')
        self.importCustomMeshBtn2.setFixedWidth(150)
        self.importCustomMeshBtn2.clicked.connect(importCustomMeshFunc2)
        self.importCustomMeshBtn2.setToolTip('Import .msh or .vtk file. The electrodes will be snapped to the closest node.')

        def meshAspectBtnFunc(state):
            if state == Qt.Checked:
                replotMesh(aspect='equal')
            else:
                replotMesh(aspect='auto')
        self.meshAspectBtn = QCheckBox('Equal aspect')
        self.meshAspectBtn.setToolTip('Check for equal aspect of axis.')
        self.meshAspectBtn.setChecked(True)
        self.meshAspectBtn.stateChanged.connect(meshAspectBtnFunc)

        def select3DRegionBtnFunc():
            self.mesh3Dplotter.clear() # clear all actors 
            self.clip = self.project.mesh.pick3Dbox(ax=self.mesh3Dplotter, darkMode=eval(resipySettings.param['dark'])) #extracts the surface and plots transparent boxed mesh
            self.select3DRegionBtn.setDisabled(True) # so button can't be clicked again 
            self.add3DRegionBtn.setDisabled(False)
            self.fin3DRegionBtn.setDisabled(False)
        self.select3DRegionBtn = QPushButton('(1) Region Mode')
        self.select3DRegionBtn.clicked.connect(select3DRegionBtnFunc)
        self.select3DRegionBtn.setVisible(False)

        def add3DRegionBtnFunc():
            clipped_mesh = self.project.mesh.addRegion3D(self.clip)
            self.regionTable.addRow()
            # self.mesh3Dplotter.clear() 
            clipped_mesh.show(ax=self.mesh3Dplotter, attr='region', color_bar=True, darkMode=eval(resipySettings.param['dark']))
        self.add3DRegionBtn = QPushButton('(2) Add region')
        self.add3DRegionBtn.clicked.connect(add3DRegionBtnFunc)
        self.add3DRegionBtn.setVisible(False)
        self.add3DRegionBtn.setDisabled(True)#start disabled 
        
        def fin3DRegionBtnFunc():#finish region selection mode 
            self.mesh3Dplotter.clear() # clear all actors 
            self.project.showMesh(ax=self.mesh3Dplotter, color_map='Spectral', attr='region', color_bar=True)
            self.select3DRegionBtn.setDisabled(False) #enable select region mode button again 
            self.add3DRegionBtn.setDisabled(True)
            self.fin3DRegionBtn.setDisabled(True)
        self.fin3DRegionBtn = QPushButton('(3) Exit mode')
        self.fin3DRegionBtn.clicked.connect(fin3DRegionBtnFunc)
        self.fin3DRegionBtn.setVisible(False)
        self.fin3DRegionBtn.setDisabled(True)#start disabled 
            
        def resetMeshBtnFunc():
            self.regionTable.reset()
            self.mwMesh.clear()
            self.project.mesh = None
            self.project.geom_input = {}
        self.resetMeshBtn = QPushButton('Reset Mesh')
        self.resetMeshBtn.setStyleSheet("color: red")
        self.resetMeshBtn.setFixedWidth(150)
        self.resetMeshBtn.clicked.connect(resetMeshBtnFunc)
        

        class RegionTable(QTableWidget):
            def __init__(self, nrow=1, ncol=4):
                super(RegionTable, self).__init__(nrow, ncol)
                self.nrow = nrow
                self.ncol = ncol
                self.setColumnCount(self.ncol)
                self.setRowCount(self.nrow)
                self.headers = ['|Z| [Ohm.m]', 'φ [mrad]', 'Zones', 'Fixed']
                self.setHorizontalHeaderLabels(self.headers)
                self.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
                self.setItem(0,0,QTableWidgetItem('100.0')) # resistivity [Ohm.m]
                self.setItem(0,1,QTableWidgetItem('0')) # phase [mrad]
                self.setItem(0,2,QTableWidgetItem('1')) # zone
                self.checkBoxWidget = QWidget()
                checkBoxLayout = QHBoxLayout()
                checkBoxLayout.setContentsMargins(5,5,5,5)
                checkBoxLayout.setAlignment(Qt.AlignCenter)
                checkBoxLayout.addWidget(QCheckBox())
                self.checkBoxWidget.setLayout(checkBoxLayout)
                self.setCellWidget(0,3, self.checkBoxWidget)
                

            def addRow(self):
                self.nrow = self.nrow + 1
                self.setRowCount(self.nrow)
                self.setItem(self.nrow-1, 0, QTableWidgetItem('100.0'))
                self.setItem(self.nrow-1, 1, QTableWidgetItem('0'))
                self.setItem(self.nrow-1, 2, QTableWidgetItem('1'))
                self.checkBoxWidget = QWidget()
                checkBoxLayout = QHBoxLayout()
                checkBoxLayout.setContentsMargins(5,5,5,5)
                checkBoxLayout.setAlignment(Qt.AlignCenter)
                checkBoxLayout.addWidget(QCheckBox())
                self.checkBoxWidget.setLayout(checkBoxLayout)
                self.setCellWidget(self.nrow-1, 3, self.checkBoxWidget)

            def getTable(self):
                res0 = np.zeros(self.nrow)
                phase0 = np.zeros(self.nrow)
                zones = np.zeros(self.nrow, dtype=int)
                fixed = np.zeros(self.nrow, dtype=bool)
                for j in range(self.nrow):
                    res0[j] = float(self.item(j,0).text())
                    phase0[j] = float(self.item(j,1).text())
                    zones[j] = int(self.item(j,2).text())
                    fixed[j] = self.cellWidget(j,3).findChildren(QCheckBox)[0].isChecked()
                return res0, phase0, zones, fixed

            def reset(self):
                self.nrow = 1
                self.setRowCount(1)


        self.instructionLabel = QLabel('Click on the buttons below to define a region.'
            ' If polyline is selected, you can close the region by using the right'
            ' click of the mouse. Note that it is possible to define a zone (so'
            ' no smoothing at the boundaries) and fix the region value for '
            'triangular mesh only.')
        self.instructionLabel.setWordWrap(True)
        
        self.instructionLabel3D = QLabel('Click on "Region Mode" to interactively'
            ' define a region in 3D using the box. When done, hit the "Add Region"'
            ' button to confirm. You can then specify different values in the '
            ' table on the side. Recreate a mesh to erase all regions.')
        self.instructionLabel3D.setWordWrap(True)
        self.instructionLabel3D.setVisible(False)

        self.mwMesh = MatplotlibWidget(navi=True, itight=False)
        self.mwMesh3D = MatplotlibWidget(threed=True, navi=True)
        if pvfound:
            self.meshFrame = QFrame()
            vlayout = QVBoxLayout()
            self.mesh3Dplotter = QtInteractor(self.meshFrame)
            vlayout.addWidget(self.mesh3Dplotter.interactor)
            self.meshFrame.setLayout(vlayout)
            self.mwMesh3D.setVisible(False)

        def meshLogTextFunc(text):
            cursor = self.meshLogText.textCursor()
            cursor.movePosition(cursor.End)
            cursor.insertText(text + '\n')
            self.meshLogText.ensureCursorVisible()
            QApplication.processEvents()
        self.meshLogText = QTextEdit()
        self.meshLogText.setReadOnly(True)
        
        self.regionTable = RegionTable()
        self.regionTable.setColumnHidden(1, True)


        # layout
        self.meshLayout = QVBoxLayout()

        self.meshOptionQuadLayout = QHBoxLayout()
        self.meshOptionQuadLayout.setAlignment(Qt.AlignVCenter)
        self.meshOptionQuadLayout.addWidget(self.nnodesLabel)
        self.meshOptionQuadLayout.addLayout(self.nnodesGrid, 1)

        self.meshOptionTrianLayout = QHBoxLayout()
        self.meshOptionTrianLayout.addWidget(self.clLabel)
        self.meshOptionTrianLayout.addLayout(self.clGrid, 1)
        self.meshOptionTrianLayout.addWidget(self.clFactorLabel, 1)
        self.meshOptionTrianLayout.addWidget(self.clFactorSld, 1)
        
        self.meshButtonTrianLayout = QHBoxLayout()
        self.meshButtonTrianLayout.addWidget(self.refineTrianCheck)
        self.meshButtonTrianLayout.addWidget(self.designModelBtn)
        self.meshButtonTrianLayout.addWidget(self.meshTrianBtn)
        
        self.importCustomMesh2dLayout = QVBoxLayout()
        self.importCustomMesh2dLayout.addWidget(self.importCustomMeshLabel2)
        self.importCustomMesh2dLayout.addWidget(self.importCustomMeshBtn2)

        self.meshOptionTetraLayout = QHBoxLayout()
        self.meshOptionTetraLayout.addWidget(self.cl3Label)
        self.meshOptionTetraLayout.addWidget(self.cl3Edit)
        self.meshOptionTetraLayout.addWidget(self.cl3FactorLabel)
        self.meshOptionTetraLayout.addWidget(self.cl3FactorEdit)
        self.meshOptionTetraLayout.addWidget(self.clnFactorLabel)
        self.meshOptionTetraLayout.addWidget(self.clnFactorEdit)
        self.meshOptionTetraLayout.addWidget(self.refineTetraCheck)

        self.meshOptionTankLayout = QHBoxLayout()
        self.meshOptionTankLayout.addWidget(self.clTankLabel)
        self.meshOptionTankLayout.addWidget(self.clTankEdit)
        self.meshOptionTankLayout.addWidget(self.originTankLabel)
        self.meshOptionTankLayout.addWidget(self.originTankEdit)
        self.meshOptionTankLayout.addWidget(self.dimensionTankLabel)
        self.meshOptionTankLayout.addWidget(self.dimensionTankEdit)
        # self.meshOptionTankLayout.addWidget(self.refineTankLabel)
        self.meshOptionTankLayout.addWidget(self.refineTankCheck)
        
        self.meshOptionCylinderLayout = QHBoxLayout()
        self.meshOptionCylinderLayout.addWidget(self.clCylinderLabel)
        self.meshOptionCylinderLayout.addWidget(self.clCylinderEdit)
        # self.meshOptionCylinderLayout.addWidget(self.originCylinderLabel)
        # self.meshOptionCylinderLayout.addWidget(self.originCylinderEdit)
        # self.meshOptionCylinderLayout.addWidget(self.dimensionCylinderLabel)
        # self.meshOptionCylinderLayout.addWidget(self.dimensionCylinderEdit)
        self.meshOptionCylinderLayout.addWidget(self.zlimCylinderLabel)
        self.meshOptionCylinderLayout.addWidget(self.zlimCylinderEdit)
        # self.meshOptionCylinderLayout.addWidget(self.refineCylinderLabel)
        self.meshOptionCylinderLayout.addWidget(self.refineCylinderCheck)

        self.importCustomMesh3dLayout = QVBoxLayout()
        self.importCustomMesh3dLayout.addWidget(self.importCustomMeshLabel)
        # self.importCustomMesh3dLayout.addWidget(self.importCustomMeshBtn)        
        
        self.meshChoiceLayout = QHBoxLayout()
        
        self.fmdLayout = QVBoxLayout()
        self.fmdLayout.setAlignment(Qt.AlignBottom | Qt.AlignCenter)
        
        self.meshQuadLayout = QVBoxLayout()
        self.meshQuadLayout.setAlignment(Qt.AlignBottom | Qt.AlignHCenter)
        self.meshTrianLayout = QVBoxLayout()
        self.meshTetraLayout = QVBoxLayout()
        self.meshTetraLayout.setAlignment(Qt.AlignHCenter)
        self.meshTankLayout = QVBoxLayout()
        self.meshTankLayout.setAlignment(Qt.AlignHCenter)
        self.meshCylinderLayout = QVBoxLayout()
        self.meshCylinderLayout.setAlignment(Qt.AlignHCenter)
        
        self.fmdGroup = QGroupBox()
        self.fmdGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        self.meshQuadGroup = QGroupBox()
        self.meshQuadGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        self.meshTrianGroup = QGroupBox()
        self.meshTrianGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        self.meshCustom2dGroup = QGroupBox()
        self.meshCustom2dGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        self.mesh3DGroup = QGroupBox()
        self.mesh3DGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        self.mesh3DGroup.setHidden(True)
        self.meshTetraGroup = QGroupBox()
        self.meshTetraGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        self.meshTetraGroup.setHidden(True)
        self.meshTankGroup = QGroupBox()
        self.meshTankGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        self.meshTankGroup.setHidden(True)
        self.meshCylinderGroup = QGroupBox()
        self.meshCylinderGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        self.meshCylinderGroup.setHidden(True)
        self.meshCustom3dGroup = QGroupBox()
        self.meshCustom3dGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        self.meshCustom3dGroup.setHidden(True)
        
        self.mesh3DComboLayout = QVBoxLayout()
        self.mesh3DComboLayout.addWidget(self.mesh3DCombo)
        self.mesh3DComboLayout.addWidget(self.mesh3DBtn)
        self.mesh3DGroup.setLayout(self.mesh3DComboLayout)
        self.meshChoiceLayout.addWidget(self.mesh3DGroup)
        
        self.fmdLayout.addWidget(self.fmdLabel)
        self.fmdLayout.addWidget(self.fmdBox)
        self.fmdGroup.setLayout(self.fmdLayout)
        self.meshChoiceLayout.addWidget(self.fmdGroup,0)
        
        self.meshQuadLayout.addLayout(self.meshOptionQuadLayout, 1)
        self.meshQuadLayout.addWidget(self.meshQuad, 0)
        self.meshQuadGroup.setLayout(self.meshQuadLayout)
        self.meshChoiceLayout.addWidget(self.meshQuadGroup,35)

        self.meshTrianLayout.addLayout(self.meshOptionTrianLayout)
        self.meshTrianLayout.addLayout(self.meshButtonTrianLayout)
        self.meshTrianGroup.setLayout(self.meshTrianLayout)
        self.meshChoiceLayout.addWidget(self.meshTrianGroup,65)

        self.meshCustom2dGroup.setLayout(self.importCustomMesh2dLayout)
        self.meshChoiceLayout.addWidget(self.meshCustom2dGroup,0)
        
        self.meshTetraLayout.addLayout(self.meshOptionTetraLayout)
        # self.meshTetraLayout.addWidget(self.meshTetraBtn)
        self.meshTetraGroup.setLayout(self.meshTetraLayout)
        self.meshChoiceLayout.addWidget(self.meshTetraGroup, 1)
        
        self.meshTankLayout.addLayout(self.meshOptionTankLayout)
        # self.meshTankLayout.addWidget(self.meshTankBtn)
        self.meshTankGroup.setLayout(self.meshTankLayout)
        self.meshChoiceLayout.addWidget(self.meshTankGroup, 1)
        
        self.meshCylinderLayout.addLayout(self.meshOptionCylinderLayout)
        # self.meshCylinderLayout.addWidget(self.meshCylinderBtn)
        self.meshCylinderGroup.setLayout(self.meshCylinderLayout)
        self.meshChoiceLayout.addWidget(self.meshCylinderGroup, 1)
        
        self.meshCustom3dGroup.setLayout(self.importCustomMesh3dLayout)
        self.meshChoiceLayout.addWidget(self.meshCustom3dGroup, 1)
        
        self.meshLayout.addLayout(self.meshChoiceLayout, 0)

        self.instructionLayout = QHBoxLayout()
        self.instructionLayout.addWidget(self.instructionLabel, 86)
        self.instructionLayout.addWidget(self.instructionLabel3D, 86)
        self.instructionLayout.addWidget(self.meshAspectBtn, 7)
        self.instructionLayout.addWidget(self.resetMeshBtn, 7)
        ## for 3D forward modelling / region selection 
        self.instructionLayout.addWidget(self.select3DRegionBtn, 7)
        self.instructionLayout.addWidget(self.add3DRegionBtn, 7)
        self.instructionLayout.addWidget(self.fin3DRegionBtn, 7)
        self.instructionLayout.addWidget(self.saveMeshBtn)
        self.meshLayout.addLayout(self.instructionLayout)
        
        # for RAM issue
        self.ramRequiredLabel = QLabel('<font color="red">Make a coarser mesh!! It is likely that <b><u>more RAM is required</u></b> for inversion!</font>')
        self.ramRequiredLabel.setAlignment(Qt.AlignCenter)
        self.ramRequiredLabel.hide()
        self.meshLayout.addWidget(self.ramRequiredLabel)
        
        self.regionLayout = QVBoxLayout()
        self.regionLayout.addWidget(self.regionTable)

        self.meshPlot = QWidget()
        self.meshPlotLayout = QHBoxLayout()
        self.meshPlotLayout.addWidget(self.mwMesh)
        self.meshPlot.setLayout(self.meshPlotLayout)

        self.meshPlot3D = QWidget()
        self.meshPlot3DLayout = QHBoxLayout()
        self.meshPlot3DLayout.addWidget(self.mwMesh3D)
        if pvfound:
            self.meshPlot3DLayout.addWidget(self.meshFrame)
        self.meshPlot3D.setLayout(self.meshPlot3DLayout)

        self.meshOutputStack = QStackedLayout()
        self.meshOutputStack.addWidget(self.meshLogText)
        self.meshOutputStack.addWidget(self.meshPlot)
        self.meshOutputStack.addWidget(self.meshPlot3D)
        self.meshOutputStack.setCurrentIndex(0)

        self.meshSubLayout = QHBoxLayout()
        self.meshSubLayout.addLayout(self.meshOutputStack, 70)
        self.meshSubLayout.addLayout(self.regionLayout, 30)
        
        self.meshLayout.addLayout(self.meshSubLayout, 1)

        self.tabMesh.setLayout(self.meshLayout)


#%% =================== Tab for forward model ===================
        self.tabForward = QWidget()
        self.tabs.addTab(self.tabForward, 'Forward model')
        self.tabs.setTabEnabled(3, False)
        
        self.fwdSlpitterLayout = QHBoxLayout() # a splitter so the graphs are easier to visualize
        self.fwdSplitter = QSplitter(Qt.Vertical)

        self.seqLabel = QLabel('Design a sequence for the forward modelling. A \
combination of multiple sequence is accepted as well as importing a custom sequence')
        self.seqLabel.setWordWrap(True)
        self.seqLabel.setAlignment(Qt.AlignTop)
        
        # alternative design
        seqData = [('dpdp1', 'Dipole-Dipole', ['a','n']),
                   ('wenner', 'Wenner', ['a']),
                   ('schlum1', 'Schlumberger', ['a','m']),
                   ('multigrad', 'Multi-Gradient', ['a','m','n']),
                   ('custSeq', 'Custom Sequence', [])]
        
        
        class RowOpt(QHBoxLayout):
            def __init__(self, parent=None):
                super(QHBoxLayout, self).__init__()
                self.combo = None
                self.rmBtn = None
                self.labels = []
                self.fields = []
                self.seq = 'dpdp1'
                self.iseq = 0
                self.importBtn = None
                self.fname = ''
                self.parent = parent
#                self.createRow() # this create floating windows
#                self.showArg()
                
            def createRow(self):
                self.combo = QComboBox()
                for row in seqData:
                    self.combo.addItem(row[1])
                self.combo.activated.connect(self.comboFunc)
                self.addWidget(self.combo)
                for a, b in zip(['a','n','m'], ['1','8','']):
                    lab = QLabel(a + '=')
                    field = QLineEdit(b)
                    field.setValidator(QIntValidator())
                    self.labels.append(lab)
                    self.fields.append(field)
                    self.addWidget(lab)
                    self.addWidget(field)
                self.importBtn = QPushButton('Import Custom Sequence')
                self.importBtn.clicked.connect(self.importFile)
                self.importBtn.setVisible(False)
                self.addWidget(self.importBtn)
                self.rmBtn = QPushButton('Remove')
                self.rmBtn.clicked.connect(self.remove)
                self.addWidget(self.rmBtn)
            
            def importFile(self):
                fname, _ = QFileDialog.getOpenFileName(self.parent.tabImportingData,'Open File', self.parent.datadir)
                if fname != '':
                    self.fname = fname
                    self.importBtn.setText(os.path.basename(fname))
                                 
            def comboFunc(self, i):
                self.seq = seqData[i][0]
                showArray(self.seq) # display help aside
                self.iseq = i
                self.showArg()
                
            def showArg(self):
                n = len(seqData[self.iseq][2])
                if self.iseq == 4: # custom sequence
                    self.importBtn.setVisible(True)
                else:
                    self.importBtn.setVisible(False)
                for i in range(3):
                    val = False
                    if i < n:
                        val = True
                    self.labels[i].setVisible(val)
                    self.fields[i].setVisible(val)
                if self.seq == 'multigrad':
                    self.labels[-1].setText('s=')
                    
            def remove(self):
                if self.seq == 'custSeq':
                    self.fname = ''
                self.combo.deleteLater()
                for w in self.fields:
                    w.deleteLater()
                for w in self.labels:
                    w.deleteLater()
                self.rmBtn.deleteLater()
                self.importBtn.deleteLater()
                self.deleteLater()
                
            
            def getData(self):
                if self.seq == 'custSeq':
                    return (self.seq, self.fname)
                else:
                    n = len(seqData[self.iseq][2])
                    vals = []
                    for i, field in enumerate(self.fields):
                        if i < n:
                            val = int(field.text()) if field.text() != '' else -9999
                            vals.append(val)
                    return (self.seq, *vals)


        self.seqRowLayout = QVBoxLayout()
        self.seqRowLayout.setAlignment(Qt.AlignTop)
        seqRows = []
        seqRow = RowOpt(parent=self)
        self.seqRowLayout.addLayout(seqRow)
        seqRows.append(seqRow)
        
        DpDp = resource_path('image/dipdip.png')
        Wenner = resource_path('image/wenner.png')
        Schlum = resource_path('image/schlum.png')
        Gradient = resource_path('image/gradient.png')
        
        seqHelp = {'dpdp1' : '<img height=140 src="%s">' % DpDp,
           'wenner': '<img height=140 src="%s">' % Wenner,
           'schlum1': '<img height=140 src="%s">' % Schlum,
           'multigrad': '<img height=140 src="%s">' % Gradient,
           'custSeq': 'Use the button to import a custom CSV file (with headers)\ncolumn1: C+, column2: C-, column3: P+, column4: P-\n' \
               'It is recommended to use a custom sequence in case of "unconventional surveys"'
            }
        
        self.arrayLabel = QLabel('Sequence help will be display here.')
        self.arrayLabel.setAlignment(Qt.AlignCenter)
        def showArray(arg):
            if arg not in seqHelp.keys():
                self.arrayLabel.setText('Sequence help not found.')
            else:
                self.arrayLabel.setText(seqHelp[arg])
        showArray('dpdp1') # default
        
        def addRowBtnFunc():
            a = RowOpt(parent=self)
            seqRows.append(a)
            self.seqRowLayout.addLayout(a)
            a.createRow()
            a.showArg()
        self.addRowBtn = QPushButton('Add sequence')
        self.addRowBtn.adjustSize()
        self.addRowBtn.clicked.connect(addRowBtnFunc)
        
        def getDataBtnFunc():
            vals = []
            for seqRow in seqRows:
                try: # because we 'deleteLater() the object it causes error
                    # would be better to check if the object exists
                    val = seqRow.getData()
                    ok = True
                    for j, a in enumerate(val):
                        if (j != 0) & (a == -9999):
                            ok = False
                    if ok is True:
                        vals.append(val)
                except:
                    pass
            return vals

        
        def seqCreateFunc():
            if self.project.elec is None:
                self.errorDump('Input electrode positions in the "Electrodes (XYZ/Topo)" tab first.')
                return
            params = getDataBtnFunc()
            if len(params) == 0:
                raise ValueError('You must specify at least one sequence.')
                return
            self.seqIdx = self.project.createSequence(params=params)
            self.writeLog('k.createSequence(params={:s})'.format(str(params)))
            self.seqOutputLabel.setText('{:d} quadrupoles in total'.format(len(self.project.sequence)))

        self.seqOutputLabel = QLabel('')
        
        # add noise possibility
        self.noiseLabel = QLabel('Resistivity noise [%]:')
        self.noiseEdit = QLineEdit('2')
        self.noiseEdit.setValidator(QDoubleValidator())

        # add IP noise
        self.noiseLabelIP = QLabel('Phase noise [mrad]:')
        self.noiseEditIP = QLineEdit('2')
        self.noiseLabelIP.hide()
        self.noiseEditIP.hide()
        self.noiseEditIP.setValidator(QDoubleValidator())

        # save sequence button
        def saveSeqBtnFunc():
            fname, _ = QFileDialog.getSaveFileName(self.tabImportingData,'Save File', self.datadir, 'Comma Separated Values (*.csv)')
            if fname != '':
                self.project.saveSequence(fname)
                self.writeLog('k.saveSequence("{:s}")'.format(fname))
        self.saveSeqBtn = QPushButton('Save Sequence')
        self.saveSeqBtn.setToolTip('This will save the sequence of the fwd modeling. Output data is already saved in <i>fwd</i> folder in the <i>working directory</i>.')
        self.saveSeqBtn.clicked.connect(saveSeqBtnFunc)

        # add a forward button
        def forwardBtnFunc():
            if self.project.mesh is None: # we need to create mesh to assign starting resistivity
                self.errorDump('Please specify a mesh and an initial model first.')
                return
            # try:
            seqCreateFunc()
            # except:
                # self.errorDump('Error in sequence generation! Use a custom sequence instead.')
                # return
            if len(self.project.sequence) == 0:
                self.errorDump('Sequence is empty, can not run forward model.')
                return
            forwardOutputStack.setCurrentIndex(0)
            self.forwardLogText.clear()
            QApplication.processEvents()
            # apply region for initial model

            x, phase0, zones, fixed = self.regionTable.getTable()
            regid = np.arange(len(x)) + 1 # region 0 doesn't exist
            pdebug('forwardBtnFunc(): with {:d} regions'.format(len(x)))
            self.project.setStartingRes(dict(zip(regid, x)),
                                   dict(zip(regid, zones)),
                                   dict(zip(regid, fixed)),
                                   dict(zip(regid, phase0)))
            self.writeLog('k.setStartingRes({:s}, {:s}, {:s}, {:s})'.format(
                                   str(dict(zip(regid, x))),
                                   str(dict(zip(regid, zones))),
                                   str(dict(zip(regid, fixed))),
                                   str(dict(zip(regid, phase0)))))
            noise = float(self.noiseEdit.text())
            noiseIP = float(self.noiseEditIP.text())
            self.project.forward(noise=noise, noiseIP=noiseIP, iplot=False, dump=forwardLogTextFunc)
            self.writeLog('k.forward(noise={:.2f}, noiseIP={:.2f})'.format(noise, noiseIP))
            if self.project.typ[-1] == '2':
                self.calcAspectRatio() # doesn't work for 3D?
            
            ### pseudo section plotting! ###
            if pvfound and self.project.typ[-1] == 't': # 3D pseudo-sections?
                self.mwFwdPseudo.hide()
                self.pseudo3Dplotter.clear() # clear all actors

                self.project.surveys[0].showPseudo(ax=self.pseudo3Dplotter, threed=True, 
                                              strIdx=self.seqIdx, darkMode=eval(resipySettings.param['dark']))
                # self.fwdContour.setDisabled(True)
                self.pseudoFramefwd.setVisible(True)
            else:
                self.mwFwdPseudo.plot(self.project.surveys[0].showPseudo, aspect='auto')
                
            self.fwdContour.setVisible(True)
            self.tabs.setTabEnabled(4, True)
            self.tabs.setTabEnabled(5, True)
            # self.tabs.setTabEnabled(6, True)
            if self.project.typ[0] == 'c':
                if pvfound and self.project.typ[-1]=='t':
                    self.mwFwdPseudoIP.hide()
                    self.pseudo3DplotterIP.clear() # clear all actors 
                    self.project.surveys[0].showPseudoIP(ax=self.pseudo3DplotterIP, threed=True, 
                                                         strIdx=self.seqIdx, darkMode=eval(resipySettings.param['dark']))
                    self.pseudoFrameIPfwd.setVisible(True)
                else:
                    self.mwFwdPseudoIP.plot(self.project.surveys[0].showPseudoIP, aspect='auto')
                
        self.forwardBtn = QPushButton('Forward Modelling')
        self.forwardBtn.setAutoDefault(True)
        self.forwardBtn.clicked.connect(forwardBtnFunc)
        self.forwardBtn.setStyleSheet('background-color: green; color:black')

        self.mwFwdPseudo = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        self.mwFwdPseudoIP = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        self.mwFwdPseudoIP.setVisible(False)
        if pvfound:
            self.pseudoFramefwd = QFrame()
            vlayout = QVBoxLayout()
            self.pseudo3Dplotter = QtInteractor(self.pseudoFramefwd)
            vlayout.addWidget(self.pseudo3Dplotter.interactor)
            self.pseudoFramefwd.setLayout(vlayout)
            self.pseudoFramefwd.setVisible(False)
            
            #IP 3D fwd pseudo section
            self.pseudoFrameIPfwd = QFrame()
            vlayout = QVBoxLayout()
            self.pseudo3DplotterIP = QtInteractor(self.pseudoFrameIPfwd)
            vlayout.addWidget(self.pseudo3DplotterIP.interactor)
            self.pseudoFrameIPfwd.setLayout(vlayout)
            self.pseudoFrameIPfwd.setVisible(False)

        self.forwardLogText = QTextEdit()
        self.forwardLogText.setReadOnly(True)
        def forwardLogTextFunc(text, end=''):
            cursor = self.forwardLogText.textCursor()
            cursor.movePosition(cursor.End)
            cursor.insertText(text + end)
            self.forwardLogText.ensureCursorVisible()
            QApplication.processEvents()
            if text == 'Forward modelling done.':
                forwardOutputStack.setCurrentIndex(1) # switch to graph
                if pvfound and self.project.typ[-1]=='t':
                    forwardOutputStack.setCurrentIndex(2)
                
        def fwdContourFunc(state):
            if state == Qt.Checked:
                contour = True
            else:
                contour = False

            if pvfound and self.project.typ[-1] == 't': 
                self.pseudo3Dplotter.clear() # clear all actors
                self.project.surveys[0].showPseudo(ax=self.pseudo3Dplotter, threed=True, contour=contour,
                                                   strIdx=self.seqIdx, darkMode=eval(resipySettings.param['dark']))
                if self.project.typ[0] == 'c':
                    self.pseudo3DplotterIP.clear() # clear all actors 
                    self.project.surveys[0].showPseudoIP(ax=self.pseudo3DplotterIP, threed=True, contour=contour,
                                                         strIdx=self.seqIdx, darkMode=eval(resipySettings.param['dark']))
            else:    
                self.mwFwdPseudo.setCallback(self.project.surveys[0].showPseudo)
                self.mwFwdPseudo.replot(aspect='auto', contour=contour)
                self.writeLog('k.showPseudo(contour={:s})'.format(str(contour)))
                if self.project.typ[0] == 'c':
                    self.mwFwdPseudoIP.setCallback(self.project.surveys[0].showPseudoIP)
                    self.mwFwdPseudoIP.replot(aspect='auto', contour=contour)
                    self.writeLog('k.showPseudoIP(contour={:s})'.format(str(contour)))
        
        self.fwdContour = QCheckBox('Contour')
        self.fwdContour.setToolTip('Check/uncheck to contour pseudo section plots (in 3D mode "Delaunay 3D" is used).')
        self.fwdContour.stateChanged.connect(fwdContourFunc)
        self.fwdContour.setVisible(False)


        # layout
        forwardLayout = QVBoxLayout()
        forwardLayout.setAlignment(Qt.AlignTop)
        seqLayout = QHBoxLayout()
        noiseLayout = QHBoxLayout()

        # top part
        seqOptionLayout = QVBoxLayout()
        seqOptionLayout.setAlignment(Qt.AlignTop)
        seqOptionLayout.addLayout(self.seqRowLayout)
        seqOptionLayout.addWidget(self.addRowBtn)
        seqLayout.addLayout(seqOptionLayout, 50)
        seqLayout.addWidget(self.arrayLabel, 50)
        
        # noise Layout
        noiseLayout.addWidget(self.noiseLabel)
        noiseLayout.addWidget(self.noiseEdit)
        noiseLayout.addWidget(self.noiseLabelIP)
        noiseLayout.addWidget(self.noiseEditIP)
        noiseLayout.addWidget(self.forwardBtn)
        noiseLayout.addWidget(self.seqOutputLabel)
        noiseLayout.addWidget(self.saveSeqBtn)
        noiseLayout.addWidget(self.fwdContour)

        # pseudo dynamic layout
        forwardPseudoLayout = QVBoxLayout()
        
        forwardPseudoLayoutBottom = QHBoxLayout()
        forwardPseudoLayoutBottom.addWidget(self.mwFwdPseudo, 50)
        forwardPseudoLayoutBottom.addWidget(self.mwFwdPseudoIP, 50)
        self.mwFwdPseudoIP.hide()
        
        if pvfound:
            self.pvFwdBottomWidget = QWidget()
            pvFwdBottomLayout = QHBoxLayout()
            pvFwdBottomLayout.addWidget(self.pseudoFramefwd, 50)
            pvFwdBottomLayout.addWidget(self.pseudoFrameIPfwd, 50)
            self.pvFwdBottomWidget.setLayout(pvFwdBottomLayout)
        
        forwardPseudoLayout.addLayout(forwardPseudoLayoutBottom)

        self.forwardPseudos = QWidget()
        self.forwardPseudos.setLayout(forwardPseudoLayout)

        forwardOutputStack = QStackedLayout()
        forwardOutputStack.addWidget(self.forwardLogText)
        forwardOutputStack.addWidget(self.forwardPseudos)
        if pvfound:
            forwardOutputStack.addWidget(self.pvFwdBottomWidget)
        forwardOutputStack.setCurrentIndex(0)
                
        # general forward layout
        forwardLayout.addWidget(self.seqLabel)
        forwardLayout.addLayout(seqLayout)
        forwardLayout.addLayout(noiseLayout)

        self.fwdTopWidget = QWidget()
        self.fwdTopWidget.setObjectName('self.fwdTopWidget')
        self.fwdTopWidget.setStyleSheet("QWidget#self.fwdTopWidget {border:1px solid rgb(185,185,185)}")
        self.fwdTopWidget.setLayout(forwardLayout)
        
        #bottom part 
        self.fwdBottomWidget = QWidget()
        self.fwdBottomWidget.setObjectName('self.fwdBottomWidget')
        self.fwdBottomWidget.setStyleSheet("QWidget#fwdBottomWidget {border:1px solid rgb(185,185,185)}")
        self.fwdBottomWidget.setLayout(forwardOutputStack)
        
        self.fwdSplitter.addWidget(self.fwdTopWidget)
        self.fwdSplitter.addWidget(self.fwdBottomWidget)
        self.fwdSplitter.setCollapsible(1, False)
        
        self.fwdSlpitterLayout.addWidget(self.fwdSplitter)
        
        # instantiate the first row of the table
        # need to do this here as the layout needs to be integrated otherwise
        # it creates floating windows at the start
        seqRow.createRow()
        seqRow.showArg()
        
        self.tabForward.setLayout(self.fwdSlpitterLayout)


        #%% tab INVERSION SETTINGS
        self.tabInversionSettings = QTabWidget()
        self.tabs.addTab(self.tabInversionSettings , 'Inversion settings')
        self.tabs.setTabEnabled(4, False)

        # general tab
        self.generalSettings = QWidget()
        generalLayout = QHBoxLayout()
        invForm = QFormLayout()

        # advanced tab
        self.advancedSettings = QWidget()
        advancedLayout = QHBoxLayout()
        advForm = QFormLayout()

        def showIpOptions(arg):
            settingsIP = [self.min_error, self.min_errorLabel,
                          ]
            settingsDC = [self.data_type, self.data_typeLabel,
                          self.reg_mode, self.reg_modeLabel,
                          ]
            if arg == True:
                self.a_wgt.setText('0.02')
                self.b_wgt.setText('2')
                self.inv_type.clear()
                self.inv_type.addItem('Pseudo Marquardt [0]')
                self.inv_type.addItem('Regularized Inversion with Linear Filtering [1]')
                self.inv_type.addItem('Blocked Linear Regularized Inversion [4]')
                self.inv_type.setCurrentIndex(1)
                [o.setVisible(True) for o in settingsIP]
                [o.setVisible(False) for o in settingsDC]
            else:
                self.a_wgt.setText('0.01')
                self.b_wgt.setText('0.02')
                self.inv_type.clear()
                self.inv_type.addItem('Pseudo Marquardt [0]')
                self.inv_type.addItem('Regularized Inversion with Linear Filtering [1]')
                self.inv_type.addItem('Regularized Inversion with Quadratic Filtering [2]')
                self.inv_type.addItem('Qualitative Solution [3]')
                self.inv_type.addItem('Blocked Linear Regularized Inversion [4]')
                self.inv_type.setCurrentIndex(1)
                [o.setVisible(False) for o in settingsIP]
                [o.setVisible(True) for o in settingsDC]
            if self.project.iForward is False and len(self.project.surveys) >= 1:
                if 'magErr' in self.project.surveys[0].df.columns:
                    self.a_wgt.setText('0.0')
                    self.b_wgt.setText('0.0')

        def show3DOptions(arg):
            settings3D = []
            settings2D = [self.flux_type, self.flux_typeLabel,
                          self.scale, self.scaleLabel,
                          self.inv_type, self.inv_typeLabel]
            if arg is True:
                [s.setVisible(True) for s in settings3D]
                [s.setVisible(False) for s in settings2D]
            else:
                [s.setVisible(False) for s in settings3D]
                [s.setVisible(True) for s in settings2D]

        # help sections
        def showHelp(arg): # for general tab
            if arg not in r2help:
                self.helpSection.setText('SORRY NOT IN HELP')
            else:
                self.helpSection.setHtml(r2help[arg])

        def showHelpAdv(arg): # for advanced tab
            if arg not in r2help:
                self.helpSection2.setText('SORRY NOT IN HELP')
            else:
                self.helpSection2.setHtml(r2help[arg])


        self.parallelLabel = QLabel('<a href="parallel">Parallel inversion</a>')
        self.parallelLabel.linkActivated.connect(showHelpAdv)
        self.parallelCheck = QCheckBox()
        advForm.addRow(self.parallelLabel, self.parallelCheck)
        
        # decide number of cpu cores to use 
        self.ncoresLabel = QLabel('<a href="ncores">Number of parallel threads</a>')
        self.ncoresLabel.linkActivated.connect(showHelpAdv)
        self.ncoresText = QLineEdit()
        self.ncoresText.setText('%i'%sysinfo['core_count'])
        advForm.addRow(self.ncoresLabel, self.ncoresText)

        self.modErrLabel = QLabel('<a href="modErr">Compute Modelling Error</a>')
        self.modErrLabel.linkActivated.connect(showHelpAdv)
        self.modErrCheck = QCheckBox()
        advForm.addRow(self.modErrLabel, self.modErrCheck)

        def notCroppingFunc(state):
            if state == Qt.Checked:
                self.iCropping = False
                if 'num_xz_poly' in self.project.param:
                    self.num_xz_poly = self.project.param['num_xz_poly'] # store value
                    self.clipCornersCheck.setChecked(False)
                    self.clipCornersCheck.setEnabled(False)
                elif 'num_xy_poly' in self.project.param:
                    self.num_xz_poly = self.project.param['num_xy_poly'] # store value
            else:
                self.iCropping = True # default
                if ('num_xz_poly' in self.project.param) and (self.num_xz_poly is not None):
                    self.project.param['num_xz_poly'] = self.num_xz_poly # restore value
                    self.writeLog('k.param["num_xz_poly"] = {:s}'.format(str(self.num_xz_poly)))
                    self.clipCornersCheck.setEnabled(True)
                elif ('num_xy_poly' in self.project.param) and (self.num_xz_poly is not None):
                    self.project.param['num_xy_poly'] = self.num_xz_poly # restore value
                    self.writeLog('k.param["num_xy_poly"] = {:s}'.format(str(self.num_xy_poly)))
        self.notCroppingLabel = QLabel('<a href="notCropping">Do not crop the output</a>')
        self.notCroppingLabel.linkActivated.connect(showHelpAdv)
        self.notCropping = QCheckBox()
        self.notCropping.stateChanged.connect(notCroppingFunc)
        advForm.addRow(self.notCroppingLabel, self.notCropping)
                    
        self.cropBelowFmdLabel = QLabel('<a href="cropBelowFmd">Crop below mesh fine region</a>')
        self.cropBelowFmdLabel.linkActivated.connect(showHelpAdv)
        self.cropBelowFmd = QCheckBox()
        self.cropBelowFmd.setChecked(True)
        # not connected but read when calling displayInvertedResults()
        advForm.addRow(self.cropBelowFmdLabel, self.cropBelowFmd)

        def modelDOIFunc(status):
            if status == Qt.Checked:
                self.doiCheck.setVisible(True)
                self.doiSensCheck.setVisible(False)
                self.displayParams['doiSens'] = False
            else:
                self.doiCheck.setVisible(False)
                self.doiSensCheck.setVisible(True)
        self.modelDOILabel = QLabel('<a href="modelDOI">Model DOI</a>')
        self.modelDOILabel.linkActivated.connect(showHelpAdv)
        self.modelDOICheck = QCheckBox()
        self.modelDOICheck.stateChanged.connect(modelDOIFunc)
        advForm.addRow(self.modelDOILabel, self.modelDOICheck)
        
        def checkTxSignFunc(state):
            if state == Qt.Checked:
                self.project.param['checkTxSign'] = True
                self.writeLog('k.param["checkTxSign"] = True')
            else:
                self.project.param['checkTxSign'] = False
                self.writeLog('k.param["checkTxSign"] = False')
        self.checkTxSignLabel = QLabel('<a href="txSign">Resistance polarity check</a>')
        self.checkTxSignLabel.linkActivated.connect(showHelpAdv)
        self.checkTxSign = QCheckBox()
        self.checkTxSign.stateChanged.connect(checkTxSignFunc)
        advForm.addRow(self.checkTxSignLabel, self.checkTxSign)

        def flux_typeFunc(index):
            if index == 0:
                self.project.param['flux_type'] = 3
                self.writeLog('k.param["flux_type"] = 3')
            else:
                self.project.param['flux_type'] = 2
                self.writeLog('k.param["flux_type"] = 2')
        self.flux_typeLabel = QLabel('<a href="flux_type">Flux Type</a>')
        self.flux_typeLabel.linkActivated.connect(showHelp)
        self.flux_type = QComboBox()
        self.flux_type.addItem('3D')
        self.flux_type.addItem('2D')
        self.flux_type.activated.connect(flux_typeFunc)
        advForm.addRow(self.flux_typeLabel, self.flux_type)

        def singular_typeFunc(state):
            if state == Qt.Checked:
                self.project.param['singular_type'] = 1
                self.writeLog('k.param["singular_type"] = 1')
            else:
                self.project.param['singular_type'] = 0
                self.writeLog('k.param["singular_type"] = 0')
        self.singular_typeLabel = QLabel('<a href="singular_type">Remove Singularity</a>')
        self.singular_typeLabel.linkActivated.connect(showHelp)
        self.singular_type = QCheckBox()
        self.singular_type.stateChanged.connect(singular_typeFunc)
        advForm.addRow(self.singular_typeLabel, self.singular_type)

        def res_matrixFunc(index):
            self.project.param['res_matrix'] = index
            self.writeLog('k.param["res_matrix"] = {:d}'.format(index))
        self.res_matrixLabel = QLabel('<a href="res_matrix">Value for <code>res_matrix</code><a/>')
        self.res_matrixLabel.linkActivated.connect(showHelpAdv)
        self.res_matrix = QComboBox()
        self.res_matrix.addItem('No sensisitivity/resolution matrix [0]')
        self.res_matrix.addItem('Sensitivity matrix [1]')
        self.res_matrix.addItem('True Resolution Matrix [2]')
        self.res_matrix.addItem('Sensitivity map [3]')
        self.res_matrix.setCurrentIndex(1)
        self.res_matrix.activated.connect(res_matrixFunc)
        advForm.addRow(self.res_matrixLabel, self.res_matrix)

        def scaleFunc():
            a = float(self.scale.text())
            self.project.param['scale'] = a
            self.writeLog('k.param["scale"] = {:.2f}'.format(a))
        self.scaleLabel = QLabel('<a href="scale"> Scale for triangular mesh</a>')
        self.scaleLabel.linkActivated.connect(showHelp)
        self.scaleLabel.setVisible(False)
        self.scale = QLineEdit()
        self.scale.setValidator(QDoubleValidator())
        self.scale.setText('1.0')
        self.scale.editingFinished.connect(scaleFunc)
        self.scale.setVisible(False)
        invForm.addRow(self.scaleLabel, self.scale)

        def patch_xFunc():
            a = int(self.patch_x.text())
            self.project.param['patch_x'] = a
            self.writeLog('k.param["patch_x"] = {:d}'.format(a))
        self.patch_xLabel = QLabel('<a href="patch">Patch size x<a/>:')
        self.patch_xLabel.linkActivated.connect(showHelpAdv)
        self.patch_x = QLineEdit()
        self.patch_x.setValidator(QIntValidator())
        self.patch_x.setText('1')
        self.patch_x.editingFinished.connect(patch_xFunc)
        advForm.addRow(self.patch_xLabel, self.patch_x)

        def patch_zFunc():
            a = int(self.patch_z.text())
            self.project.param['patch_z'] = a
            self.writeLog('k.param["patch_z"] = {:d}'.format(a))
        self.patch_zLabel = QLabel('<a href="patch">Patch size y<a/>:')
        self.patch_zLabel.linkActivated.connect(showHelpAdv)
        self.patch_z = QLineEdit()
        self.patch_z.setValidator(QIntValidator())
        self.patch_z.setText('1')
        self.patch_z.editingFinished.connect(patch_zFunc)
        advForm.addRow(self.patch_zLabel, self.patch_z)

        self.inv_typeVisible = []
        def inv_typeFunc(arg):
            index = int(self.inv_type.currentText()[-2])
            self.project.param['inverse_type'] = index
            self.writeLog('k.param["inverse_type"] = {:d}'.format(index))
            opts = [self.data_typeLabel, self.data_type,
                    self.reg_modeLabel, self.reg_mode,
                    self.toleranceLabel, self.tolerance,
                    self.max_iterationsLabel, self.max_iterations,
                    self.error_modLabel, self.error_mod,
                    self.alpha_anisoLabel, self.alpha_aniso,
                    self.a_wgtLabel, self.a_wgt,
                    self.b_wgtLabel, self.b_wgt]
            if index == 3: # qualitative solution
                self.inv_typeVisible = [o.isVisible() for o in opts]
                [o.setVisible(False) for o in opts]
            else:
                if len(self.inv_typeVisible) > 1:
                    [o.setVisible(a) for o, a in zip(opts, self.inv_typeVisible)]
                self.inv_typeVisible = []
        self.inv_typeLabel = QLabel('<a href="inverse_type">Inversion Type</a>:')
        self.inv_typeLabel.linkActivated.connect(showHelp)
        self.inv_type = QComboBox()
        self.inv_type.addItem('Pseudo Marquardt [0]')
        self.inv_type.addItem('Regularized Inversion with Linear Filtering [1]')
        self.inv_type.addItem('Regularized Inversion with Quadratic Filtering [2]')
        self.inv_type.addItem('Qualitative Solution [3]')
        self.inv_type.addItem('Blocked Linear Regularized Inversion [4]')
        self.inv_type.setCurrentIndex(1)
        self.inv_type.activated.connect(inv_typeFunc)
        invForm.addRow(self.inv_typeLabel, self.inv_type)

        def target_decreaseFunc():
            a = float(self.target_decrease.text())
            self.project.param['target_decrease'] = a
            self.writeLog('k.param["target_decrease"] = {:.2f}'.format(a))
        self.target_decreaseLabel = QLabel('<a href="target_decrease">Target decrease</a>:')
        self.target_decreaseLabel.linkActivated.connect(showHelp)
        self.target_decrease = QLineEdit()
        self.target_decrease.setValidator(QDoubleValidator())
        self.target_decrease.setText('0')
        self.target_decrease.editingFinished.connect(target_decreaseFunc)
        invForm.addRow(self.target_decreaseLabel, self.target_decrease)


        def data_typeFunc(index):
            self.project.param['data_type'] = index
            self.writeLog('k.param["data_type"] = {:d}'.format(index))
        self.data_typeLabel = QLabel('<a href="data_type">Data type</a>:')
        self.data_typeLabel.linkActivated.connect(showHelp)
        self.data_type = QComboBox()
        self.data_type.addItem('Normal [0]')
        self.data_type.addItem('Logarithmic [1]')
        self.data_type.setCurrentIndex(1)
        self.data_type.activated.connect(data_typeFunc)
        invForm.addRow(self.data_typeLabel, self.data_type)

        def reg_modeFunc(index):
            self.project.param['reg_mode'] = index
            self.writeLog('k.param["reg_mode"] = {:d}'.format(index))
            if int(index) == 1: # enable qedit for alpha_s (as only relevant to timelapse)
                self.alpha_sLabel.setVisible(True)
                self.alpha_s.setVisible(True)
            else:
                self.alpha_sLabel.setVisible(False)
                self.alpha_s.setVisible(False)
        self.reg_modeLabel = QLabel('<a href="reg_mode">Regularization mode</a>:')
        self.reg_modeLabel.linkActivated.connect(showHelp)
        self.reg_mode = QComboBox()
        self.reg_mode.addItem('Normal regularization [0]')
        self.reg_mode.addItem('Regularization from initial model [1]')
        self.reg_mode.addItem('Regularization from difference inversion [2]')
        self.reg_mode.activated.connect(reg_modeFunc)
        invForm.addRow(self.reg_modeLabel, self.reg_mode)

        def toleranceFunc():
            a = float(self.tolerance.text())
            self.project.param['tolerance'] = a
            self.writeLog('k.param["tolerance"] = {:.2e}'.format(a))
        self.toleranceLabel = QLabel('<a href="tolerance">Value for tolerance</a>:')
        self.toleranceLabel.linkActivated.connect(showHelp)
        self.tolerance = QLineEdit()
        self.tolerance.setValidator(QDoubleValidator())
        self.tolerance.setText('1.0')
        self.tolerance.editingFinished.connect(toleranceFunc)
        invForm.addRow(self.toleranceLabel, self.tolerance)

        def max_iterationsFunc():
            a = int(self.max_iterations.text())
            self.project.param['max_iter'] = a
            self.writeLog('k.param["max_iter"] = {:d}'.format(a))
        self.max_iterationsLabel = QLabel('<a href="max_iterations">Maximum number of iterations</a>:')
        self.max_iterationsLabel.linkActivated.connect(showHelp)
        self.max_iterations = QLineEdit()
        self.max_iterations.setValidator(QIntValidator())
        self.max_iterations.setText('10')
        self.max_iterations.editingFinished.connect(max_iterationsFunc)
        invForm.addRow(self.max_iterationsLabel, self.max_iterations)

        def error_modFunc(index):
            self.project.param['error_mod'] = index
            self.writeLog('k.param["error_mod"] = {:d}'.format(index))
        self.error_modLabel = QLabel('<a href="error_mod">Update the weights</a>:')
        self.error_modLabel.linkActivated.connect(showHelpAdv)
        self.error_mod = QComboBox()
        self.error_mod.addItem('Keep the same weights [0]')
        self.error_mod.addItem('Update the weights [1]')
        self.error_mod.addItem('Update the weights (recommended) [2]')
        self.error_mod.setCurrentIndex(2)
        self.error_mod.activated.connect(error_modFunc)
        advForm.addRow(self.error_modLabel, self.error_mod)

        def alpha_anisoFunc():
            a = float(self.alpha_aniso.text())
            self.project.param['alpha_aniso'] = a
            self.writeLog('k.param["alpha_aniso"] = {:.2e}'.format(a))
        self.alpha_anisoLabel = QLabel('<a href="alpha_aniso">Value for <code>alpha_aniso</code></a>:')
        self.alpha_anisoLabel.linkActivated.connect(showHelpAdv)
        self.alpha_aniso = QLineEdit()
        self.alpha_aniso.setValidator(QDoubleValidator())
        self.alpha_aniso.setText('1.0')
        self.alpha_aniso.editingFinished.connect(alpha_anisoFunc)
        advForm.addRow(self.alpha_anisoLabel, self.alpha_aniso)

        def alpha_sFunc():
            a = float(self.alpha_s.text())
            self.project.param['alpha_s'] = a
            self.writeLog('k.param["alpha_s"] = {:.2e}'.format(a))
        self.alpha_sLabel = QLabel('<a href="alpha_s"><code>alpha_s</code></a>:')
        self.alpha_sLabel.linkActivated.connect(showHelpAdv)
        self.alpha_s = QLineEdit()
        self.alpha_s.setValidator(QDoubleValidator())
        self.alpha_s.setText('1.0')
        self.alpha_s.editingFinished.connect(alpha_sFunc)
        advForm.addRow(self.alpha_sLabel, self.alpha_s)
        self.alpha_sLabel.setVisible(False)
        self.alpha_s.setVisible(False)

        def min_errorFunc():
            a = float(self.min_error.text())
            self.project.param['min_error'] = a
            self.writeLog('k.param["min_error"] = {:f}'.format(a))
        self.min_errorLabel = QLabel('<a href="errorParam"><code>min_error</code></a>:')
        self.min_errorLabel.linkActivated.connect(showHelp)
        self.min_errorLabel.setVisible(False)
        self.min_error = QLineEdit()
        self.min_error.setText('0.01')
        self.min_error.editingFinished.connect(min_errorFunc)
        self.min_error.setVisible(False)
        invForm.addRow(self.min_errorLabel, self.min_error)

        def a_wgtFunc():
            a = float(self.a_wgt.text())
            self.project.param['a_wgt'] = a
            self.writeLog('k.param["a_wgt"] = {:f}'.format(a))
        self.a_wgtLabel = QLabel('<a href="errorParam"><code>a_wgt</code></a>:')
        self.a_wgtLabel.linkActivated.connect(showHelp)
        self.a_wgt = QLineEdit()
        self.a_wgt.setValidator(QDoubleValidator())
        self.a_wgt.setText('0.01')
        self.a_wgt.editingFinished.connect(a_wgtFunc)
        invForm.addRow(self.a_wgtLabel, self.a_wgt)

        def b_wgtFunc():
            a = float(self.b_wgt.text())
            self.project.param['b_wgt'] = a
            self.writeLog('k.param["b_wgt"] = {:f}'.format(a))
        self.b_wgtLabel = QLabel('<a href="errorParam"><code>b_wgt</code></a>:')
        self.b_wgtLabel.linkActivated.connect(showHelp)
        self.b_wgt = QLineEdit()
        self.b_wgt.setValidator(QDoubleValidator())
        self.b_wgt.setText('0.02')
        self.b_wgt.editingFinished.connect(b_wgtFunc)
        invForm.addRow(self.b_wgtLabel, self.b_wgt)

        def rho_minFunc():
            a = float(self.rho_min.text())
            self.project.param['rho_min'] = a
            self.writeLog('k.param["rho_min"] = {:f}'.format(a))
        self.rho_minLabel = QLabel('<a href="rho_max">Minimum apparent resistivity</a>:')
        self.rho_minLabel.linkActivated.connect(showHelp)
        self.rho_min = QLineEdit()
        self.rho_min.setValidator(QDoubleValidator())
        self.rho_min.setText('-1000')
        self.rho_min.editingFinished.connect(rho_minFunc)
        invForm.addRow(self.rho_minLabel, self.rho_min)

        def rho_maxFunc():
            a = float(self.rho_max.text())
            self.project.param['rho_max'] = a
            self.writeLog('k.param["rho_max"] = {:f}'.format(a))
        self.rho_maxLabel = QLabel('<a href="rho_max">Maximum apparent resistivity</a>:')
        self.rho_maxLabel.linkActivated.connect(showHelp)
        self.rho_max = QLineEdit()
        self.rho_max.setValidator(QDoubleValidator())
        self.rho_max.setText('1000')
        self.rho_max.editingFinished.connect(rho_maxFunc)
        invForm.addRow(self.rho_maxLabel, self.rho_max)
        

        generalLayout.addLayout(invForm)

        self.helpSection = QTextEdit('Help will be display here')
        self.helpSection.setReadOnly(True)
        self.helpSection.setText('Click on the labels and help will be displayed here')
        generalLayout.addWidget(self.helpSection)

        self.generalSettings .setLayout(generalLayout)
        self.tabInversionSettings .addTab(self.generalSettings , 'General')


        advancedLayout.addLayout(advForm)

        self.helpSection2 = QTextEdit('Help will be display here')
        self.helpSection2.setReadOnly(True)
        self.helpSection2.setText('Click on the labels and help will be displayed here')
        advancedLayout.addWidget(self.helpSection2)

        self.advancedSettings .setLayout(advancedLayout)
        self.tabInversionSettings .addTab(self.advancedSettings , 'Advanced')


        #%% tab 5 INVERSION
        self.tabInversion = QWidget()
        self.tabs.addTab(self.tabInversion, '&Inversion')
        self.tabs.setTabEnabled(5, False)
        
        self.invtabs = QTabWidget()
        self.logTab = QWidget()
        self.showTab = QWidget()
        self.computeTab = QWidget()
        self.invtabs.addTab(self.logTab, 'Log')
        self.invtabs.addTab(self.showTab, 'Results')
        self.invtabs.addTab(self.computeTab, 'Compute attribute')
        self.invtabs.setTabEnabled(1, False)
        self.invtabs.setTabEnabled(2, False)
                
        
        def frozeUI(val=True): # when inversion is running
            n = self.tabs.count()
            if val == True: # froze them
                self.tabState = [self.tabs.isTabEnabled(i) for i in range(n)]
                for i in range(n):
                    if i != 5:
                        self.tabs.setTabEnabled(i, False)
                self.hamBtn.setEnabled(False)
            else: # unfrozing
                for i in range(n):
                    self.tabs.setTabEnabled(i, self.tabState[i])
                if self.end is True:
                    self.tabs.setTabEnabled(6, True) # post processing tab should only be activated after successful inversion
                else:
                    self.tabs.setTabEnabled(6, False)
                self.hamBtn.setEnabled(True)
                    
        def afterInversion():
            # replace the log by the R2.out
            with open(os.path.join(self.project.dirname, self.project.typ + '.out'),'r') as f:
                text = f.read()
            self.logText.setText(text)
            self.project.proc = None
            
            # check if we don't have a fatal error
            if 'FATAL' in text and not (self.iBatch or self.iTimeLapse):
                self.end = False
                self.errorDump('WARNING: Error weights too high! Use lower <b>a_wgt</b> and <b>b_wgt</b> or choose an error model.')
            
            # if fixed elements are present, the mesh will be automatically
            # sorted at writing time, meaning we need to replot it
            # in the mesh tab
            if any(self.project.mesh.df['param'] == 0):
                if (self.project.typ == 'R3t') | (self.typ == 'cR3t'):
                    if pvfound:
                        self.mesh3Dplotter.clear() # clear all actors 
                        self.project.showMesh(ax=self.mesh3Dplotter, color_map='Greys', color_bar=False)
                    else:
                        self.mwMesh3D.plot(self.project.showMesh, threed=True)
                else:
                    replotMesh() # 2D only
                    
        def afterInversionShow():
            # show results
            if self.end is True:
                self.displayInvertedResults()
                try:
                    prepareInvError() # plot error graphs
                    if self.parallelCheck.isChecked():
                        self.logText.setText(self.project.invLog) # populate inversion log correctly after parallel inversion
                except Exception as e:
                    pdebug('ERROR: ui: invertBtnFunc:' + str(e))
                    self.errorDump(e)
                    pass
            self.invertBtn.setText('Invert')
            self.invertBtn.setStyleSheet('background-color:green; color:black')
            frozeUI(False)
            
        def afterInversionPseudo():
            self.logText.setText(self.pseudo3DInvtext)
            if self.project.proc.killFlag is True: # make sure everything stops now!
                self.end = False
            
            
        # ------------------------ log sub tab
        def invertBtnFunc():
            self.end = False # will be set to True if inversion successfull
            self.invtabs.setCurrentIndex(0) # log tab
            self.invtabs.setTabEnabled(1, False)
            self.invtabs.setTabEnabled(2, False)
            
            # check if we kill or invert
            if self.invertBtn.text() == 'Invert':
                self.invertBtn.setStyleSheet('background-color:red; color:black')
                self.invertBtn.setText('Kill')
                frozeUI(True)
            else:
                print('killing...', end='')
                # self.loadingWidget('Killing in progress...')
                if self.project.proc is not None:
                    self.project.proc.kill()
                print('done')
                self.invertBtn.setStyleSheet('background-color:green; color:black')
                self.invertBtn.setText('Invert')
                # self.loadingWidget(exitflag=True)
                frozeUI(False)
                return
            
            # clear stuff and connect
            self.logText.setText('')
            self.mwRMS.clear()
            self.mwIter.clear()
            self.mwInv.clear()
            
            # reset variables for RMS graph
            self.pindex = 0
            self.rms = []
            self.rmsIndex = []
            self.rmsIP = []
            self.rmsIndexIP = []
            self.inversionOutput = ''
            
            # create default mesh is not specified
            if self.project.mesh is None:
                logTextFunc('Creating the mesh... ')
                if (self.project.typ == 'R2') | (self.project.typ == 'cR2'):
                    meshTrianFunc()
                else:
                    elec = self.elecTable.getTable()
                    if self.topoTable.useNarray:
                        topo = self.topoTable.xyz      
                    else:
                        topo = self.topoTable.getTable()[['x','y','z']].values
                    inan = ~np.isnan(topo[:,0])

                    if self.project.elec['x'].isna().sum() > 0:
                        self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                        self.invertBtn.setStyleSheet('background-color: green; color:black')
                        self.invertBtn.setText('Invert')
                        frozeUI(False)
                        return
                    elif all(elec['y'].values == 0) & all(topo[inan,1] == 0):
                        self.errorDump('For 3D meshes, Y coordinates must be supplied for topo or elec at least.')
                        self.invertBtn.setStyleSheet('background-color:green; color:black')
                        self.invertBtn.setText('Invert')
                        frozeUI(False)
                        return
                    meshTetraFunc()
                logTextFunc('done!\n')

            # don't crop the mesh if that's what we've chosen
            if self.iCropping is True:
                if self.num_xz_poly is not None:
                    self.project.param['num_xz_poly'] = self.num_xz_poly
                    self.project.param['num_xy_poly'] = self.num_xz_poly
            else:
                self.project.param['num_xz_poly'] = 0
                self.project.param['num_xy_poly'] = 0
                self.project.param['zmin'] = np.inf 
                self.project.param['zmax'] = np.inf 
                
            # check how many cores to use 
            if self.ncoresText.text() == 'all':
                ncores = None # defaults to all the cores 
            else:
                ncores = int(self.ncoresText.text())
                if ncores < 1: # then fallback to using one core 
                   ncores = 1
                   
            if self.pseudo3DCheck.isChecked(): # pseudo 3D (multiple Projects)
                self.pseudo3DInvtext = ''
                def pseudo3DInvLog(dirname, typ):
                    with open(os.path.join(dirname, typ + '.out'),'r') as f:
                        self.pseudo3DInvtext += f.read() + '\n'
                        
                invLog = pseudo3DInvLog if self.parallelCheck.isChecked() is False else logTextFunc # handling dump
            
                # in thread to not block the UI
                kwargs = {'dump': logTextFunc,
                          'iplot': False, 
                          'runParallel': self.parallelCheck.isChecked(),
                          'invLog': invLog,
                          'ncores': ncores
                          }
                
                # run inversion in different thread to not block the UI
                self.thread = QThread()
                self.worker = Worker(self.project, kwargs, pseudo=True)
                self.worker.moveToThread(self.thread)
                self.thread.started.connect(self.worker.run)
                self.worker.finished.connect(self.thread.quit)
                self.worker.finished.connect(self.worker.deleteLater)
                self.thread.finished.connect(self.thread.deleteLater)
                self.worker.progress.connect(logTextFunc)
                self.thread.start()
                
                # self.project.invertPseudo3D(iplot=False, runParallel=self.parallelCheck.isChecked(),
                #                             dump=logTextFunc, invLog=invLog)
                self.writeLog('k.invertPseudo3D(iplot=False, runParallel={:s}, invLog={:s})'.format(
                    str(kwargs['runParallel']), str(invLog)))
                
                self.thread.finished.connect(afterInversionPseudo)
                self.thread.finished.connect(afterInversionShow)
                
                
            else: # regular single Project (2D or 3D)
                # set initial model
                x, phase0, zones, fixed = self.regionTable.getTable()
                regid = np.arange(len(x)) + 1 # 1 is the background (no 0)
                self.project.setStartingRes(dict(zip(regid, x)),
                                       dict(zip(regid, zones)),
                                       dict(zip(regid, fixed)),
                                       dict(zip(regid, phase0)))
                self.writeLog('k.setStartingRes(regionValues={:s}, zoneValues={:s},'
                              ' fixedValues={:s}, ipValues={:s}) '
                              '# define region manually using k.addRegion()'.format(
                                str(dict(zip(regid, x))),
                                str(dict(zip(regid, zones))),
                                str(dict(zip(regid, fixed))),
                                str(dict(zip(regid, phase0)))))
    
                # invert
                modErr = self.modErrCheck.isChecked()
                parallel = self.parallelCheck.isChecked()
                modelDOI = self.modelDOICheck.isChecked()
                kwargs = {'dump': logTextFunc, # keyword arguments for project.invert
                          'iplot': False,
                          'modErr': modErr, 
                          'parallel': parallel,
                          'modelDOI': modelDOI,
                          'ncores': ncores
                          }
                
                # run inversion in different thread to not block the UI
                self.thread = QThread()
                self.worker = Worker(self.project, kwargs)
                self.worker.moveToThread(self.thread)
                self.thread.started.connect(self.worker.run)
                self.worker.finished.connect(self.thread.quit)
                self.worker.finished.connect(self.worker.deleteLater)
                self.thread.finished.connect(self.thread.deleteLater)
                self.worker.progress.connect(logTextFunc)
                self.thread.start()

                # self.project.invert(iplot=False, dump=logTextFunc,
                #                     modErr=modErr, parallel=parallel, modelDOI=modelDOI)
                self.writeLog('k.invert(modErr={:s}, parallel={:s}, modelDOI={:s})'.format(
                    str(modErr), str(parallel), str(modelDOI)))
                
                self.thread.finished.connect(afterInversion)
                self.thread.finished.connect(afterInversionShow)

            
        self.invertBtn = QPushButton('Invert')
        self.invertBtn.setStyleSheet('background:green;')
        self.invertBtn.setToolTip('Invert all datasets')
        self.invertBtn.setAutoDefault(True)
        self.invertBtn.clicked.connect(invertBtnFunc)
                
        # RMS parser from live R2 output (can not use R2.showRMS())
        self.pindex = 0
        self.rms = []
        self.rmsIndex = []
        self.rmsIP = []
        self.rmsIndexIP = []
        self.inversionOutput = ''
        self.end = False

        
        def parseRMS(tt):
            newFlag = False
            a = tt.split()
            if len(a) > 0:
                if a[0] == 'Initial':
                    try:
                        newFlag = True
                        self.rms.append(float(a[3]))
                        self.rmsIndex.append(self.pindex)
                    except ValueError as e:
                        print('parseRMS error:', e)
                        pass
                    if a[1] == 'Phase':
                        self.rmsIP.append(float(a[4]))
                        self.rmsIndexIP.append(self.pindex)
                if a[0] == 'Processing':
                    self.pindex = self.pindex + 1
                if a[0] == 'End':
                    self.end = True
                if self.iBatch and self.parallelCheck.isChecked(): # parallel inversion of batch surveys doesn't show "End" rather a "x/x" format 
                    try:
                        spla0 = a[0].split('/')
                        if spla0[0] == spla0[-1]:
                            self.end = True
                    except:
                        pass
                if a[0] == 'Iteration': # if initial RMS shows solution then this will be "FATAL:" instead
                    if self.typ[-1] == 't':
                        cropMaxDepth = False # this parameter doesnt make sense for 3D surveys 
                    else:
                        cropMaxDepth = self.cropBelowFmd.isChecked()
                    if self.pseudo3DCheck.isChecked():
                        showIter = self.project.projectPseudo3D.showIter
                    else:
                        showIter = self.project.showIter
                        
                    self.mwIter.plot(partial(showIter, modelDOI=self.modelDOICheck.isChecked(),
                                             cropMaxDepth=cropMaxDepth), aspect=self.plotAspect)
                if a[0] == 'FATAL:':
#                    self.invertBtn.animateClick()
                    self.errorDump('WARNING: Error weights too high! Use lower <b>a_wgt</b> and <b>b_wgt</b> or choose an error model.')
#                    raise ValueError('*** Error weights too high! lower a_wgt and b_wgt or choose an error model. ***')
                    
            return newFlag

        def plotRMS(ax):
            if len(self.rms) > 0:
                rms = np.array(self.rms)
                rmsIndex = np.array(self.rmsIndex)
                xx = np.arange(len(rms))
                iindex = np.unique(rmsIndex)
                for ii in iindex:
                    ie = rmsIndex == ii
                    ax.plot(xx[ie], rms[ie], 'o-')
            if len(self.rmsIP) > 0:
                rms = np.array(self.rmsIP)
                rmsIndex = np.array(self.rmsIndexIP)
                xx = np.arange(len(rms))
                iindex = np.unique(rmsIndex)
                ax.set_prop_cycle(None)
                for ii in iindex:
                    ie = rmsIndex == ii
                    ax.plot(xx[ie], rms[ie], '*--')
            ax.set_xticks([])
            ax.set_xticklabels([])
            ax.set_xlabel('Iterations', fontsize=8)
            ax.tick_params(axis='both', which='major', labelsize=8)
            ax.set_ylabel('RMS', fontsize=8)
            
            
        def logTextFunc(arg):
            cursor = self.logText.textCursor()
            cursor.movePosition(cursor.End)
            cursor.insertText(arg)
            self.logText.ensureCursorVisible()
            if parseRMS(arg):
                self.mwRMS.plot(plotRMS)
            QApplication.processEvents()
        self.logText = QTextEdit()
        self.logText.setReadOnly(True)
        
        self.mwRMS = MatplotlibWidget(navi=False)
        self.mwIter = MatplotlibWidget(navi=False)
        
        
        # ------------------------- show subtab
        self.mwInv = MatplotlibWidget(navi=True)
        
        #interactive vtk widget with pyvista
        self.frame = QFrame()
        vlayout = QVBoxLayout()
        self.vtkWidget = QtInteractor(self.frame)
        vlayout.addWidget(self.vtkWidget.interactor)
        self.frame.setLayout(vlayout)
        self.mSlice = None
        self.mMesh = None
        
        self.displayParams = {'index':0,'edge_color':'none',
                            'sens':True, 'attr':'Resistivity(Ohm-m)',
                            'contour':False, 'vmin':None, 'vmax':None,
                            'cmap':'viridis', 'sensPrc':0.5,
                            'doi':self.modelDOICheck.isChecked(),
                            'doiSens':False, 'clipCorners':False,
                            'pvslices':([],[],[]), 'pvthreshold':None, 'pvdelaunay3d': False,
                            'pvgrid':False, 'pvcontour':[], 'aspect':'equal'}          
                
        # reset attribute specific settings (like vmin, vmax, threshold, isosurface)
        def resetAttributeSpecificSettings():
            self.displayParams['vmin'] = None
            self.displayParams['vmax'] = None
            self.vminEdit.setText('')
            self.vmaxEdit.setText('')
            self.pvthreshMax.setText('')
            self.pvthreshMin.setText('')
            self.displayParams['pvthreshold'] = None
            self.pvcontour.setText('')
            self.displayParams['pvcontour'] = []
        
        # change survey
        def surveyComboFunc(index):
            if self.pseudo3DCheck.isChecked():
                if index == 0:
                    self.showStackedLayout.setCurrentIndex(1)
                else:
                    self.showStackedLayout.setCurrentIndex(0)
                    index -= 1 # first graph is pseudo 3D
            self.displayParams['index'] = index
            attrs = list(self.project.meshResults[index].df.keys())
            attr0 = str(self.attrCombo.currentText())
            ci = 0
            c = -1
            found = False
            for i, a in enumerate(attrs): # find same attribute or plot first one
                if a not in ['param', 'region', 'zone', 'elm_id', 'cellType', 'X', 'Y', 'Z']:
                    c = c + 1
                    if a == attr0:
                        ci = c
                        attr = attr0
                        found = True
                        break
            if found is False: # we change attr so let's reset some parameters
                resetAttributeSpecificSettings()
                attr = attrs[4] # by default, first attribute
                ci = 0
            self.displayParams['attr'] = attr
            self.attrCombo.clear()
            for attr in attrs:
                if attr not in ['param', 'region', 'zone', 'elm_id', 'cellType', 'X', 'Y', 'Z']:
                    self.attrCombo.addItem(attr)
            self.attrCombo.setCurrentIndex(ci)
            self.replotSection()
        self.surveyCombo = QComboBox()
        self.surveyCombo.activated.connect(surveyComboFunc)
        self.surveyCombo.setToolTip('Change survey or see initial model.')
        
        # change attribute
        def attrComboFunc(index):
            resetAttributeSpecificSettings()
            self.displayParams['attr'] = str(self.attrCombo.currentText())
            if self.displayParams['attr'] == 'Sigma_imag(log10)':
                sigma_imag_vals = self.project.meshResults[self.displayParams['index']].df['Sigma_imag(log10)']
                if any(val == 0 for val in sigma_imag_vals):
                    if all(val == 0 for val in sigma_imag_vals):
                        pass
                    else:
                        self.contourCheck.setChecked(True)
                        self.infoDump('Contouring data by default!')
                        return # replotting triggers by edgeCheckFunc()
            self.replotSection()
        self.attrCombo = QComboBox()
        self.attrCombo.activated.connect(attrComboFunc)
        self.attrCombo.setToolTip('Choose which attribute to plot.')
        
        # change vmin/vmax
        self.vminLabel = QLabel('Min:')
        self.vminEdit = QLineEdit()
        self.vminEdit.setToolTip('Set mininum for current scale.')
        self.vminEdit.setValidator(QDoubleValidator())
        self.vmaxLabel = QLabel('Max:')
        self.vmaxEdit = QLineEdit()
        self.vmaxEdit.setToolTip('Set maximum for current color scale.')
        self.vmaxEdit.setValidator(QDoubleValidator())
        
        def vMinMaxBtnFunc():
            vmax = self.vmaxEdit.text()
            vmin = self.vminEdit.text()
            vmax = None if vmax == '' else float(vmax)
            vmin = None if vmin == '' else float(vmin)
            self.displayParams['vmin'] = vmin
            self.displayParams['vmax'] = vmax
            if (self.contourCheck.isChecked() is True) or (self.project.typ[-1] != '2'):
                self.replotSection()
            else:
                if self.pseudo3DCheck.isChecked() and self.surveyCombo.currentIndex() == 0:
                    self.replotSection()
                else:
                    self.mwInv.setMinMax(vmin=vmin, vmax=vmax)
        self.vMinMaxBtn = QPushButton('Apply')
        self.vMinMaxBtn.setAutoDefault(True)
        self.vMinMaxBtn.clicked.connect(vMinMaxBtnFunc)
        self.vMinMaxBtn.setToolTip('Apply limits on current color scale.')
        
        def doiCheckFunc(status):
            if status == Qt.Checked:
                self.displayParams['doi'] = True
            else:
                self.displayParams['doi'] = False
            self.replotSection()
        self.doiCheck = QCheckBox('DOI')
        self.doiCheck.setVisible(False)
        self.doiCheck.stateChanged.connect(doiCheckFunc)
        self.doiCheck.setToolTip('Depth of Investigation (DOI) based on Oldenburg and Li 1999.')
        
        def doiSensCheckFunc(status):
            if status == Qt.Checked:
                self.displayParams['doiSens'] = True
            else:
                self.displayParams['doiSens'] = False
            self.replotSection()
        self.doiSensCheck = QCheckBox('DOI estimate')
        self.doiSensCheck.stateChanged.connect(doiSensCheckFunc)
        self.doiSensCheck.setToolTip('Depth of Investigation (DOI) estimated based on sensitivity.\nSee advanced inversion settings for DOI with the Oldengburg and Li method.')
        
        self.cmapComboLabel = QLabel('Colormap')
        cmaps = ['viridis','viridis_r','plasma','plasma_r','seismic','seismic_r','cividis','cividis_r','winter',
                 'winter_r','autumn','autumn_r','rainbow','rainbow_r','jet','jet_r','Spectral','Spectral_r']
        def cmapComboFunc(index):
            self.displayParams['cmap'] = cmaps[index]
            self.replotSection()
        self.cmapCombo = QComboBox()
        for cmap in cmaps:
            self.cmapCombo.addItem(cmap)
        self.cmapCombo.setCurrentIndex(0)
        self.cmapCombo.activated.connect(cmapComboFunc)
        
        def contourCheckFunc(state):
            if state == Qt.Checked:
                self.edgeCheck.setEnabled(False)
                self.displayParams['contour'] = True
            else:
                self.edgeCheck.setEnabled(True)
                self.displayParams['contour'] = False
            self.replotSection()
        self.contourCheck = QCheckBox('Contour')
        self.contourCheck.stateChanged.connect(contourCheckFunc)
        self.contourCheck.setToolTip('Grid and contour the data.')
        
        def clipCornersFunc(state):
            if state == Qt.Checked:
                self.displayParams['clipCorners'] = True
            else:
                self.displayParams['clipCorners'] = False
            self.replotSection()
        self.clipCornersCheck = QCheckBox('Crop Corners')
        self.clipCornersCheck.stateChanged.connect(clipCornersFunc)
        self.clipCornersCheck.setToolTip('Triangles from bottom corners will be cropped (only if the whole mesh is not shown).')
      
        def sensSliderFunc(val):
            if val == 0:
                self.displayParams['sens'] = False
                self.infoDump('No sensitivity overlay displayed!')
            else:
                self.displayParams['sens'] = True
                val = (val-1)/10.0
                aa = np.logspace(-6, -1, 101)
                a = aa[int(val*100)]
                self.infoDump('Overlay sensitivity threshold is set to: </i>{:.1E} X max_sensitivity<i>'.format((a))) #to remove some ambiguity with this slider for people!
                self.displayParams['sensPrc'] = val
            self.replotSection()
        self.sensWidget = QWidget()
        sensLayout = QVBoxLayout()
        sensLayout.setContentsMargins(0,9,0,9)
        sensLayout.setSpacing(2)
        self.sensSlider = QSlider(Qt.Horizontal)
        self.sensSlider.setFixedWidth(110)
        self.sensSlider.setMinimum(0)
        self.sensSlider.setMaximum(11)
        self.sensSlider.setValue(5)
        self.sensLabel = QLabel('Sensitivity overlay')
        self.sensLabel.setAlignment(Qt.AlignCenter )
        sensLayout.addWidget(self.sensLabel)
        self.sensSlider.setToolTip('Normalized sensivity threshold')
        self.sensSlider.valueChanged.connect(sensSliderFunc)
        sensLayout.addWidget(self.sensSlider)
        self.sensWidget.setLayout(sensLayout)
        
        def showEdges(status):
            if status == Qt.Checked:
                self.displayParams['edge_color'] = 'k'
            else:
                self.displayParams['edge_color'] = 'none'
            self.replotSection()
        self.edgeCheck= QCheckBox('Edges')
        self.edgeCheck.setChecked(False)
        self.edgeCheck.setToolTip('Show edges of each mesh cell.')
        self.edgeCheck.stateChanged.connect(showEdges)
        
        def aspectCheckFunc(state):
            if state == Qt.Checked:
                self.displayParams['aspect'] = 'equal'
            else:
                self.displayParams['aspect'] = 'auto'
            self.replotSection()
        self.aspectCheck = QCheckBox('Equal')
        self.aspectCheck.setChecked(True)
        self.aspectCheck.stateChanged.connect(aspectCheckFunc)
        self.aspectCheck.setToolTip('Check for equal aspect of the axis'
                                    '\nUncheck for auto aspect.')

        def saveBtnFunc():
            fdir = QFileDialog.getExistingDirectory(self.tabImportingData, 'Choose the directory to export graphs and .vtk', directory=self.datadir)
            if fdir != '':
                self.loadingWidget('Saving data...')
                if self.project.typ[-1] == '2':
                    edge_color = self.displayParams['edge_color']
                    sens = self.displayParams['sens']
                    sensPrc = self.displayParams['sensPrc']
                    attr = self.displayParams['attr']
                    contour = self.displayParams['contour']
                    clipCorners = self.displayParams['clipCorners']
                    vmin = self.displayParams['vmin']
                    vmax = self.displayParams['vmax']
                    cmap = self.displayParams['cmap']
                    self.project.saveInvPlots(outputdir=fdir, edge_color=edge_color,
                                       contour=contour, sens=sens, attr=attr,
                                       vmin=vmin, vmax=vmax, color_map=cmap,
                                       sensPrc=sensPrc, clipCorners=clipCorners)
                    self.writeLog('k.saveInvPlots(outputdir="{:s}", edge_color="{:s}",'
                                  ' contour={:s}, sens={:s}, attr="{:s}", vmin={:s},'
                                  ' vmax={:s}, color_map="{:s}", sensPrc={:.2f})'.format(
                                  fdir, edge_color, str(contour), str(sens), attr,
                                  str(vmin), str(vmax), cmap, sensPrc))
                    if self.pseudo3DCheck.isChecked():
                        fname = os.path.join(fdir, 'Pseudo_3D_result.png')
                        self.vtkWidget.screenshot(fname, transparent_background=True)
                self.project.saveVtks(fdir)
                self.writeLog('k.saveVtks("{:s}")'.format(fdir))
            self.project.saveData(fdir)
            self.writeLog('k.saveData("{:s}")'.format(fdir))
            self.loadingWidget(exitflag=True)
            self.infoDump('All data and graphs saved successfully!')
        self.saveBtn = QPushButton('Save Data')
        self.saveBtn.clicked.connect(saveBtnFunc)
        self.saveBtn.setToolTip('Save current graph to the working directory.')
        
        
        # 3D specific options for pyvista
        self.pvthreshLabel = QLabel('Threshold:')
        self.pvthreshLabel.setToolTip('Value which to keep the cells.')
        self.pvthreshMin = QLineEdit('')
        self.pvthreshMin.setPlaceholderText('Min')
        self.pvthreshMin.setValidator(QDoubleValidator())
        self.pvthreshMin.setToolTip('Minimal value above which to keep the cells.')
        
        # pvthreshMaxLabel = QLabel('Max Threshold:')
        self.pvthreshMax = QLineEdit('')
        self.pvthreshMax.setPlaceholderText('Max')
        self.pvthreshMax.setValidator(QDoubleValidator())
        self.pvthreshMax.setToolTip('Maximal value below which to keep the cells.')
        
        self.pvslicesLabel = QLabel('Axis slices:')
        self.pvslicesLabel.setToolTip('Slice the mesh normal to X, Y and/or Z axis. '
                                 'Set multiple slices on one axis by separating values with ","')
        self.pvxslices = QLineEdit('')
        self.pvxslices.setPlaceholderText('X [m]')
        self.pvxslices.setToolTip('e.g. 4, 5 to have to slices normal'
                                  'to X in 4 and 5.')
        
        # pvyslicesLabel = QLabel('Y slices:')
        self.pvyslices = QLineEdit('')
        self.pvyslices.setPlaceholderText('Y [m]')
        self.pvyslices.setToolTip('e.g. 4, 5 to have to slices normal'
                                  'to Y in 4 and 5.')
        
        # pvzslicesLabel = QLabel('Z slices:')
        self.pvzslices = QLineEdit('')
        self.pvzslices.setPlaceholderText('Z [m]')
        self.pvzslices.setToolTip('e.g. 4, 5 to have to slices normal'
                                  'to Z in 4 and 5.')
                
        self.pvcontourLabel = QLabel('Isosurfaces:')
        self.pvcontour = QLineEdit('')
        self.pvcontour.setToolTip('Values of isosurfaces (comma separated).')
        
        def pvapplyBtnFunc():
            threshMin = float(self.pvthreshMin.text()) if self.pvthreshMin.text() != '' else None
            threshMax = float(self.pvthreshMax.text()) if self.pvthreshMax.text() != '' else None
            if threshMin is None and threshMax is None:
                self.displayParams['pvthreshold'] = None
            else:
                self.displayParams['pvthreshold'] = [threshMin, threshMax]
            if self.pvxslices.text() != '':
                xslices = [float(a) for a in self.pvxslices.text().split(',')]
            else:
                xslices = []
            if self.pvyslices.text() != '':
                yslices = [float(a) for a in self.pvyslices.text().split(',')]
            else:
                yslices = []
            if self.pvzslices.text() != '':
                zslices = [float(a) for a in self.pvzslices.text().split(',')]
            else:
                zslices = []
            self.displayParams['pvslices'] = (xslices, yslices, zslices)
            if self.pvcontour.text() != '':
                pvcontour = [float(a) for a in self.pvcontour.text().split(',')]
            else:
                pvcontour = []
            self.displayParams['pvcontour'] = pvcontour
            self.replotSection()
            
        def pvgridCheckFunc(state):
            self.displayParams['pvgrid'] = self.pvgridCheck.isChecked()
            try:
                self.replotSection()
            except:
                pass
        self.pvgridCheck = QCheckBox('Grid')
        self.pvgridCheck.stateChanged.connect(pvgridCheckFunc)
        
        def pvdelaunay3dCheckFunc(state):
            self.displayParams['pvdelaunay3d'] = self.pvdelaunay3dCheck.isChecked()
            try:
                self.replotSection()
            except:
                pass
        self.pvdelaunay3dCheck = QCheckBox('Delaunay 3D')
        self.pvdelaunay3dCheck.setToolTip('Apply a Delaunay 3D filter on the plot - similar to 2D contouring.')
        self.pvdelaunay3dCheck.stateChanged.connect(pvdelaunay3dCheckFunc)
            
        self.pvapplyBtn = QPushButton('Apply 3D')
        self.pvapplyBtn.setAutoDefault(True)
        self.pvapplyBtn.clicked.connect(pvapplyBtnFunc)
        self.pvapplyBtn.setToolTip('Apply 3D options')
        
        def screenshotBtnFunc():
            fname, _ = QFileDialog.getSaveFileName(self, 'Open File', self.datadir,
                                                   'PNG (*.png);;TIFF (*.tif);;JPEG (*.jpg)')
            if fname != '':
                if fname[-3:] == 'jpg':
                    transparent_background=False
                else: 
                    transparent_background=True
                self.vtkWidget.screenshot(fname, transparent_background=transparent_background)
        self.screenshotBtn = QPushButton('Save screenshot')
        self.screenshotBtn.setAutoDefault(True)
        self.screenshotBtn.clicked.connect(screenshotBtnFunc)
            
        # opt3d = [pvthreshMinLabel, self.pvthreshMin,
        #          pvthreshMaxLabel, self.pvthreshMax,
        #          pvxslicesLabel, self.pvxslices,
        #          pvyslicesLabel, self.pvyslices,
        #          pvzslicesLabel, self.pvzslices,
        #          self.pvcontourLabel, self.pvcontour,
        #          self.pvapplyBtn,
        #          self.pvdelaunay3dCheck, self.pvgridCheck,
        #          self.screenshotBtn]
        
        opt3d = [self.pvthreshLabel, self.pvthreshMin,
                 self.pvthreshMax, self.pvslicesLabel, 
                 self.pvxslices, self.pvyslices,
                 self.pvzslices, self.pvcontourLabel, 
                 self.pvcontour, self.pvapplyBtn,
                 self.pvdelaunay3dCheck, self.pvgridCheck, self.screenshotBtn]
        
        def show3DInvOptions(a):
            [o.setVisible(a) for o in opt3d]
        show3DInvOptions(False)
        
        
        
        # subtab compute attribute
        self.evalLabel = QLabel("You can use a formula to compute new attribute. "
                                "Use x['attributeName'] to access current attribute. "
                                "Once computed, go back to 'Results' tab to view the new attribute.")
        self.evalLabel.setWordWrap(True)
        self.evalEdit = QLineEdit()
        self.evalEdit.setPlaceholderText("formula: e.g. 1/x['Resistivity']")
        self.evalNewName = QLineEdit()
        self.evalNewName.setPlaceholderText('New attribute name')
        
        def evalBtnFunc():
            self.evalLog.clear()
            self.project.computeAttribute(formula=self.evalEdit.text(),
                                     name=self.evalNewName.text(),
                                     dump=evalLogFunc)
            # udpate attribute combo box
            index = self.surveyCombo.currentIndex()
            self.attrCombo.clear()
            attrs = self.project.meshResults[index].df.keys()
            c = 0
            ci = 0
            for i, attr in enumerate(attrs):
                if attr not in ['param', 'region', 'zone', 'elm_id', 'cellType', 'X', 'Y', 'Z']:
                    self.attrCombo.addItem(attr)
                    if attr == self.displayParams['attr']:
                        ci = c
                    c += 1
            self.attrCombo.setCurrentIndex(ci)
                
        self.evalBtn = QPushButton('Compute')
        self.evalBtn.clicked.connect(evalBtnFunc)
        
        def evalLogFunc(arg):
            cursor = self.evalLog.textCursor()
            cursor.movePosition(cursor.End)
            cursor.insertText(arg)
            self.evalLog.ensureCursorVisible()
            QApplication.processEvents()
        self.evalLog = QTextEdit()
        self.evalLog.setReadOnly(True)
        
        
        # layout
        invLayout = QVBoxLayout()
        logLayout = QHBoxLayout()
        showLayout = QVBoxLayout()
        computeLayout = QVBoxLayout()
        
        logLeftLayout = QVBoxLayout()
        logLeftLayout.addWidget(self.invertBtn)
        logLeftLayout.addWidget(self.logText)
        logLayout.addLayout(logLeftLayout)
        
        logRightLayout = QVBoxLayout()
        logRightLayout.addWidget(self.mwRMS, 50)
        logRightLayout.addWidget(self.mwIter, 50)
        logLayout.addLayout(logRightLayout)
        
        optsLayout = QHBoxLayout()
        optsLayout.addWidget(self.surveyCombo, 15)
        optsLayout.addWidget(self.attrCombo, 20)
        optsLayout.addWidget(self.vminLabel)
        optsLayout.addWidget(self.vminEdit, 5)
        optsLayout.addWidget(self.vmaxLabel)
        optsLayout.addWidget(self.vmaxEdit, 5)
        optsLayout.addWidget(self.vMinMaxBtn)
        optsLayout.addWidget(self.doiCheck)
        optsLayout.addWidget(self.doiSensCheck)
        optsLayout.addWidget(self.cmapCombo)
        optsLayout.addWidget(self.contourCheck)
        optsLayout.addWidget(self.clipCornersCheck)
        optsLayout.addWidget(self.sensWidget)
        optsLayout.addWidget(self.edgeCheck)
        optsLayout.addWidget(self.aspectCheck)
        optsLayout.addWidget(self.saveBtn)
        
        opt3dLayout = QHBoxLayout()
        for o in opt3d:
            opt3dLayout.addWidget(o)
        
        showLayout.addLayout(optsLayout)
        showLayout.addLayout(opt3dLayout)
        self.showStackedLayout = QStackedLayout()
        self.showStackedLayout.addWidget(self.mwInv)
        self.showStackedLayout.addWidget(self.frame)
        showLayout.addLayout(self.showStackedLayout)
        
        computeLayout.addWidget(self.evalLabel)
        computeLayoutTop = QHBoxLayout()
        computeLayoutTop.addWidget(self.evalEdit)
        computeLayoutTop.addWidget(self.evalNewName)
        computeLayoutTop.addWidget(self.evalBtn)
        computeLayout.addLayout(computeLayoutTop)
        computeLayout.addWidget(self.evalLog)
        
        self.logTab.setLayout(logLayout)
        self.showTab.setLayout(showLayout)
        self.computeTab.setLayout(computeLayout)
        
        invLayout.addWidget(self.invtabs)
        self.tabInversion.setLayout(invLayout)
    
    

        #%% tab 6 POSTPROCESSING
        self.tabPostProcessing = QWidget()
        self.tabs.addTab(self.tabPostProcessing, 'Post-processing')
        self.tabs.setTabEnabled(6,False)
        
        self.errorGraphs = QTabWidget()
        
        self.invErrorTabLabel = QLabel('Enter error range to filter the data based on inversion error or '
                                       'manually remove points on the pseudo section and then reinvert the data.')
        
        self.invPseudoErrLabel = QLabel('Select datapoints on the pseudo section to remove.')
        
        def prepareInvError():
            names = [s.name for s in self.project.surveys]
            if len(names) > 1:
                self.invErrorComboLabel.show()
                self.invErrorCombo.show()
            else:
                self.invErrorComboLabel.hide()
                self.invErrorCombo.hide()
                
            if self.project.iTimeLapse:
                names[0] = names[0] + ' (Ref)'
                
            # self.invErrorCombo.disconnect()
            self.invErrorCombo.clear()
            self.invErrorCombo.addItem('Apply to Each')
            for name in names:
                self.invErrorCombo.addItem(name)
            self.invErrorCombo.activated.connect(invErrorComboFunc)
            invErrorComboFunc(1)
            self.invErrorIndex = 0
        
        def invErrorComboFunc(index):
            try:
                index = self.invErrorIndex if index == 0 else index-1 # for apply to each (time-lapse or batch only)
                plotInvError2(index)
                if self.iBorehole is False and self.m3DRadio.isChecked() is False:
                    self.errorGraphs.setTabEnabled(0, True)
                    # self.errorGraphs.setCurrentIndex(0)
                    plotInvError(index)
                    self.invErrorIndex = index
                else:
                    self.errorGraphs.setTabEnabled(0, False)
                    # self.errorGraphs.setCurrentIndex(1)
            except Exception as e:
                print('Could not print error: ', e)
        self.invErrorCombo = QComboBox()
        self.invErrorCombo.setMinimumWidth(250)
        self.invErrorCombo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.invErrorCombo.activated.connect(invErrorComboFunc)
        self.invErrorCombo.hide()
        self.invErrorComboLabel = QLabel('Dataset:')
        self.invErrorComboLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.invErrorComboLabel.hide()
        
        self.rangeInvErrorLabel = QLabel('Error range:')
        self.rangeInvErrorLabel.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        
        
        def filterInvError():
            if self.invErrorCombo.currentIndex() == 0: # apply to each
                index = -1
            else:
                index = self.invErrorIndex
            vminVal = self.rangeInvErrorMinInput.text()
            vmaxVal = self.rangeInvErrorMaxInput.text()
            vmin = float(vminVal) if vminVal != '' else None 
            vmax = float(vmaxVal) if vmaxVal != '' else None
            self.project.filterInvError(index=index, vmin=vmin, vmax=vmax)
            self.writeLog('k.filterInvError(index={:d}, vmin={}, vmax={})'.format(
            index, str(vmin), str(vmax)))    
            plotInvError(self.invErrorIndex)
            plotInvError2(self.invErrorIndex)
            
        def resetInvErrorFilter():
            if self.invErrorCombo.currentIndex() == 0: # apply to each
                for s in self.project.surveys:
                    s.df = s.dfInvErrOutputOrigin.copy()
            else:
                self.project.surveys[self.invErrorIndex].df = self.project.surveys[self.invErrorIndex].dfInvErrOutputOrigin.copy()
            
            self.rangeInvErrorMinInput.setText('')
            self.rangeInvErrorMaxInput.setText('')
            plotInvError(self.invErrorIndex)
            plotInvError2(self.invErrorIndex)

        
        self.rangeInvErrorMinInput = QLineEdit('')
        self.rangeInvErrorMinInput.setPlaceholderText('min')
        self.rangeInvErrorMinInput.setToolTip('Minimum accepted normalized error')
        self.rangeInvErrorMinInput.setFixedWidth(80)
        self.rangeInvErrorMinInput.setValidator(QDoubleValidator())
        
        self.rangeInvErrorMaxInput = QLineEdit('')
        self.rangeInvErrorMaxInput.setPlaceholderText('max')
        self.rangeInvErrorMaxInput.setToolTip('Maximum accepted normalized error')
        self.rangeInvErrorMaxInput.setFixedWidth(80)
        self.rangeInvErrorMaxInput.setValidator(QDoubleValidator())
        
        self.rangeInvErrorApplyBtn = QPushButton('Apply filter')
        self.rangeInvErrorApplyBtn.setToolTip('Apply range filtering to selected dataset(s)')
        self.rangeInvErrorApplyBtn.setFixedWidth(150)
        self.rangeInvErrorApplyBtn.clicked.connect(filterInvError)
        
        self.rangeInvErrorResetBtn = QPushButton('Reset filters')
        self.rangeInvErrorResetBtn.setStyleSheet("color: red")
        self.rangeInvErrorResetBtn.setToolTip('Reset plots and filters to after inversion state.')
        self.rangeInvErrorResetBtn.setFixedWidth(150)
        self.rangeInvErrorResetBtn.clicked.connect(resetInvErrorFilter)
        
        
        def plotInvError(index=0):
            pdebug('plotInvError()')
            self.mwInvErrorIP.setVisible(False) # let's show/hide it each time we invert
            self.mwInvError.setCallback(self.project.showPseudoInvError)
            self.mwInvError.replot(index=index)
            self.writeLog('k.showPseudoInvError(index={:d})'.format(index))
            if self.project.typ[0] == 'c':
                self.mwInvErrorIP.setVisible(True)
                self.mwInvErrorIP.setCallback(self.project.showPseudoInvErrorIP)
                self.mwInvErrorIP.replot(index=index)
                self.writeLog('k.showPseudoInvErrorIP(index={:d})'.format(index))
            

        self.mwInvError = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        self.mwInvErrorIP = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        self.mwInvErrorIP.setVisible(False)
        
        def invErrorFiltFunc():
            try:
                i2remove = self.project.surveys[self.invErrorIndex].iselect
                if not all(self.project.surveys[self.invErrorIndex].df['irecip'].values == 0): 
                    # as the selection is only done on normal dataset, reciprocal pair should be removed too
                    recipPaires = self.project.surveys[self.invErrorIndex].df[i2remove]['irecip'].values*-1
                    ie = np.isin(self.project.surveys[self.invErrorIndex].df['irecip'].values, recipPaires[recipPaires != 0]) 
                    i2remove = i2remove + ie
                self.project.surveys[self.invErrorIndex].filterData(~i2remove)
                plotInvError(self.invErrorIndex)
                plotInvError2(self.invErrorIndex)

            except:
                self.errorDump('Error plotting error graphs.')
                pass
        
        self.invErrorFilterBtn = QPushButton('Remove selected')
        self.invErrorFilterBtn.setFixedWidth(150)
        self.invErrorFilterBtn.clicked.connect(invErrorFiltFunc)
        
        def invErrorReinvert():
            self.tabs.setCurrentIndex(5) # jump to inversion tab
            self.invertBtn.animateClick() # invert
        
        self.invErrorReinvertBtn = QPushButton('Invert')
        self.invErrorReinvertBtn.setFixedWidth(150)
        self.invErrorReinvertBtn.clicked.connect(invErrorReinvert)
        self.invErrorReinvertBtn.setStyleSheet('background-color:green; color:black')
      

        def plotInvError2(index=0):
            pdebug('plotInvError2()')
            self.mwInvError2.setCallback(self.project.showInvError)
            self.mwInvError2.replot(index=index)
            self.writeLog('k.showInvError(index={:d})'.format(index))
        self.mwInvError2 = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        invErrorLabel = QLabel('All errors should be between +/- 3% (Binley at al. 1995). '
                               'If it\'s not the case try to fit an error model or '
                               'manually change the a_wgt and b_wgt in inversion settings.')


        # layout
        postProcessingLayout = QVBoxLayout()
        
        topInvErrorLayout = QHBoxLayout()
        # topInvErrorLayout.setAlignment(Qt.AlignRight)
        topInvErrorLayout.addWidget(self.invErrorTabLabel, 1)
        topInvErrorLayout.addWidget(self.invErrorComboLabel, 0)
        topInvErrorLayout.addWidget(self.invErrorCombo, 0)
        postProcessingLayout.addLayout(topInvErrorLayout)
        
        rangeInvErrorLayout = QHBoxLayout()
        rangeInvErrorLayoutL = QHBoxLayout()
        rangeInvErrorLayoutL.setAlignment(Qt.AlignLeft)
        rangeInvErrorLayoutL.addWidget(self.rangeInvErrorLabel)
        rangeInvErrorLayoutL.addWidget(self.rangeInvErrorMinInput)
        rangeInvErrorLayoutL.addWidget(self.rangeInvErrorMaxInput)
        rangeInvErrorLayout.addLayout(rangeInvErrorLayoutL, 0)
        
        rangeInvErrorLayoutR = QHBoxLayout()
        rangeInvErrorLayoutR.setAlignment(Qt.AlignRight)
        rangeInvErrorLayoutR.addWidget(self.rangeInvErrorApplyBtn)
        rangeInvErrorLayoutR.addWidget(self.rangeInvErrorResetBtn)
        rangeInvErrorLayoutR.addWidget(self.invErrorReinvertBtn)
        rangeInvErrorLayout.addLayout(rangeInvErrorLayoutR, 1)
        
        postProcessingLayout.addLayout(rangeInvErrorLayout)
        
        postProcessingLayout.addWidget(self.errorGraphs)
        
        self.invError = QWidget()
        self.errorGraphs.addTab(self.invError, 'Pseudo Section of Inversion Errors')

        invErrorLayout = QVBoxLayout()
        invErrorLayout.setAlignment(Qt.AlignTop)
        
        invErrorTopLayout = QHBoxLayout()
        invErrorTopLayoutL = QHBoxLayout()
        invErrorTopLayoutL.addWidget(self.invPseudoErrLabel)
        invErrorTopLayout.addLayout(invErrorTopLayoutL)
        
        invErrorTopLayoutR = QHBoxLayout()
        invErrorTopLayoutR.setAlignment(Qt.AlignRight)
        invErrorTopLayoutR.addWidget(self.invErrorFilterBtn)

        invErrorTopLayout.addLayout(invErrorTopLayoutR)
        
        invErrorLayout.addLayout(invErrorTopLayout, 0)
        
        invErrorPlotLayout = QHBoxLayout()
        invErrorPlotLayout.addWidget(self.mwInvError, 50)
        invErrorPlotLayout.addWidget(self.mwInvErrorIP, 50)
        
        invErrorLayout.addLayout(invErrorPlotLayout, 1)
        self.invError.setLayout(invErrorLayout)

        self.invError2 = QWidget()
        self.errorGraphs.addTab(self.invError2, 'Normalised Inversion Errors')
        invErrorLayout2 = QVBoxLayout()
        invErrorLayout2Plot = QVBoxLayout()

        invErrorLayout2Plot.addWidget(self.mwInvError2, Qt.AlignCenter)
        invErrorLayout2.addLayout(invErrorLayout2Plot, 1)
        invErrorLayout2.addWidget(invErrorLabel)
        self.invError2.setLayout(invErrorLayout2)
        
        self.tabPostProcessing.setLayout(postProcessingLayout)



        #%% Help tab
        self.tabHelp = QTabWidget()
        self.tabs.addTab(self.tabHelp, 'Help')

        helpLayout = QVBoxLayout()
        helpLayout.setAlignment(Qt.AlignTop)
        helpText = QTextBrowser()
        helpText.setOpenExternalLinks(True)
        helpText.setText('''
           <h1>General help</h1>\
           <p>Below are simple instructions to guide you to through the software.</p>
           <ul>
           <p><li>In the "<b>Importing</b>" tab:
           <ul>
           <li>Select if you want a 2D/3D survey, an inverse/forward solution and check if you have borehole/timelapse/batch data.</li>
           <li>Modify the default working directory if you want to keep the outputted files afterwards.</li>
           <li>Select the file type. You can choose "Custom" if you file type is not available and you will be redirected to the custom parser tab.</li>
           <ul><li>Note: Syscal files must be exported as 'Spreadsheet' files with .csv format (comma separators) from Prosys.</li>
           <li>Note: Res2DInv files are mostly supported, but it is recommended to change them in "General Array" format if your file is not recognized.</li>
           <li>Note: Res3DInv files are only supported in general 4 electrode format.</ul></li>
           <ul>
           <li>If your survey has topography, you can import it in the "Electrodes(XZY/Topo)" tab.</li>
           <ul><li><i>Pole-dipole arrays</i>: the remote electrode's X location must be exactly at 99999 or -99999 m.</li>
           <li><i>Pole-pole arrays</i>: first remote electrode's X location must be exactly at 99999 and second one at exactly -99999 m (or vice versa).</ul></li>
           <li>Then one can choose to directly invert with all the default settings or go through the other tabs on the rights.</li>
           </ul></li></p>
           <p><li>In the "<b>Pre-processing</b>" tab:
           <ul>
           <li>The first tab offers manual filtering option based on reciprocal measurements in the dataset (if any).</li>
           <li>The "Phase Filtering" tab is only enable for IP data and allows precise filtering of IP data (range filtering, removing of nested measurements, DCA, ...).</li>
           <li>The "Resistance Error Model" tab allows to fit a power-law or linear error model to resistance data.</li>
           <li>The "Phase Error Model" tab allows to fit a power-law or parabolic error model to phase data.</li>
           <br>**<b>NOTE:</b> Error modeling tabs will be disabled if there are no reciprocal measurements in the imported dataset
           </ul></li></p>
           <p><li>In the "<b>Mesh</b>" tab you can create a quadrilateral or triangular mesh (2D) or a tetrahedral mesh (3D). For 2D mesh you can specify different\
           region of given resistivity/phase and if they need to be fixed or not during inversion. For forward modelling this mesh serves as the initial model.</li></p>
           <p><li>In the "<b>Forward model</b>" tab (only available in forward mode) you can design your sequence and add noise. The resulting synthetic measurements will be\
           automatically added to as an actual survey in ResIPy and can be inverted directly.</li></p>
           <p><li>In the "<b>Inversion Settings</b>" tab, you can modify all settings for the inversion. Help is available by clicking on the label of each item. The help\
           generally refers to the document present in the R2/cR3/R3t/cR3t respective manuals.</li></p>
           <p><li>In the "<b>Inversion</b>" tab, you can invert your survey and see the output in real time. if you have selected parallel inversion in "Inversion Settings">"Advanced",\
           then nothing will be printed out until the inversion finished. When the inversion finished you will be able to see the inverted section, open it with Paraview (mainly for 3D)\
           and save the outputted .vtk file and graphs using the "Save Graphs" button.</li>
           <ul><li>Plot aspect ratio can be changed by dragging  left handle to right or left and top handle (above plot options) up and down.</li></ul></p>
           <p><li>The "<b>Post-processing</b>" tab displays the errors from the inversion. It helps to assess the quality of the inversion.</li>
           </ul></p>
           <p><b>Figure options</b> (available under each figure, the button next to save icon):</p>
           <ul>
           <li>Select the axes to edit (usually "<i>Distance [m] - Elevation [m]</i>):
           <ul><li>The axes limits/labels and plot title  can be changed in the "Axes" tab.</li>
           <li>The marker size and line width can be changed in the "Curves" tab.</li></ul>
           </ul>
           <p><b>Pseudo 3D inversion<sup>beta</sup></b>:<br>
           In 2D mode, creates a batch survey from separate 2D lines and also visualize the lines in a 3D mode.</p>
           <ul>
           <li>Electrodes must be provided as if a 3D survey is imported.
           <ul><li>Electrodes must have all <i><b>label, x, y and z</i></b> values</li>
           <li>All electrode labels must <b><i>EXACTLY</b></i> have <font color="red">"Line Number [space] Electrode Number"</font></b> format.</li>
           <li>In case of borehole surveys, each line <b>must</b> have at least one surface electrode (i.e., whole mesh problems are not supported).</li></ul>
           </ul>
           <p><b>Output files</b> (Inversion/modeling results):</p>
           <ul>
           <li>Always select a <i>Working directory</i> (in the <i>Importing > Data</i> tab), so you do not lose your inversion/modeling results.
           <li>Data in X Y Z format is saved in <i>f001_res.dat</i> for inversions and <i>*_forward.dat</i> (*R2/cR2/R3t/cR3t) for forward models.
           <li>For more on output files of ResIPy, <u>please read the readme files of R2/cR2/R3t/cR3t (links below)</u>.
           </ul>
           <p><b>Reporting an issue</b>:</p>
           <ul>
           <li>Please make sure you include the full error (copy the whole error details) when reporting it on our <a href="https://gitlab.com/hkex/resipy/issues">GitLab page</a>.
           <li>Including examples that cause the error in the report will be very helpful to debug the software.
           </ul>
           ''')
        # local manuals
        R2_help = bytearray(QUrl.fromLocalFile(resource_path('resipy/exe/R2_manual.pdf')).toEncoded()).decode()
        R3t_help = bytearray(QUrl.fromLocalFile(resource_path('resipy/exe/R3t_manual.pdf')).toEncoded()).decode()
        cR2_help = bytearray(QUrl.fromLocalFile(resource_path('resipy/exe/cR2_manual.pdf')).toEncoded()).decode()
        cR3t_help = bytearray(QUrl.fromLocalFile(resource_path('resipy/exe/cR3t_manual.pdf')).toEncoded()).decode()
        helpLabel = QLabel( 
        '''
        <p>More help for ResIPy: <a href="https://hkex.gitlab.io/resipy/">Documentation and the examples</a></p>
        <p>Video tutorials: <a href="https://www.youtube.com/channel/UCkg2drwtfaVAo_Tuyeg_z5Q">ResIPy on YouTube</a></p>
        <p>Read more on 2D resistivity inversion: <a href="%s">R2_readme.pdf</a> | Read more on 3D resistivity inversion: <a href="%s">R3t_readme.pdf</a></p>
        <p>Read more on 2D complex resistivity (IP) inversion: <a href="%s">cR2_readme.pdf</a> | Read more on 3D complex resistivity (IP) inversion: <a href="%s">cR3t_readme.pdf</a></p>
        '''% (R2_help,R3t_help, cR2_help, cR3t_help))
        helpLabel.setOpenExternalLinks(True)
        helpLayout.addWidget(helpText)
        helpLayout.addWidget(helpLabel)
        self.tabHelp.setLayout(helpLayout)


        #%% About tab

        self.tabAbout = QTabWidget()
        self.tabs.addTab(self.tabAbout, 'About')

        infoLayout = QVBoxLayout()
        aboutText = QLabel()
        aboutText.setText(
            '''
            <h1>About ResIPy </h1>
            <p><b>Version: {:s}</b></p>
            <p><i>ResIPy is a free and open source software for inversion and modeling of geoelectrical data (Resistivity and IP)</i></p>
            <p>If you encouter any issues or would like to submit a feature request, please raise an issue on our gitlab repository at:</p>
            <p><a href="https://gitlab.com/hkex/resipy/issues">https://gitlab.com/hkex/resipy/issues</a></p>
            <p>ResIPy uses 
                <a href="{}">R2</a>,
                <a href="{}">cR2</a>,
                <a href="{}">R3t</a> and 
                <a href="{}">cR3t</a> developed by Andrew Binley</p>
            <p>For generation of triangular mesh, ResIPy uses software 
                <a href="http://gmsh.info/">Gmsh</a></p>
            <p>Python packages used: 
                <a href="https://numpy.org/">numpy</a>, 
                <a href="https://pandas.pydata.org/">pandas</a>,
                <a href="https://matplotlib.org/">matplotlib</a>,
                <a href="https://scipy.org/index.html">scipy</a>,
                <a href="https://pypi.org/project/PyQt5/">PyQt5</a>,
                <a href="https://pypi.org/project/pyvista/">PyVista</a>.
            </p>
            <p><strong>ResIPy's core developers: Guillaume Blanchy, Sina Saneiyan, Jimmy Boyd and Paul McLachlan.<strong></p>
            <p>Contributors: Pedro Concha, Michael Tso</p>
            <p><b><a href="https://www.researchgate.net/project/pyR2-GUI-for-R2-family-codes">Visit our ResearchGate page!</a></b></p>
            <p><b>Citing ResIPy</b>:<br>Blanchy G., Saneiyan S., Boyd J., McLachlan P.
            and Binley A. 2020.<br>“ResIPy, an Intuitive Open Source Software for 
            Complex Geoelectrical Inversion/Modeling.”<br>Computers & Geosciences, February, 104423.
            <a href="https://doi.org/10.1016/j.cageo.2020.104423">https://doi.org/10.1016/j.cageo.2020.104423</a>.</p>
            '''.format(ResIPy_version, R2_help, cR2_help, R3t_help, cR3t_help))
        aboutText.setOpenExternalLinks(True)
        aboutText.setWordWrap(True)
        aboutText.setAlignment(Qt.AlignTop | Qt.AlignHCenter)
        infoLayout.addWidget(aboutText, 0)

        self.tabAbout.setLayout(infoLayout)

        #%% test tab
        # tabTest = QTabWidget()
        # self.tabs.addTab(tabTest, 'TEST')
        # self.tabs.setCurrentIndex(9)
        
        # self.fframe = QFrame()
        # vlayout = QVBoxLayout()
        # self.vtk_widget = QtInteractor(self.fframe)
        # vlayout.addWidget(self.vtk_widget)
        # self.fframe.setLayout(vlayout)
        
        # m = pv.read('/home/jkl/Downloads/f001.vtk')
        # self.vtk_widget.show_grid()
        # self.vtk_widget.add_axes()
        # self.vtk_widget.add_mesh(m, scalars='Resistivity')
        
        
        # tabTestLayout = QVBoxLayout()
        # tabTestLayout.addWidget(self.fframe)
        # tabTest.setLayout(tabTestLayout)
        
        
        #%% general Ctrl+Q shortcut + general tab layout

        self.layout.addWidget(self.tabs)
        self.layout.addWidget(self.errorLabel)
        self.tableWidget.setLayout(self.layout)
        self.setCentralWidget(self.tableWidget)
        self.showMaximized() # maximizing window on start
        self.show()

#%% useful methods

    def keyPressEvent(self, e):
        if (e.modifiers() == Qt.ControlModifier) & (e.key() == Qt.Key_Q):
            self.close()
        if (e.modifiers() == Qt.ControlModifier) & (e.key() == Qt.Key_W):
            self.close()
            
    def timeOut(self, timeStamp):
        if resipySettings.param['dark'] == 'True':
            timeStmpCol = 'white'
        else:
            timeStmpCol = 'black'
        self.errorLabel.setText(('<i style="color:%s">[' % timeStmpCol) +timeStamp+']: </i>')
        
    def errorDump(self, text, flag=1):
        text = str(text)
        timeStamp = time.strftime('%H:%M:%S')
        if flag == 1: # error in red
            col = 'red'
            pdebug('errorDump:', text)
        else:
            if resipySettings.param['dark'] == 'True':
                col = 'white'
            else:
                col = 'black'
            pdebug('infoDump:', text)
        self.errorLabel.setText('<i style="color:'+col+'">['+timeStamp+']: '+text+'</i>')
        self.timer.timeout.connect(partial(self.timeOut, timeStamp))
        self.timer.start(10000) # 10 secs - making sure the error/message doen't stick forever!

    def infoDump(self, text):
        self.errorDump(text, flag=0)
        
    # def writeLog(self, text):
    #     """Write a new line to the log file.
    #     """
    #     flog = os.path.join(self.project.dirname, 'log.txt')
    #     try:
    #         if os.path.exists(flog) is False:
    #             with open(flog, 'w') as f:
    #                 f.write('# -------- ResIPy log file -------\n')
    #         with open(flog, 'a') as f:
    #             f.write(text + '\n')
    #     except Exception as e:
    #         print('Error in writeLog:', e)
    
    def writeLog(self, text, dico=None):
        if dico is None:
            self.apiLog += text + '\n'
        else:
            arg = ''
            for key in dico.keys():
                val = dico[key]
                if type(val) == str:
                    arg += key + '="{:s}", '.format(val)
                elif (type(val) == int or
                    type(val) == float or
                    type(val) == tuple or
                    type(val) == list or
                    type(val) == bool or
                    val is None):
                    arg += key + '={:s}, '.format(str(val))
                elif type(val) == type(np.ndarray):
                    arg += key + '={:s}, '.format(str(list(val)))
                else:
                    pass # ignore argument
            self.apiLog += text + '(' + arg[:-2] + ')\n'
        
    def saveLog(self):
        fname, _ = QFileDialog.getSaveFileName(self)
        if fname != '':
            if fname[-3:] != '.py':
                fname = fname + '.py'
            self.apiLog = self.apiLog.replace('"None"', 'None')
            with open(fname, 'w') as f:
                f.write(self.apiLog)
            
    def loadingWidget(self, msgtxt='', exitflag=False):
        '''Shows a dialog to indicate ResIPy is working (e.g., loading large dataset, etc.)'''
        if OS != 'Linux':
            if exitflag == False:
                self.loadingDialogTxtWidget.setText(msgtxt)
                self.loadingDialog.show()            
            else:
                self.loadingDialog.accept() 
            QApplication.processEvents()
            
    def calcAspectRatio(self): # for calculating aspect ratio of long surveys
        if self.project.fmd is None:
            self.project.computeFineMeshDepth()
        surLength = np.abs(self.project.param['xz_poly_table'][0,0] - self.project.param['xz_poly_table'][1,0])
        surDepth = np.abs(self.project.param['xz_poly_table'][-1,1] - self.project.param['xz_poly_table'][-2,1])
        aspectRatio = surLength/surDepth
        if not 0.2 < aspectRatio < 5: # make sure not to get narrow plots (max aspect ratio is 5 or 1/5)
            self.plotAspect = 'auto'
            
    
    def recipOrNoRecipShow(self, recipPresence=True): # changes the reciprocal filtering tab stuff
        if recipPresence == True:
            self.tabPreProcessing.setTabText(0, 'Reciprocal Filtering')
            self.recipErrorLabel.setText('<b>Remove datapoints that have reciprocal error larger than what you prefer.</b><br>Either select (<i>click on the dots to select them</i>) the points on the pseudo section below or choose a percentage threshold or both!</br>')
            # self.recipErrorInputLabel.show()
            # self.rhoRangeInputLabel.hide()
            # self.rhoRangeMinInput.hide()
            # self.rhoRangeMaxInput.hide()
            # self.recipErrorInputLine.show()
            # self.recipErrorInputLine.setText('')
            # self.rhoRangeMinInput.setText('')
            # self.rhoRangeMaxInput.setText('')
            self.recipErrorUnpairedBtn.show()
            self.recipErrorPltBtn.setToolTip('Removes measuremtns that have either greater reciprocal error than "Percent error threshold" or are manually selected or both!')
            self.recipErrorBottomTabs.setTabEnabled(1, True)
        if recipPresence == False:
            self.tabPreProcessing.setTabText(0, 'Manual Filtering')
            self.recipErrorLabel.setText('<b>Define a range of apparent resistivity to keep and/or select datapoints to remove.</b><br><i>click on the dots to select them</i></br>')
            # self.recipErrorInputLabel.hide()
            # self.rhoRangeInputLabel.show()
            # self.rhoRangeMinInput.show()
            # self.rhoRangeMaxInput.show()
            # self.recipErrorInputLine.setText('')
            # self.rhoRangeMinInput.setText('')
            # self.rhoRangeMaxInput.setText('')
            # self.recipErrorInputLine.hide()
            self.recipErrorUnpairedBtn.hide()
            self.recipErrorPltBtn.setToolTip('Removes measuremtns that are out side of defined "Apparent Resistivity" range and/or are manually selected!')
            self.recipErrorBottomTabs.setTabEnabled(1, False)

    def plotManualFiltering(self, index=0):
        attrText = self.filterAttrCombo.currentText()
        dico = {'App. Resistivity':'app',
                'Transfer Resistance':'resist',
                'Reciprocal Error':'reciprocalErrRel',
                'Stacking Error (Dev.)': 'dev'}
        attr = dico[attrText]
        pdebug('plotManualFiltering() with attr={:s} and index={:d}'.format(attr, index))
        if attr not in self.project.surveys[index].df.columns:
            self.errorDump('Attribute {:s} not found for survey.'.format(attrText))
            return
        try:
            self.mwManualFiltering.setCallback(self.project.filterManual)
            self.mwManualFiltering.replot(index=index, attr=attr)
            self.writeLog('# k.filterManual(index={:d}, attr="{:s}" # interactive plot)'.format(index, attr))
        except ValueError as e:
            self.errorDump(e)
            self.mwManualFiltering.clear()

    def errHist(self, index=0):
        pdebug('errHist()')
        if all(self.project.surveys[index].df['irecip'].values == 0) is False:
            if self.iBatch or self.iTimeLapse:
                self.mwRecipError.setCallback(self.project.showErrorDist)
                self.mwRecipError.replot(index=index)
                self.writeLog('k.showErrorDist(index={:d})'.format(index))
            else: 
                self.mwRecipError.plot(self.project.showErrorDist)
                self.writeLog('k.showErrorDist()')
        else:
            pass
    
    def plotError(self, index=0):
        if len(self.project.surveys) == 0:
            return
        self.mwFitError.setCallback(self.project.showError)
        self.mwFitError.replot(index=index)
        self.project.err = False
        self.writeLog('k.showError(index={:d})'.format(index))
        
    def topoInterpBtnFunc(self):
        dfelec = self.elecTable.getTable()
        elec = dfelec[['x','y','z']].values
        if self.topoTable.useNarray:
            topo = self.topoTable.xyz
        else:
            topo = self.topoTable.getTable()[['x','y','z']].values
        inan = ~np.isnan(elec[:,2])
        inan2 = ~np.isnan(topo[:,2])
        points = np.r_[elec[inan,:2], topo[inan2,:2]]
        values = np.r_[elec[inan,2], topo[inan2,2]]
        xi = elec[~inan,:2]
        if self.project.typ[-1] == '2':
            if points.shape[0] < 2:
                self.errorDump('No enough known points for interpolation. Need at least two.')
                return
            from scipy.interpolate import interp1d
            func = interp1d(points[:,0], values, kind='linear', fill_value='extrapolate')
            zi = func(xi[:,0])
        else:
            from scipy.interpolate import griddata
            if points.shape[0] < 3:
                self.errorDump('No enough known points for interpolation. Need at least three.')
                return
            zi = griddata(points, values, xi, method='linear', fill_value=0)
        elec[~inan,2] = zi
        dfelec.loc[:, ['x','y','z']] = elec
        self.elecTable.setTable(dfelec)
        self.infoDump('Interpolation successful.')
        
    
    def updateElec(self):
        pdebug('updateElec()')
        try:
            elec = self.elecTable.getTable()
            ok = False
            if self.tempElec is None:
                ok = True
            else:
                if self.project.iForward is True:
                    ok = True
                elif np.sum(elec[['x','y','z']].values-self.tempElec[['x','y','z']].values) != 0:
                    ok = True
                else:
                    ok = False
            if ok:
                self.tempElec = elec
                elecList = None
                if self.pseudo3DCheck.isChecked():
                    elecList = self.project.split3DGrid(elec.copy(), changeLabel=False) # splitting lines
                    self.project.setPseudo3DElec(pd.concat(elecList, axis=0, ignore_index=True))
                    elecList = self.project._create2DLines() # convert to horizontal 2D lines
                
                self.project.setElec(elec, elecList)
                if self.project.elec['remote'].sum() > 0:
                    self.meshQuadGroup.setEnabled(False)
                else:
                    self.meshQuadGroup.setEnabled(True)
                self.writeLog('#k.setElec(elec)')
                # TODO don't know how to write to log this
                self.project.mesh = None
                self.mwMesh.clear()
                if len(self.project.surveys) > 0:
                    self.plotPseudo()
                    self.plotManualFiltering()
                    if self.project.typ[0] == 'c':
                        self.plotPseudoIP()
        except Exception as e:
            self.errorDump('Error updating pseudosection: ' + str(e))
    
    def plotPseudo(self):
        pdebug('plotPseudo()')
        if (pvfound 
            & ((self.project.typ == 'R3t') | (self.project.typ == 'cR3t')) 
            & (self.boreholeCheck.isChecked() is False)):
            self.pseudoPlotter.clear()
            self.project.showPseudo(ax=self.pseudoPlotter, **self.pParams, 
                                        darkMode=eval(resipySettings.param['dark']))
        elif self.pseudo3DCheck.isChecked():
            self.pseudoPlotter.clear()
            self.project.pseudo3DSurvey.showPseudo(ax=self.pseudoPlotter, threed=True, vmin=self.pParams['vmin'],
                                                   contour=self.psContourCheck.isChecked(),
                                                   vmax=self.pParams['vmax'], darkMode=eval(resipySettings.param['dark']))
        else:
            self.mwPseudo.setCallback(self.project.showPseudo)
            self.mwPseudo.replot(aspect='auto', **self.pParams)
        self.writeLog('k.showPseudo()')
        QApplication.processEvents()
    
    def plotPseudoIP(self):
        pdebug('plotPseudoIP()')
        if (pvfound 
            & ((self.project.typ == 'R3t') | (self.project.typ == 'cR3t')) 
            & (self.boreholeCheck.isChecked() is False)):
            self.pseudoPlotterIP.clear()
            self.project.showPseudoIP(ax=self.pseudoPlotterIP, **self.pParamsIP, 
                                          darkMode=eval(resipySettings.param['dark']))
        elif self.pseudo3DCheck.isChecked():
            self.pseudoPlotterIP.clear()
            self.project.pseudo3DSurvey.showPseudoIP(ax=self.pseudoPlotterIP, threed=True, vmin=self.pParamsIP['vmin'],
                                                     contour=self.psContourCheck.isChecked(),
                                                     vmax=self.pParamsIP['vmax'], darkMode=eval(resipySettings.param['dark']))
        else:
            self.mwPseudoIP.setCallback(self.project.showPseudoIP)
            self.mwPseudoIP.replot(aspect='auto', **self.pParamsIP)
        self.writeLog('k.showPseudoIP()')
        QApplication.processEvents()

    def activateTabs(self, val=True):
        if self.iForward is False:
            self.tabs.setTabEnabled(1,val)
            self.tabs.setTabEnabled(2,val)
            self.tabs.setTabEnabled(4,val)
            self.tabs.setTabEnabled(5,val)
            # self.tabs.setTabEnabled(6,val)
            # try:
            if self.m3DRadio.isChecked():
                self.mwManualFiltering.hide()
                # if len(self.project.surveys) > 0 and all(self.project.surveys[0].df['irecip'].values == 0):
                #     self.tabs.setTabEnabled(1, False)
                #     self.tabPreProcessing.setTabEnabled(0, False)
                # else:
                #     self.tabs.setTabEnabled(1, True)
                # if self.ipCheck.checkState() == Qt.Checked:
                #     self.tabs.setTabEnabled(1, True)
            else:
                self.mwManualFiltering.show()
            # except:
            #     pass
            
        else:
            self.tabs.setTabEnabled(2,val)
        self.meshTrianGroup.setFocus() # needed after tab activation as default is the self.fmdBox (QLineEdit) for some strange reasons...

    def errorCombosShow(self, state=False): #showing/hiding pre-processing comboboxes
        self.recipErrorfnamesCombo.setCurrentIndex(0)
        self.errFitfnamesCombo.setCurrentIndex(0)
        self.iperrFitfnamesCombo.setCurrentIndex(0)
        self.phasefiltfnamesCombo.setCurrentIndex(0)
        if state == False: 
            self.recipErrorfnamesComboLabel.hide()
            self.recipErrorfnamesCombo.hide()
            self.errFitfnamesComboLabel.hide()
            self.errFitfnamesCombo.hide()
            self.iperrFitfnamesCombo.hide()
            self.iperrFitfnamesComboLabel.hide()
            self.phasefiltfnamesComboLabel.hide()
            self.phasefiltfnamesCombo.hide()
        else:
            self.recipErrorfnamesComboLabel.show()
            self.recipErrorfnamesCombo.show()
            self.errFitfnamesComboLabel.show()
            self.errFitfnamesCombo.show()
            self.iperrFitfnamesCombo.show()
            self.iperrFitfnamesComboLabel.show()
            self.phasefiltfnamesComboLabel.show()
            self.phasefiltfnamesCombo.show()
            
            
    def importFile(self, fname):
        pdebug('importFile:', fname)
        if len(self.project.surveys) > 0:
            self.project.surveys = []
        try:                
            self.loadingWidget('Loading data, please wait...', False) # for large datasets
            self.ipCheck.setEnabled(True)
            self.psContourCheck.setEnabled(True)
            self.mergElecCheck.setEnabled(True)
            self.fname = fname
            if float(self.spacingEdit.text()) == -1:
                spacing = None
            else:
                spacing = float(self.spacingEdit.text())

            self.project.createSurvey(self.fname, ftype=self.ftype, spacing=spacing,
                                 parser=self.parser)
            self.writeLog('k.createSurvey("{:s}", ftype="{:s}")'.format(self.fname, self.ftype))
            self.settingUI()
            self.loadingWidget(exitflag=True)
            self.importDataBtn.setText(os.path.basename(fname) + ' (Press to change)')
            self.invNowBtn.setEnabled(True)
            self.activateTabs(True)
            self.infoDump(fname + ' imported successfully')
        except Exception as e:
            self.loadingWidget(exitflag=True)
            pdebug('importFile: ERROR:', e)
            self.errorDump('Importation failed. File is not being recognized. \
                      Make sure you have selected the right file type.')
            pass

            
    def settingUI(self):
        pdebug('importFile: setting up UI')
        if 'magErr' in self.project.surveys[0].df.columns:
            self.a_wgt.setText('0.0')
            self.b_wgt.setText('0.0')
        if np.sum(self.project.surveys[0].df['irecip'].values != 0) < 2:
            # we need more than a single reciprocal to fit error model and so
            listOfNoRecipBtnShow = [self.pseudo3DCheck.isChecked(), self.regular3DCheck.isChecked(), # we can't import group reciprocal files
                                    self.timeLapseCheck.isChecked(), self.batchCheck.isChecked()]
            if not np.any(listOfNoRecipBtnShow):
                self.importDataRecipBtn.show()
            self.recipOrNoRecipShow(recipPresence=False)
        else:
            self.recipOrNoRecipShow(recipPresence=True)
            self.importDataRecipBtn.hide()
            self.tabPreProcessing.setTabEnabled(2, True)
            self.filterAttrCombo.addItem('Reciprocal Error')
            self.plotError()
            self.errHist()
        if 'dev' in self.project.surveys[0].df.columns:
            self.filterAttrCombo.addItem('Stacking Error (Dev.)')

        if self.project.elec['remote'].sum() > 0:
            self.meshQuadGroup.setEnabled(False)
        else:
            self.meshQuadGroup.setEnabled(True)
        if self.boreholeCheck.isChecked() is True:
            self.project.setBorehole(True)
        else:
            self.project.setBorehole(False)
        self.plotManualFiltering()
        if self.pseudo3DCheck.isChecked() and self.project.pseudo3DSurvey is not None:
            self.elecTable.initTable(self.project.pseudo3DSurvey.elec)
        else:
            self.elecTable.initTable(self.project.elec)
        self.tabImporting.setTabEnabled(1,True)
        if 'ip' in self.project.surveys[0].df.columns:
            if np.sum(self.project.surveys[0].df['ip'].values) > 0 or np.sum(self.project.surveys[0].df['ip'].values) < 0: # np.sum(self.project.surveys[0].df['ip'].values) !=0 will result in error if all the IP values are set to NaN
                self.ipCheck.setChecked(True)
            if self.ftype == 'Syscal':
                self.dcaButton.setEnabled(True)
                self.dcaProgress.setEnabled(True)
        if np.isnan(self.project.elec[['x','y','z']].values).any(): # for users who import messed up topography files (res2dinv mostly)
            self.topoInterpBtnFunc()
            self.updateElec()
        self.plotPseudo()
        self.nbElecEdit.setText(str(self.project.elec.shape[0]))
        self.elecDxEdit.setText('{:.2f}'.format(np.abs(np.diff(self.project.elec[~self.project.elec['remote']]['x'].values[:2]))[0]))
        self.fnamesCombo.hide()
        self.fnamesComboLabel.hide()
        

    def restartFunc(self):
        pdebug('restartFunc: creating new R2 object')
        if self.loadedProjectFlag is False:
            self.project = Project(self.newwd, typ=self.typ) # create new R2 instance
            self.project.darkMode = eval(resipySettings.param['dark']) # set matplotlib theme
        self.writeLog('from resipy import Project')
        self.writeLog('k = Project("{:s}", typ="{:s}")'.format(self.newwd, self.typ))
        '''actually we don't really need to instanciate a new object each
        time but it's safe otherwise we would have to reset all attributes
        , delete all surveys and clear parameters
        '''
        # reinitiate flags
        self.project.iBatch = self.iBatch
        self.project.setBorehole(self.iBorehole)
        self.project.iTimeLapse = self.iTimeLapse
        self.project.iForward = self.iForward
        if self.iTimeLapse is True:
            self.reg_mode.setCurrentIndex(2)
        else:
            self.reg_mode.setCurrentIndex(0)
        self.activateTabs(False)
        
        # importing
        self.parser = None
        self.plotAspect = 'equal'
        self.wdBtn.setText('Working directory:' + os.path.basename(self.project.dirname))
        self.importDataBtn.setText('Import Data')
        self.ipCheck.setChecked(False)
        self.ipCheck.setEnabled(False)
        self.importDataRecipBtn.hide()
        self.fnamesCombo.hide()
        self.fnamesComboLabel.hide()
        self.psContourCheck.setEnabled(False)
        self.mergElecCheck.setEnabled(False)
        self.tabImporting.setTabEnabled(1, False)
        self.mwPseudo.clear() # clearing figure
        if pvfound:
            self.pseudoPlotter.clear()
        self.tempElec = None
        self.elecTable.initTable(np.array([['',''],['','']]))
        self.topoTable.initTable(np.array([['',''],['','']]))

        # importing - IP stuff
        self.phiConvFactor.setEnabled(True)
        self.phiConvFactorlabel.setEnabled(True)
        if self.ftype != 'Syscal':
            self.phiConvFactor.setText('1')
        else:
            self.phiConvFactor.setText('1.2')
        if self.ftype == 'ProtocolIP':
            self.phiConvFactor.setText('')
            self.phiConvFactor.setEnabled(False)
            self.phiConvFactorlabel.setEnabled(False)

        # pre-processing
        self.errorCombosShow(False)
        for combobox in self.prepFnamesComboboxes:
            combobox.blockSignals(True)
            combobox.clear()
        self.recipErrorInputLabel.hide()
        self.rhoRangeInputLabel.show()
        self.rhoRangeMinInput.show()
        self.rhoRangeMaxInput.show()
        self.recipErrorInputLine.setText('')
        self.rhoRangeMinInput.setText('')
        self.rhoRangeMaxInput.setText('')
        self.recipErrorInputLine.hide()
        self.filterAttrCombo.clear()
        self.filterAttrCombo.addItems(self.filterAttrComboItems)
        self.mwManualFiltering.clear()
        self.recipErrApplyToAll = True
        self.recipErrDataIndex = -1 # -1 is apply the function to all individually
        self.errFitApplyToAll = True
        self.errFitDataIndex = -1
        self.iperrFitDataIndex = -1
        self.errFitType.setCurrentIndex(0)
        self.mwFitError.clear()
        self.mwIPFitError.clear()
        self.iperrFitType.setCurrentIndex(0)
        self.errFitPlotIndexList = []
        self.iperrFitPlotIndexList = []
        self.phivminEdit.setText('0')
        self.phivmaxEdit.setText('25')
        self.dcaProgress.setValue(0)
        self.dcaButton.setEnabled(False)
        self.dcaProgress.setEnabled(False)
        self.phaseFiltDataIndex = -1
        self.tabPreProcessing.setTabEnabled(0, True)
        self.tabPreProcessing.setTabEnabled(1, False)
        self.tabPreProcessing.setTabEnabled(2, False)
        self.tabPreProcessing.setTabEnabled(3, False)

        # mesh
        self.mwMesh.clear()
        self.meshLogText.clear()
        self.regionTable.reset()
        if pvfound:
            self.mesh3Dplotter.clear()
        
        #forward model
        self.mwFwdPseudo.clear()
        self.mwFwdPseudoIP.clear()
        if self.iForward is True: # in case the "Restart" button is used while Forward radio is already checked
            self.tabImporting.setTabEnabled(1, True) # here because restartFunc() set it to False
            self.tabs.setTabEnabled(2, True)
            self.tabs.setTabEnabled(3, True)
            self.tabs.setTabEnabled(4, False)
            self.tabs.setTabEnabled(5, False)
            self.ftypeCombo.setEnabled(False)
            self.invRadio.setChecked(False)
            self.importDataBtn.setEnabled(False)
            self.timeLapseCheck.setEnabled(False)
            self.batchCheck.setEnabled(False)
            self.tabImporting.setTabEnabled(2,False) # no custom parser needed
            self.nbElecEdit.setEnabled(True)
            self.ipCheck.setEnabled(True)
            self.psContourCheck.setEnabled(False)
            self.mergElecCheck.setEnabled(False)
            if pvfound: # 3D pseudo-sections?
                self.pseudo3Dplotter.clear()
                self.pseudo3DplotterIP.clear()
            
        # inversion settings
        self.flux_type.setCurrentIndex(0)
        self.singular_type.setChecked(False)
        self.res_matrix.setCurrentIndex(1)
        self.scale.setText('1.0')
        self.patch_x.setText('1')
        self.patch_z.setText('1')
        self.inv_type.setCurrentIndex(1)
        self.data_type.setCurrentIndex(1)
        self.max_iterations.setText('10')
        self.error_mod.setCurrentIndex(1)
        self.alpha_aniso.setText('1.0')
        self.min_error.setText('0.01')
        self.a_wgt.setText('0.01')
        self.b_wgt.setText('0.02')
        self.rho_min.setText('-10e10')
        self.rho_max.setText('10e10')
        self.target_decrease.setText('0')
        self.helpSection2.setText('Click on the labels and help will be displayed here')

        # inversion
        self.invtabs.setCurrentIndex(0)
        self.invtabs.setTabEnabled(1, False)
        self.logText.setText('')
        self.mwRMS.clear()
        self.mwIter.clear()
        self.surveyCombo.clear()
        self.attrCombo.clear()
        self.vminEdit.setText('')
        self.vmaxEdit.setText('')
        self.pvthreshMax.setText('')
        self.pvthreshMin.setText('')
        self.pvxslices.setText('')
        self.pvyslices.setText('')
        self.pvzslices.setText('')
        self.pvcontour.setText('')
        self.pvgridCheck.setChecked(False)
        self.pvdelaunay3dCheck.setChecked(False)
        self.doiCheck.setChecked(False)
        self.doiSensCheck.setChecked(False)
        self.displayParams = {'index':0,'edge_color':'none',
                            'sens':True, 'attr':'Resistivity(Ohm-m)',
                            'contour':False, 'vmin':None, 'vmax':None,
                            'cmap':'viridis', 'sensPrc':0.5, 'clipCorners':False,
                            'doi':self.modelDOICheck.isChecked(),
                            'doiSens':False, 'pvdelaunay3d': False,
                            'pvslices':([],[],[]), 'pvthreshold':None,
                            'pvgrid':False, 'pvcontour':[], 'aspect':'equal'}
        self.mwInv.clear()
        self.tabs.setTabEnabled(6, False) # there is no post-processing after a reset!
        self.mwInvError.clear()
        self.mwInvErrorIP.clear()
        self.mwInvError2.clear()
        self.rangeInvErrorMinInput.setText('')
        self.rangeInvErrorMaxInput.setText('')

        
    def disableOptionsPseudo3D(self, flag):
        if flag:
            self.doiSensCheck.setChecked(False)
            self.doiSensCheck.setEnabled(False)
            self.aspectCheck.setEnabled(False)
            self.sensWidget.setEnabled(False)
        else:
            self.doiSensCheck.setEnabled(True)
            self.aspectCheck.setEnabled(True)
            self.sensWidget.setEnabled(True)
    
    
    def replotSection(self): # main plotting function
        pdebug('replotSection() with', self.displayParams)
        index = self.displayParams['index']
        edge_color = self.displayParams['edge_color']
        sens = self.displayParams['sens']
        attr = self.displayParams['attr']
        contour = self.displayParams['contour']
        clipCorners = self.displayParams['clipCorners']
        vmin = self.displayParams['vmin']
        vmax = self.displayParams['vmax']
        cmap = self.displayParams['cmap']
        sensPrc = self.displayParams['sensPrc']
        doi = self.displayParams['doi']
        doiSens = self.displayParams['doiSens']
        pvslices = self.displayParams['pvslices']
        pvthreshold = self.displayParams['pvthreshold']
        pvgrid = self.displayParams['pvgrid']
        pvdelaunay3d = self.displayParams['pvdelaunay3d']
        pvcontour = self.displayParams['pvcontour']
        aspect = self.displayParams['aspect']
        if self.project.typ[-1] == '2':
            if self.pseudo3DCheck.isChecked() and self.surveyCombo.currentIndex() == 0:
                self.disableOptionsPseudo3D(True)
                self.vtkWidget.clear()
                self.vtkWidget.clear_plane_widgets()
                self.project.showResults(index=-1, ax=self.vtkWidget, attr=attr, edge_color=edge_color, clipCorners=clipCorners,
                                         vmin=vmin, vmax=vmax, color_map=cmap, pvgrid=True, cropMesh=self.iCropping, pseudo3DContour=contour,
                                         background_color=(0.8,0.8,0.8), cropMaxDepth=self.cropBelowFmd.isChecked())
            else:
                self.disableOptionsPseudo3D(False)
                self.mwInv.replot(threed=False, aspect=aspect,
                                  index=index, edge_color=edge_color,
                                  contour=contour, sens=sens, clipCorners=clipCorners,
                                  attr=attr, vmin=vmin, vmax=vmax, color_map=cmap, 
                                  sensPrc=sensPrc, doi=doi, doiSens=doiSens)
                self.writeLog('k.showResults(index={:d}, edge_color="{:s}",'
                              ' contour={:s}, sens={:s}, attr="{:s}", vmin={:s}, '
                              'vmax={:s}, color_map="{:s}", sensPrc={:.2f}, doi={:s},'
                              ' doiSens={:s})'.format(index, edge_color, str(contour),
                                str(sens), attr, str(vmin), str(vmax), cmap, sensPrc, str(doi),
                                                     str(doiSens)))
        else:
            # mwInvResult3D.replot(threed=True, index=index, attr=attr,
                                  # vmin=vmin, vmax=vmax, color_map=cmap)

            # if self.project.iTimeLapse and index == 0:
                # fname = os.path.join(self.project.dirname, 'f001.vtk')
            # else:
                # fname = os.path.join(self.project.dirname, 'f{:03d}.vtk'.format(index+1))
            # m = pv.read(fname)

            self.vtkWidget.clear()
            self.vtkWidget.clear_plane_widgets()
            # self.project.meshResults[index].show3D(ax=self.vtkWidget, attr=attr,
            #                                   edge_color=edge_color, vmin=vmin,
            #                                   vmax=vmax, color_map=cmap,
            #                                   background_color=(0.8,0.8,0.8),
            #                                   pvslices=pvslices,
            #                                   pvthreshold=pvthreshold,
            #                                   pvgrid=pvgrid,
            #                                   pvcontour=pvcontour)

            self.project.showResults(index=index, ax=self.vtkWidget, attr=attr,
                                              edge_color=edge_color, vmin=vmin,
                                              vmax=vmax, color_map=cmap,
                                              background_color=(0.8,0.8,0.8),
                                              pvslices=pvslices,
                                              pvthreshold=pvthreshold,
                                              pvgrid=pvgrid, pvdelaunay3d=pvdelaunay3d,
                                              pvcontour=pvcontour)
            self.writeLog('k.showResults(index={:d}, attr="{:s}", edge_color="{:s}", vmin={:s}, '
                          'vmax={:s}, color_map="{:s}", background_color=(0.8, 0.8, 0.8),'
                          'pvslices={:s}, pvthreshold={:s}, pvdelaunay3d={:s}, pvgrid={:s}, pvcontour={:s})'.format(
                          index, attr, edge_color, str(vmin), str(vmax), cmap, str(pvslices),
                          str(pvthreshold), str(pvdelaunay3d), str(pvgrid), str(pvcontour)))
  

    def displayInvertedResults(self): # after inversion, plot graph
        self.invtabs.setCurrentIndex(1) # show tab
        self.invtabs.setTabEnabled(1, True)
        self.invtabs.setTabEnabled(2, True)
        if self.project.typ[-1] == '2':
            if self.pseudo3DCheck.isChecked():
                self.showStackedLayout.setCurrentIndex(1)
            else:
                self.showStackedLayout.setCurrentIndex(0)
            self.mwInv.setCallback(partial(self.project.showResults, cropMaxDepth=self.cropBelowFmd.isChecked()))
        else:
            self.showStackedLayout.setCurrentIndex(1)
        if self.project.typ == 'R2' or self.project.typ == 'R3t':
            self.displayParams['attr'] = 'Resistivity(log10)'
        else:
            self.displayParams['attr'] = 'Sigma_real(log10)'
        self.surveyCombo.clear()
        if self.pseudo3DCheck.isChecked():
            self.surveyCombo.addItem('Pseudo 3D display')
        for m in self.project.meshResults:
            self.surveyCombo.addItem(m.mesh_title)
        index = 0
        if self.project.iForward: # display the inverted not the initial
            self.surveyCombo.setCurrentIndex(1)
            index = 1
            self.displayParams['index'] = 1
        self.attrCombo.clear()
        attrs = self.project.meshResults[index].df.keys()
        c = 0
        ci = 0
        for i, attr in enumerate(attrs):
            if attr not in ['param', 'region', 'zone', 'elm_id', 'cellType', 'X', 'Y', 'Z']:
                self.attrCombo.addItem(attr)
                if attr == self.displayParams['attr']:
                    ci = c
                c += 1
        self.attrCombo.setCurrentIndex(ci)
        self.vminEdit.setText('')
        self.vmaxEdit.setText('')
        self.doiCheck.setChecked(self.modelDOICheck.isChecked())
        self.contourCheck.setChecked(False)
        self.clipCornersCheck.setChecked(False)
        self.edgeCheck.setChecked(False)
        self.replotSection() # this plot the results
       

#%% updater function and wine check
    # based on https://kushaldas.in/posts/pyqt5-thread-example.html
    def updateChecker(self): # check for new updates on gitlab
        version = ResIPy_version
        newChanges = ''
        try:
            versionSource = urlRequest.urlopen('https://raw.githubusercontent.com/hkexgroup/resipy/master/src/version.txt')
            versionCheck = versionSource.read().decode()
            version = versionCheck.split()[1] # assuming version number is in 2nd line of version.txt
            changeLogSource = urlRequest.urlopen('https://github.com/hkexgroup/resipy/raw/master/CHANGELOG')
            changeLogTxt = changeLogSource.read().decode()
            newChangesRaw = changeLogTxt.split('\n\n')[0].split('\n')
            newChanges = ''.join('{}<br>'*len(newChangesRaw[2:])).format(*newChangesRaw[2:])
            print('online version :', version)
        except:
            pass
        return [version, newChanges]
    
    def updateCheckerShow(self, msgInput):
        version = msgInput[0] # online version 
        newChanges = msgInput[1] # documented changes 
        ResIPy_version_int = int(ResIPy_version.replace('.',''))
        version_int = int(version.replace('.',''))
        if ResIPy_version_int < version_int: # check if version is older than online version 
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Information)
            msg.setText('''<b>ResIPy version %s is available</b>''' % (version))
            msg.setInformativeText('''Please download the latest version of ResIPy at:\
                                   <p><a href='https://gitlab.com/hkex/resipy#downloads'>https://gitlab.com/hkex/resipy</a></p>\
                                   New updates:<br>\
                                   %s''' % newChanges)
            msg.setWindowTitle("New version available")
            bttnUpY = msg.addButton(QMessageBox.Yes)
            bttnUpY.setText('Update')
            bttnUpN = msg.addButton(QMessageBox.No)
            bttnUpN.setText('Ignore')
            msg.setDefaultButton(bttnUpY)
            msg.exec_()
            if msg.clickedButton() == bttnUpY:
                webbrowser.open('https://gitlab.com/hkex/resipy#gui-for-r2-family-code') # can add download link, when we have a direct dl link
    
    def checkWine(self): # check if wine is installed, on Unix system
        #check operating system
        wineInstalled = sysinfo['wineCheck']
        return wineInstalled

    def checkWineShow(self, wineInstalled):
        print('====================system-information======================')

        for key in sysinfo.keys():
            print(key + ' = ' + str(sysinfo[key]), end='')
            if key == 'max_freq':
                print(' Mhz')
            elif key.find('Memory') != -1:
                print(' Gb')
            elif key == 'wineCheck' and sysinfo['OS'] == 'Windows':
                print(' (Native Windows)')
            else:
                print('')
         
        print('============================================================')
        #print('is wine installed?', wineInstalled)
        if OS == 'Linux':
            OpSys = 'Linux'
        elif OS == 'Darwin':
            OpSys = 'macOS'
        if wineInstalled is False:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText('''<b>No "wine" is installed on your %s</b>''' % (OpSys))
            msg.setInformativeText('''ResIPy needs "wine" to run properly,<br>without "wine", no inversion or triangular meshing is possible.<br>''')
            msg.setWindowTitle('"Wine" is not detected!')
            bttnUpY = msg.addButton(QMessageBox.Yes)
            bttnUpY.setText('Learn more')
            bttnUpN = msg.addButton(QMessageBox.No)
            bttnUpN.setText('Continue')
            msg.setDefaultButton(bttnUpY)
            msg.exec_()
            if msg.clickedButton() == bttnUpY:
                webbrowser.open('https://gitlab.com/hkex/resipy#linux-and-mac-user')
    
    def restartGUI(self):
        try:
            if 'frozen' in resipySettings.param.keys() and resipySettings.param['frozen'] == 'True': # Windows/Linux frozen package only
                exe_path = resipySettings.param['exe_path']
                Popen([exe_path], shell=False, stdout=None, stdin=None)
                sys.exit()
            else:
                os.execl(sys.executable, sys.executable, *sys.argv)
        except SystemExit:
            sys.exit() # as sys.exit is an exception already, we need this to close the app
        except: # in case something goes wrong!
            self.errorDump('Could not restart ResIPy, try manually closing and reopening!!')
    
    def darkModeFunc(self):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Question)
        msg.setText('''<b>Restart ResIPy to change the theme?</b>''')
        msg.setWindowTitle('Restart?')
        bttnY = msg.addButton(QMessageBox.Yes)
        bttnY.setText('Restart')
        bttnN = msg.addButton(QMessageBox.No)
        bttnN.setText('Cancel')
        msg.setDefaultButton(bttnN)
        msg.exec_()
        if msg.clickedButton() == bttnY:
            if resipySettings.param['dark'] == 'True':
                resipySettings.param['dark'] = 'False'
            else:
                resipySettings.param['dark'] = 'True'
                
            resipySettings.genLocalSetting()
            self.restartGUI()


if __name__ == '__main__':
    catchErrors()
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    
    import matplotlib # for dark mode should be imported here
    matplotlib.use('Qt5Agg')
    
    # loading settings
    from resipy.Settings import Settings
    resipySettings = Settings()
    localSettings = resipySettings.retLocalSetting()
    if localSettings != None: # dark mode is selected?
        if 'dark' in resipySettings.param.keys():
            if resipySettings.param['dark'] == 'True':
                # dark theme GUI
                dark_palette = QPalette()
                dark_palette.setColor(QPalette.Window, QColor(53, 53, 53))
                dark_palette.setColor(QPalette.WindowText, Qt.white)
                dark_palette.setColor(QPalette.Base, QColor(35, 35, 35))
                dark_palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
                dark_palette.setColor(QPalette.ToolTipBase, QColor(25, 25, 25))
                dark_palette.setColor(QPalette.ToolTipText, Qt.white)
                dark_palette.setColor(QPalette.Text, Qt.white)
                dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
                dark_palette.setColor(QPalette.ButtonText, Qt.white)
                dark_palette.setColor(QPalette.BrightText, Qt.red)
                dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
                dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
                dark_palette.setColor(QPalette.HighlightedText, QColor(35, 35, 35))
                dark_palette.setColor(QPalette.Active, QPalette.Button, QColor(53, 53, 53))
                dark_palette.setColor(QPalette.Disabled, QPalette.ButtonText, Qt.darkGray)
                dark_palette.setColor(QPalette.Disabled, QPalette.WindowText, Qt.darkGray)
                dark_palette.setColor(QPalette.Disabled, QPalette.Text, Qt.darkGray)
                dark_palette.setColor(QPalette.Disabled, QPalette.Light, QColor(53, 53, 53))
                app.setPalette(dark_palette)
                # dark theme matplotlib plots
                matplotlib.style.use('dark_background')
        else:
            resipySettings.param['dark'] = 'False'
            resipySettings.genLocalSetting()
            
    else:
        print('Generating local settings file...')
        resipySettings.param['dark'] = 'False'
        resipySettings.genLocalSetting()

    app.setWindowIcon(QIcon(os.path.join(bundle_dir, 'logo.png'))) # that's the true app icon
    
    splash_pix = QPixmap(os.path.join(bundle_dir, 'loadingLogo.png'))
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setWindowFlags(Qt.WindowStaysOnTopHint | Qt.FramelessWindowHint)
    splash.setEnabled(False)

    # adding progress bar
    progressBar = QProgressBar(splash)
    progressBar.setMaximum(10)
    progressBar.setGeometry(100, splash_pix.height() - 50, splash_pix.width() - 200, 20)

    from resipy import ResIPy_version, sysinfo
    global sysinfo
    splash.show()
    splash.showMessage("Loading libraries", Qt.AlignBottom | Qt.AlignCenter, Qt.black)
    app.processEvents()

    # in this section all import are made except the one for pyQt
    progressBar.setValue(1)
    app.processEvents()
    
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
    from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.figure import Figure
#    from mpl_toolkits.mplot3d import axes3d # supersede by pyvista
    progressBar.setValue(2)
    app.processEvents()

    import numpy as np
    progressBar.setValue(4)
    app.processEvents()
    import pandas as pd
    progressBar.setValue(6)
    app.processEvents()
    from datetime import datetime
    progressBar.setValue(8)
    app.processEvents()
    from matplotlib import rcParams
    rcParams.update({'font.size': 11}) # CHANGE HERE for graph font size

    # library needed for update checker + wine checker
    import platform
    OS = platform.system()
    from subprocess import PIPE, Popen
    from urllib import request as urlRequest
    import webbrowser

    try:
        import pyvista as pv
        try:
            from pyvistaqt import QtInteractor # newer version
        except:
            from pyvista import QtInteractor # older version
        pvfound = True
    except:
        pvfound = False
        print('WARNING: pyvista not found, 3D plotting capabilities will be limited.')
    
    from resipy import Project
    from resipy.r2help import r2help
    splash.showMessage("ResIPy is ready!", Qt.AlignBottom | Qt.AlignCenter, Qt.black)
    progressBar.setValue(10)
    app.processEvents()

    ex = App()
    ex.show()
    splash.hide() # hiding the splash screen when finished
    sys.exit(app.exec_())
