#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import sys
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
    QStackedLayout, QRadioButton, QGroupBox, QTextBrowser)#, QAction, QButtonGroup, QListWidget, QShortcut)
from PyQt5.QtGui import QIcon, QPixmap, QIntValidator, QDoubleValidator, QColor#, QKeySequence
from PyQt5.QtCore import QThread, pyqtSignal, QTimer#, QProcess, QSize
from PyQt5.QtCore import Qt
from functools import partial
QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True) # for high dpi display
QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

import threading
import traceback

# COMMENT below is using UI automatic testing
# from resipy.R2 import ResIPy_version
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
# try:
#     import pyvista as pv
#     pvfound = True
# except:
#     pvfound = False
#     print('WARNING: pyvista not found, 3D plotting capabilities will be limited.')
# from resipy.R2 import R2
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
    finalError =''
    for errs in errorMsg:
        finalError += errs
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText("<b>Critical error:</b>")
    msg.setInformativeText('''Please see the detailed error below.<br>You can report the errors at:<p><a href='https://gitlab.com/hkex/pyr2/issues'>https://gitlab.com/hkex/pyr2/issues</a></p><br>''')
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
        self.newwd = os.path.join(bundle_dir, 'resipy', 'invdir')

        # UI attributes (sync with self.r2 attributes when instantiated)
        self.r2 = None
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
        self.errorLabel = QLabel('<i style="color:black">Error messages will be displayed here</i>')
        QApplication.processEvents()



        #%% tab 1 importing data
        self.tabImporting = QTabWidget()
        self.tabs.addTab(self.tabImporting, 'Importing')

        self.tabImportingData = QWidget()
        self.tabImporting.addTab(self.tabImportingData, 'Data')
        tabImportingDataLayout = QVBoxLayout()

        
        self.restartBtn = QPushButton('Restart')
        self.restartBtn.setAutoDefault(True)
        self.restartBtn.clicked.connect(self.restartFunc)
        self.restartBtn.setToolTip('Press to reset all tabs and start a new survey.')

        def dimSurvey():
            if self.m2DRadio.isChecked():
                self.typ = self.typ.replace('3t','2')
                if self.r2 is not None:
                    self.r2.typ = self.r2.typ.replace('3t','2')
                # importing tab
                self.elecTable.setColumnHidden(1, True)
                self.topoTable.setColumnHidden(1, True)
                if self.r2 is not None:
                    if self.r2.elec is not None:
                        self.elecTable.initTable(self.r2.elec)
                self.elecDyEdit.setEnabled(False)
                self.fwdRadio.setEnabled(True)
                self.fwdRadio.setChecked(False)
                self.boreholeCheck.setChecked(False)
                self.boreholeCheck.setEnabled(True)
                self.regular3DCheck.setVisible(False)
                # self.ipCheck.setEnabled(True)
                
                #Pre-processing tab
                self.recipErrorBottomTabs.setTabEnabled(0, True)
                self.recipErrorBottomTabs.setCurrentIndex(0)
                self.recipErrorSaveBtn.setVisible(True)
                self.tabPreProcessing.setCurrentIndex(0)
                self.tabPreProcessing.setTabEnabled(0, True)
                try:
                    if not self.r2.surveys[0].df.empty:
                        self.tabs.setTabEnabled(1, True)
                except:
                    pass
                
                # mesh tab
                self.meshQuadGroup.setVisible(True)
                self.meshTrianGroup.setVisible(True)
                meshTetraGroup.setVisible(False)
                meshCustomGroup.setVisible(True)
                instructionLabel.setVisible(True)
                self.resetMeshBtn.setVisible(True)

                # inversion settings
                show3DOptions(False)
                if self.r2 is not None:
                    showIpOptions(self.typ[0] == 'c')

                # inversion tab
                self.contourCheck.setVisible(True)
                self.doiSensCheck.setVisible(True)
                self.sensWidget.setVisible(True)
                show3DInvOptions(False)
            else:
                self.typ = self.typ.replace('2','3t')
                if self.r2 is not None:
                    self.r2.typ = self.r2.typ.replace('2', '3t')

                # importing tab
                self.fwdRadio.setChecked(False)
                self.fwdRadio.setEnabled(False)
                self.elecTable.setColumnHidden(1, False)
                self.topoTable.setColumnHidden(1, False)
                if self.r2 is not None:
                    if self.r2.elec is not None:
                        self.elecTable.initTable(self.r2.elec)
                self.elecDyEdit.setEnabled(True)
                self.invRadio.setChecked(True)
                self.boreholeCheck.setChecked(True) # to disable pseudo-section
                self.boreholeCheck.setEnabled(False)
                self.regular3DCheck.setVisible(True)
                # self.ipCheck.setEnabled(False) # TODO disabling cR3t for now
                
                #Pre-processing tab
                self.recipErrorBottomTabs.setTabEnabled(0, False)
                self.recipErrorSaveBtn.setVisible(False)
                
                try:
                    if all(self.r2.surveys[0].df['irecip'].values == 0):
                        self.tabs.setTabEnabled(1, False)
                        self.tabPreProcessing.setTabEnabled(0, False)
                    else:
                        self.tabs.setTabEnabled(1, True)
                    if self.ipCheck.checkState() == Qt.Checked:
                        self.tabs.setTabEnabled(1, True)
                except:
                    pass

                # mesh tab
                self.meshQuadGroup.setVisible(False)
                self.meshTrianGroup.setVisible(False)
                meshTetraGroup.setVisible(True)
                meshCustomGroup.setVisible(False)
                instructionLabel.setVisible(False)
                self.resetMeshBtn.setVisible(False)

                # inversion settings
                show3DOptions(True)
                if self.r2 is not None:
                    showIpOptions(self.typ[0] == 'c')

                # inversion tab
                self.contourCheck.setVisible(False)
                self.doiSensCheck.setVisible(False)
                self.sensWidget.setVisible(False)
                show3DInvOptions(True)

        self.m2DRadio = QRadioButton('2D')
        self.m2DRadio.setChecked(True)
        self.m2DRadio.toggled.connect(dimSurvey)
        self.m3DRadio = QRadioButton('3D')
        self.m3DRadio.setChecked(False)
        self.m3DRadio.toggled.connect(dimSurvey)
        dimLayout = QHBoxLayout()
        dimLayout.addWidget(self.m2DRadio)
        dimLayout.addWidget(self.m3DRadio)
        dimLayout.setContentsMargins(0,0,0,0)
        dimGroup = QGroupBox()
        dimGroup.setLayout(dimLayout)
        dimGroup.setFlat(True)
        dimGroup.setContentsMargins(0,0,0,0)
        dimGroup.setStyleSheet('QGroupBox{border: 0px;'
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
                if self.r2 is not None: # if there is already an R2 object
                    self.restartFunc()
                    self.reg_mode.setCurrentIndex(2)
                self.importDataBtn.setText('Import multiple datasets')
                self.importDataBtn.clicked.disconnect()
                self.importDataBtn.clicked.connect(getdir)
                self.ipCheck.setEnabled(False) # timelapse IP not available for now
                self.batchCheck.setEnabled(False)
            else:
                self.iTimeLapse = False
                if self.r2 is not None:
                    self.restartFunc()
                    self.reg_mode.setCurrentIndex(0)
                self.importDataBtn.setText('Import Data')
                self.importDataBtn.clicked.disconnect()
                self.importDataBtn.clicked.connect(importDataBtnFunc)
                self.ipCheck.setEnabled(True) # timelapse IP not available for now
                self.batchCheck.setEnabled(True)
        self.timeLapseCheck = QCheckBox('Time-lapse Survey')
        self.timeLapseCheck.stateChanged.connect(timeLapseCheckFunc)
        self.timeLapseCheck.setToolTip('Check to import time-lapse datasets and enable time-lapse inversion.')

        def batchCheckFunc(state):
            if state == Qt.Checked:
                self.iBatch = True
                if self.r2 is not None:
                    self.restartFunc()
                self.importDataBtn.setText('Import multiple datasets')
                self.importDataBtn.clicked.disconnect()
                self.importDataBtn.clicked.connect(getdir)
                self.timeLapseCheck.setEnabled(False)
            else:
                self.iBatch = False
                if self.r2 is not None:
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
                if self.r2 is not None:
                    self.r2.setBorehole(True)
            else:
                self.iBorehole = False
                self.errorGraphs.setTabEnabled(0, True)
                if self.r2 is not None:
                    self.r2.setBorehole(False)
            try:
                if self.fname is not None:
                        self.plotPseudo()
                        self.plotPseudoIP()
            except:
                pass
        self.boreholeCheck = QCheckBox('Unconventional Survey')
        self.boreholeCheck.stateChanged.connect(boreholeCheckFunc)
        self.boreholeCheck.setToolTip('Check if you have an unconventional survey (e.g. boreholes).\nThis will just change the pseudo-section.')

        def regular3DFunc(state):
            if state == Qt.Checked:
                self.lineSpacing.setVisible(True)
                self.lineSpacingLabel.setVisible(True)
                self.create3DBtn.setVisible(True)
                self.importDataBtn.setVisible(False)
            else:
                self.lineSpacing.setVisible(False)
                self.lineSpacingLabel.setVisible(False)
                self.create3DBtn.setVisible(False)
                self.importDataBtn.setVisible(True)
        self.regular3DCheck = QCheckBox('3D survey from regular 2D lines')
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
            self.tabImporting.setTabEnabled(2,False) # no custom parser needed
            self.restartFunc() # let's first from previous inversion
            self.nbElecEdit.setEnabled(True)
            self.tabImporting.setTabEnabled(1, True) # here because restartFunc() set it to False
            self.ipCheck.setEnabled(True)
            self.psContourCheck.setEnabled(False)
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
            self.timeLapseCheck.setChecked(False) # not checked by default
            self.batchCheck.setChecked(False) # not checked by default
            self.activateTabs(False)
        self.invRadio = QRadioButton('Inverse')
        self.invRadio.setChecked(True)
        self.invRadio.toggled.connect(invRadioFunc)
        self.invRadio.setToolTip('To invert data that is already collected.')
        
        dimInvLayout = QHBoxLayout()
        dimInvLayout.addWidget(self.fwdRadio)
        dimInvLayout.addWidget(self.invRadio)
        dimInvLayout.setContentsMargins(0,0,0,0)
        dimInvGroup = QGroupBox()
        dimInvGroup.setLayout(dimInvLayout)
        dimInvGroup.setFlat(True)
        dimInvGroup.setContentsMargins(0,0,0,0)
        dimInvGroup.setStyleSheet('QGroupBox{border: 0px;'
                                'border-style:inset;}')


        # ask for working directory, and survey file to input
        def getwd():
            fdir = QFileDialog.getExistingDirectory(self.tabImportingData, 'Choose Working Directory')
            if fdir != '':
                self.newwd = fdir
                if self.r2 is not None:
                    self.r2.setwd(fdir)
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
                self.ftype = 'Res2Dinv'
                self.fformat = 'DAT (*.dat *.DAT)'
            elif index == 4:
                self.ftype = 'BGS Prime'
                self.fformat = 'DAT (*.dat *.DAT)'
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
        self.ftypeCombo.addItem('Res2Dinv')
        self.ftypeCombo.addItem('BGS Prime')
        self.ftypeCombo.addItem('Sting')
        self.ftypeCombo.addItem('ABEM-Lund')
        self.ftypeCombo.addItem('Lippmann')
        self.ftypeCombo.addItem('ARES (beta)')
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
#            fdir = QFileDialog.getExistingDirectory(self.tabImportingData, 'Choose the directory containing the data', directory=self.datadir)
            
            if fnames != []:
                fdir = os.path.dirname(fnames[0])
                self.restartFunc()
                self.datadir = os.path.dirname(fdir)
                try:
                    if self.r2.iBatch is False:
                        if len(fnames) < 2:
                            self.errorDump('at least two files needed for timelapse.')
                            return
                        self.r2.createTimeLapseSurvey(fnames, ftype=self.ftype, dump=self.infoDump)
                        self.ipCheck.setEnabled(False) # TODO enable IP for timelapse
                        self.infoDump('Time-lapse survey created.')
                    else:
                        self.r2.createBatchSurvey(fnames, ftype=self.ftype, dump=self.infoDump)
                        self.ipCheck.setEnabled(True)
                        self.infoDump('Batch survey created.')
                    self.fnamesCombo.clear()
                    self.psContourCheck.setEnabled(True)
                    
                    for s in self.r2.surveys:
                        self.fnamesCombo.addItem(s.name)
                        self.errFitPlotIndexList.append(0)
                        self.iperrFitPlotIndexList.append(0)
                    self.errorCombosShow(True)
                    errorCombosFill(self.prepFnamesComboboxes)
                    self.fnamesCombo.show()
                    self.fnamesComboLabel.show()
                    self.importDataBtn.setText(os.path.basename(fdir) + ' (Press to change)')
                    self.plotPseudo()
                    self.elecTable.initTable(self.r2.elec)
                    self.tabImporting.setTabEnabled(1,True)
                    self.invNowBtn.setEnabled(True)
                    self.nbElecEdit.setText(str(len(self.r2.elec)))
                    if all(self.r2.surveys[0].df['irecip'].values == 0):
                        self.recipOrNoRecipShow(recipPresence=False)
                    else:
                        self.recipOrNoRecipShow(recipPresence=True)
                        self.tabPreProcessing.setTabEnabled(2, True)
                        self.filterAttrCombo.addItem('Reciprocal Error')
                        self.plotError()
                        self.errHist()
                    self.plotManualFiltering()
                    self.activateTabs(True)
                    if 'dev' in self.r2.surveys[0].df.columns:
                        self.filterAttrCombo.addItem('Stacking Error (Dev.)')
                    if 'ip' in self.r2.surveys[0].df.columns and self.iTimeLapse is False:
                        if np.sum(self.r2.surveys[0].df['ip'].values) > 0 or np.sum(self.r2.surveys[0].df['ip'].values) < 0: # np.sum(self.r2.surveys[0].df['ip'].values) !=0 will result in error if all the IP values are set to NaN
                            self.ipCheck.setChecked(True)
                        if self.ftype == 'Syscal':
                            self.dcaButton.setEnabled(True)
                            self.dcaProgress.setEnabled(True)
                except Exception as e:
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
                self.importDataRecipBtn.setText(os.path.basename(fnameRecip))
                # if float(self.spacingEdit.text()) == -1:
                #     spacing = None
                # else:
                #     spacing = float(self.spacingEdit.text())
                self.r2.addData(fname=fnameRecip, ftype=self.ftype, parser=self.parser)
                if all(self.r2.surveys[0].df['irecip'].values == 0) is False:
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
                self.plotManualFiltering()
                self.infoDump(fnameRecip + ' imported successfully')
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

        def create3DFunc():
            fnames, _ = QFileDialog.getOpenFileNames(self.tabImportingData, 'Select file(s)', self.datadir, self.fformat)            
            if fnames != []:
                fdir = os.path.dirname(fnames[0])
                self.restartFunc()
                self.datadir = os.path.dirname(fdir)
                val = float(self.lineSpacing.text())
                try:
                    self.r2.create3DSurvey(fnames, lineSpacing=val, ftype=self.ftype, parser=self.parser)
                    self.infoDump('3D survey from regular 2D lines created.')
                    self.ipCheck.setEnabled(True)
                    self.psContourCheck.setEnabled(True)
                    self.importDataBtn.setText(os.path.basename(fdir) + ' (Press to change)')
                    self.calcAspectRatio()
                    if 'magErr' in self.r2.surveys[0].df.columns:
                        self.a_wgt.setText('0.0')
                        self.b_wgt.setText('0.0')
                    self.recipOrNoRecipShow(recipPresence = True)
                    self.importDataRecipBtn.hide()
                    self.tabPreProcessing.setTabEnabled(2, True)
                    self.plotError()
                    self.errHist()
                    self.plotManualFiltering()
                    self.elecTable.initTable(self.r2.elec)
                    self.tabImporting.setTabEnabled(1,True)
                    if 'ip' in self.r2.surveys[0].df.columns:
                        if np.sum(self.r2.surveys[0].df['ip'].values) > 0 or np.sum(self.r2.surveys[0].df['ip'].values) < 0: # np.sum(self.r2.surveys[0].df['ip'].values) !=0 will result in error if all the IP values are set to NaN
                            self.ipCheck.setChecked(True)
                        if self.ftype == 'Syscal':
                            self.dcaButton.setEnabled(True)
                            self.dcaProgress.setEnabled(True)               
                    self.plotPseudo()
                    self.invNowBtn.setEnabled(True)
                    self.activateTabs(True)
                    self.nbElecEdit.setText(str(len(self.r2.elec)))
                    self.elecDxEdit.setText('%s' %(self.r2.elec[~self.r2.iremote,:][1,0]-self.r2.elec[~self.r2.iremote,:][0,0]))
                    self.fnamesCombo.hide()
                    self.fnamesComboLabel.hide()
                    if np.isnan(self.r2.elec).any():
                        self.topoInterpBtnFunc()
                except Exception as e:
                    print('Error in create3DFunc(): ', e)
                    self.errorDump('File format is not recognized or not all files go the same number of electrodes')
        
        self.create3DBtn = QPushButton('Select 2D lines')
        self.create3DBtn.clicked.connect(create3DFunc)
        self.create3DBtn.setVisible(False)


        def invNowBtnFunc():
            self.tabs.setCurrentIndex(5) # jump to inversion tab
            self.invertBtn.animateClick() # invert
        self.invNowBtn = QPushButton('Invert')
        self.invNowBtn.setStyleSheet('background-color: green')
        self.invNowBtn.setFixedWidth(150)
        self.invNowBtn.setAutoDefault(True)
        self.invNowBtn.clicked.connect(invNowBtnFunc)
        self.invNowBtn.setEnabled(False)
        self.invNowBtn.setToolTip('Invert with default settings. This will redirect you to the inversion tab.')

        def ipCheckFunc(state):
            if state  == Qt.Checked:
                self.r2.typ = 'c' + self.r2.typ
                self.typ = 'c' + self.typ
                showIpOptions(True)
                [p.setVisible(True) for p in [self.pvminIPLabel, self.pvminIP, self.pvmaxIPLabel, self.pvmaxIP]]
#                self.timeLapseCheck.setEnabled(False)
                if self.r2.iForward == True:
                    self.mwFwdPseudoIP.setVisible(True)
                    self.noiseLabelIP.show()
                    self.noiseEditIP.show()
                else:
                    self.mwPseudoIP.setVisible(True)
                    self.plotPseudoIP()
                    self.tabPreProcessing.setTabEnabled(1, True)
                    if all(self.r2.surveys[0].df['irecip'].values == 0) is False:
                        phaseplotError()
                        self.tabPreProcessing.setTabEnabled(3, True) # no reciprocity = no IP error model
                        self.recipFiltBtn.setEnabled(True)
                    heatRaw()
    #                self.r2.surveys[0].filterDataIP_plot = self.r2.surveys[0].filterDataIP_plotOrig
                    self.r2.surveys[0].filterDataIP = self.r2.surveys[0].df
                    heatFilter()
                self.regionTable.setColumnHidden(1, False)

            else:
                self.r2.typ = self.r2.typ[1:]
                self.typ = self.typ[1:]
                showIpOptions(False)
                [p.setVisible(False) for p in [self.pvminIPLabel, self.pvminIP, self.pvmaxIPLabel, self.pvmaxIP]]
#                self.timeLapseCheck.setEnabled(True)
                self.mwPseudoIP.setVisible(False)
                self.tabPreProcessing.setTabEnabled(1, False)
                self.tabPreProcessing.setTabEnabled(3, False)
                self.regionTable.setColumnHidden(1, True)
                if self.r2.iForward == True:
                    self.mwFwdPseudoIP.setVisible(False)
                    self.noiseLabelIP.hide()
                    self.noiseEditIP.hide()
            pdebug('ipCheckFunc: mode=', self.r2.typ)

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
            if self.r2.typ[0] == 'c':
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
            if self.r2.typ[0] == 'c':
                self.plotPseudoIP()
        
        # display options for pseudo-sections
        self.psContourCheck = QCheckBox('Contour')
        self.psContourCheck.stateChanged.connect(psContourFunc)
        self.psContourCheck.setEnabled(False)
        self.psContourCheck.setToolTip('Check/uncheck to contour pseudo section plots')
        
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
        
        self.pParams = {'index':0, 'vmin':None, 'vmax':None}
        self.pParamsIP = {'index':0, 'vmin':None, 'vmax':None}
        def prescaleBtnFunc():
            if self.r2 is not None:
                self.pParams['vmin'] = float(self.pvmin.text()) if self.pvmin.text() != '' else None
                self.pParams['vmax'] = float(self.pvmax.text()) if self.pvmax.text() != '' else None
                self.pParamsIP['vmin'] = float(self.pvminIP.text()) if self.pvminIP.text() != '' else None
                self.pParamsIP['vmax'] = float(self.pvmaxIP.text()) if self.pvmaxIP.text() != '' else None    
                self.plotPseudo()
                if self.r2.typ[0] == 'c':
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

        self.mwPseudo = MatplotlibWidget(navi=True, aspect='auto', itight=True)

        self.mwPseudoIP = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        self.mwPseudoIP.setVisible(False)



        # layout
        hbox1 = QHBoxLayout()
#        hbox1.addWidget(restartBtn)
        hbox1.addWidget(dimGroup)
        hbox1.addWidget(self.titleLabel)
        hbox1.addWidget(self.titleEdit)
        hbox1.addWidget(self.dateLabel)
        hbox1.addWidget(self.dateEdit)
        hbox1.addWidget(self.restartBtn)

        hbox2 = QHBoxLayout()
        hbox2.addWidget(dimInvGroup)
        hbox2.addWidget(self.timeLapseCheck)
        hbox2.addWidget(self.batchCheck)
        hbox2.addWidget(self.boreholeCheck)
        hbox2.addWidget(self.regular3DCheck)

        hbox4 = QHBoxLayout()
        hbox4.addWidget(self.wdBtn)
        hbox4.addWidget(self.ftypeComboLabel)
        hbox4.addWidget(self.ftypeCombo)
#        hbox4.addWidget(self.spacingEdit)
        hbox4.addWidget(self.importDataBtn)
        hbox4.addWidget(self.importDataRecipBtn)
        hbox4.addWidget(self.lineSpacingLabel)
        hbox4.addWidget(self.lineSpacing)
        hbox4.addWidget(self.create3DBtn)
        hbox4.addWidget(self.invNowBtn)
        
        hbox5 = QHBoxLayout()
        hbox5.setAlignment(Qt.AlignRight)
        hbox5.addWidget(self.ipCheck, Qt.AlignLeft)
        hbox5.addWidget(self.psContourCheck)
        hbox5.addWidget(self.pvminLabel)
        hbox5.addWidget(self.pvmin)
        hbox5.addWidget(self.pvmaxLabel)
        hbox5.addWidget(self.pvmax)
        hbox5.addWidget(self.pvminIPLabel)
        hbox5.addWidget(self.pvminIP)
        hbox5.addWidget(self.pvmaxIPLabel)
        hbox5.addWidget(self.pvmaxIP)
        hbox5.addWidget(self.prescaleBtn)
        hbox5.addWidget(self.fnamesComboLabel)
        hbox5.addWidget(self.fnamesCombo)

        metaLayout = QVBoxLayout()
        metaLayout.addLayout(hbox1)
        metaLayout.addLayout(hbox2)
        metaLayout.addLayout(hbox4)
        metaLayout.addLayout(hbox5)
        tabImportingDataLayout.addLayout(metaLayout, 40)

        pseudoLayout = QHBoxLayout()
        pseudoLayout.addWidget(self.mwPseudo, 50)
        pseudoLayout.addWidget(self.mwPseudoIP, 50)
        tabImportingDataLayout.addLayout(pseudoLayout, 60)
        self.tabImportingData.setLayout(tabImportingDataLayout)


        #%% sub tab with elec and topo informations
        self.tabImportingTopo = QWidget()
        self.tabImporting.addTab(self.tabImportingTopo, 'Electrodes (XYZ/Topo)')

        # electrode table
        class ElecTable(QTableWidget):
            def __init__(self, nrow=2, headers=['x','y','z','Buried'],
                         selfInit=False, parent=None):
                """ if selfInit is true, it will automatically add rows if tt
                is bigger than the actual rows
                """
                ncol = len(headers)
                super(ElecTable, self).__init__(nrow, ncol)
                self.nrow = nrow
                self.ncol = ncol
                self.parent = parent
                self.headers = headers
                self.selfInit = selfInit
                self.initTable(np.array([['',''],['',''],['','']]))
                self.horizontalHeader().sortIndicatorChanged.connect(self.setAllBuried)
                self.setHorizontalHeaderLabels(headers)
                self.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

            def addRow(self):
                self.nrow += 1
                self.setRowCount(self.nrow)
    
            def setBuried(self, vals=None):
                if vals is None:
                    vals = np.zeros(self.nrow, dtype=bool)
                for i in range(len(vals)):
                    checkBoxWidget = QWidget()
                    checkBoxLayout = QHBoxLayout()
                    checkBoxLayout.setContentsMargins(5,5,5,5)
                    checkBoxLayout.setAlignment(Qt.AlignCenter)
                    buriedCheck = QCheckBox()
                    buriedCheck.setChecked(bool(vals[i]))
                    checkBoxLayout.addWidget(buriedCheck)
                    checkBoxWidget.setLayout(checkBoxLayout)
                    self.setCellWidget(i, 3, checkBoxWidget) # column 3 == Buried

            def setAllBuried(self, colIndex):
                if self.headers[colIndex] == 'Buried':
                    for i in range(self.nrow):
                        buriedCheck = self.cellWidget(i, 3).findChildren(QCheckBox)[0]
                        if buriedCheck.isChecked() is True:
                            buriedCheck.setChecked(False)
                        else:
                            buriedCheck.setChecked(True)
                            
            def getBuried(self):
                buried = np.zeros(self.nrow, dtype=bool)
                for i in range(self.nrow):
                    buriedCheck = self.cellWidget(i, 3).findChildren(QCheckBox)[0]
                    if buriedCheck.isChecked() is True:
                        buried[i] = True
                return buried

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
                if self.isColumnHidden(1) == True: # 'y' is hidden
                    tt = np.c_[tt[:,0], np.zeros(tt.shape[0]), tt[:,1:]]
                pdebug('elecTable.paste():', tt)
                if np.sum(tt.shape) > 0:
                    if self.selfInit is True:
                        self.initTable(tt)
                    else:
                        self.setTable(tt, c0, r0)

            def initTable(self, tt):
                pdebug('elecTable.initTable():', tt)
                self.clear() # this clear out the labels as well
                self.setRowCount(tt.shape[0])
                self.setHorizontalHeaderLabels(self.headers)
                self.nrow = tt.shape[0]
                self.setTable(tt)
                if 'Buried' in self.headers:
                    self.setBuried()
            
            def setTable(self, tt, c0=0, r0=0):
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
                ncol = self.ncol if 'Buried' not in self.headers else self.ncol - 1
                table = np.zeros((self.nrow, ncol))*np.nan
                for i in range(ncol):
                    for j in range(self.nrow):
                        if self.item(j,i) is not None:
                            item = self.item(j,i).text()
                            if item != '':
                                table[j,i] = float(item)
                return table

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
                
                pdebug('elecTable.readTable():\n', df2.head())
                tt = df2.values
                if nbElec is not None:
                    if tt.shape[0] != nbElec:
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


        self.elecTable = ElecTable(headers=['x','y','z','Buried'], parent=self)
        self.elecTable.setColumnHidden(1, True)
        self.elecLabel = QLabel('<i>Add electrode position. Use <code>Ctrl+V</code> to paste or import from CSV (x,y,z header).\
                           The last column is 1 if checked (= buried electrode) and 0 if not (=surface electrode).\
                           You can also use the form below to generate \
                           regular electrode spacing. <b>Click on the <font color="red">"Buried"</font> table header to check/unchek all</b></i>')
        self.elecLabel.setWordWrap(True)

        def importElecBtnFunc():
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingTopo,'Open File', directory=self.datadir)
            if fname != '':
                if self.iForward is False:
                    nbElec = int(self.nbElecEdit.text())
                else:
                    nbElec = None
                self.elecTable.readTable(fname, nbElec=nbElec)
        self.importElecBtn = QPushButton('Import from CSV files (x, y, z headers)')
        self.importElecBtn.setAutoDefault(True)
        self.importElecBtn.clicked.connect(importElecBtnFunc)
        
        self.nbElecEdit = QLineEdit()
        self.nbElecEdit.setValidator(QDoubleValidator())
        self.nbElecEdit.setEnabled(False)
        self.nbElecLabel = QLabel('Number of electrodes:')
        
        self.elecDxLabel = QLabel('X spacing:')
        self.elecDxEdit = QLineEdit('0.0')
        self.elecDxEdit.setValidator(QDoubleValidator())

        self.elecDyEditLabel = QLabel('Y spacing:')        
        self.elecDyEdit = QLineEdit('0.0')
        self.elecDyEdit.setValidator(QDoubleValidator())
        self.elecDyEdit.setEnabled(False)

        self.elecDzLabel = QLabel('Z spacing:')
        self.elecDzEdit = QLineEdit('0.0')
        self.elecDzEdit.setValidator(QDoubleValidator())

        def elecGenButtonFunc():
            nbElec = int(self.nbElecEdit.text())
            dx = float(self.elecDxEdit.text())
            dy = float(self.elecDyEdit.text())
            dz = float(self.elecDzEdit.text())
            if 'y' in self.elecTable.headers: # 3D case
                electrodes = np.c_[np.linspace(0.0, (nbElec-1)*dx, nbElec),
                              np.linspace(0.0, (nbElec-1)*dy, nbElec),
                              np.linspace(0.0, (nbElec-1)*dz, nbElec)]
            else:
                electrodes = np.c_[np.linspace(0.0, (nbElec-1)*dx, nbElec),
                              np.linspace(0.0, (nbElec-1)*dz, nbElec)]
            self.elecTable.initTable(electrodes)
        self.elecGenButton = QPushButton('Generate')
        self.elecGenButton.setAutoDefault(True)
        self.elecGenButton.clicked.connect(elecGenButtonFunc)

        self.topoTable = ElecTable(headers=['x','y','z'], selfInit=True, parent=self)
        self.topoTable.setColumnHidden(1, True)
        topoLabel = QLabel('<i>Add additional surface points. \
                           You can use <code>Ctrl+V</code> to paste directly \
                           into a cell.</i>')
                           
        def topoBtnFunc():
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingTopo,'Open File', directory=self.datadir)
            if fname != '':
                self.topoTable.readTable(fname)
        self.topoBtn = QPushButton('Import from CSV files (x, y, z headers)')
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
        topoLayout = QVBoxLayout()

        elecGenLayout = QHBoxLayout()
        elecGenLayout.addWidget(self.nbElecLabel)
        elecGenLayout.addWidget(self.nbElecEdit)
        elecGenLayout.addWidget(self.elecDxLabel)
        elecGenLayout.addWidget(self.elecDxEdit)
        elecGenLayout.addWidget(self.elecDyEditLabel)
        elecGenLayout.addWidget(self.elecDyEdit)
        elecGenLayout.addWidget(self.elecDzLabel)
        elecGenLayout.addWidget(self.elecDzEdit)
        elecGenLayout.addWidget(self.elecGenButton)
        topoLayout.addWidget(self.elecLabel)
        topoLayout.addLayout(elecGenLayout)
        topoLayout.addWidget(self.importElecBtn)
        topoLayout.addWidget(self.elecTable)
        
        topoLayout.addWidget(topoLabel)
        topoBtnLayout = QHBoxLayout()
        topoBtnLayout.addWidget(self.topoBtn, 70)
        topoBtnLayout.addWidget(self.topoAddRowBtn, 10)
        topoBtnLayout.addWidget(self.topoInterpBtn, 20)
        topoLayout.addLayout(topoBtnLayout)
        topoLayout.addWidget(self.topoTable)

        self.tabImportingTopo.setLayout(topoLayout)
        self.tabImporting.setTabEnabled(1, False)


        #%% sub tab for custom parser
        customParser = QWidget()
        self.tabImporting.addTab(customParser, 'Custom Parser')
        
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
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingTopo,'Open File')
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
        
        
#        elecSpacingLabel = QLabel('Electrode spacing')

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
                    df = pd.read_csv(fname, delimiter=delimiter, skiprows=skiprows, nrows=nrows)
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
                    array = df[['a','b','m','n']].values.copy()
                    arrayMin = np.min(np.unique(np.sort(array.flatten())))
                    if arrayMin != 0:
                        array -= arrayMin
                    if espacing is None:
                        espacing = np.unique(np.sort(array.flatten()))[1] - np.unique(np.sort(array.flatten()))[0]
                    array = np.round(array/espacing+1).astype(int)
                    df[['a','b','m','n']] = array
                    imax = int(np.max(array))
                    elec = np.zeros((imax,3))
                    elec[:,0] = np.arange(0,imax)*espacing
                    self.nbElecEdit.setText('%s' % (len(elec)))
    #                self.nbElecEdit.setEnabled(False)
                    self.elecDxEdit.setText('%s' % (espacing))
                    return elec, df
                except:
                    self.errorDump("Import Failed: 'nan' values must be removed before importation. Use the 'Number of rows to read or skip' to remove 'nan's.")
            self.parser = parserFunc

            # test the parser
            elec, df = parserFunc(self.fnameManual)
            pdebug('custom parserFunc: shapes = ', elec.shape, df.shape)

            if (self.r2.iTimeLapse is False) & (self.r2.iBatch is False):
                self.importFile(self.fnameManual)
            self.ftypeCombo.setCurrentIndex(9)
            self.tabImporting.setCurrentIndex(0)

        self.importBtn = QPushButton('Import Dataset')
        self.importBtn.setAutoDefault(True)
        self.importBtn.clicked.connect(importBtnFunc)


        # layout
        parserLayout = QVBoxLayout()
        parserOptions = QHBoxLayout()
        columnsAssign = QGridLayout()

        parserLayout.addWidget(self.openFileBtn)
        parserOptions.addWidget(self.delimLabel)
        parserOptions.addWidget(self.delimCombo)
        parserOptions.addWidget(self.delimiterBox)
        parserOptions.addWidget(self.skipRowsLabel)
        parserOptions.addWidget(self.skipRowsEdit)
        parserOptions.addWidget(self.nrowsLabel)
        parserOptions.addWidget(self.nrowsEdit)
        parserOptions.addWidget(self.parseBtn)
        parserLayout.addLayout(parserOptions)

        parserLayout.addWidget(self.parserTable)
        for i in range(len(boxes)):
            c = (i % 3)*2 # in 2*3 columns (with labels)
            r = int(i/3)
            columnsAssign.addWidget(boxesLabels[i], r, c, Qt.AlignRight)
            columnsAssign.addWidget(boxes[i], r, c+1)
        parserLayout.addLayout(columnsAssign)
        parserLayout.addWidget(self.importBtn)

        customParser.setLayout(parserLayout)


#%% ===================== pre-processing tab =====================
        self.tabPreProcessing = QTabWidget()
        self.tabs.addTab(self.tabPreProcessing, 'Pre-processing')
        self.tabs.setTabEnabled(1, False)

        
#%% reciprocal filtering tab (or just normal filtering if no reciprocal)
        recipErrorWidget = QWidget()
        self.tabPreProcessing.addTab(recipErrorWidget, 'Reciprocal error analysis')
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
                numElecRemoved = np.sum(self.r2.surveys[index].eselect)
                # don't need filterElec as the selected points are in iselect anyway
                numSelectRemoved = 0
                if self.r2.iBatch or self.r2.iTimeLapse:
                    if not self.recipErrApplyToAll:
                        numSelectRemoved += self.r2.surveys[index].filterData(~self.r2.surveys[index].iselect)
                    else:
                        s = self.r2.surveys[index] # index should be 0 as we select on the 1st survey only
                        quads = s.df[s.iselect][['a','b','m','n']].values
                        numSelectRemoved += self.r2._filterSimilarQuad(quads) # this will remove the same quads to all surveys
                else:
                    numSelectRemoved += self.r2.surveys[0].filterData(~self.r2.surveys[0].iselect)
                if self.recipErrorInputLine.text() != '':
                    percent = float(self.recipErrorInputLine.text())
                    if self.filterAttrCombo.currentText() == 'Reciprocal Error':
                        numRecipRemoved = self.r2.filterRecip(index=self.recipErrDataIndex, percent=percent)
                    elif self.filterAttrCombo.currentText() == 'Stacking Error (Dev.)':
                        numRecipRemoved = self.r2.filterStack(index=self.recipErrDataIndex, percent=percent)
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
                            numRhoRangeRemoved = self.r2.filterTransferRes(index=self.recipErrDataIndex, vmin=vmin, vmax=vmax)
                        elif self.filterAttrCombo.currentText() == 'App. Resistivity':
                            numRhoRangeRemoved = self.r2.filterAppResist(index=self.recipErrDataIndex, vmin=vmin, vmax=vmax)
                        rhoRangeText = '%i measurements outside of the range and ' % numRhoRangeRemoved
                    if numElecRemoved != 0:
                        self.infoDump("%s%i selected electrodes and %i measurements removed!" % (rhoRangeText, numElecRemoved, numSelectRemoved))
                    else:
                        self.infoDump("%s%i selected measurements removed!" % (rhoRangeText, numSelectRemoved))
                if self.ipCheck.checkState() == Qt.Checked:
                    for s in self.r2.surveys:
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
                for s in self.r2.surveys:
                    numRestored += len(s.dfReset) - len(s.df)
                    s.df = s.dfReset.copy()
            else:
                numRestored = len(self.r2.surveys[self.recipErrDataIndex].dfReset) - len(self.r2.surveys[self.recipErrDataIndex].df)
                self.r2.surveys[self.recipErrDataIndex].df = self.r2.surveys[self.recipErrDataIndex].dfReset.copy()
            if self.recipErrorInputLine.text() != '':
                self.errHist(self.recipErrDataIndex)
                self.recipErrorInputLine.setText('')
            if self.ipCheck.checkState() == Qt.Checked:
                if self.recipErrApplyToAll:
                    for s in self.r2.surveys:
                        s.dfPhaseReset = s.dfReset.copy()
                        s.filterDataIP = s.dfReset.copy()
                else:
                    self.r2.surveys[self.recipErrDataIndex].dfPhaseReset = self.r2.surveys[self.recipErrDataIndex].dfReset.copy()
                    self.r2.surveys[self.recipErrDataIndex].filterDataIP = self.r2.surveys[self.recipErrDataIndex].dfReset.copy()
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
        # self.rhoRangeInputLabel.hide()
        
        self.rhoRangeMinInput = QLineEdit('')
        self.rhoRangeMinInput.setPlaceholderText('min')
        self.rhoRangeMinInput.setFixedWidth(80)
        self.rhoRangeMinInput.setValidator(QDoubleValidator())
        # self.rhoRangeMinInput.hide()
        
        self.rhoRangeMaxInput = QLineEdit('')
        self.rhoRangeMaxInput.setPlaceholderText('max')
        self.rhoRangeMaxInput.setFixedWidth(80)
        self.rhoRangeMaxInput.setValidator(QDoubleValidator())
        # self.rhoRangeMaxInput.hide()
        
        
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
        
        # self.filterAttrComboItems = ['Transfer Resistance', 'App. Resistivity',
        #                              'Reciprocal Error', 'Stacking Error (Dev.)']
        
        self.filterAttrComboItems = ['Transfer Resistance', 'App. Resistivity']
        # self.filterAttrCombo.addItem('Transfer Resistance')
        # self.filterAttrCombo.addItem('App. Resistivity')
        # self.filterAttrCombo.addItem('Reciprocal Error')
        # self.filterAttrCombo.addItem('Stacking Error (Dev.)')
        self.filterAttrCombo.addItems(self.filterAttrComboItems)
        # self.filterAttrCombo.currentIndexChanged.connect(filterAttrComboFunc)
        self.filterAttrCombo.activated.connect(filterAttrComboFunc)


        def recipErrorUnpairedFunc():
            pdebug('recipErrorUnpairedFunc()')
            index = -1 if self.recipErrDataIndex < 0 else self.recipErrDataIndex
            numRemoved = self.r2.filterUnpaired(index=index)
            if self.ipCheck.checkState() == Qt.Checked:
                if self.recipErrApplyToAll:
                    for s in self.r2.surveys:
                        s.dfPhaseReset = s.dfReset.copy()
                        s.filterDataIP = s.dfReset.copy()
                else:
                    self.r2.surveys[self.recipErrDataIndex].dfPhaseReset = self.r2.surveys[self.recipErrDataIndex].dfReset.copy()
                    self.r2.surveys[self.recipErrDataIndex].filterDataIP = self.r2.surveys[self.recipErrDataIndex].dfReset.copy()
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
                elec = self.elecTable.getTable() # getting the topography info
                self.r2.param['lineTitle'] = self.titleEdit.text()
                if not (self.r2.iBatch or self.r2.iTimeLapse):
                    spacing = float(self.elecDxEdit.text())
                else:
                    spacing = None
                self.r2.saveFilteredData(fname, elec, savetyp, spacing=spacing)
                
        self.recipErrorSaveBtn = QPushButton('Save data')
        self.recipErrorSaveBtn.setStyleSheet("color: green")
        self.recipErrorSaveBtn.setToolTip('This will save the data in available formats (e.g., Res2DInv.dat)')
        self.recipErrorSaveBtn.clicked.connect(saveFilteredData)
        self.recipErrorSaveBtn.setFixedWidth(150)
       
        self.mwRecipError = MatplotlibWidget(navi=True, aspect='auto', itight=True)

        self.mwManualFiltering = MatplotlibWidget(navi=True, aspect='auto', itight=True)


        # layout
        recipErrorLayout = QVBoxLayout()

        recipErrorTopLayout = QVBoxLayout()
        recipErrorLayout.addLayout(recipErrorTopLayout, 0) # number is stretch factor

        recipErrorLabelLayout = QHBoxLayout()
        recipErrorLabelLayout.addWidget(self.recipErrorLabel, 1)
        recipErrorLabelLayout.addWidget(self.recipErrorfnamesComboLabel)
        recipErrorLabelLayout.addWidget(self.recipErrorfnamesCombo)
        recipErrorTopLayout.addLayout(recipErrorLabelLayout)
        
        recipErrorInputlayout = QHBoxLayout()
        recipErrorTopLayout.addLayout(recipErrorInputlayout)
        
        recipErrorInputLeftlayout = QHBoxLayout()
        recipErrorInputLeftlayout.setAlignment(Qt.AlignLeft)
        recipErrorInputLeftlayoutL = QHBoxLayout()
        recipErrorInputLeftlayoutL.setAlignment(Qt.AlignRight)
        recipErrorInputLeftlayoutL.addWidget(self.rhoRangeInputLabel)
        recipErrorInputLeftlayoutL.addWidget(self.recipErrorInputLabel)
        # recipErrorInputLeftlayoutL.addWidget(self.tResRangeInputLabel)
        recipErrorInputLeftlayout.addLayout(recipErrorInputLeftlayoutL)
        
        recipErrorInputLineLayout = QHBoxLayout()
        recipErrorInputLineLayout.setAlignment(Qt.AlignLeft)
        recipErrorInputLineLayout.addWidget(self.rhoRangeMinInput)
        recipErrorInputLineLayout.addWidget(self.rhoRangeMaxInput)
        recipErrorInputLineLayout.addWidget(self.recipErrorInputLine)
        # recipErrorInputLineLayout.addWidget(self.tResRangeMinInput)
        # recipErrorInputLineLayout.addWidget(self.tResRangeMaxInput)
        recipErrorInputLineLayout.addWidget(self.filterAttrCombo)
        recipErrorInputLeftlayout.addLayout(recipErrorInputLineLayout)

        recipErrorInputlayout.addLayout(recipErrorInputLeftlayout)

        recipErrorBtnLayout = QHBoxLayout()
        recipErrorBtnLayout.setAlignment(Qt.AlignRight)
        recipErrorBtnLayout.addWidget(self.recipErrorUnpairedBtn)
        recipErrorBtnLayout.addWidget(self.recipErrorPltBtn)
        recipErrorBtnLayout.addWidget(self.recipErrorResetBtn)
        recipErrorBtnLayout.addWidget(self.recipErrorSaveBtn)
        recipErrorInputlayout.addLayout(recipErrorBtnLayout, 1)
               
        #tab widgets for the graphs
        self.recipErrorBottomTabs = QTabWidget()

        recipErrorBottomLayout = QVBoxLayout()
        recipErrorPlotLayout = QVBoxLayout()
        recipErrorPlotLayout.addWidget(self.mwRecipError)

        recipErrorPseudoPlotLayout = QVBoxLayout()
        recipErrorPseudoPlotLayout.addWidget(self.mwManualFiltering)

        pseudoSectionPlotTab = QWidget()
        pseudoSectionPlotTab.setLayout(recipErrorPseudoPlotLayout)
        self.recipErrorBottomTabs.addTab(pseudoSectionPlotTab, 'Pseudo Section')

        errorHistogramPlotTab = QWidget()
        errorHistogramPlotTab.setLayout(recipErrorPlotLayout)
        self.recipErrorBottomTabs.addTab(errorHistogramPlotTab, 'Error Histogram')

        recipErrorBottomLayout.addWidget(self.recipErrorBottomTabs)
        recipErrorLayout.addLayout(recipErrorBottomLayout, 1) # '1' to keep the plot in largest strech

        recipErrorWidget.setLayout(recipErrorLayout)



#%% phase filtering tab
        ipfiltWidget = QWidget()
        self.tabPreProcessing.addTab(ipfiltWidget, 'Phase Filtering')
        self.tabPreProcessing.setTabEnabled(1, False)

        phasefiltlayout = QVBoxLayout()
        phaseLabelLayout = QHBoxLayout()
        
        self.phasefiltLabel = QLabel('<b>Filter the data based on the phase/IP measurements.</b><br>\
                                Below graphs show the status of filtered data versus raw input.')
        self.phasefiltLabel.setWordWrap(True)
        self.phasefiltLabel.setAlignment(Qt.AlignLeft)
        phaseLabelLayout.addWidget(self.phasefiltLabel, 1)
        
        self.phasefiltfnamesComboLabel = QLabel('Select a dataset:')
        self.phasefiltfnamesComboLabel.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
        phaseLabelLayout.addWidget(self.phasefiltfnamesComboLabel)
        
        def phasefiltfnamesComboFunc(index):
            if index == 0:
                self.phaseFiltDataIndex = -1
            elif index > 0:
                self.phaseFiltDataIndex = index-1
            if self.r2.surveys != []:
                heatRaw()
                heatFilter()
        
        self.phasefiltfnamesCombo = QComboBox()
        self.phasefiltfnamesCombo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.phasefiltfnamesCombo.setMinimumWidth(150)
        self.phasefiltfnamesCombo.activated.connect(phasefiltfnamesComboFunc)
        phaseLabelLayout.addWidget(self.phasefiltfnamesCombo)
        phasefiltlayout.addLayout(phaseLabelLayout)

        def phirange():
            self.r2.filterRangeIP(self.phaseFiltDataIndex,
                                  float(self.phivminEdit.text()),
                                  float(self.phivmaxEdit.text()))
            heatFilter()

        def removerecip():
            self.r2.filterRecipIP(self.phaseFiltDataIndex)
            heatFilter()

        def removenested():
            self.r2.filterNested(self.phaseFiltDataIndex)
            heatFilter()

        def convFactK():
            if self.phaseFiltDataIndex == -1:
                for s in self.r2.surveys:
                    s.kFactor = float(self.phiConvFactor.text())
            else:
                self.r2.surveys[self.phaseFiltDataIndex].kFactor = float(self.phiConvFactor.text())
            heatFilter()
            heatRaw()

        phitoplayout = QHBoxLayout()
        phitoplayoutL = QHBoxLayout()
        phitoplayoutL.setAlignment(Qt.AlignLeft)
        phitoplayoutC = QHBoxLayout()
        phitoplayoutC.setAlignment(Qt.AlignRight)
        phitoplayoutR = QHBoxLayout()
        phitoplayoutR.setAlignment(Qt.AlignRight)
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

        nestedfilt = QPushButton('Remove nested')
        nestedfilt.setFixedWidth(150)
        nestedfilt.setToolTip('Measurments where M and/or N are inbetween A and B will be removed.\nNOTE: Wenner like arrays will also be affected')
        nestedfilt.setAutoDefault(True)
        nestedfilt.clicked.connect(removenested)

        phitoplayoutL.addWidget(self.phiConvFactorlabel)
        phitoplayoutL.addWidget(self.phiConvFactor)
        phitoplayoutC.addWidget(self.rangelabel)
        phitoplayoutR.addWidget(self.phivminlabel)
        phitoplayoutR.addWidget(self.phivminEdit)
        phitoplayoutR.addWidget(self.phivmaxlabel)
        phitoplayoutR.addWidget(self.phivmaxEdit)
        phitoplayoutR.addWidget(self.rangebutton)
        phitoplayoutR.addWidget(self.recipFiltBtn)
        phitoplayoutR.addWidget(nestedfilt)
        
        phitoplayout.addLayout(phitoplayoutL, 0)
        phitoplayout.addLayout(phitoplayoutC, 1)
        phitoplayout.addLayout(phitoplayoutR, 0)
        
        phasefiltlayout.addLayout(phitoplayout,0)

        def filt_reset():
            if self.phaseFiltDataIndex == -1:
                for s in self.r2.surveys:
                    s.filterDataIP = s.dfPhaseReset.copy()
                    s.df = s.dfPhaseReset.copy()
                self.infoDump('Phase filters are now reset for all datasets!')
            else:
                self.r2.surveys[self.phaseFiltDataIndex].filterDataIP = self.r2.surveys[self.phaseFiltDataIndex].dfPhaseReset.copy()
                self.r2.surveys[self.phaseFiltDataIndex].df = self.r2.surveys[self.phaseFiltDataIndex].dfPhaseReset.copy()
                self.infoDump('Phase filters are now reset for selected dataset!')
            heatFilter()
            self.dcaProgress.setValue(0)
            

        def phiCbarRange():
            if self.phaseFiltDataIndex == -1:
                for s in self.r2.surveys:
                    s.phiCbarmin = float(phiCbarminEdit.text())
                    s.phiCbarMax = float(phiCbarMaxEdit.text())
            else:
                self.r2.surveys[self.phaseFiltDataIndex].phiCbarmin = float(phiCbarminEdit.text())
                self.r2.surveys[self.phaseFiltDataIndex].phiCbarMax = float(phiCbarMaxEdit.text())
            heatFilter()
            heatRaw()

        def phiCbarDataRange():
            if self.phaseFiltDataIndex == -1:
                for s in self.r2.surveys:
                    minDataIP = np.min(s.dfOrigin['ip'])
                    maxDataIP = np.max(s.dfOrigin['ip'])
                    if self.ftype == 'ProtocolIP':
                        s.phiCbarmin = -maxDataIP
                        s.phiCbarMax = -minDataIP
                    else:
                        s.phiCbarmin = minDataIP
                        s.phiCbarMax = maxDataIP
            else:
                minDataIP = np.min(self.r2.surveys[self.phaseFiltDataIndex].dfOrigin['ip'])
                maxDataIP = np.max(self.r2.surveys[self.phaseFiltDataIndex].dfOrigin['ip'])
                if self.ftype == 'ProtocolIP':
                    self.r2.surveys[self.phaseFiltDataIndex].phiCbarmin = -maxDataIP
                    self.r2.surveys[self.phaseFiltDataIndex].phiCbarMax = -minDataIP
                else:
                    self.r2.surveys[self.phaseFiltDataIndex].phiCbarmin = minDataIP
                    self.r2.surveys[self.phaseFiltDataIndex].phiCbarMax = maxDataIP
            heatFilter()
            heatRaw()

        resetlayout = QHBoxLayout()
        resetlayoutL = QHBoxLayout()
        resetlayoutL.setAlignment(Qt.AlignLeft)
        resetlayoutR = QHBoxLayout()
        resetlayoutR.setAlignment(Qt.AlignRight)
        
        phaseSavebtn = QPushButton('Save data')
        phaseSavebtn.setStyleSheet("color: green")
        phaseSavebtn.setToolTip('This will save the data in available formats (e.g. Res2DInv.dat)')
        phaseSavebtn.clicked.connect(saveFilteredData)
        phaseSavebtn.setFixedWidth(150)
        
        filtreset = QPushButton('Reset phase filters')
        filtreset.setStyleSheet("color: red")
        filtreset.setToolTip('Reset all the filtering.\nk factor is not affected')
        filtreset.setAutoDefault(True)
        filtreset.clicked.connect(filt_reset)
        filtreset.setFixedWidth(150)
        phiCbarminlabel = QLabel('Colorbar<sub>min</sub>: ')
        phiCbarminEdit = QLineEdit()
        phiCbarminEdit.setFixedWidth(50)
        phiCbarminEdit.setValidator(QDoubleValidator())
        phiCbarminEdit.setText('0')
        phiCbarMaxlabel = QLabel('Colorbar<sub>Max</sub>: ')
        phiCbarMaxEdit = QLineEdit()
        phiCbarMaxEdit.setFixedWidth(50)
        phiCbarMaxEdit.setValidator(QDoubleValidator())
        phiCbarMaxEdit.setText('25')
        phiCbarrangebutton = QPushButton('Apply')
        phiCbarrangebutton.setFixedWidth(100)
        phiCbarrangebutton.setToolTip('This is not a filtering step.')
        phiCbarrangebutton.setAutoDefault(True)
        phiCbarrangebutton.clicked.connect(phiCbarRange)
        phiCbarDatarangebutton = QPushButton('Raw data range')
        phiCbarDatarangebutton.setToolTip('This is not a filtering step.')
        phiCbarDatarangebutton.setAutoDefault(True)
        phiCbarDatarangebutton.clicked.connect(phiCbarDataRange)
        phiCbarDatarangebutton.setFixedWidth(150)
        resetlayoutL.addWidget(phiCbarminlabel)
        resetlayoutL.addWidget(phiCbarminEdit)
        resetlayoutL.addWidget(phiCbarMaxlabel)
        resetlayoutL.addWidget(phiCbarMaxEdit)
        resetlayoutL.addWidget(phiCbarrangebutton)
        resetlayoutL.addWidget(phiCbarDatarangebutton)
        resetlayoutR.addWidget(filtreset)
        resetlayoutR.addWidget(phaseSavebtn)
        
        resetlayout.addLayout(resetlayoutL, 0)
        resetlayout.addLayout(resetlayoutR, 1)
#        self.recipFiltBtn.clicked.connect("add function")


        ipfiltlayout = QHBoxLayout()

        def heatRaw():
            if self.phaseFiltDataIndex == -1:
                index = 0
            else:
                index = self.phaseFiltDataIndex
            self.r2.surveys[index].filt_typ = 'Raw'
            raw_hmp.setCallback(self.r2.showHeatmap)
            raw_hmp.replot(index=index)

        def heatFilter():
            if self.phaseFiltDataIndex == -1:
                index = 0
            else:
                index = self.phaseFiltDataIndex
            self.r2.surveys[index].filt_typ = 'Filtered'
            filt_hmp.setCallback(self.r2.showHeatmap)
            filt_hmp.replot(index=index)

        raw_hmp = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        filt_hmp = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        ipfiltlayout.addWidget(raw_hmp)
        ipfiltlayout.addWidget(filt_hmp)


        def dcaDump(val):
            self.dcaProgress.setValue(val)
            QApplication.processEvents()

        def dcaFiltering():
            try:
                self.dcaButton.setEnabled(False) # prevent multiple clicks!
                if len(self.r2.surveys) > 1 and self.phaseFiltDataIndex == -1: # well, DCA is long!!
                    self.infoDump('DCA will run %s times. Please be patient!' % len(self.r2.surveys))
                self.r2.filterDCA(index=self.phaseFiltDataIndex, dump=dcaDump)
                heatFilter()
                self.dcaButton.setEnabled(True)
            except:
                self.errorDump('No decay curves found or incomplete set of decay curves! Export the data from "Prosys" with M1, M2, ... , M20 and TM1 tabs enabled.')
                self.dcaButton.setEnabled(True)

        dcaLayout = QHBoxLayout()
        self.dcaButton = QPushButton('DCA filtering')
        self.dcaButton.setToolTip('Decay Curve Analysis filtering.\nFor more see: Flores Orozco, et al. (2017), Decay curve analysis for data error quantification in\ntime-domain induced polarization imaging')
        self.dcaButton.setAutoDefault(True)
        self.dcaButton.clicked.connect(dcaFiltering)
        self.dcaProgress = QProgressBar()
        self.dcaButton.setEnabled(False)
        self.dcaProgress.setEnabled(False)
        dcaLayout.addWidget(self.dcaButton)
        dcaLayout.addWidget(self.dcaProgress)

        phasefiltlayout.addLayout(dcaLayout, 1)
        phasefiltlayout.addLayout(resetlayout, 2)
        phasefiltlayout.addLayout(ipfiltlayout, 3)


        # layout
        #TODO tidy up and put all layout here
        
        ipfiltWidget.setLayout(phasefiltlayout)


#%% resistance error modelling tab
        errorWidget = QWidget()
        self.tabPreProcessing.addTab(errorWidget, 'Resistance Error Model')
        self.tabPreProcessing.setTabEnabled(2, False)
        
        errFitLabel = QLabel('Select an error model from the drop-down menu. Once\
                             fitted, the model will generate an error for each quadrupoles\
                             (even the ones with no reciprocals). This error will\
                             be written in the <code>protocol.dat</code> file \
                             and used in the inversion if both <code>a_wgt</code> and\
                             <code>b_wgt</code> are both set to 0 (see \'Inversion settings\' tab).')
        errFitLabel.setWordWrap(True)
        errFitLabel.setToolTip('In case of batch/time-lapse inversion, <i>all</i> datesets must either have an error model \
                               or not have any error models (i.e., select separate error models for each individual dataset or "Apply to all"). \
                               ResIPy can handle batch data with mixture of different error models.')
        errFitLabel.setAlignment(Qt.AlignLeft)
        
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
                self.errFitTypeFunc(self.errFitPlotIndexList[index-2])
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
            if len(self.r2.surveys) == 0:
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
                self.mwFitError.setCallback(self.r2.fitErrorLin)
                self.mwFitError.replot(index=self.errFitDataIndex)
                self.r2.err = True
            elif index == 2:
                self.mwFitError.setCallback(self.r2.fitErrorPwl)
                self.mwFitError.replot(index=self.errFitDataIndex)
                self.r2.err = True
            elif index == 3:
                self.mwFitError.setCallback(self.r2.fitErrorLME)
                self.mwFitError.replot(index=self.errFitDataIndex)
                self.r2.err = True
                
            # record the type of fit for each survey
            if self.errFitDataIndex == -1: # same model for each
                self.errFitPlotIndexList = [index]*len(self.r2.surveys)
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
                self.r2.saveErrorData(fname)                
        saveErrBtn = QPushButton('Save Error Data')
        saveErrBtn.setStyleSheet("color: green")
        saveErrBtn.setFixedWidth(150)
        saveErrBtn.clicked.connect(saveErrBtnFunc)
        saveErrBtn.setToolTip('Save error data for DC and IP (if available) as .csv')
        
        # layout
        errorLayout = QVBoxLayout()
        errorLayout.setAlignment(Qt.AlignTop)
        
        errorTopLayout = QHBoxLayout()
        errorTopLayout.addWidget(errFitLabel, 1)
        errorTopLayout.addWidget(self.errFitfnamesComboLabel)
        errorTopLayout.addWidget(self.errFitfnamesCombo)
        errorLayout.addLayout(errorTopLayout)
        
        errFitLayout = QHBoxLayout()
        errFitLayout.addWidget(self.errFitType, 70)
        errFitLayout.addWidget(saveErrBtn, 30)
        errorLayout.addLayout(errFitLayout)
        
        errorPlotLayout = QVBoxLayout()
        errorPlotLayout.addWidget(self.mwFitError)
        errorLayout.addLayout(errorPlotLayout, 1)

        errorWidget.setLayout(errorLayout)
        

#%% IP error model tab
        ipWidget = QWidget()
        self.tabPreProcessing.addTab(ipWidget, 'Phase Error Model')
        self.tabPreProcessing.setTabEnabled(3, False)

        iperrFitLabel = QLabel('Select an error model from the drop-down menu. Once\
                     fitted, the model will generate an error for each quadrupoles\
                     (even the ones with no reciprocals). This error will\
                     be written in the <code>protocol.dat</code> file \
                     and used in the inversion if both <code>a_wgt</code> and\
                     <code>b_wgt</code> are both set to 0 (see \'Inversion settings\' tab).')
        iperrFitLabel.setWordWrap(True)
        iperrFitLabel.setToolTip('In case of batch/time-lapse inversion, <i>all</i> datesets must either have an error model \
                                 or not have any error models (i.e., select separate error models for each individual dataset or "Apply to all"). \
                                 ResIPy can handle batch data with mixture of different error models.')
        iperrFitLabel.setAlignment(Qt.AlignLeft)
        
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
                self.iperrFitTypeFunc(self.iperrFitPlotIndexList[index-2])
        self.iperrFitfnamesCombo = QComboBox()
        self.iperrFitfnamesCombo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.iperrFitfnamesCombo.setMinimumWidth(150)
        self.iperrFitfnamesCombo.activated.connect(iperrFitfnamesComboFunc)
        
        self.mwIPFitError = MatplotlibWidget(navi=True, aspect='auto', itight=True)

        def phaseplotError(index=0):
            if len(self.r2.surveys) == 0:
                return
            self.mwIPFitError.setCallback(self.r2.showErrorIP)
            self.mwIPFitError.replot(index=index)
            self.r2.err = False

        def iperrFitTypeFunc(index):
            if len(self.r2.surveys) == 0:
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
                self.mwIPFitError.setCallback(self.r2.fitErrorPwlIP)
                self.mwIPFitError.replot(index=self.iperrFitDataIndex)
                self.r2.err = True
            elif index == 2:
                self.mwIPFitError.setCallback(self.r2.fitErrorParabolaIP)
                self.mwIPFitError.replot(index=self.iperrFitDataIndex)
                self.r2.err = True
                
            # record the type of fit for each survey
            if self.iperrFitDataIndex == -1: # same model for each
                self.iperrFitPlotIndexList = [index]*len(self.r2.surveys)
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
        ipLayout = QVBoxLayout()
        ipLayout.setAlignment(Qt.AlignTop)
        
        ipTopLayout = QHBoxLayout()
        ipTopLayout.addWidget(iperrFitLabel, 1)
        ipTopLayout.addWidget(self.iperrFitfnamesComboLabel)
        ipTopLayout.addWidget(self.iperrFitfnamesCombo)
        ipLayout.addLayout(ipTopLayout)

        errIPFitLayout = QHBoxLayout()
        errIPFitLayout.addWidget(self.iperrFitType, 70)
        errIPFitLayout.addWidget(self.saveIPErrBtn, 30)        
        ipLayout.addLayout(errIPFitLayout)
        
        ipErrPlotLayout = QVBoxLayout()
        ipErrPlotLayout.addWidget(self.mwIPFitError)
        ipLayout.addLayout(ipErrPlotLayout,1)

        ipWidget.setLayout(ipLayout)


#%% additional actions for pre-processing tab
        
        def errorCombosFill(comboboxes=[]): #filling pre-processing comboboxes with fnames
            for comboboxe in comboboxes:
                comboboxe.clear()
            
            for comboboxe in comboboxes:
                comboboxe.addItem('Apply to Each')
            
            for comboboxe in comboboxes[2:]: # only for error modelling
                comboboxe.addItem('Combined Fit')
            
            for s in self.r2.surveys:
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
        tabMesh= QWidget()
        self.tabs.addTab(tabMesh, 'Mesh')
        self.tabs.setTabEnabled(2, False)
        
        
        fmdToolTip = 'Boundary depth (meters) below which the coarse (or background) mesh starts.\n'\
                     'Default = 2/3 max dipole length.'
        fmdLabel = QLabel('Fine/Coarse\nboundary depth')
        fmdLabel.setToolTip(fmdToolTip)
        fmdLabel.setAlignment(Qt.AlignCenter)
        fmdBox = QLineEdit()
        fmdBox.setAlignment(Qt.AlignCenter)
        fmdBox.setPlaceholderText('[m]')
        fmdBox.setToolTip(fmdToolTip)
        fmdBox.setValidator(QDoubleValidator())
        
        
        def replotMesh():
            if self.iDesign is False:
                self.regionTable.reset()
                self.iDesign is False
            def func(ax):
                self.r2.createModel(ax=ax, addAction=self.regionTable.addRow)
            self.calcAspectRatio()
            self.mwMesh.plot(func, aspect = self.plotAspect)
            self.mwMesh.canvas.setFocusPolicy(Qt.ClickFocus) # allows the keypressevent to go to matplotlib
            self.mwMesh.canvas.setFocus() # set focus on the canvas

        
        def designModel():
            self.iDesign = True
            # read electrodes locations
            elec = self.elecTable.getTable()
            if np.sum(~np.isnan(elec[:,0])) == 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            else:
                self.r2.setElec(elec)
            
            # plot the interactive model
            self.regionTable.reset()
            def func(ax):
                self.r2.designModel(ax=ax, addAction=self.regionTable.addRow)
            self.mwMesh.plot(func, aspect = self.plotAspect)
            self.mwMesh.canvas.setFocusPolicy(Qt.ClickFocus) # allows the keypressevent to go to matplotlib
            self.mwMesh.canvas.setFocus() # set focus on the canvas
            # if self.iForward is False:
            #     self.regionTable.setColumnHidden(2, False) # show zone column
            #     self.regionTable.setColumnHidden(3, False) # show fixed column
            meshOutputStack.setCurrentIndex(1)

        self.designModelBtn = QPushButton('Design Model before meshing')
        self.designModelBtn.clicked.connect(designModel)

        def meshQuadFunc():
#            self.cropBelowFmd.setChecked(False)
            self.cropBelowFmd.setEnabled(True)
            self.regionTable.reset()
            elec = self.elecTable.getTable()
            if np.sum(~np.isnan(elec[:,0])) == 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            else:
                self.r2.setElec(elec)
            buried = self.elecTable.getBuried()
            surface = self.topoTable.getTable()
            inan = ~np.isnan(surface[:,0])
            if not np.isnan(surface[0,0]):
                surface_flag = True
            else:
                surface_flag = False
            if np.sum(~inan) == surface.shape[0]:
                surface = None
            else:
                surface = surface[inan,:]
#            if not np.isnan(surface[0,0]):
#                surface_flag = True
#            else:
#                surface_flag = False
            nnodes = nnodesSld.value()
            try:
                fmd = np.abs(float(fmdBox.text())) if fmdBox.text() != '' else None
                if surface_flag:
                    print("quad mesh + topo")
                    self.r2.createMesh(typ='quad', buried=buried, elemx=nnodes, surface=surface, fmd=fmd)
                else:
                    print("quad mesh no topo")
                    self.r2.createMesh(typ='quad', buried=buried, elemx=nnodes, fmd=fmd)
                self.scale.setVisible(False)
                self.scaleLabel.setVisible(False)
                meshOutputStack.setCurrentIndex(1)
                replotMesh()
                # if self.iForward is False:
                #     self.regionTable.setColumnHidden(2, False) # hide zone column
                #     self.regionTable.setColumnHidden(3, False) # hide fixed column
            except Exception as e:
                self.errorDump('Error creating the mesh: ' + str(e))
        self.meshQuad = QPushButton('Quadrilateral Mesh')
        self.meshQuad.setAutoDefault(True)
        self.meshQuad.clicked.connect(meshQuadFunc)
        self.meshQuad.setToolTip('Generate quadrilateral mesh.')

        def meshTrianFunc():
#            self.cropBelowFmd.setChecked(True)
            self.cropBelowFmd.setEnabled(True)
            elec = self.elecTable.getTable()
            if np.sum(~np.isnan(elec[:,0])) == 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            else:
                self.r2.setElec(elec)
            meshOutputStack.setCurrentIndex(0)
            QApplication.processEvents()
            self.meshLogText.clear()
            elecSpacing = np.sqrt((self.r2.elec[~self.r2.iremote,:][0,0]-self.r2.elec[~self.r2.iremote,:][1,0])**2+
                                  (self.r2.elec[~self.r2.iremote,:][0,2]-self.r2.elec[~self.r2.iremote,:][1,2])**2)
            cl = float(clSld.value())/10*(elecSpacing-elecSpacing/8)
            cl_factor = clFactorSld.value()
            buried = self.elecTable.getBuried()
            surface = self.topoTable.getTable()
            inan = ~np.isnan(surface[:,0])
            refine = 1 if refineTrianCheck.isChecked() else 0
            if np.sum(~inan) == surface.shape[0]:
                surface = None
            else:
                surface = surface[inan,:]
            fmd = np.abs(float(fmdBox.text())) if fmdBox.text() != '' else None
            self.r2.createModelMesh(buried=buried, surface=surface,
                                    cl=cl, cl_factor=cl_factor, show_output=True,
                                    dump=meshLogTextFunc, refine=refine, fmd=fmd)
            # if self.iForward is False:
            #     self.regionTable.setColumnHidden(2, False) # show zone column
            #     self.regionTable.setColumnHidden(3, False) # show fixed column
            self.scale.setVisible(True)
            self.scaleLabel.setVisible(True)
            meshOutputStack.setCurrentIndex(1)
            replotMesh()
        meshTrian = QPushButton('Triangular Mesh')
        meshTrian.setAutoDefault(True)
        meshTrian.setFocus()
        meshTrian.clicked.connect(meshTrianFunc)
        meshTrian.setToolTip('Generate triangular mesh.')


        def meshTetraFunc():
            self.cropBelowFmd.setChecked(False) # TODO: come back here and see if maxDepth works on 3D
            self.cropBelowFmd.setEnabled(False)
            elec = self.elecTable.getTable()
            topo = self.topoTable.getTable()
            inan = ~np.isnan(topo[:,0])
            if np.sum(~np.isnan(elec[:,0])) == 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            elif all(elec[:,1] == 0) & all(topo[inan,1] == 0):
                self.errorDump('For 3D meshes, Y coordinates must be supplied for topo or elec at least.')
                return
            else:
                self.r2.setElec(elec)
            meshOutputStack.setCurrentIndex(0)
            QApplication.processEvents()
            self.meshLogText.clear()
            cl = float(cl3Edit.text())
            cl_factor = float(cl3FactorEdit.text())
            cln_factor = float(clnFactorEdit.text()) if clnFactorEdit.text() != '' else 100
            buried = self.elecTable.getBuried()
            refine = 1 if refineTetraCheck.isChecked() else 0
            if np.sum(~inan) == topo.shape[0]:
                topo = None
            else:
                topo = topo[inan,:]
            
            fmd = np.abs(float(fmdBox.text())) if fmdBox.text() != '' else None
            self.r2.createMesh(typ='tetra', buried=buried, surface=topo, fmd=fmd,
                               cl=cl, cl_factor=cl_factor, dump=meshLogTextFunc,
                               cln_factor=cln_factor, refine=refine, show_output=True)
            if pvfound:
                mesh3Dplotter.clear() # clear all actors 
                self.r2.showMesh(ax=mesh3Dplotter, color_map='Greys', color_bar=False)
            else:
                self.mwMesh3D.plot(self.r2.showMesh, threed=True)
            meshOutputStack.setCurrentIndex(2)

        meshTetra = QPushButton('Tetrahedral Mesh')
        meshTetra.setAutoDefault(True)
        meshTetra.clicked.connect(meshTetraFunc)
        meshTetra.setToolTip('Generate tetrahedral mesh.')


        # additional options for quadrilateral mesh
        nnodesLabel = QLabel('Number of elements\nbetween electrodes')
        nnodesSld = QSlider(Qt.Horizontal)
        nnodesSld.setMinimumWidth(50)
        nnodesSld.setMinimum(1)
        nnodesSld.setMaximum(10)
        nnodesSld.setValue(4)
        nnodesGrid = QGridLayout()
        # nnodesGrid.setContentsMargins(9,0,9,0)
        nnodesGrid.setSpacing(2)
        nnodes1Label = QLabel('1')
        nnodes10Label = QLabel('10')
        nnodes10Label.setAlignment(Qt.AlignRight)#| Qt.AlignVCenter)
        nnodesGrid.addWidget(nnodesSld, 0, 0, 1, 2)
        nnodesGrid.addWidget(nnodes1Label, 1,0,1,1)
        nnodesGrid.addWidget(nnodes10Label, 1,1,1,1)
        

        # additional options for triangular mesh
        clLabel = QLabel('Characteristic Length:')
        clLabel.setToolTip('Control the number and size of elements between electrodes.')
        clGrid = QGridLayout()
        clGrid.setContentsMargins(9,0,9,0)
        clGrid.setSpacing(2)
        clFineLabel = QLabel('Fine')
        clFineLabel.setStyleSheet('font:12px;')
        clCoarseLabel = QLabel('Coarse')
        clCoarseLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        clCoarseLabel.setStyleSheet('font:12px;')
        clSld = QSlider(Qt.Horizontal)
        clSld.setMinimum(1) # this depends on electrode spacing
        clSld.setMaximum(10)
        clSld.setValue(5)
        clGrid.addWidget(clSld, 0, 0, 1, 2)
        clGrid.addWidget(clFineLabel, 1,0,1,1)
        clGrid.addWidget(clCoarseLabel, 1,1,1,1)
        clFactorLabel = QLabel('Growth factor:')
        clFactorLabel.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        clFactorLabel.setToolTip('Factor by which elements grow away from electrodes.')
        clFactorSld = QSlider(Qt.Horizontal)
        clFactorSld.setMaximumWidth(200)
        clFactorSld.setMinimum(1)
        clFactorSld.setMaximum(10)
        clFactorSld.setValue(4)
        refineTrianCheck = QCheckBox('Refine')
        refineTrianCheck.setToolTip('Refine the mesh for forward modelling'
                                    'without increasing the number of parameters.')

        # additional options for tetrahedral mesh
        cl3Label = QLabel('Characteristic Length:')
        cl3Edit = QLineEdit()
        cl3Edit.setValidator(QDoubleValidator())
        cl3Edit.setText('-1')
        cl3FactorLabel = QLabel('Growth factor Top:')
        cl3FactorEdit = QLineEdit()
        cl3FactorEdit.setValidator(QDoubleValidator())
        cl3FactorEdit.setText('8')       
        clnFactorLabel = QLabel('Growth factor Bottom:')
        clnFactorEdit = QLineEdit()
        clnFactorEdit.setValidator(QDoubleValidator())
        clnFactorEdit.setText('100')
        refineTetraCheck = QCheckBox('Refine')
        refineTetraCheck.setToolTip('Refine the mesh for forward modelling'
                                    'without increasing the number of parameters.')
        
        def saveMeshVtkBtnFunc():
            fname, _ = QFileDialog.getSaveFileName(tabMesh, 'Open File', self.datadir)
            if fname != '':
                self.r2.saveMeshVtk(fname)
                self.infoDump('Mesh saved to {:s}'.format(fname))
        saveMeshVtkBtn = QPushButton('Save Mesh as .vtk')
        saveMeshVtkBtn.clicked.connect(saveMeshVtkBtnFunc)

        def importCustomMeshFunc():
            elec = self.elecTable.getTable()
            if np.sum(~np.isnan(elec[:,0])) == 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            else:
                self.r2.setElec(elec)
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingData,'Open File', self.datadir)
            if fname != '':
                try:
                    self.r2.importMesh(fname)
                    print('mesh imported ... now displaying ... ')
                    if pvfound:
                        mesh3Dplotter.clear() # clear all actors 
                        self.r2.showMesh(ax=mesh3Dplotter, color_map='Greys', color_bar=False)
                    else:
                        self.mwMesh3D.plot(self.r2.showMesh, threed=True)
                    meshOutputStack.setCurrentIndex(2)
                except Exception as e:
                    self.errorDump('Error importing mesh' + str(e))
        self.importCustomMeshBtn = QPushButton('Import Custom Mesh')
        self.importCustomMeshBtn.setFixedWidth(150)
        self.importCustomMeshBtn.clicked.connect(importCustomMeshFunc)


        importCustomMeshLabel2 = QLabel('Import .msh or .vtk file.')
        importCustomMeshLabel2.setAlignment(Qt.AlignCenter)
        importCustomMeshLabel2.setWordWrap(True)
        def importCustomMeshFunc2():
            print('using importCustomMeshFunc2')
            elec = self.elecTable.getTable()
            if np.sum(~np.isnan(elec[:,0])) == 0:
                self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                return
            else:
                self.r2.setElec(elec)
            fname, _ = QFileDialog.getOpenFileName(self.tabImportingData,'Open File', self.datadir)
            if fname != '':
                try:
                    self.r2.importMesh(fname, mesh_type='trian')
                    if (self.r2.typ == 'R3t') or (self.r2.typ == 'cR3t'):
                        self.mwMesh.plot(self.r2.showMesh, threed=True)
                        meshOutputStack.setCurrentIndex(2)
                        print('NO WAY THIS CAN HAPPEN!')
                    else:
                        replotMesh()
                        meshOutputStack.setCurrentIndex(1)
                    self.regionTable.reset()
                    for i in range(len(np.unique(self.r2.mesh.attr_cache['region']))-1):
                        self.regionTable.addRow()
                except Exception as e:
                    self.errorDump('Error importing mesh' + str(e))
        self.importCustomMeshBtn2 = QPushButton('Import Custom Mesh')
        self.importCustomMeshBtn2.setFixedWidth(150)
        self.importCustomMeshBtn2.clicked.connect(importCustomMeshFunc2)
        self.importCustomMeshBtn2.setToolTip('Import .msh or .vtk file. The electrodes will be snapped to the closest node.')

        def resetMeshBtnFunc():
            self.regionTable.reset()
            self.mwMesh.clear()
            self.r2.mesh = None
            self.r2.geom_input = {}
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
                checkBoxWidget = QWidget()
                checkBoxLayout = QHBoxLayout()
                checkBoxLayout.setContentsMargins(5,5,5,5)
                checkBoxLayout.setAlignment(Qt.AlignCenter)
                checkBoxLayout.addWidget(QCheckBox())
                checkBoxWidget.setLayout(checkBoxLayout)
                self.setCellWidget(0,3, checkBoxWidget)
                

            def addRow(self):
                self.nrow = self.nrow + 1
                self.setRowCount(self.nrow)
                self.setItem(self.nrow-1, 0, QTableWidgetItem('100.0'))
                self.setItem(self.nrow-1, 1, QTableWidgetItem('0'))
                self.setItem(self.nrow-1, 2, QTableWidgetItem('1'))
                checkBoxWidget = QWidget()
                checkBoxLayout = QHBoxLayout()
                checkBoxLayout.setContentsMargins(5,5,5,5)
                checkBoxLayout.setAlignment(Qt.AlignCenter)
                checkBoxLayout.addWidget(QCheckBox())
                checkBoxWidget.setLayout(checkBoxLayout)
                self.setCellWidget(self.nrow-1, 3, checkBoxWidget)

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


        instructionLabel = QLabel('Click on the buttons below to define a region.'
            ' If polyline is selected, you can close the region by using the right'
            ' click of the mouse. Note that it is possible to define a zone (so'
            ' no smoothing at the boundaries) and fix the region value for '
            'triangular mesh only.')
        instructionLabel.setWordWrap(True)

        self.mwMesh = MatplotlibWidget(navi=True, itight=False)
        self.mwMesh3D = MatplotlibWidget(threed=True, navi=True)
        if pvfound:
            mesh3Dplotter = pv.QtInteractor()
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
        meshLayout = QVBoxLayout()

        meshOptionQuadLayout = QHBoxLayout()
        meshOptionQuadLayout.setAlignment(Qt.AlignVCenter)
        meshOptionQuadLayout.addWidget(nnodesLabel)
        meshOptionQuadLayout.addLayout(nnodesGrid, 1)

        meshOptionTrianLayout = QHBoxLayout()
        meshOptionTrianLayout.addWidget(clLabel)
        meshOptionTrianLayout.addLayout(clGrid, 1)
        meshOptionTrianLayout.addWidget(clFactorLabel, 1)
        meshOptionTrianLayout.addWidget(clFactorSld, 1)
        
        meshButtonTrianLayout = QHBoxLayout()
        meshButtonTrianLayout.addWidget(refineTrianCheck)
        meshButtonTrianLayout.addWidget(self.designModelBtn)
        meshButtonTrianLayout.addWidget(meshTrian)
        
        importCustomLayout = QVBoxLayout()
        importCustomLayout.addWidget(importCustomMeshLabel2)
        importCustomLayout.addWidget(self.importCustomMeshBtn2)

        meshOptionTetraLayout = QHBoxLayout()
        meshOptionTetraLayout.addWidget(cl3Label)
        meshOptionTetraLayout.addWidget(cl3Edit)
        meshOptionTetraLayout.addWidget(cl3FactorLabel)
        meshOptionTetraLayout.addWidget(cl3FactorEdit)
        meshOptionTetraLayout.addWidget(clnFactorLabel)
        meshOptionTetraLayout.addWidget(clnFactorEdit)
        meshOptionTetraLayout.addWidget(refineTetraCheck)
        meshOptionTetraLayout.addWidget(saveMeshVtkBtn)
        meshOptionTetraLayout.addWidget(self.importCustomMeshBtn)
        
        meshChoiceLayout = QHBoxLayout()
        fmdLayout = QVBoxLayout()
        fmdLayout.setAlignment(Qt.AlignBottom | Qt.AlignCenter)
        meshQuadLayout = QVBoxLayout()
        meshQuadLayout.setAlignment(Qt.AlignBottom | Qt.AlignHCenter)
        meshTrianLayout = QVBoxLayout()
        meshTetraLayout = QVBoxLayout()
        meshTetraLayout.setAlignment(Qt.AlignBottom | Qt.AlignHCenter)
        fmdGroup = QGroupBox()
        fmdGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        self.meshQuadGroup = QGroupBox()
        self.meshQuadGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        self.meshTrianGroup = QGroupBox()
        self.meshTrianGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        meshTetraGroup = QGroupBox()
        meshTetraGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        
        fmdLayout.addWidget(fmdLabel)
        fmdLayout.addWidget(fmdBox)
        fmdGroup.setLayout(fmdLayout)
        meshChoiceLayout.addWidget(fmdGroup,0)
        
        meshQuadLayout.addLayout(meshOptionQuadLayout, 1)
        meshQuadLayout.addWidget(self.meshQuad, 0)
        self.meshQuadGroup.setLayout(meshQuadLayout)
        meshChoiceLayout.addWidget(self.meshQuadGroup,35)

        meshTrianLayout.addLayout(meshOptionTrianLayout)
        meshTrianLayout.addLayout(meshButtonTrianLayout)
        self.meshTrianGroup.setLayout(meshTrianLayout)
        meshChoiceLayout.addWidget(self.meshTrianGroup,65)

        meshCustomGroup = QGroupBox()
        meshCustomGroup.setStyleSheet("QGroupBox{padding-top:1em; margin-top:-1em}")
        meshCustomGroup.setLayout(importCustomLayout)
        meshChoiceLayout.addWidget(meshCustomGroup,0)
        
        meshTetraLayout.addLayout(meshOptionTetraLayout)
        meshTetraLayout.addWidget(meshTetra)
        meshTetraGroup.setLayout(meshTetraLayout)
        meshChoiceLayout.addWidget(meshTetraGroup)
        meshTetraGroup.setHidden(True)

        meshLayout.addLayout(meshChoiceLayout, 0)

        instructionLayout = QHBoxLayout()
        instructionLayout.addWidget(instructionLabel, 92)
        instructionLayout.addWidget(self.resetMeshBtn, 8)
        meshLayout.addLayout(instructionLayout)

        regionLayout = QVBoxLayout()
        regionLayout.addWidget(self.regionTable)

        meshPlot = QWidget()
        meshPlotLayout = QHBoxLayout()
        meshPlotLayout.addWidget(self.mwMesh, 70)
        meshPlotLayout.addLayout(regionLayout, 30)
        meshPlot.setLayout(meshPlotLayout)

        meshPlot3D = QWidget()
        meshPlot3DLayout = QHBoxLayout()
        meshPlot3DLayout.addWidget(self.mwMesh3D)
        if pvfound:
            meshPlot3DLayout.addWidget(mesh3Dplotter.interactor)
        meshPlot3D.setLayout(meshPlot3DLayout)

        meshOutputStack = QStackedLayout()
        meshOutputStack.addWidget(self.meshLogText)
        meshOutputStack.addWidget(meshPlot)
        meshOutputStack.addWidget(meshPlot3D)
        meshOutputStack.setCurrentIndex(0)

        meshLayout.addLayout(meshOutputStack, 1)

        tabMesh.setLayout(meshLayout)


#%% =================== Tab for forward model ===================
        tabForward = QWidget()
        self.tabs.addTab(tabForward, 'Forward model')
        self.tabs.setTabEnabled(3, False)
        
        fwdSlpitterLayout = QHBoxLayout() # a splitter so the graphs are easier to visualize
        fwdSplitter = QSplitter(Qt.Vertical)

        seqLabel = QLabel('Design a sequence for the forward modelling. A \
combination of multiple sequence is accepted as well as importing a custom sequence')
        seqLabel.setWordWrap(True)
        seqLabel.setAlignment(Qt.AlignTop)
        
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
                self.combo.deleteLater()
                for w in self.fields:
                    w.deleteLater()
                for w in self.labels:
                    w.deleteLater()
                self.rmBtn.deleteLater()
                self.importBtn.deleteLater()
                self.deleteLater()
                if self.seq == 'custSeq':
                    self.fname = ''
            
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


        seqRowLayout = QVBoxLayout()
        seqRowLayout.setAlignment(Qt.AlignTop)
        seqRows = []
        seqRow = RowOpt(parent=self)
        seqRowLayout.addLayout(seqRow)
        seqRows.append(seqRow)
        
        DpDp = resource_path('image/dipdip.png')
        Wenner = resource_path('image/wenner.png')
        Schlum = resource_path('image/schlum.png')
        Gradient = resource_path('image/gradient.png')
        
        seqHelp = {'dpdp1' : '<img height=140 src="%s">' % DpDp,
           'wenner': '<img height=140 src="%s">' % Wenner,
           'schlum1': '<img height=140 src="%s">' % Schlum,
           'multigrad': '<img height=140 src="%s">' % Gradient,
           'custSeq': 'Use the button to import a custom CSV file (with headers)\ncolumn1: C+, column2: C-, column3: P+, column4: P-'
            }
        
        arrayLabel = QLabel('Sequence help will be display here.')
        arrayLabel.setAlignment(Qt.AlignCenter)
        def showArray(arg):
            if arg not in seqHelp.keys():
                arrayLabel.setText('Sequence help not found.')
            else:
                arrayLabel.setText(seqHelp[arg])
        showArray('dpdp1') # default
        
        def addRowBtnFunc():
            a = RowOpt(parent=self)
            seqRows.append(a)
            seqRowLayout.addLayout(a)
            a.createRow()
            a.showArg()
        addRowBtn = QPushButton('Add sequence')
        addRowBtn.adjustSize()
        addRowBtn.clicked.connect(addRowBtnFunc)
        
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
                except:# Exception as e:
                    pass
#                    print('object does not exist', e)
            print('sequences = ', vals)
            return vals
#        getDataBtn = QPushButton('Get Data')
#        getDataBtn.clicked.connect(getDataBtnFunc)
        
        def seqCreateFunc():
            if self.r2.elec is None:
                self.errorDump('Input electrode positions in the "Electrodes (XYZ/Topo)" tab first.')
                return
            else:
                self.r2.setElec(self.elecTable.getTable())
            params = getDataBtnFunc()
            if len(params) == 0:
                raise ValueError('You must specify at least one sequence.')
                return
            self.r2.createSequence(params=params)
            self.seqOutputLabel.setText('{:d} quadrupoles in total'.format(len(self.r2.sequence)))

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
                self.r2.saveSequence(fname)
        self.saveSeqBtn = QPushButton('Save Sequence')
        self.saveSeqBtn.setToolTip('This will save the sequence of the fwd modeling. Output data is already saved in <i>fwd</i> folder in the <i>working directory</i>.')
        self.saveSeqBtn.clicked.connect(saveSeqBtnFunc)

        # add a forward button
        def forwardBtnFunc():
            if self.r2.mesh is None: # we need to create mesh to assign starting resistivity
                self.errorDump('Please specify a mesh and an initial model first.')
                return
            try:
                seqCreateFunc()
            except:
                self.errorDump('Error in sequence generation')
                return
            if len(self.r2.sequence) == 0:
                self.errorDump('Sequence is empty, can not run forward model.')
                return
            forwardOutputStack.setCurrentIndex(0)
            self.forwardLogText.clear()
            QApplication.processEvents()
            # apply region for initial model

            x, phase0, zones, fixed = self.regionTable.getTable()
            regid = np.arange(len(x)) + 1 # region 0 doesn't exist
            pdebug('forwardBtnFunc(): with {:d} regions'.format(len(x)))
            self.r2.setStartingRes(dict(zip(regid, x)),
                                   dict(zip(regid, zones)),
                                   dict(zip(regid, fixed)),
                                   dict(zip(regid, phase0)))
            noise = float(self.noiseEdit.text())
            noiseIP = float(self.noiseEditIP.text())
            self.r2.forward(noise=noise, noiseIP=noiseIP, iplot=False, dump=forwardLogTextFunc)
            self.calcAspectRatio()
            self.mwFwdPseudo.plot(self.r2.surveys[0].showPseudo, aspect='auto')
            self.fwdContour.setVisible(True)
            self.tabs.setTabEnabled(4, True)
            self.tabs.setTabEnabled(5, True)
            self.tabs.setTabEnabled(6, True)
            if self.r2.typ[0] == 'c':
                self.mwFwdPseudoIP.plot(self.r2.surveys[0].showPseudoIP, aspect='auto')
        self.forwardBtn = QPushButton('Forward Modelling')
        self.forwardBtn.setAutoDefault(True)
        self.forwardBtn.clicked.connect(forwardBtnFunc)
        self.forwardBtn.setStyleSheet('background-color: green')

        self.mwFwdPseudo = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        self.mwFwdPseudoIP = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        self.mwFwdPseudoIP.setVisible(False)

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
                
        def fwdContourFunc(state):
            if state == Qt.Checked:
                contour = True
            else:
                contour = False
                
            self.mwFwdPseudo.setCallback(self.r2.surveys[0].showPseudo)
            self.mwFwdPseudo.replot(aspect='auto', contour=contour)
            if self.r2.typ[0] == 'c':
                self.mwFwdPseudoIP.setCallback(self.r2.surveys[0].showPseudoIP)
                self.mwFwdPseudoIP.replot(aspect='auto', contour=contour)
        
        self.fwdContour = QCheckBox('Contour')
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
        seqOptionLayout.addLayout(seqRowLayout)
        seqOptionLayout.addWidget(addRowBtn)
        seqLayout.addLayout(seqOptionLayout, 50)
        seqLayout.addWidget(arrayLabel, 50)
        
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
        forwardPseudoLayoutBottom.addWidget(self.mwFwdPseudo)
        forwardPseudoLayoutBottom.addWidget(self.mwFwdPseudoIP)
        self.mwFwdPseudoIP.hide()
        
        forwardPseudoLayout.addLayout(forwardPseudoLayoutBottom)

        forwardPseudos = QWidget()
        forwardPseudos.setLayout(forwardPseudoLayout)

        forwardOutputStack = QStackedLayout()
        forwardOutputStack.addWidget(self.forwardLogText)
        forwardOutputStack.addWidget(forwardPseudos)
        forwardOutputStack.setCurrentIndex(0)
                
        # general forward layout
        forwardLayout.addWidget(seqLabel)
        forwardLayout.addLayout(seqLayout)
        forwardLayout.addLayout(noiseLayout)

        fwdTopWidget = QWidget()
        fwdTopWidget.setObjectName('fwdTopWidget')
        fwdTopWidget.setStyleSheet("QWidget#fwdTopWidget {border:1px solid rgb(185,185,185)}")
        fwdTopWidget.setLayout(forwardLayout)
        
        fwdBottomWidget = QWidget()
        fwdBottomWidget.setObjectName('fwdBottomWidget')
        fwdBottomWidget.setStyleSheet("QWidget#fwdBottomWidget {border:1px solid rgb(185,185,185)}")
        fwdBottomWidget.setLayout(forwardOutputStack)
        
        fwdSplitter.addWidget(fwdTopWidget)
        fwdSplitter.addWidget(fwdBottomWidget)
        fwdSplitter.setCollapsible(1, False)
        
        fwdSlpitterLayout.addWidget(fwdSplitter)
        
        # instantiate the first row of the table
        # need to do this here as the layout needs to be integrated otherwise
        # it creates floating windows at the start
        seqRow.createRow()
        seqRow.showArg()
        
        tabForward.setLayout(fwdSlpitterLayout)


        #%% tab INVERSION SETTINGS
        tabInversionSettings = QTabWidget()
        self.tabs.addTab(tabInversionSettings, 'Inversion settings')
        self.tabs.setTabEnabled(4, False)

        # general tab
        generalSettings = QWidget()
        generalLayout = QHBoxLayout()
        invForm = QFormLayout()

        # advanced tab
        advancedSettings = QWidget()
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
            if self.r2.iForward is False and len(self.r2.surveys) >= 1:
                if 'magErr' in self.r2.surveys[0].df.columns:
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
                helpSection.setText('SORRY NOT IN HELP')
            else:
                helpSection.setHtml(r2help[arg])

        def showHelpAdv(arg): # for advanced tab
            if arg not in r2help:
                self.helpSection2.setText('SORRY NOT IN HELP')
            else:
                self.helpSection2.setHtml(r2help[arg])


        self.parallelLabel = QLabel('<a href="parallel">Parallel inversion</a>')
        self.parallelLabel.linkActivated.connect(showHelpAdv)
        self.parallelCheck = QCheckBox()
        advForm.addRow(self.parallelLabel, self.parallelCheck)

        self.modErrLabel = QLabel('<a href="modErr">Compute Modelling Error</a>')
        self.modErrLabel.linkActivated.connect(showHelpAdv)
        self.modErrCheck = QCheckBox()
        advForm.addRow(self.modErrLabel, self.modErrCheck)

        def notCroppingFunc(state):
            if state == Qt.Checked:
                self.iCropping = False
                if 'num_xz_poly' in self.r2.param:
                    self.num_xz_poly = self.r2.param['num_xz_poly'] # store value
                elif 'num_xy_poly' in self.r2.param:
                    self.num_xz_poly = self.r2.param['num_xy_poly'] # store value
            else:
                self.iCropping = True # default
                if ('num_xz_poly' in self.r2.param) and (self.num_xz_poly is not None):
                    self.r2.param['num_xz_poly'] = self.num_xz_poly # restore value
                elif ('num_xy_poly' in self.r2.param) and (self.num_xz_poly is not None):
                    self.r2.param['num_xy_poly'] = self.num_xz_poly # restore value
        self.notCroppingLabel = QLabel('<a href="notCropping">Do not crop the output</a>')
        self.notCroppingLabel.linkActivated.connect(showHelpAdv)
        self.notCropping = QCheckBox()
        self.notCropping.stateChanged.connect(notCroppingFunc)
        advForm.addRow(self.notCroppingLabel, self.notCropping)
                    
        cropBelowFmdLabel = QLabel('<a href="cropBelowFmd">Crop below mesh fine region</a>')
        cropBelowFmdLabel.linkActivated.connect(showHelpAdv)
        self.cropBelowFmd = QCheckBox()
        self.cropBelowFmd.setChecked(True)
        self.cropBelowFmd.stateChanged.connect(notCroppingFunc)
        advForm.addRow(cropBelowFmdLabel, self.cropBelowFmd)

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
                self.r2.param['checkTxSign'] = True
            else:
                self.r2.param['checkTxSign'] = False
        self.checkTxSignLabel = QLabel('<a href="txSign">Resistance polarity check</a>')
        self.checkTxSignLabel.linkActivated.connect(showHelpAdv)
        self.checkTxSign = QCheckBox()
        self.checkTxSign.stateChanged.connect(checkTxSignFunc)
        advForm.addRow(self.checkTxSignLabel, self.checkTxSign)

        def flux_typeFunc(index):
            if index == 0:
                self.r2.param['flux_type'] = 3
            else:
                self.r2.param['flux_type'] = 2
        self.flux_typeLabel = QLabel('<a href="flux_type">Flux Type</a>')
        self.flux_typeLabel.linkActivated.connect(showHelp)
        self.flux_type = QComboBox()
        self.flux_type.addItem('3D')
        self.flux_type.addItem('2D')
        self.flux_type.activated.connect(flux_typeFunc)
        advForm.addRow(self.flux_typeLabel, self.flux_type)

        def singular_typeFunc(state):
            if state == Qt.Checked:
                self.r2.param['singular_type'] = 1
            else:
                self.r2.param['singular_type'] = 0
        self.singular_typeLabel = QLabel('<a href="singular_type">Remove Singularity</a>')
        self.singular_typeLabel.linkActivated.connect(showHelp)
        self.singular_type = QCheckBox()
        self.singular_type.stateChanged.connect(singular_typeFunc)
        advForm.addRow(self.singular_typeLabel, self.singular_type)

        def res_matrixFunc(index):
            self.r2.param['res_matrix'] = index
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
            self.r2.param['scale'] = float(self.scale.text())
        self.scaleLabel = QLabel('<a href="scale"> Scale for triangular mesh</a>')
        self.scaleLabel.linkActivated.connect(showHelp)
        self.scaleLabel.setVisible(False)
        self.scale = QLineEdit()
        self.scale.setValidator(QDoubleValidator())
        self.scale.setText('1.0')
        self.scale.editingFinished.connect(scaleFunc)
        self.scale.setVisible(False)
        invForm.addRow(self.scaleLabel, self.scale)

        def patch_size_xFunc():
            self.r2.param['patch_size_x'] = int(self.patch_size_x.text())
        self.patch_size_xLabel = QLabel('<a href="patch">Patch size x<a/>:')
        self.patch_size_xLabel.linkActivated.connect(showHelpAdv)
        self.patch_size_x = QLineEdit()
        self.patch_size_x.setValidator(QIntValidator())
        self.patch_size_x.setText('1')
        self.patch_size_x.editingFinished.connect(patch_size_xFunc)
        advForm.addRow(self.patch_size_xLabel, self.patch_size_x)

        def patch_size_yFunc():
            self.r2.param['patch_size_y'] = int(self.patch_size_y.text())
        self.patch_size_yLabel = QLabel('<a href="patch">Patch size y<a/>:')
        self.patch_size_yLabel.linkActivated.connect(showHelpAdv)
        self.patch_size_y = QLineEdit()
        self.patch_size_y.setValidator(QIntValidator())
        self.patch_size_y.setText('1')
        self.patch_size_y.editingFinished.connect(patch_size_yFunc)
        advForm.addRow(self.patch_size_yLabel, self.patch_size_y)

        self.inv_typeVisible = []
        def inv_typeFunc(arg):
            index = int(self.inv_type.currentText()[-2])
            self.r2.param['inverse_type'] = index
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
            self.r2.param['target_decrease'] = float(self.target_decrease.text())
        self.target_decreaseLabel = QLabel('<a href="target_decrease">Target decrease</a>:')
        self.target_decreaseLabel.linkActivated.connect(showHelp)
        self.target_decrease = QLineEdit()
        self.target_decrease.setValidator(QDoubleValidator())
        self.target_decrease.setText('0')
        self.target_decrease.editingFinished.connect(target_decreaseFunc)
        invForm.addRow(self.target_decreaseLabel, self.target_decrease)


        def data_typeFunc(index):
            self.r2.param['data_type'] = index
        self.data_typeLabel = QLabel('<a href="data_type">Data type</a>:')
        self.data_typeLabel.linkActivated.connect(showHelp)
        self.data_type = QComboBox()
        self.data_type.addItem('Normal [0]')
        self.data_type.addItem('Logarithmic [1]')
        self.data_type.setCurrentIndex(1)
        self.data_type.activated.connect(data_typeFunc)
        invForm.addRow(self.data_typeLabel, self.data_type)

        def reg_modeFunc(index):
            self.r2.param['reg_mode'] = index
        self.reg_modeLabel = QLabel('<a href="reg_mode">Regularization mode</a>:')
        self.reg_modeLabel.linkActivated.connect(showHelp)
        self.reg_mode = QComboBox()
        self.reg_mode.addItem('Normal regularization [0]')
        self.reg_mode.addItem('Regularization from initial model [1]')
        self.reg_mode.addItem('Regularization from difference inversion [2]')
        self.reg_mode.activated.connect(reg_modeFunc)
        invForm.addRow(self.reg_modeLabel, self.reg_mode)

        def toleranceFunc():
            self.r2.param['tolerance'] = float(self.tolerance.text())
        self.toleranceLabel = QLabel('<a href="tolerance">Value for tolerance</a>:')
        self.toleranceLabel.linkActivated.connect(showHelp)
        self.tolerance = QLineEdit()
        self.tolerance.setValidator(QDoubleValidator())
        self.tolerance.setText('1.0')
        self.tolerance.editingFinished.connect(toleranceFunc)
        invForm.addRow(self.toleranceLabel, self.tolerance)

        def max_iterationsFunc():
            self.r2.param['max_iter'] = int(self.max_iterations.text())
        self.max_iterationsLabel = QLabel('<a href="max_iterations">Maximum number of iterations</a>:')
        self.max_iterationsLabel.linkActivated.connect(showHelp)
        self.max_iterations = QLineEdit()
        self.max_iterations.setValidator(QIntValidator())
        self.max_iterations.setText('10')
        self.max_iterations.editingFinished.connect(max_iterationsFunc)
        invForm.addRow(self.max_iterationsLabel, self.max_iterations)

        def error_modFunc(index):
            self.r2.param['error_mod'] = index
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
            self.r2.param['alpha_aniso'] = float(self.alpha_aniso.text())
        self.alpha_anisoLabel = QLabel('<a href="alpha_aniso">Value for <code>alpha_aniso</code></a>:')
        self.alpha_anisoLabel.linkActivated.connect(showHelpAdv)
        self.alpha_aniso = QLineEdit()
        self.alpha_aniso.setValidator(QDoubleValidator())
        self.alpha_aniso.setText('1.0')
        self.alpha_aniso.editingFinished.connect(alpha_anisoFunc)
        advForm.addRow(self.alpha_anisoLabel, self.alpha_aniso)

        def alpha_sFunc():
            self.r2.param['alpha_s'] = float(alpha_s.text())
        alpha_sLabel = QLabel('<a href="alpha_s"><code>alpha_s</code></a>:')
        alpha_sLabel.linkActivated.connect(showHelpAdv)
        alpha_s = QLineEdit()
        alpha_s.setValidator(QDoubleValidator())
        alpha_s.setText('1.0')
        alpha_s.editingFinished.connect(alpha_sFunc)
        advForm.addRow(alpha_sLabel, alpha_s)
        alpha_sLabel.setVisible(False)
        alpha_s.setVisible(False)

        def min_errorFunc():
            self.r2.param['min_error'] = float(self.min_error.text())
        self.min_errorLabel = QLabel('<a href="errorParam"><code>min_error</code></a>:')
        self.min_errorLabel.linkActivated.connect(showHelp)
        self.min_errorLabel.setVisible(False)
        self.min_error = QLineEdit()
        self.min_error.setText('0.01')
        self.min_error.editingFinished.connect(min_errorFunc)
        self.min_error.setVisible(False)
        invForm.addRow(self.min_errorLabel, self.min_error)

        def a_wgtFunc():
            self.r2.param['a_wgt'] = float(self.a_wgt.text())
        self.a_wgtLabel = QLabel('<a href="errorParam"><code>a_wgt</code></a>:')
        self.a_wgtLabel.linkActivated.connect(showHelp)
        self.a_wgt = QLineEdit()
        self.a_wgt.setValidator(QDoubleValidator())
        self.a_wgt.setText('0.01')
        self.a_wgt.editingFinished.connect(a_wgtFunc)
        invForm.addRow(self.a_wgtLabel, self.a_wgt)

        def b_wgtFunc():
            self.r2.param['b_wgt'] = float(self.b_wgt.text())
        self.b_wgtLabel = QLabel('<a href="errorParam"><code>b_wgt</code></a>:')
        self.b_wgtLabel.linkActivated.connect(showHelp)
        self.b_wgt = QLineEdit()
        self.b_wgt.setValidator(QDoubleValidator())
        self.b_wgt.setText('0.02')
        self.b_wgt.editingFinished.connect(b_wgtFunc)
        invForm.addRow(self.b_wgtLabel, self.b_wgt)

        def rho_minFunc():
            self.r2.param['rho_min'] = float(self.rho_min.text())
        self.rho_minLabel = QLabel('<a href="rho_max">Minimum apparent resistivity</a>:')
        self.rho_minLabel.linkActivated.connect(showHelp)
        self.rho_min = QLineEdit()
        self.rho_min.setValidator(QDoubleValidator())
        self.rho_min.setText('-1000')
        self.rho_min.editingFinished.connect(rho_minFunc)
        invForm.addRow(self.rho_minLabel, self.rho_min)

        def rho_maxFunc():
            self.r2.param['rho_max'] = float(self.rho_max.text())
        self.rho_maxLabel = QLabel('<a href="rho_max">Maximum apparent resistivity</a>:')
        self.rho_maxLabel.linkActivated.connect(showHelp)
        self.rho_max = QLineEdit()
        self.rho_max.setValidator(QDoubleValidator())
        self.rho_max.setText('1000')
        self.rho_max.editingFinished.connect(rho_maxFunc)
        invForm.addRow(self.rho_maxLabel, self.rho_max)
        

        generalLayout.addLayout(invForm)

        helpSection = QTextEdit('Help will be display here')
        helpSection.setReadOnly(True)
        helpSection.setText('Click on the labels and help will be displayed here')
        generalLayout.addWidget(helpSection)

        generalSettings.setLayout(generalLayout)
        tabInversionSettings.addTab(generalSettings, 'General')


        advancedLayout.addLayout(advForm)

        self.helpSection2 = QTextEdit('Help will be display here')
        self.helpSection2.setReadOnly(True)
        self.helpSection2.setText('Click on the labels and help will be displayed here')
        advancedLayout.addWidget(self.helpSection2)

        advancedSettings.setLayout(advancedLayout)
        tabInversionSettings.addTab(advancedSettings, 'Advanced')


        #%% tab 5 INVERSION
        tabInversion = QWidget()
        self.tabs.addTab(tabInversion, '&Inversion')
        self.tabs.setTabEnabled(5, False)
        
        invtabs = QTabWidget()
        logTab = QWidget()
        showTab = QWidget()
        invtabs.addTab(logTab, 'Log')
        invtabs.addTab(showTab, 'Results')
        invtabs.setTabEnabled(1, False)
                
        
        def frozeUI(val=True): # when inversion is running
            n = self.tabs.count()
            if val == True: # froze them
                self.tabState = [self.tabs.isTabEnabled(i) for i in range(n)]
                for i in range(n):
                    if i != 5:
                        self.tabs.setTabEnabled(i, False)
            else: # unfrozing
                for i in range(n):
                    self.tabs.setTabEnabled(i, self.tabState[i])
                    
            
        # ------------------------ log sub tab
        def invertBtnFunc():
            self.end = False # will be set to True if inversion successfull
            invtabs.setCurrentIndex(0) # log tab
            invtabs.setTabEnabled(1, False)
            
            # check if we kill or invert
            if self.invertBtn.text() == 'Invert':
                self.invertBtn.setStyleSheet('background-color: red')
                self.invertBtn.setText('Kill')
                frozeUI(True)
            else:
                print('killing...', end='')
                # self.loadingWidget('Killing in progress...')
                if self.r2.proc is not None:
                    self.r2.proc.kill()
                print('done')
                self.invertBtn.setStyleSheet('background-color: green')
                self.invertBtn.setText('Invert')
                # self.loadingWidget(exitflag=True)
                frozeUI(False)
                return
            
            # clear stuff and connect
            self.logText.setText('')
            self.mwRMS.clear()
            self.mwIter.clear()
            self.mwInv.clear()
            
            # create default mesh is not specified
            if self.r2.mesh is None:
                logTextFunc('Creating the mesh... ')
                if (self.r2.typ == 'R2') | (self.r2.typ == 'cR2'):
                    meshTrianFunc()
                else:
                    elec = self.elecTable.getTable()
                    topo = self.topoTable.getTable()
                    inan = ~np.isnan(topo[:,0])
                    if np.sum(~np.isnan(elec[:,0])) == 0:
                        self.errorDump('Please first import data or specify electrodes in the "Electrodes (XYZ/Topo)" tab.')
                        self.invertBtn.setStyleSheet('background-color: green')
                        self.invertBtn.setText('Invert')
                        frozeUI(False)
                        return
                    elif all(elec[:,1] == 0) & all(topo[inan,1] == 0):
                        self.errorDump('For 3D meshes, Y coordinates must be supplied for topo or elec at least.')
                        self.invertBtn.setStyleSheet('background-color: green')
                        self.invertBtn.setText('Invert')
                        frozeUI(False)
                        return
                    meshTetraFunc()
                logTextFunc('done!\n')

            # don't crop the mesh if that's what we've chosen
            if self.iCropping is True:
                if self.num_xz_poly is not None:
                    self.r2.param['num_xz_poly'] = self.num_xz_poly
                    self.r2.param['num_xy_poly'] = self.num_xz_poly
            else:
                self.r2.param['num_xz_poly'] = 0
                self.r2.param['num_xy_poly'] = 0
            
            # set initial model
            x, phase0, zones, fixed = self.regionTable.getTable()
            regid = np.arange(len(x)) + 1 # 1 is the background (no 0)
            self.r2.setStartingRes(dict(zip(regid, x)),
                                   dict(zip(regid, zones)),
                                   dict(zip(regid, fixed)),
                                   dict(zip(regid, phase0)))

            # invert
            # TODO run inversion in different thread
            self.r2.invert(iplot=False, dump=logTextFunc,
                           modErr=self.modErrCheck.isChecked(),
                           parallel=self.parallelCheck.isChecked(),
                           modelDOI=self.modelDOICheck.isChecked())

            # replace the log by the R2.out
            with open(os.path.join(self.r2.dirname, self.r2.typ + '.out'),'r') as f:
                text = f.read()
            self.logText.setText(text)
            self.r2.proc = None
            
            if any(self.r2.mesh.attr_cache['param'] == 0): # if fixed element are present, the mesh
            # will be sorted, meaning we need to replot it
                self.mwMesh.replot()
            
            # show results
            if self.end is True:
                # try:
                invtabs.setCurrentIndex(1) # show tab
                invtabs.setTabEnabled(1, True)
                if self.r2.typ[-1] == '2':
                    showStackedLayout.setCurrentIndex(0)
#                    maxDepth = self.r2.fmd if self.cropBelowFmd.isChecked() else None
                    self.mwInv.setCallback(partial(self.r2.showResults, cropMaxDepth=self.cropBelowFmd.isChecked()))# maxDepth=maxDepth))
                else:
                    showStackedLayout.setCurrentIndex(1)
                
                if self.r2.typ == 'R2' or self.r2.typ == 'R3t':
                    self.displayParams['attr'] = 'Resistivity(log10)'
                else:
                    self.displayParams['attr'] = 'Sigma_real(log10)'
                self.surveyCombo.clear()
                for m in self.r2.meshResults:
                    self.surveyCombo.addItem(m.mesh_title)
                index = 0
                if self.r2.iForward: # display the inverted not the initial
                    self.surveyCombo.setCurrentIndex(1)
                    index = 1
                    self.displayParams['index'] = 1
                self.attrCombo.clear()
                attrs = self.r2.meshResults[index].attr_cache.keys()
                c = 0
                ci = 0
                for i, attr in enumerate(attrs):
                    if attr not in ['param', 'region', 'zone']:
                        self.attrCombo.addItem(attr)
                        if attr == self.displayParams['attr']:
                            ci = c
                        c += 1
                self.attrCombo.setCurrentIndex(ci)
                self.vminEdit.setText('')
                self.vmaxEdit.setText('')
                self.doiCheck.setChecked(self.modelDOICheck.isChecked())
                self.contourCheck.setChecked(False)
                self.edgeCheck.setChecked(False)
                replotSection() # this plot the results
                try:
                    prepareInvError() # plot error graphs
                except Exception as e:
                    print('ERROR: ui: invertBtnFunc:', e)
                    self.errorDump(e)
                    pass
            self.invertBtn.setText('Invert')
            self.invertBtn.setStyleSheet('background-color: green')
            frozeUI(False)
            
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
                if a[0] == 'Iteration':
                    if self.typ[-1]=='t':
                        cropMaxDepth=False # this parameter doesnt make sense for 3D surveys 
                    else:
                        cropMaxDepth=self.cropBelowFmd.isChecked()
                    self.mwIter.plot(partial(self.r2.showIter, modelDOI=self.modelDOICheck.isChecked(), 
                                             cropMaxDepth=cropMaxDepth), aspect=self.plotAspect)
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
            ax.set_xticklabels([],[])
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
        self.vtkWidget = pv.QtInteractor(self.frame)
        vlayout.addWidget(self.vtkWidget.interactor)
        self.frame.setLayout(vlayout)
        self.mSlice = None
        self.mMesh = None
        
        self.displayParams = {'index':0,'edge_color':'none',
                            'sens':True, 'attr':'Resistivity(Ohm-m)',
                            'contour':False, 'vmin':None, 'vmax':None,
                            'cmap':'viridis', 'sensPrc':0.5,
                            'doi':self.modelDOICheck.isChecked(),
                            'doiSens':False,
                            'pvslices':([],[],[]), 'pvthreshold':None,
                            'pvgrid':False, 'pvcontour':[]}
        

        def replotSection(): # main plotting function
            pdebug('replotSection() with', self.displayParams)
            index = self.displayParams['index']
            edge_color = self.displayParams['edge_color']
            sens = self.displayParams['sens']
            attr = self.displayParams['attr']
            contour = self.displayParams['contour']
            vmin = self.displayParams['vmin']
            vmax = self.displayParams['vmax']
            cmap = self.displayParams['cmap']
            sensPrc = self.displayParams['sensPrc']
            doi = self.displayParams['doi']
            doiSens = self.displayParams['doiSens']
            pvslices = self.displayParams['pvslices']
            pvthreshold = self.displayParams['pvthreshold']
            pvgrid = self.displayParams['pvgrid']
            pvcontour = self.displayParams['pvcontour']
            if self.r2.typ[-1] == '2':
                self.mwInv.replot(threed=False, aspect=self.plotAspect,
                                  index=index, edge_color=edge_color,
                                  contour=contour, sens=sens, attr=attr,
                                  vmin=vmin, vmax=vmax, color_map=cmap, 
                                  sensPrc=sensPrc, doi=doi, doiSens=doiSens)
            else:
                # mwInvResult3D.replot(threed=True, index=index, attr=attr,
                                      # vmin=vmin, vmax=vmax, color_map=cmap)
                
                # if self.r2.iTimeLapse and index == 0:
                    # fname = os.path.join(self.r2.dirname, 'f001.vtk')
                # else:
                    # fname = os.path.join(self.r2.dirname, 'f{:03d}.vtk'.format(index+1))
                # m = pv.read(fname)
                
                self.vtkWidget.clear()
                self.vtkWidget.clear_plane_widgets()
                self.r2.meshResults[index].show3D(ax=self.vtkWidget, attr=attr,
                                                  edge_color=edge_color, vmin=vmin,
                                                  vmax=vmax, color_map=cmap,
                                                  background_color=(0.8,0.8,0.8),
                                                  pvslices=pvslices,
                                                  pvthreshold=pvthreshold,
                                                  pvgrid=pvgrid,
                                                  pvcontour=pvcontour)
                
                
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
            self.displayParams['index'] = index
            attrs = list(self.r2.meshResults[index].attr_cache.keys())
            attr0 = str(self.attrCombo.currentText())
            ci = 0
            c = -1
            found = False
            for i, a in enumerate(attrs): # find same attribute or plot first one
                if a not in ['param', 'region', 'zone']:
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
                if attr not in ['param', 'region', 'zone']:
                    self.attrCombo.addItem(attr)
            self.attrCombo.setCurrentIndex(ci)
            replotSection()
        self.surveyCombo = QComboBox()
        self.surveyCombo.activated.connect(surveyComboFunc)
        self.surveyCombo.setToolTip('Change survey or see initial model.')
        
        # change attribute
        def attrComboFunc(index):
            resetAttributeSpecificSettings()
            self.displayParams['attr'] = str(self.attrCombo.currentText())
            if self.displayParams['attr'] == 'Sigma_imag(log10)':
                sigma_imag_vals = self.r2.meshResults[self.displayParams['index']].attr_cache['Sigma_imag(log10)']
                if any(val == 0 for val in sigma_imag_vals):
                    if all(val == 0 for val in sigma_imag_vals):
                        pass
                    else:
                        self.contourCheck.setChecked(True)
                        self.infoDump('Contouring data by default!')
                        return # replotting triggers by edgeCheckFunc()
            replotSection()
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
            if (self.contourCheck.isChecked() is True) or (self.r2.typ[-1] != '2'):
                replotSection()
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
            replotSection()
        self.doiCheck = QCheckBox('DOI')
        self.doiCheck.setVisible(False)
        self.doiCheck.stateChanged.connect(doiCheckFunc)
        self.doiCheck.setToolTip('Depth of Investigation (DOI) based on Oldenburg and Li 1999.')
        
        def doiSensCheckFunc(status):
            if status == Qt.Checked:
                self.displayParams['doiSens'] = True
            else:
                self.displayParams['doiSens'] = False
            replotSection()
        self.doiSensCheck = QCheckBox('DOI estimate')
        self.doiSensCheck.stateChanged.connect(doiSensCheckFunc)
        self.doiSensCheck.setToolTip('Depth of Investigation (DOI) estimated based on sensitivity.\nSee advanced inversion settings for DOI with the Oldengburg and Li method.')
        
        self.cmapComboLabel = QLabel('Colormap')
        cmaps = ['viridis','plasma','seismic','cividis','winter','autumn','rainbow','jet','Spectral']
        def cmapComboFunc(index):
            self.displayParams['cmap'] = cmaps[index]
            replotSection()
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
            replotSection()
        self.contourCheck = QCheckBox('Contour')
        self.contourCheck.stateChanged.connect(contourCheckFunc)
        self.contourCheck.setToolTip('Grid and contour the data.')
      
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
            replotSection()
        self.sensWidget = QWidget()
        sensLayout = QVBoxLayout()
        sensLayout.setContentsMargins(0,9,0,9)
        sensLayout.setSpacing(2)
        self.sensSlider = QSlider(Qt.Horizontal)
        self.sensSlider.setFixedWidth(110)
        self.sensSlider.setMinimum(0)
        self.sensSlider.setMaximum(11)
        self.sensSlider.setValue(5)
        sensLabel = QLabel('Sensitivity overlay')
        sensLabel.setAlignment(Qt.AlignCenter )
        sensLayout.addWidget(sensLabel)
        self.sensSlider.setToolTip('Normalized sensivity threshold')
        self.sensSlider.valueChanged.connect(sensSliderFunc)
        sensLayout.addWidget(self.sensSlider)
        self.sensWidget.setLayout(sensLayout)
        
        def showEdges(status):
            if status == Qt.Checked:
                self.displayParams['edge_color'] = 'k'
            else:
                self.displayParams['edge_color'] = 'none'
            replotSection()
        self.edgeCheck= QCheckBox('Edges')
        self.edgeCheck.setChecked(False)
        self.edgeCheck.setToolTip('Show edges of each mesh cell.')
        self.edgeCheck.stateChanged.connect(showEdges)

        def saveBtnFunc():
            fdir = QFileDialog.getExistingDirectory(self.tabImportingData, 'Choose the directory to export graphs and .vtk', directory=self.datadir)
            if fdir != '':
                if self.r2.typ[-1] == '2':
                    edge_color = self.displayParams['edge_color']
                    sens = self.displayParams['sens']
                    sensPrc = self.displayParams['sensPrc']
                    attr = self.displayParams['attr']
                    contour = self.displayParams['contour']
                    vmin = self.displayParams['vmin']
                    vmax = self.displayParams['vmax']
                    cmap = self.displayParams['cmap']
                    self.r2.saveInvPlots(outputdir=fdir, edge_color=edge_color,
                                       contour=contour, sens=sens, attr=attr,
                                       vmin=vmin, vmax=vmax, color_map=cmap,
                                       sensPrc=sensPrc)
                self.r2.saveVtks(fdir)
            self.r2.saveData(fdir)
            self.infoDump('All graphs saved successfully in the working directory.')
        self.saveBtn = QPushButton('Save Data')
        self.saveBtn.clicked.connect(saveBtnFunc)
        self.saveBtn.setToolTip('Save current graph to the working directory.')
        
        
        # 3D specific options for pyvista
        pvthreshLabel = QLabel('Threshold:')
        pvthreshLabel.setToolTip('Value which to keep the cells.')
        self.pvthreshMin = QLineEdit('')
        self.pvthreshMin.setPlaceholderText('Min')
        self.pvthreshMin.setValidator(QDoubleValidator())
        self.pvthreshMin.setToolTip('Minimal value above which to keep the cells.')
        
        # pvthreshMaxLabel = QLabel('Max Threshold:')
        self.pvthreshMax = QLineEdit('')
        self.pvthreshMax.setPlaceholderText('Max')
        self.pvthreshMax.setValidator(QDoubleValidator())
        self.pvthreshMax.setToolTip('Maximal value below which to keep the cells.')
        
        pvslicesLabel = QLabel('Axis slices:')
        pvslicesLabel.setToolTip('Slice the mesh normal to X, Y and/or Z axis. '
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
                
        pvcontourLabel = QLabel('Isosurfaces:')
        self.pvcontour = QLineEdit('')
        self.pvcontour.setToolTip('Values of isosurfaces (comma separated).')
        
        def pvapplyBtnFunc():
            threshMin = float(self.pvthreshMin.text()) if self.pvthreshMin.text() != '' else None
            threshMax = float(self.pvthreshMax.text()) if self.pvthreshMax.text() != '' else None
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
            replotSection()
            
        def pvgridCheckFunc(state):
            self.displayParams['pvgrid'] = self.pvgridCheck.isChecked()
            replotSection()
        self.pvgridCheck = QCheckBox('Grid')
        self.pvgridCheck.stateChanged.connect(pvgridCheckFunc)
            
        self.pvapplyBtn = QPushButton('Apply 3D')
        self.pvapplyBtn.setAutoDefault(True)
        self.pvapplyBtn.clicked.connect(pvapplyBtnFunc)
        self.pvapplyBtn.setToolTip('Apply 3D options')
        
        def screenshotBtnFunc():
            fname, _ = QFileDialog.getSaveFileName(self, 'Open File', self.datadir,
                                                   'PNG (*.png);;TIFF (*.tif)')
            if fname != '':
                self.vtkWidget.screenshot(fname, transparent_background=True)
        self.screenshotBtn = QPushButton('Save screenshot')
        self.screenshotBtn.setAutoDefault(True)
        self.screenshotBtn.clicked.connect(screenshotBtnFunc)
            
        # opt3d = [pvthreshMinLabel, self.pvthreshMin,
        #          pvthreshMaxLabel, self.pvthreshMax,
        #          pvxslicesLabel, self.pvxslices,
        #          pvyslicesLabel, self.pvyslices,
        #          pvzslicesLabel, self.pvzslices,
        #          pvcontourLabel, self.pvcontour,
        #          self.pvapplyBtn,
        #          self.pvgridCheck,
        #          self.screenshotBtn]
        
        opt3d = [pvthreshLabel, self.pvthreshMin,
                 self.pvthreshMax, pvslicesLabel, 
                 self.pvxslices, self.pvyslices,
                 self.pvzslices, pvcontourLabel, 
                 self.pvcontour, self.pvapplyBtn,
                 self.pvgridCheck, self.screenshotBtn]
        
        def show3DInvOptions(a):
            [o.setVisible(a) for o in opt3d]
        show3DInvOptions(False)
        
        
        # layout
        invLayout = QVBoxLayout()
        logLayout = QHBoxLayout()
        showLayout = QVBoxLayout()
        
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
        optsLayout.addWidget(self.vminEdit, 10)
        optsLayout.addWidget(self.vmaxLabel)
        optsLayout.addWidget(self.vmaxEdit, 10)
        optsLayout.addWidget(self.vMinMaxBtn)
        optsLayout.addWidget(self.doiCheck)
        optsLayout.addWidget(self.doiSensCheck)
        optsLayout.addWidget(self.cmapCombo)
        optsLayout.addWidget(self.contourCheck)
        optsLayout.addWidget(self.sensWidget)
        optsLayout.addWidget(self.edgeCheck)
        optsLayout.addWidget(self.saveBtn)
        
        opt3dLayout = QHBoxLayout()
        for o in opt3d:
            opt3dLayout.addWidget(o)
        
        showLayout.addLayout(optsLayout)
        showLayout.addLayout(opt3dLayout)
        showStackedLayout = QStackedLayout()
        showStackedLayout.addWidget(self.mwInv)
        showStackedLayout.addWidget(self.frame)
        showLayout.addLayout(showStackedLayout)
        
        logTab.setLayout(logLayout)
        showTab.setLayout(showLayout)
        
        invLayout.addWidget(invtabs)
        tabInversion.setLayout(invLayout)
    
    

        #%% tab 6 POSTPROCESSING
        tabPostProcessing = QWidget()
        self.tabs.addTab(tabPostProcessing, 'Post-processing')
        self.tabs.setTabEnabled(6,False)
        
        self.errorGraphs = QTabWidget()
        
        self.invPseudoErrLabel = QLabel('Select datapoints on the pseudo section to remove and reinvert the data.')

        def prepareInvError():
            names = [s.name for s in self.r2.surveys]
            if len(names) > 1:
                self.invErrorComboLabel.show()
                self.invErrorCombo.show()
            else:
                self.invErrorComboLabel.hide()
                self.invErrorCombo.hide()
                
            if self.r2.iTimeLapse:
                names[0] = names[0] + ' (Ref)'
                
            # self.invErrorCombo.disconnect()
            self.invErrorCombo.clear()
            for name in names:
                self.invErrorCombo.addItem(name)
            self.invErrorCombo.activated.connect(invErrorComboFunc)
            invErrorComboFunc(0)
        
        def invErrorComboFunc(index):
            try:
                plotInvError2(index)
                if self.iBorehole is False:
                    plotInvError(index)
                    self.invErrorIndex = index
            except Exception as e:
                print('Could not print error: ', e)
        self.invErrorCombo = QComboBox()
        self.invErrorCombo.setMinimumWidth(250)
        self.invErrorCombo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.invErrorCombo.activated.connect(invErrorComboFunc)
        self.invErrorCombo.hide()
        self.invErrorComboLabel = QLabel('Choose a dataset to plot the error:')
        self.invErrorComboLabel.hide()
        
        def plotInvError(index=0):
            pdebug('plotInvError()')
            self.mwInvError.setCallback(self.r2.showPseudoInvError)
            self.mwInvError.replot(index=index)

        self.mwInvError = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        
        def invErrorFiltFunc():
            try:
                i2remove = self.r2.surveys[self.invErrorIndex].iselect
                if not all(self.r2.surveys[self.invErrorIndex].df['irecip'].values == 0): 
                    # as the selection is only done on normal dataset, reciprocal pair should be removed too
                    recipPaires = self.r2.surveys[self.invErrorIndex].df[i2remove]['irecip'].values*-1
                    ie = np.isin(self.r2.surveys[self.invErrorIndex].df['irecip'].values, recipPaires[recipPaires != 0]) 
                    i2remove = i2remove + ie
                self.r2.surveys[self.invErrorIndex].filterData(~i2remove)
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
        self.invErrorReinvertBtn.setStyleSheet('background-color: green')
      

        def plotInvError2(index=0):
            pdebug('plotInvError2()')
            self.mwInvError2.setCallback(self.r2.showInvError)
            self.mwInvError2.replot(index=index)
        self.mwInvError2 = MatplotlibWidget(navi=True, aspect='auto', itight=True)
        invErrorLabel = QLabel('All errors should be between +/- 3% (Binley at al. 1995). '
                               'If it\'s not the case try to fit an error model or '
                               'manually change the a_wgt and b_wgt in inversion settings.')


        # layout
        postProcessingLayout = QVBoxLayout()
        
        topInvErrorLayout = QHBoxLayout()
        topInvErrorLayout.addWidget(self.invErrorComboLabel, 1)
        topInvErrorLayout.addWidget(self.invErrorCombo, 0)
        postProcessingLayout.addLayout(topInvErrorLayout)
        postProcessingLayout.addWidget(self.errorGraphs)
        
        invError = QWidget()
        self.errorGraphs.addTab(invError, 'Pseudo Section of Inversion Errors')

        invErrorLayout = QVBoxLayout()
        invErrorLayout.setAlignment(Qt.AlignTop)
        
        invErrorTopLayout = QHBoxLayout()
        invErrorTopLayoutL = QHBoxLayout()
        invErrorTopLayoutL.addWidget(self.invPseudoErrLabel)
        invErrorTopLayout.addLayout(invErrorTopLayoutL)
        
        invErrorTopLayoutR = QHBoxLayout()
        invErrorTopLayoutR.setAlignment(Qt.AlignRight)
        invErrorTopLayoutR.addWidget(self.invErrorFilterBtn)
        invErrorTopLayoutR.addWidget(self.invErrorReinvertBtn)
        invErrorTopLayout.addLayout(invErrorTopLayoutR)
        
        invErrorLayout.addLayout(invErrorTopLayout, 0)
        
        invErrorLayout.addWidget(self.mwInvError, Qt.AlignCenter)
        invError.setLayout(invErrorLayout)

        invError2 = QWidget()
        self.errorGraphs.addTab(invError2, 'Normalised Inversion Errors')
        invErrorLayout2 = QVBoxLayout()
        invErrorLayout2Plot = QVBoxLayout()

        invErrorLayout2Plot.addWidget(self.mwInvError2, Qt.AlignCenter)
        invErrorLayout2.addLayout(invErrorLayout2Plot, 1)
        invErrorLayout2.addWidget(invErrorLabel)
        invError2.setLayout(invErrorLayout2)
        
        tabPostProcessing.setLayout(postProcessingLayout)



        #%% Help tab
        tabHelp = QTabWidget()
        self.tabs.addTab(tabHelp, 'Help')

        helpLayout = QVBoxLayout()
        helpLayout.setAlignment(Qt.AlignTop)
        helpText = QTextBrowser()
        helpText.setReadOnly(True)
        helpText.setOpenExternalLinks(True)
        helpText.setText('''
           <h1>General help</h1>\
           <p>Below are simple instructions to guide you to through the software.</p>
           <ul>
           <p><li>In the "<b>Importing</b>" tab:
           <ul>
           <li>Select if you want a 2D/3D survey, an inverse/forward solution and check if you have borehole/timelapse/batch data.</li>
           <li>Modify the default working directory if you want to keep the outputed files afterwards.</li>
           <li>Select the file type. You can choose "Custom" if you file type is not available and you will be redirected to the custom parser tab.</li>
           <ul><li>Note: Syscal files must be exported as 'Spreadsheet' files with .csv format (comma separators) from Prosys.</li>
           <li>Note: Res2DInv files are mostly supported, but it is recommended to change them in "General Array" format if your file is not recognized.</ul></li>
           <ul>
           <li>If your survey has topography, you can import it in the "Electrodes(XZY/Topo)" tab.</li>
           <ul><li><i>Pole-dipole arrays</i>: the remote electrode's X location must be exactly at 99999 or -99999 m.</li>
           <li><i>Pole-pole arrays</i>: first remote electrode's X location must be exactly at 99999 and second one at exactly -99999 m (or vice versa).</ul></li>
           <li>Then one can choose to directly invert with all the default settings or go through the other tabs on the rights.</li>
           </ul></li></p>
           <p><li>In the "<b>Pre-processing</b>" tab:
           <ul>
           <li>The first tab offers manual filtering option based on reciprocal measurements in the dataset (if any).</li>
           <li>The "Phase Filtering" tab is only enable for IP data and allows precise filtering of IP data (range filtering, removing of nested measuremetns, DCA, ...).</li>
           <li>The "Resistance Error Model" tab allows to fit a power-law or linear error model to resistance data.</li>
           <li>The "Phase Error Model" tab allows to fit a power-law or parabolic error model to phase data.</li>
           </ul></li></p>
           <p><li>In the "<b>Mesh</b>" tab you can create a quadrilateral or triangular mesh (2D) or a tetrahedral mesh (3D). For 2D mesh you can specify different\
           region of given resistivity/phase and if they need to be fixed or not during inversion. For forward modelling this mesh serves as the initial model.</li></p>
           <p><li>In the "<b>Forward model</b>" tab (only available in forward mode) you can design your sequence and add noise. The resulting synthetic measurements will be\
           automatically added to as an actual survey in ResIPy and can be inverted directly.</li></p>
           <p><li>In the "<b>Inversion Settings</b>" tab, you can modify all settings for the inversion. Help is available by clicking on the label of each item. The help\
           generally refers to the document present in the R2/cR3/R3t/cR3t respective manuals.</li></p>
           <p><li>In the "<b>Inversion</b>" tab, you can invert your survey and see the output in real time. if you have selected parallel inversion in "Inversion Settings">"Advanced",\
           then nothing will be printed out until the inversion finished. When the inversion finished you will be able to see the inverted section, open it with Paraview (mainly for 3D)\
           and save the outputed .vtk file and graphs using the "Save Graphs" button.</li>
           <ul><li>Plot aspect ratio can be changed by dragging  left handle to right or left and top handle (above plot options) up and down.</li></ul></p>
           <p><li>The "<b>Post-processing</b>" tab displays the errors from the invesrion. It helps to assess the quality of the inversion.</li>
           </ul></p>
           <p><b>Figure options</b> (available under each figure, the button next to save icon):</p>
           <ul>
           <li>Select the axes to edit (usually "<i>Distance [m] - Elevation [m]</i>):
           <ul><li>The axes limits/labels and plot title  can be changed in the "Axes" tab.</li>
           <li>The marker size and line width can be changed in the "Curves" tab.</li></ul>
           </ul>
           <p>More help for ResIPy: <a href="https://hkex.gitlab.io/pyr2/">https://hkex.gitlab.io/pyr2</a></p>
           <p>Read more on 2D resistivity inversion: <a href="http://www.es.lancs.ac.uk/people/amb/Freeware/R2/R2_readme.pdf">R2_readme.pdf</a></p>
           <p>Read more on 3D resistivity inversion: <a href="http://www.es.lancs.ac.uk/people/amb/Freeware/R3t/R3t_readme.pdf">R3t_readme.pdf</a></p>
           <p>Read more on 2D complex resistivity (IP) inversion: <a href="http://www.es.lancs.ac.uk/people/amb/Freeware/cR2/cR2_readme.pdf">cR2_readme.pdf</a></p>
        ''')
        helpLayout.addWidget(helpText)
        tabHelp.setLayout(helpLayout)


        #%% About tab

        tabAbout = QTabWidget()
        self.tabs.addTab(tabAbout, 'About')

        infoLayout = QVBoxLayout()
        aboutText = QLabel()
        aboutText.setText('''<h1>About ResIPy </h1> \
                          <p><b>Version: %s</b></p> \
                          <p><i>ResIPy is a free and open source software for inversion and modeling of geoelectrical data (Resistivity and IP)</i></p> \
                          <p>If you encouter any issues or would like to submit a feature request, please raise an issue on our gitlab repository at:</p> \
                          <p><a href="https://gitlab.com/hkex/pyr2/issues">https://gitlab.com/hkex/pyr2/issues</a></p> \
                          <p>ResIPy uses 
                              <a href="http://www.es.lancs.ac.uk/people/amb/Freeware/R2/R2.htm">R2</a>,
                              <a href="http://www.es.lancs.ac.uk/people/amb/Freeware/cR2/cR2.htm">cR2</a>,
                              <a href="http://www.es.lancs.ac.uk/people/amb/Freeware/R3t/R3t.htm">R3t</a> and 
                              <a href="http://www.es.lancs.ac.uk/people/amb/Freeware/cR3t/cR3t.htm">cR3t</a> developed by Andrew Binley</p> \
                          <p>For generation of triangular mesh, ResIPy uses software 
                              <a href="http://gmsh.info/">Gmsh</a></p>\
                          <p>Python packages used: 
                              <a href="https://numpy.org/">numpy</a>, 
                              <a href="https://pandas.pydata.org/">pandas</a>,
                              <a href="https://matplotlib.org/">matplotlib</a>,
                              <a href="https://scipy.org/index.html">scipy</a>,
                              <a href="https://pypi.org/project/PyQt5/">PyQt5</a>.
                          </p>
<p><strong>ResIPy's core developers: Guillaume Blanchy, Sina Saneiyan, Jimmy Boyd and Paul McLachlan.<strong></p>
<p>Contributors: Pedro Concha, Michael Tso</p>
<p><b><a href="https://www.researchgate.net/project/pyR2-GUI-for-R2-family-codes">Visit our ResearchGate page!</a></b></p>
<p><b>Citing ResIPy</b>:<br>Blanchy G., Saneiyan S., Boyd J., McLachlan P. and Binley A. 2020.<br>“ResIPy, an Intuitive Open Source Software for Complex Geoelectrical Inversion/Modeling.”<br>Computers & Geosciences, February, 104423. <a href="https://doi.org/10.1016/j.cageo.2020.104423">https://doi.org/10.1016/j.cageo.2020.104423</a>.</p>

'''%ResIPy_version)
#        aboutText.setText('''<h1>About ResIPy</h1> \
#                          <p><b>Version: %s</b></p> \
#                          <p><i>ResIPy is a free and open source software for inversion of geoelectrical data (Resistivity and IP)</i></p> \
#                          <p>If you encouter any issues or would like to submit a feature request, please raise an issue on our gitlab repository at:</p> \
#                          <p><a href="https://gitlab.com/hkex/pyr2/issues">https://gitlab.com/hkex/pyr2/issues</a></p> \
#                          <p>pyR2 uses R2 and cR2 codes developed by Andrew Binley:</p> \
#                          <p><a href="http://www.es.lancs.ac.uk/people/amb/Freeware/R2/R2.htm">http://www.es.lancs.ac.uk/people/amb/Freeware/R2/R2.htm</a></p> \
#                          <p>For generation of triangular mesh, pyR2 uses "Gmsh" software:</p> \
#                          <p><a href="http://gmsh.info/">http://gmsh.info/</a></p>\
#                          <p>Python packages used: scipy, numpy, pandas, matplotlib.
#<ul>
#<li>Jones E, Oliphant E, Peterson P, <em>et al.</em>
#<strong>SciPy: Open Source Scientific Tools for Python</strong>, 2001-,
#<a class="reference external" href="http://www.scipy.org/">http://www.scipy.org/</a> [Online; accessed 2018-10-02].
#</li>
#<li>
# Wes McKinney.
#<strong>Data Structures for Statistical Computing in Python</strong>,
#Proceedings of the 9th Python in Science Conference, 51-56 (2010)
#(<a class="reference external" href="http://conference.scipy.org/proceedings/scipy2010/mckinney.html">publisher link</a>)
#</li>
#<li>
#John D. Hunter.
#<strong>Matplotlib: A 2D Graphics Environment</strong>,
#Computing in Science &amp; Engineering, <strong>9</strong>, 90-95 (2007),
#<a class="reference external" href="https://doi.org/10.1109/MCSE.2007.55">DOI:10.1109/MCSE.2007.55</a>
#</li>
#<li>Travis E, Oliphant. <strong>A guide to NumPy</strong>,
#USA: Trelgol Publishing, (2006).
#</li>
#</ul>
#</p>
        aboutText.setOpenExternalLinks(True)
        aboutText.setWordWrap(True)
        aboutText.setAlignment(Qt.AlignTop | Qt.AlignHCenter)
        infoLayout.addWidget(aboutText)

        tabAbout.setLayout(infoLayout)

        #%% test tab
        # tabTest = QTabWidget()
        # self.tabs.addTab(tabTest, 'TEST')
        # self.tabs.setCurrentIndex(9)
        
        # self.fframe = QFrame()
        # vlayout = QVBoxLayout()
        # self.vtk_widget = pv.QtInteractor(self.fframe)
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
        self.errorLabel.setText('<i style="color:black">['+timeStamp+']: </i>')
        
    def errorDump(self, text, flag=1):
        text = str(text)
        timeStamp = time.strftime('%H:%M:%S')
        if flag == 1: # error in red
            col = 'red'
            pdebug('errorDump:', text)
        else:
            col = 'black'
            pdebug('infoDump:', text)
        self.errorLabel.setText('<i style="color:'+col+'">['+timeStamp+']: '+text+'</i>')
        self.timer.timeout.connect(partial(self.timeOut, timeStamp))
        self.timer.start(10000) # 10 secs - making sure the error/message doen't stick forever!

    def infoDump(self, text):
        self.errorDump(text, flag=0)
        
    def loadingWidget(self, msgtxt='', exitflag=False):
        '''Shows a dialog to indicate ResIPy is working (e.g., loading large dataset, etc.)'''
        if exitflag == False:
            self.loadingDialogTxtWidget.setText(msgtxt)
            self.loadingDialog.show()
            QApplication.processEvents()
        
        else:
            self.loadingDialog.accept()  
            
    def calcAspectRatio(self): # for calculating aspect ratio of long surveys
        self.r2.computeFineMeshDepth()
        surLength = np.abs(self.r2.param['xz_poly_table'][0,0] - self.r2.param['xz_poly_table'][1,0])
        surDepth = np.abs(self.r2.param['xz_poly_table'][-1,1] - self.r2.param['xz_poly_table'][-2,1])
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
        if attr not in self.r2.surveys[index].df.columns:
            self.errorDump('Attribute {:s} not found for survey.'.format(attrText))
            return
        try:
            self.mwManualFiltering.setCallback(self.r2.filterManual)
            self.mwManualFiltering.replot(index=index, attr=attr)
        except ValueError as e:
            self.errorDump(e)
            self.mwManualFiltering.clear()

    def errHist(self, index=0):
        pdebug('errHist()')
        if all(self.r2.surveys[index].df['irecip'].values == 0) is False:
            if self.iBatch or self.iTimeLapse:
                self.mwRecipError.setCallback(self.r2.showErrorDist)
                self.mwRecipError.replot(index=index)
            else: 
                self.mwRecipError.plot(self.r2.showErrorDist)
        else:
            pass
    
    def plotError(self, index=0):
        if len(self.r2.surveys) == 0:
            return
        self.mwFitError.setCallback(self.r2.showError)
        self.mwFitError.replot(index=index)
        self.r2.err = False
        
    def topoInterpBtnFunc(self):
        elec = self.elecTable.getTable()
        topo = self.topoTable.getTable()
        inan = ~np.isnan(elec[:,2])
        inan2 = ~np.isnan(topo[:,2])
        points = np.r_[elec[inan,:2], topo[inan2,:2]]
        values = np.r_[elec[inan,2], topo[inan2,2]]
        xi = elec[~inan,:2]
        if self.r2.typ[-1] == '2':
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
        self.elecTable.setTable(elec)
        self.infoDump('Interpolation successful.')
        
    def updateElec(self):
        try:
            elec = self.elecTable.getTable()
            buried = self.elecTable.getBuried()
            if self.tempElec is None or np.sum(elec-self.tempElec) != 0:
                self.tempElec = elec
                self.r2.setElec(elec)
                self.r2.buried = buried
                self.r2.mesh = None
                self.mwMesh.clear()
                if len(self.r2.surveys) > 0:
                    self.plotPseudo()
                    self.plotManualFiltering()
                    if self.r2.typ[0] == 'c':
                        self.plotPseudoIP()
        except Exception as e:
            self.errorDump('Error updating pseudosection: ' + e)
    
    def plotPseudo(self):
        pdebug('plotPseudo()')
        self.mwPseudo.setCallback(self.r2.showPseudo)
        if (self.r2.typ == 'R3t') | (self.r2.typ == 'cR3t'):
            self.mwPseudo.replot(aspect='auto', **self.pParams)
        else:
            self.mwPseudo.replot(aspect='auto', **self.pParams)
    
    def plotPseudoIP(self):
        pdebug('plotPseudoIP()')
        self.mwPseudoIP.setCallback(self.r2.showPseudoIP)
        self.mwPseudoIP.replot(aspect='auto', **self.pParamsIP)

    def activateTabs(self, val=True):
        if self.iForward is False:
            self.tabs.setTabEnabled(1,val)
            self.tabs.setTabEnabled(2,val)
            self.tabs.setTabEnabled(4,val)
            self.tabs.setTabEnabled(5,val)
            self.tabs.setTabEnabled(6,val)
            try:
                if self.m3DRadio.isChecked():
                    if all(self.r2.surveys[0].df['irecip'].values == 0):
                        self.tabs.setTabEnabled(1, False)
                        self.tabPreProcessing.setTabEnabled(0, False)
                    else:
                        self.tabs.setTabEnabled(1, True)
                    if self.ipCheck.checkState() == Qt.Checked:
                        self.tabs.setTabEnabled(1, True)
            except:
                pass
        else:
            self.tabs.setTabEnabled(2,val)
        self.meshTrianGroup.setFocus() # needed after tab activation as default is the fmdBox (QLineEdit) for some strange reasons...

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
            
    def importFile(self, fname): # TODO to test the UI automatically this needs to be a methods
        pdebug('importFile:', fname)
        if len(self.r2.surveys) > 0:
            self.r2.surveys = []
        try:                
            self.loadingWidget('Loading data, please wait...', False) # for large datasets
            self.ipCheck.setEnabled(True)
            self.psContourCheck.setEnabled(True)
            self.fname = fname
            self.importDataBtn.setText(os.path.basename(self.fname) + ' (Press to change)')
            if float(self.spacingEdit.text()) == -1:
                spacing = None
            else:
                spacing = float(self.spacingEdit.text())
            try:
                self.r2.createSurvey(self.fname, ftype=self.ftype, spacing=spacing,
                                     parser=self.parser)
                # self.calcAspectRatio()
                if 'magErr' in self.r2.surveys[0].df.columns:
                    self.a_wgt.setText('0.0')
                    self.b_wgt.setText('0.0')
            except:
                self.errorDump('File is not recognized.')
                pass
            pdebug('importFile: setting up UI')
            if all(self.r2.surveys[0].df['irecip'].values == 0):
                self.importDataRecipBtn.show()
                self.recipOrNoRecipShow(recipPresence = False)
            else:
                self.recipOrNoRecipShow(recipPresence = True)
                self.importDataRecipBtn.hide()
                self.tabPreProcessing.setTabEnabled(2, True)
                self.filterAttrCombo.addItem('Reciprocal Error')
                self.plotError()
                self.errHist()
            if 'dev' in self.r2.surveys[0].df.columns:
                self.filterAttrCombo.addItem('Stacking Error (Dev.)')
                
            if self.r2.iremote is not None:
                if np.sum(self.r2.iremote) > 0:
                    self.meshQuadGroup.setEnabled(False)
                else:
                    self.meshQuadGroup.setEnabled(True)
            if self.boreholeCheck.isChecked() is True:
                self.r2.setBorehole(True)
            else:
                self.r2.setBorehole(False)
            self.plotManualFiltering()
            self.elecTable.initTable(self.r2.elec)
            self.tabImporting.setTabEnabled(1,True)
            if 'ip' in self.r2.surveys[0].df.columns:
                if np.sum(self.r2.surveys[0].df['ip'].values) > 0 or np.sum(self.r2.surveys[0].df['ip'].values) < 0: # np.sum(self.r2.surveys[0].df['ip'].values) !=0 will result in error if all the IP values are set to NaN
                    self.ipCheck.setChecked(True)
                if self.ftype == 'Syscal':
                    self.dcaButton.setEnabled(True)
                    self.dcaProgress.setEnabled(True)
            if np.isnan(self.r2.elec).any(): # for users who import messed up topography files (res2dinv mostly)
                self.topoInterpBtnFunc()
                self.updateElec()
            self.plotPseudo()
            self.loadingWidget(exitflag=True) #close the loading widget
            self.infoDump(fname + ' imported successfully')
            self.invNowBtn.setEnabled(True)
            self.activateTabs(True)
            self.nbElecEdit.setText(str(len(self.r2.elec)))
            self.elecDxEdit.setText('%s' %(self.r2.elec[~self.r2.iremote,:][1,0]-self.r2.elec[~self.r2.iremote,:][0,0]))
            self.fnamesCombo.hide()
            self.fnamesComboLabel.hide()
            
        except Exception as e:
            self.loadingWidget(exitflag=True)
            pdebug('importFile: ERROR:', e)
            self.errorDump('Importation failed. File is not being recognized. \
                      Make sure you have selected the right file type.')
            pass

    def restartFunc(self):
        pdebug('restartFunc: creating new R2 object')
        self.r2 = R2(self.newwd, typ=self.typ) # create new R2 instance
        '''actually we don't really need to instanciate a new object each
        time but it's safe otherwise we would have to reset all attributes
        , delete all surveys and clear parameters
        '''
        # reinitiate flags
        self.r2.iBatch = self.iBatch
        self.r2.setBorehole(self.iBorehole)
        self.r2.iTimeLapse = self.iTimeLapse
        self.r2.iForward = self.iForward
        if self.iTimeLapse is True:
            self.reg_mode.setCurrentIndex(2)
        else:
            self.reg_mode.setCurrentIndex(0)
        self.activateTabs(False)
        
        # importing
        self.parser = None
        self.plotAspect = 'equal'
        self.wdBtn.setText('Working directory:' + os.path.basename(self.r2.dirname))
        self.importDataBtn.setText('Import Data')
        self.ipCheck.setChecked(False)
        self.ipCheck.setEnabled(False)
        self.importDataRecipBtn.hide()
        self.fnamesCombo.hide()
        self.fnamesComboLabel.hide()
        self.psContourCheck.setEnabled(False)
        self.tabImporting.setTabEnabled(1, False)
        self.mwPseudo.clear() # clearing figure
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
        self.regionTable.reset()
        
        #forward model
        self.mwFwdPseudo.clear()
        self.mwFwdPseudoIP.clear()

        # inversion settings
        self.flux_type.setCurrentIndex(0)
        self.singular_type.setChecked(False)
        self.res_matrix.setCurrentIndex(1)
        self.scale.setText('1.0')
        self.patch_size_x.setText('1')
        self.patch_size_y.setText('1')
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
        self.logText.setText('')
        self.mwRMS.clear()
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
        self.doiCheck.setChecked(False)
        self.doiSensCheck.setChecked(False)
        self.displayParams = {'index':0,'edge_color':'none',
                            'sens':True, 'attr':'Resistivity(Ohm-m)',
                            'contour':False, 'vmin':None, 'vmax':None,
                            'cmap':'viridis', 'sensPrc':0.5,
                            'doi':self.modelDOICheck.isChecked(),
                            'doiSens':False,
                            'pvslices':([],[],[]), 'pvthreshold':None,
                            'pvgrid':False, 'pvcontour':[]}
        self.mwInv.clear()
        self.mwInvError.clear()
        self.mwInvError2.clear()

  
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
                                   <p><a href='https://gitlab.com/hkex/pyr2#downloads'>https://gitlab.com/hkex/pyr2</a></p>\
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
                webbrowser.open('https://gitlab.com/hkex/pyr2#gui-for-r2-family-code') # can add download link, when we have a direct dl link
    
    def checkWine(self): # check if wine is installed, on Unix system
        #check operating system
        OpSys=platform.system()
        wineInstalled = True
        #detect wine
        if OpSys == 'Linux':
            p = Popen("wine --version", stdout=PIPE, shell=True)
            is_wine = str(p.stdout.readline())
            if is_wine.find("wine") == -1:
                wineInstalled = False
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
                wineInstalled = False
        return wineInstalled


    def checkWineShow(self, wineInstalled):
        print('is wine installed?', wineInstalled)
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
                webbrowser.open('https://www.winehq.org/')
#                webbrowser.open('https://gitlab.com/hkex/pyr2#linux-and-mac-user')

    


if __name__ == '__main__':
    catchErrors()
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    app.setWindowIcon(QIcon(os.path.join(bundle_dir, 'logo.png'))) # that's the true app icon
    
    splash_pix = QPixmap(os.path.join(bundle_dir, 'loadingLogo.png'))
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setWindowFlags(Qt.WindowStaysOnTopHint | Qt.FramelessWindowHint)
    splash.setEnabled(False)

    # adding progress bar
    progressBar = QProgressBar(splash)
    progressBar.setMaximum(10)
    progressBar.setGeometry(100, splash_pix.height() - 50, splash_pix.width() - 200, 20)

    from resipy.R2 import ResIPy_version
    splash.show()
    splash.showMessage("Loading libraries", Qt.AlignBottom | Qt.AlignCenter, Qt.black)
    app.processEvents()

    # in this section all import are made except the one for pyQt
    progressBar.setValue(1)
    app.processEvents()

    import matplotlib
    matplotlib.use('Qt5Agg')
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
        pvfound = True
    except:
        pvfound = False
        print('WARNING: pyvista not found, 3D plotting capabilities will be limited.')

    from resipy.R2 import R2
    from resipy.r2help import r2help
    splash.showMessage("ResIPy is ready!", Qt.AlignBottom | Qt.AlignCenter, Qt.black)
    progressBar.setValue(10)
    app.processEvents()

    ex = App()
    ex.show()
    splash.hide() # hiding the splash screen when finished
    
    sys.exit(app.exec_())
