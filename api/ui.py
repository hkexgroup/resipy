#!/usr/bin/python3
# -*- coding: utf-8 -*-
from PyQt5.QtWidgets import (QMainWindow, QApplication, QPushButton, QWidget, 
    QAction, QTabWidget,QVBoxLayout, QGridLayout, QLabel, QLineEdit, QMessageBox,
    QListWidget, QFileDialog, QCheckBox, QComboBox, QTextEdit, QSlider, QHBoxLayout,
    QTableWidget, QFormLayout)
from PyQt5.QtGui import QIcon, QPixmap, QIntValidator, QDoubleValidator
from PyQt5.QtCore import QThread, pyqtSignal, QProcess
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import os
import sys
import shutil
import platform
OS = platform.system()
from datetime import datetime

from matplotlib import rcParams
rcParams.update({'font.size': 13}) # CHANGE HERE for graph font size

sys.path.append(os.path.relpath('../api')) # not needed anymore

#from mesh_class import mesh_obj
#import meshTools as mt
#from meshTools import Mesh_obj
from R2 import R2
from r2help import r2help

# small code to see where are all the directories
frozen = 'not'
if getattr(sys, 'frozen', False):
        # we are running in a bundle
        frozen = 'ever so'
        bundle_dir = sys._MEIPASS
else:
        # we are running in a normal Python environment
        bundle_dir = os.path.dirname(os.path.abspath(__file__))
print( 'we are',frozen,'frozen')
print( 'bundle dir is', bundle_dir )
print( 'sys.argv[0] is', sys.argv[0] )
print( 'sys.executable is', sys.executable )
print( 'os.getcwd is', os.getcwd() )


class MatplotlibWidget(QWidget):
    def __init__(self, parent=None, figure=None, navi=False):
        super(MatplotlibWidget, self).__init__(parent) # we can pass a figure but we can replot on it when
        # pushing on a button (I didn't find a way to do it) while with the axes, you can still clear it and
        # plot again on them
        if figure is None:
            figure = Figure()
            axes = figure.add_subplot(111)
        else:
            axes = figure.get_axes()
        self.axis = axes
        self.figure = figure
        self.canvas = FigureCanvasQTAgg(self.figure)
        
        self.layoutVertical = QVBoxLayout(self)
        self.layoutVertical.addWidget(self.canvas)

        if navi is True:
            self.navi_toolbar = NavigationToolbar(self.canvas, self)
            self.navi_toolbar.setMaximumHeight(30)
            self.layoutVertical.addWidget(self.navi_toolbar)
    
    
    def setMinMax(self, vmin, vmax):
        coll = self.axis.collections[0]
        oldLimits = coll.get_clim()
        if vmin == '':
            vmin = oldLimits[0]
        else:
            vmin = int(vmin)
        if vmax == '':
            vmax = oldLimits[0]
        else:
            vmax = int(vmax)
        coll.set_clim(vmin, vmax)
        self.canvas.draw()

    
    def plot(self, callback):
        ''' call a callback plot function and give it the ax to plot to
        '''
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        self.axis = ax
        self.callback = callback
        callback(ax=ax)
        ax.set_aspect('auto')
        self.figure.tight_layout()
        self.canvas.draw()
    
    def replot(self, **kwargs):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        self.axis = ax
        callback = self.callback
        callback(ax=ax, **kwargs)
        ax.set_aspect('auto')
        self.figure.tight_layout()
        self.canvas.draw()
        
#
#    def draw(self, fig):
#        print('creating new figure')
#        self.figure = fig
#        self.figure.suptitle('hello')
#        self.canvas = FigureCanvasQTAgg(self.figure)
#        self.canvas.draw()
#        self.figure.canvas.draw()
#    
#    def drawAxes(self, ax):
#        print('draw Awes')
#        self.figure.clf()
#        self.figure.axes.append(ax)
##        self.canvas = FigureCanvasQTAgg(self.figure)
#        self.figure.canvas.draw()
#        self.canvas.draw()
    
#class ThreadSample(QThread):
#    newSample = pyqtSignal(list)
#    
#    def __init__(self, parent=None):
#        super(ThreadSample, self).__init__(parent)
#    
#    def run(self):
#        randomSample = random.sample(range(0, 10), 10)
#        
#        self.newSample.emit(randomSample)


class App(QMainWindow):
 
    def __init__(self, parent=None):
        super().__init__()
        self.setWindowTitle('R2 GUI')
        self.setGeometry(100,100,1000,600)
        self.r2 = R2(os.path.join(bundle_dir, 'invdir'))
        
        self.table_widget = QWidget()
        layout = QVBoxLayout()
        tabs = QTabWidget()
        
        #%% tab 1 importing data
        tab1 = QWidget()
        tabs.addTab(tab1, 'Importing')
        gridImport = QGridLayout()
        topLayout = QHBoxLayout()
        
        # meta data (title and date of survey)
        title = QLabel('Title')
        titleEdit = QLineEdit()
        titleEdit.setText('My beautiful survey')
        
        date = QLabel('Date')
        dateEdit = QLineEdit()
        dateEdit.setText(datetime.now().strftime('%Y-%m-%d')) # get today date
        
        hbox1 = QHBoxLayout()
        hbox1.addWidget(title)
        hbox1.addWidget(titleEdit)
        
        hbox2 = QHBoxLayout()
        hbox2.addWidget(date)
        hbox2.addWidget(dateEdit)
        
        # ask for working directory, and survey file to input
        def getwd():
            fdir = QFileDialog.getExistingDirectory(tab1, 'Choose Working Directory')
            self.r2.setwd(fdir)
            wd.setText(fdir)
            
        wd = QPushButton('Working directory:' + self.r2.dirname + ' (Press to change)')
        wd.clicked.connect(getwd)
        

        def getfile():
            fname, _ = QFileDialog.getOpenFileName(tab1,'Open File', directory=self.r2.dirname)
            if len(self.r2.surveys) > 0: # will need to change this for time-lapse
                self.r2.surveys = []
            if fname != '':
                self.fname = fname
                buttonf.setText(self.fname)
                self.r2.createSurvey(self.fname)
                if all(self.r2.surveys[0].df['irecip'].values == 0):
                    hbox4.addWidget(buttonfr)
                else:
                    tabPreProcessing.setTabEnabled(1, True)
                    plotError()
                generateMesh()
                plotPseudo()
        
        buttonf = QPushButton('Import Data') 
        buttonf.clicked.connect(getfile)
        
        def getfileR():
            fnameRecip, _ = QFileDialog.getOpenFileName(tab1,'Open File', directory=self.r2.dirname)
            buttonfr.setText(fnameRecip)
            self.r2.surveys[0].addData(fnameRecip)
            if all(self.r2.surveys[0].df['irecip'].values == 0) is False:
                tabPreProcessing.setTabEnabled(1, True) # no point in doing error processing if there is no reciprocal
                plotError()

        buttonfr = QPushButton('If you have reciprocals upload them here') 
        buttonfr.clicked.connect(getfileR)
        
        hbox4 = QHBoxLayout()
        hbox4.addWidget(buttonf)
        
        def diplayTopo(state):
            if state  == Qt.Checked:
                elecTable.setVisible(True)
            else:
                elecTable.setVisible(False)
        
        def diplayPseudoIP(state):
            if state  == Qt.Checked:
                self.r2.typ = 'cR2'
                plotPseudoIP()
                phaseplotError()
                showIpOptions(True)
                mwPseudoIP.setVisible(True)
                tabPreProcessing.setTabEnabled(2, True)
            else:
                self.r2.typ = 'R2'
                showIpOptions(False)
                mwPseudoIP.setVisible(False)
                tabPreProcessing.setTabEnabled(2, False)

        ipCheck = QCheckBox('Induced Polarization')
        ipCheck.stateChanged.connect(diplayPseudoIP)
        topoCheck = QCheckBox('Topography')
        topoCheck.stateChanged.connect(diplayTopo)
        hbox5 = QHBoxLayout()
        hbox5.addWidget(ipCheck)
        hbox5.addWidget(topoCheck)
        
        metaLayout = QVBoxLayout()
        metaLayout.addLayout(hbox1)
        metaLayout.addLayout(hbox2)
        metaLayout.addWidget(wd)
        metaLayout.addLayout(hbox4)
        metaLayout.addLayout(hbox5)
        topLayout.addLayout(metaLayout)

        # electrode table
        # TODO add copy paste functionnality
        elecTable = QTableWidget()
        elecTable.setRowCount(10)
        elecTable.setColumnCount(3)
        elecTable.setVisible(False)
        topLayout.addWidget(elecTable)
        gridImport.addLayout(topLayout, 0, 0)        
        
        def plotPseudo():
            mwPseudo.plot(self.r2.surveys[0].pseudo)
        
        def plotPseudoIP():
            mwPseudoIP.plot(self.r2.surveys[0].pseudoIP)
        
        pseudoLayout = QHBoxLayout()

        mwPseudo = MatplotlibWidget(navi=True)
        pseudoLayout.addWidget(mwPseudo)
                
        mwPseudoIP = MatplotlibWidget(navi=True)
        mwPseudoIP.setVisible(False)
        pseudoLayout.addWidget(mwPseudoIP)
        
        gridImport.addLayout(pseudoLayout, 1, 0)
        
#        def plotError():
#            mwError.plot(self.r2.surveys[0].plotError)
#            
#        mwError = MatplotlibWidget()
#        grid.addWidget(mwError, 4, 1)
        
        tab1.setLayout(gridImport)
        
    
        
        #%% tab 2 PRE PROCESSING
        tabPreProcessing = QTabWidget()
        tabs.addTab(tabPreProcessing, 'Pre-processing')
        
        manualLayout = QVBoxLayout()
        
        def plotManualFiltering():
            mwManualFiltering.plot(self.r2.surveys[0].manualFilter)
        
        notice = QLabel('Press "Start" and then click on the dots to select them. Press "Done" to remove them.')
        manualLayout.addWidget(notice)
        
        btnLayout = QHBoxLayout()
        btnStart = QPushButton('Start')
        btnStart.clicked.connect(plotManualFiltering)
        btnLayout.addWidget(btnStart)
        btnDone = QPushButton('Done')
        btnDone.clicked.connect(lambda x: print('TO BE IMPLEMENTED'))
        btnLayout.addWidget(btnDone)
        manualLayout.addLayout(btnLayout)
        
        
        mwManualFiltering = MatplotlibWidget(navi=True)
        manualLayout.addWidget(mwManualFiltering)
        
        
        errorLayout = QVBoxLayout()
        
        def plotError():
            mwFitError.plot(self.r2.plotError)
            self.r2.errTyp = 'obs'
            
        def fitLinError():
            mwFitError.plot(self.r2.linfit)
            self.r2.errTyp = 'lin'
        
        def fitLmeError():
            print('NOT READY YET')
            mwFitError.plot(self.r2.lmefit)
            self.r2.errTyp = 'lme'
        
        def fitpwl():
            mwFitError.plot(self.r2.pwlfit)
            self.r2.errTyp = 'pwl'

        def fitModel(index):
            print(index)
            if index == 0:
                plotError()
            elif index == 1:
                fitLinError()
            elif index == 2:
                fitpwl()
            elif index == 3:
                fitLmeError()
            else:
                print('NOT IMPLEMENTED YET')
        
        errFitType = QComboBox()
        errFitType.addItem('Observed Errors')
        errFitType.addItem('Linear')
        errFitType.addItem('Power-law')
        errFitType.addItem('Exponential')
        errFitType.addItem('Linear Mixed Effect')
        errFitType.currentIndexChanged.connect(fitModel)
        errorLayout.addWidget(errFitType)
        
        mwFitError = MatplotlibWidget(navi=True)
        errorLayout.addWidget(mwFitError)
       
        
        ipLayout = QVBoxLayout()
        
        def phaseplotError():
            mwIPFiltering.plot(self.r2.phaseplotError)
        
        def phasePLerr():
            mwIPFiltering.plot(self.r2.plotIPFit)
            
        def ipfitModel(index):
            print(index)
            if index == 0:
                phaseplotError()
            elif index == 1:
                phasePLerr()
            else:
                print('NOT IMPLEMENTED YET')
            
        iperrFitType = QComboBox()
        iperrFitType.addItem('Observed discrepancies') ##### BY default does not show!! should be selected after the power law (don't know why!!!)
        iperrFitType.addItem('Power law')
        iperrFitType.currentIndexChanged.connect(ipfitModel)
        ipLayout.addWidget(iperrFitType)
        
        mwIPFiltering = MatplotlibWidget(navi=True)
        ipLayout.addWidget(mwIPFiltering)
        
        def dcaFiltering():
            self.r2.surveys[0].dca(dump=dcaProgress.setText)
            
        dcaLayout = QHBoxLayout()
        dcaButton = QPushButton('DCA filtering')
        dcaButton.clicked.connect(dcaFiltering)
        dcaProgress = QLineEdit('0%')
        dcaProgress.setReadOnly(True)
        dcaLayout.addWidget(dcaButton)
        dcaLayout.addWidget(dcaProgress)
        ipLayout.addLayout(dcaLayout)
        
        
        manualWidget = QWidget()
        manualWidget.setLayout(manualLayout)
        tabPreProcessing.addTab(manualWidget, 'Manual Filtering')
        errorWidget = QWidget()
        errorWidget.setLayout(errorLayout)
        tabPreProcessing.addTab(errorWidget, 'Resistance Error Model')
        ipWidget = QWidget()
        ipWidget.setVisible(False)
        ipWidget.setLayout(ipLayout)
        tabPreProcessing.addTab(ipWidget, 'Phase Error Model')
        tabPreProcessing.setTabEnabled(2, False)
        tabPreProcessing.setTabEnabled(1, False)
        
        
        #%% tab MESH
        tabMesh= QWidget()
        tabs.addTab(tabMesh, 'Mesh')
        meshLayout = QVBoxLayout()
                
        def callback2(ax):
            ax.plot(np.random.randn(20,5), '+--')
            ax.set_title('Random data nnnnndfghdfh')

        def generateMesh(index=0):
            if index == 0:
                self.r2.createMesh(typ='quad')
                scale.setVisible(False)
                scaleLabel.setVisible(False)
            elif index == 1:
                scale.setVisible(True)
                scaleLabel.setVisible(True)
                self.r2.createMesh(typ='trian')
                # TODO to implemente the triangular mesh
            else:
                print('NOT IMPLEMENTED')
            print(self.r2.mesh.summary())
            mwMesh.plot(self.r2.mesh.show) # mesh.show() is the callback function to be called with ax
        
        
        meshType = QComboBox()
        meshType.addItem('Quadrilateral Mesh')
        meshType.addItem('Triangular Mesh')
        meshType.currentIndexChanged.connect(generateMesh)
        meshLayout.addWidget(meshType)
        
        # TODO EVENTUALLY SHOW MESH OPTION HERE
        
        mwMesh = MatplotlibWidget(navi=True)
        meshLayout.addWidget(mwMesh)
        
        '''
        def changeValue(value):
            print(value)
            self.r2.createMesh(elemy=value)
            plotQuadMesh()
            
        def plotQuadMesh():
            print('plotQuadMesh')
            mwqm.plot(self.r2.mesh.show)
        
        def genQuadMesh():
            self.r2.createMesh()
            
        btn = QPushButton('QuadMesh')
        btn.clicked.connect(genQuadMesh)
        grid.addWidget(btn, 3, 0)
        
        sld = QSlider(Qt.Horizontal)
        sld.setGeometry(30, 40, 100, 30)
        sld.valueChanged[int].connect(changeValue)
        grid.addWidget(sld, 4, 0)
        
        mwqm = MatplotlibWidget(navi=True)
        grid.addWidget(mwqm, 5, 0)
        '''
       
        tabMesh.setLayout(meshLayout)
        
        
        #%% tab INVERSION SETTINGS
        tabInversionSettings = QTabWidget()
        tabs.addTab(tabInversionSettings, 'Inversion settings')

        # general tab
        generalSettings = QWidget()
        generalLayout = QHBoxLayout()
        invForm = QFormLayout()
        
        # advanced tab
        advancedSettings = QWidget()
        advancedLayout = QHBoxLayout()
        advForm = QFormLayout()
        
        def showIpOptions(arg):
            opts = [c_wgt, c_wgtLabel, d_wgt, d_wgtLabel,
                    singular_type, singular_typeLabel,
                    res_matrix, res_matrixLabel]
            iopts = arg*[1, 1, 1, 1, 0, 0, 0, 0]
            [o.setVisible(io) for o, io in zip(opts, iopts)]
        
        # help sections
        def showHelp(arg):
            if arg not in r2help:
                helpSection.setText('SORRY NOT IN HELP')
            else:
                helpSection.setHtml(r2help[arg])
        
        def showHelp2(arg):
            if arg not in r2help:
                helpSection2.setText('SORRY NOT IN HELP')
            else:
                helpSection2.setHtml(r2help[arg])
        
        
        def job_typeFunc(index):
            self.r2.param['job_type'] = index
        job_type = QComboBox()
        job_type.addItem('Inversion [1]')
        job_type.addItem('Forward [0]')
        job_type.currentIndexChanged.connect(job_typeFunc)
        invForm.addRow(QLabel('Job Type:'), job_type)
        
        def flux_typeFunc(index):
            if index == 0:
                self.r2.param['flux_type'] = 3
            else:
                self.r2.param['flux_type'] = 2
        flux_typeLabel = QLabel('<a href="flux_type">Flux Type</a>:')
        flux_typeLabel.linkActivated.connect(showHelp)
        flux_type = QComboBox()
        flux_type.addItem('3D')
        flux_type.addItem('2D')
        flux_type.currentIndexChanged.connect(flux_typeFunc)
        invForm.addRow(flux_typeLabel, flux_type)
        
        def singular_typeFunc(state):
            if state == Qt.Checked:
                self.r2.param['singular_type'] = 1
            else:
                self.r2.param['singular_type'] = 0
        singular_typeLabel = QLabel('<a href="singular_type">Remove Singularity</a>')
        singular_typeLabel.linkActivated.connect(showHelp)
        singular_type = QCheckBox()
        singular_type.stateChanged.connect(singular_typeFunc)
        invForm.addRow(singular_typeLabel, singular_type)
        
        def res_matrixFunc(index):
            self.r2.param['res_matrix'] = index
        res_matrixLabel = QLabel('<a href="res_matrix">Value for <code>res_matrix</code><a/>')
        res_matrixLabel.linkActivated.connect(showHelp2)
        res_matrix = QComboBox()
        res_matrix.addItem('No sensisitivity/resolution matrix [0]')
        res_matrix.addItem('Sensitivity matrix [1]')
        res_matrix.addItem('True Resolution Matrix [2]')
        res_matrix.addItem('Sensitivity map [3]')
        res_matrix.setCurrentIndex(1)
        res_matrix.currentIndexChanged.connect(res_matrixFunc)
        advForm.addRow(res_matrixLabel, res_matrix)
        
        def scaleFunc():
            self.r2.param['scale'] = float(scale.text())
        scaleLabel = QLabel('<a href="scale"> Scale for triangular mesh</a>:')
        scaleLabel.linkActivated.connect(showHelp)
        scaleLabel.setVisible(False)
        scale = QLineEdit()
        scale.setValidator(QDoubleValidator())
        scale.setText('1.0')
        scale.editingFinished.connect(scaleFunc)
        scale.setVisible(False)
        invForm.addRow(scaleLabel, scale)
        
        # put in adv
        def patch_size_xFunc():
            self.r2.param['patch_size_x'] = int(patch_size_x.text())
        patch_size_xLabel = QLabel('<a href="patch">Patch size x<a/>:')
        patch_size_xLabel.linkActivated.connect(showHelp2)
        patch_size_x = QLineEdit()
        patch_size_x.setValidator(QIntValidator())
        patch_size_x.setText('1')
        patch_size_x.editingFinished.connect(patch_size_xFunc)
        advForm.addRow(patch_size_xLabel, patch_size_x)

        # put in adv
        def patch_size_yFunc():
            self.r2.param['patch_size_y'] = int(patch_size_y.text())
        patch_size_yLabel = QLabel('<a href="patch">Patch size y<a/>:')
        patch_size_yLabel.linkActivated.connect(showHelp2)
        patch_size_y = QLineEdit()
        patch_size_y.setValidator(QIntValidator())
        patch_size_y.setText('1')
        patch_size_y.editingFinished.connect(patch_size_yFunc)
        advForm.addRow(patch_size_yLabel, patch_size_y)
        
        def inv_typeFunc(index):
            self.r2.param['inversion_type'] = index
            opts = [data_typeLabel, data_type,
                    reg_modeLabel, reg_mode,
                    toleranceLabel, tolerance,
                    max_iterationsLabel, max_iterations,
                    error_modLabel, error_mod,
                    alpha_anisoLabel, alpha_aniso,
                    a_wgtLabel, a_wgt,
                    b_wgtLabel, b_wgt]
            if index == 3: # qualitative solution
                [o.setVisible(False) for o in opts]
            else:
                [o.setVisible(True) for o in opts]
        inv_typeLabel = QLabel('<a href="inverse_type">Inversion Type</a>:')
        inv_typeLabel.linkActivated.connect(showHelp)
        inv_type = QComboBox()
        inv_type.addItem('Pseudo Marquardt [0]')
        inv_type.addItem('Regularized Inversion with Linear Filtering [1]')
        inv_type.addItem('Regularized Inversion with Quadratic Filtering [2]')
        inv_type.addItem('Qualitative Solution [3]')
        inv_type.addItem('Blocked Linear Regularized Inversion [4]')
        inv_type.setCurrentIndex(1)
        inv_type.currentIndexChanged.connect(inv_typeFunc)
        invForm.addRow(inv_typeLabel, inv_type)
        
        def data_typeFunc(index):
            self.r2.param['data_type'] = index
        data_typeLabel = QLabel('<a href="data_type">Data type</a>:')
        data_typeLabel.linkActivated.connect(showHelp)
        data_type = QComboBox()
        data_type.addItem('Normal [0]')
        data_type.addItem('Logarithmic [1]')
        data_type.currentIndexChanged.connect(data_typeFunc)
        invForm.addRow(data_typeLabel, data_type)
        
        def reg_modeFunc(index):
            self.r2.param['reg_mode'] = index
        reg_modeLabel = QLabel('<a href="reg_mode">Regularization mode</a>:')
        reg_modeLabel.linkActivated.connect(showHelp)
        reg_mode = QComboBox()
        reg_mode.addItem('Normal regularization [0]')
        reg_mode.addItem('Regularization from initial model [1]')
        reg_mode.addItem('Regularization from difference inversion [2]')
        reg_mode.currentIndexChanged.connect(reg_modeFunc)
        invForm.addRow(reg_modeLabel, reg_mode)
        
        def toleranceFunc():
            self.r2.param['tolerance'] = float(tolerance.text())
        toleranceLabel = QLabel('<a href="tolerance">Value for tolerance</a>:')
        toleranceLabel.linkActivated.connect(showHelp)
        tolerance = QLineEdit()
        tolerance.setValidator(QDoubleValidator())
        tolerance.setText('1.0')
        tolerance.editingFinished.connect(toleranceFunc)
        invForm.addRow(toleranceLabel, tolerance)
        
        def max_iterationsFunc():
            self.r2.param['max_iterations'] = int(max_iterations.text())
        max_iterationsLabel = QLabel('<a href="max_iterations">Maximum number of iterations</a>:')
        max_iterationsLabel.linkActivated.connect(showHelp)
        max_iterations = QLineEdit()
        max_iterations.setValidator(QIntValidator())
        max_iterations.setText('10')
        max_iterations.editingFinished.connect(max_iterationsFunc)
        invForm.addRow(max_iterationsLabel, max_iterations)
        
        def error_modFunc(index):
            self.r2.param['error_mod'] = index
        error_modLabel = QLabel('<a href="error_mod">Update the weights</a>:')
        error_modLabel.linkActivated.connect(showHelp2)
        error_mod = QComboBox()
        error_mod.addItem('Keep the same weights [0]')
        error_mod.addItem('Update the weights [1]')
        error_mod.addItem('Update the weights (recommended) [2]')
        error_mod.setCurrentIndex(1)
        error_mod.currentIndexChanged.connect(error_modFunc)
        advForm.addRow(error_modLabel, error_mod)
        
        
#        def createLabel(helpTag, title): # doesn't seem to work
#            ql = QLabel('<a href="' + helpTag + '>' + title + '</a>:')
#            ql.linkActivated.connect(showHelp)
#            return ql

        def alpha_anisoFunc():
            self.r2.param['alpha_aniso'] = float(alpha_aniso.text())
        alpha_anisoLabel = QLabel('<a href="alpha_aniso">Value for <code>alpha_aniso</code></a>:')
        alpha_anisoLabel.linkActivated.connect(showHelp2)
        alpha_aniso = QLineEdit()
        alpha_aniso.setValidator(QDoubleValidator())
        alpha_aniso.setText('1.0')
        alpha_aniso.editingFinished.connect(alpha_anisoFunc)
        advForm.addRow(alpha_anisoLabel, alpha_aniso)
        
        def a_wgtFunc():
            self.r2.param['a_wgt'] = float(a_wgt.text())
        a_wgtLabel = QLabel('<a href="errorParam"><code>a_wgt</code></a>:')
        a_wgtLabel.linkActivated.connect(showHelp)
        a_wgt = QLineEdit()
        a_wgt.setValidator(QDoubleValidator())
        a_wgt.setText('0.0')
        a_wgt.editingFinished.connect(a_wgtFunc)
        invForm.addRow(a_wgtLabel, a_wgt)
        
        def b_wgtFunc():
            self.r2.param['b_wgt'] = float(b_wgt.text())
        b_wgtLabel = QLabel('<a href="errorParam"><code>b_wgt</code></a>:')
        b_wgtLabel.linkActivated.connect(showHelp)
        b_wgt = QLineEdit()
        b_wgt.setValidator(QDoubleValidator())
        b_wgt.setText('0.0')
        b_wgt.editingFinished.connect(b_wgtFunc)
        invForm.addRow(b_wgtLabel, b_wgt)
        
        def c_wgtFunc():
            self.r2.param['b_wgt'] = float(c_wgt.text())
        c_wgtLabel = QLabel('<a href="errorParam"><code>c_wgt</code></a>:')
        c_wgtLabel.linkActivated.connect(showHelp)
        c_wgtLabel.setVisible(False)
        c_wgt = QLineEdit()
        c_wgt.setValidator(QDoubleValidator())
        c_wgt.setText('2')
        c_wgt.editingFinished.connect(c_wgtFunc)
        c_wgt.setVisible(False)
        invForm.addRow(c_wgtLabel, c_wgt)
        
        def d_wgtFunc():
            self.r2.param['d_wgt'] = float(d_wgt.text())
        d_wgtLabel = QLabel('<a href="errorParam"><code>b_wgt</code></a>:')
        d_wgtLabel.linkActivated.connect(showHelp)
        d_wgtLabel.setVisible(False)
        d_wgt = QLineEdit()
        d_wgt.setValidator(QDoubleValidator())
        d_wgt.setText('1')
        d_wgt.editingFinished.connect(b_wgtFunc)
        d_wgt.setVisible(False)
        invForm.addRow(d_wgtLabel, d_wgt)
        
        def rho_minFunc():
            self.r2.param['rho_min'] = float(rho_min.text())
        rho_minLabel = QLabel('<a href="errorParam">Minimum apparent resistivity</a>:')
        rho_minLabel.linkActivated.connect(showHelp)
        rho_min = QLineEdit()
        rho_min.setValidator(QDoubleValidator())
        rho_min.setText('-1000')
        rho_min.editingFinished.connect(rho_minFunc)
        invForm.addRow(rho_minLabel, rho_min)

        def rho_maxFunc():
            self.r2.param['rho_max'] = float(rho_max.text())
        rho_maxLabel = QLabel('<a href="errorParam">Maximum apparent resistivity</a>:')
        rho_maxLabel.linkActivated.connect(showHelp)
        rho_max = QLineEdit()
        rho_max.setValidator(QDoubleValidator())
        rho_max.setText('1000')
        rho_max.editingFinished.connect(rho_maxFunc)
        invForm.addRow(rho_maxLabel, rho_max)        
        
        def target_decreaseFunc():
            self.r2.param['target_decrease'] = float(target_decrease.text())
        target_decreaseLabel = QLabel('<a href="target_decrease">Target decrease</a>:')
        target_decreaseLabel.linkActivated.connect(showHelp)
        target_decrease = QLineEdit()
        target_decrease.setValidator(QDoubleValidator())
        target_decrease.setText('1.0')
        target_decrease.editingFinished.connect(target_decreaseFunc)
        invForm.addRow(target_decreaseLabel, target_decrease)
                
        
        generalLayout.addLayout(invForm)
        
        helpSection = QTextEdit('Help will be display here')
        helpSection.setReadOnly(True)
        helpSection.setText('Click on the labels and help will be displayed here')
        generalLayout.addWidget(helpSection)
        
        generalSettings.setLayout(generalLayout)
        tabInversionSettings.addTab(generalSettings, 'General')
        
        
        advancedLayout.addLayout(advForm)
        
        helpSection2 = QTextEdit('Help will be display here')
        helpSection2.setReadOnly(True)
        helpSection2.setText('Click on the labels and help will be displayed here')
        advancedLayout.addWidget(helpSection2)
        
        advancedSettings.setLayout(advancedLayout)
        tabInversionSettings.addTab(advancedSettings, 'Advanced')

        
        #%% tab 5 INVERSION
        
        tabInversion = QWidget()
        tabs.addTab(tabInversion, 'Inversion')
        invLayout = QVBoxLayout()

#        def callback(ax):
#            ax.plot(np.random.randn(20,5), '*-')
#            ax.set_title('Random data')
#        
#        def updateGraph():
#            mw.plot(callback)
#
#        btn = QPushButton('Draw random sample')
#        btn.clicked.connect(updateGraph)
#        mw = MatplotlibWidget()
#        
#        grid.addWidget(btn, 1, 0)
#        grid.addWidget(mw, 2, 0)
#        
        
        # run R2
        def dataReady():
            cursor = logText.textCursor()
            cursor.movePosition(cursor.End)
#            cursor.insertText(text)
            cursor.insertText(str(self.process.readAll(), 'utf-8'))
            logText.ensureCursorVisible()
            
        def logInversion():
#            self.r2.invert(callback=dataReady)
#            dataReady('kk\n')
            param = {} # TODO to be set in the previous tab
            if 'mesh' not in self.r2.param:
                generateMesh() # that will call mesh creation
        
            # write configuration file
#            if self.r2.configFile == '':
            self.r2.write2in(param=param)
        
            self.r2.surveys[0].write2protocol(os.path.join(self.r2.dirname, 'protocol.dat'))
            exeName = self.r2.typ + '.exe'
            
            if frozen == 'not':
                shutil.copy(os.path.join('./exe', exeName),
                    os.path.join(self.r2.dirname, exeName))
            else:
                shutil.copy(os.path.join(bundle_dir, 'exe', exeName),
                    os.path.join(self.r2.dirname, exeName))
            
            self.process = QProcess(self)
            self.process.setWorkingDirectory(self.r2.dirname)
            # QProcess emits `readyRead` when there is data to be read
            self.process.readyRead.connect(dataReady)
            if OS == 'Linux':
                self.process.start('wine R2.exe')
            else:
                wdpath = "\"" + os.path.join(self.r2.dirname, 'R2.exe').replace('\\','/') + "\""
                self.process.start(wdpath) # need absolute path and escape quotes (if space in the path)
            self.process.finished.connect(plotSection)
        
        def plotSection():
            mwInvResult.plot(self.r2.showResults)
            plotInvError()
            # TODO if we want to plot different attribute
#            mwInvResult.plot(self.r2.showSection)
        
        def replotSection(index=-1):
            print(index)
            mwInvResult.plot(self.r2.showResults) # TODO PASS ATTRIBUTE
        
        def setCBarLimit():
            print(vmaxEdit.text())
            mwInvResult.setMinMax(vmaxEdit.text(), vminEdit.text())
            
        def showEdges(status):
            if status == Qt.Checked:
                mwInvResult.replot(edge_color='k')
            else:
                mwInvResult.replot(edge_color='none')
                

        btn = QPushButton('Invert')
        btn.clicked.connect(logInversion)
        invLayout.addWidget(btn)
        
        logText = QTextEdit()
        logText.setReadOnly(True)
        invLayout.addWidget(logText)
        
        # option for display
        displayOptions = QHBoxLayout()
        attributeName = QComboBox()
        attributeName.addItem('Log(Resistivity)')
        attributeName.addItem('Resistivity')
        attributeName.addItem('Sensitivity')
        attributeName.addItem('Phase')
        attributeName.currentIndexChanged.connect(plotSection)
        displayOptions.addWidget(attributeName)
        
        vminLabel = QLabel('Min:')
        vminEdit = QLineEdit()
        vminEdit.setValidator(QIntValidator())
        vminEdit.editingFinished.connect(setCBarLimit)
        vmaxLabel = QLabel('Max:')
        vmaxEdit = QLineEdit()
        vmaxEdit.editingFinished.connect(setCBarLimit)
        vmaxEdit.setValidator(QIntValidator())
        displayOptions.addWidget(vminLabel)
        displayOptions.addWidget(vminEdit)
        displayOptions.addWidget(vmaxLabel)
        displayOptions.addWidget(vmaxEdit)
        
        showEdge = QCheckBox('Show edges')
        showEdge.setChecked(True)
        showEdge.stateChanged.connect(showEdges)
        displayOptions.addWidget(showEdge)
        
        contour = QCheckBox('Contour')
        contour.stateChanged.connect(lambda x : print(x))
        displayOptions.addWidget(contour)
        
        showSens = QCheckBox('Sensitivity overlay')
        showSens.stateChanged.connect(lambda x : print(x))
        displayOptions.addWidget(showSens)
        
        invLayout.addLayout(displayOptions)
        
        mwInvResult = MatplotlibWidget(navi=True)
        invLayout.addWidget(mwInvResult)
        
        
        tabInversion.setLayout(invLayout)
                    
        
        #%% tab 6 POSTPROCESSING
        tabPostProcessing = QTabWidget()
        tabs.addTab(tabPostProcessing, 'Post-processing')
        
        invError = QWidget()
        tabPostProcessing.addTab(invError, 'Inversion Errors')
        invErrorLayout = QVBoxLayout()
        
        def plotInvError():
            mwInvError.plot(self.r2.pseudoError)
            
        mwInvError = MatplotlibWidget(navi=True)
        invErrorLayout.addWidget(mwInvError)
        
        
        invError.setLayout(invErrorLayout)
        
        
        # add Jimmy graph
        
        
        
        #%%
        layout.addWidget(tabs)
        self.table_widget.setLayout(layout)
        self.setCentralWidget(self.table_widget)
        self.show()

''' 
class MyTableWidget(QWidget):        
 
    def __init__(self, parent):   
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)
 
        # Initialize tab screen
        self.tabs = QTabWidget()
        self.tab1 = QWidget()	
        self.tab2 = QWidget()
        self.tabs.resize(300,200) 
 
        # Add tabs
        self.tabs.addTab(self.tab1,"Tab 1")
        self.tabs.addTab(self.tab2,"Tab 2")
 
        # Create first tab
        self.tab1.layout = QVBoxLayout(self)
        self.pushButton1 = QPushButton("PyQt5 button")
        self.tab1.layout.addWidget(self.pushButton1)
        self.tab1.setLayout(self.tab1.layout)
 
        # Add tabs to widget        
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
 
'''

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())
