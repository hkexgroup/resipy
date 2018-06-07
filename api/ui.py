#!/usr/bin/python3
# -*- coding: utf-8 -*-
from PyQt5.QtWidgets import (QMainWindow, QApplication, QPushButton, QWidget, 
    QAction, QTabWidget,QVBoxLayout, QGridLayout, QLabel, QLineEdit, QMessageBox,
    QListWidget, QFileDialog, QCheckBox, QComboBox, QTextEdit, QSlider, QHBoxLayout,
    QTableWidget)
from PyQt5.QtGui import QIcon, QPixmap
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

sys.path.append(os.path.relpath('../api')) # not needed anymore

#from mesh_class import mesh_obj
#import meshTools as mt
#from meshTools import Mesh_obj
from R2 import R2


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
            

    
    def plot(self, callback):
        ''' call a callback plot function and give it the ax to plot to
        '''
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        callback(ax=ax)
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
        self.r2 = R2('/media/jkl/data/phd/tmp/r2gui/api/test')
        
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
            
        wd = QPushButton('Press here to set working directory')
        wd.clicked.connect(getwd)
        

        def getfile():
            fname, _ = QFileDialog.getOpenFileName(tab1,'Open File', directory=self.r2.dirname)
            if fname != '':
                self.fname = fname
                buttonf.setText(self.fname)
                self.r2.createSurvey(self.fname)
                if all(self.r2.surveys[0].df['irecip'].values == 0):
                    hbox4.addWidget(buttonfr)
                else:
                    plotError()
                    generateMesh()
                plotPseudo()
        
        buttonf = QPushButton('Import Data') 
        buttonf.clicked.connect(getfile)
        
        def getfileR():
            fnameRecip, _ = QFileDialog.getOpenFileName(tab1,'Open File', directory=self.r2.dirname)
            buttonfr.setText(fnameRecip)
            self.r2.surveys[0].addData(fnameRecip)
        
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
                plotPseudoIP()
                mwPseudoIP.setVisible(True)
                tabPreProcessing.setTabEnabled(2, True)
            else:
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
            mwPseudoIP.plot(self.r2.surveys[0].pseudo)
        
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
        tabs.addTab(tabPreProcessing, 'PreProcessing')
        
        manualLayout = QVBoxLayout()
        
        mwManualFiltering = MatplotlibWidget(navi=True)
        manualLayout.addWidget(mwManualFiltering)
        
        
        errorLayout = QVBoxLayout()
        
        def plotError():
            mwFitError.plot(self.r2.plotError)
            
        def fitLinError():
            mwFitError.plot(self.r2.linfit)
        
        def fitLmeError():
            print('NOT READY YET')
            mwFitError.plot(self.r2.lmefit)
        
        def fitModel(index):
            print(index)
            if index == 0:
                plotError()
            elif index == 1:
                fitLinError()
            elif index == 3:
                fitLmeError()
            else:
                print('NOT IMPLEMENTED YET')
        
        errFitType = QComboBox()
        errFitType.addItem('Observed Errors')
        errFitType.addItem('Linear')
        errFitType.addItem('Exponential')
        errFitType.addItem('Linear Mixed Effect')
        errFitType.currentIndexChanged.connect(fitModel)
        errorLayout.addWidget(errFitType)
        
        mwFitError = MatplotlibWidget(navi=True)
        errorLayout.addWidget(mwFitError)
       
        
        ipLayout = QVBoxLayout()
        
        mwIPFiltering = MatplotlibWidget(navi=True)
        ipLayout.addWidget(mwIPFiltering)
        
        
        manualWidget = QWidget()
        manualWidget.setLayout(manualLayout)
        tabPreProcessing.addTab(manualWidget, 'Manual Filtering')
        errorWidget = QWidget()
        errorWidget.setLayout(errorLayout)
        tabPreProcessing.addTab(errorWidget, 'Error-based filtering')
        ipWidget = QWidget()
        ipWidget.setVisible(False)
        ipWidget.setLayout(ipLayout)
        tabPreProcessing.addTab(ipWidget, 'IP filtering')
        tabPreProcessing.setTabEnabled(2, False)
        
        
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
            elif index == 1:
                self.r2.createMesh(typ='quad')
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
        grid = QGridLayout()
            
        singular_type = QCheckBox('Singularity Removal')
        grid.addWidget(singular_type, 0, 1)
        
        button = QPushButton('hello')
        button.clicked.connect(QMessageBox)
        grid.addWidget(button, 3, 0)
        
        # chose inversion type
        inv_type = QComboBox()
        inv_type.addItem('Regular One')
        inv_type.addItem('New one')
        inv_type.addItem('No a good one')
        grid.addWidget(inv_type, 4, 0)

        inv_type = QListWidget()        
        inv_type.addItem('Regular One')
        inv_type.addItem('New one')
        inv_type.addItem('No a good one')
        grid.addWidget(inv_type, 5, 0)
        
        
        generalSettings.setLayout(grid)
        tabInversionSettings.addTab(generalSettings, 'General')
        
        
        # advanced settings
        advancedSettings = QWidget()
        gridAdv = QGridLayout()
        
        btn = QPushButton('Press me')
        gridAdv.addWidget(btn)
        
        advancedSettings.setLayout(gridAdv)
        tabInversionSettings.addTab(advancedSettings, 'Advanced')

        
        
        
        #%% tab 5 INVERSION
        
        tabInversion = QWidget()
        tabs.addTab(tabInversion, 'Inversion')
        grid = QGridLayout()

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
            if self.r2.configFile == '':
                self.r2.write2in(param=param)
        
            self.r2.surveys[0].write2protocol(os.path.join(self.r2.dirname, 'protocol.dat'))
        
            if frozen == 'not':
                shutil.copy('R2.exe',
                    os.path.join(self.r2.dirname, 'R2.exe'))
            else:
                shutil.copy(os.path.join(bundle_dir, 'R2.exe'),
                    os.path.join(self.r2.dirname, 'R2.exe'))
            
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
#            mwInvResult.plot(self.r2.showResults)
            mwInvResult.plot(self.r2.showSection)


        btn = QPushButton('Invert')
        btn.clicked.connect(logInversion)
        
        grid.addWidget(btn, 0, 0)
        
        logText = QTextEdit()
        grid.addWidget(logText, 1, 0)
        
        mwInvResult = MatplotlibWidget(navi=True)
        grid.addWidget(mwInvResult, 2, 0)        
        
        
        tabInversion.setLayout(grid)
                    
        
        #%% tab 6 POSTPROCESSING
        tabPostProcessing = QWidget()
        tabs.addTab(tabPostProcessing, 'Post Processing')
        grid = QGridLayout()
        
        grid.addWidget(QLabel('TO BE IMPLEMENTED'), 0, 0)
        
        tabPostProcessing.setLayout(grid)
        
        
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
