#!/usr/bin/python3
# -*- coding: utf-8 -*-
from PyQt5.QtWidgets import (QMainWindow, QApplication, QPushButton, QWidget, 
    QAction, QTabWidget,QVBoxLayout, QGridLayout, QLabel, QLineEdit, QMessageBox,
    QListWidget, QFileDialog, QCheckBox, QComboBox, QTextEdit)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import QThread, pyqtSignal, QProcess


from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import os
import sys
import platform
OS = platform.system()

sys.path.append(os.path.relpath('../api'))

#from mesh_class import mesh_obj
import meshTools as mt
from meshTools import Mesh_obj
from R2 import R2


class MatplotlibWidget(QWidget):
    def __init__(self, parent=None, figure=None):
        super(MatplotlibWidget, self).__init__(parent)
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
    
    def plot(self, callback):
        ''' call a callback plot function and give it the ax to plot to
        '''
#        print('plot is called')
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        callback(ax=ax)
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
        grid = QGridLayout()
        
        title = QLabel('Title')
        
        def getwd():
            fdir = QFileDialog.getExistingDirectory(tab1, 'Choose Working Directory')
            self.r2.setwd(fdir)
            wdEdit.setText(fdir)
            
        wd = QPushButton('Working directory')
        wd.clicked.connect(getwd)
        
        titleEdit = QLineEdit()
        titleEdit.setText('My beautiful survey')
        wdEdit = QLineEdit()
                
        grid.addWidget(title, 1, 0)
        grid.addWidget(titleEdit, 1, 1)
        
        grid.addWidget(wd, 2, 0)
        grid.addWidget(wdEdit, 2, 1)
        

        def getfile():
            self.fname, _ = QFileDialog.getOpenFileName(tab1,'Open File', directory=self.r2.dirname)
            buttonf.setText(self.fname)
            self.r2.createSurvey(self.fname)
            plotError()
        
        buttonf = QPushButton('Import Data') 
        buttonf.clicked.connect(getfile)
        grid.addWidget(buttonf, 3, 0)
        
        
        def plotPseudo():
            mwPseudo.plot(self.r2.surveys[0].pseudo)
        
        btn = QPushButton('Plot Pseudo')
        btn.clicked.connect(plotPseudo)
        grid.addWidget(btn, 3, 1)

        mwPseudo = MatplotlibWidget()
        grid.addWidget(mwPseudo, 4, 0)
        
        def plotError():
            mwError.plot(self.r2.surveys[0].plotError)
            
        mwError = MatplotlibWidget()
        grid.addWidget(mwError, 4, 1)
        
        tab1.setLayout(grid)
        
    
        
        #%% tab 2 PRE PROCESSING
        tabPreProcessing = QWidget()
        tabs.addTab(tabPreProcessing, 'PreProcessing')
        grid = QGridLayout()
 
        def fitLinError():
            mwFitError.plot(self.r2.surveys[0].linfit)
        
        def fitLmeError():
            print('NOT READY YET')
#            mwFitError.plot(self.r2.surveys[0].lmefit)
        
        btn = QPushButton('Fit Linear model')
        btn.clicked.connect(fitLinError)
        grid.addWidget(btn, 0, 0)
        
        btn = QPushButton('Fit LME model')
        btn.clicked.connect(fitLmeError)
        grid.addWidget(btn, 0, 1)
               
        mwFitError = MatplotlibWidget()
        grid.addWidget(mwFitError, 2, 1)
        
        tabPreProcessing.setLayout(grid)
        
        
        #%% tab MESH
        tab3 = QWidget()
        tabs.addTab(tab3, 'Mesh')
        grid = QGridLayout()

                
        def callback2(ax):
            ax.plot(np.random.randn(20,5), '+--')
            ax.set_title('Random data nnnnndfghdfh')
        
#        def updateGraph2():
#            mw1.plot(callback2)
            
        def generateMesh():
#            mesh_dict=mt.vtk_import(self.fname)#makes a dictionary of a mesh 
#            mesh = mesh_obj.mesh_dict2obj(mesh_dict)# this is a mesh_obj class instance 
            self.r2.createMesh(typ='quad')
            print(self.r2.mesh.summary())
            mw1.plot(self.r2.mesh.show) # mesh.show() is the callback function to be called with ax
            
        btn = QPushButton('Generate Quadrilateral Mesh')
        btn.clicked.connect(generateMesh)
        grid.addWidget(btn, 1, 0)

        mw1 = MatplotlibWidget()
        grid.addWidget(mw1, 2, 0)
       
        tab3.setLayout(grid)
        
        
        #%% tab INVERSION SETTINGS
        tabInversionSettings = QWidget()
        tabs.addTab(tabInversionSettings, 'Inversion settings')
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
        
        
        tabInversionSettings.setLayout(grid)
        
        
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
                self.r2.createMesh()
        
            # write configuration file
            if self.r2.configFile == '':
                self.r2.write2in(param=param)
        
            self.r2.surveys[0].write2protocol(os.path.join(self.r2.dirname, 'protocol.dat'))
        
        
            self.process = QProcess(self)
            self.process.setWorkingDirectory(self.r2.dirname)
            # QProcess emits `readyRead` when there is data to be read
            self.process.readyRead.connect(dataReady)
            if OS == 'Linux':
                self.process.start('wine R2.exe')
            else:
                self.process.start(os.path.join(self.r2.dirname, 'R2.exe')) # need absolute path
            self.process.finished.connect(plotSection)
        
        def plotSection():
#            mwInvResult.plot(self.r2.showResults)
            mwInvResult.plot(self.r2.showSection)


        btn = QPushButton('Invert')
        btn.clicked.connect(logInversion)
        
        grid.addWidget(btn, 0, 0)
        
        logText = QTextEdit()
        grid.addWidget(logText, 1, 0)
        
        mwInvResult = MatplotlibWidget()
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
