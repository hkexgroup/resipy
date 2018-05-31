#!/usr/bin/python3
# -*- coding: utf-8 -*-
import sys
import random

from PyQt5.QtWidgets import (QMainWindow, QApplication, QPushButton, QWidget, 
    QAction, QTabWidget,QVBoxLayout, QGridLayout, QLabel, QLineEdit, QMessageBox,
    QListWidget, QFileDialog, QCheckBox, QComboBox)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import QThread, pyqtSignal


from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import  numpy as np


from mesh_class import mesh_obj
import meshTools as mt # need to put in API and create modules


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
        print('plot is called')
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        callback(ax=ax)
        self.canvas.draw()
        

    def draw(self, fig):
        print('creating new figure')
        self.figure = fig
        self.figure.suptitle('hello')
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.draw()
        self.figure.canvas.draw()
    
    def drawAxes(self, ax):
        print('draw Awes')
        self.figure.clf()
        self.figure.axes.append(ax)
#        self.canvas = FigureCanvasQTAgg(self.figure)
        self.figure.canvas.draw()
        self.canvas.draw()
#    def update(self, fig):
        
        
        
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
        self.setGeometry(300,300,600,400)
        
        self.table_widget = QWidget()
        layout = QVBoxLayout()
        tabs = QTabWidget()
        tab1 = QWidget()
        tab2 = QWidget()
        tab3 = QWidget()
        
        #tab 1
        grid = QGridLayout()
        
        title = QLabel('Title')
        author = QLabel('Author')
        
        titleEdit = QLineEdit()
        titleEdit.setText('My beautiful survey')
        authorEdit = QLineEdit()
        
        grid.addWidget(title, 1, 0)
        grid.addWidget(titleEdit, 1, 1)
        
        grid.addWidget(author, 2, 0)
        grid.addWidget(authorEdit, 2, 1)
        
#        openFile = QAction(QIcon('open.png'), 'Open', self)
#        openFile.setShortcut('Ctrl+O')
#        openFile.setStatusTip('Open new File')
#        openFile.triggered.connect(self.showDialog)
#        
        
        button = QPushButton('Import Data')
        def getfile():
            self.fname, _ = QFileDialog.getOpenFileName(tab1,'Open File')
            authorEdit.setText(self.fname)
            
        button.clicked.connect(getfile)
        grid.addWidget(button, 3, 0)
        
        button = QPushButton('hello')
        button.clicked.connect(QMessageBox)
        grid.addWidget(button, 4, 0)
        
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
        
        
        tab1.setLayout(grid)
        
        
        
        
        # tab 2
        
        grid = QGridLayout()
#        self.fname = r'C:\Users\blanchy\Downloads\pyScripts\r2gui\f001_res.vtk'

                
        def callback2(ax):
            ax.plot(np.random.randn(20,5), '+--')
            ax.set_title('Random data nnnnndfghdfh')
        
        def updateGraph2():
            mw1.plot(callback2)
                
            
        def plotMesh():
            print(self.fname)
            mesh_dict=mt.vtk_import(self.fname)#makes a dictionary of a mesh 
            mesh = mesh_obj.mesh_dict2obj(mesh_dict)# this is a mesh_obj class instance 
    
#            mesh.summary()
            print('ready to draw')
            mw1.plot(mesh.show) # mesh.show() is the callback function to be called with ax
#            mw.plot(callback)
            
        btn = QPushButton('Plot Mesh')
        btn.clicked.connect(plotMesh)
        grid.addWidget(btn, 1, 0)

        mw1 = MatplotlibWidget()
        grid.addWidget(mw1, 2, 0)
       
        tab2.setLayout(grid)
        
        # tab 3
        
        grid = QGridLayout()
        
        singular_type = QCheckBox('Singularity Removal')
        grid.addWidget(singular_type, 0, 1)
        
        
        tab3.setLayout(grid)
        
        
        # tab 4
        
        tab4 = QWidget()

        grid = QGridLayout()

#        def drawSample():
#            mw.axis.clear()
#            mw.axis.plot(random.sample(range(0,10),10))
#            mw.canvas.draw()
        
        def callback(ax):
            ax.plot(np.random.randn(20,5), '*-')
            ax.set_title('Random data')
        
        def updateGraph():
            mw.plot(callback)
            
            
#        def newGraph():
#            fig, ax = plt.subplots()
#            ax.plot(np.random.randn(100,3))
#            mw.draw(fig)
##            mw.drawAxes(ax)
        

            
        btn = QPushButton('Draw random sample')
#        btn.clicked.connect(drawSample)
#        btn.clicked.connect(newGraph)
        btn.clicked.connect(updateGraph)
        mw = MatplotlibWidget()
        
        
        grid.addWidget(btn, 1, 0)
        grid.addWidget(mw, 2, 0)
        
        tab4.setLayout(grid)
                    
        
        
        tabs.addTab(tab1, 'Processing')
        tabs.addTab(tab2, 'Mesh')
        tabs.addTab(tab3, 'Inversion')
        tabs.addTab(tab4, 'Results')
        
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
