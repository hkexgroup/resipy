#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 15:14:25 2019

@author: jkl
"""
# to run use: python3 -m pytest test-ui.py


from PyQt5 import QtCore

from ui import App
import time


testdir = 'examples/'
sleepTimeDefault = 1 # s

#%%

def test_dc_2d_syscal(qtbot, qapp):
    app = App() 
    qtbot.addWidget(app)
    def proc(sleepTime=None):
        if sleepTime is None:
            sleepTime = sleepTimeDefault
        qapp.processEvents()
        time.sleep(sleepTime)
        
    # importing
    app.restartFunc()
    app.ftype = 'Syscal'
    app.importFile('examples/dc-2d/syscal.csv')
    proc()
    app.psContourCheck.setChecked(True)
    app.pvmin.setText('10')
    app.pvmax.setText('50')
    proc()
    
    app.tabs.setCurrentIndex(1)
    proc()
    
    app.tabs.setCurrentIndex(2)
    proc()
    
    app.tabs.setCurrentIndex(4)
    proc()
    
    app.tabs.setCurrentIndex(5)
    qtbot.mouseClick(app.invertBtn, QtCore.Qt.LeftButton)
    proc(5)
    
    # app.vminEdit.setText('0')
    # app.vmaxEdit.setText('50')
    # qtbot.mouseClick(app.keepApplyBtn, QtCore.Qt.LeftButton)
    # qtbot.keyClick(app.coilCombo, QtCore.Qt.Key_Up)
    # qtbot.keyClick(app.coilCombo, QtCore.Qt.Key_Up)
    # qtbot.keyClick(app.coilCombo, QtCore.Qt.Key_Up)
    # qtbot.keyClick(app.coilCombo, QtCore.Qt.Key_Up)
    # qtbot.keyClick(app.coilCombo, QtCore.Qt.Key_Up)
    # qtbot.mouseClick(app.rollingBtn, QtCore.Qt.LeftButton)
    # qtbot.mouseClick(app.mapRadio, QtCore.Qt.LeftButton)
    # qtbot.keyClick(app.coilCombo, QtCore.Qt.Key_Down)
    # qtbot.keyClick(app.coilCombo, QtCore.Qt.Key_Down)
    # qtbot.keyClick(app.coilCombo, QtCore.Qt.Key_Down)
    # qtbot.keyClick(app.coilCombo, QtCore.Qt.Key_Down)
    # qtbot.keyClick(app.coilCombo, QtCore.Qt.Key_Down)
    # qtbot.mouseClick(app.contourCheck, QtCore.Qt.LeftButton)
    # qtbot.mouseClick(app.ptsCheck, QtCore.Qt.LeftButton)
    # app.vminfEdit.setText('0')
    # app.vmaxfEdit.setText('30')
    # qtbot.mouseClick(app.applyBtn, QtCore.Qt.LeftButton)
    # proc()

    # # calibration
    # app.tabs.setCurrentIndex(2)
    # app.fnameEC = 'examples/calib/dfec.csv'
    # app.fnameECa = 'examples/calib/dfeca.csv'
    # qtbot.mouseClick(app.fitCalibBtn, QtCore.Qt.LeftButton)
    # qtbot.mouseClick(app.applyCalibBtn, QtCore.Qt.LeftButton)
    # proc()
    
    # # inversion setting
    # app.tabs.setCurrentIndex(4)
    # app.nLayerEdit.setText('4')
    # app.thicknessEdit.setText('0.3')
    # app.startingEdit.setText('30')
    # qtbot.mouseClick(app.createModelBtn, QtCore.Qt.LeftButton)
    # qtbot.mouseClick(app.lcurveBtn, QtCore.Qt.LeftButton)
    # proc()
    
    # # inversion
    # app.tabs.setCurrentIndex(5)
    # qtbot.keyClick(app.forwardCombo, QtCore.Qt.Key_Down)
    # qtbot.keyClick(app.methodCombo, QtCore.Qt.Key_Down)
    # qtbot.mouseClick(app.invertBtn, QtCore.Qt.LeftButton)
    # proc(2)
    
    # # misfit
    # app.tabs.setCurrentIndex(6)
    # proc(2)
    
    # # about
    # app.tabs.setCurrentIndex(7)
    # proc()
    

