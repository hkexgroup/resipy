#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 29 15:38:17 2018

@author: jkl
"""

import sys
from PyQt5.QtWidgets import QApplication, QWidget, QTabWidget, QPushButton, \
    QMessageBox

if __name__ == '__main__':

    app = QApplication(sys.argv)

    w = QWidget()
    w.resize(250, 150)
    w.move(300, 300)
    w.setWindowTitle('Example windows')
#
#tab = QTabWidget(w)
#tab1 = tab.addTab(w,'Importing data')
#tab2 = tab.addTab(w,'Processing')
#tab3 = tab.addTab(w,'Inversion')

    submit = QPushButton('hello',w)
    bquit = QPushButton('quit app',w)

#bquit.clicked.connect(app.quit)

    w.show()

    app.exec()   
    sys.exit(app.exec_())

