"""
Splash screen example

Eli Bendersky (eliben@gmail.com)
License: this code is in the public domain
Last modified: 09.05.2009
"""
from PyQt5.QtWidgets import QSplashScreen, QApplication, QProgressBar
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import *
import zipfile
from subprocess import Popen
import os
import sys
import time
import shutil

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



if __name__ == "__main__":

    app = QApplication(sys.argv)

    splash_pix = QPixmap('logo.png')
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setWindowFlags(Qt.WindowStaysOnTopHint | Qt.FramelessWindowHint)
    splash.setEnabled(False)
    # splash = QSplashScreen(splash_pix)
    # adding progress bar
    progressBar = QProgressBar(splash)
    progressBar.setMaximum(10)
    progressBar.setGeometry(0, splash_pix.height() - 50, splash_pix.width(), 20)

    # splash.setMask(splash_pix.mask())

    splash.show()
    splash.showMessage("Uncompressin app", Qt.AlignBottom, Qt.white)
    app.processEvents()
#    app.processEvents()    
#    progressBar.setValue(1)
#    app.processEvents()
#    app.processEvents()
#    app.processEvents()
#    time.sleep(4)
##    progressBar.setValue(8)
#    app.processEvents()
#    time.sleep(4) # here insert the code to unzip the file
#    progressBar.setValue(10)
#    app.processEvents()
    
#    for i in range(1, 11):
#        progressBar.setValue(i)
#        t = time.time()
#        app.processEvents()
#    t = time.time()
#    while time.time() < t + 4:
#       app.processEvents()

    zf = zipfile.ZipFile(os.path.join(bundle_dir, 'pyR2.zip'),'r')
    extractDir = os.path.join(bundle_dir, 'pyR2')
    if os.path.exists(extractDir):
        shutil.rmtree(extractDir)
    os.mkdir(extractDir)
    uncompress_size = sum((file.file_size for file in zf.infolist()))
    extracted_size = 0

    for file in zf.infolist():
        extracted_size += file.file_size
        percentage = extracted_size * 100/uncompress_size
        progressBar.setValue(percentage/10)
        app.processEvents()
        zf.extract(file, extractDir)
    zf.close()

    
    splash.hide()
#    os.chdir(os.path.join(bundle_dir, 'pyR2', 'ui'))
#    Popen(['python3', 'ui.py']) # or the exe file is compiled 
    print(os.path.join(bundle_dir, 'pyR2', 'ui', 'ui'))
    Popen(os.path.join(bundle_dir, 'pyR2', 'ui', 'ui')) # permission denied because the zip
    # doesn't keep the permission

#    sys.exit()
#    sys.exit(app.exec_()) # TODO solve that !