"""
Splash screen example

Eli Bendersky (eliben@gmail.com)
License: this code is in the public domain
Last modified: 09.05.2009
"""
from PyQt5.QtWidgets import QSplashScreen, QApplication
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import *
import zipfile
import os
import sys


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
    splash.show()

    app.processEvents()
    
    zip_ref = zipfile.ZipFile(os.path.join(bundle_dir, 'pyR2.zip'),'r')
#    os.mkdir(os.path.join(bundle_dir, 'pyR2'))
    zip_ref.extractall(os.path.join(bundle_dir,'pyR2'))
    zip_ref.close()
    
    splash.finish()
#    os.popen(os.path.join(bundle_dir, 'pyR2', 'ui.exe'))
    
    app.exec_()