"""
Splash screen example

Eli Bendersky (eliben@gmail.com)
License: this code is in the public domain
Last modified: 09.05.2009
"""
from PyQt5.QtWidgets import QSplashScreen, QApplication, QProgressBar
from PyQt5.QtGui import QPixmap, QIcon
from PyQt5.QtCore import Qt
from zipfile import ZipFile, ZipInfo
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


""" PERMISSION ISSUE WITH ZIPFILE MODULE
https://stackoverflow.com/questions/39296101/python-zipfile-removes-execute-permissions-from-binaries
by default zipfile does not umpack any binary bit to say executable or not
Below is a way to do it
"""

class MyZipFile(ZipFile):
    def extract(self, member, path=None, pwd=None):
        if not isinstance(member, ZipInfo):
            member = self.getinfo(member)

        if path is None:
            path = os.getcwd()

        ret_val = self._extract_member(member, path, pwd)
        attr = member.external_attr >> 16
        os.chmod(ret_val, attr)
        return ret_val


if __name__ == "__main__":

    app = QApplication(sys.argv)
#    app.setWindowIcon(QIcon(os.path.join(bundle_dir, 'logo.png')))

    splash_pix = QPixmap(os.path.join(bundle_dir, 'logo.png'))
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
    splash.showMessage("Expanding app", Qt.AlignBottom, Qt.white)
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

    zf = MyZipFile(os.path.join(bundle_dir, 'ui.zip'),'r')
    extractDir = os.path.join(bundle_dir, 'ui')
    if os.path.exists(extractDir):
#        print('overwritting pyR2 dir')
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
    print('finished unzipping')
    
    splash.hide()
    splash.close()
    appDir = os.path.join(bundle_dir, 'ui', 'ui') # zip always putting a double dir ... don't know why
    print('Main app will be run in appDir = ', appDir)
    os.chdir(appDir)
#    os.system(['python3', 'ui.py']) # this work fine
    os.system(os.path.join(appDir, 'pyR2')) # this works now as well !

#  need to comment the following lines as the exit signal is given by the main app
#    sys.exit()
#    sys.exit(app.exec_()) 
    
""" NOTE
This approach increase significantly the size of the package from 150 to 210 MB
Another approach would be to load all modules in this script and just unzip the
sources and run `python ui.py`.
"""

