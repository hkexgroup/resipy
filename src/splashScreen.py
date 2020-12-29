""" splashScreen

"""
from PyQt5.QtWidgets import QSplashScreen, QApplication, QProgressBar
from PyQt5.QtGui import QPixmap, QIcon, QMovie
from PyQt5.QtCore import Qt
from zipfile import ZipFile, ZipInfo
from subprocess import Popen, PIPE
import os, sys, shutil, platform, time 
from os.path import expanduser
QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True) # for high dpi display
QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

ctime = time.time()
OS = platform.system()     
if OS == 'Windows':
    import winreg      

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

#workaround to deal with removing old _MEI folders on windows (works when compile WITH console=False)
try:
    if OS == 'Windows':
        week_seconds = 7*24*60*60
        cmd = "echo %LOCALAPPDATA%" # get local storage 
        p = Popen(cmd,shell=True,stdout=PIPE)
        appdata = p.stdout.readline().rstrip().decode()
        temp_path = os.path.join(appdata,'Temp')
        files = sorted(os.listdir(temp_path))
        print('Checking for old _MEI directories in %s'%temp_path)
        for f in files:
            mtime = os.stat(os.path.join(temp_path,f)).st_mtime#modifcation time
            dtime = ctime - mtime # delta time in seconds 
            if f.find('_MEI')==0 and dtime>week_seconds:
                print('removing %s ...'%f,end='')
                try:
                    cmd = "RMDIR {:s} /q /s".format(os.path.join(temp_path,f))
                    p = Popen(cmd,shell=True)
                    print('done.')
                except:
                    print('ERROR')
except:
    pass


""" PERMISSION ISSUE WITH ZIPFILE MODULE
https://stackoverflow.com/questions/39296101/python-zipfile-removes-execute-permissions-from-binaries
by default zipfile does not umpack any binary bit to say executable or not
Below is a way to do it
"""
# Windows/Linux bundle workaround for restarting the app and saving splash settings
class SplashSettings(object):
    def __init__(self):
        """For help see resipy.Settings()"""
        self.param = {} # Dict of settings
    def setSetting(self):
        if self.param != {}:
            return '\n'.join('{} = {}'.format(str(key), str(value)) for key, value in self.param.items())
    def readSetting(self, settings):
        settingList = settings.split('\n')
        for val in settingList:
            key = val.split(' = ')[0]
            value = val.split(' = ')[-1]
            self.param[key] = value
    # for windows only
    def setReg(self, name, value, REG_PATH=r'SOFTWARE\ResIPy'):
        try:
            winreg.CreateKeyEx(winreg.HKEY_CURRENT_USER, REG_PATH, 0)
            registry_key = winreg.OpenKey(winreg.HKEY_CURRENT_USER, REG_PATH, 0, 
                                           winreg.KEY_WRITE)
            winreg.SetValueEx(registry_key, name, 0, winreg.REG_SZ, value)
            winreg.CloseKey(registry_key)
        except WindowsError:
            return None
    def getReg(self, name, REG_PATH=r'SOFTWARE\ResIPy'):
        try:
            registry_key = winreg.OpenKey(winreg.HKEY_CURRENT_USER, REG_PATH, 0)
            value, regtype = winreg.QueryValueEx(registry_key, name)
            winreg.CloseKey(registry_key)
            return value
        except WindowsError:
            return None     
    # for all OS platforms
    def genLocalSetting(self):
        try:
            settingRaw = self.setSetting()
            if OS == 'Windows':
                self.setReg('settings', settingRaw)             
            elif OS == 'Darwin':
                userDir = expanduser('~/Library/Preferences/')
                with open(os.path.join(userDir,'com.resipy.plist'), 'w') as settings_file:
                    settings_file.write(settingRaw)
            elif OS == 'Linux':
                userDir = expanduser('~/.local/share/')
                with open(os.path.join(userDir,'resipy.stt'), 'w') as settings_file:
                    settings_file.write(settingRaw)            
        except:
            pass
    def retLocalSetting(self):
        try:
            if OS == 'Windows':
                settingRaw = self.getReg('settings')
            elif OS == 'Darwin':
                userDir = expanduser('~/Library/Preferences/')
                with open(os.path.join(userDir,'com.resipy.plist'), 'r') as settings_file:
                    settingRaw = settings_file.read()
            elif OS == 'Linux':
                userDir = expanduser('~/.local/share/')
                with open(os.path.join(userDir,'resipy.stt'), 'r') as settings_file:
                    settingRaw = settings_file.read()
            self.readSetting(settingRaw)
            return True
        except:
            return None


if frozen == 'ever so': # only for frozen packages we need this
    settings = SplashSettings()
    settings.retLocalSetting() # in case there are settings in there already
    settings.param['exe_path'] = os.path.abspath(sys.executable)
    settings.param['frozen'] = True
    settings.genLocalSetting()

class MySplashScreen(QSplashScreen):
    def __init__(self, animation, flags):
        # run event dispatching in another thread
        QSplashScreen.__init__(self, QPixmap(), flags)
        self.movie = QMovie(animation)
        self.movie.frameChanged.connect(self.onNextFrame)
        #self.connect(self.movie, SIGNAL('frameChanged(int)'), SLOT('onNextFrame()'))
        self.movie.start()

    def onNextFrame(self):
        pixmap = self.movie.currentPixmap()
        self.setPixmap(pixmap)
        self.setMask(pixmap.mask())
        
        
class MyZipFile(ZipFile):
    def extract(self, member, path=None, pwd=None):
        if not isinstance(member, ZipInfo):
            member = self.getinfo(member)

        if path is None:
            path = os.getcwd()

        ret_val = self._extract_member(member, path, pwd)
        attr = member.external_attr >> 16
        if OS != 'Windows':
            os.chmod(ret_val, attr) # IMPORTANT this line needs to be commented otherwise we got the _MEIxxxxx issue
        # changing the permission of somes files makes them unremovable by the splascreen bootloader when the program finished
        # this leads to accumulation of _MEIxxxxx temporary files in C:\Users\User\AppData\Local\Temp\
        # this issue is windows specific, on Linux, the temporary folder in /tmp is removed even when we uncomment this line
        return ret_val


if __name__ == "__main__":

    app = QApplication(sys.argv)
    app.setStyle('Fusion')
#    app.setWindowIcon(QIcon(os.path.join(bundle_dir, 'logo.png')))

    splash_pix = QPixmap(os.path.join(bundle_dir, 'loadingLogo.png'))
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
#    splash = MySplashScreen('chicken.gif', Qt.WindowStaysOnTopHint)
    splash.setWindowFlags(Qt.WindowStaysOnTopHint | Qt.FramelessWindowHint | Qt.Tool)
    splash.setEnabled(False)
    # splash = QSplashScreen(splash_pix)
    # adding progress bar
    progressBar = QProgressBar(splash)
#    progressBar.setMaximum(10)
    progressBar.setGeometry(100, splash_pix.height() - 50, splash_pix.width() - 200, 20)
#    progressBar.setGeometry(150, 320, 200, 18)
    # splash.setMask(splash_pix.mask())

    splash.show()
    splash.showMessage("Expanding app", Qt.AlignBottom | Qt.AlignCenter, Qt.black)
    app.processEvents()
    
#    initLoop = Qt.QEventLoop()
#    pool = Pool(processes=1)
#    pool.apply_async(longInitialization, [2], callback=lambda exitCode: initLoop.exit(exitCode))
#    initLoop.exec_()
    
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
#        print('overwritting ResIPy dir')
        shutil.rmtree(extractDir)
    os.mkdir(extractDir)
    uncompress_size = sum((file.file_size for file in zf.infolist()))
    extracted_size = 0

    for file in zf.infolist():
        extracted_size += file.file_size
        percentage = extracted_size/uncompress_size*100
        progressBar.setValue(percentage)
        if percentage > 50 and percentage < 70:
            splash.showMessage("Copying temp files", Qt.AlignBottom | Qt.AlignCenter, Qt.black)
        if percentage >= 70 and percentage < 80:
            splash.showMessage("Checking files", Qt.AlignBottom | Qt.AlignCenter, Qt.black)
        if percentage >= 80 and percentage < 90:
            splash.showMessage("Loading PyQt", Qt.AlignBottom | Qt.AlignCenter, Qt.black)
        if percentage >= 90 and percentage < 98:
            splash.showMessage("Loading App", Qt.AlignBottom | Qt.AlignCenter, Qt.black)
        if percentage >= 98:
            splash.showMessage("Almost there!", Qt.AlignBottom | Qt.AlignCenter, Qt.black)
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
    if OS == 'Linux':
        os.system(os.path.join(appDir, 'ResIPy'))
    else:
        p = Popen(os.path.join(appDir, 'ResIPy.exe'), shell=False, stdout=None, stdin=None) # this works now as well !
        p.communicate() # block and wait for the main program to finish
    # this last one doesn't work on linux WHEN COMPILED and I don't know why
    
    
    # clearing splash screen settings in case file is relocated/updated/etc.
    settings.retLocalSetting() # in case there are settings in there already
    del settings.param['exe_path']
    del settings.param['frozen']
    settings.genLocalSetting()
    
    print('splashScreen is exiting')
    sys.exit(0) # send the SIGTERM signal -> works
#    sys.exit(app.exec_()) # doesn't work
    
""" NOTE
This approach increase significantly the size of the package from 150 to 210 MB
Another approach would be to load all modules in this script and just unzip the
sources and run `python ui.py`.
"""
