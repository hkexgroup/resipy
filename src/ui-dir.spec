# -*- mode: python -*-
# this file is used by pyinstaller to generate a zip file that would 
# then be uncompressed by the splashScreen.spec
# nb needs to be a python 3.7 environment for the final executable to work 
import platform, os, sys
sys.setrecursionlimit(5000)
block_cipher = None

OS = platform.system()

# https://stackoverflow.com/questions/11322538/including-a-directory-using-pyinstaller  
datas=[('./resipy/exe/*','./resipy/exe'),
        ('./logo.png', '.'),
        ('./loadingLogo.png', '.'),
        ('./image/dipdip.png', './image'),
        ('./image/schlum.png', './image'),
        ('./image/wenner.png', './image'),
        ('./image/gradient.png', './image'),
        ]

if OS == 'Linux':
    datas.append(('./resipy/cext/*gnu.so','./resipy/cext'))
elif OS == 'Windows':
    datas.append(('./resipy/cext/*.pyd','./resipy/cext'))
#raise error if not linux or windows? 

def extra_datas(mydir):
    def rec_glob(p, files):
        import os
        import glob
        for d in glob.glob(p):
            if os.path.isfile(d):
                files.append(d)
            rec_glob("%s/*" % d, files)
    files = []
    rec_glob("%s/*" % mydir, files)
    extra_datas = []
    for f in files:
        extra_datas.append((f, os.path.dirname(os.path.join('resipy',f))))

    return extra_datas

datas += extra_datas('examples')


a = Analysis(['ui.py'],
             pathex=[],
             binaries=[],
             datas = datas,
             hiddenimports=['pkg_resources.py2_warn','vtkmodules.all','vtk.util.numpy_support','vtkmodules.util.numpy_support',
			 'vtkmodules.qt.QVTKRenderWindowInteractor','scipy._lib.array_api_compat.numpy.fft','scipy.special._special_ufuncs'],
             hookspath=[],
			 hooksconfig={"matplotlib": {"backends": "all"}},
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
           
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='ResIPy',
          debug=False,
          strip=False,
          upx=True,
          console=False,
          version='Version.details',
          icon='logo.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='ui')
