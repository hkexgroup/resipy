# -*- mode: python -*-
# this file is used by pyinstaller to generate a zip file that would 
# then be uncompressed by the splashScreen.spec

block_cipher = None

# https://stackoverflow.com/questions/11322538/including-a-directory-using-pyinstaller  
datas=[('./resipy/exe/R2.exe','./resipy/exe'),
        ('./resipy/exe/gmsh.exe','./resipy/exe'),
        ('./resipy/exe/cR2.exe', './resipy/exe'),
        ('./resipy/exe/R3t.exe', './resipy/exe'),
        ('./resipy/exe/cR3t.exe', './resipy/exe'),
        ('./logo.png', '.'),
        ('./loadingLogo.png', '.'),
        ('./image/dipdip.png', './image'),
        ('./image/schlum.png', './image'),
        ('./image/wenner.png', './image'),
        ('./image/gradient.png', './image')
        ]

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
        extra_datas.append((f, os.path.dirname(os.path.join('resipy', 'invdir',f))))

    return extra_datas

datas += extra_datas('examples')


a = Analysis(['ui.py'],
             pathex=[],
             binaries=[],
             datas = datas,
             hiddenimports=['numpy.core._dtype_ctypes'],
             hookspath=[],
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
          console=False)
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='ui')
