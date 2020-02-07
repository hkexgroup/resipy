## -*- mode: python -*-
import sys
sys.setrecursionlimit(1500)
block_cipher = None


a = Analysis(['ui.py'],
             pathex=[],
             binaries=[],
             datas=[('./resipy/exe/R2.exe','./resipy/exe'),
                    ('./resipy/exe/gmsh.exe','./resipy/exe'),
                    ('./resipy/exe/cR2.exe', './resipy/exe'),
                    ('./resipy/exe/R3t.exe', './resipy/exe'),
                    ('./resipy/exe/cR3t.exe', './resipy/exe'),
                    ('./resipy/test/testTimelapse/*', './resipy/test/testTimelapse'),
                    ('./logo.png', '.'),
                    ('./logo.ico', '.')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)

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


pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='ResIPy',
          debug=True,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=False,
		  icon='logo.ico')

