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
                    ('./resipy/test/*','./resipy/test'),
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
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='ResIPy',
          debug=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=False,
		  icon='logo.ico')

