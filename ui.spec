## -*- mode: python -*-
import sys
sys.setrecursionlimit(1500)
block_cipher = None


a = Analysis(['ui.py'],
             pathex=[],
             binaries=[],
             datas=[('./api/exe/R2.exe','./api/exe'),
                    ('./api/exe/gmsh.exe','./api/exe'),
                    ('./api/exe/cR2.exe', './api/exe'),
                    ('./logo.png', '.'),
                    ('./r2icon.ico', '.')],
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
          name='pyR2',
          debug=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=True,
		  icon='r2icon.ico')

