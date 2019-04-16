# -*- mode: python -*-
# this file is used by pyinstaller to generate a zip file that would 
# then be uncompressed by the splashScreen.spec

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
					('./resipy/test/IP/*','./resipy/test/IP'),
                    ('./resipy/test/testTimelapse/*', './resipy/test/testTimelapse'),
                    ('./logo.png', '.'),
                    ('./loadingLogo.png', '.'),
                    ('./image/dipdip.png', './image'),
                    ('./image/schlum.png', './image'),
                    ('./image/wenner.png', './image'),
                    ('./image/gradient.png', './image')
                    ],
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
          console=False )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='ui')
