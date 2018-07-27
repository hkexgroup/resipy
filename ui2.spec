# -*- mode: python -*-

block_cipher = None


a = Analysis(['ui.py'],
             pathex=['C:\\Users\\blanchy\\Downloads\\r2gui-master'],
             binaries=[],
             datas=[('./api/exe/R2.exe','./api/exe'),
					('./api/exe/gmsh.exe','./api/exe'),
					('./api/exe/cR2.exe', './api/exe')],
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
          exclude_binaries=True,
          name='ui',
          debug=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='ui')
