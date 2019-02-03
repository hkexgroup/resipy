# -*- mode: python -*-

block_cipher = None


a = Analysis(['ui.py'],
             pathex=['/Users/jkl/Downloads/r2gui'],
             binaries=[],
             datas=[('./api/exe/R2.exe','./api/exe'),
                    ('./api/exe/gmsh.exe','./api/exe'),
                    ('./api/exe/cR2.exe', './api/exe'),
                    ('./api/exe/R3t.exe', './api/exe'),
                    ('./api/exe/cR3t.exe', './api/exe'),
                    ('./api/test/*','./api/test'),
                    ('./api/test/testTimelapse/*', './api/test/testTimelapse'),
                    ('./logo.png', '.'),
                    ('./r2icon.ico', '.')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='ui',
          debug=False,
          bootloader_ignore_signals=False,
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
app = BUNDLE(coll,
             name='pyR2.app',
             icon='icon.icns',
             bundle_identifier=None,
             info_plist={'NSHighResolutionCapable': 'True'})
