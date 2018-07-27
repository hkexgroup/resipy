# -*- mode: python -*-

block_cipher = None


a = Analysis(['splashScreen.py'],
             pathex=['C:\\Users\\blanchy\\Downloads\\r2gui-master'],
             binaries=[],
             datas=[('./dist/ui.zip','./pyR2.zip'),
                     ('./logo.png','./logo.png')],
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
          name='splashScreen',
          debug=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=True )
