# -*- mode: python -*-

block_cipher = None


datas=[('./resipy/exe/*','./resipy/exe'),
       ('./resipy/cext/*', './resipy/cext'),
       ('./logo.png', '.'),
       ('./logo.ico', '.'),
       ('./loadingLogo.png', '.'),
       ('./image/dipdip.png', './image'),
       ('./image/schlum.png', './image'),
       ('./image/wenner.png', './image'),
       ('./image/gradient.png', './image')]


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
             datas=datas,
             hiddenimports=['pkg_resources.py2_warn', 'vtkmodules', 'vtkmodules.all', 
             'vtkmodules.util.numpy_support', 'vtkmodules.numpy_interface', 'vtkmodules.numpy_interface.dataset_adapter'],
             hookspath=[],
			 hooksconfig={"matplotlib": {"backends": "all"}},
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
             name='ResIPy.app',
             icon='logo.icns',
             bundle_identifier=None,
             info_plist={'NSHighResolutionCapable': 'True'})
