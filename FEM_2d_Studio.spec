# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_all

datas = [('C:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\build\\bin\\Release\\fem_run.exe', 'bin'), ('C:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\example_continuous_beam.txt', '.'), ('C:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\example_multi_story_frame.txt', '.'), ('C:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\example_truss_bridge.txt', '.')]
binaries = []
hiddenimports = []
tmp_ret = collect_all('PyQt6')
datas += tmp_ret[0]; binaries += tmp_ret[1]; hiddenimports += tmp_ret[2]


a = Analysis(
    ['C:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\python\\fem_gui.py'],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='FEM_2d_Studio',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
