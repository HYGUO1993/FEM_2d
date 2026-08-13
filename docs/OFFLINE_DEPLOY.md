# FEM_2d 线下部署指南（离线 / 本地）

本指南说明如何**不依赖网络**在本地把 FEM_2d 编译、运行并启动 GUI。已在本机（Windows 11 + MSVC 19.44）实测通过。

---

## A. 求解内核（C++）编译

项目只需一个源文件 `barsystem.cpp`，无第三方 C/C++ 依赖。两种本地编译方式任选其一。

### 方式 1：CMake（推荐，跨平台，产出路径与 GUI 一致）

前置：安装 **Visual Studio 2022（含“使用 C++ 的桌面开发”）** 或 **MinGW-w64（g++）**。

```bash
# 在仓库根目录
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
# 产物：build/bin/Release/fem_run.exe  （GUI 直接识别此路径）

# 运行单元测试
ctest --test-dir build --build-config Release
# 或直接：build/bin/Release/unit_tests.exe
```

> Linux/macOS 把生成器换成 `Unix Makefiles` / `Xcode`，产物为 `build/bin/Release/fem_run`。

### 方式 2：直接用 `cl.exe`（无需打开 VS 命令行，本机实测）

仓库已附 `build_local.bat`，它绕过 `cmd` 限制直接用 `cl.exe` 并手动设置 `INCLUDE/LIB`。若你机器上 VS 路径不同，按下面环境变量改 `build_local.bat` 即可：

```bat
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
cl /EHsc /utf-8 /O2 barsystem.cpp /Fobuild\obj\barsystem.obj /Febuild\bin\fem_run.exe
cl /EHsc /utf-8 /O2 /DUNIT_TEST tests\test_basic.cpp barsystem.cpp /Fobuild\obj\ /Febuild\bin\unit_tests.exe
```

> 注：直接 `cl` 产出在 `build/bin/`，而 GUI 期望 `build/bin/Release/`。请额外复制一份：
> `mkdir build\bin\Release && copy build\bin\fem_run.exe build\bin\Release\`

---

## B. 运行求解（命令行）

```bash
# 悬臂梁（frame）
build/bin/Release/fem_run.exe --no-stiff --input test05.txt --output build/Results_frame.dat

# 简支梁（三点载荷）
build/bin/Release/fem_run.exe --no-stiff --input test_beam.txt --output build/Results_beam.dat

# 可视化（需要 Python + matplotlib，见下）
python scripts/visualize_results.py --results build/Results_frame.dat --scale 1000 --out build/plot_frame.png
```

输出 `Results_*.dat` 含三段：`Node Displacements` / `Element End Forces` / `Support Reactions`。

> ✅ 自 2026-08-12 起，`Element End Forces` 已改为 **LOCAL 坐标**（`N=轴向力, V=剪力, M=弯矩`），表头为 `N_i V_i M_i N_j V_j M_j`，倾斜杆件可直接读轴力，无需投影。均布荷载（`LATERAL_UNIFORM_PRESSURE`）等单元荷载已接入固端力法。

> ⚠️ 约束文件语义：**`0`=自由 DOF、`-1`=约束 DOF**（与常见约定相反）。简支梁两端应写 `-1 -1 0`（X、Y 约束，R 自由）；若误写 `0 -1 0` 会形成轴向刚体模态，程序会报「刚度矩阵奇异或病态」而非输出假解。

---

## C. GUI（Python + PyQt6）部署

GUI 通过子进程调用上面的 `fem_run.exe` 真正求解并可视化。

### 1. 安装依赖（一次性，需网络一次；之后离线可用）

```bash
pip install PyQt6 matplotlib numpy qtawesome
```

> 若本机无法联网装包，可在一台能联网的机器上用 `pip download` 把这几个包下载成离线 whl，再 `pip install *.whl` 离线安装。

### 2. 启动 GUI

```bash
python quick_start_gui.py      # 推荐：会自动检查依赖、定位 fem_run.exe
# 或
python python/fem_gui.py        # 需从仓库根目录运行（路径已处理）
```

界面：左侧模型树/编辑器/视图控制，右侧画板（变形图、轴力色图、支座反力箭头），F5 求解。

> 启动前确保 `build/bin/Release/fem_run.exe` 已存在（见 A），否则 GUI 会提示“找不到求解器”。

---

## D. 本机已交付的产物

| 文件 | 说明 |
|---|---|
| `build/bin/fem_run.exe` | MSVC 19.44 编译的求解器（已验证与解析解吻合） |
| `build/bin/unit_tests.exe` | 单元测试（运行输出 `all tests passed`） |
| `build/bin/Release/fem_run.exe` | GUI 期望路径，已就位 |
| `build_local.bat` | 离线编译脚本 |
| `verify/verify_solver.py` | 独立 Python 参考求解器，用于交叉验证内核正确性 |
| `test_inclined.txt` | 45° 桁架复现用例（证明端力坐标系问题） |

---

## E. 常见问题

- **GUI 提示找不到求解器**：把 `fem_run.exe` 放到 `build/bin/Release/`（方式 1 自动满足；方式 2 需手动复制）。
- **DLL load failed（PyQt6）**：通常是 Python 与 PyQt6 位数不一致（都需 64 位），重装对应版本即可。
- **求解结果位移异常大 / NaN**：多为约束不足（机构）或误用 TRUSS 分析梁问题（梁需 FRAME 且 `Iz>0`）。✅ 自 2026-08-12 起奇异/病态会明确报错 `ERROR: 刚度矩阵奇异或病态` 并返回退出码 -2，不再掩盖。
- **`python python/fem_gui.py` 报 `No module named gui`**：改用 `quick_start_gui.py` 启动（改进报告 §11）。
