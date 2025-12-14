# FEM_2d 项目说明
本项目是一个用于教学/实验的二维杆/桁架/刚架有限元程序，支持读取结构与载荷、组装整体刚度、求解位移、输出杆端力与支座反力，并提供 Python 可视化脚本生成变形图与反力箭头。

- 语言与平台：C++（C++11），Windows（VS 2022/MSVC），Python 3.12 可视化
- 可执行文件： build\\bin\\Debug\\fem_run.exe
- 单元测试： build\\bin\\Debug\\unit_tests.exe
- 可视化脚本： scripts\\visualize_results.py
快速开始

- 构建与测试（推荐 VS 2022/MSBuild）：
  - & 'C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\MSBuild\\Current\\Bin\\MSBuild.exe' .\\build\\FEM_2d.sln /p:Configuration=Debug /m
  - & .\\build\\bin\\Debug\\unit_tests.exe
- 运行示例：
  - 框架悬臂示例（梁）： & .\\build\\bin\\Debug\\fem_run.exe --no-stiff --input test05.txt --output build\\Results_frame.dat
  - 简支梁三点载荷： & .\\build\\bin\\Debug\\fem_run.exe --no-stiff --input test_beam.txt --output build\\Results_beam.dat
- 可视化：
  - python scripts\\visualize_results.py --results build\\Results_frame.dat --scale 1000 --out build\\plot_frame.png
  - python scripts\\visualize_results.py --results build\\Results_beam.dat --scale 2000 --out build\\plot_beam.png
目录结构

- barsystem.cpp 主程序与核心有限元函数
- barsystem.h 数据结构与函数声明
- tests\\test_basic.cpp 轻量单元测试
- CMakeLists.txt 构建配置（Windows 生成 VS 解决方案）
- scripts\\visualize_results.py Python 可视化脚本
- test05.txt 、 test_beam.txt 示例输入
- build\\ 构建产物（可执行、结果与图像）
程序入口与主流程

- 入口： barsystem.cpp:12-39 解析命令行参数 --input/--output/--quiet/--no-stiff/--stiff
- 核心计算流程（关键调用位置： barsystem.cpp:280-339 ）
  - 读取全局与实体数据（节点/约束/元素/材料/截面/载荷）
  - 计算单元长度与方向余弦： barsystem.cpp:400-416
  - 自由度编号与带宽地址： barsystem.cpp:767-835 、 barsystem.cpp:872-918
  - 组装整体刚度矩阵 Skyline 下三角： barsystem.cpp:523-614
  - 组装载荷向量（目前支持节点集中力 FORCE_ON_NODE ）： barsystem.cpp:1035-1050
  - LDLT 求解自由度位移： barsystem.cpp:631-680
  - 杆端内力计算（按单元刚度乘以局部位移）： barsystem.cpp:1052-1101
  - 支座反力： barsystem.cpp:1103-1113 （约束 DOF 处的反力回写到 pLoadVect ）
  - 结果输出（位移/内力/反力）： barsystem.cpp:987-1033
数据格式（输入文件） 第 1 行为总控数据：

- TotalNodes ConstrainedNodes TotalElements MaterialTypes SectionTypes LoadCount
- 示例： 2 1 1 1 1 1 表示 2 节点、1 约束节点、1 单元、1 材料、1 截面、1 载荷
后续块：

- 节点（每行）： NodeType X Y ， NodeType=1 桁架节点、 NodeType=2 刚架节点
- 约束（每行）： Node X Y R ，值 -1 表示该 DOF 受约束， 0 表示自由
- 元素（每行）： ElemType StartNode EndNode SectionIndex [MaterialIndex可选]
  - ElemType=1 桁架单元， ElemType=2 刚架单元
- 材料（每行）： E Mu Alpha
- 截面（每行）： A Iz H （梁问题必须有 Iz>0 ）
- 载荷（每行）： Type Dir Value LoadedElem LoadedNode Position T0 T1
  - 节点集中力： Type=1 ，方向 Dir=0(X)/1(Y)/2(R) ，示例： 1 1 -1000 -1 1 0 0 0 表示在节点 1 施加竖向 -1000N

数据字段详解

- 节点类型 `NodeType`：
  - `1`（TRUSS）：节点仅使用 `Ux, Uy` 两个 DOF，旋转 `Rz` 不参与
  - `2`（FRAME）：节点使用 `Ux, Uy, Rz` 三个 DOF（平面刚架）
- 约束行 `Node X Y R`：
  - `-1` 表示该自由度被约束（位移或转角为 0），`0` 表示自由
  - 约束节点编号为 `Node`，对应节点的 `Ux(U)`、`Uy(V)`、`Rz(θ)` 三个分量
- 元素行：
  - TRUSS 单元使用 4 自由度局部向量 `[Ux_i, Uy_i, Ux_j, Uy_j]`
  - FRAME 单元使用 6 自由度局部向量 `[Ux_i, Uy_i, Rz_i, Ux_j, Uy_j, Rz_j]`
- 材料与截面：
  - 材料 `E`（弹性模量）、`Mu`（泊松比）、`Alpha`（热膨胀系数，当前未用于计算）
  - 截面 `A`（面积）、`Iz`（关于 Z 的惯性矩，FRAME 弯曲刚度依赖 `Iz`）、`H`（辅助几何参数，当前不参与刚度）

坐标系与自由度

- 节点自由度编号：`Node.iaDOFIndex[0..2]` 分别对应 `Ux, Uy, Rz`（见位移输出： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:987-1001）
- 单元局部坐标系沿杆轴方向定义，长度与方向余弦由 `LengthSinCosCalcu` 计算（ c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:540 ）
- 元端力输出顺序：`Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j`（ c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:1003-1016 ）
- 桁架（TRUSS）不输出端弯矩（`Mz_i, Mz_j = 0`），刚架（FRAME）输出三向端力与端弯矩

Skyline 存储与访问

- 总刚度矩阵以 Skyline 下三角形式存储与装配（ c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:523-614 ）
- 单元刚度装配函数：`GKAssembly`（ c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:528-619 ）
- 访问函数：`GetElementInGK` 使用索引公式 `idx = pDiag[iRow] - iRow + iCol`（下三角）读取条目（ c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:621-626 ）

示例文件：

- 悬臂梁（框架） test05.txt ：固定左端三向约束，右端节点竖向力
- 简支梁三点载荷 test_beam.txt ：两端约束 Y，跨中三节点竖向力

示例输入片段

- `test05.txt`（悬臂框架）片段：
  - `2 1 1 1 1 1`
  - `2 0.0 0.0`
  - `2 2.0 0.0`
  - `0 0 -1 -1`（左端三向约束）
  - `2 0 1 0`
  - `210000000000 0.3 0.0`
  - `0.01 1e-05 0.2`
  - `1 1 -1000 -1 1 0 0 0`（右端竖向-1000N）
- `test_beam.txt`（简支梁）片段：
  - `5 2 4 1 1 3`
  - `2 0.0 0.0` … `2 2.0 0.0`（5 节点沿 X 布置）
  - `0 0 -1 0`、`4 0 -1 0`（两端 Y 约束）
  - `2 0 1 0`、`2 1 2 0`、`2 2 3 0`、`2 3 4 0`（四个元素）
  - `210000000000 0.3 0.0`；`0.01 1e-05 0.2`（材料与截面，梁需 `Iz>0`）
  - `1 1 -1000 -1 1 0 0 0`、`1 1 -1000 -1 2 0 0 0`、`1 1 -1000 -1 3 0 0 0`（中跨三节点竖向力）
输出文件（Results.dat/Results_*.dat） 包含以下块：

- 全局与实体数据回显（节点/约束/元素/材料/截面/载荷）
- 节点位移（按 DOF 序）： Node Displacements （ barsystem.cpp:987-1001 ）
- 杆端力（局部端力 6 分量）： Element End Forces （ barsystem.cpp:1003-1016 ）
- 支座反力（约束节点处）： Support Reactions （ barsystem.cpp:1018-1033 ）
注意：

- 若输入为水平桁架对竖向节点力，体系在竖向接近奇异，位移可能极大或被数值裁剪，可视化需要使用 FRAME 或斜桁架。
可视化脚本 脚本： scripts\\visualize_results.py

- 解析 Results.dat 并绘制：
  - 原始几何（灰色）
  - 位移后的变形图（红色）
  - 支座反力箭头（蓝色）
  - 元素颜色随轴力大小渐变（plasma）
- 用法：
  - python scripts\\visualize_results.py --results build\\Results_frame.dat --scale 1000 --out build\\plot_frame.png
- 可选参数：
  - --hide-reactions 关闭支座反力箭头
  - --hide-forces 关闭元素轴力着色
  - --reac-scale <float> 调整箭头长度比例
- 解析稳健性：
  - 对于数值列相邻无空格的情况，使用正则提取数字，避免解析失败（见 parse_results： scripts\\visualize_results.py:11-108 ）
  - 变形倍率自动裁剪，避免极端值导致图形异常（见 plot_structure： scripts\\visualize_results.py:118-129 ）
  - 轴力着色采用 `plasma`，归一化范围按当前结果自动估计（见 plot_structure： scripts\\visualize_results.py:156-169 ）

命令行参数（fem_run.exe）

- `--input <path>` 指定输入文件，默认 `test05.txt`
- `--output <path>` 指定结果输出文件，默认 `Results.dat`
- `--quiet` 静默模式（减少控制台输出）
- `--stiff <path>` 输出单元刚度矩阵到文件（文本），默认路径 `ElemStiff.dat`
- `--no-stiff` 关闭单元刚度矩阵输出
- 示例：
  - `& .\\build\\bin\\Debug\\fem_run.exe --input test_beam.txt --output build\\Results_beam.dat --no-stiff`
  - `& .\\build\\bin\\Debug\\fem_run.exe --input test05.txt --output build\\Results_frame.dat --stiff build\\ElemStiff.txt`
核心函数参考

- 程序入口与参数解析： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:12-39
- 刚度组装 Skyline 下三角： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:523-614
- LDLT 求解器： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:631-680 （回代： 669-673 ）
- 位移输出： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:987-1001
- 杆端力输出： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:1003-1016
- 反力输出： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:1018-1033
- 载荷装配（节点力）： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:1035-1050
- 杆端内力计算（Truss/Frame）： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:1052-1101
已知限制与建议

- 载荷装配目前仅支持节点集中力；分布荷载、温度荷载等可通过扩展 LoadVectorAssembly 与 FixedEndForceCalcu 实现。
- 对近奇异体系（如水平桁架承受竖向力），建议改为刚架单元或提供适当约束。
- 支座反力采用 A·u - 已施加载荷 ，如需严格处理，应在装配时将边界条件与载荷一致性考虑完整。
- 数值稳定性：LDLT 对零主元采用微小对角保护（ 1e-12 ），在物理问题上应通过合理模型避免奇异。

故障诊断

- 位移几乎为零：检查是否错误使用 `TRUSS` 分析梁问题；梁需 `FRAME` 且 `Iz>0`
- 位移过大或 `NaN`：检查约束是否充分；脚本已进行裁剪，但建议修正模型
- 反力看起来不平衡：确认载荷方向与单位，支座反力计算为 `A·u - P`
- 可视化没有颜色变化：轴力着色依赖 `Element End Forces`；确认内力输出正常
常用命令速查

- 构建：
  - & 'C:\\Program Files\\Microsoft Visual Studio\\2022\\Community\\MSBuild\\Current\\Bin\\MSBuild.exe' .\\build\\FEM_2d.sln /p:Configuration=Debug /m
- 测试：
  - & .\\build\\bin\\Debug\\unit_tests.exe
- 运行：
  - & .\\build\\bin\\Debug\\fem_run.exe --no-stiff --input test05.txt --output build\\Results_frame.dat
  - & .\\build\\bin\\Debug\\fem_run.exe --no-stiff --input test_beam.txt --output build\\Results_beam.dat
- 可视化：
  - python scripts\\visualize_results.py --results build\\Results_frame.dat --scale 1000 --out build\\plot_frame.png
  - python scripts\\visualize_results.py --results build\\Results_beam.dat --scale 2000 --out build\\plot_beam.png
 - 一键构建（CMake）：
   - PowerShell：`& scripts\\build_and_test.ps1`
   - Bash：`bash scripts/build_and_test.sh`
后续可扩展方向

- 完善 LoadVectorAssembly 支持分布荷载与温度荷载（ barsystem.h 已声明多个 Load 类型）
- 在 InternalForceCalcu 中加入固端力叠加与截面力（剪力/弯矩）更丰富的输出
- 增加简支梁、连续梁、桁架网架等更多示例输入及对应图像
- 为 Python 可视化添加数值标签、剪力/弯矩图、变形倍率自动估计

附加说明

- 函数补充：
  - 杆端力数组初始化： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:978-985
  - 反力计算： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:1103-1113
  - 载荷装配（节点力）： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:1035-1050
- 结果字段含义：
  - Node Displacements：每行 `Node Ux Uy Rz`，单位与输入一致
  - Element End Forces：每行 `Elem Fx_i Fy_i Mz_i Fx_j Fy_j Mz_j`，为局部坐标端力
  - Support Reactions：每行 `Node Rx Ry Rz`，为约束自由度处的反力
- 建模与使用建议：
  - 梁/刚架问题请使用 `FRAME` 单元且截面 `Iz>0`，否则无弯曲刚度
  - 简支梁：两端约束 `Y`，不约束 `X/Rz`，中跨施加竖向集中力
  - 位移倍率 `--scale` 需与几何尺度匹配；脚本已对过大位移进行裁剪
  - 若出现奇异或位移异常，请检查约束是否充分以及单元类型是否正确
