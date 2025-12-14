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
  - 组装载荷向量（目前支持节点集中力 FORCE_ON_NODE ）： barsystem.cpp:1026-1031
  - LDLT 求解自由度位移： barsystem.cpp:631-680
  - 杆端内力计算（按单元刚度乘以局部位移）： barsystem.cpp:1047-1069
  - 支座反力： barsystem.cpp:1047-1051 （约束 DOF 处的反力回写到 pLoadVect ）
  - 结果输出（位移/内力/反力）： barsystem.cpp:978-1024
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
示例文件：

- 悬臂梁（框架） test05.txt ：固定左端三向约束，右端节点竖向力
- 简支梁三点载荷 test_beam.txt ：两端约束 Y，跨中三节点竖向力
输出文件（Results.dat/Results_*.dat） 包含以下块：

- 全局与实体数据回显（节点/约束/元素/材料/截面/载荷）
- 节点位移（按 DOF 序）： Node Displacements （ barsystem.cpp:978-992 ）
- 杆端力（局部端力 6 分量）： Element End Forces （ barsystem.cpp:994-1007 ）
- 支座反力（约束节点处）： Support Reactions （ barsystem.cpp:1009-1024 ）
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
  - 对于数值列相邻无空格的情况，使用正则提取数字，避免解析失败（见 scripts\\visualize_results.py:93-122 ）
核心函数参考

- 程序入口与参数解析： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:12-39
- 刚度组装 Skyline 下三角： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:523-614
- LDLT 求解器： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:631-680 （回代： 669-673 ）
- 位移输出： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:978-992
- 杆端力输出： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:994-1007
- 反力输出： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:1009-1024
- 载荷装配（节点力）： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:1026-1031
- 杆端内力计算（Truss/Frame）： c:\\Users\\guoho\\Documents\\GitHub\\FEM_2d\\barsystem.cpp:1047-1069
已知限制与建议

- 载荷装配目前仅支持节点集中力；分布荷载、温度荷载等可通过扩展 LoadVectorAssembly 与 FixedEndForceCalcu 实现。
- 对近奇异体系（如水平桁架承受竖向力），建议改为刚架单元或提供适当约束。
- 支座反力采用 A·u - 已施加载荷 ，如需严格处理，应在装配时将边界条件与载荷一致性考虑完整。
- 数值稳定性：LDLT 对零主元采用微小对角保护（ 1e-12 ），在物理问题上应通过合理模型避免奇异。
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
后续可扩展方向

- 完善 LoadVectorAssembly 支持分布荷载与温度荷载（ barsystem.h 已声明多个 Load 类型）
- 在 InternalForceCalcu 中加入固端力叠加与截面力（剪力/弯矩）更丰富的输出
- 增加简支梁、连续梁、桁架网架等更多示例输入及对应图像
- 为 Python 可视化添加数值标签、剪力/弯矩图、变形倍率自动估计