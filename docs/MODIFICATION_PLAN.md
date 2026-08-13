# FEM_2d 详细修改计划（可执行版）

> 前置结论（已实测）：求解内核在标准工况正确；缺陷集中在端力坐标系、载荷未实现、奇异掩盖、Skyline 未省内存等。本计划按「先防错、再修错、后扩展」推进，每阶段给出**改动位置 → 代码 → 验收标准 → 回归测试 → 风险**。
> 行号以当前 `barsystem.cpp`（1131 行）为准。

---

## Phase 0｜建立回归基线与 CI（0.5 天）

**目的**：任何重构前先把"现在是正确的"固化为黄金测试，防止改坏。

| 步骤 | 内容 | 验收标准 |
|---|---|---|
| 0.1 | 在 `tests/` 增加 `golden/`：用 `verify/verify_solver.py` 生成悬臂/简支梁/45° 桁架的**解析基准值**（位移、端力、反力） | 生成 `golden_*.json` |
| 0.2 | `tests/test_basic.cpp` 增加黄金值回归：跑 3 个算例，与基准比对（容差 1e-9 相对） | `unit_tests` 全绿 |
| 0.3 | `.github/workflows/ci.yml`：配置 MSVC + MinGW + Linux 三平台矩阵构建、跑 `ctest` | CI 通过 |

**黄金基准（已实测，可直接抄）**
- 悬臂 test05：`Uy(node1)=-1.9047619048e-4`，`Mz(固定端)=1000.0`
- 简支梁 test_beam：`Ry(node0)=Ry(node4)=1500.0`
- 45° 桁架 test_inclined：整体端力 `(707.10678, 707.10678)`，投影轴力 `N=1000.0`

---

## Phase 1｜P0 修复（正确性 / 可用性）

### 1.1 奇异与病态诊断（0.5 天）
**改动位置**：`LDLTSolve`（`barsystem.cpp:653-702`）+ `main` 错误分支（`:321-336`）。

**问题**：零主元用 `1e-12` 兜底且恒返回 `true`（`:671,:677,:688`），`main` 的 `if(!LDLTSolve)` 永不触发。

**改法**（分解时判奇异）：
```cpp
// 分解主对角元处（原 :676）
double sumd = A[i][i];
for (int k = 0; k < i; ++k) sumd -= L[i][k]*D[k]*L[i][k];
double diag0 = fabs(A[i][i]);
if (diag0 <= 1e-30) return false;            // 结构奇异（约束不足/机构）
if (fabs(sumd) <= 1e-12 * diag0) return false; // 相对容差判病态
D[i] = sumd;
```
并在 `main` 输出：`"刚度矩阵奇异/病态：请检查约束是否充分（自由DOF=%d）"`。

**验收**：构造"机构"算例（桁架缺约束）→ 程序返回非零并给出诊断，而不是输出 NaN。
**风险**：极低；不影响良态问题。

---

### 1.2 单元杆端力改为局部坐标（0.5~1 天）
**改动位置**：`InternalForceCalcu`（`:1069-1118`）+ `EndInternalForceOutput`（`:1020-1033`）。

**问题**：`K_global × u_global` 得到**整体坐标**端力，却标注"局部坐标端力"（45° 桁架证实）。

**改法**：算出整体端力后，用变换矩阵转回局部。对 FRAME（6 维），局部端力 `f_local = T^T · f_global`，`T` 复用 `FrameElemStiffCalcu` 里 `:517-521` 的定义。对 TRUSS（4 维）：
```cpp
double c = pElem[i].dCos, s = pElem[i].dSin;
double N_i =  c*f[0] + s*f[1];      // 轴力
double V_i = -s*f[0] + c*f[1];      // 剪力
// ...同法得 j 端；FRAME 的 Mz 转动分量不变
```
**输出表头改为**：`Elem  N_i  V_i  M_i  N_j  V_j  M_j`（工程师直接可用）。
**同步**：`README` / 可视化脚本 / GUI 的"轴力色图"直接用 `N`，无需再投影。

**验收**：45° 桁架 `N_i=N_j=1000`；悬臂梁 `V_i=1000, M_i=1000`。
**风险**：需同步改 `scripts/visualize_results.py` 与 `fem_parser.py` 的列解析（向后兼容可加 `--legacy`）。

---

### 1.3 固端力框架 + 均布/温度荷载（2~3 天，核心功能解锁）
**改动位置**：`FixedEndForceCalcu`（`:714-777` 死代码）、`LoadVectorAssembly`（`:1052-1067`）、`InternalForceCalcu`。

**步骤**：
1. **修死代码**：让 `FixedEndForceCalcu` 真正写入 `pFixedEndF[6]`（当前结果全丢弃）。
2. **核对公式**（FRAME 局部坐标，构件长 L，荷载距 i 端 a，b=L−a）：
   - 均布荷载 q（全跨）：`N=0; V_i=V_j=qL/2; M_i=-qL²/12; M_j=+qL²/12`
   - 跨中集中力 P（距 i 端 a）：`V_i=P b²(3a+b)/L³; M_i=-P a b²/L²; V_j=P a²(a+3b)/L³; M_j=P a² b/L²`
   - 温度梯度（上下面 T0/T1，膨胀系数 α，梁高 h）：`N=EαA(T0+T1)/2; M=EαIz(T0-T1)/h`
3. **接入求解**：
   - 等效节点荷载 `P_eq = -FixedEndForce`，在 `LoadVectorAssembly` 中对作用单元叠加（需转整体坐标再装到对应 DOF）；
   - 求出位移后，真实杆端力 `= k·u + FixedEndForce`（`InternalForceCalcu` 里加回）。
4. **载荷类型枚举**：`barsystem.h` 已声明 10 种，本次先实现 `LATERAL_UNIFORM_PRESSURE`(3)、`LATERAL_FORCE`(2)、`TEMPERATURE`(9)、`MOMENT_ON_A_POINT`(4)，其余留接口。

**验收**：
- 悬臂梁满跨均布荷载 q：`V(固定端)=qL, M(固定端)=qL²/2, 挠度=qL⁴/(8EI)` 对照解析；
- 简支梁均布荷载：跨中挠度 `=5qL⁴/(384EI)`；
- 两端固定梁温度梯度：自应力弯矩对照。
**风险**：中。等效节点荷载的**符号与坐标变换**最易错，必须用 1.2 的同一套 T 并配解析回归。

---

## Phase 2｜P1 架构与健壮性

### 2.1 真正的剖面(Profile) LDLT，取代稠密重建（2~3 天）
**改动位置**：`LDLTSolve`（`:653-702`）。

**问题**：内部 `TwoArrayDoubAlloc(nRow,nRow)` 重建 n×n 稠密 A 与 L（`:655,:662`），等于放弃 Skyline 的全部省内存收益，数千 DOF 即 O(n²)。

**改法**：直接在剖面数组 `pGK` + 主元地址 `pDiag` 上做原地 LDLᵀ（标准 profile 算法）：
```
for i in 0..n:
  lo = i - (pDiag[i]-pDiag[i-1]) + 1          // 该行最早非零列
  for j in lo..i-1:
    s = K(i,j) - Σ_{k=max(lo,loj)}^{j-1} L(i,k) D(k) L(j,k)
    L(i,j) = s / D(j)
  D(i) = K(i,i) - Σ_{k=lo}^{i-1} L(i,k)² D(k)
```
前代/回代同样只遍历剖面内非零。**内存与运算量 → O(剖面长度)**。
**验收**：黄金算例数值与旧实现一致（容差内）；用随机稠密/稀疏 SPD 对比 `numpy.linalg`；基准测试 1 万 DOF 内存下降一个量级。
**风险**：中高。建议**保留旧稠密实现做对照**（`#ifdef` 或两个函数），逐算例 diff。

### 2.2 内存与资源管理（1~2 天）
- 全局 `new[]/delete[]` → `std::vector<double/int>` / `std::unique_ptr`（`main` 的 10 处裸分配、各 `TwoArray*`）；
- `TwoArrayDoubAlloc/IntAlloc` 失败 `exit(-1)`（`:954,:962`）改为返回/抛异常；
- 一次性消除 `main` 里失败分支的重复 `delete[]`（`:323-335` vs `:345-354`），用 RAII 自然回收。
**验收**：`-fsanitize=address`（Linux/MinGW）或 MSVC `/RTC1` 无泄漏、无越界。
**风险**：低-中；纯重构，数值不变。

### 2.3 输入校验（0.5 天）
**位置**：`main` 解析段（`:75-113`）。
- 单元 `iSection/iMaterial`、载荷 `iLoadedNode/iLoadedElem`、约束 `iNode` 越界检查并给出行号错误；
- 耦合标志由 `iaDOFIndex[j] >= 1e4`（浮点阈值，`:840,:849`）改为显式整型位/独立字段。
**验收**：喂入坏索引 → 明确报错而非崩溃/垃圾结果。

---

## Phase 3｜P2 工程化与文档

| # | 项 | 位置 | 说明 | 工作量 |
|---|---|---|---|---|
| 3.1 | 单元刚度缓存 | `GKAssembly`/`InternalForceCalcu` | 一次计算存 `Element::daKe[36]`，装配与内力共用，消除"两套可能分歧的刚度" | 0.5 天 |
| 3.2 | `ElemStiff.dat` 语义 | `Truss/FrameElemStiffCalcu`（`:456-459,:512-514`） | 当前写的是**未旋转**局部刚度；改为旋转后写，或明确标注 | 0.5 天 |
| 3.3 | 去除重复 `LengthSinCosCalcu` | `main:300` 与 `GKAssembly:555` | 只算一次 | 0.5 天 |
| 3.4 | 文档同步 | `README` | test05 参数、"局部坐标端力"说法、坐标系约定小节 | 0.5 天 |
| 3.5 | 仓库卫生 | `.gitignore` | 移除提交进仓库的 Linux 二进制 `build/bin/fem_run`/`unit_tests`；`build/` 整体忽略；产物走 Release/CI | 0.5 天 |
| 3.6 | GUI 入口 | `fem_gui.py:22` / `quick_start_gui.py` | 统一入口；`fem_gui.py` 顶部加 `sys.path.insert(0, os.path.dirname(__file__))` | 0.5 天 |
| 3.7 | 测试扩充 | `tests/` | 端力/反力/奇异/边界用例 + 黄金回归 | 1 天 |

---

## 验收总清单（Definition of Done）
- [x] `unit_tests` 全绿（含黄金回归 + 奇异用例 + 端力/反力）
- [x] 45° 桁架轴力直接输出 1000（局部坐标）
- [x] 悬臂均布荷载对照解析通过（`δ=qL⁴/8EI`、`V=M=qL²/2`，2026-08-12 实测）
- [x] 机构/缺约束算例明确报错（test_inclined 验证，2026-08-12）
- [x] 回归测试抓出 test_beam.txt 输入缺陷并修正（`0 -1 0` → `-1 -1 0`，2026-08-12）
- [ ] 简支均布荷载对照解析通过（可作 1.3 补充验收）
- [ ] 1 万 DOF 剖面求解内存/时间达标（剖面法，Phase 2.1）
- [ ] 静态分析无泄漏（Phase 2.2）
- [ ] CI 三平台通过；README 与代码一致（Phase 3）

**总工作量估计**：约 9~13 人天（不含 Phase 4 拓展）。

---

## 建议落地顺序
```
Phase 0（基线）→ 1.1 防崩 → 1.2 修坐标系 → 1.3 载荷（核心功能）
→ 2.1 剖面求解 → 2.2/2.3 健壮性 → Phase 3 收尾
```
先保证"不错、不崩"，再解锁"均布/温度荷载"这一最大功能缺口，最后上规模。

---

# Phase 4（新增 2026-08-12）: 架构重构 — JSON 统一模型 + femcli CLI

## 背景
用户反馈: GUI 打包难受(PyInstaller + venv 依赖)、项目粗糙、想走"CLI + GUI + LLM"现代工具形态(类比 OpenCode/WorkBuddy)。

## 已完成

| # | 项 | 状态 |
|---|----|------|
| 4.1 | `FemSolveModel()` 纯求解接口(不依赖 main/文件) | ✅ |
| 4.2 | 引入 `third_party/nlohmann/json.hpp`(v3.11.3 单头文件) | ✅ |
| 4.3 | `femcli.cpp`: JSON 输入/输出 + txt 兼容(`--legacy-txt`) + `validate` 子命令 | ✅ |
| 4.4 | CMake 增加 `femcore` 静态库 + `femcli` target;`barsystem.cpp` 用 `FEM_NO_MAIN` 宏控制 | ✅ |
| 4.5 | `docs/FEM_MODEL_SCHEMA.md`: 统一输入/输出 Schema(全字段映射) | ✅ |
| 4.6 | `scripts/txt2json.py`: 旧 txt → JSON 无损转换 | ✅ |
| 4.7 | `scripts/llm_run.py` + `docs/LLM_INTEGRATION.md`: LLM 自然语言 → 案例 JSON → 求解 | ✅ |
| 4.8 | `verify/verify_femcli.py`: 4 个黄金用例回归(位移/反力严格对比) | ✅ 全 PASS |

## 验证结果
- `verify/verify_femcli.py`: cantilever_frame / simply_supported_beam / triangle_truss / cantilever_udl **4/4 PASS,maxdiff=0**
- JSON 与 txt 输入求解结果逐位一致(maxdiff=0)
- 转换工具 txt2json 求解与手写 JSON 一致(maxdiff=0)

## Phase 5: Tauri GUI 重写(进行中)

### 已完成
- [x] 选型 Tauri v2(Rust 壳 + 前端静态资源,体积 ~12MB)
- [x] 前端纯原生 HTML/CSS/JS(无 Node 依赖),直接读写 JSON 模型
- [x] 后端 `solve_model` Tauri 命令:JSON 模型 → 临时文件 → femcli.exe → 结果 JSON
- [x] `withGlobalTauri: true` + `frontendDist: "../src"`(白屏修复,见下)
- [x] dev.bat / build_gui.bat 一键脚本(femcli 自动复制到产物目录)

### 白屏问题排查记录(2026-08-13)
**根因一:devUrl 依赖外部 dev server。** 旧配置 `devUrl: http://localhost:1420` + `beforeDevCommand` 起 python http.server。直接跑 debug exe 时若 server 未启动 → 白屏。
修复:配置改为仅 `frontendDist: "../src"`,静态资源由 Tauri 内置加载。

**根因二:withGlobalTauri 缺失。** 前端 `window.__TAURI__.core.invoke` 未注入,IPC 不可用 → 即使显示也无法求解。
修复:tauri.conf.json 加 `"withGlobalTauri": true`。

**根因三:构建锁文件 os error 5。** cargo/tauri dev 反复报 `.cargo-build-lock` 拒绝访问。
- 排查:非沙箱 cargo build 成功,tauri dev 失败 → tauri dev 内部 watcher 与 cargo 的锁冲突
- 修复:直接用 `cargo build` 产物跑 exe(绕过 tauri dev),或清理锁文件后重试
- 注:本项目在沙箱环境构建时,Windows 文件共享锁会导致 target 下已存在文件无法覆写;用全新 CARGO_TARGET_DIR 可绕过

**根因四(疑似):多实例窗口最小化。** 窗口 rect 显示 (-32000,-32000)(Windows 最小化位置),多个实例竞争时窗口可能被移到屏幕外 → 用户看到"白屏"。确保单实例运行。

### 验证结果
- [x] 窗口正常创建渲染(1296x839,非白屏,UI 内容分布正常)
- [x] 点击"求解"后结果区域内容增加 ~4 倍(IPC → femcli → 结果渲染全链路成功)
- [x] femcli 求解结果正确(位移/端力/反力)

### 待办
- [ ] 画布可视化(结构图/变形图)未做
- [ ] 正式打包(tauri build / NSIS / MSI),bundle.icon 配置
- [ ] 单实例保护(single-instance plugin)
- [ ] 接 LLM: 聊天式建模 + 案例库
