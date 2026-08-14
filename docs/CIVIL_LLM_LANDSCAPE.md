# Civil × LLM Agent 同行调研（2024-2026）

> 调研目的：对照同行进展，校准 FemLab Studio 的"确定性内核 + LLM 编排"战略定位。
> 调研时间：2026-08-14，信息来自 arXiv / 期刊 / 公开新闻。

---

## 一、同行地图

### 1. AIstructure-Copilot — 清华陆新征团队（生成式结构设计，商业化）
- **论文**：Qin S, Liao W, et al. *AIstructure-Copilot: assistant for generative AI-driven intelligent design of building structures*. Smart Constr. 2024. https://doi.org/10.55092/sc20240001
- **产品**：本地-云协同 CAD 插件（2025-12 与中望联合发布"中望版 AI-structure Copilot"，适配国产 CAD）
- **方法**：GAN / GNN / Diffusion 三种生成模型做剪力墙智能布置 + 材料用量预测 + 结构模型参数化构建
- **效果**：剪力墙/梁布置 2 小时 → 3~10 分钟；流程"初步布置 → 模型计算 → 方案优化"
- **本质**：**生成式设计**（生成模型直接出方案），不是 LLM 工具编排；聚焦住宅剪力墙结构设计，面向**生产工程师**
- 团队背景：陆新征（清华土木教授，防灾减灾，高被引学者）；覃思中、廖文杰等

### 2. "LLM-Empowered Agent for Reliable and Robust Structural Analysis"（arXiv 2507.02938, 2025-07）
- **核心实验**：Llama-3.3 70B 直接做梁结构分析 → **定性理解好、定量不可靠**（8 题基准，重复运行）
- **方案**：把结构分析**重构为代码生成任务**——CoT + few-shot 生成 OpenSeesPy 代码，自动执行取结果 → **准确率 >99%**
- **意义**：学术级证明了"LLM 不能直接算数，必须转确定性执行"——与我们的核心判断一致
- 作者：Liu, Geng 等（注：**不是**用户记忆中的 Chen Minhui，见附录）

### 3. Lightweight Multi-Agent System for 2D Frame（arXiv 2510.05414, 2025-10）
- **架构**：5 个专职 agent 流水线：问题分析 → 几何建模（**专家规则**推节点/单元）→ 代码翻译（OpenSeesPy）→ 模型校验（一致性检查）→ 荷载施加
- **结果**：20 基准 ×10 次重复，准确率 >80%，**超过 Gemini-2.5 Pro 与 ChatGPT-4o**；用轻量 Llama-3.3 70B
- **与我们最接近的同行**（同为 2D 框架分析）

### 4. Multi-Agent Architecture to Reduce Hallucinations（arXiv 2603.07728, 2026-03）
- **架构**：规划 agent 制定分步建模计划 → 节点/单元 agent **并行**组装几何 → 荷载 agent → 代码翻译 agent
- **结果**：20 框架 ×10 次重复，**18 例 100%、2 例 90%**——当前公开最佳成绩
- **要点**：任务分解 + 并行 agent + 一致性校验，显著减少长序列操作中的幻觉累积

### 5. Agentic LLMs for 3D Frame Systems（arXiv 2606.06525, 2026-06）
- **架构**：3D 框架投影为 2D 平面网格表示 + 楼层分解 agent + 节点/梁/板/柱 agent + SAP2000 脚本生成
- **结果**：10 个 3D 框架平均准确率 90%
- **信号**：同行已从 2D 走向 3D

### 6. 其他相关
- Liang-Zhou 等 *Automating Structural Engineering Workflows with LLM Agents*（Semantic Scholar）
- SSRN 综述 *LLM Applications in Civil and Structural Engineering: Review*
- 马智亮（清华土木教授）2025-08 观点：**"通用大模型在建筑领域离实用差很远"**——垂直化是共识
- CIC（香港建造业议会）Webinar：生成式 AI 工程设计应用

---

## 二、行业共识（已形成的三条）

1. **LLM 直接产出数值结果不可靠**——必须转"代码/工具 + 确定性执行"（2507.02938 用基准证明了这一点）
2. **任务分解 + 多 agent + 校验层**是提升长任务可靠性的标准解法（2603.07728 达 90-100%）
3. **轻量开源模型（Llama-3.3 70B 级）足够**做编排——不依赖顶级闭源模型即可达成 80-100%

## 三、对照分析：我们的位置

### 我们已有的差异化优势
| 维度 | 同行 | FemLab Studio |
|---|---|---|
| 确定性内核 | **生成 OpenSeesPy/SAP2000 脚本**（依赖外部软件，脚本层有额外脆弱性） | **原生 C++ FEM 内核**（femcli），零外部依赖 |
| 数值级互证 | 一致性/结构校验为主 | **独立参考求解器逐位对照**（13 案例 ≤1e-11）+ golden + 理论值 |
| 产品形态 | 论文原型 / 商业 CAD 插件 | 完整桌面应用：交互画布 + 精确内力图 + 发版闭环 |
| 定位 | 自动化分析 / 生产设计 | **教学**（交互实验、理论对照 5qL⁴/384EI、可视化反弯点）——同行无此定位 |

### 我们的差距（同行已验证而我们没有的）
1. **无量化基准**：同行都有公开 benchmark（8/20/20/10 题）+ 重复试验统计；我们只有案例级验证
2. **单 agent 循环**：我们仍是"单 agent 工具循环 8 轮"，同行已验证**多 agent 分解 + 并行**更稳（90-100% vs 我们"看运气"）
3. **无 3D**：同行已到 3D 框架（投影表示法可借鉴）
4. **无生成式设计**：AIstructure-Copilot 的 GAN/GNN/diffusion 布置能力我们完全没有（教学场景可暂缓）
5. **模型可复现性**：我们用闭源 API（deepseek-v4-flash）；同行用开源 Llama 可复现（教学场景其实无碍）

## 四、对路线图的启示（按优先级）

1. **建 Agent benchmark**：10-20 个框架/梁问题集 + 10 次重复试验，量化当前成功率
   ——没有数字就没有产品级可信度，也是论文/宣传的硬通货
2. **升级多 agent 架构**：规划 agent（分解任务）→ 几何 agent（专家规则，借鉴 2510.05414）→ 校验 agent（物理+结构校验）→ 求解/解读 agent；目标 90%+
3. **校验层强化**：在 schema 校验基础上加物理校验（机构判定、平衡检查）——同行已有"模型校验 agent"先例
4. **确定性内核服务化**（femcli 常驻/REST）：守住"原生内核"优势，同时暴露为可编程 API
5. **教学差异化深挖**：同行全部面向"自动化/生产"，**教学闭环（出题-求解-判分-讲解）是无人占据的蓝海**，是我们的主场

## 五、附录：关于"Chen Minhui"的核实

- 检索到 **Minghui Chen**（IEEE Xplore 作者页，https://ieeexplore.ieee.org/author/37716807800），未能从其公开页面确认土木教学 LLM 文章。
- 检索到的"LLM → 土木教学"方向未发现高引用论文；**最接近用户记忆的可能是 arXiv 2507.02938（Liu & Geng）**或教育向综述（EduLLMs 等通用教育方向）。
- 建议：若用户有更具体的线索（论文标题/会议/年份），可补充核实。**无论该文内容如何，我们判断"教学闭环 + 数值互证"的产品深度在同类中领先。**

## 参考链接

- https://ar5iv.labs.arxiv.org/html/2507.02938 （可靠鲁棒结构分析 agent）
- https://ar5iv.labs.arxiv.org/html/2510.05414 （2D 框架轻量多 agent）
- https://ar5iv.labs.arxiv.org/html/2603.07728 （多 agent 抑制幻觉）
- https://ar5iv.labs.arxiv.org/html/2606.06525 （3D 框架 agentic LLM）
- https://www.zwsoft.cn/news/294-15512.html （中望×清华 AI-structure Copilot）
- https://doi.org/10.55092/sc20240001 （AIstructure-Copilot 论文）
- https://m.mp.oeeee.com/oe/BAAFRD0000202508221115236.html （马智亮观点）
