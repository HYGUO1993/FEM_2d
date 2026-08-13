# FEM_2d + LLM 接入指南

> 目标: 用自然语言描述结构,LLM 生成符合 `FEM_MODEL_SCHEMA.md` 的 JSON 案例,
> 交给 `femcli` 求解。这是 CLI/GUI/LLM 三入口中的 LLM 入口。

## 1. 工作流

```
用户自然语言描述 ("一个 2 米简支梁,E=210GPa,截面 10cm²,跨中受 10kN 集中力")
        │
        ▼
┌──────────────────────────┐
│  LLM(注入 Schema 提示词) │  ← 关键: 把 Schema 文档片段注入 system prompt
└──────────────────────────┘
        │  产出 model.json
        ▼
  femcli solve model.json -o result.json
        │
        ▼
  返回结构化解算结果(位移/端力/反力)
```

## 2. LLM 提示词模板

### System Prompt(建议注入)

```
你是一个结构力学有限元建模助手。用户会用自然语言描述一个二维杆系结构,
你需要输出一个符合以下 JSON Schema 的模型文件(直接输出 JSON,不要多余文字):

- nodes: 节点数组,每项 {id, type: "truss"|"frame", x, y}
- constraints: 约束数组,每项 {node, dofs: ["ux","uy","rz"] 子集}
  (ux=水平位移, uy=竖向位移, rz=转角; 固定支座给全部三个)
- elements: 单元数组,每项 {id, type: "truss"|"frame", nodeI, nodeJ, section, material}
- materials: 材料数组,每项 {id, E(Pa), mu(泊松比), alpha(线膨胀系数)}
- sections: 截面数组,每项 {id, A(m²), Iz(m⁴), height(m)}
- loads: 荷载数组,每项:
  {type: "nodalForce", direction: "y", value: 负值向下, node}
  {type: "lateralUniformPressure", direction: "y", value: 均布压力, element, position: 分布长度}
  等等(见 FEM_MODEL_SCHEMA.md)

物理约定:
- 长度单位米,力单位牛,弹性模量 Pa。
- 重力/向下荷载 value 为负(方向由 direction 决定,y 向下为负)。
- 若用户只说"钢",E=2.1e11 Pa, mu=0.3; 混凝土 E=3e10, mu=0.2。
- 不要猜测不存在的边界条件; 约束不足会报奇异,提示用户补充。

用户: <自然语言描述>
```

### 示例对话

**用户**: "一个 2 米长的简支梁,跨中节点受 10kN 向下的集中力,E=210GPa,截面 A=0.01m², I=1e-5 m⁴"

**LLM 输出**(节选):
```json
{
  "schemaVersion": "1.0",
  "title": "2m 简支梁跨中集中力",
  "nodes": [
    { "id": 0, "type": "frame", "x": 0.0, "y": 0.0 },
    { "id": 1, "type": "frame", "x": 1.0, "y": 0.0 },
    { "id": 2, "type": "frame", "x": 2.0, "y": 0.0 }
  ],
  "constraints": [
    { "node": 0, "dofs": ["ux", "uy"] },
    { "node": 2, "dofs": ["uy"] }
  ],
  "elements": [
    { "id": 0, "type": "frame", "nodeI": 0, "nodeJ": 1, "section": 0, "material": 0 },
    { "id": 1, "type": "frame", "nodeI": 1, "nodeJ": 2, "section": 0, "material": 0 }
  ],
  "materials": [ { "id": 0, "E": 210000000000.0, "mu": 0.3, "alpha": 0.0 } ],
  "sections": [ { "id": 0, "A": 0.01, "Iz": 1e-5, "height": 0.1 } ],
  "loads": [ { "type": "nodalForce", "direction": "y", "value": -10000.0, "node": 1 } ]
}
```

## 3. 自动校验 + 求解脚本

见 `scripts/llm_run.py`:

```bash
# 把 LLM 输出存为 model.json 后
python scripts/llm_run.py model.json
# 输出: 校验结果 + femcli 求解结果摘要
```

## 4. 与 GUI 的关系

- 未来 Tauri/Electron GUI 直接读写同一 JSON 格式,模型可互相导入。
- LLM 生成的案例可直接在 GUI 中打开/求解/可视化。
- 每个案例即一个 .json 文件,天然适合 Git 管理、分享、教材复用。

## 5. 进阶: 批量案例生成

给 LLM 一个"案例库"system prompt(含多个已生成的例子做 few-shot),
可批量产出教学案例:简支梁、悬臂梁、三角桁架、多层框架、连续梁等。
