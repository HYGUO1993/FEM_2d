# femcli 命令行使用指南

> femcli 是 FemLab Studio 的 JSON 模型命令行求解器（`femcli.cpp`），适合批处理、脚本集成和自动化验证。

## 快速开始

### 1. 构建

```bash
cmake -S . -B build
cmake --build build --parallel
```

构建产物：`build/bin/femcli`（或 `femcli.exe`）。

### 2. 校验模型

```bash
build/bin/femcli.exe validate examples/beam_simple.json
# 输出: 模型校验通过: 5 节点, 4 单元, 3 荷载
```

### 3. 求解

```bash
build/bin/femcli.exe solve examples/beam_simple.json -o result.json
```

`result.json` 包含：
- `displacements` — 节点位移（ux / uy / rz）
- `endForces` — 单元杆端内力（nodeI/nodeJ 的 N / V / M，局部坐标）
- `reactions` — 支座反力
- `stats` — 节点/单元/DOF 统计
- `status` / `message` — 求解状态

### 4. 批量求解脚本示例

```bash
# 循环求解 examples/ 下所有模型
for f in examples/*.json; do
  build/bin/femcli.exe solve "$f" -o "results/$(basename "$f" .json)_result.json"
done
```

## 命令

| 命令 | 说明 |
|---|---|
| `femcli validate <model.json>` | 校验模型合法性 |
| `femcli solve <model.json> -o <out.json>` | 求解并输出结果 |
| `femcli solve -i <legacy.txt> -o <out.json> --legacy-txt` | 兼容旧 txt 输入格式 |

## 模型格式

JSON 模型（见 [FEM_MODEL_SCHEMA.md](../FEM_MODEL_SCHEMA.md)）：

```json
{
  "schemaVersion": "1.0",
  "nodes": [{ "id": 0, "type": "frame", "x": 0.0, "y": 0.0 }],
  "constraints": [{ "node": 0, "dofs": ["ux", "uy"] }],
  "elements": [{ "id": 0, "type": "frame", "nodeI": 0, "nodeJ": 1, "section": 0, "material": 0 }],
  "materials": [{ "id": 0, "E": 210000000000.0, "mu": 0.3, "alpha": 0.0 }],
  "sections": [{ "id": 0, "A": 0.01, "Iz": 1e-5, "height": 0.1 }],
  "loads": [{ "type": "nodalForce", "direction": "y", "value": -10000.0, "node": 1 }]
}
```

字段速查：
- `nodes[].type`：`frame`（刚架，3 DOF）| `truss`（桁架，2 DOF）
- `constraints[].dofs`：`["ux","uy","rz"]` 表示该自由度被约束（位移=0）
- `loads[].direction`：`x` | `y` | `rz`，`value` 带符号（向下为负 y）

## 结果字段说明

| 字段 | 含义 |
|---|---|
| `displacements[].ux/uy/rz` | 节点位移（单位与模型一致：米 / 弧度） |
| `endForces[].nodeI/nodeJ` | 单元端部内力：`N` 轴力(N)、`V` 剪力(N)、`M` 弯矩(N·m)，**局部坐标** |
| `reactions[].rx/ry/rz` | 约束节点支座反力 |

## 退出码

- `0` — 成功（校验通过 / 求解完成）
- `非 0` — 失败（模型不合法 / 求解器报错），`stderr` 有错误信息

## 常见问题

**Q: 提示找不到模型文件？**
A: 路径用相对项目根目录或绝对路径均可。

**Q: 求解结果位移特别大？**
A: 检查约束是否充分（刚体位移未约束），或单元类型是否合适（梁用 `frame` 且 `Iz>0`）。

**Q: 与旧版 txt 输入兼容？**
A: 用 `--legacy-txt` 参数，输入旧格式 txt 文件。
