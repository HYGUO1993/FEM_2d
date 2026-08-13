# FEM_2d 统一模型 Schema (JSON)

> 目标:让 **GUI / CLI / LLM 三个入口共享同一个模型格式**。
> 本 Schema 是唯一事实源。C++ 求解内核通过 `femcli` 读取本格式,输出结构化的 JSON 结果。

## 1. 设计原则

1. **人类可读 + LLM 可生成**:字段用全名、枚举用字符串,不用魔法数字。
2. **与 txt 一一对应**:旧的 `.txt` 输入可以无损转成 JSON(见 §4 对照表)。
3. **向后兼容**:`femcli` 同时支持 JSON 输入与 txt 输入;txt 模式行为与 `fem_run.exe` 完全一致。
4. **可扩展**:根对象带 `solver` 字段,为将来接入 OpenSees 等外部求解器预留。

## 2. 输入 Schema

```json
{
  "schemaVersion": "1.0",
  "title": "简支梁 - 均布节点荷载(示例)",
  "solver": "builtin",

  "nodes": [
    { "id": 0, "type": "frame", "x": 0.0, "y": 0.0 },
    { "id": 1, "type": "frame", "x": 0.5, "y": 0.0 },
    { "id": 2, "type": "frame", "x": 1.0, "y": 0.0 },
    { "id": 3, "type": "frame", "x": 1.5, "y": 0.0 },
    { "id": 4, "type": "frame", "x": 2.0, "y": 0.0 }
  ],

  "constraints": [
    { "node": 0, "dofs": ["ux", "uy"] },
    { "node": 4, "dofs": ["uy"] }
  ],

  "elements": [
    { "id": 0, "type": "frame", "nodeI": 0, "nodeJ": 1, "section": 0, "material": 0 },
    { "id": 1, "type": "frame", "nodeI": 1, "nodeJ": 2, "section": 0, "material": 0 },
    { "id": 2, "type": "frame", "nodeI": 2, "nodeJ": 3, "section": 0, "material": 0 },
    { "id": 3, "type": "frame", "nodeI": 3, "nodeJ": 4, "section": 0, "material": 0 }
  ],

  "materials": [
    { "id": 0, "E": 210000000000.0, "mu": 0.3, "alpha": 0.0 }
  ],

  "sections": [
    { "id": 0, "A": 0.01, "Iz": 1.0e-5, "height": 0.2 }
  ],

  "loads": [
    { "type": "nodalForce", "direction": "y", "value": -1000.0, "node": 1 },
    { "type": "nodalForce", "direction": "y", "value": -1000.0, "node": 2 },
    { "type": "nodalForce", "direction": "y", "value": -1000.0, "node": 3 }
  ]
}
```

### 2.1 字段说明

| 字段 | 类型 | 必填 | 说明 |
|------|------|------|------|
| `schemaVersion` | string | 是 | 固定 `"1.0"`,将来升级用 |
| `title` | string | 否 | 案例名称,便于 LLM 命名 |
| `solver` | string | 否 | `"builtin"`(默认)/ 预留 `"opensees"` |
| `nodes[].type` | string | 是 | `"truss"` / `"frame"`(对应 TRUSS_NODE/FRAME_NODE) |
| `constraints[].dofs` | string[] | 是 | `"ux"` / `"uy"` / `"rz"` 子集 |
| `elements[].type` | string | 是 | `"truss"` / `"frame"`(对应 TRUSS/FRAME) |
| `loads[].type` | string | 是 | 见 §3 荷载类型表 |

### 2.2 荷载类型枚举(字符串 ↔ txt 数值)

| JSON `type` | txt 数值 | 说明 |
|-------------|----------|------|
| `nodalForce` | 1 | 节点力 |
| `lateralForce` | 2 | 单元横向集中力 |
| `lateralUniformPressure` | 3 | 单元横向均布压力 |
| `momentOnPoint` | 4 | 节点集中弯矩 |
| `lateralLinearlyPressure` | 5 | 单元横向线性分布压力 |
| `axialPressure` | 6 | 单元轴向均布压力 |
| `axialForce` | 7 | 单元轴向集中力 |
| `momentOnBeam` | 8 | 梁上弯矩 |
| `temperature` | 9 | 温度荷载 |
| `supportMove` | 10 | 支座位移 |

`loads[]` 各字段与 txt 行 `iType iDirect dValue iLoadedElem iLoadedNode dPosition dT0 dT1` 一一对应:

| JSON 字段 | txt 列 | 说明 |
|-----------|--------|------|
| `type` | iType | 上表字符串枚举 |
| `direction` | iDirect | `"x"` / `"y"` / `"rz"`(对应 DIRECT_X/Y/R) |
| `value` | dValue | 荷载值(注意方向由 direction 决定,正负号需 LLM 理解) |
| `element` | iLoadedElem | 荷载作用的单元 id,无则 `-1` |
| `node` | iLoadedNode | 荷载作用的节点 id,无则 `-1` |
| `position` | dPosition | 作用位置或分布长度 |
| `T0` / `T1` | dT0 / dT1 | 温度上下表面值,无则 `0.0` |

## 3. 输出 Schema

```json
{
  "schemaVersion": "1.0",
  "solver": "builtin",
  "status": "ok",
  "message": "求解完成",
  "stats": {
    "nodeCount": 5, "elementCount": 4, "totalDOF": 15, "freeDOF": 10
  },
  "displacements": [
    { "node": 0, "ux": 0.0, "uy": 0.0, "rz": 0.0 },
    { "node": 1, "ux": 0.0, "uy": -0.0007143, "rz": -0.0006122 },
    { "node": 2, "ux": 0.0, "uy": -0.0008571, "rz": 0.0 },
    { "node": 3, "ux": 0.0, "uy": -0.0007143, "rz": 0.0006122 },
    { "node": 4, "ux": 0.0, "uy": 0.0, "rz": 0.0 }
  ],
  "endForces": [
    { "element": 0, "type": "frame", "N": -0.5, "V": 1500.0, "M": -200.0, "nodeJ": { "N": 0.5, "V": -1500.0, "M": 100.0 } },
    { "element": 1, "type": "frame", "N": 0.0, "V": 1000.0, "M": 100.0, "nodeJ": { "N": 0.0, "V": -1000.0, "M": 0.0 } }
  ],
  "reactions": [
    { "node": 0, "ux": 0.0, "uy": 1500.0, "rz": 0.0 },
    { "node": 4, "ux": 0.0, "uy": 1500.0, "rz": 0.0 }
  ]
}
```

> 说明:上表数值为示意,实际以 femcli 求解为准。
> `endForces` 的 `N/V/M` 为**局部坐标系**杆端力(与现有 Results.dat 表头 `Elem N_i V_i M_i N_j V_j M_j` 一致)。
> 出错时 `status` 为 `"error"`,`message` 含错误详情(如"刚度矩阵奇异或病态")。

## 4. txt ↔ JSON 对照表

| txt 块 | JSON 块 |
|--------|---------|
| 首行 `nNode nCons nElem nMat nSect nLoad` | 各数组的 length(隐含) |
| 节点行 `iType X Y` | `nodes[]` |
| 约束行 `iNode dx dy dr` | `constraints[]`(`0`=自由,`-1`=约束) |
| 单元行 `iType i j iSect iMat` | `elements[]` |
| 材料行 `E mu alpha` | `materials[]` |
| 截面行 `A Iz H` | `sections[]` |
| 荷载行 `type dir value elem node pos T0 T1` | `loads[]` |

## 5. 用法示例

```bash
# JSON 模式(推荐,LLM/GUI/CI 用)
femcli solve model.json -o result.json

# txt 模式(向后兼容)
femcli solve -i test_beam.txt -o result.json --legacy-txt

# 校验模型(不求解,快速检查 Schema)
femcli validate model.json
```
