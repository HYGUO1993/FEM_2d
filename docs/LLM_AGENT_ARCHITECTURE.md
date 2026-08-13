# LLM Agent 架构设计（参考 Claude Code / OpenCode）

> 目标：把「自然语言 → JSON」的一次性调用，升级为**工具调用式 Agent**——
> LLM 通过函数调用驱动 FEM 应用（建模/求解/解读/修改），多轮迭代直到任务完成。
> 本文是设计稿，用于对齐架构方向，再分阶段实现。

---

## 1. 现状与差距

| 维度 | 现状（单轮 JSON 生成） | Agent 目标 |
|---|---|---|
| 交互模型 | 一次性：用户输入 → 返回模型 | 多轮循环：思考 → 调用工具 → 观察结果 → 再思考 |
| 能力边界 | 只会"生成模型 JSON" | 建模 + 求解 + 解读 + 修改 + 项目管理 |
| 上下文 | 当前模型 + 最近 12 条聊天 | 结构化会话状态 + 结果摘要 + 工具轨迹 |
| 错误处理 | 解析失败就放弃 | 自动重试 / 自修复 |
| 用户参与 | 手动点「应用+求解」 | 按需确认（危险操作才打断） |

---

## 2. 目标架构

```
┌─────────────────────────────────────────────────────────┐
│  前端 (React) — Agent 运行时 (orchestrator)               │
│                                                         │
│  ChatPanel ──► agentLoop()                               │
│                   │                                     │
│                   ▼                                     │
│  ┌──────────────────────────────┐                       │
│  │ 循环: 调 LLM → 解析 tool_call │                       │
│  │   → 执行工具 → 结果回填 → 再调  │  (最多 N 轮)          │
│  └──────────────────────────────┘                       │
│         │            │                                  │
│   ┌─────▼─────┐ ┌────▼──────────┐                       │
│   │ 工具注册表  │ │ 上下文管理器    │                       │
│   │ getModel   │ │ 会话消息       │                       │
│   │ applyModel │ │ 结果摘要压缩    │                       │
│   │ solve      │ │ 工具轨迹       │                       │
│   │ validate   │ └───────────────┘                       │
│   │ projects…  │                                        │
│   └─────┬─────┘                                        │
└─────────┼───────────────────────────────────────────────┘
          │ invoke
┌─────────▼───────────────────────────────────────────────┐
│  后端 (Rust)                                            │
│  llm_chat_tools: 转发 tools 参数到 OpenAI 兼容 API       │
│  solve_model / save_project / list_projects ... (已有)  │
└─────────────────────────────────────────────────────────┘
```

**关键决策：Agent 循环放在前端**（OpenCode/Claude Code 同为 client-side loop）。
理由：
- 前端已拥有全部"工具"（store actions + ipc 命令），零后端重写
- 中间步骤可实时可视化（工具卡片、流式打字）
- 用户确认机制天然在前端（apply_model 前弹确认）

---

## 3. 工具协议（OpenAI 兼容 function calling）

### 3.1 后端扩展
新增 `llm_chat_tools(config, messages, tools)` 命令：转发 `tools` + `tool_choice` 到
`/chat/completions`，返回完整 `message`（含 `tool_calls` 时原样返回）。现有 `llm_chat` 保留作
纯文本回退。

### 3.2 工具清单（v1 核心 6 个）

| 工具名 | 参数 | 作用 | 副作用 |
|---|---|---|---|
| `get_current_model` | - | 返回当前画布模型 JSON | 无 |
| `apply_model` | `model` (JSON 字符串) | 应用到画布（**需用户确认**）| 覆盖模型 |
| `solve` | - | 求解当前模型，返回结果 JSON | 无（只读计算）|
| `validate_model` | `model` | 校验模型合法性 | 无 |
| `get_result_summary` | - | 返回求解结果的压缩摘要（max\|N\|/max\|V\|/max\|M\|、跨中挠度等）| 无 |
| `list_projects` | - | 列出项目 | 无 |

v2 可加：`save_project`、`load_project`、`new_project`、`set_llm_config`。

### 3.3 工具定义示例（发给 LLM 的 JSON Schema）
```json
{
  "type": "function",
  "function": {
    "name": "apply_model",
    "description": "将给定的 FEM 模型 JSON 应用到画布（覆盖当前模型）。调用前先用 validate_model 校验。",
    "parameters": {
      "type": "object",
      "properties": {
        "model": { "type": "string", "description": "完整模型 JSON 字符串" }
      },
      "required": ["model"]
    }
  }
}
```

---

## 4. Agent 循环（agentLoop）

```
用户消息 M
  ─► messages = [system, ...history, M]
  ─► 调 llm_chat_tools(messages, tools)
  ─► 响应 R
     ├─ R 含 tool_calls?
     │    ├─ 对每个 tool_call:
     │    │    ├─ 危险工具(apply_model) → 弹确认卡片, 等用户点头
     │    │    ├─ 执行工具 → 结果 result
     │    │    └─ 把 {role:"tool", tool_call_id, content:result} 追加进 messages
     │    ├─ 循环上限(默认 8 轮) 或工具结果提示"任务完成" → 停止
     │    └─ 回到「调 llm_chat_tools」
     └─ R 无 tool_calls → 这是最终回复 → 展示
```

**每轮工具执行都渲染成步骤卡片**（类似 OpenCode 的 tool 块）：
```
⚙ 正在求解…        →  ✔ 求解完成 (V_max=10000N, M_max=30000N·m)
⚙ 校验模型…        →  ✔ 模型有效
📋 等待确认: 应用模型 (3节点/2单元)  [应用] [取消]
```

---

## 5. 上下文管理

### 5.1 消息结构（OpenAI 格式）
```
system   角色 + FEM 知识 + 工具清单 + 当前模型摘要
user     「建一个3x3框架, 每层10kN」
assistant tool_call: solve
tool     solve 结果(完整 JSON 仅保留给工具, 不再进 LLM)
assistant tool_call: get_result_summary
tool     摘要文本
assistant 最终自然语言回复
```

### 5.2 结果摘要压缩（关键！防 token 爆炸）
solve 返回的完整 `endForces`/`displacements` 太大，喂回 LLM 前必须压缩：
```json
{
  "status": "ok",
  "maxAbsN": 12345, "maxAbsN_elem": 3,
  "maxAbsV": 10000, "maxAbsV_elem": 0,
  "maxAbsM": 30000, "maxAbsM_elem": 0,
  "maxDisp_mm": 42.9, "maxDisp_node": 3,
  "reactions": [{"node":0,"rx":0,"ry":10000,"rz":30000}]
}
```
只有 LLM **主动要求**完整数据（如"列出所有节点位移"）时才传全量。

### 5.3 会话状态
- 前端维护 `sessionMessages`（含 tool 轨迹，最多保留最近 ~40 条）
- 会话可「新建/清空」；历史会话暂存 store（`chatSessions`），v2 持久化到磁盘

---

## 6. 用户确认机制（对齐 Claude Code 的权限模型）

| 工具 | 策略 |
|---|---|
| `get_current_model` / `solve` / `validate_model` / `get_result_summary` / `list_projects` | 自动执行，无需确认 |
| `apply_model`（覆盖画布） | **必须确认**：显示模型预览卡片（节点/单元/荷载数），用户点「应用」才执行 |

好处：LLM 可以自由"试算/试探"，只有真正改动用户数据时才打断。

---

## 7. 系统提示词（Agent 版）

```
你是 FemLab Studio 的结构力学建模 Agent。
你可以通过工具完成: 查看当前模型、生成/修改模型、求解、解读结果。

工作流程:
1. 用户提出需求 → 先用 get_current_model 了解现状
2. 生成或修改模型 → 先用 validate_model 校验, 再 apply_model
3. 求解 → solve → get_result_summary 解读(与理论值对比, 如 PL、qL²/8、5qL⁴/384EI)
4. 用自然语言向用户汇报结论

规则:
- 新结构必须按用户尺寸/跨度/层数完整生成, 禁止回显当前模型
- 单位换算: mm→m, kN→N
- 求解失败(约束不足等) → 自动修复模型后重试, 最多2次
- 最终回复用中文/英文(跟随用户语言), 简明扼要, 可含关键数值
```

---

## 8. 分阶段实施

### Phase 1 — 工具循环骨架（本次）
- [ ] 后端 `llm_chat_tools` 命令（透传 tools 参数）
- [ ] 前端 `agentLoop()` + 工具注册表（6 个核心工具）
- [ ] 步骤卡片渲染 + apply_model 确认卡片
- [ ] 系统提示词 Agent 版 + 结果摘要压缩
- [ ] 循环上限/超时/停止按钮

### Phase 2 — 解读与教学
- [ ] 结果解读模板（对比理论值、判断合理性）
- [ ] 结构类型诊断（静定/超静定次数）
- [ ] 方案对比（简支 vs 连续并排）

### Phase 3 — 自修复与持久化
- [ ] 求解失败自动修复循环
- [ ] 会话持久化（磁盘）+ 多会话管理
- [ ] 流式输出（SSE）提升体验

---

## 9. 风险与对策

| 风险 | 对策 |
|---|---|
| 工具调用格式不标准（部分模型不支持 function calling）| 后端探测：响应无 tool_calls 但有 JSON 动作块时，前端解析兜底 |
| 循环失控（LLM 反复调用）| 硬上限 8 轮 + 每轮后判断"任务完成"关键词 + 停止按钮 |
| token 膨胀 | 结果摘要压缩 + 历史裁剪（只保留 tool 轨迹摘要） |
| apply_model 误覆盖 | 确认卡片 + 撤销可用（store 历史栈已支持） |
