// FemLab Agent: 工具调用式 LLM 循环 (参考 Claude Code / OpenCode client-side loop)
// agentLoop(messages) → 调 llm_chat_tools → 执行 tool_calls → 结果回填 → 再调, 直到无工具调用

import * as ipc from "./ipc.js";
import { validateModel } from "./model.js";
import { useStore } from "./store.js";
import { buildAgentSystemPrompt } from "./llmSystemPrompt.js";

const MAX_ROUNDS = 8; // 循环硬上限, 防止失控

// —— 工具定义 (OpenAI function calling schema) ——
export const AGENT_TOOLS = [
  {
    type: "function",
    function: {
      name: "get_current_model",
      description: "返回当前画布上的 FEM 模型 JSON (字符串)。用于了解当前结构。",
      parameters: { type: "object", properties: {}, required: [] },
    },
  },
  {
    type: "function",
    function: {
      name: "validate_model",
      description: "校验一个模型 JSON 是否合法 (节点数/单元数/约束/引用)。返回错误消息或 OK。",
      parameters: {
        type: "object",
        properties: {
          model: { type: "string", description: "完整模型 JSON 字符串" },
        },
        required: ["model"],
      },
    },
  },
  {
    type: "function",
    function: {
      name: "solve",
      description: "对当前画布模型执行有限元求解。返回完整结果 JSON。",
      parameters: { type: "object", properties: {}, required: [] },
    },
  },
  {
    type: "function",
    function: {
      name: "get_result_summary",
      description: "返回当前求解结果的压缩摘要 (max|N|/max|V|/max|M|/最大位移/反力)，用于解读。",
      parameters: { type: "object", properties: {}, required: [] },
    },
  },
  {
    type: "function",
    function: {
      name: "apply_model",
      description: "将给定模型 JSON 应用到画布 (覆盖当前模型)。应用前请先 validate_model。",
      parameters: {
        type: "object",
        properties: {
          model: { type: "string", description: "完整模型 JSON 字符串" },
        },
        required: ["model"],
      },
    },
  },
];

// —— 结果摘要压缩 (防 token 爆炸) ——
export function summarizeResults(results) {
  if (!results) return "暂无求解结果";
  if (results.status !== "ok") return `求解失败: ${results.message || "未知错误"}`;
  const ef = results.endForces || [];
  let maxN = 0, maxV = 0, maxM = 0, nE = "-", vE = "-", mE = "-";
  for (const e of ef) {
    for (const side of ["nodeI", "nodeJ"]) {
      const f = e[side];
      if (!f) continue;
      if (Math.abs(f.N) > maxN) { maxN = Math.abs(f.N); nE = e.element; }
      if (Math.abs(f.V) > maxV) { maxV = Math.abs(f.V); vE = e.element; }
      if (Math.abs(f.M) > maxM) { maxM = Math.abs(f.M); mE = e.element; }
    }
  }
  let maxDisp = 0, dispNode = "-";
  for (const d of results.displacements || []) {
    const mag = Math.hypot(d.ux || 0, d.uy || 0);
    if (mag > maxDisp) { maxDisp = mag; dispNode = d.node; }
  }
  return {
    status: "ok",
    maxAbsN: maxN, maxAbsN_elem: nE,
    maxAbsV: maxV, maxAbsV_elem: vE,
    maxAbsM: maxM, maxAbsM_elem: mE,
    maxDisp_m: maxDisp, maxDisp_node: dispNode,
    reactions: (results.reactions || []).slice(0, 6),
    nodeCount: results.stats?.nodeCount,
    elementCount: results.stats?.elementCount,
  };
}

// —— JSON 容错解析: 修复 LLM 输出的常见小错 (markdown 包裹/尾逗号/缺尾括号/截断) ——
export function parseModelJson(text) {
  if (typeof text !== "string") return null;
  let s = text.trim();
  // 1. 剥 markdown 代码块 ```json ... ```
  const fence = s.match(/```(?:json)?\s*([\s\S]*?)```/);
  if (fence) s = fence[1].trim();
  // 2. 尝试直接解析
  try {
    return JSON.parse(s);
  } catch {
    /* 继续容错 */
  }
  // 3. 去掉首尾非 JSON 内容 (取第一个 { 到最后一个 })
  const start = s.indexOf("{");
  const end = s.lastIndexOf("}");
  if (start >= 0 && end > start) s = s.slice(start, end + 1);
  // 4. 修尾逗号 (},] 前多余的逗号)
  s = s.replace(/,\s*([}\]])/g, "$1");
  // 5. 补缺的右括号 (数左/右括号差)
  for (let i = 0; i < 8; i++) {
    try {
      return JSON.parse(s);
    } catch {
      /* 尝试补一个括号 */
    }
    const open = (s.match(/{/g) || []).length;
    const close = (s.match(/}/g) || []).length;
    if (open > close) s += "}";
    else if (open < close) break;
    else break;
  }
  try {
    return JSON.parse(s);
  } catch {
    return null;
  }
}

// —— 工具执行器 (前端 store 直接操作) ——
async function executeTool(name, args, st, onStep) {
  switch (name) {
    case "get_current_model":
      return JSON.stringify(st.model);
    case "validate_model": {
      const m = parseModelJson(args.model);
      if (!m) return `模型 JSON 解析失败: 不是合法的 JSON (请输出完整无省略的 JSON, 不要用 markdown 代码块或省略号)`;
      const err = validateModel(m);
      if (err) return `校验失败: ${err} (生成后请核对节点/单元数量是否与公式一致)`;
      return `OK: 模型有效 (${m.nodes.length} 节点 / ${m.elements.length} 单元 / ${(m.constraints || []).length} 约束 / ${(m.loads || []).length} 荷载)`;
    }
    case "solve": {
      const res = await st.solve();
      return JSON.stringify(res);
    }
    case "get_result_summary":
      return JSON.stringify(summarizeResults(st.results));
    case "apply_model": {
      const m = parseModelJson(args.model);
      if (!m) return `模型 JSON 解析失败: 不是合法的 JSON`;
      st.setModel(m);
      st.setTool("select");
      st.resetView();
      return `OK: 模型已应用 (${m.nodes?.length ?? 0} 节点 / ${m.elements?.length ?? 0} 单元)`;
    }
    default:
      return `错误: 未知工具 ${name}`;
  }
}

/**
 * Agent 循环: 从初始 messages 开始, 执行工具直到 LLM 给出最终回复或达上限
 * onStep: (kind, payload) 回调用于 UI 渲染 (step/tool/done/error)
 * 返回最终回复文本
 */
export async function agentLoop({ config, messages, onStep }) {
  const st = () => useStore.getState();
  let msgs = [...messages];
  for (let round = 0; round < MAX_ROUNDS; round++) {
    onStep?.("thinking", null);
    let reply;
    try {
      reply = await ipc.llmChatTools(config, msgs, AGENT_TOOLS);
    } catch (e) {
      const msg = e?.message || String(e);
      onStep?.("error", msg);
      return `请求失败: ${msg}`;
    }
    const content = reply.content || "";
    const toolCalls = reply.tool_calls;

    // 无工具调用 → 最终回复
    if (!toolCalls || toolCalls.length === 0) {
      onStep?.("done", content);
      return content;
    }

    // 有工具调用 → 逐条执行
    // 清洗 assistant 消息: 只保留 OpenAI 兼容字段 (剥离 reasoning_content 等推理模型特有字段,
    // 部分 API 不接受回传; content 补空串防 null)
    msgs.push({
      role: reply.role || "assistant",
      content: typeof reply.content === "string" ? reply.content : "",
      tool_calls: reply.tool_calls,
    });
    for (const tc of toolCalls) {
      const fn = tc.function || {};
      const name = fn.name || "";
      let args = {};
      try { args = JSON.parse(fn.arguments || "{}"); } catch { args = {}; }
      onStep?.("tool", { name, args });
      let result;
      try {
        result = await executeTool(name, args, st(), onStep);
      } catch (e) {
        result = `工具执行错误: ${e?.message || String(e)}`;
      }
      msgs.push({ role: "tool", tool_call_id: tc.id || "", content: result });
    }
  }
  const msg = `已达到最大循环次数 (${MAX_ROUNDS} 轮)，请简化请求。`;
  onStep?.("error", msg);
  return msg;
}
