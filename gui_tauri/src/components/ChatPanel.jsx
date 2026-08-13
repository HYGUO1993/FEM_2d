import { useState } from "react";
import { useStore } from "../store.js";
import * as ipc from "../ipc.js";
import { buildSystemPrompt } from "../llmSystemPrompt.js";
import { normalizeModel } from "../model.js";
import { t } from "../i18n.js";

/** 从 LLM 回复提取 JSON 模型: 取第一个 { 到最后一个 }, 校验是合法模型对象 */
function extractModelJson(text) {
  const start = text.indexOf("{");
  const end = text.lastIndexOf("}");
  if (start === -1 || end === -1 || end <= start) return null;
  try {
    const obj = JSON.parse(text.slice(start, end + 1));
    if (obj && Array.isArray(obj.nodes) && Array.isArray(obj.elements)) {
      return normalizeModel(obj);
    }
    return null;
  } catch {
    return null;
  }
}

/** 模型是否与当前模型完全一致(说明 LLM 没按需求生成, 只是回显) */
function isSameModel(a, b) {
  try {
    return JSON.stringify(a) === JSON.stringify(b);
  } catch {
    return false;
  }
}

export default function ChatPanel() {
  const chatMessages = useStore((s) => s.chatMessages);
  const pendingLlmModel = useStore((s) => s.pendingLlmModel);
  const [input, setInput] = useState("");
  const [showConfig, setShowConfig] = useState(false);
  const [sending, setSending] = useState(false);

  const llmConfig = useStore((s) => s.llmConfig);
  const setLlmConfig = useStore((s) => s.setLlmConfig);

  async function send() {
    const text = input.trim();
    if (!text || sending) return;
    const st = useStore.getState();
    st.pushChat("user", text);
    setInput("");

    const cfg = st.llmConfig;
    if (!cfg.base_url || !cfg.api_key || !cfg.model) {
      st.pushChat("assistant", "请先点击「设置」配置 base_url / api_key / model。");
      return;
    }

    setSending(true);
    const messages = [
      { role: "system", content: buildSystemPrompt(st.model) },
      ...st.chatMessages.slice(-12),
    ];
    try {
      const reply = await ipc.llmChat(
        { base_url: cfg.base_url, api_key: cfg.api_key, model: cfg.model },
        messages
      );
      st.pushChat("assistant", reply);
      const parsed = extractModelJson(reply);
      if (parsed) {
        if (isSameModel(parsed, st.model)) {
          st.pushChat(
            "assistant",
            "⚠️ 模型似乎与当前画布完全相同 —— LLM 可能没有按你的要求生成新模型，请重新描述需求，或先新建一个空白项目再提问。"
          );
        } else {
          st.setPendingLlmModel(parsed);
        }
      }
    } catch (e) {
      st.pushChat("assistant", "请求失败: " + (e.message || String(e)));
    } finally {
      setSending(false);
    }
  }

  function applyModel() {
    const st = useStore.getState();
    const m = st.pendingLlmModel;
    if (!m) return;
    st.setModel(m);
    st.setPendingLlmModel(null);
    st.setTool("select");
    st.resetView();
    st.pushChat(
      "assistant",
      `已应用模型（${m.nodes.length} 节点 / ${m.elements.length} 单元），正在求解…`
    );
    st.solve();
  }

  function saveConfig() {
    const st = useStore.getState();
    const cfg = {
      base_url: st.llmConfig.base_url.trim(),
      api_key: st.llmConfig.api_key.trim(),
      model: st.llmConfig.model.trim(),
    };
    ipc
      .setLlmConfig(cfg)
      .then(() => {
        st.setLlmConfig(cfg);
        setShowConfig(false);
        st.setToast("LLM 配置已保存");
      })
      .catch((e) => st.setToast("保存配置失败: " + e.message, true));
  }

  return (
    <div className="pane chat-pane">
      <div className="pane-head">
        <h3>{t("llmAssistant")}</h3>
        <button className="btn ghost small" onClick={() => setShowConfig((v) => !v)}>
          设置
        </button>
      </div>

      {showConfig && (
        <div className="llm-config">
          <label className="field">
            <span>base_url</span>
            <input
              type="text"
              placeholder="https://api.deepseek.com/v1"
              value={llmConfig.base_url}
              onChange={(e) => setLlmConfig({ ...llmConfig, base_url: e.target.value })}
            />
          </label>
          <label className="field">
            <span>api_key</span>
            <input
              type="password"
              placeholder="sk-..."
              value={llmConfig.api_key}
              onChange={(e) => setLlmConfig({ ...llmConfig, api_key: e.target.value })}
            />
          </label>
          <label className="field">
            <span>model</span>
            <input
              type="text"
              placeholder="deepseek-chat"
              value={llmConfig.model}
              onChange={(e) => setLlmConfig({ ...llmConfig, model: e.target.value })}
            />
          </label>
          <button className="btn primary small block" onClick={saveConfig}>
            保存配置
          </button>
        </div>
      )}

      <div className="chat-messages">
        {chatMessages.length === 0 && (
          <p className="placeholder">
            用自然语言描述结构模型, 例如: 「一根 2m 的简支梁, 跨中施加 10kN 竖向向下荷载」。
          </p>
        )}
        {chatMessages.map((m, i) => (
          <div key={i} className={`chat-msg ${m.role}`}>
            {m.content}
          </div>
        ))}

        {pendingLlmModel && (
          <div className="llm-model-card">
            <span>
              已生成模型: {pendingLlmModel.nodes.length} 节点 / {pendingLlmModel.elements.length}{" "}
              单元
            </span>
            <button className="btn primary small" onClick={applyModel}>
              应用到画布 + 求解
            </button>
          </div>
        )}
      </div>

      <div className="chat-input-row">
        <input
          className="text-input"
          placeholder="描述你的结构模型…"
          value={input}
          onChange={(e) => setInput(e.target.value)}
          onKeyDown={(e) => {
            if (e.key === "Enter") send();
          }}
        />
        <button className="btn primary" onClick={send} disabled={sending || !input.trim()}>
          {sending ? "…" : "发送"}
        </button>
      </div>
    </div>
  );
}
