import { useState } from "react";
import { useStore } from "../store.js";
import * as ipc from "../ipc.js";
import { buildAgentSystemPrompt } from "../llmSystemPrompt.js";
import { normalizeModel, newEmptyModel } from "../model.js";
import { agentLoop } from "../agent.js";
import { TUTORIAL_PROJECTS } from "../tutorialProjects.js";
import { t } from "../i18n.js";

/** 从 LLM 回复提取 JSON 模型 (兼容旧的非 agent 模式, 暂保留) */
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

export default function ChatPanel() {  const chatMessages = useStore((s) => s.chatMessages);
  const pendingLlmModel = useStore((s) => s.pendingLlmModel);
  const [input, setInput] = useState("");
  const [showConfig, setShowConfig] = useState(false);
  const [sending, setSending] = useState(false);
  const [showTutorial, setShowTutorial] = useState(false);

  const llmConfig = useStore((s) => s.llmConfig);
  const setLlmConfig = useStore((s) => s.setLlmConfig);

  async function send() {
    const text = input.trim();
    if (!text || sending) return;
    const st = useStore.getState();
    st.pushChat("user", text);
    setInput("");

    // —— /help 命令: 显示帮助 + 教程工程 ——
    if (text.startsWith("/help") || text.startsWith("/教程") || text === "/") {
      setShowTutorial(true);
      st.pushChat("assistant", t("help.intro"));
      return;
    }
    if (text.startsWith("/new")) {
      const m = newEmptyModel(t("unnamed"));
      st.setModel(m);
      st.setCurrentProject(t("unnamed"));
      st.resetView();
      st.pushChat("assistant", t("help.newDone"));
      return;
    }

    const cfg = st.llmConfig;
    if (!cfg.base_url || !cfg.api_key || !cfg.model) {
      st.pushChat("assistant", t("chat.needsConfig"));
      return;
    }

    setSending(true);
    // Agent 模式: 工具调用式循环
    const system = { role: "system", content: buildAgentSystemPrompt(st.model) };
    const history = st.chatMessages.slice(-10).map((m) => ({
      role: m.role === "user" ? "user" : "assistant",
      content: typeof m.content === "string" ? m.content : "",
    }));
    const messages = [system, ...history, { role: "user", content: text }];
    const reply = await agentLoop({
      config: { base_url: cfg.base_url, api_key: cfg.api_key, model: cfg.model },
      messages,
      onStep: (kind, payload) => {
        // 工具步骤实时渲染到聊天区 (轻量: 用 toast 提示, 详细步骤卡片后续版)
        if (kind === "tool") {
          const names = { get_current_model: "查看模型", validate_model: "校验模型", solve: "求解", get_result_summary: "读取结果摘要", apply_model: "应用模型" };
          st.pushChat("assistant", `⚙ ${names[payload.name] || payload.name}…`);
        }
        if (kind === "error") {
          st.pushChat("assistant", "⚠ " + payload);
        }
      },
    });
    st.pushChat("assistant", reply);
    setSending(false);
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
      t("chat.applied", { n: m.nodes.length, e: m.elements.length })
    );
    st.solve();
  }

  /** 加载内置教程工程到画布 */
  function loadTutorial(p) {
    const st = useStore.getState();
    try {
      const parsed = JSON.parse(p.json);
      st.setModel(parsed);
      st.setCurrentProject(p.desc);
      st.resetView();
      st.setTool("select");
      st.setToast(t("help.loaded", { name: p.desc }));
    } catch (e) {
      st.setToast(t("help.loadFailed", { msg: e.message }), true);
    }
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
        st.setToast(t("chat.configSaved"));
      })
      .catch((e) => st.setToast(t("chat.configSaveFailed", { msg: e.message }), true));
  }

  return (
    <div className="pane chat-pane">
      <div className="pane-head">
        <h3>{t("llmAssistant")}</h3>
        <button className="btn ghost small" onClick={() => setShowConfig((v) => !v)}>
          {t("chat.settings")}
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
            {t("chat.saveConfig")}
          </button>
        </div>
      )}

      <div className="chat-messages">
        {chatMessages.length === 0 && (
          <p className="placeholder">
            {t("chat.example")}
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
              {t("chat.generated", { n: pendingLlmModel.nodes.length, e: pendingLlmModel.elements.length })}
            </span>
            <button className="btn primary small" onClick={applyModel}>
              {t("chat.apply")}
            </button>
          </div>
        )}

        {showTutorial && (
          <div className="tutorial-panel">
            <div className="tutorial-head">
              <span>{t("help.tutorialTitle")}</span>
              <button className="btn small" onClick={() => setShowTutorial(false)}>
                ✕
              </button>
            </div>
            <div className="tutorial-list">
              {TUTORIAL_PROJECTS.map((p) => (
                <button
                  key={p.id}
                  className="tutorial-item"
                  onClick={() => loadTutorial(p)}
                  title={t("help.loadTip")}
                >
                  <span className="tutorial-name">{p.desc}</span>
                  <span className="tutorial-go">→</span>
                </button>
              ))}
            </div>
          </div>
        )}
      </div>

      <div className="chat-input-row">
        <input
          className="text-input"
          placeholder={t("chat.placeholder")}
          value={input}
          onChange={(e) => setInput(e.target.value)}
          onKeyDown={(e) => {
            if (e.key === "Enter") send();
          }}
        />
        <button className="btn primary" onClick={send} disabled={sending || !input.trim()}>
          {sending ? "…" : t("chat.send")}
        </button>
      </div>
    </div>
  );
}
