import { useEffect, useRef, useState } from "react";
import Sidebar from "./components/Sidebar.jsx";
import CanvasView from "./components/CanvasView.jsx";
import Toolbar from "./components/Toolbar.jsx";
import Inspector from "./components/Inspector.jsx";
import ResultsPanel from "./components/ResultsPanel.jsx";
import ChatPanel from "./components/ChatPanel.jsx";
import { useStore } from "./store.js";
import * as ipc from "./ipc.js";
import { defaultModel, validateModel, modelSummary } from "./model.js";
import { t } from "./i18n.js";

// 拖拽分隔条: 通过 onDelta(dx) 回调上报鼠标横向位移
function SplitHandle({ onDelta }) {
  const dragging = useRef(false);
  const lastX = useRef(0);

  function onMouseDown(e) {
    e.preventDefault();
    dragging.current = true;
    lastX.current = e.clientX;
    document.body.style.cursor = "col-resize";
    document.body.style.userSelect = "none";
    window.addEventListener("mousemove", onMouseMove);
    window.addEventListener("mouseup", onMouseUp);
  }
  function onMouseMove(e) {
    if (!dragging.current) return;
    const dx = e.clientX - lastX.current;
    lastX.current = e.clientX;
    onDelta(dx);
  }
  function onMouseUp() {
    dragging.current = false;
    document.body.style.cursor = "";
    document.body.style.userSelect = "";
    window.removeEventListener("mousemove", onMouseMove);
    window.removeEventListener("mouseup", onMouseUp);
  }
  useEffect(() => {
    return () => {
      window.removeEventListener("mousemove", onMouseMove);
      window.removeEventListener("mouseup", onMouseUp);
      document.body.style.cursor = "";
      document.body.style.userSelect = "";
    };
  }, []);
  return <div className="split-handle" onMouseDown={onMouseDown} />;
}

export default function App() {
  const currentProject = useStore((s) => s.currentProject);
  const model = useStore((s) => s.model);
  const solving = useStore((s) => s.solving);
  const solveTime = useStore((s) => s.solveTime);
  const toast = useStore((s) => s.toast);
  const solve = useStore((s) => s.solve);
  const sidebarWidth = useStore((s) => s.sidebarWidth);
  const rightWidth = useStore((s) => s.rightWidth);
  const llmPosition = useStore((s) => s.llmPosition);
  const setSidebarWidth = useStore((s) => s.setSidebarWidth);
  const setRightWidth = useStore((s) => s.setRightWidth);
  const setLlmPosition = useStore((s) => s.setLlmPosition);
  const canUndo = useStore((s) => s.past.length > 0);
  const canRedo = useStore((s) => s.future.length > 0);
  const theme = useStore((s) => s.theme);
  const lang = useStore((s) => s.lang);
  const displayOptions = useStore((s) => s.displayOptions);
  const setTheme = useStore((s) => s.setTheme);
  const setLang = useStore((s) => s.setLang);
  const setDisplayOption = useStore((s) => s.setDisplayOption);
  const [showSettings, setShowSettings] = useState(false);

  // 启动初始化: 拉项目列表 + LLM 配置 + 加载最近项目
  useEffect(() => {
    (async () => {
      try {
        const st = useStore.getState();
        const projects = await ipc.listProjects();
        st.setProjects(projects || []);
        try {
          const cfg = await ipc.getLlmConfig();
          if (cfg) {
            st.setLlmConfig({
              base_url: cfg.base_url || "",
              api_key: cfg.api_key || "",
              model: cfg.model || "",
            });
          }
        } catch {
          /* LLM 配置读取失败不阻塞启动 */
        }
        if (projects && projects.length) {
          const first = projects[0];
          const json = await ipc.loadProject(first);
          const parsed = JSON.parse(json);
          // 防御: 历史版本可能存成双重编码(字符串), 非对象则回退默认模型, 避免 .length 崩溃
          st.setModel(parsed && typeof parsed === "object" ? parsed : defaultModel());
          st.setCurrentProject(first);
        } else {
          st.setModel(defaultModel());
          st.setCurrentProject("未命名项目");
        }
        st.resetView();
      } catch (e) {
        console.error("初始化失败", e);
        const msg = (e && e.message) || (typeof e === "string" ? e : JSON.stringify(e));
        useStore.setState({ toast: { msg: "初始化失败: " + msg, isError: true } });
      }
    })();
  }, []);

  // 自动保存: model 变更后 debounce 500ms
  useEffect(() => {
    if (!currentProject) return;
    const t = setTimeout(() => {
      ipc.saveProject(currentProject, model).catch((e) =>
        console.error("自动保存失败", e)
      );
    }, 500);
    return () => clearTimeout(t);
  }, [model, currentProject]);

  // 快捷键: Ctrl+Enter 求解 / Ctrl+Z 撤销 / Ctrl+Y+Ctrl+Shift+Z 重做 / Delete 删除
  useEffect(() => {
    const h = (e) => {
      const st = useStore.getState();
      // 输入框内不触发全局快捷键 (除了 Ctrl+Enter)
      const tag = e.target?.tagName;
      const inInput = tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT";
      if ((e.ctrlKey || e.metaKey) && e.key === "Enter") {
        e.preventDefault();
        st.solve();
        return;
      }
      if (inInput) return;
      if ((e.ctrlKey || e.metaKey) && !e.shiftKey && e.key.toLowerCase() === "z") {
        e.preventDefault();
        st.undo();
        return;
      }
      if ((e.ctrlKey || e.metaKey) && (e.key.toLowerCase() === "y" || (e.shiftKey && e.key.toLowerCase() === "z"))) {
        e.preventDefault();
        st.redo();
        return;
      }
      if ((e.key === "Delete" || e.key === "Backspace") && st.selection) {
        e.preventDefault();
        const sel = st.selection;
        if (sel.type === "node") st.deleteNode(sel.id);
        else if (sel.type === "element") st.deleteElement(sel.id);
        else if (sel.type === "load") st.deleteLoad(sel.id);
        else if (sel.type === "constraint") st.removeConstraint(sel.id);
      }
      // 数字键 1-6 切换工具 (select/node/element/load/constraint/erase)
      const toolMap = {
        "1": "select", "2": "node", "3": "element",
        "4": "load", "5": "constraint", "6": "erase",
      };
      if (toolMap[e.key]) {
        st.setTool(toolMap[e.key]);
      }
    };
    window.addEventListener("keydown", h);
    return () => window.removeEventListener("keydown", h);
  }, []);

  const err = validateModel(model);
  const summary = modelSummary(model);

  return (
    <div className="app">
      <header className="app-header">
        <div className="brand">
          <span className="logo">FemLab</span>
          <div>
            <h1>FemLab Studio</h1>
            <p className="subtitle">二维杆系有限元 · 画布建模 · LLM 辅助</p>
          </div>
        </div>
        <div className="header-actions">
          <button
            className={`btn ghost ${showSettings ? "active" : ""}`}
            onClick={() => setShowSettings((v) => !v)}
            title={t("display")}
          >
            ⚙ {t("display")}
          </button>
          <button
            className="btn ghost"
            onClick={() => setLlmPosition(llmPosition === "right" ? "bottom" : "right")}
            title="切换 LLM 对话位置"
          >
            {llmPosition === "right" ? "LLM 移到底部" : "LLM 移到右侧"}
          </button>
          <button
            className="btn ghost"
            onClick={() => useStore.getState().undo()}
            disabled={!canUndo}
            title="撤销 (Ctrl+Z)"
          >
            ↩ 撤销
          </button>
          <button
            className="btn ghost"
            onClick={() => useStore.getState().redo()}
            disabled={!canRedo}
            title="重做 (Ctrl+Y)"
          >
            ↪ 重做
          </button>
          <button className="btn ghost" onClick={() => useStore.getState().resetView()}>
            重置视图
          </button>
          <button
            className="btn primary"
            onClick={() => solve()}
            disabled={solving || !currentProject}
          >
            {solving ? "求解中…" : "求解"} <span className="shortcut">Ctrl+Enter</span>
          </button>
        </div>
      </header>

      <main
        className="app-body"
        style={{ gridTemplateColumns: `${sidebarWidth}px 7px 1fr 7px ${rightWidth}px` }}
      >
        <Sidebar />
        <SplitHandle
          onDelta={(dx) => {
            const cur = useStore.getState().sidebarWidth;
            setSidebarWidth(Math.max(160, Math.min(480, cur + dx)));
          }}
        />
        <CanvasView />
        <SplitHandle
          onDelta={(dx) => {
            const cur = useStore.getState().rightWidth;
            setRightWidth(Math.max(280, Math.min(640, cur - dx)));
          }}
        />
        <aside className="right-pane">
          <Toolbar />
          <Inspector />
          <ResultsPanel />
          {llmPosition === "right" && <ChatPanel />}
        </aside>
      </main>

      {llmPosition === "bottom" && (
        <div className="llm-bottom-pane">
          <ChatPanel />
        </div>
      )}

      <footer className="status-bar">
        <span className="status-project">{currentProject || "未命名项目"}</span>
        <span>
          {summary.nodes} 节点 · {summary.elements} 单元 · {summary.dofs} DOF
        </span>
        {solveTime && <span>求解耗时 {solveTime}s</span>}
        <span className={err ? "badge err" : "badge ok"}>{err ? "✘ " + err : "✔ 模型有效"}</span>
      </footer>

      {toast && <div className={`toast ${toast.isError ? "error" : ""}`}>{toast.msg}</div>}

      {showSettings && (
        <div className="settings-panel">
          <div className="settings-group">
            <div className="settings-label">{t("theme")}</div>
            <div className="settings-options">
              {[
                ["dark", t("theme.dark")],
                ["light", t("theme.light")],
                ["ocean", t("theme.ocean")],
              ].map(([id, label]) => (
                <button
                  key={id}
                  className={`btn small ${theme === id ? "active" : ""}`}
                  onClick={() => setTheme(id)}
                >
                  {label}
                </button>
              ))}
            </div>
          </div>
          <div className="settings-group">
            <div className="settings-label">{t("language")}</div>
            <div className="settings-options">
              {[
                ["zh", "中文"],
                ["en", "English"],
              ].map(([id, label]) => (
                <button
                  key={id}
                  className={`btn small ${lang === id ? "active" : ""}`}
                  onClick={() => setLang(id)}
                >
                  {label}
                </button>
              ))}
            </div>
          </div>
          <div className="settings-group">
            <div className="settings-label">{t("display")}</div>
            {[
              ["nodeLabels", t("display.nodeLabels"), true],
              ["elementLabels", t("display.elementLabels"), false],
              ["loads", t("display.loads"), true],
              ["constraints", t("display.constraints"), true],
            ].map(([key, label]) => (
              <label key={key} className="check-row">
                <input
                  type="checkbox"
                  checked={displayOptions[key]}
                  onChange={(e) => setDisplayOption(key, e.target.checked)}
                />
                <span>{label}</span>
              </label>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
