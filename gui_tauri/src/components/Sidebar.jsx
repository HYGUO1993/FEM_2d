import { useState } from "react";
import { useStore } from "../store.js";
import * as ipc from "../ipc.js";
import { defaultModel, newEmptyModel } from "../model.js";

export default function Sidebar() {
  const projects = useStore((s) => s.projects);
  const currentProject = useStore((s) => s.currentProject);
  const [showNew, setShowNew] = useState(false);
  const [newName, setNewName] = useState("");

  // 导出当前模型到本地 JSON 文件
  async function handleExport() {
    const st = useStore.getState();
    const name = st.currentProject || "model";
    try {
      const savedPath = await ipc.exportModelFile(name, JSON.stringify(st.model, null, 2));
      st.setToast("已导出: " + savedPath);
    } catch (e) {
      const msg = (e && e.message) || (typeof e === "string" ? e : JSON.stringify(e));
      if (msg !== "用户取消了保存") st.setToast("导出失败: " + msg, true);
    }
  }

  // 从本地 JSON 文件导入模型
  async function handleImport() {
    const st = useStore.getState();
    try {
      const content = await ipc.importModelFile();
      const parsed = JSON.parse(content);
      if (!parsed || typeof parsed !== "object") {
        st.setToast("文件内容不是有效的模型 JSON", true);
        return;
      }
      // 兼容字段名差异: 后端 schema 用 model_json 传参, 但文件本身应直接是模型对象
      const model = parsed;
      if (!Array.isArray(model.nodes)) {
        st.setToast("模型缺少 nodes 数组", true);
        return;
      }
      st.setModel(model);
      st.setCurrentProject(model.title || "导入模型");
      st.resetView();
      st.setToast("已导入: " + (model.title || "未命名模型"));
    } catch (e) {
      const msg = (e && e.message) || (typeof e === "string" ? e : JSON.stringify(e));
      if (msg !== "用户取消了选择") st.setToast("导入失败: " + msg, true);
    }
  }

  async function handleNew() {
    const name = (newName || `项目 ${projects.length + 1}`).trim();
    if (!name) return;
    const st = useStore.getState();
    if (!projects.includes(name)) st.setProjects([...st.projects, name]);
    const m = newEmptyModel(name);
    st.setModel(m);
    st.setCurrentProject(name);
    st.resetView();
    try {
      await ipc.saveProject(name, m);
    } catch (e) {
      console.error("新建项目保存失败", e);
    }
    st.setToast("已创建项目: " + name);
    setShowNew(false);
    setNewName("");
  }

  async function handleSwitch(name) {
    if (name === currentProject) return;
    const st = useStore.getState();
    try {
      await ipc.saveProject(currentProject, st.model); // 切换前自动保存当前项目
      const json = await ipc.loadProject(name);
      st.setModel(JSON.parse(json));
      st.setCurrentProject(name);
      st.resetView();
    } catch (e) {
      st.setToast("切换项目失败: " + e.message, true);
    }
  }

  async function handleDelete(name) {
    const st = useStore.getState();
    try {
      await ipc.deleteProject(name);
      const remaining = st.projects.filter((p) => p !== name);
      st.setProjects(remaining);
      if (name === currentProject) {
        if (remaining.length) {
          const next = remaining[0];
          const json = await ipc.loadProject(next);
          st.setModel(JSON.parse(json));
          st.setCurrentProject(next);
          st.resetView();
        } else {
          st.setModel(defaultModel());
          st.setCurrentProject("");
        }
      }
      st.setToast("已删除项目: " + name);
    } catch (e) {
      st.setToast("删除失败: " + e.message, true);
    }
  }

  return (
    <aside className="sidebar">
      <div className="sidebar-header">
        <button
          className="btn primary block"
          onClick={() => setShowNew((v) => !v)}
        >
          + 新建项目
        </button>
        {showNew && (
          <div className="new-project-row">
            <input
              className="text-input"
              placeholder="项目名称"
              value={newName}
              onChange={(e) => setNewName(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter") handleNew();
                if (e.key === "Escape") setShowNew(false);
              }}
              autoFocus
            />
            <button className="btn small" onClick={handleNew}>
              创建
            </button>
          </div>
        )}
        <div className="file-actions">
          <button className="btn block" onClick={handleExport} title="将当前模型保存为本地 JSON 文件">
            保存到本地文件
          </button>
          <button className="btn block" onClick={handleImport} title="从本地 JSON 文件导入模型">
            导入本地文件
          </button>
        </div>
      </div>
      <div className="sidebar-list">
        {projects.length === 0 && <div className="sidebar-empty">暂无项目</div>}
        {projects.map((name) => (
          <div
            key={name}
            className={`project-item ${name === currentProject ? "active" : ""}`}
            onClick={() => handleSwitch(name)}
            title={name}
          >
            <span className="project-name">{name}</span>
            <button
              className="project-del"
              title="删除项目"
              onClick={(e) => {
                e.stopPropagation();
                handleDelete(name);
              }}
            >
              ✕
            </button>
          </div>
        ))}
      </div>
    </aside>
  );
}
