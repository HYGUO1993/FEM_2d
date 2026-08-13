import { useState } from "react";
import { useStore } from "../store.js";
import * as ipc from "../ipc.js";
import { defaultModel, newEmptyModel } from "../model.js";
import { t } from "../i18n.js";

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
      st.setToast(t("toast.exported", { path: savedPath }));
    } catch (e) {
      const msg = (e && e.message) || (typeof e === "string" ? e : JSON.stringify(e));
      if (msg !== t("toast.canceledSave")) st.setToast(t("toast.exportFailed", { msg }), true);
    }
  }

  // 从本地 JSON 文件导入模型
  async function handleImport() {
    const st = useStore.getState();
    try {
      const content = await ipc.importModelFile();
      const parsed = JSON.parse(content);
      if (!parsed || typeof parsed !== "object") {
        st.setToast(t("toast.invalidJson"), true);
        return;
      }
      // 兼容字段名差异: 后端 schema 用 model_json 传参, 但文件本身应直接是模型对象
      const model = parsed;
      if (!Array.isArray(model.nodes)) {
        st.setToast(t("toast.missingNodes"), true);
        return;
      }
      st.setModel(model);
      st.setCurrentProject(model.title || t("importedModel"));
      st.resetView();
      st.setToast(t("toast.imported", { name: model.title || t("unnamed") }));
    } catch (e) {
      const msg = (e && e.message) || (typeof e === "string" ? e : JSON.stringify(e));
      if (msg !== t("toast.canceledPick")) st.setToast(t("toast.importFailed", { msg }), true);
    }
  }

  async function handleNew() {
    const name = (newName || t("projectDefault", { n: projects.length + 1 })).trim();
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
    st.setToast(t("toast.created", { name }));
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
      st.setToast(t("toast.switchFailed", { msg: e.message }), true);
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
      st.setToast(t("toast.deleted", { name }));
    } catch (e) {
      st.setToast(t("toast.deleteFailed", { msg: e.message }), true);
    }
  }

  // 复制项目: 读取原模型 → 以新名字保存
  async function handleDuplicate(name) {
    const st = useStore.getState();
    try {
      const json = await ipc.loadProject(name);
      let newName = name + t("dupeSuffix");
      let i = 1;
      while (st.projects.includes(newName)) {
        newName = `${name}${t("dupeSuffixN", { n: i })}`;
        i++;
      }
      await ipc.saveProject(newName, json);
      st.setProjects([...st.projects, newName]);
      st.setToast(t("toast.duplicated", { name: newName }));
    } catch (e) {
      st.setToast(t("toast.duplicateFailed", { msg: e.message }), true);
    }
  }

  // 改名: 保存新名 + 删除旧名 (后端无 rename, 用 save+delete 实现)
  const [renaming, setRenaming] = useState(null);
  const [renameVal, setRenameVal] = useState("");
  async function handleRename(oldName) {
    const newName = renameVal.trim();
    if (!newName || newName === oldName) {
      setRenaming(null);
      return;
    }
    if (projects.includes(newName)) {
      st.setToast(t("toast.renameExists"), true);
      return;
    }
    const st = useStore.getState();
    try {
      const json = await ipc.loadProject(oldName);
      await ipc.saveProject(newName, json);
      await ipc.deleteProject(oldName);
      const list = st.projects.map((p) => (p === oldName ? newName : p));
      st.setProjects(list);
      if (st.currentProject === oldName) st.setCurrentProject(newName);
      st.setToast(t("toast.renamed", { name: newName }));
    } catch (e) {
      st.setToast(t("toast.renameFailed", { msg: e.message }), true);
    }
    setRenaming(null);
  }

  return (
    <aside className="sidebar">
      <div className="sidebar-header">
        <button
          className="btn primary block"
          onClick={() => setShowNew((v) => !v)}
        >
          {t("newProject")}
        </button>
        {showNew && (
          <div className="new-project-row">
            <input
              className="text-input"
              placeholder={t("projectName")}
              value={newName}
              onChange={(e) => setNewName(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter") handleNew();
                if (e.key === "Escape") setShowNew(false);
              }}
              autoFocus
            />
            <button className="btn small" onClick={handleNew}>
              {t("create")}
            </button>
          </div>
        )}
        <div className="file-actions">
          <button className="btn block" onClick={handleExport} title={t("saveLocalTip")}>
            {t("saveLocal")}
          </button>
          <button className="btn block" onClick={handleImport} title={t("importLocalTip")}>
            {t("importLocal")}
          </button>
        </div>
      </div>
      <div className="sidebar-list">
        {projects.length === 0 && <div className="sidebar-empty">{t("noProjects")}</div>}
        {projects.map((name) => (
          <div
            key={name}
            className={`project-item ${name === currentProject ? "active" : ""}`}
            onClick={() => handleSwitch(name)}
            title={name}
          >
            {renaming === name ? (
              <input
                className="text-input rename-input"
                value={renameVal}
                autoFocus
                onClick={(e) => e.stopPropagation()}
                onChange={(e) => setRenameVal(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === "Enter") handleRename(name);
                  if (e.key === "Escape") setRenaming(null);
                }}
                onBlur={() => handleRename(name)}
              />
            ) : (
              <>
                <span className="project-name">{name}</span>
                <span className="project-ops" onClick={(e) => e.stopPropagation()}>
                  <button
                    className="project-op"
                    title={t("rename")}
                    onClick={() => {
                      setRenaming(name);
                      setRenameVal(name);
                    }}
                  >
                    ✎
                  </button>
                  <button
                    className="project-op"
                    title={t("duplicate")}
                    onClick={() => handleDuplicate(name)}
                  >
                    ⧉
                  </button>
                  <button
                    className="project-op del"
                    title={t("delete")}
                    onClick={() => handleDelete(name)}
                  >
                    ✕
                  </button>
                </span>
              </>
            )}
          </div>
        ))}
      </div>
    </aside>
  );
}
