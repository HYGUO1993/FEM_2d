// Tauri IPC 封装（GUI_REDESIGN_PLAN §4.5）
// 命令名与参数形状必须与后端 §3 完全一致:
//   solve_model({modelJson}) / llm_chat({config, messages})
//   save_project({name, modelJson}) / list_projects() / load_project({name})
//   delete_project({name}) / get_llm_config() / set_llm_config({config})
import { t } from "./i18n.js";

const inTauri =
  typeof window !== "undefined" && typeof window.__TAURI_INTERNALS__ !== "undefined";

export function isTauri() {
  return inTauri;
}

export async function invoke(cmd, args) {
  if (!inTauri) {
    throw new Error(t("ipc.notInTauri"));
  }
  return window.__TAURI__.core.invoke(cmd, args || {});
}

export const solve = (model) => invoke("solve_model", { modelJson: JSON.stringify(model) });
export const llmChat = (config, messages) => invoke("llm_chat", { config, messages });
export const llmChatTools = (config, messages, tools) =>
  invoke("llm_chat_tools", { config, messages, tools });
export const saveProject = (name, model) =>
  invoke("save_project", { name, modelJson: JSON.stringify(model) });
export const listProjects = () => invoke("list_projects");
export const loadProject = (name) => invoke("load_project", { name });
export const deleteProject = (name) => invoke("delete_project", { name });
export const getLlmConfig = () => invoke("get_llm_config");
export const setLlmConfig = (config) => invoke("set_llm_config", { config });
export const exportModelFile = (defaultName, modelJson) =>
  invoke("export_model_file", { defaultName, modelJson });
export const importModelFile = () => invoke("import_model_file");
