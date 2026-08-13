//! FemLab Studio Tauri 后端
//! 提供 `solve_model` 命令: 接收 JSON 模型 → 调 femcli 求解 → 返回 JSON 结果。

use std::process::Command;
use serde::{Deserialize, Serialize};
use tauri::Manager;

/// 求解模型
/// `model_json`: 符合 FEM_MODEL_SCHEMA 的输入 JSON 字符串
/// 返回: femcli 输出的结果 JSON 字符串
#[tauri::command]
fn solve_model(model_json: String, app: tauri::AppHandle) -> Result<String, String> {
    let femcli_path = locate_femcli(&app)?;

    // 把模型 JSON 写到临时文件,再调 femcli
    let tmp_dir = std::env::temp_dir().join("fem2d_tauri");
    std::fs::create_dir_all(&tmp_dir).map_err(|e| format!("创建临时目录失败: {e}"))?;
    let input_path = tmp_dir.join("model.json");
    let output_path = tmp_dir.join("result.json");
    std::fs::write(&input_path, model_json).map_err(|e| format!("写模型文件失败: {e}"))?;

    // 调 femcli: solve <input> -o <output>
    let output = Command::new(&femcli_path)
        .arg("solve")
        .arg(&input_path)
        .arg("-o")
        .arg(&output_path)
        .output()
        .map_err(|e| format!("无法启动 femcli: {e}"))?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(format!("femcli 求解失败 (exit={}): {}", output.status.code().unwrap_or(-1), stderr.trim()));
    }

    let result = std::fs::read_to_string(&output_path)
        .map_err(|e| format!("读取结果失败: {e}"))?;
    Ok(result)
}

/// 定位 femcli.exe: 优先随应用打包,其次向上遍历父目录找项目根 build/bin
fn locate_femcli(app: &tauri::AppHandle) -> Result<std::path::PathBuf, String> {
    // 1. 打包形态: exe 同目录
    if let Some(exe_dir) = std::env::current_exe().ok().and_then(|p| p.parent().map(|p| p.to_path_buf())) {
        for cand in ["femcli.exe", "femcli"] {
            let p = exe_dir.join(cand);
            if p.exists() {
                return Ok(p);
            }
        }
        // 2. 开发形态: 从 exe 目录向上遍历, 找 <项目根>/build/bin/femcli.exe
        //    例: ...\FEM_2d\gui_tauri\src-tauri\target\debug\femlab-studio.exe
        //        → 逐级上溯直到 FEM_2d\build\bin\femcli.exe
        let mut dir = exe_dir.as_path().parent();
        while let Some(d) = dir {
            for cand in ["build/bin/femcli.exe", "build/bin/femcli"] {
                let p = d.join(cand);
                if p.exists() {
                    return Ok(p);
                }
            }
            dir = d.parent();
        }
    }
    // 3. resource_dir / 当前工作目录
    let resource_dir = app.path().resource_dir().unwrap_or_default();
    for base in [resource_dir, std::env::current_dir().unwrap_or_default()] {
        for cand in ["build/bin/femcli.exe", "build/bin/femcli"] {
            let p = base.join(cand);
            if p.exists() {
                return Ok(p);
            }
        }
    }
    Err("未找到 femcli.exe。请先在 FEM_2d 项目根目录构建求解器 (build/bin/femcli.exe)".into())
}

/// 前端 ping 用,验证 IPC 正常
#[tauri::command]
fn ping() -> String {
    "pong".into()
}

/// LLM 配置(OpenAI 兼容 API)
#[derive(Serialize, Deserialize, Clone, Default)]
struct LlmConfig {
    #[serde(default)] base_url: String,   // 例 "https://api.deepseek.com/v1"
    #[serde(default)] api_key: String,
    #[serde(default)] model: String,      // 例 "deepseek-chat"
}

/// 对话消息
#[derive(Serialize, Deserialize)]
struct ChatMessage {
    role: String,     // "system" | "user" | "assistant"
    content: String,
}

/// 获取应用数据目录
fn app_data_dir(app: &tauri::AppHandle) -> Result<std::path::PathBuf, String> {
    app.path()
        .app_data_dir()
        .map_err(|e| format!("获取应用数据目录失败: {e}"))
}

/// 获取项目存储目录(<app_data_dir>/projects)
fn projects_dir(app: &tauri::AppHandle) -> Result<std::path::PathBuf, String> {
    Ok(app_data_dir(app)?.join("projects"))
}

/// 清洗项目名,替换路径分隔符等危险字符,防止路径穿越
fn sanitize_name(name: &str) -> String {
    name.chars()
        .map(|c| match c {
            '\\' | '/' | ':' | '*' | '?' | '"' | '<' | '>' | '|' => '_',
            _ => c,
        })
        .collect()
}

/// LLM 对话: POST {base_url}/chat/completions,返回 choices[0].message.content
#[tauri::command]
fn llm_chat(config: LlmConfig, messages: Vec<ChatMessage>) -> Result<String, String> {
    let url = format!("{}/chat/completions", config.base_url.trim_end_matches('/'));
    let body = serde_json::json!({
        "model": config.model,
        "messages": messages,
        "temperature": 0.2,
        "stream": false,
    });

    let response = ureq::post(&url)
        .set("Authorization", &format!("Bearer {}", config.api_key))
        .set("Content-Type", "application/json")
        .send_json(&body)
        .map_err(|e| format!("LLM 请求失败: {e}"))?;

    let value: serde_json::Value = response
        .into_json()
        .map_err(|e| format!("LLM 响应解析失败: {e}"))?;

    value
        .pointer("/choices/0/message/content")
        .and_then(|c| c.as_str())
        .map(|s| s.to_string())
        .ok_or_else(|| "LLM 响应缺少 choices[0].message.content".to_string())
}

/// 保存项目: 写 <app_data_dir>/projects/<name>.json
#[tauri::command]
fn save_project(app: tauri::AppHandle, name: String, model_json: String) -> Result<(), String> {
    let safe_name = sanitize_name(&name);
    let dir = projects_dir(&app)?;
    std::fs::create_dir_all(&dir).map_err(|e| format!("创建项目目录失败: {e}"))?;
    let path = dir.join(format!("{safe_name}.json"));
    std::fs::write(&path, model_json).map_err(|e| format!("保存项目失败: {e}"))
}

/// 列出项目: 返回 <app_data_dir>/projects/ 下 *.json 文件名(去后缀,按名排序)
#[tauri::command]
fn list_projects(app: tauri::AppHandle) -> Result<Vec<String>, String> {
    let dir = projects_dir(&app)?;
    if !dir.exists() {
        return Ok(Vec::new());
    }
    let mut names: Vec<String> = Vec::new();
    let entries = std::fs::read_dir(&dir).map_err(|e| format!("读取项目目录失败: {e}"))?;
    for entry in entries.flatten() {
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) == Some("json") {
            if let Some(stem) = path.file_stem().and_then(|s| s.to_str()) {
                names.push(stem.to_string());
            }
        }
    }
    names.sort();
    Ok(names)
}

/// 加载项目: 读 <app_data_dir>/projects/<name>.json 返回原文
#[tauri::command]
fn load_project(app: tauri::AppHandle, name: String) -> Result<String, String> {
    let safe_name = sanitize_name(&name);
    let path = projects_dir(&app)?.join(format!("{safe_name}.json"));
    std::fs::read_to_string(&path).map_err(|e| format!("读取项目失败: {e}"))
}

/// 删除项目: 删对应文件; 不存在则 Ok(幂等)
#[tauri::command]
fn delete_project(app: tauri::AppHandle, name: String) -> Result<(), String> {
    let safe_name = sanitize_name(&name);
    let path = projects_dir(&app)?.join(format!("{safe_name}.json"));
    match std::fs::remove_file(&path) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(format!("删除项目失败: {e}")),
    }
}

/// 读取 LLM 配置; 不存在返回默认(空)
#[tauri::command]
fn get_llm_config(app: tauri::AppHandle) -> Result<LlmConfig, String> {
    let path = app_data_dir(&app)?.join("llm_config.json");
    if !path.exists() {
        return Ok(LlmConfig::default());
    }
    let content = std::fs::read_to_string(&path).map_err(|e| format!("读取 LLM 配置失败: {e}"))?;
    serde_json::from_str(&content).map_err(|e| format!("解析 LLM 配置失败: {e}"))
}

/// 保存 LLM 配置到 <app_data_dir>/llm_config.json
#[tauri::command]
fn set_llm_config(app: tauri::AppHandle, config: LlmConfig) -> Result<(), String> {
    let dir = app_data_dir(&app)?;
    std::fs::create_dir_all(&dir).map_err(|e| format!("创建配置目录失败: {e}"))?;
    let path = dir.join("llm_config.json");
    let content =
        serde_json::to_string_pretty(&config).map_err(|e| format!("序列化 LLM 配置失败: {e}"))?;
    std::fs::write(&path, content).map_err(|e| format!("保存 LLM 配置失败: {e}"))
}

/// 导出模型到本地文件: 弹出保存对话框, 写 model_json, 返回保存路径
#[tauri::command]
fn export_model_file(default_name: String, model_json: String) -> Result<String, String> {
    use std::io::Write;
    let mut default_name = default_name.trim().to_string();
    if default_name.is_empty() {
        default_name = "model.json".to_string();
    }
    if !default_name.to_lowercase().ends_with(".json") {
        default_name.push_str(".json");
    }
    let path = rfd::FileDialog::new()
        .set_file_name(&default_name)
        .add_filter("JSON 模型", &["json"])
        .save_file()
        .ok_or_else(|| "用户取消了保存".to_string())?;
    let mut f = std::fs::File::create(&path).map_err(|e| format!("创建文件失败: {e}"))?;
    f.write_all(model_json.as_bytes())
        .map_err(|e| format!("写入文件失败: {e}"))?;
    Ok(path.display().to_string())
}

/// 导入模型文件: 弹出打开对话框, 读取文件内容返回字符串
#[tauri::command]
fn import_model_file() -> Result<String, String> {
    let path = rfd::FileDialog::new()
        .add_filter("JSON 模型", &["json"])
        .pick_file()
        .ok_or_else(|| "用户取消了选择".to_string())?;
    std::fs::read_to_string(&path).map_err(|e| format!("读取文件失败: {e}"))
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![
            solve_model, ping,
            llm_chat, save_project, list_projects, load_project, delete_project,
            get_llm_config, set_llm_config,
            export_model_file, import_model_file
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
