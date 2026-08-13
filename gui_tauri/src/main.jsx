import React from "react";
import ReactDOM from "react-dom/client";
import App from "./App.jsx";
import "./styles.css";
import { useStore } from "./store.js";
import { setLang } from "./i18n.js";

// 应用持久化的主题与语言 (store 初始化已读取 localStorage)
const st0 = useStore.getState();
document.documentElement.setAttribute("data-theme", st0.theme || "dark");
setLang(st0.lang || "zh");

ReactDOM.createRoot(document.getElementById("root")).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
);
