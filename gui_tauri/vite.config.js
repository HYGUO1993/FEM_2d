import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [
    react(),
    // Tauri custom-protocol 下不允许 crossorigin，否则 WebView2 拒绝加载 JS 模块 → 黑屏
    { name: "strip-crossorigin", transformIndexHtml(html) { return html.replace(/ crossorigin/g, ""); } },
  ],
  clearScreen: false,
  // 关键: Tauri custom-protocol 下必须用相对路径, 否则 /assets/xxx.js 加载不到 → 黑屏
  base: "./",
  server: { port: 5173, strictPort: true },
  build: {
    outDir: "dist", target: "chrome105", sourcemap: false,
    // Tauri custom-protocol 下 module preload + crossorigin 会导致资源加载失败
    modulePreload: { polyfill: false },
  },
});
