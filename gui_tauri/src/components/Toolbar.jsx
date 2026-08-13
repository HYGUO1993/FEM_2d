import { useStore } from "../store.js";

const TOOLS = [
  { id: "select", label: "选择", icon: "◈" },
  { id: "node", label: "节点", icon: "●" },
  { id: "element", label: "杆件", icon: "╱" },
  { id: "load", label: "荷载", icon: "↑" },
  { id: "constraint", label: "约束", icon: "▽" },
  { id: "erase", label: "删除", icon: "✕" },
];

export default function Toolbar() {
  const tool = useStore((s) => s.tool);
  const setTool = useStore((s) => s.setTool);

  return (
    <div className="toolbar">
      {TOOLS.map((t) => (
        <button
          key={t.id}
          className={`tool-btn ${tool === t.id ? "active" : ""}`}
          onClick={() => setTool(t.id)}
          title={t.label}
        >
          <span className="tool-icon">{t.icon}</span>
          <span className="tool-label">{t.label}</span>
        </button>
      ))}
    </div>
  );
}
