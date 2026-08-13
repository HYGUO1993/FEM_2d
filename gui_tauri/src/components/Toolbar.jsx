import { useStore } from "../store.js";
import { t } from "../i18n.js";

const TOOLS = [
  { id: "select", label: "tool.select", icon: "◈" },
  { id: "node", label: "tool.node", icon: "●" },
  { id: "element", label: "tool.element", icon: "╱" },
  { id: "load", label: "tool.load", icon: "↑" },
  { id: "constraint", label: "tool.constraint", icon: "▽" },
  { id: "erase", label: "tool.erase", icon: "✕" },
];

export default function Toolbar() {
  const tool = useStore((s) => s.tool);
  const setTool = useStore((s) => s.setTool);

  return (
    <div className="toolbar">
      {TOOLS.map((tb) => (
        <button
          key={tb.id}
          className={`tool-btn ${tool === tb.id ? "active" : ""}`}
          onClick={() => setTool(tb.id)}
          title={t(tb.label)}
        >
          <span className="tool-icon">{tb.icon}</span>
          <span className="tool-label">{t(tb.label)}</span>
        </button>
      ))}
    </div>
  );
}
