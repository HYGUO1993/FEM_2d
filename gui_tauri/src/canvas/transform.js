// 画布坐标变换 + pan/zoom 数学（GUI_REDESIGN_PLAN §4.6）
// 模型坐标 y 向上, 屏幕 y 向下 → worldToScreen 翻转 y。
// view = { scale, panX, panY }
//   sx = x*scale + panX
//   sy = -y*scale + panY

export function worldToScreen(x, y, view) {
  if (!view) return { x, y: -y };
  return { x: x * view.scale + view.panX, y: -y * view.scale + view.panY };
}

export function screenToWorld(sx, sy, view) {
  if (!view) return { x: sx, y: -sy };
  return { x: (sx - view.panX) / view.scale, y: (view.panY - sy) / view.scale };
}

export function dist(x1, y1, x2, y2) {
  return Math.hypot(x2 - x1, y2 - y1);
}

/** 点到线段距离（屏幕坐标） */
export function pointSegDist(px, py, x1, y1, x2, y2) {
  const dx = x2 - x1;
  const dy = y2 - y1;
  const len2 = dx * dx + dy * dy;
  let t = len2 === 0 ? 0 : ((px - x1) * dx + (py - y1) * dy) / len2;
  t = Math.max(0, Math.min(1, t));
  return dist(px, py, x1 + t * dx, y1 + t * dy);
}

export function clamp(v, lo, hi) {
  return Math.max(lo, Math.min(hi, v));
}
