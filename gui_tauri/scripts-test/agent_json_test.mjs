// parseModelJson 容错单测
import { parseModelJson } from "../src/agent.js";

let fails = 0;
function check(name, cond, detail) {
  if (!cond) fails++;
  console.log(`[${cond ? "PASS" : "FAIL"}] ${name}${detail ? "  " + detail : ""}`);
}

// 1. 正常 JSON
check("正常 JSON", parseModelJson('{"a":1,"b":[1,2]}')?.a === 1);
// 2. markdown 包裹
check("markdown 代码块", parseModelJson('```json\n{"a":1}\n```')?.a === 1);
// 3. 尾逗号
check("尾逗号", parseModelJson('{"a":1,"b":[1,2,]}')?.b?.length === 2);
// 4. 缺右括号
check("缺右括号", parseModelJson('{"a":{"b":1}')?.a?.b === 1);
// 5. 前后杂音
check("前后杂音", parseModelJson('好的, 模型如下: {"a":1} 以上')?.a === 1);
// 6. 截断 (缺多个括号)
check("截断补括号", parseModelJson('{"a":{"b":{"c":1')?.a?.b?.c === 1);
// 7. 完全非法
check("完全非法返回 null", parseModelJson('not json at all') === null);
// 8. 空串
check("空串返回 null", parseModelJson('') === null);

console.log(fails === 0 ? "全部通过!" : `${fails} 项失败!`);
process.exit(fails === 0 ? 0 : 1);
