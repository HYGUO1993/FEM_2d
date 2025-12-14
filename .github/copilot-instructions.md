**项目快速概览（已更新）**

- **构建与测试**: 仓库已添加 `CMakeLists.txt`、`tests/test_basic.cpp`、以及 CI 工作流 `.github/workflows/ci.yml`。本地可用：

```bash
# Linux/macOS
./scripts/build_and_test.sh
# Windows (PowerShell)
./scripts/build_and_test.ps1
```

- **有用的目标**:
  - `fem_run`：主程序（使用 `test05.txt` 作为默认输入）
  - `unit_tests`：小型单元测试（封装为 CTest）

- **注意**: 我在 `barsystem.cpp` 中把 `main()` 用 `#ifndef UNIT_TEST` 包起，`unit_tests` 通过定义 `UNIT_TEST` 来编译并排除 `main`，以便将函数链接到测试程序。

- **CI/自动化**: 已添加 GitHub Actions 以在 `ubuntu-latest` 上安装 `cmake`/`g++`、构建并运行 `ctest`。如果你期望其他平台（例如 MSVC）加入 CI，请告诉我。

其他内容见仓库根目录的 `README.md` 和评分注释（代码中已注出几个高优先级缺陷）。
