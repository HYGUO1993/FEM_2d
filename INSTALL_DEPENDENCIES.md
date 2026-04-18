# 依赖安装故障排除

如果启动脚本无法自动安装依赖，请按照本指南手动安装。

---

## 问题：启动脚本提示"安装失败"

### 原因可能是：
1. 网络连接问题
2. pip 权限不足
3. Python 环境配置问题

---

## 解决方案

### 方案 1：使用 `--user` 标志（最可能成功！）

```bash
python -m pip install --user PySide6 matplotlib numpy
```

**适用于：** 大多数情况，特别是权限受限时

---

### 方案 2：使用 conda（如果已安装 Anaconda）

```bash
conda install -c conda-forge PySide6 matplotlib numpy
```

**适用于：** Anaconda/Miniconda 用户

---

### 方案 3：一个包一个包安装

如果上述方案失败，尝试逐个安装：

```bash
# 首先升级 pip
python -m pip install --upgrade pip

# 然后逐个安装
python -m pip install --user PySide6
python -m pip install --user matplotlib
python -m pip install --user numpy
```

---

### 方案 4：检查网络后重试

```bash
# 先测试网络连接
ping pypi.org

# 然后使用不同的 PyPI 源
python -m pip install -i https://pypi.org/simple PySide6 matplotlib numpy
```

---

### 方案 5：使用国内镜像（中国用户）

```bash
# 阿里云镜像
python -m pip install -i https://mirrors.aliyun.com/pypi/simple PySide6 matplotlib numpy

# 或清华大学镜像
python -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple PySide6 matplotlib numpy
```

---

## 验证安装

安装后验证是否成功：

```bash
# 逐个测试
python -c "import PySide6; print('PySide6: OK')"
python -c "import matplotlib; print('matplotlib: OK')"
python -c "import numpy; print('numpy: OK')"

# 或一次全部测试
python -c "import PySide6, matplotlib, numpy; print('All OK!')"
```

---

## 如果仍然失败

### 检查 Python 环境

```bash
# 查看 Python 位置
python --version
which python          # Linux/macOS
where python          # Windows

# 查看 pip 位置
pip --version
which pip             # Linux/macOS
where pip             # Windows
```

### 检查 pip 配置

```bash
# 查看 pip 配置
pip config list

# 重置为默认
pip config unset global.index-url
```

### 清理缓存后重试

```bash
python -m pip cache purge
python -m pip install --user --no-cache-dir PySide6 matplotlib numpy
```

---

## 备选方案：不安装依赖直接运行

如果实在无法安装依赖，可以使用命令行工具：

```bash
# 使用命令行求解（无 GUI）
./build/bin/Release/fem_run.exe --input test05.txt --output build/Results.dat

# 生成可视化（无 GUI）
python scripts/visualize_results.py --results build/Results.dat --out build/plot.png
```

---

## 获取帮助

如果以上方案都不奏效：

1. 检查 Python 是否安装正确（`python --version` 应显示 3.8+）
2. 尝试在不同的命令行工具中运行（PowerShell、CMD、git bash）
3. 查看完整的 GUI 使用指南：`GUI_USAGE.md`
4. 查看项目 README：`README.md`

---

## 常见错误及解决

### 错误：`No module named 'pip'`

```bash
# 重新安装 pip
python -m ensurepip --upgrade
```

### 错误：`permission denied` / `Access is denied`

使用 `--user` 标志：
```bash
python -m pip install --user PySide6 matplotlib numpy
```

### 错误：`Collecting ... ERROR: Could not find a version`

网络问题或包名错误，尝试不同的源。

### 错误：`SSL: CERTIFICATE_VERIFY_FAILED`

```bash
# 临时忽略 SSL 验证
python -m pip install --trusted-host pypi.python.org --user PySide6 matplotlib numpy
```

---

**最后更新：** 2026年4月18日
