# FEM_2d v1.0 - 快速开始指南

## 立即使用（3步）

### 1️⃣ 克隆并切换到发布分支
```bash
git clone https://github.com/HYGUO1993/FEM_2d.git
cd FEM_2d
git checkout claude/improve-visualization-for-engineers
```

### 2️⃣ 编译
```bash
mkdir build && cd build
cmake ..
make
```

### 3️⃣ 运行示例
```bash
# 从项目根目录运行
./build/bin/fem_run --input test05.txt --output build/Results.dat

# 查看结果
cat build/Results.dat
```

## 输出结果示例

```
Node Displacements:
    Node          Ux          Uy          Rz
       0           0           0           0
       1           0-0.000190484-0.000285726

Element End Forces:
    Elem        Fx_i        Fy_i        Mz_i        Fx_j        Fy_j        Mz_j
       0           0        1000        1000           0       -1000           0

Support Reactions:
    Node          Rx          Ry          Rz
       0           0        1000        1000
```

## 可视化（可选）

```bash
# 安装依赖
pip install matplotlib numpy

# 生成图片
python scripts/visualize_results.py --results build/Results.dat --scale 1000 --out plot.png
```

## 完整文档

- **使用手册**: [USAGE_GUIDE.md](USAGE_GUIDE.md)
- **发布说明**: [RELEASE_NOTES.md](RELEASE_NOTES.md)
- **开发文档**: [DEVELOPMENT.md](DEVELOPMENT.md)
- **项目说明**: [README.md](README.md)

## ✅ 验证状态

- ✅ 编译成功
- ✅ 单元测试通过 (ctest: 100% passed)
- ✅ 悬臂梁示例运行正常
- ✅ 结果输出正确

**版本**: v1.0  
**分支**: claude/improve-visualization-for-engineers  
**日期**: 2026-04-18  
**测试平台**: Linux (Ubuntu latest)
