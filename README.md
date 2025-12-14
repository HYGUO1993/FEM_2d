# FEM_2d

一个用于教学/实验的二维杆/桁架/刚架有限元程序。

Quick start

```bash
# Linux / macOS
./scripts/build_and_test.sh

# Windows (PowerShell)
./scripts/build_and_test.ps1
```

Notes

- The program reads `test05.txt` as input and writes `Results.dat` as output by default.
- There is a minimal test suite under `tests/` and CI configured under `.github/workflows/ci.yml`.
