# NucMutInt

[![CI](https://github.com/Scorpion1020/NucMutInt/actions/workflows/ci.yml/badge.svg)](https://github.com/Scorpion1020/NucMutInt/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C.svg)](https://en.cppreference.com/w/cpp/17)

[English](README.md) | **简体中文**

**NucMutInt** 是一个轻量级 C++ 命令行工具，用于**核酸原位碱基突变（in-silico mutation）**与**蛋白–核酸相互作用差异分析**。

输入蛋白–DNA/RNA 复合物的 PDB 结构后，NucMutInt 会在保留糖磷酸骨架的前提下原位替换指定碱基，并自动调用 [PLIP](https://github.com/pharmai/plip) 对比突变前后的非共价相互作用变化，导出结构化的 Markdown 分析报告。

## 主要特性

- **原位 DNA/RNA 碱基替换**：基于刚体几何变换替换目标碱基（`A`、`C`、`G`、`T`、`U`）：
  1. 计算原始碱基平面法向量 $n_1$ 与模板碱基平面法向量 $n_2$，通过旋转矩阵 $R_1$ 对齐碱基平面；
  2. 将糖苷键向量投影至碱基平面，绕法向量 $n_1$ 旋转 $R_2$ 对齐糖苷键方向；
  3. 以原 `C1'` 原子为锚点平移新碱基，完整保留糖磷酸骨架原子。
- **多链与多区域支持**：单条命令即可同时对一条或多条核酸链的单个/连续位点进行突变。
- **非标准 `HETATM` 残基支持**：支持修饰核苷酸突变（交互指定糖苷键对齐原子），并自动构建临时链以兼容 PLIP 对非标准残基的解析，最终报告中自动过滤临时链冗余项。
- **相互作用自动对比**：解析 PLIP 生成的 XML 报告，按以下 6 类相互作用统计缺失（Lost）与新增（Gained）的相互作用：
  - 氢键（`hydrogen_bonds`）
  - $\pi$-堆积（`pi_stacks`）
  - $\pi$-阳离子相互作用（`pi_cation_interactions`）
  - 盐桥（`salt_bridges`）
  - 疏水相互作用（`hydrophobic_interactions`）
  - 水桥（`water_bridges`）
- **空间位阻（Steric Clash）检测**：突变完成后自动基于 Bondi 范德华半径扫描新碱基重原子与周边蛋白/核酸原子的距离（$d < 0.70 \times (R_1 + R_2)$），在终端和报告中提示潜在碰撞。
- **Markdown 表格化报告**：同时在终端输出对比结果，并生成 `<mutated_pdb>_outputFolder/report.md`，以表格形式展示相互作用距离、角度及位阻检测结果。

## 目录结构

```text
NucMutInt/
├── CMakeLists.txt
├── LICENSE
├── Makefile
├── README.md                  # 英文说明
├── README_CN.md               # 中文说明
├── data/
│   └── nucleicacid.dat        # 标准碱基重原子及 C1' 模板坐标（A/C/G/T/U）
├── examples/
│   ├── 1emh.pdb               # 公开 RCSB PDB 测试样例（Uracil-DNA 糖基化酶 + DNA）
│   └── sample_report.md       # 相互作用差异报告样例
├── src/
│   ├── HelixFromTmpl.cpp      # PDB 解析、刚体碱基替换、命令行入口
│   └── PNI_Calculate.hpp      # PLIP 调用、XML 解析、差异对比与 Markdown 导出
└── third_party/
    └── tinyxml2/              # 内置 TinyXML-2（zlib 协议）
```

## 环境依赖

- **C++17 编译器**（`g++ >= 8` 或 `clang++ >= 7`）
- **Eigen 3**（`>= 3.3`，仅头文件线性代数库）
- **PLIP**（外部运行时工具，需在 `PATH` 中可用；未安装 PLIP 时仍可单独完成碱基突变并输出 PDB）

### 安装依赖（Ubuntu / Conda）

```bash
# Ubuntu / Debian
sudo apt-get install g++ cmake libeigen3-dev

# 通过 Conda 安装 PLIP（推荐）
conda create -n plip python=3.9
conda activate plip
conda install -c conda-forge openbabel plip
```

## 编译方式

### 方式一：CMake（推荐）

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### 方式二：Makefile

```bash
make EIGEN_INC=/usr/include/eigen3
```

## 使用方法

```bash
./build/NucMutInt <pdb_file> (<chainID> <resID_begin> <resID_end> <sequence>) ...
```

### 运行示例

将 `examples/1emh.pdb` 中 `B` 链第 2–6 位碱基突变为 `ACGAT`：

```bash
./build/NucMutInt examples/1emh.pdb B 2 6 ACGAT
```

单次命令同时突变多条链：

```bash
./build/NucMutInt complex.pdb A 1 8 ACGTACGT B 10 12 CGT
```

### 输出文件

- `<stem>_<chain>_<begin>_<end>_<seq>.pdb`：突变后的复合物 PDB 结构。
- `<mutated_pdb>_outputFolder/original/report.xml`：野生型结构的 PLIP 分析结果。
- `<mutated_pdb>_outputFolder/mutated/report.xml`：突变型结构的 PLIP 分析结果。
- `<mutated_pdb>_outputFolder/report.md`：相互作用变化汇总报告（样例见 [examples/sample_report.md](examples/sample_report.md)）。

## 开源协议与第三方声明

本项目采用 **MIT License** 开源，详见 [LICENSE](LICENSE)。

### 第三方依赖与致谢

- **PLIP (Protein-Ligand Interaction Profiler)** — 采用 **GPL-2.0** 协议。本项目仅通过命令行调用用户本地安装的 PLIP，不包含或分发 PLIP 源码。使用 PLIP 请参考并引用：
  > Adasme MF, Linnemann KL, Bolz SN, Kaiser F, Salentin S, Haupt VJ, Schroeder M. *PLIP 2021: expanding the scope of the protein-ligand interaction profiler to DNA and RNA.* Nucleic Acids Research, 2021, 49(W1): W530–W534.
- **Eigen** — 采用 **MPL-2.0** 协议，作为外部头文件构建依赖引入。
- **TinyXML-2** — 采用 **zlib License**，内置于 `third_party/tinyxml2/` 并保留原始版权声明。
- **测试结构（`examples/1emh.pdb`）** — 来自公开数据库 [RCSB Protein Data Bank](https://www.rcsb.org/structure/1EMH)（PDB ID: `1EMH`）。
