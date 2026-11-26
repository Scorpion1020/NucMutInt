# NucMutInt

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C.svg)](https://en.cppreference.com/w/cpp/17)

**English** | [简体中文](README_CN.md)

**NucMutInt** is a lightweight C++ command-line tool for **in-silico nucleic acid base mutation** and **protein–nucleic acid interaction difference analysis**.

Given a protein–DNA/RNA complex in PDB format, NucMutInt replaces specified nucleotide bases in place while preserving the sugar-phosphate backbone, then invokes [PLIP](https://github.com/pharmai/plip) to compare non-covalent interactions before and after mutation and exports a structured Markdown report.

## Features

- **In-place DNA/RNA base mutation**: Replaces target bases (`A`, `C`, `G`, `T`, `U`) using rigid-body geometric alignment:
  1. Aligns the template base plane normal ($n_2$) to the original base plane normal ($n_1$) via rotation $R_1$.
  2. Projects the glycosidic bond vectors onto the base plane and applies an in-plane rotation $R_2$ around $n_1$.
  3. Translates the rotated base to the original `C1'` anchor while keeping backbone atoms intact.
- **Multi-chain & multi-region support**: Mutate single positions or contiguous segments across one or more chains in a single command.
- **Non-standard `HETATM` residue handling**: Supports mutating modified nucleotides (prompts for glycosidic alignment atom names) and automatically handles temporary chain duplication so PLIP can analyze modified residues without polluting the final diff.
- **Automated interaction diff**: Parses PLIP XML reports and summarizes lost/gained interactions across 6 categories:
  - Hydrogen bonds (`hydrogen_bonds`)
  - $\pi$-stacking (`pi_stacks`)
  - $\pi$-cation interactions (`pi_cation_interactions`)
  - Salt bridges (`salt_bridges`)
  - Hydrophobic interactions (`hydrophobic_interactions`)
  - Water bridges (`water_bridges`)
- **Markdown report export**: Outputs both console summary and `<mutated_pdb>_outputFolder/report.md`.

## Directory Structure

```text
NucMutInt/
├── CMakeLists.txt
├── LICENSE
├── Makefile
├── README.md
├── data/
│   └── nucleicacid.dat        # Standard base heavy-atom & C1' coordinates (A/C/G/T/U)
├── examples/
│   ├── 1emh.pdb               # Public RCSB PDB example (Uracil-DNA glycosylase + DNA)
│   └── sample_report.md       # Sample interaction diff report
├── src/
│   ├── HelixFromTmpl.cpp      # PDB parser, geometric base replacement, CLI entry point
│   └── PNI_Calculate.hpp      # PLIP runner, XML parser, interaction diff & Markdown writer
└── third_party/
    └── tinyxml2/              # Bundled TinyXML-2 (zlib license)
```

## Dependencies

- **C++17 compiler** (`g++ >= 8` or `clang++ >= 7`)
- **Eigen 3** (`>= 3.3`, header-only linear algebra library)
- **PLIP** (external runtime tool on `PATH`, required for interaction analysis; base mutation still runs without PLIP)

### Installing Dependencies (Ubuntu / Conda)

```bash
# Ubuntu / Debian
sudo apt-get install g++ cmake libeigen3-dev

# PLIP via Conda (recommended)
conda create -n plip python=3.9
conda activate plip
conda install -c conda-forge openbabel plip
```

## Build

### Option 1: CMake (Recommended)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### Option 2: Makefile

```bash
make EIGEN_INC=/usr/include/eigen3
```

## Usage

```bash
./build/NucMutInt <pdb_file> (<chainID> <resID_begin> <resID_end> <sequence>) ...
```

### Examples

Mutate residues 2–6 on chain `B` of `examples/1emh.pdb` to `ACGAT`:

```bash
./build/NucMutInt examples/1emh.pdb B 2 6 ACGAT
```

Mutate multiple chains in one run:

```bash
./build/NucMutInt complex.pdb A 1 8 ACGTACGT B 10 12 CGT
```

### Outputs

- `<stem>_<chain>_<begin>_<end>_<seq>.pdb`: Mutated complex structure.
- `<mutated_pdb>_outputFolder/original/report.xml`: PLIP report for the wild-type structure.
- `<mutated_pdb>_outputFolder/mutated/report.xml`: PLIP report for the mutated structure.
- `<mutated_pdb>_outputFolder/report.md`: Markdown summary of lost and gained interactions (see [examples/sample_report.md](examples/sample_report.md)).

## License & Third-Party Notices

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

### Third-Party Dependencies & Acknowledgments

- **PLIP (Protein-Ligand Interaction Profiler)** — Licensed under **GPL-2.0**. PLIP is invoked as an external command-line dependency and is not distributed in this repository. Please install PLIP separately and cite:
  > Adasme MF, Linnemann KL, Bolz SN, Kaiser F, Salentin S, Haupt VJ, Schroeder M. *PLIP 2021: expanding the scope of the protein-ligand interaction profiler to DNA and RNA.* Nucleic Acids Research, 2021, 49(W1): W530–W534.
- **Eigen** — Licensed under **MPL-2.0**. Required at build time as an external header-only dependency.
- **TinyXML-2** — Licensed under the **zlib License**. Bundled under `third_party/tinyxml2/` with its original license notice preserved.
- **Example Structure (`examples/1emh.pdb`)** — Publicly available crystal structure from the [RCSB Protein Data Bank](https://www.rcsb.org/structure/1EMH) (PDB ID: `1EMH`).
