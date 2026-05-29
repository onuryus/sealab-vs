# sealab-vs

**Ligand-based virtual screening engine in C++.** Given a target SMILES, screens a large `.smi` library by fingerprint similarity, layered with drug-likeness rules, structural alerts, QED scoring, BOILED-Egg permeability, and a two-stage workflow for million-scale libraries. Multicore via OpenMP.

A Sealab Project.

---

## Table of Contents

- [What is sealab-vs?](#what-is-sealab-vs)
- [Features](#features)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Parameter reference](#parameter-reference)
- [Screening presets](#screening-presets)
- [Output format](#output-format)
- [Performance](#performance)
- [Documentation](#documentation)
- [License](#license)

---

## What is sealab-vs?

In structure-based drug discovery, **docking** is the expensive step — it can take seconds to minutes per molecule, and a typical library has millions. You don't want to dock everything; you want to dock the *interesting subset*.

**sealab-vs is that pre-filter.** It takes:

1. A **query molecule** (your target — a known active, a fragment, a reference compound).
2. A **library of candidate molecules** in `.smi` format (ZINC, ChEMBL, your in-house collection, …).

…and produces a **ranked CSV** of library molecules sorted by similarity to the query, with every drug-likeness and toxicity check pre-computed as a column. You decide which subset goes into docking based on that CSV.

It is **standalone**. No GUI, no database, no network. Input: one file + one SMILES. Output: one file. Stream through it.

### Why C++?

The original prototype was a single-threaded Python script and could not keep up with modern library sizes. The C++ rewrite uses RDKit C++, OpenMP, and a custom streaming reader to process ~2 000 mol/s on 16 cores for the full pipeline, or ~60 000 mol/s for a fast prefilter pass — about 10–100× the Python version, with deterministic output.

---

## Features

### Fingerprint similarity
- **Methods**: Morgan (ECFP-like), AtomPair, RDKit path-based
- **Similarity coefficients**: Tanimoto, Dice, Sokal, McConnaughey, Kulczynski, Cosine
- **Configurable** bit length and (Morgan) radius

### Drug-likeness rules (boolean gates)
Lipinski • Veber • Ghose • Egan • Muegge • Rule-of-Three • Lead-like • REOS • GSK 4/400 • Pfizer 3/75

### Continuous-score gates
- **QED** (Bickerton 2012) — single composite drug-likeness score
- **Fsp3** threshold (Lovering 2009) — escape-from-flatland
- **Stereocenter** count cap

### Structural alert catalogs
- **PAINS** — pan-assay interference (Baell & Holloway 2010)
- **BRENK** — toxic/reactive groups (Brenk 2008)
- **Lilly MedChem Rules** (Bruns & Watson 2012)
- **Inpharmatica** alerts (Sutherland 2003)
- **Ames** mutagenicity — Benigni–Bossa SMARTS subset
- **Aggregators** — Shoichet-style SMARTS subset
- **NIH / ZINC / CHEMBL** extra bundle

### Permeability
- **BOILED-Egg** GIA / BBB regions (Daina & Zoete 2016)
- **BBB-MPO** score threshold (Wager 2010)

### Performance
- **OpenMP** parallel processing — scales linearly up to ~6.5× on 16 cores
- **Streaming reader** — never loads the whole library into RAM
- **Two-stage screening** — fast prefilter pass for million-scale libraries
- **Live progress bar** with rate and ETA

---

## Installation

A single conda environment carries everything. Full step-by-step in [`SETUP.txt`](SETUP.txt) (bilingual EN/TR). TL;DR:

```bash
conda create -n sealab-vs -c conda-forge \
    librdkit-dev boost-cpp cmake cxx-compiler make
conda activate sealab-vs

git clone https://github.com/<user>/sealab-vs.git
cd sealab-vs
cmake -S . -B build
cmake --build build -j
```

You get `build/sealab-vs`. The conda lib path is baked in as rpath, so the binary runs even without the env activated.

---

## Quick start

```bash
./build/sealab-vs \
    --lib examples/example.smi \
    --target 'CC(=O)Oc1ccccc1C(=O)O' \
    -o results.csv
```

This screens `example.smi` for similarity to **aspirin**, using defaults (Morgan / Tanimoto / 2048 bits / radius 2 / Lipinski rule / all CPU cores). Results go to `results.csv`, sorted by descending similarity.

For the full option list at any time:
```bash
./build/sealab-vs --help
```

---

## Parameter reference

Every flag with what it does and a worked example.

### Required input

#### `--lib <path>`
Input library file in `.smi` format. Each line: `SMILES [whitespace] OPTIONAL_NAME`.

```bash
--lib /data/chembl_subset.smi
```

#### `--target <SMILES>`
The query molecule. Similarity is computed against this. Single SMILES string in quotes.

```bash
--target 'CC(=O)Oc1ccccc1C(=O)O'           # aspirin
--target 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'    # caffeine
```

### Output

#### `-o, --out <path>`
Output CSV path. Default: `results.csv`.

```bash
-o /tmp/screen_run_42.csv
```

### Fingerprint settings

#### `--method morgan | aptt | rd`
Fingerprint algorithm. Default: `morgan`.

| Value | Algorithm | Use case |
|---|---|---|
| `morgan` | ECFP-like circular | Most common; balanced |
| `aptt` | Atom-pair / topological torsion | Captures longer-range features |
| `rd` | RDKit path-based | Substructure-focused |

```bash
--method morgan
--method aptt
```

#### `--similarity Tanimoto | Dice | Sokal | McConnaughey | Kulczynski | Cosine`
Bit-vector similarity metric. Default: `Tanimoto`.

```bash
--similarity Tanimoto    # the field standard
--similarity Dice         # gives slightly more weight to common bits
```

#### `--fp-bits <N>`
Fingerprint bit length. Default: `2048`. Higher = more discriminative, slower similarity. Halving → ~2× faster.

```bash
--fp-bits 1024     # fast prefilter
--fp-bits 4096     # high precision for close-analog search
```

#### `--fp-radius <N>`
Morgan radius (ignored for `aptt`/`rd`). Default: `2` (≈ ECFP4).

```bash
--fp-radius 2     # ECFP4 — standard
--fp-radius 3     # ECFP6 — more specific scaffold
```

### Drug-likeness rules

#### `--rules <comma-list>`
Comma-separated list of rule sets. **ALL listed rules must pass.** Default: `lipinski`.

| Rule | Formula |
|---|---|
| `lipinski` | MW ≤ 500, HBD ≤ 5, HBA ≤ 10, logP ≤ 5 |
| `veber` | rotb ≤ 10, TPSA ≤ 140 |
| `ghose` | 160 ≤ MW ≤ 480, –0.4 ≤ logP ≤ 5.6, 40 ≤ MR ≤ 130 |
| `egan` | logP ≤ 5.88, TPSA ≤ 131.6 |
| `muegge` | 200 ≤ MW ≤ 600, –2 ≤ logP ≤ 5, … (full Muegge filter) |
| `ro3` | Rule of Three (fragments): MW ≤ 300, logP ≤ 3, … |
| `leadlike` | MW ≤ 350, logP ≤ 3.5 |
| `reos` | REOS (Walters 2002) |
| `gsk` | GSK 4/400: logP < 4, MW < 400 |
| `pfizer` | Pfizer 3/75: TPSA > 75, logP < 3 |

```bash
--rules lipinski                    # default
--rules lipinski,veber              # both must pass
--rules lipinski,veber,gsk          # toxicity-aware Lipinski-plus
--rules ro3                         # fragment-only screen
```

### Continuous-score gates

#### `--qed-min <x>`
Minimum QED score (Bickerton 2012). Range 0–1. Typical threshold: 0.5.

```bash
--qed-min 0.5      # standard drug-likeness cutoff
--qed-min 0.7      # strict
```

#### `--fsp3-min <x>`
Minimum Fsp3 (fraction of sp3 carbons). Lovering 2009 — correlates with clinical success.

```bash
--fsp3-min 0.42    # the Lovering threshold
```

#### `--stereo-max <n>`
Maximum number of stereocenters. Too many → harder synthesis.

```bash
--stereo-max 4
```

### Structural alerts

#### `--alerts <comma-list>`
Comma-separated alert catalogs. **Each listed catalog must be CLEAN** (no substructure match).

| Value | Catalog |
|---|---|
| `pains` | PAINS — pan-assay interference |
| `brenk` | BRENK — toxic/reactive groups |
| `lilly` | Eli Lilly MedChem Rules (CHEMBL_LINT) |
| `inpharmatica` | Inpharmatica alerts |
| `ames` | Ames mutagenicity (Benigni–Bossa subset) |
| `aggregators` | Shoichet aggregator-prone scaffolds |
| `extra` | NIH + ZINC + CHEMBL bundled |

```bash
--alerts pains,brenk                            # the medkim baseline
--alerts pains,brenk,ames,aggregators           # HTS-clean
--alerts pains,brenk,lilly,inpharmatica,ames    # paranoid
```

### Permeability gates

#### `--permeability bbb | gia | bbb_and_gia`
BOILED-Egg region the molecule must fall into.

| Value | Meaning |
|---|---|
| `gia` | Inside gastrointestinal absorption ellipse |
| `bbb` | Inside blood-brain barrier ellipse |
| `bbb_and_gia` | Inside both |

```bash
--permeability gia            # oral drug
--permeability bbb            # CNS drug
--permeability bbb_and_gia    # both
```

#### `--bbb-threshold <x>`
Minimum BBB-MPO score (Wager 2010-style). Typical range 0–6.

```bash
--bbb-threshold 3.0           # moderate CNS preference
--bbb-threshold 4.5           # strong CNS preference
```

### Two-stage workflow

#### `--prefilter <K>`
Activates two-stage screening. Stage 1 = fast similarity-only pass on the whole library (no standardization, no filters). Top-K kept. Stage 2 = full pipeline on those K only.

**Rule of thumb:** pick K = 5–10× the number of hits you actually want.

```bash
--prefilter 1000     # keep top 1000 by similarity → full filtering
--prefilter 10000    # for million-scale libraries
```

### Performance / behavior

#### `--threads <n>`
Number of OpenMP worker threads. `0` = all cores. Default: `0`.

```bash
--threads 0       # use everything
--threads 8       # leave headroom
```

#### `--no-standardize`
Skip the entire standardization pipeline (largest-fragment, normalize, uncharge, tautomer). Faster, but processes raw SMILES — salts/charges stay as input.

```bash
--no-standardize     # ~3× faster, less accurate
```

#### `--no-tautomer`
Skip **only** the tautomer canonicalization step (the slowest part of standardization). Keeps largest-fragment, normalize, and uncharge.

```bash
--no-tautomer        # ~2× faster, near-identical accuracy
```

---

## Screening presets

Common scenarios, ready to paste.

### Cheap & fast — exploratory scan

For a quick first look: minimal filters, just rank by similarity. Good for "is there anything close to my target at all?"

```bash
./build/sealab-vs \
    --lib library.smi \
    --target '<your SMILES>' \
    --rules lipinski \
    --no-tautomer \
    --threads 0 \
    -o results_cheap.csv
```

Expectation: ~3 000 mol/s on 16 cores; minimal filtering means many hits.

### Medium — drug-discovery baseline

Standard "med-chem clean" screen — drug-like, no PAINS/BRENK noise, oral-absorbable.

```bash
./build/sealab-vs \
    --lib library.smi \
    --target '<your SMILES>' \
    --rules lipinski,veber \
    --alerts pains,brenk \
    --qed-min 0.5 \
    --permeability gia \
    --threads 0 \
    -o results_medium.csv
```

Expectation: a sane shortlist suitable for direct docking input.

### Detailed — strict shortlist for high-cost docking

When docking is expensive and you want only the highest-confidence candidates.

```bash
./build/sealab-vs \
    --lib library.smi \
    --target '<your SMILES>' \
    --rules lipinski,veber,gsk,pfizer \
    --alerts pains,brenk,lilly,inpharmatica,ames,aggregators \
    --qed-min 0.6 \
    --fsp3-min 0.42 \
    --stereo-max 4 \
    --permeability gia \
    --bbb-threshold 2.5 \
    --threads 0 \
    -o results_detailed.csv
```

Expectation: very strict; many libraries collapse to dozens of hits.

### Massive library — two-stage screening

For ZINC, Enamine REAL, or any library > 1 M molecules.

```bash
./build/sealab-vs \
    --lib enamine_real.smi \
    --target '<your SMILES>' \
    --prefilter 10000 \
    --rules lipinski,veber \
    --alerts pains,brenk \
    --qed-min 0.5 \
    --threads 0 \
    -o results_massive.csv
```

Expectation: stage 1 ~60 000 mol/s, stage 2 only on the top 10 000. A 10 M-molecule library finishes in minutes.

### CNS-focused — brain-penetrant candidates

For targets in the central nervous system.

```bash
./build/sealab-vs \
    --lib library.smi \
    --target '<your SMILES>' \
    --rules lipinski,veber \
    --alerts pains,brenk \
    --qed-min 0.5 \
    --permeability bbb \
    --bbb-threshold 4.0 \
    --threads 0 \
    -o results_cns.csv
```

### Fragment-based — Rule-of-Three

For fragment libraries where Lipinski is too permissive.

```bash
./build/sealab-vs \
    --lib fragments.smi \
    --target '<your SMILES>' \
    --rules ro3 \
    --alerts pains,brenk \
    --threads 0 \
    -o results_fragments.csv
```

### Highest precision — close-analog search

Tight fingerprint settings, narrow filters.

```bash
./build/sealab-vs \
    --lib library.smi \
    --target '<your SMILES>' \
    --fp-bits 4096 \
    --fp-radius 3 \
    --similarity Tanimoto \
    --rules lipinski \
    --alerts pains,brenk \
    --threads 0 \
    -o results_analog.csv
```

---

## Output format

The CSV is sorted by descending similarity. One row per molecule that passes every selected gate. Columns:

| Column | Meaning |
|---|---|
| `rank` | 1-based rank by similarity |
| `smiles` | Canonical SMILES after standardization |
| `similarity` | Similarity to the target (0–1) |
| `mw, logp, tpsa, hbd, hba, rotb, arom, fsp3, heavy, stereo` | Physicochemical descriptors |
| `qed` | QED score (Bickerton 2012; 0–1) |
| `bbb_mpo` | BBB-MPO score (Wager-style) |
| `gia, bbb` | 1 if inside the corresponding BOILED-Egg ellipse |
| `pains_free, brenk_free, lilly_free, inpharmatica_free, ames_free, aggregators_free, extra_free` | 1 = no alert match for that catalog |
| `lipinski, veber, ghose, egan, muegge, ro3, leadlike, reos, gsk, pfizer` | 1 = passed that rule |

Every column is populated for every surviving row, so post-hoc filtering in pandas / Excel / awk is straightforward.

---

## Performance

Reference: 16-core Intel/AMD CPU.

| Mode | Throughput |
|---|---|
| Full pipeline (single stage) | ~1 500 – 2 500 mol/s |
| Full pipeline + `--no-tautomer` | ~3 000 – 5 000 mol/s |
| Stage 1 of `--prefilter` | ~50 000 – 80 000 mol/s |

For a 10 M-molecule ZINC subset on 16 threads:
- Single stage, full filtering: ~1 hour
- Two-stage with `--prefilter 10000`: ~3 minutes

Scales sub-linearly with cores (typically 6.5× on 16 threads) due to heap allocation contention.

---

## Documentation

- [`SETUP.txt`](SETUP.txt) — bilingual (EN/TR) full installation guide for a fresh machine
- [`ARCHITECTURE.txt`](ARCHITECTURE.txt) — bilingual (EN/TR) detailed design, algorithms, formulas, references — intended as the methods reference for a manuscript

---

## License

MIT — see [`LICENSE`](LICENSE).

Created by Hacı Aslan Onur İşcil, Ecem Bulut Turhan, Saliha Ece Acuner.
