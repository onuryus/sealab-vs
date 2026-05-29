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
- [How `--prefilter` works](#how---prefilter-works)
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

The original prototype was a single-threaded Python script and could not keep up with modern library sizes. The C++ rewrite uses RDKit C++, OpenMP, and a custom streaming reader. On a 16-core CPU it sustains **~4 600 mol/s on the full pipeline** or **~130 000 mol/s on the prefilter pass** (measured — see [Performance](#performance)), about an order of magnitude or more above the Python version, with deterministic output.

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

git clone https://github.com/onuryus/sealab-vs.git 
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

## How `--prefilter` works

Two-stage screening (the `--prefilter K` flag) is the single biggest
performance lever in sealab-vs. The trade-off is real, though, so here is
how it works in detail.

### The cost problem

Without `--prefilter`, every molecule in the library walks through the
full pipeline:

```
parse → standardize → descriptors → rules → PAINS/BRENK → QED →
optional alerts → BOILED-Egg → fingerprint → similarity → CSV
```

Most of that work is wasted on molecules that turn out to be irrelevant.
If you are searching a 1 M-molecule library for analogs of aspirin, you
do not really need to run PAINS substructure matching on methane.

### The two-stage idea

`--prefilter K` splits the work in two:

**Stage 1 — fast similarity-only sweep over the entire library:**

```
parse → fingerprint → similarity → keep (similarity, raw SMILES)
```

No standardization, no descriptors, no rules, no alerts, no BOILED-Egg.
Three cheap steps per molecule, ~130 000 mol/s on 16 cores.

After Stage 1, `std::nth_element` picks the **top K** records by
similarity in O(N) time — much cheaper than a full sort.

**Stage 2 — full pipeline on those K molecules only:**

The K molecules retained from Stage 1 are pushed through the complete
pipeline (standardize, descriptors, rules, alerts, BOILED-Egg, …). All
the expensive work happens here, but only on K records instead of N.

### Why it is fast

| Library size | Single stage | `--prefilter 10000` |
|---|---:|---:|
| 100 K mol | 22 s | < 1 s |
| 1 M mol | 3.6 min | ~8 s |
| 10 M mol | 36 min | ~80 s |

For a 1 M library with K = 1 000, the speedup is around **25×** —
almost the entire pipeline cost is replaced by a cheap similarity sweep.

### The trade-off you must understand

Stage 1 applies **no filters** — it ranks purely by similarity. This has
two consequences:

1. If the top-K molecules all happen to fail your filters
   (e.g. they all match PAINS), Stage 2 produces zero hits even though
   acceptable molecules might have existed further down the similarity
   ranking.
2. Molecules that look very similar to the target but fail an alert
   filter still occupy slots in the top-K, displacing slightly less
   similar molecules that would have passed.

**Rule of thumb:** pick K to be **5–10× the number of hits you actually
expect**. If you want ~100 candidates for docking, ask for `--prefilter
1000` to leave slack.

### When to use it

| Situation | Decision |
|---|---|
| Library larger than ~100 K | **Use it** |
| You want the highest-similarity candidates for docking | **Use it** |
| Library smaller than ~10 K | Skip — single stage is already fast |
| You want every molecule that passes your filters, regardless of rank | Skip — filters must run first |
| Scaffold diversity matters more than similarity rank | Skip — top-K clusters near the target |

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

Measured on a 16-core CPU with a 100 000-molecule test library, default settings (Morgan/2048/radius 2, Tanimoto, Lipinski + Veber rules, PAINS + BRENK alerts, full standardization):

| Configuration | Time (100 k mol) | Throughput |
|---|---:|---:|
| Single thread (baseline) | 196.6 s | 508 mol/s |
| 16 threads, full pipeline | **21.8 s** | **4 593 mol/s** |
| 16 threads, `--no-tautomer` | 19.5 s | 5 126 mol/s |
| 16 threads, `--no-standardize` | 17.6 s | 5 692 mol/s |
| 16 threads, `--prefilter 1000` | **0.76 s** | **132 100 mol/s** |

**Multi-core scaling:** 9.0× on 16 threads (508 → 4 593 mol/s). Sub-linear because of heap-allocator contention — RDKit allocates 5–6 ROMol objects per molecule.

**Extrapolations** (same hardware, same defaults):

| Library size | Single stage | Two-stage `--prefilter 10000` |
|---|---:|---:|
| 100 K mol | 22 s | < 1 s |
| 1 M mol | 3.6 min | ~8 s |
| 10 M mol | 36 min | ~80 s |

Real diverse libraries (ZINC, Enamine REAL, ChEMBL subsets) behave very close to these synthetic numbers because parsing, standardization, and fingerprint cost dominate — not cache locality of unique SMILES.

---

## Documentation

- [`SETUP.txt`](SETUP.txt) — bilingual (EN/TR) full installation guide for a fresh machine
- [`ARCHITECTURE.txt`](ARCHITECTURE.txt) — bilingual (EN/TR) detailed design, algorithms, formulas, references — intended as the methods reference for a manuscript

---

## Acknowledgments & citations

sealab-vs implements algorithms and rule sets from a number of openly
published works. If you use the program in a publication, please cite the
relevant primary references:

- **RDKit** — Landrum, G. *RDKit: Open-source cheminformatics*. https://www.rdkit.org (BSD-3 license)
- **QED** — Bickerton, G. R. et al. (2012) *Nat. Chem.* **4**, 90–98
- **Lipinski Rule of Five** — Lipinski, C. A. et al. (1997) *Adv. Drug Deliv. Rev.* **23**, 3–25
- **Veber rule** — Veber, D. F. et al. (2002) *J. Med. Chem.* **45**, 2615–2623
- **Ghose filter** — Ghose, A. K. et al. (1999) *J. Comb. Chem.* **1**, 55–68
- **Egan filter** — Egan, W. J. et al. (2000) *J. Med. Chem.* **43**, 3867–3877
- **Muegge filter** — Muegge, I. et al. (2001) *J. Med. Chem.* **44**, 1841–1846
- **Rule of Three** — Congreve, M. et al. (2003) *Drug Discov. Today* **8**, 876–877
- **Lead-likeness** — Teague, S. J. et al. (1999) *Angew. Chem. Int. Ed.* **38**, 3743–3748
- **REOS** — Walters, W. P. & Murcko, M. A. (2002) *Drug Discov. Today* **7**, S40–S47
- **GSK 4/400 rule** — Gleeson, M. P. (2008) *J. Med. Chem.* **51**, 817–834
- **Pfizer 3/75 rule** — Hughes, J. D. et al. (2008) *Bioorg. Med. Chem. Lett.* **18**, 4872–4875
- **Fsp3** — Lovering, F. et al. (2009) *J. Med. Chem.* **52**, 6752–6756
- **PAINS** — Baell, J. B. & Holloway, G. A. (2010) *J. Med. Chem.* **53**, 2719–2740
- **BRENK alerts** — Brenk, R. et al. (2008) *ChemMedChem* **3**, 435–444
- **Eli Lilly MedChem Rules** — Bruns, R. F. & Watson, I. A. (2012) *J. Med. Chem.* **55**, 9763–9772
- **Inpharmatica alerts** — Sutherland, J. J. et al. (2003) *J. Chem. Inf. Comput. Sci.* **43**, 1666–1673
- **Ames mutagenicity SMARTS** — Benigni, R. & Bossa, C. (2008) *Mutat. Res.* **659**, 248–261
- **BOILED-Egg** — Daina, A. & Zoete, V. (2016) *ChemMedChem* **11**, 1117–1121
- **BBB-MPO** — Wager, T. T. et al. (2010) *ACS Chem. Neurosci.* **1**, 435–449

The SMARTS rule sets bundled by RDKit (PAINS, BRENK, Lilly MedChem Rules,
Inpharmatica, NIH, ZINC, CHEMBL) are accessed via the RDKit FilterCatalog
API; they originate from the respective open-access publications cited
above and are not redistributed by this project independently of RDKit.

The Ames and Aggregators SMARTS lists shipped in `src/chem.cpp` are
small, hand-curated subsets of the patterns published in the
corresponding peer-reviewed papers; users are welcome to extend them.

## License

MIT — see [`LICENSE`](LICENSE).

Created by Hacı Aslan Onur İşcil, Ecem Bulut Turhan, Saliha Ece Acuner.
