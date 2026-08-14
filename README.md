# MErlin

**M**ethylation-driven **E**xpression & **R**egulation **Link**age in **I**nteracting **N**uclear-domains — a prokaryotic multi-omics toolkit, shipped as a **Claude skill** with a full command-line tool underneath.

Point it at bacterial DNA methylation data and it will cross that methylation with expression, motifs, methyltransferases, genome position and Hi-C — running only the modules your inputs actually support.

> **New here?** The fastest path is [Quick start](#3-quick-start): create the conda environment, install the package, install the skill, then just describe your data in plain language.

---

## Table of contents

- [1. What this repository contains](#1-what-this-repository-contains)
- [2. The two halves: the skill and the tool](#2-the-two-halves-the-skill-and-the-tool)
  - [2.1 The skill (the recommended way in)](#21-the-skill-the-recommended-way-in)
  - [2.2 The tool (the `merlin` command)](#22-the-tool-the-merlin-command)
- [3. Quick start](#3-quick-start)
- [4. Step 1 · Create the conda environment](#4-step-1--create-the-conda-environment)
  - [4.1 With the provided `environment.yml`](#41-with-the-provided-environmentyml)
  - [4.2 By hand](#42-by-hand)
  - [4.3 Without conda](#43-without-conda)
- [5. Step 2 · Install the MErlin package](#5-step-2--install-the-merlin-package)
- [6. Step 3 · Install the skill](#6-step-3--install-the-skill)
  - [6.1 Build the `.skill` file](#61-build-the-skill-file)
  - [6.2 Claude desktop app or claude.ai](#62-claude-desktop-app-or-claudeai)
  - [6.3 Claude Code](#63-claude-code)
  - [6.4 Check that it took](#64-check-that-it-took)
- [7. Using the skill](#7-using-the-skill)
- [8. Using it as a command-line tool](#8-using-it-as-a-command-line-tool)
  - [8.1 The modules](#81-the-modules)
  - [8.2 One module at a time](#82-one-module-at-a-time)
  - [8.3 A whole pipeline with `merlin run`](#83-a-whole-pipeline-with-merlin-run)
  - [8.4 On a cluster](#84-on-a-cluster)
- [9. Docker](#9-docker)
- [10. Upgrading and removing](#10-upgrading-and-removing)
- [11. Troubleshooting](#11-troubleshooting)
- [12. Documentation map](#12-documentation-map)
- [13. Citation and licence](#13-citation-and-licence)

---

## 1. What this repository contains

```text
MErlin/
├── README.md                  ← you are here
├── environment.yml            ← the conda environment
├── LICENSE
└── merlin-multiomics/         ← THE SKILL — zip this folder to get a .skill file
    ├── SKILL.md               ← what Claude reads to drive MErlin correctly
    ├── SETUP.md               ← the long-form install guide
    ├── references/            ← tools, statistics, preflight, interpretation, troubleshooting
    ├── scripts/               ← install_merlin.sh, check_env.sh, roundtable.R
    └── assets/
        └── merlin_multiomics-2.1.0-py3-none-any.whl   ← THE TOOL
```

Everything you need is in the repository — the skill folder and the installable wheel. Nothing is downloaded from PyPI.

> ⚠️ **Do not rename the `merlin-multiomics/` folder.** Its name has to match the `name:` field inside `SKILL.md`, and a mismatch is the most common cause of a rejected skill upload.

---

## 2. The two halves: the skill and the tool

MErlin is one codebase with two front doors. They are independent: the skill without the package lets Claude explain MErlin but not run it; the package without the skill works fine from a terminal. **Install both** — that is what this README walks you through.

### 2.1 The skill (the recommended way in)

A *skill* is a folder of instructions that teaches Claude how to use a tool properly. `merlin-multiomics` is what turns "I have three WT and three knockout bedMethyl files and a Prokka GFF3" into a correct, defensible analysis, because it encodes the things that are easy to get wrong:

| The skill makes Claude… | …instead of |
| --- | --- |
| run `merlin doctor`, then `merlin preflight`, **before** any analysis | running first and discovering the seqids never matched |
| ask which condition is the treatment, and what a positive log2FC means | guessing, and inverting every directional conclusion |
| pass replicates as a list so the beta-binomial LRT is licensed | concatenating them into one fake sample |
| report numbers straight from MErlin's output tables, naming the file | recomputing statistics in pandas and mixing two provenances |
| say "association", and reserve "regulates" for what a knockout licenses | overclaiming mechanism from a correlation |
| treat a methylation-only run as a complete run | telling you your project is incomplete because you have no Hi-C |

In short: the skill is the **methodological guardrail**. It is why this repository is skill-first rather than a bare CLI.

### 2.2 The tool (the `merlin` command)

Underneath sits a single `merlin` CLI with **seventeen subcommands in three groups**. Every module is independent and runs on whatever inputs exist — that modularity is the design, not a fallback.

```text
analysis        pileup  xam  dixam  spacem  mematch  diem  hicam  phase  discover  mtase
integration     harmonise  evidence  report
orchestration   preflight  doctor  run  export
```

You can drive it entirely from a terminal, with no Claude involved — see [§8](#8-using-it-as-a-command-line-tool).

---

## 3. Quick start

Four commands, from a fresh clone to a working install:

```bash
git clone https://github.com/IacopoPasseri/MErlin.git
cd MErlin
conda env create -f environment.yml     # step 1 — environment + MErlin package
conda activate merlin
merlin doctor                           # step 3 — confirm it works
```

Then build and upload the skill:

```bash
zip -r merlin-multiomics.skill merlin-multiomics -x '*.DS_Store' '*__pycache__*'
```

…and drop that file into **Customize → Skills** in the Claude app ([§6.2](#62-claude-desktop-app-or-claudeai)).

The three steps are explained in full below.

---

## 4. Step 1 · Create the conda environment

**Prerequisites:** Python **3.10 or newer**, and [conda](https://docs.conda.io/projects/miniconda/) or [mamba](https://mamba.readthedocs.io/). Nothing else is mandatory.

### 4.1 With the provided `environment.yml`

Run this **from the repository root**, because the file installs the bundled wheel by relative path:

```bash
conda env create -f environment.yml
conda activate merlin
```

That single command gives you Python 3.11, every MErlin dependency, the optional external binaries (`modkit`, `samtools`, `minimap2`) and MErlin itself. `mamba env create -f environment.yml` does the same, faster.

To use a different environment name:

```bash
conda env create -f environment.yml -n merlin-dev
```

### 4.2 By hand

If you would rather build it yourself, or you want to skip the optional binaries:

```bash
conda create -n merlin -c conda-forge -c bioconda python=3.11 pip
conda activate merlin

# optional external binaries — each one has a fallback inside MErlin
conda install -c bioconda modkit samtools minimap2
```

| Binary | Needed for | Without it |
| --- | --- | --- |
| `modkit` | the reference modBAM pileup | MErlin's own pysam pileup runs instead |
| `minimap2` + `samtools` | aligning an *unaligned* modBAM | `merlin pileup` refuses and prints the command to run by hand |
| R + `optparse`, `data.table`, `circlize`, `Biostrings` | the `roundtable.R` circos figure | `merlin spacem --circular` draws one without R |

R is only ever needed for that one figure. If you want it:

```bash
conda install -c conda-forge -c bioconda r-base r-optparse r-data.table r-circlize bioconductor-biostrings
```

### 4.3 Without conda

A plain virtual environment works just as well:

```bash
python3 -m venv ~/.venvs/merlin
source ~/.venvs/merlin/bin/activate
```

> **Important, either way:** for Claude to run MErlin for you, `merlin` has to be findable on `PATH`. Activate the environment **before** starting your Claude session, or tell Claude the absolute path to the executable (`which merlin`).

---

## 5. Step 2 · Install the MErlin package

Skip this if you used `environment.yml` in [§4.1](#41-with-the-provided-environmentyml) — it already installed the wheel.

**The one-liner** (installs the bundled wheel with all extras, then runs `merlin doctor`):

```bash
sh merlin-multiomics/scripts/install_merlin.sh
```

It is idempotent: if `merlin` already runs, it only reports and checks.

**By hand, from the bundled wheel:**

```bash
python3 -m pip install "merlin-multiomics/assets/merlin_multiomics-2.1.0-py3-none-any.whl[all]"
```

**By hand, from a source tree** (what you want if you are developing MErlin — changes take effect without reinstalling):

```bash
python3 -m pip install -e "/path/to/merlin-2.1.0[all]"
# or: sh merlin-multiomics/scripts/install_merlin.sh /path/to/merlin-2.1.0
```

**What `[all]` buys you.** It is the sensible default; individually:

| Extra | Packages | Without it |
| --- | --- | --- |
| *(core)* | numpy, pandas, scipy, matplotlib | nothing runs |
| `ml` | scikit-learn, statsmodels, seaborn | `diem` loses its classifier and logistic regression |
| `store` / `hic` | h5py | no `harmonise` store; no `.cool` reading |
| `cooler` | cooler | h5py's native `.cool` reader is used instead — fine for most files |
| `reads` | pysam | `phase` and `pileup` cannot run |
| `config` | pyyaml | `merlin run` accepts JSON configs only |

**Verify:**

```bash
merlin --version        # merlin 2.1.0
merlin doctor --check-r
```

`doctor` prints each dependency, then the external binaries, then a module-availability table — `ready`, `degraded` or `unavailable`, with the reason. A `warn` on `cooler`, `modkit` or `Rscript` is **normal**: each has a fallback inside MErlin.

---

## 6. Step 3 · Install the skill

### 6.1 Build the `.skill` file

A `.skill` file is simply a zip archive **of the folder**, not of its contents:

```bash
zip -r merlin-multiomics.skill merlin-multiomics -x '*.DS_Store' '*__pycache__*'
```

`SKILL.md` must sit one level down inside the archive (`merlin-multiomics/SKILL.md`). If the upload dialog only accepts `.zip`, rename the extension — the contents are identical.

### 6.2 Claude desktop app or claude.ai

1. Open **Customize → Skills**.
2. Click **+** → **Create skill** → **Upload a skill**.
3. Select `merlin-multiomics.skill`.
4. The skill appears in the list with a toggle. **Turn it on.**

Personal skills stay private to your account. On Team and Enterprise plans you can share one organisation-wide from the same screen.

### 6.3 Claude Code

Unzip it straight into your skills directory:

```bash
mkdir -p ~/.claude/skills
unzip merlin-multiomics.skill -d ~/.claude/skills/
ls ~/.claude/skills/merlin-multiomics/SKILL.md    # must exist
```

| Scope | Path | Available in |
| --- | --- | --- |
| Personal | `~/.claude/skills/merlin-multiomics/SKILL.md` | all your projects |
| Project | `.claude/skills/merlin-multiomics/SKILL.md` | that project only |

A skill enabled on claude.ai also syncs down to Claude Code, so you generally only need one of §6.2 and §6.3. A local copy takes precedence over the synced one if both exist.

### 6.4 Check that it took

Start a new session and ask something the skill should own:

> I have bedMethyl files for a *Sinorhizobium* WT and a *dam* knockout, three replicates each, plus a Prokka GFF3. What can MErlin tell me?

A working install answers with `merlin doctor` → `merlin preflight` → a replicated `dixam` design, and **asks which condition is the treatment** before running anything. If it instead offers to write a pandas script, the skill did not load — check that the folder name matches and the toggle is on.

---

## 7. Using the skill

Once the skill is on, you do not memorise flags. You describe your data:

> *"I have a basecalled modBAM with MM/ML tags, a Bakta GFF3 and the reference FASTA. No motif list. Find the motifs and tell me where methylation sits relative to oriC."*

> *"Three WT and three ∆dam bedMethyl files, plus a DESeq2 table. Positive log2FC means up in the knockout. Is methylation associated with expression?"*

Claude will then, in order:

1. **`merlin doctor`** — check the install and which modules are runnable.
2. **Inventory your inputs** and propose a plan (`merlin run --dry-run` prints exactly which modules run, and for each skipped one, the missing input). A plan with three skipped modules is a normal methylation-only run, not a broken configuration.
3. **`merlin preflight`** — catch the silent killers before they cost you a day. Both of MErlin's joins (on seqid and on locus_tag) fail *silently*: a seqid namespace mismatch returns zero overlaps, which is indistinguishable from a genome with no methylation.
4. **Ask you what the files cannot say** — contrast direction, which condition is which, replicate structure, which replicons matter, the knockout gene.
5. **Run the modules**, reporting progress, and stopping rather than feeding a failed output downstream.
6. **Integrate and report** — `harmonise` → `evidence` → a self-contained HTML report, delivered to you.

Two habits worth knowing about, because they change what your results mean:

- **`xam --include-zero`** before `diem`. Without it, features with no methylation call are *absent* rather than zero, and every correlation describes the methylated subset only.
- **Replicates as a list** (`--methylation a b c`), never concatenated. Three files licenses the beta-binomial LRT; `cat a b c` is one sample and licenses only a descriptive comparison.

---

## 8. Using it as a command-line tool

The skill is optional. Everything below runs in a plain terminal with the conda environment activated.

### 8.1 The modules

| Command | Needs | Produces |
| --- | --- | --- |
| `pileup` | basecalled modBAM (MM/ML tags) | bedMethyl — where a raw BAM becomes an input |
| `xam` | annotation + methylation | methylation per feature (per replicate + pooled) |
| `dixam` | annotation + 2 conditions | differential methylation, replicate-aware |
| `spacem` | annotation + methylation + FASTA | oriC/ter spatial distribution |
| `mematch` | + motif list | motif × methylation association, enrichment |
| `discover` | methylation + FASTA | motifs found *de novo* (no motif list needed) |
| `mtase` | annotation + motifs (+ knockout) | motif → MTase gene attribution |
| `diem` | expression table + `xam` table | expression × methylation, confounder-adjusted |
| `hicam` | `.cool`/`.mcool` (+ any other layer) | compartments, CIDs, contact-weighted co-methylation |
| `phase` | modBAM with MM/ML tags | hemimethylation, bistability, cis co-methylation |
| `harmonise` | a results directory | one HDF5 store, everything joined on feature and bin |
| `evidence` | that store | genes ranked by combined evidence across layers |
| `report` | that results directory | one self-contained HTML report |
| `preflight` `doctor` `run` `export` | — | validation, diagnostics, orchestration, workflow export |

Every module takes `--help`, which is the authority if the docs and the binary ever disagree. Every module also writes `<prefix>.provenance.txt` (command line, versions, SHA-256 of each input) and `<prefix>.warnings.txt` next to its tables.

### 8.2 One module at a time

**Always start with preflight.** Pass whatever exists; exit code `0` is clean, `1` warnings, `2` blocking.

```bash
merlin preflight --gff ANN.gff3 --fasta REF.fasta --methylation METH.bed \
    --motifs motifs.txt --geneexp deseq2.csv --hic contacts.cool
```

Then:

```bash
merlin pileup   --modbam calls.bam --fasta ref.fasta -o out/pileup
merlin xam      -g ann.gff3 --methylation r1.bed r2.bed r3.bed -o out --include-zero
merlin dixam    -g ann.gff3 --methylation wt1.bed wt2.bed wt3.bed \
                            --methylation2 ko1.bed ko2.bed ko3.bed \
                            --labels WT dam_ko -o out
merlin discover --fasta ref.fasta --methylation calls.bed -o out
```

> ⚠️ `--labels` is **positional** and defaults to `control treatment` regardless of file contents. A swapped pair inverts every differential call and looks completely normal.

### 8.3 A whole pipeline with `merlin run`

For more than one or two modules, write a config. `merlin run` resolves paths relative to the config, orders the modules correctly (`discover` before `mematch`, `xam` before `diem`), records the plan, and skips modules with a stated reason.

```yaml
# merlin.yaml
outdir: results
prefix: SmBL225C
annotation: annotation.gff3
reference:  reference.fasta
methylation:                                                   # 'modbam:' takes the
  WT:     [wt_rep1.bed.gz, wt_rep2.bed.gz, wt_rep3.bed.gz]      # same shape; with no
  dam_kd: [ko_rep1.bed.gz, ko_rep2.bed.gz, ko_rep3.bed.gz]      # 'methylation:' a
motifs:     motifs.txt        # omit it and 'discover' finds them   pileup runs first
expression: deseq2_results.csv
hic:        contacts.cool
modbam:     reads.modbam.bam
resolution: 5000
knockout_gene: dam
min_frac: 0.5
options:
  dixam: {features: CDS, regions: [upstream:-300..50, body]}
  diem:  {gff: annotation.gff3, fasta: reference.fasta, n-perm: 5000}
```

```bash
merlin run --config merlin.yaml --dry-run     # show the plan first
merlin run --config merlin.yaml
```

Then integrate — worth doing whenever more than one layer ran:

```bash
merlin harmonise --results results --gff ann.gff3 -o results/harmonised
merlin evidence  --store results/harmonised/merlin.h5 -o results/evidence
merlin report    --results results --store results/harmonised/merlin.h5 -o results
```

`evidence` combines per-layer statistics per gene (Stouffer + robust rank aggregation, BH across genes). It **ranks**; it does not test a new hypothesis. A gene high in the ranking is a gene worth looking at, not a gene shown to be regulated by methylation.

### 8.4 On a cluster

```bash
merlin export --config merlin.yaml --workflow nextflow    # or snakemake, or json
```

---

## 9. Docker

If you would rather not install anything, the source tree ships a Dockerfile with R, its packages, `samtools` and `minimap2` already in place:

```bash
docker build -t merlin:2.1.0 .
docker run --rm -v "$PWD":/data merlin:2.1.0 doctor
docker run --rm -v "$PWD":/data merlin:2.1.0 run --config merlin.yaml
```

The entrypoint is `merlin`, so pass subcommands directly. `/data` is the working directory inside the container — mount your data there and use relative paths.

To bake `modkit` into the image, pass the release tarball URL for your architecture:

```bash
docker build --build-arg MODKIT_URL=https://github.com/nanoporetech/modkit/releases/download/vX.Y.Z/modkit_vX.Y.Z_<target>.tar.gz -t merlin:2.1.0 .
```

---

## 10. Upgrading and removing

```bash
# package
python3 -m pip install --upgrade "merlin-multiomics/assets/merlin_multiomics-<new>-py3-none-any.whl[all]"
python3 -m pip uninstall merlin-multiomics

# environment
conda env remove -n merlin

# local skill
rm -rf ~/.claude/skills/merlin-multiomics
```

For a skill uploaded to claude.ai, upload the new version from **Customize → Skills** and it replaces the old one; the toggle there also removes it.

---

## 11. Troubleshooting

| Symptom | Cause and fix |
| --- | --- |
| Upload rejected: "folder name doesn't match" | the zip's top-level folder must be `merlin-multiomics` |
| Upload rejected: "missing SKILL.md" | zip the *folder*, not its contents — `SKILL.md` must sit one level down |
| Claude ignores the skill | the toggle is off, or a local `.claude/skills/merlin-multiomics` is shadowing the synced copy |
| `merlin: command not found` | not installed, wrong environment, or the user script directory is not on `PATH` |
| `ModuleNotFoundError: h5py` | installed without `[all]` — reinstall with the extras ([§5](#5-step-2--install-the-merlin-package)) |
| `pileup` says modkit is missing | expected; the internal backend runs. Install `modkit` only if you need modkit-identical numbers |
| `pip` refuses: externally-managed-environment | use the conda environment, a venv, or `--break-system-packages` |
| `merlin doctor` says a module is `unavailable` | the reason is printed next to it, and it names the package to install |
| A module returns zero overlaps | almost always a seqid namespace mismatch — run `merlin preflight`, do **not** interpret the zero |

Two quick diagnostics:

```bash
sh merlin-multiomics/scripts/check_env.sh    # can this machine run MErlin?
merlin doctor --check-r                      # what is available, and why not
```

If `merlin` installed but is not on `PATH`:

```bash
python3 -m merlin.cli doctor
python3 -m site --user-base    # then add .../bin to PATH
```

---

## 12. Documentation map

Everything below lives inside `merlin-multiomics/`. Claude reads these on demand; you can read them too.

| File | What it covers |
| --- | --- |
| `SKILL.md` | the workflow Claude follows, and the rules it works under |
| `SETUP.md` | the long-form installation guide (skill, package, Docker, R) |
| `references/tools.md` | every module's real arguments, outputs and columns |
| `references/statistics.md` | which model each experimental design licenses, and how the tests are calibrated |
| `references/run_config.md` | the `merlin run` YAML schema, with worked configs |
| `references/preflight.md` | every preflight finding and its fix |
| `references/interpretation.md` | what MErlin's statistics do and do not license |
| `references/troubleshooting.md` | analysis failures and their causes |

**The rules that matter most**, if you read nothing else:

- **Never present a zero result before preflight passes.** "No methylation found" and "the seqids didn't match" produce identical-looking output.
- **Match the claim to the design.** With one sample per condition, the p-value speaks to sampling of reads, not to biological reproducibility. Replicates change that; nothing else does.
- **Association, not causation.** MErlin quantifies co-variation, and significance does not convert a correlation into a mechanism. The narrow exception is `mtase` with a knockout, which licenses "motif X is lost when gene Y is deleted".
- **Never loosen a threshold to obtain a result**, and never edit an input file to make a join work. Fix the flag, not the data.

---

## 13. Citation and licence

TO ADD
