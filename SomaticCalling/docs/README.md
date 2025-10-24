# SomaticCalling: Parabricks + GATK Somatic (Tumor/Normal) – Nextflow

> **End‑to‑end GPU‑accelerated somatic short‑read pipeline**: `fq2bam → Picard RG → BQSR (table+apply) → Mutect → FilterMutectCalls`

---

## ✨ What this pipeline does

- **Aligns** Tumor/Normal FASTQs with \*\*Parabricks \*\*`` (GPU)
- **Adds read groups** with **Picard** `AddOrReplaceReadGroups`
- **Recalibrates base quality** with **Parabricks** `bqsr` + `applybqsr` (GPU)
- **Calls somatic variants** with **Parabricks MutectCaller** (GPU)
- **Filters calls** with \*\*GATK \*\*`` (requires `*.stats` + VCF index)
- Handles **multiple samples** and both **custom** and **GIAB‑style** filenames
- Supports **multi‑lane** samples (optionally via `--in-fq-list`)

> The pipeline is contained under `SomaticCalling/` in this repository; run it from the **repo root**.

---

## 🗺️ Contents

- [Requirements](#-requirements)
- [Folder layout](#-folder-layout)
- [Input naming](#-input-naming)
- [Reference data](#-reference-data)
- [Quick start](#-quick-start)
- [Parameters](#-parameters)
- [Outputs](#-outputs)
- [Reproducibility](#-reproducibility)
- [Example dataset (GIAB HG008)](#-example-dataset-giab-hg008)
- [Troubleshooting](#-troubleshooting)
- [Roadmap](#-roadmap)
- [Cite & License](#-cite--license)

---

## ✅ Requirements

- **Nextflow** ≥ 25.04.1
- **Docker** with **NVIDIA GPU** support (nvidia‑container‑toolkit)
- Access to the container images:
  - `nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1`
  - `broadinstitute/picard:2.27.5`
  - `broadinstitute/gatk:4.4.0.0`

> 💡 The tools run **inside containers**; your local Conda env is optional and only needed for Nextflow.

---

## 📁 Folder layout

### Default repo layout (recommended)

```
TACoGenomics/
├─ SomaticCalling/
│  ├─ main.nf                 # the pipeline
│  ├─ conf/
│  │  ├─ base.config          # parameters & defaults (expects paths under SomaticCalling/)
│  │  └─ docker_gpu.config    # container pins + --gpus all
│  ├─ docs/
│  │  ├─ README.md            # pipeline‑specific usage
│  ├─ scripts/                # optional helpers
│  ├─ InputData/              # FASTQs (not tracked)
│  ├─ GenomeRef/              # reference + indexes (not tracked)
│  └─ results/                # published outputs
└─ (repo root README & LICENSE)
```

---

## 📦 Default layout & path assumptions

**When you clone this repository and run from the repo root**, the pipeline assumes the following **by default** (from `SomaticCalling/conf/base.config`:

| Item          | Default path (relative to repo root)                    | Must exist?  | How to override                        |
| ------------- | ------------------------------------------------------- | ------------ | -------------------------------------- |
| Input folder  | `SomaticCalling/InputData`                              | **Yes**      | `--input_dir /path/to/fastqs`          |
| Output folder | `SomaticCalling/results`                                | No (created) | `--outdir /path/to/out`                |
| Reference FA  | `SomaticCalling/GenomeRef/GRCh38.d1.vd1.fa`             | **Yes**      | `--ref /path/to/ref.fa`                |
| Known sites   | `SomaticCalling/GenomeRef/*{dbsnp,mills,known_indels}*` | **Yes**      | `--dbsnp … --mills … --known_indels …` |

**Key points**

- The **input folder must exist** and contain your `*.fastq.gz` files.
- The **output folder will be created** automatically if it does not exist.
- All **reference sidecar indexes** (`.fai`, `.dict`, BWA indexes) must sit next to the FASTA (or supply absolute paths).
- If you prefer a different layout, **override paths at run time** (examples below).

---

## 🧾 Input naming

The pipeline auto‑pairs Tumor/Normal and R1/R2 using either style:

**Custom style**

```
Sample1_T_L001_001.R1.fastq.gz  Sample1_T_L001_001.R2.fastq.gz
Sample1_N_L001_001.R1.fastq.gz  Sample1_N_L001_001.R2.fastq.gz
```

**GIAB style**

```
HG001-T_..._L001_001.R{1,2}.fastq.gz   (Tumor)
HG001-N-..._L001_001.R{1,2}.fastq.gz   (Normal)
```

- The **sample token** (`Sample1`, `HG001`, etc.) must match between T & N.
- Multi‑lane support: lanes (e.g., `L001`, `L002`) can be merged during `fq2bam` using a list file (see Advanced note below).

> 🔎 The pairing regex is implemented in `main.nf` and supports both styles. If you adopt different names, update the regex accordingly.

---

## 🧬 Reference data

Place these in the reference folder (layout A: `SomaticCalling/GenomeRef/`, layout B: `GenomeRef/`). The filenames below mirror the tested setup you used:

```
GRCh38.d1.vd1.fa
GRCh38.d1.vd1.fa.fai
GRCh38.d1.vd1.dict
GRCh38.d1.vd1.fa.amb
GRCh38.d1.vd1.fa.ann
GRCh38.d1.vd1.fa.bwt
GRCh38.d1.vd1.fa.pac
GRCh38.d1.vd1.fa.sa
Homo_sapiens_assembly38.dbsnp138.vcf
Homo_sapiens_assembly38.dbsnp138.vcf.idx
Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
Homo_sapiens_assembly38.known_indels.vcf.gz
Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
```

> You may also see downloaded archives like `GRCh38.d1.vd1.fa.tar.gz` or `GRCh38.d1.vd1_GATK_indices.tar.gz`; those are **not used at run time** once extracted.

**Checklist before running**

- FASTA present **and** indexed (`.fai`) with a sequence dictionary (`.dict`).
- BWA index files present (`.amb/.ann/.bwt/.pac/.sa`).
- Known-sites VCFs present, with indexes: `.idx` for the plain VCF, `.tbi` for bgzipped VCFs.

**Build indexes (if needed):**

```bash
# samtools faidx
samtools faidx <ref_dir>/GRCh38.d1.vd1.fa

# GATK sequence dictionary
gatk CreateSequenceDictionary \
  -R <ref_dir>/GRCh38.d1.vd1.fa \
  -O <ref_dir>/GRCh38.d1.vd1.dict

# BWA index
bwa index <ref_dir>/GRCh38.d1.vd1.fa
```

\-------------------------------------------------- | --------------------------------- | | `GRCh38.d1.vd1.fa`                                 | Reference FASTA                   | | `GRCh38.d1.vd1.fa.fai`                             | FASTA index (samtools)            | | `GRCh38.d1.vd1.dict`                               | Sequence dictionary (GATK/Picard) | | `GRCh38.d1.vd1.fa.{amb,ann,bwt,pac,sa}`            | BWA indexes                       | | `Homo_sapiens_assembly38.dbsnp138.vcf[.gz]`        | Known sites                       | | `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | Known sites                       | | `Homo_sapiens_assembly38.known_indels.vcf.gz`      | Known sites                       |

**Build indexes (if needed):**

```bash
# samtools faidx GRCh38.d1.vd1.fa
samtools faidx SomaticCalling/GenomeRef/GRCh38.d1.vd1.fa

# GATK create sequence dictionary
gatk CreateSequenceDictionary \
  -R SomaticCalling/GenomeRef/GRCh38.d1.vd1.fa \
  -O SomaticCalling/GenomeRef/GRCh38.d1.vd1.dict

# BWA index
bwa index SomaticCalling/GenomeRef/GRCh38.d1.vd1.fa
```

> Keep filenames consistent; the pipeline points to the `.fa` and expects its sidecar indexes alongside it.

---

## 🚀 Quick start

### Option A — Run with defaults (layout A: inside `SomaticCalling/`)

> Run this from the **repo root** so relative paths resolve correctly.

```bash
nextflow run SomaticCalling \
  -c SomaticCalling/conf/base.config \
  -c SomaticCalling/conf/docker_gpu.config \
  -resume -ansi-log false
```

This expects `SomaticCalling/InputData/` and `SomaticCalling/GenomeRef/` to exist.

---

### Option B — Use explicit paths (works for both layouts)

```bash
nextflow run SomaticCalling \
  -c SomaticCalling/conf/docker_gpu.config \
  --input_dir   /ABS/PATH/TO/InputData \
  --outdir      /ABS/PATH/TO/results \
  --ref         /ABS/PATH/TO/GenomeRef/GRCh38.d1.vd1.fa \
  --dbsnp       /ABS/PATH/TO/GenomeRef/Homo_sapiens_assembly38.dbsnp138.vcf \
  --mills       /ABS/PATH/TO/GenomeRef/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known_indels /ABS/PATH/TO/GenomeRef/Homo_sapiens_assembly38.known_indels.vcf.gz \
  -resume -ansi-log false
```

Rules:

- `--input_dir` must exist and contain FASTQs; `--outdir` will be created.
- All reference sidecars must be present (see checklist above).

---

### Option C — Run from inside `SomaticCalling/` (layout A) with absolute paths

If you prefer `cd SomaticCalling`, keep **absolute** `--input_dir/--outdir/--ref` paths to avoid resolving to `SomaticCalling/SomaticCalling/...`.

---

### Option D — Flat layout helper (layout B) using a small config

Create `SomaticCalling/conf/flat_root.config` with:

```groovy
params {
  input_dir    = 'InputData'
  outdir       = 'results'
  ref          = 'GenomeRef/GRCh38.d1.vd1.fa'
  dbsnp        = 'GenomeRef/Homo_sapiens_assembly38.dbsnp138.vcf'
  mills        = 'GenomeRef/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
  known_indels = 'GenomeRef/Homo_sapiens_assembly38.known_indels.vcf.gz'
}
```

Then run from the repo root:

```bash
nextflow run SomaticCalling \
  -c SomaticCalling/conf/flat_root.config \
  -c SomaticCalling/conf/docker_gpu.config \
  -resume -ansi-log false
```

This matches a repo that has `InputData/`, `GenomeRef/`, and `results/` at the **repo root**, exactly like your quick tests.

---

### Option B — Use your own folders (explicit paths)

```bash
nextflow run SomaticCalling \
  -c SomaticCalling/conf/docker_gpu.config \
  --input_dir   /data/fastqs \
  --outdir      /results/somatic_run1 \
  --ref         /refs/GRCh38.d1.vd1.fa \
  --dbsnp       /refs/Homo_sapiens_assembly38.dbsnp138.vcf \
  --mills       /refs/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known_indels /refs/Homo_sapiens_assembly38.known_indels.vcf.gz \
  -resume -ansi-log false
```

**Rules**

- `--input_dir` must exist and contain FASTQs.
- `--outdir` will be created if missing.
- All reference sidecars must be present (or provide their absolute paths if you relocate them).

---

### Option C — Run from inside `SomaticCalling/`

If you want to `cd SomaticCalling` and run there, either:

1. keep using **absolute paths**, or
2. add a small alt‑config where defaults are `InputData`, `results`, `GenomeRef` relative to the subfolder.

> Default configs here assume you run from the **repo root**. Running inside the subfolder without overrides will point to the wrong relative paths.

---

## 🔧 Parameters

Defined in `SomaticCalling/conf/base.config` (overridable via `--param=value`):

| Param          | Default                                                                     | Description                     |
| -------------- | --------------------------------------------------------------------------- | ------------------------------- |
| `input_dir`    | `SomaticCalling/InputData`                                                  | Folder with FASTQs              |
| `outdir`       | `SomaticCalling/results`                                                    | Publish directory               |
| `ref`          | `SomaticCalling/GenomeRef/GRCh38.d1.vd1.fa`                                 | Reference FASTA                 |
| `dbsnp`        | `SomaticCalling/GenomeRef/Homo_sapiens_assembly38.dbsnp138.vcf`             | Known sites                     |
| `mills`        | `SomaticCalling/GenomeRef/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | Known sites                     |
| `known_indels` | `SomaticCalling/GenomeRef/Homo_sapiens_assembly38.known_indels.vcf.gz`      | Known sites                     |
| `bwa_nstreams` | `1`                                                                         | BWA threads (per‑GPU streams)   |
| `low_memory`   | `true`                                                                      | Pass `--low-memory` to `fq2bam` |

**Containers** are pinned in `SomaticCalling/conf/docker_gpu.config`.

---

## 📤 Outputs

```
SomaticCalling/results/
├─ 01_fq2bam/               SampleX_Tumor.raw.bam, SampleX_Normal.raw.bam
├─ 02_readgroups/           SampleX_Tumor.rg.bam, SampleX_Normal.rg.bam
├─ 03_bqsr_tables/          SampleX_Tumor.recal.txt, SampleX_Normal.recal.txt
├─ 04_bqsr_bams/            SampleX_Tumor.bqsr.bam, SampleX_Normal.bqsr.bam
├─ 05_mutect/               SampleX.mutect.vcf.gz(.tbi) + .stats
└─ 06_filtered_vcf/         SampleX.filtered.vcf.gz(.tbi)
```

---

## 🔁 Reproducibility

- **Containers** pin exact versions (can be further pinned by **digest**).
- **Nextflow** version requirement: `>= 25.04.1`.
- `-resume` ensures cached results are reused if code/inputs unchanged.

**GPU notes**

- Docker needs `--gpus all` (set in `docker_gpu.config`).
- Validate with `nvidia-smi` inside containers if troubleshooting.
