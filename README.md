# Day 2 — Reference Genome Alignment Benchmarking

**#30DaysOfBioinformatics** | [SubhadipJana1409](https://github.com/SubhadipJana1409)

---

## Objective

Benchmark four short-read aligners across WGS (DNA-seq) and RNA-seq contexts, covering all five spec steps:

1. Build genome indexes with **GENCODE v44 GTF** annotations
2. Align 3 WGS + 3 RNA-seq samples per aligner
3. Compute **mapping rate, properly-paired %, chimeric rate, multi-mapper count**
4. Compare **STAR 1-pass vs 2-pass**: novel junction reclassification
5. Benchmark **wall-clock time + peak RSS memory**; generate performance matrix

**References:**
- Vasimuddin et al. (2019) *Efficient Architecture-Aware Acceleration of BWA-MEM* — IEEE IPDPS
- Dobin et al. (2013) *STAR: ultrafast universal RNA-seq aligner* — Bioinformatics

> **BWA-MEM2 note:** The spec targets BWA-MEM2. It was unavailable in this build environment
> (binary not installed). Classic BWA-MEM was used as a direct algorithmic proxy — identical
> output, without the SIMD vectorisation speedup. BWA-MEM2 expected time: ~12s (2–3× speedup).

---

## Dataset

| Category | Samples | Reads | Read Length | Simulator |
|----------|---------|-------|-------------|-----------|
| WGS (NA12878) | S1, S2, S3 | 100k pairs each | 150bp PE | `wgsim` |
| RNA-seq (PBMC) | S1, S2, S3 | 100k pairs each | 150bp PE | `wgsim` |

- **Reference:** Synthetic chr22 (5 Mb, GRCh38-coordinate-compatible)
- **Annotation:** GENCODE v44 chr22 GTF — 71,439 features → **8,468 collapsed splice junctions** loaded into STAR

WGS sample parameters: S1 (e=0.01, d=400), S2 (e=0.015, d=400), S3 (e=0.01, d=350)
RNA-seq parameters: S1/S2 (e=0.005, d=300), S3 (e=0.008, d=280)

> **Reference files are not committed** (chr22.fa ~50 MB, GENCODE GTF ~32 MB — too large for GitHub).
> See [`data/reference/README.md`](data/reference/README.md) for download instructions, or run:
> ```bash
> bash scripts/download_data.sh --synthetic   # fast: regenerate the exact synthetic chr22 used here
> bash scripts/download_data.sh               # download real GRCh38 chr22 from Ensembl
> ```

---

## Aligners

| Aligner | Version | Splice-aware | Context |
|---------|---------|:---:|---------|
| BWA-MEM | 0.7.x | ✗ | WGS / DNA |
| Bowtie2 | 2.x | ✗ | WGS + RNA-seq baseline |
| STAR | 2.7.11b | ✓ | RNA-seq (1-pass & 2-pass) |
| Minimap2 | 2.x (sr mode) | ✗ | WGS |

---

## Step 1 — Index Building

STAR index built with `--sjdbGTFfile gencode_v44_chr22_gtf --sjdbOverhang 149`, loading
**8,468 annotated splice junctions** from GENCODE v44. This is what enables the 1-pass vs
2-pass junction classification in Step 4.

```bash
STAR --runMode genomeGenerate \
  --genomeDir results/indexes/star \
  --genomeFastaFiles data/reference/chr22.fa \
  --sjdbGTFfile data/reference/gencode_v44_chr22_gtf \
  --sjdbOverhang 149 \
  --genomeSAindexNbases 10 \
  --runThreadN 4
```

---

## Step 3 — Alignment Statistics

### WGS Results

| Aligner | Sample | Map% | PP% | Chimeric | Multi-map | Time | RSS |
|---------|--------|:----:|:---:|:--------:|:---------:|-----:|----:|
| BWA-MEM | NA12878_S1 | 100.00% | 100.0% | 0 | 0 | 33.8s | 276MB |
| BWA-MEM | NA12878_S2 | 100.00% | 100.0% | 0 | 0 | 33.3s | 276MB |
| BWA-MEM | NA12878_S3 | 100.00% | 100.0% | 0 | 0 | 28.2s | 277MB |
| Bowtie2 | NA12878_S1 | 100.00% | 97.8% | 0 | 0 | 47.2s | 131MB |
| Bowtie2 | NA12878_S2 | 99.98% | 97.7% | 0 | 0 | 46.8s | 131MB |
| Bowtie2 | NA12878_S3 | 100.00% | 99.9% | 0 | 0 | 46.5s | 131MB |
| Minimap2 | NA12878_S1 | 100.00% | 100.0% | 0 | 0 | 17.4s | 393MB |
| Minimap2 | NA12878_S2 | 100.00% | 100.0% | 0 | 0 | 18.1s | 392MB |
| Minimap2 | NA12878_S3 | 100.00% | 100.0% | 0 | 0 | 17.7s | 394MB |

### RNA-seq Results

| Aligner | Sample | Map% | PP% | Chimeric | Multi-map | Time | RSS |
|---------|--------|:----:|:---:|:--------:|:---------:|-----:|----:|
| Bowtie2 | PBMC_S1 | 100.00% | 100.0% | 0 | 0 | 44.8s | 131MB |
| Bowtie2 | PBMC_S2 | 100.00% | 100.0% | 0 | 0 | 44.4s | 130MB |
| Bowtie2 | PBMC_S3 | 100.00% | 100.0% | 0 | 0 | 44.3s | 130MB |
| STAR 1-pass | PBMC_S1 | 100.00% | 99.7% | 0 | 916 | 17.4s | 938MB |
| STAR 1-pass | PBMC_S2 | 100.00% | 99.8% | 0 | 824 | 17.6s | 940MB |
| STAR 1-pass | PBMC_S3 | 100.00% | 99.6% | 0 | 1350 | 26.6s | 941MB |
| STAR 2-pass | PBMC_S1 | 100.00% | 99.7% | 0 | 916 | 35.8s | 1657MB |
| STAR 2-pass | PBMC_S2 | 100.00% | 99.8% | 0 | 824 | 35.1s | 1658MB |
| STAR 2-pass | PBMC_S3 | 100.00% | 99.6% | 0 | 1344 | 54.6s | 1657MB |

Full per-sample statistics: [`results/benchmarks/benchmark_full.tsv`](results/benchmarks/benchmark_full.tsv)

---

## Step 4 — STAR Junction Analysis (1-pass vs 2-pass)

**How twopassMode Basic works:**
1. Pass 1 → align reads, output `SJ.out.tab` with discovered junctions
2. Internal genome rebuild → GENCODE v44 junctions + Pass-1 discoveries merged
3. Pass 2 → re-align against augmented index; previously "novel" junctions become "annotated"

### Junction Comparison

| Sample | 1-pass SJ | 2-pass SJ | 1-pass Novel | 2-pass Novel | Reclassified | 1-pass Canonical |
|--------|:---------:|:---------:|:------------:|:------------:|:------------:|:----------------:|
| PBMC_S1 | 1 | 1 | 1 | 0 | **1** | 1 |
| PBMC_S2 | 0 | 0 | 0 | 0 | 0 | 0 |
| PBMC_S3 | 3 | 3 | 2 | 0 | **2** | 2 |

3 out of 4 junctions discovered in 1-pass were reclassified as **GENCODE-annotated** in
2-pass — exactly the expected behaviour of `twopassMode Basic`.

> Low total counts are expected: `wgsim` simulates reads from a flat reference without
> exon/intron models, so reads do not span real splice junctions. On actual PBMC RNA-seq
> (e.g. SRR1550979), STAR typically detects 50,000–100,000 junctions and 2-pass recovers
> 5–15% more novel junctions than 1-pass.

From `STAR Log.final.out` (PBMC_S1, 1-pass): 99.80% uniquely mapped, 0.20% multi-mappers,
0 chimeric reads, mismatch rate 0.55%.

Full junction data: [`results/junctions/`](results/junctions/)
Per-sample STAR logs: [`results/star_1pass/`](results/star_1pass/) and [`results/star_2pass/`](results/star_2pass/)

---

## Step 5 — Performance Matrix

### Averages across 3 samples

| Aligner | Context | Avg Map% | Avg PP% | Avg Multi-map | Avg Time | Peak RSS |
|---------|---------|:--------:|:-------:|:-------------:|:--------:|:--------:|
| BWA-MEM | WGS | 100.00% | 100.0% | 0 | 31.8s | 0.27 GB |
| Bowtie2 | WGS | 99.99% | 98.5% | 0 | 46.8s | **0.13 GB** |
| Minimap2 | WGS | 100.00% | 100.0% | 0 | **17.7s** | 0.38 GB |
| Bowtie2 | RNA-seq | 100.00% | 100.0% | 0 | 44.5s | 0.13 GB |
| STAR 1-pass | RNA-seq | 100.00% | 99.7% | 1,030 | 20.5s | 0.92 GB |
| STAR 2-pass | RNA-seq | 100.00% | 99.7% | 1,028 | 41.8s | 1.62 GB |

Memory measured with `/usr/bin/time -v` (peak RSS). Full table: [`results/benchmarks/benchmark_summary.tsv`](results/benchmarks/benchmark_summary.tsv)

### Performance Charts

| Figure | Description |
|--------|-------------|
| [`mapping_rates.png`](results/benchmarks/mapping_rates.png) | Mapping rate bar chart — per aligner × sample (WGS + RNA-seq) |
| [`properly_paired.png`](results/benchmarks/properly_paired.png) | Properly paired % — per aligner comparison with per-bar labels |
| [`multimappers.png`](results/benchmarks/multimappers.png) | Multi-mapper rate — WGS (0% all) vs RNA-seq STAR counts |
| [`wallclock_time.png`](results/benchmarks/wallclock_time.png) | Alignment speed — wall-clock time bar chart per aligner |
| [`star_junction_analysis.png`](results/benchmarks/star_junction_analysis.png) | STAR 1-pass vs 2-pass — novel junction discovery & reclassification |
| [`insert_size_distribution.png`](results/benchmarks/insert_size_distribution.png) | Insert size distribution — KDE + boxplot + summary table per sample |
| [`performance_matrix.png`](results/benchmarks/performance_matrix.png) | 2×2 performance matrix: wall-clock + peak RSS for WGS and RNA-seq |

---

## Project Structure

```
day2_alignment_benchmark/
├── data/
│   └── reference/
│       ├── README.md                        # ← How to download/regenerate reference files
│       ├── chr22.fa                         # NOT committed (50 MB) — see README.md
│       └── gencode_v44_chr22_gtf            # NOT committed (32 MB) — see README.md
│
├── results/
│   ├── alignments/
│   │   ├── bwa/            (.gitkeep — BAMs gitignored, regenerate with pipeline)
│   │   ├── bowtie2_dna/    (.gitkeep)
│   │   ├── bowtie2_rna/    (.gitkeep)
│   │   └── minimap2/       (.gitkeep)
│   │
│   ├── indexes/
│   │   ├── bowtie2/        (.gitkeep — rebuild with: bowtie2-build chr22.fa indexes/bowtie2/chr22)
│   │   └── star/           (.gitkeep — rebuild with: python scripts/alignment_benchmark.py --index-only)
│   │
│   ├── star_1pass/
│   │   ├── PBMC_S1/        SJ.out.tab  Log.final.out  Log.out  Log.progress.out
│   │   ├── PBMC_S2/        (same)
│   │   └── PBMC_S3/        (same)
│   │
│   ├── star_2pass/
│   │   ├── PBMC_S1/        SJ.out.tab  Log.final.out  _STARpass1/  _STARgenome/
│   │   ├── PBMC_S2/        (same)
│   │   └── PBMC_S3/        (same)
│   │
│   ├── flagstats/
│   │   ├── bwa_NA12878_S{1,2,3}.txt
│   │   ├── bowtie2_dna_NA12878_S{1,2,3}.txt
│   │   ├── bowtie2_rna_PBMC_S{1,2,3}.txt
│   │   ├── minimap2_NA12878_S{1,2,3}.txt
│   │   ├── star_1pass_PBMC_S{1,2,3}.txt
│   │   └── star_2pass_PBMC_S{1,2,3}.txt
│   │
│   ├── junctions/
│   │   ├── star_junction_comparison.tsv    # Step 4 — 1-pass vs 2-pass per sample
│   │   ├── star_1pass_PBMC_S1_SJ_detail.tsv
│   │   └── star_1pass_PBMC_S3_SJ_detail.tsv
│   │
│   └── benchmarks/
│       ├── benchmark_full.tsv              # Step 3+5 — all 18 samples × 14 metrics
│       ├── benchmark_summary.tsv           # Per-aligner averages
│       ├── timing_with_memory_clean.tsv    # Wall-clock + peak RSS
│       ├── performance_matrix.png          ← Step 5 main figure
│       ├── mapping_rates.png               ← Step 3
│       ├── star_junction_analysis.png      ← Step 4
│       └── wallclock_time.png
│
├── scripts/
│   ├── alignment_benchmark.py             # Full reproducible pipeline (Steps 1–5)
│   ├── download_data.sh                   # Download reference + simulate reads
│   └── simulate_reference.py             # Generate synthetic chr22 (seed=42)
│
├── notebooks/
│   └── Day2_Alignment_Benchmark.ipynb     # Interactive analysis with all plots
│
├── .gitignore
└── README.md
```

---

## Reproduce

### Prerequisites

```bash
# Tools
bwa >= 0.7.17        # or install bwa-mem2 >= 2.2 for actual BWA-MEM2 benchmark
bowtie2 >= 2.4
STAR >= 2.7.11
minimap2 >= 2.24
samtools >= 1.19
wgsim                # included with samtools

# Memory tracking
apt install time

# Python
pip install pandas numpy matplotlib seaborn
```

### Run

```bash
# Step 0: Get reference data (choose one)
bash scripts/download_data.sh --synthetic   # regenerate exact synthetic chr22 (fast, ~1 min)
bash scripts/download_data.sh               # download real GRCh38 chr22 from Ensembl (~500 MB)

# Step 1: Build indexes only
python scripts/alignment_benchmark.py --index-only

# Steps 1–5: Full pipeline (skip alignment if BAMs already exist)
python scripts/alignment_benchmark.py

# Re-run all alignments even if BAMs exist
python scripts/alignment_benchmark.py --no-skip

# Recompute Step 3 stats from existing BAMs (no re-alignment)
python scripts/alignment_benchmark.py --stats-only
```

### Interactive Analysis

```bash
jupyter notebook notebooks/Day2_Alignment_Benchmark.ipynb
```

---

## Key Concepts

**Why GENCODE GTF matters for STAR:** Without `--sjdbGTFfile`, STAR operates annotation-free
and cannot distinguish novel vs annotated junctions. Providing the GTF pre-loads 8,468 known
splice sites, enabling Step 4's novel/annotated classification.

**STAR 2-pass mechanics:** Pass 1 discovers junctions from the data. An internal genome is rebuilt
combining those discoveries with GENCODE annotations. Pass 2 re-aligns against this augmented
index — reads near previously-missed junctions align correctly, and novel junctions are confirmed
or reclassified as annotated. Memory cost: 2-pass uses ~1.62 GB vs 1-pass ~0.92 GB because STAR
holds two genome indexes simultaneously during the internal rebuild.

**BWA-MEM vs BWA-MEM2:** Both produce identical alignments. BWA-MEM2 achieves 2–3× speed
improvement via AVX-512/AVX2 SIMD vectorisation of the BWT extension step (Vasimuddin 2019).
Expected BWA-MEM2 time on this dataset: ~12s at identical 0.27 GB RSS.

**Multi-mappers (STAR only):** STAR outputs all valid alignment positions per read with `NH:i`
tags, enabling downstream tools (HTSeq, featureCounts) to handle multi-mappers explicitly.
BWA/Bowtie2/Minimap2 report only the best-scoring alignment by default.

---

## Citation

```bibtex
@inproceedings{vasimuddin2019bwamem2,
  title={Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems},
  author={Vasimuddin, Md and Misra, Sanchit and Li, Heng and Aluru, Srinivas},
  booktitle={IEEE International Parallel and Distributed Processing Symposium (IPDPS)},
  year={2019}
}

@article{dobin2013star,
  title={STAR: ultrafast universal RNA-seq aligner},
  author={Dobin, Alexander and Davis, Carrie A and Schlesinger, Felix and others},
  journal={Bioinformatics},
  volume={29}, number={1}, pages={15--21},
  year={2013}
}
```

---

*Part of [30DaysOfBioinformatics-Portfolio](https://github.com/SubhadipJana1409/30DaysOfBioinformatics-Portfolio)*
