#!/usr/bin/env python3
"""
Day 2 — Reference Genome Alignment Benchmarking
BWA-MEM vs Bowtie2 vs STAR (1-pass + 2-pass) vs Minimap2

Steps from spec:
  1. Build indexes: BWA, Bowtie2, STAR with GENCODE v44 GTF
  2. Align 3 WGS (BWA-MEM, Bowtie2, Minimap2) + 3 RNA-seq samples (STAR 1p/2p, Bowtie2)
  3. Compute: mapping rate, properly paired %, chimeric rate, multi-mappers
  4. STAR 2-pass vs 1-pass: novel junction discovery & splice site reclassification
  5. Benchmark wall-clock time + peak RSS; plot performance matrix

Reference: Vasimuddin et al. (2019) BWA-MEM2, IEEE IPDPS
           Dobin et al. (2013) STAR, Bioinformatics

#30DaysOfBioinformatics | SubhadipJana1409
"""

import os, re, subprocess, sys, time, csv
import argparse
from pathlib import Path

# ─────────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────────
PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR    = PROJECT_DIR / "data"
REF_DIR     = DATA_DIR / "reference"
WGS_DIR     = DATA_DIR / "wgs"
RNA_DIR     = DATA_DIR / "rnaseq"
RES_DIR     = PROJECT_DIR / "results"
BENCH_DIR   = RES_DIR / "benchmarks"
FLAG_DIR    = RES_DIR / "flagstats"
ALIGN_DIR   = RES_DIR / "alignments"

REF_FA      = REF_DIR / "chr22.fa"
GTF_FILE    = REF_DIR / "gencode_v44_chr22_gtf"
BT2_IDX     = RES_DIR / "indexes" / "bowtie2" / "chr22"
STAR_IDX    = RES_DIR / "indexes" / "star"

WGS_SAMPLES = ["NA12878_S1", "NA12878_S2", "NA12878_S3"]
RNA_SAMPLES = ["PBMC_S1", "PBMC_S2", "PBMC_S3"]
THREADS = 4


# ─────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────
def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r.stdout, r.stderr, r.returncode


def run_timed(cmd):
    """Wrap with /usr/bin/time -v; return (stdout, wall_s, rss_mb, rc)."""
    full = f"/usr/bin/time -v sh -c '{cmd}'"
    r = subprocess.run(full, shell=True, capture_output=True, text=True)
    wall_s = _parse_wall(r.stderr)
    rss_mb = _parse_rss(r.stderr)
    return r.stdout, wall_s, rss_mb, r.returncode


def _parse_wall(txt):
    m = re.search(r"Elapsed.*?: ([\d:]+\.?\d*)", txt)
    if not m:
        return 0.0
    parts = m.group(1).split(":")
    return round(float(parts[0]) * 60 + float(parts[1]), 1)


def _parse_rss(txt):
    m = re.search(r"Maximum resident.*?(\d+)", txt)
    return round(int(m.group(1)) / 1024, 1) if m else 0.0


def samtools_index(bam):
    run(f"samtools index {bam}")


def flagstat(bam):
    """Parse samtools flagstat into dict."""
    out, _, _ = run(f"samtools flagstat {bam}")
    d = {}
    for line in out.splitlines():
        if "in total" in line:
            d["total"] = int(line.split()[0])
        elif "mapped (" in line and "primary" not in line:
            d["mapped"]  = int(line.split()[0])
            d["map_pct"] = float(line.split("(")[1].split("%")[0])
        elif "properly paired" in line:
            d["properly_paired"] = int(line.split()[0])
        elif "singletons" in line:
            d["singletons"] = int(line.split()[0])
        elif "with mate mapped to a different chr" in line and "mapQ" not in line:
            d["chimeric"] = int(line.split()[0])
    d.setdefault("chimeric", 0)
    d.setdefault("map_pct", 0.0)
    return d


def count_multimappers(bam):
    """Count primary aligned reads with NH:i > 1."""
    out, _, _ = run(
        f"samtools view -F 4 {bam} "
        f"| grep -oP 'NH:i:\\K\\d+' "
        f"| awk '$1>1{{c++}} END{{print c+0}}'"
    )
    try:
        return int(out.strip())
    except ValueError:
        return 0


def make_dirs():
    for p in [
        ALIGN_DIR / "bwa", ALIGN_DIR / "bowtie2_dna", ALIGN_DIR / "bowtie2_rna",
        ALIGN_DIR / "minimap2",
        BENCH_DIR, FLAG_DIR,
        RES_DIR / "indexes" / "bowtie2",
        RES_DIR / "indexes" / "star",
        RES_DIR / "star_1pass",
        RES_DIR / "star_2pass",
        RES_DIR / "junctions",
    ]:
        p.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────
# STEP 1 — INDEXES
# ─────────────────────────────────────────────────────────────
def build_indexes():
    print("\n[STEP 1] Building aligner indexes")
    print("─" * 60)

    # BWA (proxy for BWA-MEM2 — see README)
    if not (REF_DIR / "chr22.fa.amb").exists():
        print("  BWA index...")
        t0 = time.time()
        run(f"bwa index {REF_FA}")
        print(f"  Done in {time.time()-t0:.1f}s")
    else:
        print("  BWA index: already built")

    # Bowtie2
    if not Path(str(BT2_IDX) + ".1.bt2").exists():
        print("  Bowtie2 index...")
        t0 = time.time()
        run(f"bowtie2-build {REF_FA} {BT2_IDX} 2>/dev/null")
        print(f"  Done in {time.time()-t0:.1f}s")
    else:
        print("  Bowtie2 index: already built")

    # STAR with GENCODE v44
    if not (STAR_IDX / "Genome").exists():
        gtf_arg = ""
        if GTF_FILE.exists():
            n_lines, _, _ = run(f"grep -v '^#' {GTF_FILE} | wc -l")
            print(f"  STAR index + GENCODE v44 GTF ({n_lines.strip()} annotation lines, sjdbOverhang=149)...")
            gtf_arg = f"--sjdbGTFfile {GTF_FILE} --sjdbOverhang 149"
        else:
            print("  STAR index (no GTF found — annotation-free mode)...")
        t0 = time.time()
        run(
            f"STAR --runMode genomeGenerate "
            f"--genomeDir {STAR_IDX} "
            f"--genomeFastaFiles {REF_FA} "
            f"--genomeSAindexNbases 10 "
            f"--runThreadN {THREADS} {gtf_arg} 2>/dev/null"
        )
        log = (STAR_IDX / "Log.out").read_text() if (STAR_IDX / "Log.out").exists() else ""
        m = re.search(r"(\d+) collapsed junctions", log)
        n_junc = m.group(1) if m else "0"
        print(f"  Done in {time.time()-t0:.1f}s — {n_junc} splice junctions loaded from GTF")
    else:
        log = (STAR_IDX / "Log.out").read_text() if (STAR_IDX / "Log.out").exists() else ""
        m = re.search(r"(\d+) collapsed junctions", log)
        n_junc = m.group(1) if m else "?"
        print(f"  STAR index: already built ({n_junc} annotated junctions)")

    print()


# ─────────────────────────────────────────────────────────────
# STEP 2 — ALIGNMENTS
# ─────────────────────────────────────────────────────────────
def align_bwa(sample, r1, r2):
    out = ALIGN_DIR / "bwa" / f"{sample}.bam"
    rg = f"@RG\\\\tID:{sample}\\\\tSM:{sample}\\\\tPL:ILLUMINA"
    cmd = (f"bwa mem -t {THREADS} -R '{rg}' {REF_FA} {r1} {r2} 2>/dev/null "
           f"| samtools sort -@{THREADS} -o {out}")
    _, wall, rss, _ = run_timed(cmd)
    samtools_index(out)
    return out, wall, rss


def align_bowtie2(sample, r1, r2, subdir):
    out = ALIGN_DIR / subdir / f"{sample}.bam"
    cmd = (f"bowtie2 -x {BT2_IDX} -1 {r1} -2 {r2} -p {THREADS} 2>/dev/null "
           f"| samtools sort -o {out}")
    _, wall, rss, _ = run_timed(cmd)
    samtools_index(out)
    return out, wall, rss


def align_minimap2(sample, r1, r2):
    out = ALIGN_DIR / "minimap2" / f"{sample}.bam"
    rg = f"@RG\\\\tID:{sample}\\\\tSM:{sample}"
    cmd = (f"minimap2 -ax sr -t {THREADS} -R '{rg}' {REF_FA} {r1} {r2} 2>/dev/null "
           f"| samtools sort -o {out}")
    _, wall, rss, _ = run_timed(cmd)
    samtools_index(out)
    return out, wall, rss


def align_star(sample, r1, r2, twopass=False):
    mode = "star_2pass" if twopass else "star_1pass"
    out_dir = RES_DIR / mode / sample
    out_dir.mkdir(parents=True, exist_ok=True)
    tp_flag = "--twopassMode Basic" if twopass else ""
    cmd = (
        f"STAR --runMode alignReads "
        f"--genomeDir {STAR_IDX} "
        f"--readFilesIn {r1} {r2} "
        f"--readFilesCommand zcat "
        f"--outSAMtype BAM SortedByCoordinate "
        f"--outFileNamePrefix {out_dir}/ "
        f"--runThreadN {THREADS} "
        f"--outSAMattributes NH HI AS NM MD "
        f"--outSAMstrandField intronMotif "
        f"--outFilterIntronMotifs RemoveNoncanonical "
        f"{tp_flag} 2>/dev/null"
    )
    _, wall, rss, _ = run_timed(cmd)
    bam = out_dir / "Aligned.sortedByCoord.out.bam"
    samtools_index(bam)
    return bam, wall, rss, out_dir


# ─────────────────────────────────────────────────────────────
# STEP 3 — STATS
# ─────────────────────────────────────────────────────────────
def collect_stats(label, sample, dtype, bam, wall, rss):
    fs  = flagstat(bam)
    mm  = count_multimappers(bam)
    tot = fs.get("total", 0)
    pp_pct   = round(fs.get("properly_paired", 0) * 100 / tot, 2) if tot else 0.0
    chim_pct = round(fs.get("chimeric",         0) * 100 / tot, 4) if tot else 0.0
    mm_pct   = round(mm * 100 / tot, 2)                             if tot else 0.0

    run(f"samtools flagstat {bam} > {FLAG_DIR}/{label}_{sample}.txt")

    rec = {
        "aligner": label, "sample": sample, "type": dtype,
        "wall_s": wall, "peak_rss_mb": rss,
        "total_reads": tot,
        "mapped": fs.get("mapped", 0),
        "map_pct": fs.get("map_pct", 0.0),
        "properly_paired": fs.get("properly_paired", 0),
        "pp_pct": pp_pct,
        "chimeric": fs.get("chimeric", 0),
        "chimeric_pct": chim_pct,
        "multi_mappers": mm,
        "mm_pct": mm_pct,
    }
    print(
        f"    {sample}: map={rec['map_pct']:.2f}%  PP={pp_pct:.1f}%  "
        f"chimeric={rec['chimeric']}  MM={mm}  "
        f"{wall:.1f}s  {rss:.0f}MB"
    )
    return rec


# ─────────────────────────────────────────────────────────────
# STEP 4 — JUNCTION ANALYSIS
# ─────────────────────────────────────────────────────────────
def parse_sj_tab(sj_path):
    """
    STAR SJ.out.tab columns:
      0:chrom 1:intron_start 2:intron_end 3:strand 4:intron_motif
      5:annotated(0=novel,1=annotated) 6:n_uniq 7:n_mm 8:max_overhang
    Motif: 1=GT/AG 2=CT/AC(GT/AG antisense) 3=GC/AG 4=CT/GC 5=AT/AC 6=GT/AT 0=non-canonical
    """
    stats = {"total": 0, "novel": 0, "annotated": 0, "canonical": 0, "noncanonical": 0, "rows": []}
    if not Path(sj_path).exists():
        return stats
    with open(sj_path) as f:
        for line in f:
            c = line.strip().split("\t")
            if len(c) < 9:
                continue
            stats["total"] += 1
            is_novel = c[5] == "0"
            motif    = int(c[4])
            is_canon = motif in (1, 2)
            stats["novel"]       += int(is_novel)
            stats["annotated"]   += int(not is_novel)
            stats["canonical"]   += int(is_canon)
            stats["noncanonical"] += int(not is_canon)
            stats["rows"].append({
                "chrom": c[0], "start": c[1], "end": c[2],
                "motif": motif, "annotated": int(c[5]),
                "novel": is_novel, "canonical": is_canon,
                "uniq_reads": int(c[6]), "mm_reads": int(c[7])
            })
    return stats


def junction_analysis(sj1, sj2):
    print("\n[STEP 4] STAR Junction Analysis — 1-pass vs 2-pass")
    print("─" * 70)
    hdr = f"  {'Sample':<14} {'1p_SJ':<8} {'2p_SJ':<8} {'1p_novel':<10} {'2p_novel':<10} {'reclassified':<14} {'canonical'}"
    print(hdr)
    print("  " + "-" * 66)

    rows = []
    for s in RNA_SAMPLES:
        s1 = sj1.get(s, {})
        s2 = sj2.get(s, {})
        reclassified = max(0, s1.get("novel", 0) - s2.get("novel", 0))
        row = {
            "sample": s,
            "pass1_total":    s1.get("total", 0),
            "pass2_total":    s2.get("total", 0),
            "pass1_novel":    s1.get("novel", 0),
            "pass2_novel":    s2.get("novel", 0),
            "pass1_annotated":  s1.get("annotated", 0),
            "pass2_annotated":  s2.get("annotated", 0),
            "reclassified":   reclassified,
            "pass1_canonical": s1.get("canonical", 0),
            "pass1_noncanonical": s1.get("noncanonical", 0),
        }
        rows.append(row)
        print(
            f"  {s:<14} {row['pass1_total']:<8} {row['pass2_total']:<8} "
            f"{row['pass1_novel']:<10} {row['pass2_novel']:<10} "
            f"{reclassified:<14} {row['pass1_canonical']}"
        )

    # Save junction comparison TSV
    junc_tsv = RES_DIR / "junctions" / "star_junction_comparison.tsv"
    with open(junc_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # Save per-sample SJ detail TSVs (for downstream inspection)
    for s in RNA_SAMPLES:
        detail_rows = sj1.get(s, {}).get("rows", [])
        if detail_rows:
            out_f = RES_DIR / "junctions" / f"star_1pass_{s}_SJ_detail.tsv"
            with open(out_f, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=detail_rows[0].keys(), delimiter="\t")
                w.writeheader()
                w.writerows(detail_rows)

    print(f"\n  Saved → {junc_tsv}")
    print(
        "\n  Note: 2-pass reclassifies 1-pass novel junctions as annotated by\n"
        "  re-indexing with GENCODE v44 + discovered junctions before pass 2.\n"
        "  On real PBMC RNA-seq: ~50k-100k junctions, 5-15% novel→annotated gain."
    )
    return rows


# ─────────────────────────────────────────────────────────────
# STEP 5 — BENCHMARK TABLES
# ─────────────────────────────────────────────────────────────
def save_benchmark(records):
    fields = [
        "aligner", "sample", "type", "wall_s", "peak_rss_mb",
        "total_reads", "mapped", "map_pct",
        "properly_paired", "pp_pct",
        "chimeric", "chimeric_pct",
        "multi_mappers", "mm_pct",
    ]
    # Full per-sample TSV
    tsv = BENCH_DIR / "benchmark_full.tsv"
    with open(tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(records)

    # Per-aligner summary
    from collections import defaultdict
    by_aligner = defaultdict(list)
    for r in records:
        by_aligner[(r["aligner"], r["type"])].append(r)

    summary_tsv = BENCH_DIR / "benchmark_summary.tsv"
    with open(summary_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["aligner","type","n","avg_map_pct","avg_pp_pct",
                    "avg_chimeric","avg_multimappers","avg_wall_s","avg_rss_mb","avg_rss_gb"])
        for (aligner, dtype), rows in sorted(by_aligner.items()):
            n = len(rows)
            avg_rss = round(sum(r["peak_rss_mb"] for r in rows) / n, 1)
            w.writerow([
                aligner, dtype, n,
                round(sum(r["map_pct"]      for r in rows) / n, 2),
                round(sum(r["pp_pct"]       for r in rows) / n, 2),
                round(sum(r["chimeric"]     for r in rows) / n, 1),
                round(sum(r["multi_mappers"]for r in rows) / n, 1),
                round(sum(r["wall_s"]       for r in rows) / n, 1),
                avg_rss,
                round(avg_rss / 1024, 3),
            ])

    print(f"\n[STEP 5] Benchmark tables saved:")
    print(f"  Full      → {tsv}")
    print(f"  Summary   → {summary_tsv}")
    return tsv, summary_tsv


def print_table(records):
    print("\n" + "=" * 90)
    print(f"{'BENCHMARK RESULTS — STEPS 3 + 5':^90}")
    print("=" * 90)
    print(
        f"  {'Aligner':<14} {'Sample':<14} {'Type':<9}"
        f" {'Map%':<8} {'PP%':<7} {'Chim':<7} {'MultiMap':<10}"
        f" {'Time(s)':<9} {'RSS(MB)'}"
    )
    print("  " + "-" * 86)
    for r in records:
        print(
            f"  {r['aligner']:<14} {r['sample']:<14} {r['type']:<9}"
            f" {r['map_pct']:<8.2f} {r['pp_pct']:<7.1f} {r['chimeric']:<7}"
            f" {r['multi_mappers']:<10} {r['wall_s']:<9.1f} {r['peak_rss_mb']:.0f}"
        )


# ─────────────────────────────────────────────────────────────
# PIPELINE
# ─────────────────────────────────────────────────────────────
def run_pipeline(skip_existing=True):
    make_dirs()
    build_indexes()

    records   = []
    sj1_data  = {}
    sj2_data  = {}

    # WGS aligners
    for label, fn, hint_dir in [
        ("bwa_mem",  lambda s: align_bwa(s, WGS_DIR/f"{s}_R1.fastq.gz", WGS_DIR/f"{s}_R2.fastq.gz"), "bwa"),
        ("bowtie2",  lambda s: align_bowtie2(s, WGS_DIR/f"{s}_R1.fastq.gz", WGS_DIR/f"{s}_R2.fastq.gz", "bowtie2_dna"), "bowtie2_dna"),
        ("minimap2", lambda s: align_minimap2(s, WGS_DIR/f"{s}_R1.fastq.gz", WGS_DIR/f"{s}_R2.fastq.gz"), "minimap2"),
    ]:
        print(f"\n[{label.upper()}] WGS")
        for s in WGS_SAMPLES:
            hint = ALIGN_DIR / hint_dir / f"{s}.bam"
            if skip_existing and hint.exists():
                bam, wall, rss = hint, 0.0, 0.0
            else:
                bam, wall, rss = fn(s)
            records.append(collect_stats(label, s, "WGS", bam, wall, rss))

    # Bowtie2 RNA-seq (splice-unaware baseline)
    print("\n[BOWTIE2_RNA] RNA-seq splice-unaware baseline")
    for s in RNA_SAMPLES:
        hint = ALIGN_DIR / "bowtie2_rna" / f"{s}.bam"
        if skip_existing and hint.exists():
            bam, wall, rss = hint, 0.0, 0.0
        else:
            bam, wall, rss = align_bowtie2(s, RNA_DIR/f"{s}_R1.fastq.gz", RNA_DIR/f"{s}_R2.fastq.gz", "bowtie2_rna")
        records.append(collect_stats("bowtie2_rna", s, "RNA-seq", bam, wall, rss))

    # STAR 1-pass
    print("\n[STAR_1PASS] GTF-informed, splice-aware")
    for s in RNA_SAMPLES:
        hint = RES_DIR / "star_1pass" / s / "Aligned.sortedByCoord.out.bam"
        if skip_existing and hint.exists():
            bam, wall, rss, out_dir = hint, 0.0, 0.0, RES_DIR / "star_1pass" / s
        else:
            bam, wall, rss, out_dir = align_star(s, RNA_DIR/f"{s}_R1.fastq.gz", RNA_DIR/f"{s}_R2.fastq.gz", twopass=False)
        records.append(collect_stats("star_1pass", s, "RNA-seq", bam, wall, rss))
        sj1_data[s] = parse_sj_tab(out_dir / "SJ.out.tab")

    # STAR 2-pass
    print("\n[STAR_2PASS] twopassMode Basic — novel junction recovery")
    for s in RNA_SAMPLES:
        hint = RES_DIR / "star_2pass" / s / "Aligned.sortedByCoord.out.bam"
        if skip_existing and hint.exists():
            bam, wall, rss, out_dir = hint, 0.0, 0.0, RES_DIR / "star_2pass" / s
        else:
            bam, wall, rss, out_dir = align_star(s, RNA_DIR/f"{s}_R1.fastq.gz", RNA_DIR/f"{s}_R2.fastq.gz", twopass=True)
        records.append(collect_stats("star_2pass", s, "RNA-seq", bam, wall, rss))
        sj2_data[s] = parse_sj_tab(out_dir / "SJ.out.tab")

    junction_analysis(sj1_data, sj2_data)
    print_table(records)
    save_benchmark(records)
    return records


# ─────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Day 2: Alignment Benchmark (Steps 1-5)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Usage examples:
  python alignment_benchmark.py              # full run (skip if BAMs exist)
  python alignment_benchmark.py --no-skip    # re-run all alignments
  python alignment_benchmark.py --index-only # build indexes only
  python alignment_benchmark.py --stats-only # recompute Step3 stats from existing BAMs
""")
    p.add_argument("--index-only",  action="store_true")
    p.add_argument("--no-skip",     action="store_true")
    p.add_argument("--stats-only",  action="store_true")
    args = p.parse_args()

    if args.index_only:
        make_dirs(); build_indexes()
    elif args.stats_only:
        make_dirs()
        bam_map = [
            ("bwa_mem",     "NA12878_S1", "WGS",     ALIGN_DIR/"bwa"/"NA12878_S1.bam"),
            ("bwa_mem",     "NA12878_S2", "WGS",     ALIGN_DIR/"bwa"/"NA12878_S2.bam"),
            ("bwa_mem",     "NA12878_S3", "WGS",     ALIGN_DIR/"bwa"/"NA12878_S3.bam"),
            ("bowtie2",     "NA12878_S1", "WGS",     ALIGN_DIR/"bowtie2_dna"/"NA12878_S1.bam"),
            ("bowtie2",     "NA12878_S2", "WGS",     ALIGN_DIR/"bowtie2_dna"/"NA12878_S2.bam"),
            ("bowtie2",     "NA12878_S3", "WGS",     ALIGN_DIR/"bowtie2_dna"/"NA12878_S3.bam"),
            ("minimap2",    "NA12878_S1", "WGS",     ALIGN_DIR/"minimap2"/"NA12878_S1.bam"),
            ("minimap2",    "NA12878_S2", "WGS",     ALIGN_DIR/"minimap2"/"NA12878_S2.bam"),
            ("minimap2",    "NA12878_S3", "WGS",     ALIGN_DIR/"minimap2"/"NA12878_S3.bam"),
            ("bowtie2_rna", "PBMC_S1", "RNA-seq",   ALIGN_DIR/"bowtie2_rna"/"PBMC_S1.bam"),
            ("bowtie2_rna", "PBMC_S2", "RNA-seq",   ALIGN_DIR/"bowtie2_rna"/"PBMC_S2.bam"),
            ("bowtie2_rna", "PBMC_S3", "RNA-seq",   ALIGN_DIR/"bowtie2_rna"/"PBMC_S3.bam"),
            ("star_1pass",  "PBMC_S1", "RNA-seq",   RES_DIR/"star_1pass"/"PBMC_S1"/"Aligned.sortedByCoord.out.bam"),
            ("star_1pass",  "PBMC_S2", "RNA-seq",   RES_DIR/"star_1pass"/"PBMC_S2"/"Aligned.sortedByCoord.out.bam"),
            ("star_1pass",  "PBMC_S3", "RNA-seq",   RES_DIR/"star_1pass"/"PBMC_S3"/"Aligned.sortedByCoord.out.bam"),
            ("star_2pass",  "PBMC_S1", "RNA-seq",   RES_DIR/"star_2pass"/"PBMC_S1"/"Aligned.sortedByCoord.out.bam"),
            ("star_2pass",  "PBMC_S2", "RNA-seq",   RES_DIR/"star_2pass"/"PBMC_S2"/"Aligned.sortedByCoord.out.bam"),
            ("star_2pass",  "PBMC_S3", "RNA-seq",   RES_DIR/"star_2pass"/"PBMC_S3"/"Aligned.sortedByCoord.out.bam"),
        ]
        recs = []
        for label, sample, dtype, bam in bam_map:
            if Path(bam).exists():
                recs.append(collect_stats(label, sample, dtype, bam, 0.0, 0.0))
        print_table(recs)
        save_benchmark(recs)
    else:
        run_pipeline(skip_existing=not args.no_skip)
