#!/usr/bin/env python3
"""
simulate_reference.py — Generate synthetic chr22 benchmark reference

Produces the exact 5 Mb chr22.fa used in the Day 2 alignment benchmark.
Uses random.seed(42) for full reproducibility.

Output: data/reference/chr22.fa

Usage:
    python scripts/simulate_reference.py
    python scripts/simulate_reference.py --length 5000000 --out data/reference/chr22.fa

#30DaysOfBioinformatics | SubhadipJana1409
"""

import random
import argparse
from pathlib import Path

def simulate_chr22(length: int = 5_000_000, n_fraction: float = 0.05, seed: int = 42) -> str:
    """
    Generate a synthetic chr22 sequence.

    Args:
        length:     Total sequence length in bp (default 5 Mb)
        n_fraction: Fraction of bases to set as N (centromere/telomere simulation)
        seed:       Random seed for reproducibility

    Returns:
        DNA sequence string
    """
    random.seed(seed)
    bases = ["A", "T", "G", "C"]

    # Generate base sequence
    seq = [random.choice(bases) for _ in range(length)]

    # Insert N-regions to simulate centromere/heterochromatic regions
    # (~5% of chr22 is N in GRCh38)
    n_bases = int(length * n_fraction)
    n_positions = random.sample(range(length), n_bases)
    for pos in n_positions:
        seq[pos] = "N"

    return "".join(seq)


def write_fasta(seq: str, out_path: Path, header: str = "chr22", line_width: int = 60):
    """Write sequence to FASTA format."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f">{header} Homo sapiens chromosome 22 (synthetic benchmark, GRCh38-coord-compatible)\n")
        for i in range(0, len(seq), line_width):
            f.write(seq[i:i + line_width] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic chr22 benchmark reference genome"
    )
    parser.add_argument(
        "--length", type=int, default=5_000_000,
        help="Sequence length in bp (default: 5,000,000)"
    )
    parser.add_argument(
        "--n-fraction", type=float, default=0.05,
        help="Fraction of N bases (default: 0.05)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed (default: 42)"
    )
    parser.add_argument(
        "--out", type=str,
        default=str(Path(__file__).resolve().parent.parent / "data" / "reference" / "chr22.fa"),
        help="Output FASTA path"
    )
    args = parser.parse_args()

    out_path = Path(args.out)
    print(f"Generating synthetic chr22...")
    print(f"  Length:     {args.length:,} bp")
    print(f"  N fraction: {args.n_fraction:.1%}")
    print(f"  Seed:       {args.seed}")
    print(f"  Output:     {out_path}")

    seq = simulate_chr22(
        length=args.length,
        n_fraction=args.n_fraction,
        seed=args.seed
    )

    write_fasta(seq, out_path)

    # Report composition
    from collections import Counter
    counts = Counter(seq)
    total = len(seq)
    print(f"\nSequence composition:")
    for base in ["A", "T", "G", "C", "N"]:
        pct = counts[base] / total * 100
        print(f"  {base}: {counts[base]:>8,}  ({pct:.2f}%)")

    size_mb = out_path.stat().st_size / 1024 / 1024
    print(f"\nDone → {out_path} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
