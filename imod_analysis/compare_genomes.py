#!/usr/bin/env python3
"""
compare_assemblies.py

Given two genome assemblies (FASTA) and their annotation files (GFF),
compute several quick similarity statistics:

1. Assembly statistics – size, #contigs, N50, GC%.
2. Average nucleotide identity (ANI) using Mash, FastANI, or dnadiff if available.
3. Gene content overlap – Jaccard index of gene IDs / product names.

Usage:
    python compare_assemblies.py A.fna A.gff B.fna B.gff [-o report.txt]

External dependencies (optional but strongly recommended):
    - Mash  (https://mash.readthedocs.io)  OR
    - fastANI (https://github.com/ParBLiSS/fastANI) OR
    - MUMmer's dnadiff (https://mummer4.github.io/)
Install at least one of the above and make sure the executable is in $PATH.

Author: ChatGPT
Date: 2025‑06‑04
"""
import argparse, shutil, subprocess, tempfile, textwrap, os, sys
from collections import Counter
from statistics import mean, median

from typing import List, Tuple

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("Biopython is required: pip install biopython")

import pandas as pd


def seq_stats(fasta: str):
    lengths = []
    gcs = []
    for rec in SeqIO.parse(fasta, "fasta"):
        seq = rec.seq.upper()
        lengths.append(len(seq))
        gcs.append((seq.count("G") + seq.count("C")) / len(seq) * 100)
    total = sum(lengths)
    n_contigs = len(lengths)
    gc = mean(gcs) if gcs else 0
    # N50
    sorted_lens = sorted(lengths, reverse=True)
    csum = 0
    n50 = 0
    for l in sorted_lens:
        csum += l
        if csum >= total / 2:
            n50 = l
            break
    return {"size": total, "contigs": n_contigs, "N50": n50, "GC": gc}


def mash_ani(fasta1: str, fasta2: str):
    if shutil.which("mash") is None:
        return None
    tempdir = tempfile.mkdtemp()
    try:
        sketch1 = os.path.join(tempdir, "1.msh")
        sketch2 = os.path.join(tempdir, "2.msh")
        # Create sketches
        subprocess.run(["mash", "sketch", "-o", sketch1[:-4], fasta1], check=True, stdout=subprocess.PIPE)
        subprocess.run(["mash", "sketch", "-o", sketch2[:-4], fasta2], check=True, stdout=subprocess.PIPE)
        # Mash dist
        res = subprocess.run(["mash", "dist", sketch1, sketch2], check=True, capture_output=True, text=True)
        line = res.stdout.strip().split("\t")
        if len(line) >= 3:
            mash_dist = float(line[2])  # Mash distance ~ 1-ANI
            ani = (1 - mash_dist) * 100
            return ani
    except subprocess.CalledProcessError as e:
        print("Mash failed:", e, file=sys.stderr)
    finally:
        shutil.rmtree(tempdir, ignore_errors=True)
    return None


def fastani_ani(fasta1: str, fasta2: str):
    if shutil.which("fastani") is None:
        return None
    temp_out = tempfile.NamedTemporaryFile(delete=False)
    temp_out.close()
    try:
        res = subprocess.run(
            ["fastANI", "-q", fasta1, "-r", fasta2, "-o", temp_out.name],
            check=True,
            capture_output=True,
            text=True,
        )
        with open(temp_out.name) as f:
            line = f.readline().strip().split("\t")
            if len(line) >= 3:
                ani = float(line[2])
                return ani
    except subprocess.CalledProcessError as e:
        print("fastANI failed:", e, file=sys.stderr)
    finally:
        os.unlink(temp_out.name)
    return None


def parse_gff_genes(gff: str) -> List[str]:
    genes = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2].lower()
            if feature not in {"gene", "cds"}:
                continue
            attr = parts[8]
            attrs = {}
            for field in attr.split(";"):
                if "=" in field:
                    k, v = field.split("=", 1)
                    attrs[k] = v
            gene_id = attrs.get("ID") or attrs.get("locus_tag") or attrs.get("Name") or ""
            product = attrs.get("product", "")
            key = gene_id or product
            if key:
                genes.append(key)
    return genes


def gene_content_similarity(gff1: str, gff2: str):
    genes1 = set(parse_gff_genes(gff1))
    genes2 = set(parse_gff_genes(gff2))
    if not genes1 or not genes2:
        return None
    inter = len(genes1 & genes2)
    union = len(genes1 | genes2)
    return inter / union * 100


def reciprocal_best_hits(proteins_a: str, proteins_b: str, out_prefix: str = "rbh", evalue: float = 1e-5) -> pd.DataFrame:
    """
    Map genes between two genomes using reciprocal best BLASTP hits (RBH).
    proteins_a: FASTA file of proteins from genome A
    proteins_b: FASTA file of proteins from genome B
    out_prefix: prefix for temporary BLAST output files
    evalue: BLASTP e-value threshold
    Returns: DataFrame with columns [gene_a, gene_b, score_a2b, score_b2a]
    """
    import tempfile
    import os
    import subprocess
    from Bio import SeqIO

    with tempfile.TemporaryDirectory() as tmpdir:
        db_b = os.path.join(tmpdir, "b_db")
        db_a = os.path.join(tmpdir, "a_db")
        blast_a2b = os.path.join(tmpdir, f"{out_prefix}_a2b.tsv")
        blast_b2a = os.path.join(tmpdir, f"{out_prefix}_b2a.tsv")
        # Make BLAST databases
        subprocess.run(["makeblastdb", "-in", proteins_b, "-dbtype", "prot", "-out", db_b], check=True)
        subprocess.run(["makeblastdb", "-in", proteins_a, "-dbtype", "prot", "-out", db_a], check=True)
        # Run BLASTP: A vs B
        subprocess.run([
            "blastp", "-query", proteins_a, "-db", db_b, "-outfmt", "6 qseqid sseqid bitscore evalue", "-evalue", str(evalue), "-max_target_seqs", "1", "-out", blast_a2b
        ], check=True)
        # Run BLASTP: B vs A
        subprocess.run([
            "blastp", "-query", proteins_b, "-db", db_a, "-outfmt", "6 qseqid sseqid bitscore evalue", "-evalue", str(evalue), "-max_target_seqs", "1", "-out", blast_b2a
        ], check=True)
        # Parse BLAST outputs
        a2b = pd.read_csv(blast_a2b, sep="\t", names=["gene_a", "gene_b", "bitscore", "evalue"])
        b2a = pd.read_csv(blast_b2a, sep="\t", names=["gene_b", "gene_a", "bitscore", "evalue"])
        # Merge to find reciprocal best hits
        merged = pd.merge(a2b, b2a, on=["gene_a", "gene_b"])
        merged = merged[["gene_a", "gene_b", "bitscore_x", "bitscore_y"]]
        merged.columns = ["gene_a", "gene_b", "score_a2b", "score_b2a"]
        return merged


def main():
    p = argparse.ArgumentParser(description="Compare two genome assemblies (FASTA+GFF).")
    p.add_argument("fasta1")
    p.add_argument("gff1")
    p.add_argument("fasta2")
    p.add_argument("gff2")
    p.add_argument("-o", "--output", help="Write report to file (default: stdout)")
    args = p.parse_args()

    stats1 = seq_stats(args.fasta1)
    stats2 = seq_stats(args.fasta2)

    ani = mash_ani(args.fasta1, args.fasta2) or fastani_ani(args.fasta1, args.fasta2)

    gene_sim = gene_content_similarity(args.gff1, args.gff2)

    report = []
    report.append("Assembly comparison report\n")
    report.append("File A: {}\n".format(args.fasta1))
    report.append("File B: {}\n".format(args.fasta2))
    report.append("\n--- Assembly stats ---")
    report.append(f"\nA: Size={stats1['size']:,} bp  Contigs={stats1['contigs']}  N50={stats1['N50']:,}  GC%={stats1['GC']:.2f}")
    report.append(f"\nB: Size={stats2['size']:,} bp  Contigs={stats2['contigs']}  N50={stats2['N50']:,}  GC%={stats2['GC']:.2f}")

    if ani is not None:
        report.append(f"\n\n--- Nucleotide similarity ---\nAverage Nucleotide Identity (ANI): {ani:.2f}%")
    else:
        report.append("\n\n--- Nucleotide similarity ---\nANI: Could not compute (Mash or fastANI not found)")

    if gene_sim is not None:
        report.append(f"\n\n--- Gene content similarity ---\nJaccard similarity: {gene_sim:.2f}% ({gene_sim/100:.3f})")
    else:
        report.append("\n\n--- Gene content similarity ---\nCould not compute (no genes parsed)")

    text = "\n".join(report)
    if args.output:
        with open(args.output, "w") as out:
            out.write(text)
    else:
        print(text)


if __name__ == "__main__":
    main()