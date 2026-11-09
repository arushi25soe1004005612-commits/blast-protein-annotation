#!/usr/bin/env python3
"""Build an annotation table from an NCBI BLAST hit table (CSV, outfmt 10).

Usage: python annotate_hits.py blast_hit_table.csv [--qlen 427] [--top 10] [--out annotation_table.tsv]

Columns of the hit table (no header): query, subject, %identity, aln_len,
mismatches, gap_opens, q_start, q_end, s_start, s_end, evalue, bitscore, %positives.
Query coverage is computed here as (q_end - q_start + 1) / qlen.
The hit table has no protein names, so the source database is inferred from the
accession prefix. Stdlib only.
"""
import argparse, csv, re
from collections import Counter

def source(acc):
    if acc.startswith("WP_"):
        return "RefSeq non-redundant protein (WP_): one record shared by many identical bacterial/phage entries"
    if re.match(r"^[0-9][A-Za-z0-9]{3}_", acc):
        return "PDB structure (experimentally solved 3D structure)"
    if re.match(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]", acc):
        return "UniProt/Swiss-Prot style accession (curated protein entry)"
    return "GenBank/ENA translated CDS (submitted sequence)"

def interpret(ident, cov, evalue):
    if evalue > 1e-5 or cov < 30:
        return "weak: not significant, no homology inferred"
    if ident >= 99 and cov >= 99:
        return "essentially identical full-length match: the same protein (phiX174-type capsid F)"
    if ident >= 95 and cov >= 99:
        return "near-identical full-length match: same protein, close natural variant"
    if ident >= 95:
        return "near-identical but partial (truncated record): same protein"
    if ident >= 60 and cov >= 70:
        return "close homolog: likely same function"
    if ident >= 30 and cov >= 50:
        return "distant homolog: shared family, function unconfirmed"
    return "partial match: only a region shared"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv"); ap.add_argument("--qlen", type=int, default=427)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--out", default="annotation_table.tsv")
    a = ap.parse_args()
    hits = []
    with open(a.csv, newline="") as f:
        for r in csv.reader(f):
            if len(r) < 12: continue
            hits.append(dict(acc=r[1], ident=float(r[2]), alen=int(r[3]), mism=int(r[4]),
                             qs=int(r[6]), qe=int(r[7]), ev=float(r[10]), bits=float(r[11]),
                             cov=100.0 * (int(r[7]) - int(r[6]) + 1) / a.qlen))
    cols = ["rank","accession","identity_%","query_cover_%","mismatches","evalue","bitscore","source","interpretation"]
    lines = ["\t".join(cols)]
    for i, h in enumerate(hits[:a.top], 1):
        lines.append("\t".join(map(str, [i, h["acc"], f'{h["ident"]:.2f}', f'{h["cov"]:.1f}', h["mism"],
                     h["ev"], h["bits"], source(h["acc"]), interpret(h["ident"], h["cov"], h["ev"])])))
    open(a.out, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    idents = [h["ident"] for h in hits]
    print(f"\n# hits: {len(hits)}  identity range: {min(idents):.2f}-{max(idents):.2f}%  "
          f"max e-value: {max(h['ev'] for h in hits)}")
    print("# hits by source:", dict(Counter(source(h["acc"]).split(" (")[0].split(":")[0] for h in hits)))
    print("# partial-length hits (cover < 99%):", [(h["acc"], round(h["cov"],1)) for h in hits if h["cov"] < 99])

if __name__ == "__main__":
    main()
