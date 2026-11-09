# Task 02: Database search and annotation (BLAST)

## What I did
Searched the NCBI non-redundant protein database with a query protein and annotated the best matches.

- **Query:** `query.fasta`, the phiX174 major capsid protein (gene F, 427 aa), translated from the Task 01 genome (NC_001422.1, nt 1001 onward until the first stop codon).
- **Why this query:** a real, known protein from a public genome I had already analysed, so I can check whether the hits make biological sense.

## Search parameters (from the BLAST "Search Summary", Request ID F359PR81016)
| Parameter | Value |
|---|---|
| Program | blastp |
| Database | nr (non-redundant protein sequences), posted 4 Nov 2026 |
| Database size | 1,207,322,105 sequences, 453,381,857,326 letters |
| Expect value | 0.05 |
| Max hits (hitlist size) | 100 |
| Word size | 5 |
| Matrix | BLOSUM62 |
| Gap costs | existence 11, extension 1 |
| Composition-based stats | 2 |
| Low-complexity filter | off (filter string F) |
| Window size / threshold | 40 / 0 |

## Files
| File | Purpose |
|---|---|
| `query.fasta` | query protein |
| `blast_hit_table.csv` | raw BLAST output (hit table, 100 hits, no header) |
| `annotate_hits.py` | script that builds the annotation table |
| `annotation_table.tsv` | annotation of the top 10 hits |
| `annotation_summary.txt` | table plus summary statistics printed by the script |

## How to rerun
```bash
python annotate_hits.py blast_hit_table.csv --qlen 427 --top 10 --out annotation_table.tsv
```
Python 3 only, no packages. To redo the search itself, paste `query.fasta` into NCBI Protein BLAST (blastp, nr) with the parameters above and download the Hit Table (CSV).

## Annotation table (top 10)
See `annotation_table.tsv` for the full table with the source and interpretation columns.

| Rank | Accession | Identity % | Query cover % | Mismatches | E-value | Bit score |
|---|---|---|---|---|---|---|
| 1 | WP_000065572.1 | 100.00 | 100 | 0 | 0.0 | 895 |
| 2 | WP_000065548.1 | 99.77 | 100 | 1 | 0.0 | 894 |
| 3 | ADC55535.1 | 99.77 | 100 | 1 | 0.0 | 892 |
| 4 | AXC43224.1 | 99.53 | 100 | 2 | 0.0 | 891 |
| 5 | ACO94986.1 | 99.53 | 100 | 2 | 0.0 | 891 |
| 6 | QIC52958.1 | 99.53 | 100 | 2 | 0.0 | 890 |
| 7 | 1AL0_F | 99.77 | 99.8 | 1 | 0.0 | 890 |
| 8 | ACY07120.1 | 99.53 | 100 | 2 | 0.0 | 890 |
| 9 | ADM34310.1 | 99.53 | 100 | 2 | 0.0 | 890 |
| 10 | ADM34695.1 | 99.30 | 100 | 3 | 0.0 | 889 |

## Interpretation
- **The best hit is the query itself.** WP_000065572.1 is 100% identical over the full 427 residues (bit score 895), as expected because the query came from phiX174. This also works as a sanity check that the pipeline is correct.
- **Every hit is the same protein.** All 100 hits have identity between 95.55% and 100% and an E-value of 0.0. At this level the E-value is not a useful way to rank hits, so I ranked by bit score and identity. The small differences (0 to 19 mismatches out of 427) are natural variants between phiX174-type phage isolates and strains, not different proteins.
- **Sources of the hits.** Of the top 100, 93 are GenBank submissions (many separate sequencing records of phiX174-like phages), 3 are RefSeq non-redundant entries (WP_), 3 are PDB structures (for example 1AL0_F, 1CD3_F, 9K3M_A) and 1 is a UniProt-style entry (P07931.1). The PDB hits show the capsid protein has a solved 3D structure, which supports the annotation as a structural capsid protein.
- **Partial hits.** EKY23402.1 and EPJ19293.1 (start at residue 11, 97.7% cover) and WP_068857377.1 (361 residues, 84.5% cover) are 99.7 to 100% identical over the stretch they cover. Their shorter coverage most likely means incomplete or truncated records, not real differences in the protein.
- **Limit of this search.** The 100-hit cap was filled entirely by near-identical sequences, so distant relatives (other microvirus capsid proteins, with perhaps 30 to 60% identity) did not appear. This is a consequence of the setting, not proof that none exist. To look for them I could exclude phiX174 organisms in the BLAST filter, or raise the hit limit.

## Decisions and limits
- I chose `nr` because it is the broadest protein database and the default.
- I left all other parameters at their defaults so that anyone can reproduce the search.
- The downloaded file was the Hit Table, which has no protein names or organisms. The script therefore classifies each hit by its accession type (RefSeq, PDB, UniProt, GenBank) and does not state organism names. Query coverage is computed by the script as aligned length divided by 427.
- The interpretation in the last column of `annotation_table.tsv` is produced by simple threshold rules (identity, coverage, E-value) written in `annotate_hits.py`. I read the results myself and the rules agree with the discussion above.
