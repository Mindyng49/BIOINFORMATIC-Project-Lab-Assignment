# Bioinformatics Lab 2 - Pairwise Sequence Alignment

Streamlit application that pairs Lab 1 (UniProt retrieval + protein analysis)
with Lab 2 (Needleman-Wunsch global alignment and Smith-Waterman local
alignment, implemented from scratch).

## Run

```bash
streamlit run lab2-Mindy.py
```

Requires `streamlit`, `biopython` and `pandas`. Only Biopython is used for the
Lab 1 data retrieval - the alignment code imports nothing outside the standard
library.

## Files

| File | Purpose |
| --- | --- |
| `lab2-Mindy.py` | Streamlit UI: navigation, inputs, scoring, result rendering |
| `pairwise_alignment.py` | DP matrix, trace-back, aligned strings, text export (pure Python) |
| `lab1_core.py` | Lab 1 functions (`retrieve_data`, `get_basic_analysis`) reused here |
| `lab1-Mindy.py` | Original Lab 1 app, left untouched |

## What the app shows

* alignment score (global = `D[n][m]`, local = highest cell in the matrix)
* the aligned sequences, colour-coded, in 60-column blocks
* the full dynamic-programming matrix with heat colours and the trace-back
  path highlighted (green = path start, amber = path, red = path end / max)
* the trace-back steps (cell, score, move, aligned column) in a table
* identity, matches, mismatches, gaps, alignment length
* downloads: alignment report (.txt) and complete matrix (.csv)

## Usage notes

* Sequences can be pasted raw or as a FASTA record; headers, digits and spaces
  are ignored.
* Global vs local and the match / mismatch / gap scores are set below the
  inputs, then **Run alignment**.
* In Lab 1, after retrieving a protein, **Use as Sequence 1 / 2** sends it
  straight to Lab 2 (also available from the *From Lab 1* dropdown there).
