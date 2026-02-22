# LIR Motif Finder

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Find **LC3-Interacting Region (LIR) motifs** within **intrinsically disordered regions (IDRs)** of protein sequences.  
Disorder prediction is powered by [metapredict V2](https://github.com/idptools/metapredict).

---

## Background

LIR motifs (also called AIM motifs in yeast) mediate selective autophagy by binding ATG8-family proteins (LC3/GABARAP).  
The canonical LIR core is **[W/F/Y]-x-x-[L/I/V]**, and functional LIRs are typically found in intrinsically disordered regions of proteins.

This tool:
1. Predicts per-residue disorder scores with **metapredict**
2. Identifies contiguous disordered regions above a user-defined threshold
3. Scans those regions for three classes of LIR motifs
4. Outputs a structured CSV/TSV table and a plain-text summary

### Motif classes

| Type | Pattern | Description |
|------|---------|-------------|
| `basic` | `[WFY]xx[LIV]` | Canonical LIR core |
| `extended` | `[DE][WFY]xx[LIV]` | With single upstream acidic residue |
| `acidic_extended` | `[DE]{2}[WFY]xx[LIV]` | With two upstream acidic residues |

---

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/lirmotiffinder.git
cd lirmotiffinder

# Install (editable mode recommended for development)
pip install -e .

# Or install directly from GitHub
pip install git+https://github.com/yourusername/lirmotiffinder.git
```

**Dependencies** (installed automatically):
- `biopython >= 1.79`
- `metapredict >= 2.0`

---

## Quick start

### Command line

```bash
# Run with defaults (disorder threshold 0.5, min region 5 residues)
lirmotiffinder proteins.fasta

# Custom thresholds, write to a custom output directory
lirmotiffinder proteins.fasta --threshold 0.6 --min-disorder 10 --outdir my_results

# Search only basic LIR motifs, output as TSV
lirmotiffinder proteins.fasta --motif-types basic --format tsv

# Suppress per-hit stdout print (only write files)
lirmotiffinder proteins.fasta --quiet
```

Output files in `results/` (or `--outdir`):

| File | Contents |
|------|----------|
| `lir_hits.csv` | One row per motif hit |
| `summary.txt` | Plain-text summary statistics |

### Python API

```python
from lirmotiffinder import run_analysis

results = run_analysis(
    fasta_path="proteins.fasta",
    threshold=0.5,          # metapredict disorder threshold
    min_disorder_length=5,  # minimum IDR length (residues)
    motif_types=None,       # None → search all types
)

for result in results:
    print(f"{result.protein_id}: {len(result.motif_hits)} hit(s)")
    for hit in result.motif_hits:
        print(
            f"  [{hit.motif_type}] {hit.motif_sequence} "
            f"@ {hit.motif_start_protein}-{hit.motif_end_protein} "
            f"(IDR {hit.disorder_start_protein}-{hit.disorder_end_protein})"
        )
```

---

## Output columns (CSV/TSV)

| Column | Description |
|--------|-------------|
| `protein_id` | FASTA record ID |
| `protein_description` | Full FASTA header |
| `protein_length` | Total sequence length |
| `motif_type` | `basic`, `extended`, or `acidic_extended` |
| `motif_sequence` | Matched motif sequence |
| `motif_start_protein` | 1-based start position in full protein |
| `motif_end_protein` | 1-based end position in full protein |
| `disorder_start_protein` | 1-based start of containing IDR |
| `disorder_end_protein` | 1-based end of containing IDR |
| `disorder_length` | Length of containing IDR (residues) |
| `mean_disorder_score` | Mean metapredict score over motif residues |

---

## CLI reference

```
usage: lirmotiffinder [-h] [--threshold FLOAT] [--min-disorder INT]
                      [--motif-types TYPE [TYPE ...]]
                      [--outdir DIR] [--format {csv,tsv}]
                      [--no-summary] [--quiet]
                      FASTA

positional arguments:
  FASTA                 Input FASTA file (protein sequences)

options:
  --threshold FLOAT     Disorder score threshold (default: 0.5)
  --min-disorder INT    Minimum IDR length in residues (default: 5)
  --motif-types TYPE    Space-separated list from: basic extended acidic_extended
  --outdir DIR          Output directory (default: ./results/)
  --format {csv,tsv}    Output format (default: csv)
  --no-summary          Skip writing summary.txt
  --quiet               Suppress per-hit stdout output
```

---

## Running tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

---

## Project structure

```
lirmotiffinder/
├── lirmotiffinder/
│   ├── __init__.py      # Public API
│   ├── core.py          # Disorder prediction + motif search logic
│   └── cli.py           # Command-line interface
├── tests/
│   ├── conftest.py      # Shared fixtures (metapredict mock)
│   └── test_core.py     # Unit tests
├── pyproject.toml
└── README.md
```

---

## Citation

If you use this tool in your research, please cite metapredict:

> Emenecker R.J., Griffith D., Holehouse A.S. (2021).  
> Metapredict: a fast, accurate, and easy-to-use predictor of consensus disorder and structure.  
> *Biophysical Journal*, 120(20), 4312–4319.  
> https://doi.org/10.1016/j.bpj.2021.08.039

---

## License

MIT
