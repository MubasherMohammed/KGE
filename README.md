# KGE — MAGeCK with Interactive HTML Reports

Enhanced [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) (v0.5.9.5) with interactive Plotly-based HTML reports, pathway enrichment analysis, and FLUTE-style gene classification.

## Features

- **Interactive Plotly plots** — zoom, pan, hover, and export all QC and results plots
- **Dynamic gene tables** — sort, filter, search, and paginate across all genes (not just top 20)
- **GO/Hallmark pathway enrichment** — ORA via decoupler with interactive dot plots
- **FLUTE-style nine-square plot** — 4-group gene classification for MLE treatment analysis
- **sgRNA line plots** — click genes in tables to visualize sgRNA-level fold changes
- **Sample description section** — dataset overview with design matrix and detected files
- **Download section** — one-click CSV downloads for all gene summary files

## Quick Start

### Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda

### Install

```bash
git clone https://github.com/<your-username>/KGE.git
cd KGE
bash install.sh
```

This creates a conda environment `mageckenv` with all dependencies and installs the modified MAGeCK package.

To use a custom environment name:

```bash
bash install.sh my_custom_env
```

### Activate

```bash
conda activate mageckenv
```

### Usage

All standard MAGeCK commands work as before, with the added `--html-report` flag and the new `report` subcommand:

```bash
# Count sgRNAs from FASTQ files
mageck count -l library.txt -n experiment \
  --sample-label Treat,Control \
  --fastq treat.fastq control.fastq \
  --html-report

# RRA test
mageck test -k experiment.count.txt \
  -t Treat -c Control -n experiment \
  --html-report

# MLE analysis
mageck mle -k counts.csv -d design.txt \
  -n experiment --html-report

# Standalone report generation
mageck report -n experiment -k counts.csv
```

### Report Subcommand Options

```
mageck report -n NAME -k COUNT_TABLE [options]

Required:
  -n, --output-prefix    Output prefix (used to find gene_summary files)
  -k, --count-table      Read count table file

Optional:
  --gene-summary FILE [FILE ...]   Additional gene summary files
  --skip-enrichment                Skip GO/Hallmark enrichment analysis
  --enrichment-top-n N             Number of top genes for enrichment (default: 50)
  --enrichment-fdr F               FDR cutoff for enrichment (default: 0.05)
  --organism ORGANISM              Organism for pathway DB (default: human)
```

## What's Modified

These files contain modifications from the original MAGeCK v0.5.9.5:

| File | Changes |
|------|---------|
| `mageck/htmlReport.py` | New file — complete interactive HTML report generator |
| `mageck/argsParser.py` | Added `report` subcommand and enrichment arguments |
| `mageck/crisprFunction.py` | Added `mageck_report_main` import and hook |
| `mageck/mlemageck.py` | Fixed numpy array isinstance check for design matrix |
| `mageck/cnv_normalization.py` | Fixed numpy.str_ decode and np.float deprecation |
| `mageck/cli.py` | New entry point with report subcommand support |

## Environment

The conda environment includes:

- **Python 3.11**
- **MAGeCK 0.5.9.5** (base from bioconda)
- **Plotly 6.x** — interactive plots
- **decoupler 2.x** — pathway enrichment (ORA)
- **scanpy, anndata** — single-cell ecosystem
- **scikit-learn, scipy, statsmodels** — statistical analysis
- **matplotlib, seaborn** — static plotting (used by original MAGeCK)
- **PyDESeq2** — differential expression

## License

MAGeCK is distributed under the BSD License. See the original [MAGeCK repository](https://bitbucket.org/liulab/mageck/src/master/) for details. The HTML report modifications in this repository are provided under the same license.
