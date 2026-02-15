"""MAGeCK HTML Report Generator v2 — Interactive Plotly + Decoupler Enrichment

Generates a self-contained HTML report with interactive Plotly charts,
dynamic gene tables, pathway enrichment, and sgRNA line plots.

Only one external JS dependency: Plotly.js v2.27.0 (loaded from CDN).
Everything else is vanilla HTML/CSS/JavaScript.

Copyright (c) 2024 MAGeCK contributors
BSD License
"""

from __future__ import print_function
import os
import sys
import glob
import base64
import json
import logging
from datetime import datetime

import numpy as np
import pandas as pd

try:
    from sklearn.decomposition import PCA
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

try:
    import decoupler as dc
    HAS_DECOUPLER = True
except ImportError:
    HAS_DECOUPLER = False

from mageck.version import __version__

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_QUAL_COLORS = [
    '#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A',
    '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52',
]
_PLOTLY_CDN = 'https://cdn.plot.ly/plotly-2.27.0.min.js'
_PLOTLY_CONFIG = json.dumps({
    'displayModeBar': True, 'responsive': True,
    'displaylogo': False,
    'toImageButtonOptions': {'format': 'svg'},
})

# Color scales for enrichment dot plots
_CSCALE_DEPLETED = [[0, '#D1E5F0'], [0.25, '#92C5DE'], [0.5, '#4393C3'],
                     [0.75, '#2166AC'], [1, '#053061']]
_CSCALE_ENRICHED = [[0, '#FDDBC7'], [0.25, '#F4A582'], [0.5, '#D6604D'],
                     [0.75, '#B2182B'], [1, '#67001F']]

# Core essential genes for beta-score normalization (cell-cycle, ribosome,
# proteasome, spliceosome, DNA replication — from DepMap/Hart 2015).
# Works for both human (uppercase) and mouse (Yusa) libraries.
_CORE_ESSENTIAL_GENES = {
    'RPL3', 'RPL4', 'RPL5', 'RPL6', 'RPL7', 'RPL8', 'RPL9', 'RPL10', 'RPL11',
    'RPL13', 'RPL14', 'RPL15', 'RPL17', 'RPL18', 'RPL19', 'RPL21', 'RPL23',
    'RPL24', 'RPL26', 'RPL27', 'RPL30', 'RPL31', 'RPL32', 'RPL34', 'RPL35',
    'RPL36', 'RPL37', 'RPL38', 'RPS2', 'RPS3', 'RPS5', 'RPS6', 'RPS7',
    'RPS8', 'RPS9', 'RPS10', 'RPS11', 'RPS12', 'RPS13', 'RPS14', 'RPS15',
    'RPS16', 'RPS17', 'RPS18', 'RPS19', 'RPS20', 'RPS23', 'RPS24', 'RPS25',
    'PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7',
    'PSMB1', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7',
    'PSMC1', 'PSMC2', 'PSMC3', 'PSMC4', 'PSMC5', 'PSMC6',
    'PSMD1', 'PSMD2', 'PSMD3', 'PSMD4', 'PSMD6', 'PSMD7', 'PSMD11', 'PSMD12', 'PSMD14',
    'SF3B1', 'SF3B2', 'SF3B3', 'SF3A1', 'SF3A2', 'SF3A3',
    'SNRPD1', 'SNRPD2', 'SNRPD3', 'SNRPE', 'SNRPF',
    'CDK1', 'CDK2', 'CDK7', 'CDK9', 'CCNA2', 'CCNB1', 'CCND1', 'CCNE1', 'CCNH',
    'CDC5L', 'CDC16', 'CDC20', 'CDC23', 'CDC27',
    'PCNA', 'MCM2', 'MCM3', 'MCM4', 'MCM5', 'MCM6', 'MCM7',
    'RFC1', 'RFC2', 'RFC3', 'RFC4', 'RFC5',
    'POLR2A', 'POLR2B', 'POLR2C', 'POLR2D', 'POLR2E', 'POLR2F', 'POLR2G', 'POLR2H',
    'EEF2', 'EIF3A', 'EIF3B', 'EIF3C', 'EIF3D',
    'NUP93', 'NUP107', 'NUP133', 'NUP155', 'NUP160', 'SUPT5H', 'SUPT6H',
    'UBA1', 'UBB', 'RAN', 'KPNB1', 'XPO1', 'NHP2', 'DKC1', 'GAR1', 'NOP10',
}

# Plot descriptions (professional bioinformatics captions)
_DESC = {
    'lib_size': ('Total read counts per sample. Uniform library sizes across '
                 'conditions indicate comparable sequencing depth, which is '
                 'essential for reliable differential analysis.'),
    'count_dist': ('Distribution of median-normalized read counts across all samples. '
                   'Each box shows the interquartile range; whiskers extend to 1.5\u00d7IQR. '
                   'Comparable distributions indicate effective normalization.'),
    'density_norm': ('Distribution of median-normalized read counts '
                     '(log\u2082-transformed). The Y-axis shows the total frequency '
                     '(number of sgRNAs). A shift toward lower values in treated '
                     'samples may indicate library-wide depletion or selection pressure.'),
    'zero_count': ('sgRNAs with zero reads per sample. Elevated zero-counts '
                   'in treatment conditions suggest guide dropout under '
                   'selection, while high baseline zeros may reflect '
                   'insufficient library coverage.'),
    'gini': ('Gini index quantifies guide representation evenness: 0 = '
             'perfectly uniform, 1 = maximally skewed. Typical CRISPR '
             'library screens range 0.1\u20130.3; higher values may indicate '
             'strong selection pressure or sequencing bias.'),
    'pca': ('Principal Component Analysis of normalized sgRNA counts. Samples '
            'should cluster by biological condition along the leading '
            'components. Outlier samples may warrant quality review.'),
    'corr': ('Pearson correlation of log\u2082-normalized sgRNA counts. '
             'Replicates typically show r > 0.9; lower correlation between '
             'conditions reflects differential selection.'),
    'rra_volcano': ('Volcano plot of gene-level log\u2082 fold-change '
                    '(x-axis) vs. \u2212log\u2081\u2080(FDR) (y-axis) from '
                    'the Robust Rank Aggregation (RRA) algorithm. Blue = '
                    'depleted, Red = enriched. Use the controls above to '
                    'adjust FDR and LFC cutoffs.'),
    'mle_volcano': ('Volcano plot of MLE beta scores (x-axis) vs. '
                    '\u2212log\u2081\u2080(FDR) (y-axis). '
                    'Negative beta = depleted (dropout); '
                    'positive beta = enriched. Use the controls '
                    'above to adjust FDR and Beta score cutoffs.'),
    'mle_scatter': ('Beta score comparison between two MLE conditions. Genes '
                    'along the diagonal show concordant effects; off-diagonal '
                    'genes exhibit condition-specific phenotypes. Top 10 '
                    'depleted and enriched genes are labeled.'),
    'enrich_depl': ('MSigDB Hallmark pathway enrichment (ORA) of the top '
                    'depleted genes. Dot size reflects statistical '
                    'significance (\u2212log\u2081\u2080 adj. p-value); '
                    'color encodes ORA enrichment score. '
                    'These pathways are under negative selection.'),
    'enrich_enr': ('MSigDB Hallmark pathway enrichment (ORA) of the top '
                   'enriched genes. Dot size reflects statistical '
                   'significance (\u2212log\u2081\u2080 adj. p-value); '
                   'color encodes ORA enrichment score. '
                   'These pathways are under positive selection.'),
    'enrich_sens': ('MSigDB Hallmark pathway enrichment (ORA) of top '
                    'depleted genes (negative beta). Dot size = '
                    '\u2212log\u2081\u2080(adj. p-value); color = ORA score. '
                    'These pathways are under negative selection.'),
    'enrich_res': ('MSigDB Hallmark pathway enrichment (ORA) of top '
                   'enriched genes (positive beta). Dot size = '
                   '\u2212log\u2081\u2080(adj. p-value); color = ORA score. '
                   'These pathways are under positive selection.'),
    'sgrna_line': ('sgRNA-level log\u2082 fold-change relative to the '
                   'sample mean, centered at 0. By default, the top 5 '
                   'enriched and top 5 depleted genes are displayed. '
                   'Click gene rows in the table to add or remove traces. '
                   'Each line represents one sgRNA.'),
    'gene_table': ('Full gene-level summary table. Sort by any column '
                   '(click header), filter by FDR threshold, or search '
                   'by gene name. Click a row to display its sgRNA counts '
                   'in the line plot below.'),
    'mapping_stats': ('Mapping statistics from MAGeCK count. Total reads, '
                      'mapped reads, and mapping rate per sample. Samples '
                      'with low mapping rates may indicate library mismatch '
                      'or sequencing issues.'),
    'mle_beta_density': ('Distribution of gene-level beta scores for each MLE '
                         'condition. Comparing density shapes reveals whether '
                         'treatment induces broader selection effects relative '
                         'to control. A wider distribution in the treatment '
                         'condition indicates stronger selective pressure.'),
    'mle_consistency': ('Scatterplot of gene beta scores between treatment and '
                        'control conditions (ConsistencyView). The red line '
                        'shows the linear regression fit; the grey dashed '
                        'line represents y = x (perfect concordance). Genes '
                        'along the diagonal show consistent effects across '
                        'conditions.'),
    'mle_selection_scatter': ('Scatterplot of beta scores for treatment vs. control. '
                              'Use the diagonal cutoff slider to classify genes: '
                              'pink dots = enriched (+), genes with increased beta '
                              'after treatment; blue dots = depleted (\u2212), genes '
                              'with decreased beta after treatment; grey dots = '
                              'genes within the cutoff range. Diagonal lines show '
                              'the \u00b1cutoff boundary relative to y\u2009=\u2009x.'),
    'mle_nine_square': ('Treatment-associated gene classification (nine-square). '
                        '<b>Group 1 (green)</b>: strongly negatively selected in control '
                        'but weakly selected in treatment \u2014 potential drug pathway targets. '
                        '<b>Group 2 (orange)</b>: weakly selected in control but strongly '
                        'positively selected in treatment \u2014 possible drug resistance genes. '
                        '<b>Group 3 (blue)</b>: strongly positively selected in control '
                        'but weakly selected in treatment. '
                        '<b>Group 4 (purple)</b>: weakly selected in control but strongly '
                        'negatively selected in treatment \u2014 possibly synthetically lethal '
                        'with treatment. Adjust the beta cutoff to change classification '
                        'stringency.'),
}


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _json_safe(obj):
    """Recursively convert numpy types and NaN/Inf to JSON-safe values."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _json_safe(obj.tolist())
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return None if (np.isnan(v) or np.isinf(v)) else v
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, float):
        return None if (np.isnan(obj) or np.isinf(obj)) else obj
    return obj


def _jdumps(obj):
    """JSON-serialize with numpy and NaN handling, compact output."""
    return json.dumps(_json_safe(obj), separators=(',', ':'))


def _calculate_gini(x):
    """Compute Gini index — exact replica of MAGeCK's mageckcount_gini().

    MAGeCK computes Gini on log(count+1) transformed values, so callers
    should pass log-transformed values to match MAGeCK output.
    Formula: 1.0 - 2.0*(n - gssum/ysum)/(n-1)
    Reference: http://en.wikipedia.org/wiki/Gini_coefficient
    """
    arr = np.array(x, dtype=float)
    arr = arr[~np.isnan(arr)]
    if len(arr) == 0:
        return 0.0
    arr = np.sort(arr)
    n = len(arr)
    total = arr.sum()
    if total == 0.0:
        total = 1.0
    index = np.arange(1.0, n + 1.0)
    gssum = np.sum(index * arr)
    return 1.0 - 2.0 * (n - gssum / total) / (n - 1)


def _median_normalize(count_df):
    """Median-ratio normalization (DESeq2-style) with total-count fallback.

    Only rows where ALL samples have non-zero counts are used to compute
    size factors.  If the result is degenerate (all factors ~1 despite
    different library sizes), falls back to total-count normalization.
    """
    sample_cols = [c for c in count_df.columns if c not in ('sgRNA', 'Gene')]
    data = count_df[sample_cols].values.astype(float)
    totals = data.sum(axis=0)
    mean_total = totals.mean()

    # Try median-ratio on rows where all samples are non-zero
    nonzero_mask = np.all(data > 0, axis=1)
    use_total = True
    if nonzero_mask.sum() >= 100:
        sub = data[nonzero_mask]
        log_sub = np.log(sub)
        geo_mean = np.exp(log_sub.mean(axis=1))
        ratios = sub / geo_mean[:, None]
        size_factors = np.median(ratios, axis=0)
        # Check for degenerate results (all ~1 when library sizes differ)
        sf_range = size_factors.max() - size_factors.min()
        total_cv = np.std(totals) / mean_total if mean_total > 0 else 0
        if sf_range > 0.01 or total_cv < 0.01:
            use_total = False

    if use_total:
        # Total-count normalization: scale to mean library size
        size_factors = totals / mean_total

    size_factors[size_factors == 0] = 1.0
    normalized = data / size_factors[None, :]
    result = count_df.copy()
    result[sample_cols] = normalized
    return result


def _desc_html(key):
    """Return a plot description paragraph for the given key."""
    text = _DESC.get(key, '')
    if not text:
        return ''
    return '<p class="plot-desc">' + text + '</p>\n'


def _trunc(name, maxlen=42):
    """Truncate a string, adding ellipsis if needed."""
    if len(name) <= maxlen:
        return name
    return name[:maxlen - 1] + '\u2026'


# ---------------------------------------------------------------------------
# Main report class
# ---------------------------------------------------------------------------

class MAGeCKHTMLReport:
    """Interactive HTML report generator for MAGeCK CRISPR screen analysis."""

    def __init__(self, output_prefix, count_table=None, title=None,
                 gene_summary_files=None, sgrna_summary_files=None,
                 design_matrix=None, fdr_threshold=0.05, top_n=20,
                 norm_method='median', skip_qc=False, skip_results=False,
                 enrichment_top_n=50, organism='human',
                 skip_enrichment=False, enrichment_fdr=0.05):
        self.output_prefix = output_prefix
        self.count_table_path = count_table
        self.title = title or 'CRISPR Screen Report'
        self.gene_summary_files = gene_summary_files or []
        self.sgrna_summary_files = sgrna_summary_files or []
        self.design_matrix_path = design_matrix
        self.fdr_threshold = fdr_threshold
        self.top_n = top_n
        self.norm_method = norm_method
        self.skip_qc = skip_qc
        self.skip_results = skip_results
        self.enrichment_top_n = enrichment_top_n
        self.organism = organism
        self.skip_enrichment = skip_enrichment
        self.enrichment_fdr = enrichment_fdr

        # Data containers
        self.count_df = None
        self.normalized_df = None
        self.sample_names = []
        self.rra_gene_summaries = {}
        self.rra_sgrna_summaries = {}
        self.mle_gene_summary = None
        self.mle_sgrna_summary = None
        self.mle_conditions = []
        self.design_matrix_df = None
        self.countsummary_df = None
        self._hallmark_net = None

        # Enrichment results for download
        self._enrichment_counter = 0
        self._enrichment_store = {}  # key -> {source, score, padj} list

        # Default gene selections for line plots
        self._default_genes = {}  # lp_id -> [gene_names]

        # Downloadable files (detected in detect_outputs)
        self._log_files = []  # list of (label, path) tuples

        # HTML assembly
        self.sections = []
        self._div_counter = 0
        self._table_ids = []

    # ------------------------------------------------------------------
    # File detection
    # ------------------------------------------------------------------

    def detect_outputs(self):
        """Auto-detect MAGeCK output files in the working directory."""
        base_dir = os.path.dirname(self.output_prefix) or '.'
        prefix_base = os.path.basename(self.output_prefix)

        if self.count_table_path and os.path.isfile(self.count_table_path):
            self._load_count_table(self.count_table_path)
        else:
            for candidate in [self.output_prefix + '.count.txt',
                              self.output_prefix + '.count_normalized.txt']:
                if os.path.isfile(candidate):
                    self._load_count_table(candidate)
                    break

        # Try to load MAGeCK's pre-computed normalized count table.
        # This is more accurate than re-computing normalization ourselves.
        self._load_normalized_count_file(base_dir)

        if self.gene_summary_files:
            gs_files = self.gene_summary_files
        else:
            gs_files = sorted(glob.glob(os.path.join(base_dir, '*.gene_summary.txt')))

        for gsf in gs_files:
            if not os.path.isfile(gsf):
                continue
            if self._is_mle_gene_summary(gsf):
                self.mle_gene_summary = pd.read_csv(gsf, sep='\t')
                self.mle_conditions = self._extract_mle_conditions(self.mle_gene_summary)
                logging.info('Detected MLE gene summary: %s (conditions: %s)',
                             gsf, ', '.join(self.mle_conditions))
            else:
                label = self._extract_comparison_label(gsf)
                self.rra_gene_summaries[label] = pd.read_csv(gsf, sep='\t')
                logging.info('Detected RRA gene summary: %s (%s)', gsf, label)

        if self.sgrna_summary_files:
            sg_files = self.sgrna_summary_files
        else:
            sg_files = sorted(glob.glob(os.path.join(base_dir, '*.sgrna_summary.txt')))

        for sgf in sg_files:
            if not os.path.isfile(sgf):
                continue
            header = open(sgf).readline().strip()
            if 'control_count' in header:
                label = self._extract_comparison_label(sgf)
                self.rra_sgrna_summaries[label] = pd.read_csv(sgf, sep='\t')
            elif self.mle_sgrna_summary is None:
                self.mle_sgrna_summary = pd.read_csv(sgf, sep='\t')

        if self.design_matrix_path and os.path.isfile(self.design_matrix_path):
            self.design_matrix_df = pd.read_csv(self.design_matrix_path, sep='\t')
        else:
            # Search for common design matrix naming patterns
            for dm_name in ['design_matrix.txt', 'designmatrix.txt']:
                dm_candidate = os.path.join(base_dir, dm_name)
                if os.path.isfile(dm_candidate):
                    self.design_matrix_df = pd.read_csv(dm_candidate, sep='\t')
                    logging.info('Loaded design matrix: %s', dm_candidate)
                    break

        # Detect countsummary from mageck count (mapping statistics)
        cs_candidate = self.output_prefix + '.countsummary.txt'
        if not os.path.isfile(cs_candidate):
            cs_files = sorted(glob.glob(os.path.join(base_dir, '*.countsummary.txt')))
            if cs_files:
                cs_candidate = cs_files[0]
        if os.path.isfile(cs_candidate):
            self.countsummary_df = pd.read_csv(cs_candidate, sep='\t')
            logging.info('Loaded count summary: %s', cs_candidate)

        # Detect MAGeCK log files
        log_candidate = self.output_prefix + '.log'
        if os.path.isfile(log_candidate):
            self._log_files.append((os.path.basename(log_candidate), log_candidate))
        else:
            log_files = sorted(glob.glob(os.path.join(base_dir, '*.log')))
            for lf in log_files:
                self._log_files.append((os.path.basename(lf), lf))
        if self._log_files:
            logging.info('Detected %d log file(s)', len(self._log_files))

    def _load_count_table(self, path):
        sep = ',' if path.endswith('.csv') else '\t'
        self.count_df = pd.read_csv(path, sep=sep)
        cols = list(self.count_df.columns)
        if cols[0].lower() in ('sgrna', 'shrna', 'guide'):
            cols[0] = 'sgRNA'
        if len(cols) > 1 and cols[1].lower() == 'gene':
            cols[1] = 'Gene'
        self.count_df.columns = cols
        self.sample_names = [c for c in cols if c not in ('sgRNA', 'Gene')]
        self.count_table_path = path
        logging.info('Loaded count table: %s (%d sgRNAs, %d samples)',
                     path, len(self.count_df), len(self.sample_names))

    def _load_normalized_count_file(self, base_dir):
        """Load MAGeCK's pre-computed normalized count file if available."""
        candidates = []
        # Derive from count_table_path: foo.count.txt -> foo.count_normalized.txt
        if self.count_table_path:
            ct = self.count_table_path
            if ct.endswith('.count.txt'):
                candidates.append(ct.replace('.count.txt', '.count_normalized.txt'))
        # Also check output_prefix variants
        candidates.append(self.output_prefix + '.count_normalized.txt')
        # Search directory for any normalized file
        candidates.extend(sorted(glob.glob(os.path.join(base_dir, '*.count_normalized.txt'))))
        for cand in candidates:
            if os.path.isfile(cand):
                try:
                    sep = ',' if cand.endswith('.csv') else '\t'
                    ndf = pd.read_csv(cand, sep=sep)
                    ncols = list(ndf.columns)
                    if ncols[0].lower() in ('sgrna', 'shrna', 'guide'):
                        ncols[0] = 'sgRNA'
                    if len(ncols) > 1 and ncols[1].lower() == 'gene':
                        ncols[1] = 'Gene'
                    ndf.columns = ncols
                    self.normalized_df = ndf
                    logging.info('Loaded pre-computed normalized counts: %s', cand)
                    return
                except Exception as e:
                    logging.warning('Failed to load normalized file %s: %s', cand, e)

    @staticmethod
    def _is_mle_gene_summary(filepath):
        with open(filepath) as fh:
            header = fh.readline()
        return '|beta' in header

    @staticmethod
    def _extract_mle_conditions(df):
        conditions = []
        for col in df.columns:
            if col.endswith('|beta'):
                conditions.append(col.replace('|beta', ''))
        return conditions

    @staticmethod
    def _extract_comparison_label(filepath):
        base = os.path.basename(filepath)
        return base.replace('.gene_summary.txt', '').replace('.sgrna_summary.txt', '')

    def _condition_to_sample(self, cond):
        """Map MLE condition name to sample name using design matrix."""
        if self.design_matrix_df is not None:
            sample_col = self.design_matrix_df.columns[0]
            if cond in self.design_matrix_df.columns:
                mask = self.design_matrix_df[cond] == 1
                if mask.any():
                    return str(self.design_matrix_df.loc[mask, sample_col].values[0])
        return cond

    # ------------------------------------------------------------------
    # Normalization
    # ------------------------------------------------------------------

    def normalize_counts(self):
        if self.count_df is None:
            return
        # If we already loaded MAGeCK's pre-computed normalized file, use it
        if self.normalized_df is not None:
            logging.info('Using pre-computed normalized counts (from file)')
            return
        if self.norm_method == 'median':
            self.normalized_df = _median_normalize(self.count_df)
        elif self.norm_method == 'total':
            sample_cols = self.sample_names
            totals = self.count_df[sample_cols].sum()
            mean_total = totals.mean()
            self.normalized_df = self.count_df.copy()
            for s in sample_cols:
                self.normalized_df[s] = self.count_df[s] * mean_total / totals[s]
        else:
            self.normalized_df = self.count_df.copy()

    def normalize_betas(self):
        """Normalize MLE beta scores using core essential genes (FLUTE-style).

        For each condition, divides all betas by the absolute median beta
        of core essential genes.  This rescales so that the typical
        essential-gene depletion is -1, making conditions comparable.
        """
        if self.mle_gene_summary is None or not self.mle_conditions:
            return
        df = self.mle_gene_summary
        gene_upper = df['Gene'].str.upper()
        present = gene_upper.isin(_CORE_ESSENTIAL_GENES)
        n_present = present.sum()
        if n_present < 10:
            logging.warning('Only %d core essential genes found — skipping '
                            'beta normalization.', n_present)
            return
        for cond in self.mle_conditions:
            beta_col = cond + '|beta'
            if beta_col not in df.columns:
                continue
            cc_betas = df.loc[present, beta_col].values.astype(float)
            cc_median = float(np.median(cc_betas[np.isfinite(cc_betas)]))
            if abs(cc_median) < 1e-6:
                continue
            df[beta_col] = df[beta_col] / abs(cc_median)
        logging.info('Beta scores normalized using %d essential genes '
                     '(FLUTE cell-cycle style)', n_present)

    def _get_baseline_name(self):
        """Return the baseline sample name from the design matrix."""
        if self.design_matrix_df is not None:
            cols = list(self.design_matrix_df.columns)
            sample_col = cols[0]
            if 'baseline' in cols:
                mask = self.design_matrix_df['baseline'] == 1
                if mask.any():
                    return str(self.design_matrix_df.loc[mask, sample_col].values[0])
        return 'baseline'

    # ------------------------------------------------------------------
    # Plotly helpers
    # ------------------------------------------------------------------

    def _next_id(self, prefix='plt'):
        self._div_counter += 1
        return '%s-%d' % (prefix, self._div_counter)

    def _plotly_div(self, traces, layout, div_id=None):
        """Generate an HTML div + script for a Plotly chart."""
        div_id = div_id or self._next_id()
        t = _jdumps(traces)
        l = _jdumps(layout)
        return ('<div id="' + div_id + '" class="plotly-chart"></div>\n'
                '<script>Plotly.newPlot("' + div_id + '",' + t + ',' + l + ',' + _PLOTLY_CONFIG + ');</script>\n')

    # ------------------------------------------------------------------
    # Section 1: QC Plots (compact, automargin)
    # ------------------------------------------------------------------

    def plot_mapping_stats(self):
        """Bar plot + table of mapping statistics from mageck count."""
        cs = self.countsummary_df
        if cs is None:
            return None, None
        labels = cs['Label'].values.astype(str) if 'Label' in cs.columns else cs.iloc[:, 1].values.astype(str)
        total = cs['Reads'].values.astype(float) if 'Reads' in cs.columns else None
        mapped = cs['Mapped'].values.astype(float) if 'Mapped' in cs.columns else None
        if total is None or mapped is None:
            return None, None
        pct = cs['Percentage'].values.astype(float) if 'Percentage' in cs.columns else mapped / np.maximum(total, 1)
        unmapped = total - mapped
        # Bar plot: stacked mapped + unmapped
        traces = [
            {'type': 'bar', 'x': labels.tolist(), 'y': mapped.tolist(),
             'name': 'Mapped', 'marker': {'color': '#2980b9'},
             'hovertemplate': 'Sample: %{x}<br>Mapped: %{y:,.0f}<extra></extra>'},
            {'type': 'bar', 'x': labels.tolist(), 'y': unmapped.tolist(),
             'name': 'Unmapped', 'marker': {'color': '#e74c3c'},
             'hovertemplate': 'Sample: %{x}<br>Unmapped: %{y:,.0f}<extra></extra>'},
        ]
        layout = {
            'title': {'text': 'Mapping Statistics', 'font': {'size': 13}},
            'barmode': 'stack',
            'xaxis': {'title': 'Sample', 'automargin': True, 'tickfont': {'size': 10}},
            'yaxis': {'title': 'Read Count', 'automargin': True, 'tickfont': {'size': 10}},
            'legend': {'font': {'size': 10}},
            'height': 350, 'margin': {'t': 45, 'b': 45, 'l': 65, 'r': 15}}
        bar_html = self._plotly_div(traces, layout)
        # Table
        tbl = '<table class="simple-table"><thead><tr>'
        tbl += '<th>Sample</th><th>Total Reads</th><th>Mapped Reads</th><th>Mapping Rate</th>'
        tbl += '</tr></thead><tbody>'
        for i in range(len(labels)):
            rate_str = '{:.1%}'.format(pct[i]) if pct[i] <= 1 else '{:.1%}'.format(pct[i] / 100)
            tbl += '<tr><td>' + labels[i] + '</td>'
            tbl += '<td>' + '{:,.0f}'.format(total[i]) + '</td>'
            tbl += '<td>' + '{:,.0f}'.format(mapped[i]) + '</td>'
            tbl += '<td>' + rate_str + '</td></tr>'
        tbl += '</tbody></table>'
        return bar_html, tbl

    def plot_library_size(self):
        totals = self.count_df[self.sample_names].sum()
        colors = _QUAL_COLORS[:len(self.sample_names)]
        traces = [{'type': 'bar', 'x': self.sample_names,
                   'y': totals.values.tolist(),
                   'text': ['%s' % '{:,.0f}'.format(v) for v in totals.values],
                   'textposition': 'outside', 'cliponaxis': False,
                   'textfont': {'size': 10},
                   'marker': {'color': colors}}]
        layout = {'title': {'text': 'Library Size (Total Reads)', 'font': {'size': 13}},
                  'xaxis': {'automargin': True, 'tickfont': {'size': 10}},
                  'yaxis': {'title': 'Total Read Counts', 'automargin': True,
                            'tickfont': {'size': 10}},
                  'height': 340, 'margin': {'t': 45, 'b': 40, 'l': 65, 'r': 15}}
        return self._plotly_div(traces, layout)

    def plot_count_distribution(self):
        """Box plot + density plot of median-normalized read counts."""
        df_norm = self.normalized_df if self.normalized_df is not None else self.count_df
        colors = _QUAL_COLORS[:len(self.sample_names)]

        # Left (1st): Box plot with outlier dots (MAGeCK style)
        traces_box = []
        for i, s in enumerate(self.sample_names):
            vals = np.log2(df_norm[s].values.astype(float) + 1)
            vals = vals[np.isfinite(vals)]
            traces_box.append({
                'type': 'box', 'y': np.round(vals, 2).tolist(),
                'name': s,
                'marker': {'color': colors[i], 'size': 4, 'opacity': 0.6},
                'boxpoints': 'outliers',
                'jitter': 0.3,
                'pointpos': 0})
        layout_box = {'title': {'text': 'Normalized Read Count Distribution', 'font': {'size': 13}},
                      'yaxis': {'title': 'log\u2082(read counts)', 'automargin': True,
                                'tickfont': {'size': 10}},
                      'xaxis': {'automargin': True, 'tickfont': {'size': 10}},
                      'legend': {'font': {'size': 10}},
                      'height': 340, 'margin': {'t': 45, 'b': 40, 'l': 50, 'r': 15}}
        box_html = self._plotly_div(traces_box, layout_box)

        # Right (2nd): Histogram of read count distribution (total frequency)
        traces_dens = []
        for i, s in enumerate(self.sample_names):
            vals = np.log2(df_norm[s].values.astype(float) + 1)
            vals = vals[np.isfinite(vals)]
            counts_d, edges_d = np.histogram(vals, bins=80, density=False)
            centers_d = ((edges_d[:-1] + edges_d[1:]) / 2).tolist()
            traces_dens.append({
                'type': 'scatter', 'mode': 'lines', 'fill': 'tozeroy',
                'x': centers_d, 'y': counts_d.tolist(),
                'name': s, 'opacity': 0.5,
                'line': {'color': colors[i], 'width': 2}})
        layout_dens = {'title': {'text': 'Distribution of Read Counts', 'font': {'size': 13}},
                       'xaxis': {'title': 'log\u2082(read counts)', 'automargin': True,
                                 'tickfont': {'size': 9},
                                 'dtick': 4, 'tick0': 0},
                       'yaxis': {'title': 'Frequency', 'automargin': True,
                                 'tickfont': {'size': 10},
                                 'rangemode': 'tozero'},
                       'legend': {'font': {'size': 10}},
                       'height': 340, 'margin': {'t': 45, 'b': 45, 'l': 60, 'r': 15}}
        dens_html = self._plotly_div(traces_dens, layout_dens)
        return box_html, dens_html

    def plot_zero_count_guides(self):
        zeros = [(self.count_df[s] == 0).sum() for s in self.sample_names]
        total = len(self.count_df)
        pcts = [z / total * 100 for z in zeros]
        colors = _QUAL_COLORS[:len(self.sample_names)]
        traces = [{'type': 'bar', 'x': self.sample_names, 'y': zeros,
                   'text': ['{:,}'.format(z) for z in zeros],
                   'textposition': 'outside', 'cliponaxis': False,
                   'textfont': {'size': 10},
                   'marker': {'color': colors}}]
        layout = {'title': {'text': 'Zero-Count sgRNAs', 'font': {'size': 13}},
                  'xaxis': {'automargin': True, 'tickfont': {'size': 10}},
                  'yaxis': {'title': 'Number of Zero-Count sgRNAs', 'automargin': True,
                            'tickfont': {'size': 10}},
                  'height': 340, 'margin': {'t': 45, 'b': 40, 'l': 65, 'r': 15}}
        return self._plotly_div(traces, layout)

    def plot_gini_index(self):
        # Use countsummary GiniIndex if available (most accurate — from mageck count)
        ginis = None
        if self.countsummary_df is not None and 'GiniIndex' in self.countsummary_df.columns:
            cs = self.countsummary_df
            label_col = 'Label' if 'Label' in cs.columns else cs.columns[1]
            cs_map = dict(zip(cs[label_col].astype(str), cs['GiniIndex'].astype(float)))
            if all(s in cs_map for s in self.sample_names):
                ginis = [cs_map[s] for s in self.sample_names]
                logging.info('Using Gini index from countsummary file')
        if ginis is None:
            # Compute matching MAGeCK: Gini on log(count + 1) values
            ginis = [_calculate_gini(np.log(self.count_df[s].values.astype(float) + 1.0))
                     for s in self.sample_names]
        colors = _QUAL_COLORS[:len(self.sample_names)]
        traces = [{'type': 'bar', 'x': self.sample_names, 'y': ginis,
                   'text': ['{:.3f}'.format(g) for g in ginis],
                   'textposition': 'outside', 'cliponaxis': False,
                   'textfont': {'size': 10},
                   'marker': {'color': colors}}]
        layout = {'title': {'text': 'Gini Index', 'font': {'size': 13}},
                  'xaxis': {'automargin': True, 'tickfont': {'size': 10}},
                  'yaxis': {'title': 'Gini Index', 'automargin': True,
                            'range': [0, max(ginis) * 1.3 if max(ginis) > 0 else 1],
                            'tickfont': {'size': 10}},
                  'height': 340, 'margin': {'t': 45, 'b': 40, 'l': 60, 'r': 15}}
        return self._plotly_div(traces, layout)

    def plot_pca(self):
        df = self.normalized_df if self.normalized_df is not None else self.count_df
        data = np.log2(df[self.sample_names].values.astype(float) + 1).T
        if HAS_SKLEARN:
            pca = PCA(n_components=min(len(self.sample_names), 3))
            coords = pca.fit_transform(data)
            var_exp = pca.explained_variance_ratio_ * 100
        else:
            centered = data - data.mean(axis=0)
            U, S, Vt = np.linalg.svd(centered, full_matrices=False)
            coords = U[:, :3] * S[:3]
            total_var = (S ** 2).sum()
            var_exp = (S[:3] ** 2) / total_var * 100
        colors = _QUAL_COLORS[:len(self.sample_names)]
        traces = [{'type': 'scatter', 'mode': 'markers+text',
                   'x': coords[:, 0].tolist(), 'y': coords[:, 1].tolist(),
                   'text': self.sample_names, 'textposition': 'top right',
                   'textfont': {'size': 10},
                   'marker': {'color': colors, 'size': 12,
                              'line': {'width': 1, 'color': '#333'}},
                   'hovertemplate': '%{text}<br>PC1: %{x:.2f}<br>PC2: %{y:.2f}'}]
        layout = {'title': {'text': 'PCA', 'font': {'size': 13}},
                  'xaxis': {'title': 'PC1 ({:.1f}% var)'.format(var_exp[0]),
                            'automargin': True, 'tickfont': {'size': 10}},
                  'yaxis': {'title': 'PC2 ({:.1f}% var)'.format(var_exp[1]),
                            'automargin': True, 'tickfont': {'size': 10}},
                  'height': 340, 'margin': {'t': 45, 'b': 45, 'l': 60, 'r': 15},
                  'showlegend': False}
        return self._plotly_div(traces, layout)

    def plot_sample_correlation(self):
        df = self.normalized_df if self.normalized_df is not None else self.count_df
        log_data = np.log2(df[self.sample_names].values.astype(float) + 1)
        corr = np.corrcoef(log_data.T)

        # Hierarchical clustering to reorder rows/columns
        sample_order = list(range(len(self.sample_names)))
        Z = None
        try:
            from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram as scipy_dendrogram
            from scipy.spatial.distance import squareform
            dist = 1.0 - corr
            np.fill_diagonal(dist, 0)
            dist = np.clip(dist, 0, None)
            condensed = squareform(dist)
            Z = linkage(condensed, method='average')
            sample_order = leaves_list(Z).tolist()
        except ImportError:
            pass  # scipy not available — use original order

        ordered_names = [self.sample_names[i] for i in sample_order]
        corr_ordered = corr[np.ix_(sample_order, sample_order)]

        text_vals = [['{:.3f}'.format(corr_ordered[i][j])
                      for j in range(len(ordered_names))]
                     for i in range(len(ordered_names))]

        n = len(self.sample_names)

        # Build dendrograms if scipy is available
        if Z is not None:
            # Get dendrogram coordinates (no_plot=True to avoid matplotlib)
            dendro = scipy_dendrogram(Z, no_plot=True, labels=ordered_names)

            # Top dendrogram traces
            dendro_top_traces = []
            for i in range(len(dendro['icoord'])):
                dendro_top_traces.append({
                    'type': 'scatter', 'mode': 'lines',
                    'x': dendro['icoord'][i],
                    'y': dendro['dcoord'][i],
                    'line': {'color': '#2c3e50', 'width': 1.2},
                    'hoverinfo': 'skip', 'showlegend': False,
                    'xaxis': 'x2', 'yaxis': 'y2'})

            # Left dendrogram traces (rotated: swap x/y, flip y)
            dendro_left_traces = []
            for i in range(len(dendro['icoord'])):
                dendro_left_traces.append({
                    'type': 'scatter', 'mode': 'lines',
                    'x': dendro['dcoord'][i],
                    'y': dendro['icoord'][i],
                    'line': {'color': '#2c3e50', 'width': 1.2},
                    'hoverinfo': 'skip', 'showlegend': False,
                    'xaxis': 'x3', 'yaxis': 'y3'})

            # Heatmap trace on its own axes
            heatmap_trace = {
                'type': 'heatmap',
                'z': corr_ordered.tolist(),
                'x': list(range(n)), 'y': list(range(n)),
                'text': text_vals, 'texttemplate': '%{text}',
                'textfont': {'size': 11},
                'colorscale': 'RdYlBu', 'reversescale': True,
                'zmin': 0.9, 'zmax': 1.0,
                'hovertemplate': '%{customdata[0]} vs %{customdata[1]}<br>r = %{z:.4f}',
                'customdata': [[[ordered_names[i], ordered_names[j]]
                               for j in range(n)] for i in range(n)],
                'xaxis': 'x', 'yaxis': 'y'}

            traces = [heatmap_trace] + dendro_top_traces + dendro_left_traces

            # Dendrogram tick positions: 5, 15, 25, ... (default scipy spacing)
            tick_vals = [5 + 10 * i for i in range(n)]
            hm_tick_vals = list(range(n))

            plot_h = max(420, n * 70 + 150)
            dendro_frac = 0.15  # fraction of plot for dendrograms

            layout = {
                'title': {'text': 'Sample Correlation (Pearson, clustered)', 'font': {'size': 13}},
                'height': plot_h,
                'margin': {'t': 45, 'b': 50, 'l': 15, 'r': 15},
                'showlegend': False,
                # Main heatmap axes
                'xaxis': {'domain': [dendro_frac, 1.0],
                          'tickvals': hm_tick_vals, 'ticktext': ordered_names,
                          'automargin': True, 'tickfont': {'size': 10},
                          'side': 'bottom'},
                'yaxis': {'domain': [0, 1.0 - dendro_frac],
                          'autorange': 'reversed',
                          'tickvals': hm_tick_vals, 'ticktext': ordered_names,
                          'automargin': True, 'tickfont': {'size': 10}},
                # Top dendrogram
                'xaxis2': {'domain': [dendro_frac, 1.0],
                           'showgrid': False, 'zeroline': False,
                           'showticklabels': False,
                           'range': [min(d for ic in dendro['icoord'] for d in ic) - 2,
                                     max(d for ic in dendro['icoord'] for d in ic) + 2],
                           'anchor': 'y2'},
                'yaxis2': {'domain': [1.0 - dendro_frac, 1.0],
                           'showgrid': False, 'zeroline': False,
                           'showticklabels': False, 'anchor': 'x2'},
                # Left dendrogram
                'xaxis3': {'domain': [0, dendro_frac - 0.02],
                           'showgrid': False, 'zeroline': False,
                           'showticklabels': False,
                           'autorange': 'reversed', 'anchor': 'y3'},
                'yaxis3': {'domain': [0, 1.0 - dendro_frac],
                           'showgrid': False, 'zeroline': False,
                           'showticklabels': False,
                           'range': [min(d for ic in dendro['icoord'] for d in ic) - 2,
                                     max(d for ic in dendro['icoord'] for d in ic) + 2],
                           'anchor': 'x3'},
            }
        else:
            # Fallback: simple heatmap without dendrograms
            traces = [{'type': 'heatmap',
                       'z': corr_ordered.tolist(),
                       'x': ordered_names, 'y': ordered_names,
                       'text': text_vals, 'texttemplate': '%{text}',
                       'textfont': {'size': 11},
                       'colorscale': 'RdYlBu', 'reversescale': True,
                       'zmin': 0.9, 'zmax': 1.0,
                       'hovertemplate': '%{y} vs %{x}<br>r = %{z:.4f}'}]
            layout = {'title': {'text': 'Sample Correlation (Pearson, clustered)', 'font': {'size': 13}},
                      'height': max(320, n * 70 + 80),
                      'margin': {'t': 45, 'b': 50, 'l': 90, 'r': 15},
                      'xaxis': {'automargin': True, 'tickfont': {'size': 10}},
                      'yaxis': {'autorange': 'reversed', 'automargin': True,
                                'tickfont': {'size': 10}}}
        return self._plotly_div(traces, layout)

    # ------------------------------------------------------------------
    # Section 2: Results Plots
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Interactive volcano helper
    # ------------------------------------------------------------------

    def _volcano_html(self, div_id, vdata, title, x_label, y_label, eff_label,
                      table_id='', enr_dep_div='', enr_enr_div=''):
        """Generate interactive volcano with FDR and effect-size cutoff controls."""
        html = []
        html.append('<div class="vol-ctrl">')
        html.append('  FDR \u2264 <input type="number" value="0.05" step="0.01" min="0" max="1" '
                     'id="vfdr-' + div_id + '" oninput="drawVolcano(\'' + div_id + '\')">')
        html.append('  |' + eff_label + '| \u2265 <input type="number" value="0" step="0.1" min="0" '
                     'id="veff-' + div_id + '" oninput="drawVolcano(\'' + div_id + '\')">')
        html.append('</div>')
        html.append('<div id="' + div_id + '" class="plotly-chart"></div>')
        vdata['title'] = title
        vdata['xl'] = x_label
        vdata['yl'] = y_label
        vdata['tid'] = table_id
        vdata['ediv'] = enr_dep_div
        vdata['eenr'] = enr_enr_div
        html.append('<script>window["VD_' + div_id + '"]=' + _jdumps(vdata)
                    + ';document.addEventListener("DOMContentLoaded",function(){drawVolcano("'
                    + div_id + '");});</script>')
        return '\n'.join(html)

    def plot_rra_volcano(self, label, df, table_id='', enr_dep_div='', enr_enr_div=''):
        """Combined RRA volcano: neg+pos in one plot, interactive cutoffs."""
        div_id = self._next_id('vol')
        gene_col = 'id' if 'id' in df.columns else df.columns[0]
        genes = df[gene_col].values.astype(str)
        lfc = df['neg|lfc'].values.astype(float) if 'neg|lfc' in df.columns else np.zeros(len(df))
        neg_fdr = df['neg|fdr'].values.astype(float) if 'neg|fdr' in df.columns else np.ones(len(df))
        pos_fdr = df['pos|fdr'].values.astype(float) if 'pos|fdr' in df.columns else np.ones(len(df))
        y_fdr = np.where(lfc < 0, neg_fdr, pos_fdr)
        y = -np.log10(np.clip(y_fdr, 1e-300, 1.0))
        pos_mask = lfc > 0
        neg_mask = lfc < 0
        top_enr = np.where(pos_mask)[0][np.argsort(pos_fdr[pos_mask])][:20].tolist() if pos_mask.any() else []
        top_dep = np.where(neg_mask)[0][np.argsort(neg_fdr[neg_mask])][:20].tolist() if neg_mask.any() else []
        vdata = {'g': genes.tolist(),
                 'lfc': np.round(lfc, 2).tolist(),
                 'y': np.round(y, 1).tolist(),
                 'nf': np.round(neg_fdr, 3).tolist(),
                 'pf': np.round(pos_fdr, 3).tolist(),
                 'top_enr': top_enr, 'top_dep': top_dep}
        return self._volcano_html(div_id, vdata, label + ': Volcano',
                                  'Log2 Fold Change', '-log10(FDR)', 'LFC',
                                  table_id=table_id, enr_dep_div=enr_dep_div,
                                  enr_enr_div=enr_enr_div)

    def plot_mle_volcano(self, condition, table_id='', enr_dep_div='', enr_enr_div=''):
        """MLE volcano with -log10(FDR) y-axis and interactive cutoffs."""
        div_id = self._next_id('vol')
        df = self.mle_gene_summary
        beta_col = condition + '|beta'
        fdr_col = condition + '|wald-fdr'
        if fdr_col not in df.columns:
            fdr_col = condition + '|fdr'
        if beta_col not in df.columns or fdr_col not in df.columns:
            return None
        sample_name = self._condition_to_sample(condition)
        genes = df['Gene'].values.astype(str)
        beta = df[beta_col].values.astype(float)
        fdr = np.clip(df[fdr_col].values.astype(float), 1e-300, 1.0)
        y = -np.log10(fdr)
        pos_mask = beta > 0
        neg_mask = beta < 0
        top_enr = np.where(pos_mask)[0][np.argsort(fdr[pos_mask])][:20].tolist() if pos_mask.any() else []
        top_dep = np.where(neg_mask)[0][np.argsort(fdr[neg_mask])][:20].tolist() if neg_mask.any() else []
        vdata = {'g': genes.tolist(),
                 'lfc': np.round(beta, 2).tolist(),
                 'y': np.round(y, 1).tolist(),
                 'nf': np.round(fdr, 3).tolist(),
                 'pf': np.round(fdr, 3).tolist(),
                 'top_enr': top_enr, 'top_dep': top_dep}
        baseline_name = self._get_baseline_name()
        title = condition + ' vs ' + baseline_name + ' (Normalized Beta)'
        return self._volcano_html(div_id, vdata, title,
                                  'Normalized Beta Score', '\u2212log\u2081\u2080(FDR)', 'Beta',
                                  table_id=table_id, enr_dep_div=enr_dep_div,
                                  enr_enr_div=enr_enr_div)

    def plot_mle_beta_scatter(self, cond_x, cond_y):
        """Beta scatter with enriched/depleted coloring and top 10 labels."""
        df = self.mle_gene_summary
        bx = cond_x + '|beta'
        by = cond_y + '|beta'
        fx = cond_x + '|wald-fdr'
        fy = cond_y + '|wald-fdr'
        if bx not in df.columns or by not in df.columns:
            return None
        if fx not in df.columns:
            fx = cond_x + '|fdr'
        if fy not in df.columns:
            fy = cond_y + '|fdr'
        name_x = self._condition_to_sample(cond_x)
        name_y = self._condition_to_sample(cond_y)
        genes = df['Gene'].values.astype(str)
        vx = df[bx].values.astype(float)
        vy = df[by].values.astype(float)
        # Significance: either condition
        sig = np.zeros(len(df), dtype=bool)
        if fx in df.columns:
            sig |= (df[fx].values.astype(float) < self.fdr_threshold)
        if fy in df.columns:
            sig |= (df[fy].values.astype(float) < self.fdr_threshold)
        # Categories based on treatment (cond_y) beta direction
        depleted = sig & (vy < 0)
        enriched = sig & (vy > 0)
        ns = ~(depleted | enriched)
        vx = np.round(vx, 2)
        vy = np.round(vy, 2)
        traces = []
        ht = 'Gene: %{text}<br>' + name_x + ': %{x:.3f}<br>' + name_y + ': %{y:.3f}'
        if ns.any():
            traces.append({'type': 'scatter', 'mode': 'markers',
                           'x': vx[ns].tolist(), 'y': vy[ns].tolist(),
                           'text': genes[ns].tolist(),
                           'marker': {'color': 'rgba(187,187,187,0.3)', 'size': 4},
                           'hovertemplate': ht, 'name': 'Not significant'})
        if depleted.any():
            traces.append({'type': 'scatter', 'mode': 'markers',
                           'x': vx[depleted].tolist(), 'y': vy[depleted].tolist(),
                           'text': genes[depleted].tolist(),
                           'marker': {'color': '#3498DB', 'size': 7},
                           'hovertemplate': ht, 'name': 'Depleted'})
        if enriched.any():
            traces.append({'type': 'scatter', 'mode': 'markers',
                           'x': vx[enriched].tolist(), 'y': vy[enriched].tolist(),
                           'text': genes[enriched].tolist(),
                           'marker': {'color': '#E74C3C', 'size': 7},
                           'hovertemplate': ht, 'name': 'Enriched'})
        # Annotations: top 10 depleted and top 10 enriched
        annotations = []
        dep_idx = np.where(depleted)[0]
        if len(dep_idx) > 0:
            sorted_dep = dep_idx[np.argsort(vy[dep_idx])][:10]
            for k, idx in enumerate(sorted_dep):
                annotations.append({
                    'x': float(vx[idx]), 'y': float(vy[idx]),
                    'text': genes[idx], 'showarrow': True,
                    'arrowhead': 0, 'ax': 25, 'ay': -8 - (k % 3) * 12,
                    'font': {'size': 9, 'color': '#2166AC'}})
        enr_idx = np.where(enriched)[0]
        if len(enr_idx) > 0:
            sorted_enr = enr_idx[np.argsort(-vy[enr_idx])][:10]
            for k, idx in enumerate(sorted_enr):
                annotations.append({
                    'x': float(vx[idx]), 'y': float(vy[idx]),
                    'text': genes[idx], 'showarrow': True,
                    'arrowhead': 0, 'ax': -25, 'ay': -8 - (k % 3) * 12,
                    'font': {'size': 9, 'color': '#B2182B'}})
        lims = [min(vx.min(), vy.min()) - 0.2, max(vx.max(), vy.max()) + 0.2]
        layout = {'title': {'text': name_x + ' vs ' + name_y + ' Beta Scores', 'font': {'size': 13}},
                  'xaxis': {'title': name_x + ' (beta)', 'automargin': True,
                            'tickfont': {'size': 10}},
                  'yaxis': {'title': name_y + ' (beta)', 'automargin': True,
                            'tickfont': {'size': 10}},
                  'legend': {'font': {'size': 10}},
                  'height': 450, 'margin': {'t': 45, 'b': 45, 'l': 55, 'r': 15},
                  'annotations': annotations,
                  'shapes': [{'type': 'line', 'x0': lims[0], 'x1': lims[1],
                              'y0': lims[0], 'y1': lims[1],
                              'line': {'dash': 'dash', 'color': 'grey', 'width': 1}}]}
        return self._plotly_div(traces, layout)

    # ------------------------------------------------------------------
    # FLUTE-like MLE plots
    # ------------------------------------------------------------------

    def plot_mle_beta_density(self):
        """Plot 1: Density distribution of beta scores for all MLE conditions."""
        df = self.mle_gene_summary
        if df is None or not self.mle_conditions:
            return None
        colors = ['#F8766D', '#00BFC4', '#7CAE00', '#C77CFF', '#FF61CC']
        traces = []
        for i, cond in enumerate(self.mle_conditions):
            beta_col = cond + '|beta'
            if beta_col not in df.columns:
                continue
            vals = df[beta_col].dropna().values.astype(float)
            counts, edges = np.histogram(vals, bins=80, density=True)
            centers = ((edges[:-1] + edges[1:]) / 2).tolist()
            color = colors[i % len(colors)]
            sample_name = self._condition_to_sample(cond)
            traces.append({
                'type': 'scatter', 'mode': 'lines',
                'x': centers, 'y': counts.tolist(),
                'name': sample_name,
                'line': {'color': color, 'width': 2.5}})
        if not traces:
            return None
        layout = {
            'title': {'text': 'Distribution of Gene Beta Scores', 'font': {'size': 13}},
            'xaxis': {'title': 'Beta Score', 'automargin': True, 'tickfont': {'size': 10}},
            'yaxis': {'title': 'Density', 'automargin': True, 'tickfont': {'size': 10}},
            'legend': {'font': {'size': 11}},
            'height': 360, 'margin': {'t': 45, 'b': 40, 'l': 50, 'r': 10}}
        return self._plotly_div(traces, layout)

    def plot_mle_consistency_scatter(self, cond_ctrl, cond_treat):
        """Plot 2: ConsistencyView scatter — treatment vs control with regression."""
        df = self.mle_gene_summary
        bx_col = cond_ctrl + '|beta'
        by_col = cond_treat + '|beta'
        if bx_col not in df.columns or by_col not in df.columns:
            return None
        ctrl_name = self._condition_to_sample(cond_ctrl)
        treat_name = self._condition_to_sample(cond_treat)
        genes = df['Gene'].values.astype(str)
        vx = df[bx_col].values.astype(float)
        vy = df[by_col].values.astype(float)
        mask = np.isfinite(vx) & np.isfinite(vy)
        vx_f, vy_f = vx[mask], vy[mask]
        # Linear regression
        if len(vx_f) > 2:
            coeffs = np.polyfit(vx_f, vy_f, 1)
            slope, intercept = coeffs[0], coeffs[1]
        else:
            slope, intercept = 1.0, 0.0
        vx_f = np.round(vx_f, 2)
        vy_f = np.round(vy_f, 2)
        lims = [min(vx_f.min(), vy_f.min()) - 0.2, max(vx_f.max(), vy_f.max()) + 0.2]
        reg_x = [lims[0], lims[1]]
        reg_y = [slope * lims[0] + intercept, slope * lims[1] + intercept]
        ht = 'Gene: %{text}<br>' + ctrl_name + ': %{x:.3f}<br>' + treat_name + ': %{y:.3f}'
        traces = [
            {'type': 'scatter', 'mode': 'markers',
             'x': vx_f.tolist(), 'y': vy_f.tolist(),
             'text': genes[mask].tolist(),
             'marker': {'color': '#4A90D9', 'size': 4, 'opacity': 0.5},
             'hovertemplate': ht, 'name': 'Genes', 'showlegend': False},
            {'type': 'scatter', 'mode': 'lines',
             'x': reg_x, 'y': reg_y,
             'line': {'color': '#E74C3C', 'width': 2},
             'name': 'Regression (slope={:.3f})'.format(slope)},
        ]
        layout = {
            'title': {'text': 'ConsistencyView: ' + ctrl_name + ' vs ' + treat_name, 'font': {'size': 13}},
            'xaxis': {'title': ctrl_name + ' (beta)', 'automargin': True, 'tickfont': {'size': 10}},
            'yaxis': {'title': treat_name + ' (beta)', 'automargin': True, 'tickfont': {'size': 10}},
            'legend': {'font': {'size': 10}},
            'height': 360, 'margin': {'t': 45, 'b': 40, 'l': 50, 'r': 10},
            'shapes': [{'type': 'line', 'x0': lims[0], 'x1': lims[1],
                        'y0': lims[0], 'y1': lims[1],
                        'line': {'dash': 'dash', 'color': 'grey', 'width': 1}}]}
        return self._plotly_div(traces, layout)

    def plot_mle_selection_scatter(self, cond_ctrl, cond_treat):
        """Plot 3: Interactive scatter with user-adjustable diagonal cutoff."""
        df = self.mle_gene_summary
        bx_col = cond_ctrl + '|beta'
        by_col = cond_treat + '|beta'
        if bx_col not in df.columns or by_col not in df.columns:
            return None
        ctrl_name = self._condition_to_sample(cond_ctrl)
        treat_name = self._condition_to_sample(cond_treat)
        genes = df['Gene'].values.astype(str)
        vx = df[bx_col].values.astype(float)
        vy = df[by_col].values.astype(float)
        mask = np.isfinite(vx) & np.isfinite(vy)
        diff_std = float(np.nanstd(vy[mask] - vx[mask]))
        default_cut = round(0.5 * diff_std, 2) if diff_std > 0 else 0.1
        div_id = self._next_id('sel')
        sdata = {
            'g': genes[mask].tolist(),
            'vx': np.round(vx[mask], 2).tolist(),
            'vy': np.round(vy[mask], 2).tolist(),
            'ctrl': ctrl_name, 'treat': treat_name,
        }
        html = []
        html.append('<div class="vol-ctrl">')
        html.append('  Diagonal cutoff \u2265 <input type="number" value="' + str(default_cut)
                     + '" step="0.05" min="0" id="scut-' + div_id
                     + '" oninput="drawSelection(\'' + div_id + '\')">')
        html.append('</div>')
        html.append('<div id="' + div_id + '" class="plotly-chart"></div>')
        html.append('<script>window["SD_' + div_id + '"]=' + _jdumps(sdata)
                    + ';document.addEventListener("DOMContentLoaded",function(){drawSelection("'
                    + div_id + '");});</script>')
        return '\n'.join(html)

    def plot_mle_nine_square(self, cond_ctrl, cond_treat):
        """Plot 4: FLUTE nine-square treatment gene classification."""
        df = self.mle_gene_summary
        bx_col = cond_ctrl + '|beta'
        by_col = cond_treat + '|beta'
        if bx_col not in df.columns or by_col not in df.columns:
            return None
        ctrl_name = self._condition_to_sample(cond_ctrl)
        treat_name = self._condition_to_sample(cond_treat)
        genes = df['Gene'].values.astype(str)
        vx = df[bx_col].values.astype(float)
        vy = df[by_col].values.astype(float)
        mask = np.isfinite(vx) & np.isfinite(vy)
        # Default cutoff: use the standard deviation of absolute beta values
        beta_abs = np.abs(np.concatenate([vx[mask], vy[mask]]))
        default_cut = round(float(np.percentile(beta_abs[beta_abs > 0], 75)), 2) if len(beta_abs[beta_abs > 0]) > 0 else 0.5
        div_id = self._next_id('nsq')
        ndata = {
            'g': genes[mask].tolist(),
            'vx': np.round(vx[mask], 2).tolist(),
            'vy': np.round(vy[mask], 2).tolist(),
            'ctrl': ctrl_name, 'treat': treat_name,
        }
        html = []
        html.append('<div class="vol-ctrl">')
        html.append('  Beta cutoff |&beta;| \u2265 <input type="number" value="' + str(default_cut)
                     + '" step="0.05" min="0" id="ncut-' + div_id
                     + '" oninput="drawNineSquare(\'' + div_id + '\')">')
        html.append('</div>')
        html.append('<div id="' + div_id + '" class="plotly-chart"></div>')
        html.append('<script>window["ND_' + div_id + '"]=' + _jdumps(ndata)
                    + ';document.addEventListener("DOMContentLoaded",function(){drawNineSquare("'
                    + div_id + '");});</script>')
        return '\n'.join(html)

    # ------------------------------------------------------------------
    # Enrichment (decoupler)
    # ------------------------------------------------------------------

    def _fetch_hallmark_net(self):
        if self._hallmark_net is not None:
            return self._hallmark_net
        if not HAS_DECOUPLER:
            logging.warning('decoupler not installed; skipping enrichment.')
            return None
        try:
            self._hallmark_net = dc.op.hallmark(organism=self.organism)
            logging.info('Fetched %d hallmark gene-set entries.', len(self._hallmark_net))
            return self._hallmark_net
        except Exception as e:
            logging.warning('Could not fetch hallmark gene sets: %s', e)
            return None

    def _run_ora(self, gene_list, all_genes):
        """Run ORA on gene_list against hallmark gene sets."""
        net = self._fetch_hallmark_net()
        if net is None or len(gene_list) == 0:
            return None
        sig_set = set(gene_list)
        values = [1.0 if g in sig_set else 0.0 for g in all_genes]
        mat = pd.DataFrame([values], columns=list(all_genes), index=['sample'])
        try:
            result = dc.mt.ora(mat, net, tmin=1)
            if isinstance(result, tuple) and len(result) >= 2:
                scores_df, padjs_df = result[0], result[1]
            else:
                scores_df = result
                padjs_df = None
            res = pd.DataFrame({'source': scores_df.columns,
                                'score': scores_df.iloc[0].values})
            if padjs_df is not None:
                res['padj'] = padjs_df.iloc[0].values
            else:
                res['padj'] = np.nan
            res = res.dropna(subset=['padj']).sort_values('padj').reset_index(drop=True)
            return res
        except Exception as e:
            logging.warning('ORA enrichment failed: %s', e)
            return None

    def plot_enrichment_dotplot(self, ora_df, title, direction='depleted', div_id=None):
        """Enrichment dot plot with direction-specific coloring and download.

        direction: 'depleted' or 'enriched'
        Returns (html_string, div_id) tuple.
        """
        if ora_df is None or len(ora_df) == 0:
            return ('<p class="plot-desc" style="color:#999;">No enrichment results available for this gene set.</p>', '')
        sig = ora_df[ora_df['padj'] < self.enrichment_fdr]
        show = sig.sort_values('score', ascending=False).head(15) if len(sig) > 0 else ora_df.sort_values('score', ascending=False).head(15)
        if len(show) == 0:
            return ('<p class="plot-desc" style="color:#999;">No significant pathways found (adj. p-value &lt; ' + str(self.enrichment_fdr) + ').</p>', '')
        show = show.iloc[::-1]  # reverse for bottom-to-top (highest score at top in Plotly)

        neg_log_p = (-np.log10(show['padj'].clip(1e-300, 1.0))).tolist()
        max_nlp = max(neg_log_p) if neg_log_p else 1
        min_nlp = min(neg_log_p) if neg_log_p else 0
        sizes = [max(8, v / max_nlp * 28) for v in neg_log_p]

        # Direction-specific color scheme
        is_cold = direction == 'depleted'
        cscale = _CSCALE_DEPLETED if is_cold else _CSCALE_ENRICHED
        marker_line_color = '#2166AC' if is_cold else '#B2182B'

        # Truncate long pathway names for readability
        y_labels = [_trunc(n) for n in show['source'].tolist()]

        traces = [{'type': 'scatter', 'mode': 'markers',
                   'x': show['score'].tolist(),
                   'y': y_labels,
                   'marker': {'size': sizes,
                              'color': show['score'].tolist(),
                              'colorscale': cscale,
                              'colorbar': {'title': {'text': 'ORA Score', 'font': {'size': 9}},
                                           'len': 0.4, 'thickness': 12,
                                           'tickfont': {'size': 8},
                                           'x': 0.82, 'xpad': 3,
                                           'y': 0.8, 'yanchor': 'middle'},
                              'line': {'width': 1.2, 'color': marker_line_color},
                              'opacity': 0.85},
                   'text': ['-log10(padj)={:.2f}'.format(v) for v in neg_log_p],
                   'hovertemplate': '%{y}<br>ORA Score: %{x:.3f}<br>%{text}',
                   'showlegend': False}]

        # Size legend: reference dots as paper-coordinate annotations (right side, near colorbar)
        ref_vals = []
        if max_nlp > 2:
            ref_vals = [round(min_nlp + (max_nlp - min_nlp) * f, 1) for f in [0.0, 0.5, 1.0]]
        elif max_nlp > 0:
            ref_vals = [round(min_nlp, 1), round(max_nlp, 1)]
        seen_rv = set()
        unique_rv = []
        for rv in ref_vals:
            if rv not in seen_rv:
                seen_rv.add(rv)
                unique_rv.append(rv)
        ref_vals = unique_rv

        # Build size legend as scatter on xaxis2 placed to the right (near colorbar)
        if ref_vals:
            ref_sizes = [max(8, v / max_nlp * 28) for v in ref_vals]
            y_pos = list(range(len(ref_vals)))
            traces.append({
                'type': 'scatter', 'mode': 'markers+text',
                'x': [0] * len(ref_vals), 'y': y_pos,
                'marker': {'size': ref_sizes,
                           'color': marker_line_color,
                           'opacity': 0.5,
                           'line': {'width': 1, 'color': marker_line_color}},
                'text': ['%.1f' % v for v in ref_vals],
                'textposition': 'middle right',
                'textfont': {'size': 8, 'color': '#555'},
                'hoverinfo': 'skip',
                'showlegend': False,
                'xaxis': 'x2', 'yaxis': 'y2',
            })

        annotations = []
        if ref_vals:
            annotations.append({
                'text': '<b>-log\u2081\u2080(p)</b>',
                'xref': 'x2 domain', 'yref': 'y2 domain',
                'x': 0.5, 'y': 1.15, 'showarrow': False,
                'font': {'size': 8, 'color': '#555'}})

        n_pathways = len(show)
        plot_h = max(270, n_pathways * 24 + 100)

        layout = {'title': {'text': title, 'font': {'size': 12},
                            'y': 0.98, 'yanchor': 'top', 'pad': {'t': 10}},
                  'xaxis': {'title': 'ORA Enrichment Score', 'automargin': True,
                            'tickfont': {'size': 9},
                            'domain': [0, 0.78]},
                  'yaxis': {'automargin': True, 'tickfont': {'size': 9},
                            'domain': [0, 1.0]},
                  'xaxis2': {'anchor': 'y2', 'domain': [0.88, 0.98],
                             'showgrid': False, 'zeroline': False,
                             'showticklabels': False, 'fixedrange': True},
                  'yaxis2': {'anchor': 'x2', 'domain': [0.05, 0.45],
                             'showgrid': False, 'zeroline': False,
                             'showticklabels': False, 'fixedrange': True},
                  'height': plot_h,
                  'margin': {'l': 10, 't': 50, 'b': 30, 'r': 80},
                  'annotations': annotations}

        if div_id is None:
            div_id = self._next_id('enr')
        plot_html = self._plotly_div(traces, layout, div_id=div_id)

        # Store enrichment data for download
        self._enrichment_counter += 1
        ekey = 'enr_%d' % self._enrichment_counter
        store_rows = []
        for _, row in ora_df.iterrows():
            store_rows.append({
                'Pathway': row['source'],
                'ORA_Score': float(row['score']) if pd.notna(row['score']) else None,
                'Adjusted_P_Value': float(row['padj']) if pd.notna(row['padj']) else None,
            })
        self._enrichment_store[ekey] = store_rows
        safe_fn = title.replace(' ', '_').replace(':', '').replace('|', '_')
        dl_btn = ('<div style="text-align:right;margin-top:-5px;margin-bottom:10px;">'
                  '<button class="download-btn-sm" onclick="dlEnr(\'' + ekey + '\',\'' + safe_fn + '.csv\')">'
                  'Download enrichment data</button></div>\n')

        return (plot_html + dl_btn, div_id)

    # ------------------------------------------------------------------
    # Dynamic tables
    # ------------------------------------------------------------------

    def _generate_dynamic_table(self, df, table_id, gene_col='id',
                                fdr_cols=None, lineplot_id='',
                                eff_col=None):
        """Generate a dynamic filterable/sortable table with all genes.

        eff_col: column key used for effect-size filtering (|LFC| or |beta|).
        """
        self._table_ids.append(table_id)
        columns = []
        # Detect an effect-size column automatically if not given
        auto_eff = None
        for col in df.columns:
            fmt = 's'
            is_fdr = False
            if col == gene_col:
                fmt = 's'
            elif col in ('num', 'sgRNA') or 'rank' in col.lower() or 'goodsgrna' in col.lower():
                fmt = 'd'
            elif 'p-value' in col or 'fdr' in col or 'score' in col:
                fmt = '.2e'
                if 'fdr' in col:
                    is_fdr = True
            elif 'lfc' in col or 'beta' in col or '|z' in col:
                fmt = '.4f'
                if auto_eff is None and ('lfc' in col or 'beta' in col):
                    auto_eff = col
            label = col
            label = label.replace('neg|', 'Neg ').replace('pos|', 'Pos ')
            for old, new in [('|beta', ' Beta'), ('|z', ' Z'), ('|p-value', ' P'),
                             ('|fdr', ' FDR'), ('|wald-p-value', ' Wald-P'),
                             ('|wald-fdr', ' Wald-FDR')]:
                label = label.replace(old, new)
            if col == 'id':
                label = 'Gene'
            elif col == 'num':
                label = '#sgRNAs'
            columns.append({'key': col, 'label': label, 'fmt': fmt, 'isFdr': is_fdr})

        if fdr_cols is None:
            fdr_cols = [c['key'] for c in columns if c['isFdr']]
        if eff_col is None:
            eff_col = auto_eff or ''

        # Compact column-oriented format: {keys:[], vals:[[row1],[row2],...]}
        # Much smaller than records format since column names aren't repeated
        df_slim = df.copy()
        for col in df_slim.columns:
            if df_slim[col].dtype in ('float64', 'float32'):
                df_slim[col] = df_slim[col].round(3)
        col_keys = [c['key'] for c in columns]
        col_vals = []
        for _, row in df_slim.iterrows():
            col_vals.append([_json_safe(row[k]) for k in col_keys])
        data_json = json.dumps({'k': col_keys, 'v': col_vals}, separators=(',', ':'), default=str)
        cols_json = json.dumps(columns, separators=(',', ':'))
        fdr_keys_json = json.dumps(fdr_cols, separators=(',', ':'))

        thead = '<tr>' + ''.join('<th>' + c['label'] + '</th>' for c in columns) + '</tr>'

        html_parts = []
        html_parts.append('<div class="table-container">')
        html_parts.append('<div class="table-controls">')
        html_parts.append('  <label>FDR &le; <input type="number" id="fdr-' + table_id + '" value="1" step="0.01" min="0" max="1"></label>')
        html_parts.append('  <label>|Effect| &ge; <input type="number" id="eff-' + table_id + '" value="0" step="0.1" min="0"></label>')
        html_parts.append('  <input type="text" id="search-' + table_id + '" placeholder="Search gene...">')
        html_parts.append('  <span id="count-' + table_id + '"></span>')
        html_parts.append('  <select id="psize-' + table_id + '">')
        html_parts.append('    <option value="25">25/page</option><option value="50" selected>50/page</option>')
        html_parts.append('    <option value="100">100/page</option><option value="0">All</option>')
        html_parts.append('  </select>')
        html_parts.append('  <button class="download-btn-sm" onclick="downloadCSV(\'' + table_id + '\',\'' + table_id + '.csv\')">Export All CSV</button>')
        html_parts.append('  <button class="download-btn-sm" onclick="downloadFilteredCSV(\'' + table_id + '\',\'' + table_id + '_filtered.csv\')">Export Filtered CSV</button>')
        html_parts.append('</div>')
        html_parts.append('<div class="table-scroll">')
        html_parts.append('<table id="tbl-' + table_id + '"><thead>' + thead + '</thead><tbody></tbody></table>')
        html_parts.append('</div>')
        html_parts.append('<div class="table-pag">')
        html_parts.append('  <button id="prev-' + table_id + '">Prev</button>')
        html_parts.append('  <span id="pinfo-' + table_id + '"></span>')
        html_parts.append('  <button id="next-' + table_id + '">Next</button>')
        html_parts.append('</div>')
        html_parts.append('</div>')

        js = _TABLE_JS_IIFE % {
            'TID': json.dumps(table_id),
            'GK': json.dumps(gene_col),
            'LP': json.dumps(lineplot_id),
            'EK': json.dumps(eff_col),
            'COLS': cols_json,
            'DATA': data_json,
            'FKEYS': fdr_keys_json,
        }
        html_parts.append('<script>' + js + '</script>')
        return '\n'.join(html_parts)

    @staticmethod
    def _df_to_html_table(df):
        """Convert a small DataFrame to a simple HTML table."""
        html = '<table class="simple-table">\n<thead><tr>'
        for col in df.columns:
            html += '<th>' + str(col) + '</th>'
        html += '</tr></thead>\n<tbody>\n'
        for _, row in df.iterrows():
            html += '<tr>'
            for col in df.columns:
                html += '<td>' + str(row[col]) + '</td>'
            html += '</tr>\n'
        html += '</tbody></table>'
        return html

    # ------------------------------------------------------------------
    # Default gene selection helpers
    # ------------------------------------------------------------------

    def _get_default_genes_rra(self, df):
        """Return top 5 depleted + top 5 enriched genes for RRA."""
        gene_col = 'id' if 'id' in df.columns else df.columns[0]
        result = []
        for direction in ['neg', 'pos']:
            rank_col = direction + '|rank'
            fdr_col = direction + '|fdr'
            if rank_col in df.columns:
                top = df.sort_values(rank_col).head(5)[gene_col].tolist()
            elif fdr_col in df.columns:
                top = df.sort_values(fdr_col).head(5)[gene_col].tolist()
            else:
                top = []
            result.extend(top)
        # deduplicate while preserving order
        seen = set()
        unique = []
        for g in result:
            if g not in seen:
                seen.add(g)
                unique.append(g)
        return unique

    def _get_default_genes_mle(self, condition=None):
        """Return top 5 depleted + top 5 enriched genes for MLE."""
        df = self.mle_gene_summary
        if df is None:
            return []
        if condition is None and self.mle_conditions:
            condition = self.mle_conditions[0]
        if condition is None:
            return []
        beta_col = condition + '|beta'
        fdr_col = condition + '|wald-fdr'
        if fdr_col not in df.columns:
            fdr_col = condition + '|fdr'
        if beta_col not in df.columns:
            return []
        result = []
        # Top 5 depleted (most negative beta, lowest FDR)
        neg = df[df[beta_col] < 0]
        if fdr_col in df.columns and len(neg) > 0:
            result.extend(neg.sort_values(fdr_col).head(5)['Gene'].tolist())
        # Top 5 enriched (most positive beta, lowest FDR)
        pos = df[df[beta_col] > 0]
        if fdr_col in df.columns and len(pos) > 0:
            result.extend(pos.sort_values(fdr_col).head(5)['Gene'].tolist())
        # deduplicate
        seen = set()
        unique = []
        for g in result:
            if g not in seen:
                seen.add(g)
                unique.append(g)
        return unique

    # ------------------------------------------------------------------
    # Section builders
    # ------------------------------------------------------------------

    def _build_sample_section(self):
        html = []
        now_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

        # Report metadata banner
        html.append('<div class="report-meta">')
        html.append('<table class="info-table">')
        html.append('<tr><td>Report Generated:</td><td>' + now_str + '</td></tr>')
        html.append('<tr><td>MAGeCK Version:</td><td>' + __version__ + '</td></tr>')
        html.append('<tr><td>Analysis Prefix:</td><td>' + os.path.basename(self.output_prefix) + '</td></tr>')
        html.append('<tr><td>Normalization:</td><td>' + self.norm_method.capitalize() + '-ratio</td></tr>')
        html.append('<tr><td>FDR Threshold:</td><td>' + str(self.fdr_threshold) + '</td></tr>')
        if not self.skip_enrichment:
            html.append('<tr><td>Enrichment:</td><td>MSigDB Hallmark (' + self.organism + '), top ' + str(self.enrichment_top_n) + ' genes, ORA adj. p &lt; ' + str(self.enrichment_fdr) + '</td></tr>')
        html.append('</table></div>')

        if self.count_df is not None:
            n_sgrna = len(self.count_df)
            n_genes = self.count_df['Gene'].nunique() if 'Gene' in self.count_df.columns else 'N/A'
            html.append('<h3>Library Overview</h3>')
            html.append('<div class="summary-box">')
            html.append('<table class="info-table">')
            html.append('<tr><td>Total sgRNAs:</td><td>{:,}</td></tr>'.format(n_sgrna))
            html.append('<tr><td>Total Genes:</td><td>{:,}</td></tr>'.format(n_genes))
            html.append('<tr><td>Samples:</td><td>{} ({})</td></tr>'.format(
                len(self.sample_names), ', '.join(self.sample_names)))
            if self.count_table_path:
                html.append('<tr><td>Count Table:</td><td>{}</td></tr>'.format(
                    os.path.basename(self.count_table_path)))
            html.append('</table></div>')

        if self.rra_gene_summaries:
            html.append('<h3>RRA Comparisons Detected</h3><ul>')
            for label in sorted(self.rra_gene_summaries.keys()):
                df = self.rra_gene_summaries[label]
                n_g = len(df)
                html.append('<li><strong>' + label + '</strong> (' + str(n_g) + ' genes)</li>')
            html.append('</ul>')

        if self.mle_gene_summary is not None:
            html.append('<h3>MLE Conditions Detected</h3>')
            html.append('<ul>')
            for cond in self.mle_conditions:
                html.append('<li><strong>' + cond + '</strong></li>')
            html.append('</ul>')

        if html:
            self.sections.append(('sample-desc', 'Sample Description', '\n'.join(html)))

    def _build_qc_section(self):
        if self.count_df is None:
            return
        html = []

        # Mapping statistics (from mageck count)
        if self.countsummary_df is not None:
            html.append('<h3>1.0 Mapping Statistics</h3>')
            html.append(_desc_html('mapping_stats'))
            bar_html, tbl_html = self.plot_mapping_stats()
            if bar_html and tbl_html:
                html.append('<div class="plot-grid">')
                html.append('<div>' + bar_html + '</div>')
                html.append('<div>' + tbl_html + '</div>')
                html.append('</div>')

        # Row 1: Library size + Zero-count (side by side)
        html.append('<h3>1.1 Library Size &amp; Zero-Count sgRNAs</h3>')
        html.append('<div class="plot-grid">')
        html.append('<div>')
        html.append(_desc_html('lib_size'))
        html.append(self.plot_library_size())
        html.append('</div>')
        html.append('<div>')
        html.append(_desc_html('zero_count'))
        html.append(self.plot_zero_count_guides())
        html.append('</div>')
        html.append('</div>')

        # Row 2: Count distribution + Box plot (side by side)
        html.append('<h3>1.2 Count Distribution</h3>')
        hist_html, box_html = self.plot_count_distribution()
        html.append('<div class="plot-grid">')
        html.append('<div>')
        html.append(_desc_html('count_dist'))
        html.append(hist_html)
        html.append('</div>')
        html.append('<div>')
        html.append(_desc_html('density_norm'))
        html.append(box_html)
        html.append('</div>')
        html.append('</div>')

        # Row 3: Gini + PCA (side by side)
        if len(self.sample_names) >= 2:
            html.append('<h3>1.3 Library Evenness &amp; Sample Clustering</h3>')
            html.append('<div class="plot-grid">')
            html.append('<div>')
            html.append(_desc_html('gini'))
            html.append(self.plot_gini_index())
            html.append('</div>')
            html.append('<div>')
            html.append(_desc_html('pca'))
            html.append(self.plot_pca())
            html.append('</div>')
            html.append('</div>')

            # Row 4: Correlation (full width but compact)
            html.append('<h3>1.4 Sample Correlation</h3>')
            html.append(_desc_html('corr'))
            html.append(self.plot_sample_correlation())
        else:
            html.append('<h3>1.3 Library Evenness</h3>')
            html.append(_desc_html('gini'))
            html.append(self.plot_gini_index())

        self.sections.append(('qc', 'Quality Control', '\n'.join(html)))

    def _build_results_section(self):
        html = []
        has_counts = self.count_df is not None

        # --- RRA ---
        if self.rra_gene_summaries:
            html.append('<h3>2.1 Rank Aggregation (RRA) Results</h3>')
            for label, df in sorted(self.rra_gene_summaries.items()):
                html.append('<h4>' + label + '</h4>')

                # FDR summary
                for direction, dlabel in [('neg', 'Depleted'), ('pos', 'Enriched')]:
                    fdr_col = direction + '|fdr'
                    if fdr_col in df.columns:
                        fdrs = df[fdr_col].values.astype(float)
                        html.append(
                            '<div class="summary-box"><strong>' + dlabel + ':</strong> '
                            'FDR &lt; 1%: ' + str((fdrs < 0.01).sum()) + ' | '
                            'FDR &lt; 5%: ' + str((fdrs < 0.05).sum()) + ' | '
                            'FDR &lt; 25%: ' + str((fdrs < 0.25).sum()) + '</div>')

                # Full gene table
                tid = 'rra-' + label.replace(' ', '_')
                lp_id = 'lineplot-rra-' + label.replace(' ', '_') if has_counts else ''
                html.append('<h5>Gene Summary Table</h5>')
                html.append(_desc_html('gene_table'))
                html.append(self._generate_dynamic_table(df, tid, gene_col='id', lineplot_id=lp_id))

                # Pre-create enrichment div IDs for volcano linkage
                enr_dep_div = self._next_id('enr')
                enr_enr_div = self._next_id('enr')

                # Volcano + sgRNA line plot side by side
                html.append('<div class="plot-grid">')
                html.append('<div>')
                html.append('<h5>Volcano Plot</h5>')
                html.append(_desc_html('rra_volcano'))
                v = self.plot_rra_volcano(label, df, table_id=tid,
                                          enr_dep_div=enr_dep_div, enr_enr_div=enr_enr_div)
                if v:
                    html.append(v)
                html.append('</div>')
                if has_counts:
                    default_genes = self._get_default_genes_rra(df)
                    if lp_id:
                        self._default_genes[lp_id] = default_genes
                    html.append('<div>')
                    html.append('<h5>sgRNA log2 Fold Change '
                                '<button class="download-btn-sm" onclick="clearLP(\'' + lp_id + '\')" '
                                'style="margin-left:10px;">Clear Selection</button></h5>')
                    html.append(_desc_html('sgrna_line'))
                    html.append('<div id="' + lp_id + '" class="plotly-chart" style="height:380px;"></div>')
                    html.append('<script>Plotly.newPlot("' + lp_id + '",[],{title:{text:"Loading default genes..."},xaxis:{title:"Sample"},yaxis:{title:"log2 Fold Change",zeroline:true},height:380},' + _PLOTLY_CONFIG + ');</script>')
                    html.append('</div>')
                html.append('</div>')

                # Enrichment (side by side: depleted + enriched)
                if not self.skip_enrichment and HAS_DECOUPLER:
                    gene_col = 'id' if 'id' in df.columns else df.columns[0]
                    all_genes = df[gene_col].tolist()
                    html.append('<h5>Hallmark Pathway Enrichment</h5>')
                    html.append('<div class="plot-grid">')
                    for direction, dtitle, edir, ediv in [
                        ('neg', 'Depleted Genes', 'depleted', enr_dep_div),
                        ('pos', 'Enriched Genes', 'enriched', enr_enr_div)]:
                        rank_col = direction + '|rank'
                        if rank_col in df.columns:
                            top_genes = df.sort_values(rank_col).head(self.enrichment_top_n)[gene_col].tolist()
                        else:
                            fdr_c = direction + '|fdr'
                            top_genes = df.sort_values(fdr_c).head(self.enrichment_top_n)[gene_col].tolist()
                        ora_result = self._run_ora(top_genes, all_genes)
                        full_title = label + ': ' + dtitle
                        html.append('<div>')
                        html.append(_desc_html('enrich_depl' if direction == 'neg' else 'enrich_enr'))
                        enr_html, _ = self.plot_enrichment_dotplot(ora_result, full_title,
                                                                    direction=edir, div_id=ediv)
                        html.append(enr_html)
                        html.append('</div>')
                    html.append('</div>')

        # --- MLE ---
        if self.mle_gene_summary is not None and self.mle_conditions:
            baseline_name = self._get_baseline_name()
            html.append('<h3>2.2 Maximum Likelihood Estimation (MLE) Results</h3>')
            html.append('<div class="summary-box">Beta scores normalized using core '
                        'essential genes (FLUTE cell-cycle method). Baseline: <strong>'
                        + baseline_name + '</strong>. Conditions: <strong>'
                        + ', '.join(self.mle_conditions) + '</strong>.</div>')
            df = self.mle_gene_summary

            # Design matrix display
            if self.design_matrix_df is not None:
                html.append('<h4>Design Matrix</h4>')
                html.append('<p class="plot-desc">The design matrix used for MLE '
                            'estimation. Rows represent samples and columns represent '
                            'conditions being tested.</p>')
                html.append(self._df_to_html_table(self.design_matrix_df))

            # MLE Gene Summary Table FIRST (harmonized with RRA order)
            tid = 'mle-' + os.path.basename(self.output_prefix).replace(' ', '_')
            mle_lp_ids = []
            if has_counts:
                for cond in self.mle_conditions:
                    mle_lp_ids.append('lineplot-mle-' + cond.replace(' ', '_'))
                default_genes = self._get_default_genes_mle()
                for _lp in mle_lp_ids:
                    self._default_genes[_lp] = default_genes
            combined_lp = ','.join(mle_lp_ids)

            html.append('<h5>MLE Gene Summary Table</h5>')
            html.append(_desc_html('gene_table'))
            html.append(self._generate_dynamic_table(df, tid, gene_col='Gene', lineplot_id=combined_lp))

            # Per-condition: FDR summary + volcano + line plot + enrichment
            for ci_idx, cond in enumerate(self.mle_conditions):
                sample_name = self._condition_to_sample(cond)
                cond_label = cond + ' vs ' + baseline_name
                html.append('<h4>' + cond_label + '</h4>')
                fdr_col = cond + '|wald-fdr'
                if fdr_col not in df.columns:
                    fdr_col = cond + '|fdr'
                if fdr_col in df.columns:
                    fdrs = df[fdr_col].values.astype(float)
                    n_sig = (fdrs < self.fdr_threshold).sum()
                    html.append('<div class="summary-box"><strong>' + cond_label + ':</strong> '
                                + str(n_sig) + ' genes with FDR &lt; ' + str(self.fdr_threshold) + '</div>')

                # Pre-create enrichment div IDs
                enr_dep_div = self._next_id('enr')
                enr_enr_div = self._next_id('enr')

                # Volcano + sgRNA line plot side by side
                html.append('<div class="plot-grid">')
                html.append('<div>')
                html.append('<h5>Volcano Plot</h5>')
                html.append(_desc_html('mle_volcano'))
                v = self.plot_mle_volcano(cond, table_id=tid,
                                          enr_dep_div=enr_dep_div, enr_enr_div=enr_enr_div)
                if v:
                    html.append(v)
                html.append('</div>')
                if has_counts and ci_idx < len(mle_lp_ids):
                    cond_lp = mle_lp_ids[ci_idx]
                    html.append('<div>')
                    html.append('<h5>sgRNA log2 Fold Change '
                                '<button class="download-btn-sm" onclick="clearLP(\'' + combined_lp + '\')" '
                                'style="margin-left:10px;">Clear Selection</button></h5>')
                    html.append(_desc_html('sgrna_line'))
                    html.append('<div id="' + cond_lp + '" class="plotly-chart" style="height:380px;"></div>')
                    html.append('<script>Plotly.newPlot("' + cond_lp + '",[],{title:{text:"Loading default genes..."},xaxis:{title:"Sample"},yaxis:{title:"log2 Fold Change",zeroline:true},height:380},' + _PLOTLY_CONFIG + ');</script>')
                    html.append('</div>')
                html.append('</div>')

                # Enrichment (side by side: depleted + enriched)
                if not self.skip_enrichment and HAS_DECOUPLER:
                    all_genes = df['Gene'].tolist()
                    beta_col = cond + '|beta'
                    if beta_col in df.columns and fdr_col in df.columns:
                        html.append('<h5>Pathway Enrichment: ' + cond_label + '</h5>')
                        html.append('<div class="plot-grid">')
                        # Depleted (negative beta)
                        neg_df = df[df[beta_col] < 0].sort_values(fdr_col).head(self.enrichment_top_n)
                        ora_neg = self._run_ora(neg_df['Gene'].tolist(), all_genes)
                        html.append('<div>')
                        html.append(_desc_html('enrich_sens'))
                        enr_html, _ = self.plot_enrichment_dotplot(
                            ora_neg, cond_label + ': Depleted Genes',
                            direction='depleted', div_id=enr_dep_div)
                        html.append(enr_html)
                        html.append('</div>')
                        # Enriched (positive beta)
                        pos_df = df[df[beta_col] > 0].sort_values(fdr_col).head(self.enrichment_top_n)
                        ora_pos = self._run_ora(pos_df['Gene'].tolist(), all_genes)
                        html.append('<div>')
                        html.append(_desc_html('enrich_res'))
                        enr_html, _ = self.plot_enrichment_dotplot(
                            ora_pos, cond_label + ': Enriched Genes',
                            direction='enriched', div_id=enr_enr_div)
                        html.append(enr_html)
                        html.append('</div>')
                        html.append('</div>')

            # FLUTE-like plots: Beta density, consistency, selection, nine-square (2-column grid)
            if len(self.mle_conditions) >= 1:
                html.append('<h4>Beta Score Analysis (FLUTE-style)</h4>')

            if len(self.mle_conditions) >= 2:
                # Identify treatment/control pair
                cond_ctrl = None
                cond_treat = None
                for c in self.mle_conditions:
                    cl = c.lower()
                    if 'dmso' in cl or 'control' in cl or 'baseline' in cl or 'untreated' in cl:
                        cond_ctrl = c
                if cond_ctrl is None:
                    cond_ctrl = self.mle_conditions[0]
                for c in self.mle_conditions:
                    if c != cond_ctrl:
                        cond_treat = c
                        break

                if cond_ctrl and cond_treat:
                    # Row 1: Density + ConsistencyView
                    html.append('<div class="plot-grid">')
                    html.append('<div>')
                    html.append(_desc_html('mle_beta_density'))
                    v = self.plot_mle_beta_density()
                    if v:
                        html.append(v)
                    html.append('</div>')
                    html.append('<div>')
                    html.append(_desc_html('mle_consistency'))
                    v = self.plot_mle_consistency_scatter(cond_ctrl, cond_treat)
                    if v:
                        html.append(v)
                    html.append('</div>')
                    html.append('</div>')

                    # Row 2: Selection scatter + Nine-square
                    html.append('<div class="plot-grid">')
                    html.append('<div>')
                    html.append(_desc_html('mle_selection_scatter'))
                    v = self.plot_mle_selection_scatter(cond_ctrl, cond_treat)
                    if v:
                        html.append(v)
                    html.append('</div>')
                    html.append('<div>')
                    html.append(_desc_html('mle_nine_square'))
                    v = self.plot_mle_nine_square(cond_ctrl, cond_treat)
                    if v:
                        html.append(v)
                    html.append('</div>')
                    html.append('</div>')

            elif len(self.mle_conditions) == 1:
                html.append(_desc_html('mle_beta_density'))
                v = self.plot_mle_beta_density()
                if v:
                    html.append(v)

        if html:
            self.sections.append(('results', 'Analysis Results', '\n'.join(html)))

    # ------------------------------------------------------------------
    # sgRNA data embedding (log2 fold-change)
    # ------------------------------------------------------------------

    def _embed_sgrna_data(self):
        """Embed sgRNA data as log2 fold-change (centered at 0) for line plots."""
        if self.count_df is None:
            return ''
        df = self.normalized_df if self.normalized_df is not None else self.count_df
        genes_data = {}
        for _, row in df.iterrows():
            gene = str(row.get('Gene', 'Unknown'))
            sgrna = str(row.get('sgRNA', ''))
            vals = np.array([float(row[s]) if pd.notna(row[s]) else 0.0 for s in self.sample_names])
            log_vals = np.log2(vals + 1)
            mean_log = log_vals.mean()
            lfc_vals = np.round(log_vals - mean_log, 1).tolist()
            entry = [sgrna] + lfc_vals
            if gene not in genes_data:
                genes_data[gene] = []
            genes_data[gene].append(entry)
        return ('<script>\n'
                'var SGRNA_DATA=' + json.dumps(_json_safe(genes_data), separators=(',', ':')) + ';\n'
                'var SAMPLE_NAMES=' + json.dumps(self.sample_names) + ';\n'
                'var selGenes={};\n'
                '</script>\n')

    def _embed_enrichment_data(self):
        """Embed enrichment ORA results as JSON for download."""
        if not self._enrichment_store:
            return ''
        return ('<script>\nvar ENR_DATA=' + json.dumps(self._enrichment_store) + ';\n</script>\n')

    def _embed_default_genes(self):
        """Embed default gene selections for line plots."""
        if not self._default_genes:
            return ''
        return ('<script>\nvar DEFAULT_GENES=' + json.dumps(self._default_genes) + ';\n</script>\n')

    def _embed_hallmark_data(self):
        """Embed hallmark gene sets for client-side enrichment recomputation."""
        if self.skip_enrichment or not HAS_DECOUPLER:
            return ''
        net = self._fetch_hallmark_net()
        if net is None:
            return ''
        # Build pathway -> gene list mapping
        pathway_genes = {}
        for _, row in net.iterrows():
            src = str(row.get('source', ''))
            tgt = str(row.get('target', ''))
            if src and tgt:
                if src not in pathway_genes:
                    pathway_genes[src] = []
                pathway_genes[src].append(tgt)
        if not pathway_genes:
            return ''
        return ('<script>\nvar HALLMARK_SETS=' + json.dumps(pathway_genes, separators=(',', ':')) + ';\n</script>\n')

    # ------------------------------------------------------------------
    # Downloads section
    # ------------------------------------------------------------------

    def _build_downloads_section(self):
        """Build a downloads section with links to count table and log files."""
        html = []

        # Embed count table as downloadable CSV
        if self.count_df is not None and self.count_table_path:
            ct_name = os.path.basename(self.count_table_path)
            try:
                with open(self.count_table_path, 'r') as f:
                    ct_content = f.read()
                ct_b64 = base64.b64encode(ct_content.encode('utf-8')).decode('ascii')
                html.append('<h3>Count Table</h3>')
                html.append('<p class="plot-desc">The sgRNA count table used for this analysis.</p>')
                html.append('<div class="download-section">')
                html.append('<a class="download-btn" href="data:text/plain;base64,'
                            + ct_b64 + '" download="' + ct_name
                            + '">Download ' + ct_name + '</a>')
                html.append('</div>')
            except Exception:
                logging.warning('Could not embed count table for download')

        # Embed log files
        if self._log_files:
            html.append('<h3>MAGeCK Log Files</h3>')
            html.append('<p class="plot-desc">Log files containing the MAGeCK '
                        'run commands and processing details.</p>')
            html.append('<div class="download-section">')
            for label, path in self._log_files:
                try:
                    with open(path, 'r') as f:
                        log_content = f.read()
                    log_b64 = base64.b64encode(log_content.encode('utf-8')).decode('ascii')
                    html.append('<a class="download-btn" href="data:text/plain;base64,'
                                + log_b64 + '" download="' + label
                                + '">Download ' + label + '</a>')
                except Exception:
                    logging.warning('Could not embed log file: %s', path)
            html.append('</div>')

        # Embed countsummary if available
        if self.countsummary_df is not None:
            cs_name = os.path.basename(self.output_prefix) + '.countsummary.txt'
            try:
                cs_text = self.countsummary_df.to_csv(sep='\t', index=False)
                cs_b64 = base64.b64encode(cs_text.encode('utf-8')).decode('ascii')
                html.append('<h3>Count Summary</h3>')
                html.append('<p class="plot-desc">Mapping statistics from MAGeCK count.</p>')
                html.append('<div class="download-section">')
                html.append('<a class="download-btn" href="data:text/plain;base64,'
                            + cs_b64 + '" download="' + cs_name
                            + '">Download ' + cs_name + '</a>')
                html.append('</div>')
            except Exception:
                logging.warning('Could not embed count summary for download')

        if html:
            self.sections.append(('downloads', 'Downloads', '\n'.join(html)))

    # ------------------------------------------------------------------
    # HTML assembly
    # ------------------------------------------------------------------

    def _build_html(self):
        nav_links = ' '.join(
            '<a href="#' + sid + '">' + stitle + '</a>' for sid, stitle, _ in self.sections)
        section_html = ''
        for sid, stitle, content in self.sections:
            section_html += ('<div class="section" id="' + sid + '">'
                             '<h2>' + stitle + '</h2>' + content + '</div>\n')
        sgrna_script = self._embed_sgrna_data()
        enr_script = self._embed_enrichment_data()
        default_script = self._embed_default_genes()
        hallmark_script = self._embed_hallmark_data()
        html = HTML_TEMPLATE
        html = html.replace('{{TITLE}}', self.title)
        html = html.replace('{{VERSION}}', __version__)
        html = html.replace('{{DATE}}', datetime.now().strftime('%Y-%m-%d %H:%M'))
        html = html.replace('{{PREFIX}}', os.path.basename(self.output_prefix))
        html = html.replace('{{NAV_LINKS}}', nav_links)
        html = html.replace('{{SGRNA_DATA}}', sgrna_script)
        html = html.replace('{{ENR_DATA}}', enr_script)
        html = html.replace('{{DEFAULT_GENES}}', default_script)
        html = html.replace('{{HALLMARK_DATA}}', hallmark_script)
        html = html.replace('{{SECTIONS}}', section_html)
        html = html.replace('{{PLOTLY_CDN}}', _PLOTLY_CDN)
        n_rra = len(self.rra_gene_summaries)
        n_mle = len(self.mle_conditions)
        html = html.replace('{{N_RRA}}', str(n_rra))
        html = html.replace('{{N_MLE}}', str(n_mle))
        n_genes = self.count_df['Gene'].nunique() if self.count_df is not None and 'Gene' in self.count_df.columns else 0
        n_sgrna = len(self.count_df) if self.count_df is not None else 0
        html = html.replace('{{N_GENES}}', '{:,}'.format(n_genes))
        html = html.replace('{{N_SGRNA}}', '{:,}'.format(n_sgrna))
        html = html.replace('{{N_SAMPLES}}', str(len(self.sample_names)))
        return html

    # ------------------------------------------------------------------
    # Main entry
    # ------------------------------------------------------------------

    def generate_report(self, output_file=None):
        logging.info('MAGeCK HTML Report v2 starting...')
        self.detect_outputs()
        if self.count_df is not None:
            self.normalize_counts()
        self.normalize_betas()
        self._build_sample_section()
        if not self.skip_qc:
            self._build_qc_section()
        if not self.skip_results:
            self._build_results_section()
        self._build_downloads_section()
        if not self.sections:
            logging.warning('No data found. Check output prefix and file paths.')
            return
        html = self._build_html()
        if output_file is None:
            output_file = self.output_prefix + '.report.html'
        with open(output_file, 'w') as fh:
            fh.write(html)
        logging.info('HTML report written to: %s', output_file)
        print('Report generated: ' + output_file)


# ---------------------------------------------------------------------------
# Per-table JavaScript IIFE template (using %s substitution for data)
# ---------------------------------------------------------------------------

_TABLE_JS_IIFE = """(function(){
var TID=%(TID)s,GK=%(GK)s,LP=%(LP)s,EK=%(EK)s;
var COLS=%(COLS)s;
var _RAW=%(DATA)s;
var FKEYS=%(FKEYS)s;
var _lpPrimary=LP&&LP.indexOf(",")>=0?LP.split(",")[0]:LP;
/* Reconstruct records from compact column format */
var DATA=[];
for(var _i=0;_i<_RAW.v.length;_i++){var _o={};for(var _j=0;_j<_RAW.k.length;_j++)_o[_RAW.k[_j]]=_RAW.v[_i][_j];DATA.push(_o);}
_RAW=null; /* Free memory */
window["TD_"+TID]=DATA;window["TC_"+TID]=COLS;
var flt=DATA.slice(),sKey="",sAsc=true,pg=0,ps=50;
function fv(v,f){if(v==null||v===undefined)return"";if(f=="s")return v;if(f=="d")return Math.round(Number(v)).toLocaleString();if(f==".2e")return Number(v).toExponential(2);if(f==".4f")return Number(v).toFixed(4);return String(v);}
function ren(){
var tb=document.querySelector("#tbl-"+TID+" tbody");
var s2=ps>0?pg*ps:0,e2=ps>0?Math.min(s2+ps,flt.length):flt.length,h="";
for(var i=s2;i<e2;i++){var r=flt[i];
var sel=window.selGenes&&window.selGenes[_lpPrimary]&&window.selGenes[_lpPrimary].has(r[GK]);
h+='<tr'+(sel?' class="gene-sel"':'')+' data-g="'+r[GK]+'" onclick="tglGene(this,\\''+LP+'\\')">';
for(var j=0;j<COLS.length;j++)h+="<td>"+fv(r[COLS[j].key],COLS[j].fmt)+"</td>";
h+="</tr>";}
tb.innerHTML=h;
document.getElementById("count-"+TID).textContent=flt.length+"/"+DATA.length+" genes";
var tp=ps>0?Math.ceil(flt.length/ps):1;
document.getElementById("pinfo-"+TID).textContent="Page "+(pg+1)+"/"+tp;
}
function af(){
var fv2=parseFloat(document.getElementById("fdr-"+TID).value);if(isNaN(fv2))fv2=1;
var ev2=parseFloat(document.getElementById("eff-"+TID).value);if(isNaN(ev2))ev2=0;
var sv=document.getElementById("search-"+TID).value.toLowerCase();
ps=parseInt(document.getElementById("psize-"+TID).value)||0;
flt=DATA.filter(function(r){
var pf=true;if(fv2<1){pf=false;for(var k=0;k<FKEYS.length;k++){if(r[FKEYS[k]]!=null&&r[FKEYS[k]]<=fv2){pf=true;break;}}}
if(ev2>0&&EK&&r[EK]!=null){if(Math.abs(r[EK])<ev2)pf=false;}
var ps2=!sv||!r[GK]||r[GK].toLowerCase().indexOf(sv)>=0;
return pf&&ps2;});
if(sKey)ds();pg=0;ren();}
function ds(){flt.sort(function(a,b){var va=a[sKey],vb=b[sKey];if(va==null)return 1;if(vb==null)return-1;var d=typeof va==="number"?va-vb:String(va).localeCompare(String(vb));return sAsc?d:-d;});}
/* Expose for external sync from volcano */
window["AF_"+TID]=af;
window["FLT_"+TID]=function(){return flt;};
window["SYNC_"+TID]=function(fdr,eff){
 document.getElementById("fdr-"+TID).value=fdr;
 if(eff!==undefined)document.getElementById("eff-"+TID).value=eff;
 af();
};
document.getElementById("fdr-"+TID).addEventListener("input",af);
document.getElementById("eff-"+TID).addEventListener("input",af);
document.getElementById("search-"+TID).addEventListener("input",af);
document.getElementById("psize-"+TID).addEventListener("change",af);
document.getElementById("prev-"+TID).addEventListener("click",function(){if(pg>0){pg--;ren();}});
document.getElementById("next-"+TID).addEventListener("click",function(){var tp=ps>0?Math.ceil(flt.length/ps):1;if(pg<tp-1){pg++;ren();}});
var ths=document.querySelectorAll("#tbl-"+TID+" thead th");
for(var i=0;i<ths.length;i++)(function(idx){ths[idx].addEventListener("click",function(){var k=COLS[idx].key;if(sKey===k)sAsc=!sAsc;else{sKey=k;sAsc=true;}ds();pg=0;ren();for(var x=0;x<ths.length;x++)ths[x].textContent=ths[x].textContent.replace(/ [\\u25B2\\u25BC]/,"");this.textContent+=sAsc?"\\u0020\\u25B2":"\\u0020\\u25BC";});})(i);
ren();
})();"""


# ---------------------------------------------------------------------------
# HTML Template
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{{TITLE}}</title>
<script src="{{PLOTLY_CDN}}" charset="utf-8"></script>
<style>
*{box-sizing:border-box;margin:0;padding:0;}
body{font-family:'Segoe UI',Tahoma,Geneva,Verdana,sans-serif;color:#333;background:#f4f6f8;line-height:1.6;}
/* Header */
.header{background:linear-gradient(135deg,#1a2332 0%,#2c3e50 40%,#2980b9 100%);color:white;padding:28px 40px;position:relative;overflow:hidden;}
.header::after{content:'';position:absolute;top:0;right:0;width:300px;height:100%;background:url("data:image/svg+xml,%3Csvg width='300' height='100' xmlns='http://www.w3.org/2000/svg'%3E%3Ctext x='50%25' y='50%25' font-family='monospace' font-size='60' fill='rgba(255,255,255,0.04)' text-anchor='middle' dominant-baseline='middle'%3ECRISPR%3C/text%3E%3C/svg%3E") no-repeat center center;}
.header h1{font-size:22px;margin-bottom:3px;font-weight:600;letter-spacing:0.5px;}
.header .subtitle{color:#a0b4c8;font-size:12px;letter-spacing:0.3px;}
.header-stats{display:flex;gap:25px;margin-top:12px;flex-wrap:wrap;}
.header-stat{background:rgba(255,255,255,0.1);border-radius:6px;padding:6px 14px;font-size:12px;backdrop-filter:blur(4px);}
.header-stat strong{font-size:16px;display:block;color:#fff;margin-bottom:1px;}
.header-stat span{color:#a0b4c8;font-size:10px;text-transform:uppercase;letter-spacing:0.5px;}
/* Navigation */
.nav{background:#222d3a;padding:8px 40px;position:sticky;top:0;z-index:100;box-shadow:0 2px 8px rgba(0,0,0,0.25);white-space:nowrap;overflow-x:auto;}
.nav a{color:#c8d6e5;text-decoration:none;margin-right:22px;font-size:13px;padding:5px 0;display:inline-block;border-bottom:2px solid transparent;transition:all 0.2s;}
.nav a:hover{color:#3498db;border-bottom-color:#3498db;}
/* Container — wider */
.container{max-width:1600px;margin:18px auto;padding:0 24px;}
/* Sections */
.section{margin:18px 0;padding:22px 28px;background:white;border:1px solid #e2e8f0;border-radius:8px;box-shadow:0 1px 4px rgba(0,0,0,0.06);}
.section h2{color:#1a2332;border-bottom:3px solid #2980b9;padding-bottom:8px;margin-bottom:18px;font-size:18px;font-weight:600;}
.section h3{color:#2c3e50;margin-top:22px;margin-bottom:8px;font-size:15px;font-weight:600;}
.section h4{color:#445;margin-top:18px;margin-bottom:8px;font-size:14px;border-left:4px solid #2980b9;padding-left:10px;font-weight:600;}
.section h5{color:#556;margin-top:14px;margin-bottom:6px;font-size:13px;font-weight:600;}
/* Plot containers */
.plotly-chart{margin:8px 0;}
.plot-grid{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin:10px 0;}
.plot-grid>div{min-width:0;}
.plot-desc{color:#5a6a7a;font-size:11.5px;margin:4px 0 8px;padding:6px 10px;border-left:3px solid #cbd5e0;background:#f8fafc;line-height:1.5;border-radius:0 4px 4px 0;}
/* Summary and info */
.summary-box{background:#f0f4f8;border-left:4px solid #2980b9;padding:10px 16px;margin:10px 0;font-size:13px;border-radius:0 4px 4px 0;}
.report-meta{background:#f8fafc;border:1px solid #e2e8f0;border-radius:6px;padding:14px 18px;margin-bottom:16px;}
.report-meta .info-table td:first-child{color:#2980b9;font-weight:600;min-width:150px;}
.metadata{color:#7f8c8d;font-size:12px;margin-bottom:8px;font-style:italic;}
.info-table{border-collapse:collapse;width:auto;margin:6px 0;font-size:13px;}
.info-table td{padding:3px 14px 3px 0;border:none;}
.info-table td:first-child{font-weight:600;color:#445;}
.simple-table{border-collapse:collapse;width:auto;margin:10px 0;font-size:12px;}
.simple-table th{background:#2c3e50;color:white;padding:7px 11px;text-align:left;font-weight:500;}
.simple-table td{padding:5px 11px;border-bottom:1px solid #eee;}
.simple-table tr:nth-child(even){background:#f8fafc;}
/* Dynamic table styles */
.table-container{margin:12px 0;}
.table-controls{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin-bottom:6px;font-size:12px;}
.table-controls input[type="number"]{width:65px;padding:4px 6px;border:1px solid #ccc;border-radius:3px;}
.table-controls input[type="text"]{width:170px;padding:4px 6px;border:1px solid #ccc;border-radius:3px;}
.table-controls select{padding:4px 6px;border:1px solid #ccc;border-radius:3px;}
.table-controls span{color:#7f8c8d;font-size:11px;}
.table-scroll{max-height:500px;overflow:auto;border:1px solid #ddd;border-radius:4px;}
.table-scroll table{border-collapse:collapse;width:100%;font-size:11.5px;}
.table-scroll th{background:#2c3e50;color:white;padding:7px 9px;text-align:left;cursor:pointer;user-select:none;white-space:nowrap;position:sticky;top:0;z-index:1;font-weight:500;font-size:11px;}
.table-scroll th:hover{background:#3d5a73;}
.table-scroll td{padding:5px 9px;border-bottom:1px solid #eee;white-space:nowrap;}
.table-scroll tr:nth-child(even){background:#f8fafc;}
.table-scroll tr:hover{background:#eef5fd;}
.table-scroll tr.gene-sel{background:#d4edda !important;}
.table-scroll tr.gene-sel:hover{background:#c3e6cb !important;}
.table-pag{display:flex;gap:8px;align-items:center;margin-top:6px;font-size:12px;}
.table-pag button{padding:3px 10px;cursor:pointer;border:1px solid #ccc;border-radius:3px;background:white;font-size:12px;}
.table-pag button:hover{background:#f0f0f0;}
/* Volcano controls */
.vol-ctrl{display:flex;gap:12px;align-items:center;padding:6px 10px;background:#f8fafc;border:1px solid #e2e8f0;border-radius:4px;margin:6px 0 4px;font-size:12px;color:#445;}
.vol-ctrl input[type="number"]{width:65px;padding:3px 6px;border:1px solid #ccc;border-radius:3px;font-size:12px;}
/* Downloads */
.download-section{display:flex;flex-wrap:wrap;gap:10px;margin:12px 0;}
.download-btn{display:inline-block;padding:8px 16px;background:#2980b9;color:white;border:none;border-radius:5px;font-size:12px;cursor:pointer;font-weight:500;transition:background 0.2s;}
.download-btn:hover{background:#1a6da0;}
.download-btn-sm{display:inline-block;padding:4px 10px;background:#e2e8f0;color:#445;border:1px solid #cbd5e0;border-radius:3px;font-size:10.5px;cursor:pointer;transition:all 0.2s;}
.download-btn-sm:hover{background:#cbd5e0;color:#222;}
/* Footer */
.footer{text-align:center;color:#7f8c8d;font-size:11px;padding:18px 20px;border-top:1px solid #e2e8f0;margin-top:25px;background:white;line-height:1.8;}
.footer a{color:#2980b9;text-decoration:none;}
.footer .footer-line{margin-bottom:2px;}
/* Print & responsive */
@media print{.nav{display:none;}.section{break-inside:avoid;box-shadow:none;border:1px solid #ccc;}body{background:white;}.plot-grid{grid-template-columns:1fr;}}
@media(max-width:900px){.plot-grid{grid-template-columns:1fr;}.container{padding:0 8px;}.section{padding:14px;}.table-controls{flex-direction:column;align-items:flex-start;}.header-stats{flex-direction:column;gap:8px;}}
</style>
</head>
<body>
<div class="header">
 <h1>{{TITLE}}</h1>
 <div class="subtitle">MAGeCK v{{VERSION}} | Generated {{DATE}} | Prefix: {{PREFIX}}</div>
 <div class="header-stats">
  <div class="header-stat"><strong>{{N_SGRNA}}</strong><span>sgRNAs</span></div>
  <div class="header-stat"><strong>{{N_GENES}}</strong><span>Genes</span></div>
  <div class="header-stat"><strong>{{N_SAMPLES}}</strong><span>Samples</span></div>
  <div class="header-stat"><strong>{{N_RRA}}</strong><span>RRA Comparisons</span></div>
  <div class="header-stat"><strong>{{N_MLE}}</strong><span>MLE Conditions</span></div>
 </div>
</div>
<nav class="nav">{{NAV_LINKS}}</nav>
<div class="container">
{{SGRNA_DATA}}
{{ENR_DATA}}
{{DEFAULT_GENES}}
{{HALLMARK_DATA}}
{{SECTIONS}}
</div>
<div class="footer">
 <div class="footer-line">Package used: <strong>MAGeCK v{{VERSION}}</strong></div>
 <div class="footer-line">Generated by <strong>CRISPR Functional Genomics, Karolinska Institute</strong></div>
</div>
<script>
var CONFIG=''' + _PLOTLY_CONFIG + r''';
/* Interactive volcano redraw */
function drawVolcano(divId){
 var D=window["VD_"+divId];if(!D)return;
 var fc=parseFloat(document.getElementById("vfdr-"+divId).value);if(isNaN(fc))fc=0.05;
 var ec=parseFloat(document.getElementById("veff-"+divId).value)||0;
 var ht="Gene: %{text}<br>"+D.xl+": %{x:.3f}<br>"+D.yl+": %{y:.2f}";
 var ns={x:[],y:[],text:[],type:"scattergl",mode:"markers",marker:{color:"rgba(187,187,187,0.5)",size:4},name:"NS",hovertemplate:ht};
 var sd={x:[],y:[],text:[],type:"scattergl",mode:"markers",marker:{color:"#3498DB",size:7},name:"Depleted (\u2212)",hovertemplate:ht};
 var se={x:[],y:[],text:[],type:"scattergl",mode:"markers",marker:{color:"#E74C3C",size:7},name:"Enriched (+)",hovertemplate:ht};
 var depGenes=[],enrGenes=[];
 for(var i=0;i<D.g.length;i++){
  var fdr=D.lfc[i]<0?D.nf[i]:D.pf[i];
  if(fdr<fc&&Math.abs(D.lfc[i])>=ec){if(D.lfc[i]<0){sd.x.push(D.lfc[i]);sd.y.push(D.y[i]);sd.text.push(D.g[i]);depGenes.push(D.g[i]);}else{se.x.push(D.lfc[i]);se.y.push(D.y[i]);se.text.push(D.g[i]);enrGenes.push(D.g[i]);}}
  else{ns.x.push(D.lfc[i]);ns.y.push(D.y[i]);ns.text.push(D.g[i]);}
 }
 var shapes=[{type:"line",x0:0,x1:1,xref:"paper",y0:-Math.log10(Math.max(fc,1e-10)),y1:-Math.log10(Math.max(fc,1e-10)),line:{dash:"dash",color:"grey",width:1}},{type:"line",y0:0,y1:1,yref:"paper",x0:0,x1:0,line:{dash:"solid",color:"rgba(150,150,150,0.3)",width:1}}];
 if(ec>0){shapes.push({type:"line",y0:0,y1:1,yref:"paper",x0:-ec,x1:-ec,line:{dash:"dot",color:"#aaa",width:1}});shapes.push({type:"line",y0:0,y1:1,yref:"paper",x0:ec,x1:ec,line:{dash:"dot",color:"#aaa",width:1}});}
 /* Top 20 enriched/depleted gene labels */
 var anns=[];
 if(D.top_enr){for(var j=0;j<D.top_enr.length;j++){var k=D.top_enr[j];anns.push({x:D.lfc[k],y:D.y[k],text:D.g[k],showarrow:true,arrowhead:0,arrowcolor:"rgba(231,76,60,0.4)",ax:18+(j%3)*12,ay:-8-(Math.floor(j/3)%7)*11,font:{size:7,color:"#C0392B"},bgcolor:"rgba(255,255,255,0.7)"});}}
 if(D.top_dep){for(var j=0;j<D.top_dep.length;j++){var k=D.top_dep[j];anns.push({x:D.lfc[k],y:D.y[k],text:D.g[k],showarrow:true,arrowhead:0,arrowcolor:"rgba(52,152,219,0.4)",ax:-18-(j%3)*12,ay:-8-(Math.floor(j/3)%7)*11,font:{size:7,color:"#2166AC"},bgcolor:"rgba(255,255,255,0.7)"});}}
 Plotly.react(divId,[ns,sd,se],{title:{text:D.title,font:{size:13}},xaxis:{title:D.xl,automargin:true,tickfont:{size:10}},yaxis:{title:D.yl,automargin:true,tickfont:{size:10}},legend:{font:{size:10}},height:450,margin:{t:45,b:45,l:55,r:15},shapes:shapes,annotations:anns},CONFIG);
 /* Table has its own independent FDR filter (no sync from volcano) */
 /* Update linked enrichment dot plots */
 if(D.ediv&&typeof HALLMARK_SETS!=="undefined"){updateEnrichment(D.ediv,depGenes,D.g.length,"depleted");}
 if(D.eenr&&typeof HALLMARK_SETS!=="undefined"){updateEnrichment(D.eenr,enrGenes,D.g.length,"enriched");}
}
/* Client-side enrichment recomputation */
function updateEnrichment(divId,geneList,totalGenes,direction){
 if(!document.getElementById(divId))return;
 if(!geneList||geneList.length===0){Plotly.react(divId,[],{title:{text:"No significant genes at current cutoffs"},height:270,margin:{t:40,b:20,l:10,r:10}},CONFIG);return;}
 var gs=new Set(geneList),nSig=geneList.length,results=[];
 var pathways=Object.keys(HALLMARK_SETS);
 for(var p=0;p<pathways.length;p++){
  var pw=pathways[p],members=HALLMARK_SETS[pw],K=members.length,overlap=0;
  for(var m=0;m<members.length;m++){if(gs.has(members[m]))overlap++;}
  if(overlap===0)continue;
  var expected=nSig*K/totalGenes;
  var fold=overlap/Math.max(expected,0.001);
  /* Approx -log10 p-value via normal approx to hypergeometric */
  var variance=nSig*(K/totalGenes)*(1-K/totalGenes)*((totalGenes-nSig)/(totalGenes-1));
  var z=variance>0?(overlap-expected)/Math.sqrt(variance):0;
  var nlp=z>0?Math.min(z*z/4.6,20):0;  /* rough -log10(p) approx */
  results.push({name:pw,overlap:overlap,fold:fold,nlp:nlp});
 }
 results.sort(function(a,b){return b.fold-a.fold;});
 var show=results.slice(0,15).reverse();
 if(show.length===0){Plotly.react(divId,[],{title:{text:"No pathway overlap at current cutoffs"},height:270,margin:{t:40,b:20,l:10,r:10}},CONFIG);return;}
 var maxNlp=Math.max.apply(null,show.map(function(r){return r.nlp;}));
 var minNlp=Math.min.apply(null,show.map(function(r){return r.nlp;}));
 var isCold=direction==="depleted";
 var lineCol=isCold?"#2166AC":"#B2182B";
 var traces=[{type:"scatter",mode:"markers",
  x:show.map(function(r){return r.fold;}),
  y:show.map(function(r){return r.name.length>42?r.name.substring(0,41)+"\u2026":r.name;}),
  marker:{size:show.map(function(r){return Math.max(8,r.nlp/(maxNlp||1)*28);}),
   color:show.map(function(r){return r.fold;}),
   colorscale:isCold?[[0,"#D1E5F0"],[0.5,"#4393C3"],[1,"#053061"]]:[[0,"#FDDBC7"],[0.5,"#D6604D"],[1,"#67001F"]],
   colorbar:{title:{text:"Fold Enrich.",font:{size:9}},len:0.4,thickness:12,tickfont:{size:8},x:0.82,xpad:3,y:0.8,yanchor:"middle"},
   line:{width:1.2,color:lineCol},opacity:0.85},
  text:show.map(function(r){return"Overlap: "+r.overlap+", -log10(p)\u2248"+r.nlp.toFixed(1);}),
  hovertemplate:"%{y}<br>Fold Enrichment: %{x:.2f}<br>%{text}",showlegend:false}];
 /* Size legend reference dots — placed on right side near colorbar */
 var refV=[];
 if(maxNlp>2){refV=[Math.round(minNlp*10)/10,Math.round((minNlp+(maxNlp-minNlp)*0.5)*10)/10,Math.round(maxNlp*10)/10];}
 else if(maxNlp>0){refV=[Math.round(minNlp*10)/10,Math.round(maxNlp*10)/10];}
 var uniqR=[],seenR={};refV.forEach(function(v){if(!seenR[v]){seenR[v]=1;uniqR.push(v);}});refV=uniqR;
 if(refV.length>0){
  traces.push({type:"scatter",mode:"markers+text",
   x:refV.map(function(){return 0;}),y:refV.map(function(v,i){return i;}),
   marker:{size:refV.map(function(v){return Math.max(8,v/(maxNlp||1)*28);}),color:lineCol,opacity:0.5,line:{width:1,color:lineCol}},
   text:refV.map(function(v){return v.toFixed(1);}),textposition:"middle right",textfont:{size:8,color:"#555"},
   hoverinfo:"skip",showlegend:false,xaxis:"x2",yaxis:"y2"});
 }
 var pH=Math.max(270,show.length*24+100);
 var anns2=[];
 if(refV.length>0){anns2.push({text:"<b>-log\u2081\u2080(p)</b>",xref:"x2 domain",yref:"y2 domain",x:0.5,y:1.15,showarrow:false,font:{size:8,color:"#555"}});}
 var ttl=(direction==="depleted"?"Depleted":"Enriched")+" Genes: Hallmark Pathways (n="+geneList.length+")";
 Plotly.react(divId,traces,{title:{text:ttl,font:{size:12},y:0.98,yanchor:"top",pad:{t:10}},
  xaxis:{title:"Fold Enrichment",automargin:true,tickfont:{size:9},domain:[0,0.78]},
  yaxis:{automargin:true,tickfont:{size:9},domain:[0,1.0]},
  xaxis2:{anchor:"y2",domain:[0.88,0.98],showgrid:false,zeroline:false,showticklabels:false,fixedrange:true},
  yaxis2:{anchor:"x2",domain:[0.05,0.45],showgrid:false,zeroline:false,showticklabels:false,fixedrange:true},
  height:pH,margin:{l:10,t:50,b:30,r:80},annotations:anns2},CONFIG);
}
/* Interactive selection scatter */
function drawSelection(divId){
 var D=window["SD_"+divId];if(!D)return;
 var cut=parseFloat(document.getElementById("scut-"+divId).value)||0;
 var ht="Gene: %{text}<br>"+D.ctrl+": %{x:.3f}<br>"+D.treat+": %{y:.3f}";
 var grey={x:[],y:[],text:[]},pink={x:[],y:[],text:[]},blue={x:[],y:[],text:[]};
 for(var i=0;i<D.g.length;i++){
  var diff=D.vy[i]-D.vx[i];
  if(diff>cut){pink.x.push(D.vx[i]);pink.y.push(D.vy[i]);pink.text.push(D.g[i]);}
  else if(diff<-cut){blue.x.push(D.vx[i]);blue.y.push(D.vy[i]);blue.text.push(D.g[i]);}
  else{grey.x.push(D.vx[i]);grey.y.push(D.vy[i]);grey.text.push(D.g[i]);}
 }
 var traces=[];
 if(grey.x.length)traces.push({type:"scatter",mode:"markers",x:grey.x,y:grey.y,text:grey.text,marker:{color:"rgba(180,180,180,0.4)",size:4},hovertemplate:ht,name:"Others ("+grey.x.length+")"});
 if(pink.x.length)traces.push({type:"scatter",mode:"markers",x:pink.x,y:pink.y,text:pink.text,marker:{color:"#FF69B4",size:6,opacity:0.7},hovertemplate:ht,name:"Enriched + ("+pink.x.length+")"});
 if(blue.x.length)traces.push({type:"scatter",mode:"markers",x:blue.x,y:blue.y,text:blue.text,marker:{color:"#6495ED",size:6,opacity:0.7},hovertemplate:ht,name:"Depleted \u2212 ("+blue.x.length+")"});
 var allx=D.vx,ally=D.vy;
 var mn=Math.min(Math.min.apply(null,allx),Math.min.apply(null,ally))-0.2;
 var mx=Math.max(Math.max.apply(null,allx),Math.max.apply(null,ally))+0.2;
 var shapes=[{type:"line",x0:mn,x1:mx,y0:mn,y1:mx,line:{dash:"dot",color:"#333",width:1}}];
 if(cut>0){
  shapes.push({type:"line",x0:mn,x1:mx,y0:mn+cut,y1:mx+cut,line:{dash:"dash",color:"rgba(255,105,180,0.5)",width:1}});
  shapes.push({type:"line",x0:mn,x1:mx,y0:mn-cut,y1:mx-cut,line:{dash:"dash",color:"rgba(100,149,237,0.5)",width:1}});
 }
 Plotly.react(divId,traces,{title:{text:"Selection: "+D.treat+" vs "+D.ctrl,font:{size:13}},
  xaxis:{title:D.ctrl+" (beta)",automargin:true,tickfont:{size:10}},
  yaxis:{title:D.treat+" (beta)",automargin:true,tickfont:{size:10}},
  legend:{font:{size:10}},height:360,margin:{t:45,b:40,l:50,r:10},shapes:shapes},CONFIG);
}
/* FLUTE nine-square treatment gene classification */
function drawNineSquare(divId){
 var D=window["ND_"+divId];if(!D)return;
 var cut=parseFloat(document.getElementById("ncut-"+divId).value)||0;
 var ht="Gene: %{text}<br>"+D.ctrl+" (\u03b2): %{x:.3f}<br>"+D.treat+" (\u03b2): %{y:.3f}";
 /* 4 treatment groups + others */
 var g1={x:[],y:[],text:[]}, /* green: strong neg ctrl, weak treat -> drug targets */
     g2={x:[],y:[],text:[]}, /* orange: weak ctrl, strong pos treat -> resistance */
     g3={x:[],y:[],text:[]}, /* blue: strong pos ctrl, weak treat */
     g4={x:[],y:[],text:[]}, /* purple: weak ctrl, strong neg treat -> synth lethal */
     grey={x:[],y:[],text:[]};
 for(var i=0;i<D.g.length;i++){
  var x=D.vx[i],y=D.vy[i],ax=Math.abs(x),ay=Math.abs(y);
  var xStrong=ax>=cut,yStrong=ay>=cut;
  if(cut>0 && x<-cut && !yStrong){g1.x.push(x);g1.y.push(y);g1.text.push(D.g[i]);}
  else if(cut>0 && !xStrong && y>cut){g2.x.push(x);g2.y.push(y);g2.text.push(D.g[i]);}
  else if(cut>0 && x>cut && !yStrong){g3.x.push(x);g3.y.push(y);g3.text.push(D.g[i]);}
  else if(cut>0 && !xStrong && y<-cut){g4.x.push(x);g4.y.push(y);g4.text.push(D.g[i]);}
  else{grey.x.push(x);grey.y.push(y);grey.text.push(D.g[i]);}
 }
 var traces=[];
 if(grey.x.length)traces.push({type:"scattergl",mode:"markers",x:grey.x,y:grey.y,text:grey.text,marker:{color:"rgba(180,180,180,0.3)",size:4},hovertemplate:ht,name:"Others ("+grey.x.length+")"});
 if(g1.x.length)traces.push({type:"scattergl",mode:"markers",x:g1.x,y:g1.y,text:g1.text,marker:{color:"#27AE60",size:6,opacity:0.8},hovertemplate:ht,name:"G1: Drug Targets ("+g1.x.length+")"});
 if(g2.x.length)traces.push({type:"scattergl",mode:"markers",x:g2.x,y:g2.y,text:g2.text,marker:{color:"#E67E22",size:6,opacity:0.8},hovertemplate:ht,name:"G2: Resistance ("+g2.x.length+")"});
 if(g3.x.length)traces.push({type:"scattergl",mode:"markers",x:g3.x,y:g3.y,text:g3.text,marker:{color:"#3498DB",size:6,opacity:0.8},hovertemplate:ht,name:"G3: Pos. Ctrl Only ("+g3.x.length+")"});
 if(g4.x.length)traces.push({type:"scattergl",mode:"markers",x:g4.x,y:g4.y,text:g4.text,marker:{color:"#8E44AD",size:6,opacity:0.8},hovertemplate:ht,name:"G4: Synth. Lethal ("+g4.x.length+")"});
 var allx=D.vx,ally=D.vy;
 var mn=Math.min(Math.min.apply(null,allx),Math.min.apply(null,ally))-0.3;
 var mx=Math.max(Math.max.apply(null,allx),Math.max.apply(null,ally))+0.3;
 /* Grid: diagonal y=x, horizontal y=+/-cut, vertical x=+/-cut */
 var shapes=[{type:"line",x0:mn,x1:mx,y0:mn,y1:mx,line:{dash:"dot",color:"#555",width:1}}];
 if(cut>0){
  /* Vertical lines at x = +/- cut */
  shapes.push({type:"line",y0:mn,y1:mx,x0:-cut,x1:-cut,line:{dash:"dash",color:"rgba(100,100,100,0.4)",width:1}});
  shapes.push({type:"line",y0:mn,y1:mx,x0:cut,x1:cut,line:{dash:"dash",color:"rgba(100,100,100,0.4)",width:1}});
  /* Horizontal lines at y = +/- cut */
  shapes.push({type:"line",x0:mn,x1:mx,y0:-cut,y1:-cut,line:{dash:"dash",color:"rgba(100,100,100,0.4)",width:1}});
  shapes.push({type:"line",x0:mn,x1:mx,y0:cut,y1:cut,line:{dash:"dash",color:"rgba(100,100,100,0.4)",width:1}});
  /* Colored region highlights (transparent rectangles) */
  shapes.push({type:"rect",x0:mn,x1:-cut,y0:-cut,y1:cut,fillcolor:"rgba(39,174,96,0.06)",line:{width:0}});  /* G1 green */
  shapes.push({type:"rect",x0:-cut,x1:cut,y0:cut,y1:mx,fillcolor:"rgba(230,126,34,0.06)",line:{width:0}});  /* G2 orange */
  shapes.push({type:"rect",x0:cut,x1:mx,y0:-cut,y1:cut,fillcolor:"rgba(52,152,219,0.06)",line:{width:0}});  /* G3 blue */
  shapes.push({type:"rect",x0:-cut,x1:cut,y0:mn,y1:-cut,fillcolor:"rgba(142,68,173,0.06)",line:{width:0}}); /* G4 purple */
 }
 /* Group labels and gene counts */
 var anns=[];
 if(cut>0){
  anns.push({x:mn+0.15,y:0,text:"<b>G1</b> ("+g1.x.length+")",showarrow:false,font:{size:10,color:"#27AE60"},bgcolor:"rgba(255,255,255,0.85)"});
  anns.push({x:0,y:mx-0.15,text:"<b>G2</b> ("+g2.x.length+")",showarrow:false,font:{size:10,color:"#E67E22"},bgcolor:"rgba(255,255,255,0.85)"});
  anns.push({x:mx-0.15,y:0,text:"<b>G3</b> ("+g3.x.length+")",showarrow:false,font:{size:10,color:"#3498DB"},bgcolor:"rgba(255,255,255,0.85)"});
  anns.push({x:0,y:mn+0.15,text:"<b>G4</b> ("+g4.x.length+")",showarrow:false,font:{size:10,color:"#8E44AD"},bgcolor:"rgba(255,255,255,0.85)"});
 }
 /* Label top 5 genes in each group */
 function labelTop(arr,col,axOff){
  if(arr.x.length===0)return;
  var sc=[];for(var i=0;i<arr.x.length;i++)sc.push({i:i,s:Math.abs(arr.x[i])+Math.abs(arr.y[i])});
  sc.sort(function(a,b){return b.s-a.s;});
  for(var j=0;j<Math.min(5,sc.length);j++){var idx=sc[j].i;
   anns.push({x:arr.x[idx],y:arr.y[idx],text:arr.text[idx],showarrow:true,arrowhead:0,arrowcolor:col,ax:axOff,ay:-8-(j%3)*12,font:{size:7.5,color:col},bgcolor:"rgba(255,255,255,0.7)"});
  }
 }
 labelTop(g1,"#1E8449",-22);labelTop(g2,"#D35400",0);labelTop(g3,"#2471A3",22);labelTop(g4,"#6C3483",0);
 Plotly.react(divId,traces,{title:{text:"Nine-Square Treatment Classification: "+D.treat+" vs "+D.ctrl,font:{size:13}},
  xaxis:{title:D.ctrl+" (normalized \u03b2)",automargin:true,tickfont:{size:10},zeroline:true,zerolinecolor:"rgba(0,0,0,0.15)",zerolinewidth:1},
  yaxis:{title:D.treat+" (normalized \u03b2)",automargin:true,tickfont:{size:10},zeroline:true,zerolinecolor:"rgba(0,0,0,0.15)",zerolinewidth:1},
  legend:{font:{size:9},x:1.02,y:1,xanchor:"left"},height:480,margin:{t:45,b:45,l:55,r:160},shapes:shapes,annotations:anns},CONFIG);
}
/* Gene selection -> line plot */
function tglGene(tr,lpId){
 if(!lpId)return;
 var gene=tr.getAttribute("data-g");
 var ids=lpId.split(","),primary=ids[0];
 if(!window.selGenes[primary])window.selGenes[primary]=new Set();
 var s=window.selGenes[primary];
 if(s.has(gene)){s.delete(gene);tr.classList.remove("gene-sel");}
 else{if(s.size>=50){alert("Max 50 genes");return;}s.add(gene);tr.classList.add("gene-sel");}
 for(var i=0;i<ids.length;i++){window.selGenes[ids[i]]=s;updLP(ids[i]);}
}
function updLP(lpId){
 var gs=window.selGenes[lpId];
 if(!gs||gs.size===0){Plotly.react(lpId,[],{title:{text:"Click genes in the table to show sgRNA log2 FC"},xaxis:{title:"Sample"},yaxis:{title:"log2 Fold Change",zeroline:true},height:380});return;}
 if(typeof SGRNA_DATA==="undefined")return;
 var traces=[],colors=["#E74C3C","#3498DB","#2ECC71","#9B59B6","#F39C12","#1ABC9C","#E67E22","#34495E","#16A085","#C0392B"],ci=0;
 gs.forEach(function(gene){
  var sgs=SGRNA_DATA[gene];if(!sgs)return;
  for(var k=0;k<sgs.length;k++){
   var sg=sgs[k],y=[];
   for(var m=1;m<sg.length;m++)y.push(sg[m]);
   traces.push({x:SAMPLE_NAMES,y:y,mode:"lines+markers",name:gene+" - "+sg[0],line:{color:colors[ci%colors.length]},marker:{size:5},legendgroup:gene,
    hovertemplate:"Gene: "+gene+"<br>sgRNA: "+sg[0]+"<br>Sample: %{x}<br>LFC: %{y:.2f}"});
  }
  ci++;
 });
 Plotly.react(lpId,traces,{title:{text:"sgRNA log2 Fold Change for Selected Genes"},xaxis:{title:"Sample"},yaxis:{title:"log2 Fold Change",zeroline:true,zerolinecolor:"rgba(0,0,0,0.3)",zerolinewidth:1.5},height:380,showlegend:true,legend:{font:{size:10}}});
}
/* Clear all gene selections for a line plot */
function clearLP(lpId){
 var ids=lpId.split(","),primary=ids[0];
 if(!window.selGenes[primary])return;
 window.selGenes[primary].clear();
 document.querySelectorAll("[data-g]").forEach(function(tr){tr.classList.remove("gene-sel");});
 for(var i=0;i<ids.length;i++){window.selGenes[ids[i]]=window.selGenes[primary];updLP(ids[i]);}
}
/* Auto-select default genes on page load */
function initDefaults(){
 if(typeof DEFAULT_GENES==="undefined")return;
 /* Build shared Set for LPs that have identical default gene lists */
 var sharedSets={};
 Object.keys(DEFAULT_GENES).forEach(function(lpId){
  var genes=DEFAULT_GENES[lpId];
  if(!genes||genes.length===0)return;
  var key=genes.slice().sort().join(",");
  if(!sharedSets[key]){sharedSets[key]=new Set();genes.forEach(function(g){sharedSets[key].add(g);});}
  window.selGenes[lpId]=sharedSets[key];
  updLP(lpId);
 });
 /* Highlight default genes in visible table rows */
 setTimeout(function(){
  document.querySelectorAll("[data-g]").forEach(function(tr){
   var g=tr.getAttribute("data-g");
   Object.keys(DEFAULT_GENES).forEach(function(lpId){
    if(window.selGenes[lpId]&&window.selGenes[lpId].has(g))
     tr.classList.add("gene-sel");
   });
  });
 },200);
}
initDefaults();
/* CSV download from table data (all rows) */
function downloadCSV(tableId,filename){
 var data=window["TD_"+tableId],cols=window["TC_"+tableId];
 if(!data||!cols)return;
 var lines=[cols.map(function(c){return'"'+c.label+'"';}).join(",")];
 data.forEach(function(r){lines.push(cols.map(function(c){var v=r[c.key];if(v==null)return"";if(typeof v==="string")return'"'+v.replace(/"/g,'""')+'"';return String(v);}).join(","));});
 var blob=new Blob([lines.join("\n")],{type:"text/csv;charset=utf-8;"});
 var url=URL.createObjectURL(blob);
 var a=document.createElement("a");a.href=url;a.download=filename;document.body.appendChild(a);a.click();document.body.removeChild(a);URL.revokeObjectURL(url);
}
/* CSV download from filtered table data */
function downloadFilteredCSV(tableId,filename){
 var fn=window["FLT_"+tableId];
 if(!fn)return downloadCSV(tableId,filename);
 var data=fn(),cols=window["TC_"+tableId];
 if(!data||!cols)return;
 var lines=[cols.map(function(c){return'"'+c.label+'"';}).join(",")];
 data.forEach(function(r){lines.push(cols.map(function(c){var v=r[c.key];if(v==null)return"";if(typeof v==="string")return'"'+v.replace(/"/g,'""')+'"';return String(v);}).join(","));});
 var blob=new Blob([lines.join("\n")],{type:"text/csv;charset=utf-8;"});
 var url=URL.createObjectURL(blob);
 var a=document.createElement("a");a.href=url;a.download=filename;document.body.appendChild(a);a.click();document.body.removeChild(a);URL.revokeObjectURL(url);
}
/* Enrichment CSV download */
function dlEnr(key,filename){
 if(typeof ENR_DATA==="undefined"||!ENR_DATA[key])return;
 var rows=ENR_DATA[key];
 var lines=["Pathway,ORA_Score,Adjusted_P_Value"];
 rows.forEach(function(r){
  lines.push('"'+r.Pathway+'",'+(r.ORA_Score!=null?r.ORA_Score:"")+","+(r.Adjusted_P_Value!=null?r.Adjusted_P_Value:""));
 });
 var blob=new Blob([lines.join("\n")],{type:"text/csv;charset=utf-8;"});
 var url=URL.createObjectURL(blob);
 var a=document.createElement("a");a.href=url;a.download=filename;document.body.appendChild(a);a.click();document.body.removeChild(a);URL.revokeObjectURL(url);
}
</script>
</body>
</html>'''


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def mageck_report_main(args):
    """Entry point for 'mageck report' subcommand."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
        datefmt='%a, %d %b %Y %H:%M:%S',
        stream=sys.stderr,
    )
    report = MAGeCKHTMLReport(
        output_prefix=getattr(args, 'output_prefix', 'sample1'),
        count_table=getattr(args, 'count_table', None),
        title=getattr(args, 'title', None),
        gene_summary_files=getattr(args, 'gene_summary', None),
        sgrna_summary_files=getattr(args, 'sgrna_summary', None),
        design_matrix=getattr(args, 'design_matrix', None),
        fdr_threshold=getattr(args, 'fdr_threshold', 0.05),
        top_n=getattr(args, 'top_n', 20),
        norm_method=getattr(args, 'norm_method', 'median'),
        skip_qc=getattr(args, 'skip_qc', False),
        skip_results=getattr(args, 'skip_results', False),
        enrichment_top_n=getattr(args, 'enrichment_top_n', 50),
        organism=getattr(args, 'organism', 'human'),
        skip_enrichment=getattr(args, 'skip_enrichment', False),
        enrichment_fdr=getattr(args, 'enrichment_fdr', 0.05),
    )
    report.generate_report(output_file=getattr(args, 'output_file', None))
