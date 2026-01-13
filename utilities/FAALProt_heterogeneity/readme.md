# FAALProt Heterogeneity Analysis

**Complete documentation for the FAALProt heterogeneity CLI (`faalprot_heterogeneity_cli_v2.py`) and its mathematical background.**  

---

## Table of Contents

- [1. Overview](#1-overview)
- [2. Approach](#2-approach)
- [3. Input Data & Metadata](#3-input-data--metadata)
- [4. Software & Environment](#4-software--environment)
- [5. Mathematical Definitions](#5-mathematical-definitions)
  - [5.1 Aligned Sequences](#51-aligned-sequences)
  - [5.2 Pairwise Identity](#52-pairwise-identity)
  - [5.3 Dissimilarity](#53-dissimilarity)
  - [5.4 Word2Vec k-mer Embeddings](#54-word2vec-kmer-embeddings)
  - [5.5 Identity Bins and Pair Sampling](#55-identity-bins-and-pair-sampling)
  - [5.6 Correlation: Identity vs Cosine Similarity](#56-correlation-identity-vs-cosine-similarity)
- [6. Pipeline Stages](#6-pipeline-stages)
  - [6.1 Multiple Sequence Alignment (MAFFT)](#61-multiple-sequence-alignment-mafft)
  - [6.2 Pairwise Identities & Dissimilarity Matrix](#62-pairwise-identities--dissimilarity-matrix)
  - [6.3 Subsampling of Sequences](#63-subsampling-of-sequences)
  - [6.4 MMseqs2 Clustering Across Identity Cutoffs](#64-mmseqs2-clustering-across-identity-cutoffs)
  - [6.5 Word2Vec Grid Search (Part B)](#65-word2vec-grid-search-part-b)
  - [6.6 Visualizations and Figures](#66-visualizations-and-figures)
- [7. Color Policy (Distinct Phylum Colors)](#7-color-policy-distinct-phylum-colors)
- [8. Outputs & Directory Layout](#8-outputs--directory-layout)
- [9. Reproducibility, Performance & Practical Notes](#9-reproducibility-performance--practical-notes)
- [10. How to Run — CLI Version](#10-how-to-run--cli-version)
- [11. FAQ](#11-faq)
- [12. References](#12-references)

---

## 1. Overview

This repository provides a complete, publication-grade pipeline to quantify and visualize **heterogeneity in FAAL proteins** across thousands of sequences. It integrates:

- **MAFFT** for multiple sequence alignment (MSA)
- **Streaming pairwise identity** estimation with deterministic subsampling of the upper triangle
- **UMAP (precomputed metric)** on *dissimilarity* derived from alignment identity
- **Hierarchical clustering** and **dendrograms** with **branch colors by Phylum**
- **MMseqs2 clustering** over fixed identity thresholds (10–90%) and diagnostic plots
- **(Part B) Word2Vec** (k-mer embeddings from aligned sequences without gaps):
  - cosine-vs-identity for **intra**- and **inter-phylum** pairs,
  - UMAP/t-SNE on embeddings,
  - and an embedding-based dendrogram.

The CLI `faalprot_heterogeneity_cli_v2.py` adds a **Word2Vec grid search** with:

- multiple **embedding dimensions** (default: 100 and 390; configurable), and  
- multiple **epoch values** (e.g. 200, 500, 1500, 2500; configurable),

evaluated across **4 subsets of the data**: **25%, 50%, 75% and 100%** of the sequences.  
For each combination `(subset, dimension, epochs)` the pipeline:

1. aligns sequences,
2. generates k-mers,
3. trains Word2Vec with a consistent configuration,
4. computes sequence embeddings,
5. recomputes similarity statistics and all W2V-dependent figures.

All figures are exported as **PNG (900 dpi)**, **TIFF (900 dpi)** and **SVG**, without titles, and with legends **below** the plots.

> **Note**: Any Phylum string equal to `Actinobacteria` is normalized to `Actinomycetota`. Phyla labelled as `Unknown` / `uncultured` / `uncultivated` are excluded from color-coded plots and legends.

---

## 2. Approach

High-level steps:

1. **Read FASTA** and an optional **Table_S2** metadata file.
2. Derive **Phylum** from `Lineage` (if provided) and normalize names.
3. Optionally **subsample sequences** (25/50/75/100%) for MSA and analyses.
4. Run **MAFFT** → produce a **multiple sequence alignment**.
5. From the aligned sequences, compute **pairwise identities** and derive **dissimilarity**.
6. Visualize:
   - **UMAP** on dissimilarity,
   - a **dendrogram** from dissimilarity,
   - **identity distribution** over bins.
7. Run **MMseqs2 clustering** across identity thresholds, and visualize intra-cluster diagnostics.
8. *(Part B)* Run **Word2Vec** k-mer embedding analysis with grid search over dimensions and epochs, re-generating all W2V-dependent figures and a summary correlation table (identity vs cosine similarity).

---

## 3. Input Data & Metadata

- **FASTA**: protein sequences (headers’ first token used as `sequence_id`).
- **Table_S2.tsv (optional)**: a tab-separated table with at least:
  - `Protein Accession` (or compatible ID that matches FASTA headers),
  - `Lineage` (e.g., `Bacteria; Pseudomonadota; Gammaproteobacteria; ...`).

Phylum is extracted from `Lineage` **position 2** if present; otherwise it falls back to a single-token lineage that is **not** a domain-level token (e.g. not `Bacteria`, `Archaea`, `Eukaryota`).

---

## 4. Software & Environment

- **MAFFT** (e.g., `/usr/bin/mafft`)
- **MMseqs2** (e.g., `<conda-env>/bin/mmseqs`)
- **Python 3.9+** with libraries:
  - `numpy`, `pandas`, `matplotlib`
  - `umap-learn` (UMAP)
  - `scikit-learn` (t-SNE, pairwise distances, metrics)
  - `scipy` (hierarchical clustering/dendrogram)
  - `biopython` (FASTA parsing; recommended)
  - `gensim` (Word2Vec; used in Part B)

The CLI tries to auto-locate MAFFT and MMseqs2 via `PATH`, `CONDA_PREFIX`, and `whereis`. You can force paths with `--mafft-bin` and `--mmseqs-bin`.

---

## 5. Mathematical Definitions

### 5.1 Aligned Sequences

Let the aligned sequence of protein $i$ be

$$
S_i = [s_{i,1}, s_{i,2}, \dots, s_{i,L}], \quad s_{i,t} \in \{\text{amino acids} \cup \{-\}\}.
$$

Here, `-` denotes a gap introduced by the MSA.  
For **Word2Vec** we first remove gaps to obtain a gap-free sequence

$$
\tilde{S}_i = [\tilde{s}_{i,1}, \tilde{s}_{i,2}, \dots, \tilde{s}_{i,\tilde{L}}], \quad \tilde{s}_{i,u} \in \{\text{amino acids}\}.
$$

### 5.2 Pairwise Identity

For sequences $i$ and $j$, the residue at aligned position $t$ is $A_i(t)$ and $A_j(t)$. Using `-` for gaps, the **pairwise identity** (in %) is

$$
\mathrm{Id}(i,j) = 100 \times
\frac{\sum_{t} \mathbf{1}\big[\,A_i(t)=A_j(t) \wedge A_i(t)\neq - \wedge A_j(t)\neq -\,\big]}
     {\sum_{t} \mathbf{1}\big[\,A_i(t)\neq - \wedge A_j(t)\neq -\,\big]}\,.
$$

Only positions where both sequences are non-gapped contribute to the denominator.

### 5.3 Dissimilarity

From identity we define a **dissimilarity**:

$$
d_{ij} = 1 - \frac{\mathrm{Id}(i,j)}{100}\,.
$$

This is a symmetric, non-negative quantity in $[0,1]$.  
For UMAP based on identity, we pass a **precomputed** matrix $D = [d_{ij}]$.

### 5.4 Word2Vec k-mer Embeddings

From each **gap-free** sequence $\tilde{S}_i$ we define overlapping k-mers (default $k = 3$) with stride $s$ (default $s = 1$):

$$
k_{i,u} = \tilde{s}_{i,u} \, \tilde{s}_{i,u+1} \dots \tilde{s}_{i,u+k-1},
\qquad u = 1,\dots, (\tilde{L}-k+1).
$$

Each k-mer acts as a **token** in the Word2Vec vocabulary. Let $M_i$ be the number of k-mers in sequence $i$. Across the training set we define

$$
m_{\min} = \min_i M_i
$$

the minimum number of k-mers across all sequences used to **standardize** sequence embeddings.

For each sequence:

- If $M_i \ge m_{\min}$, we keep only the first $m_{\min}$ k-mers.
- If $M_i < m_{\min}$, we pad with zero vectors so the effective length is $m_{\min}$.

Let $\mathbf{w}(k_{i,u}) \in \mathbb{R}^d$ denote the Word2Vec embedding of the $u$-th k-mer (dimension $d$ is `vector_size`, e.g. 100 or 390). After truncation/padding each sequence has **exactly** $m_{\min}$ k-mer embeddings.

We then build a **sequence-level** embedding. In this project, we use **mean aggregation**:

$$
\mathbf{v}_i = \frac{1}{m_{\min}} \sum_{u=1}^{m_{\min}} \mathbf{w}(k_{i,u}) \in \mathbb{R}^d.
$$

These $\mathbf{v}_i$ are the W2V sequence embeddings used in cosine-similarity analyses and W2V-based dendrograms.

Cosine similarity between two embeddings $\mathbf{v}_i, \mathbf{v}_j \in \mathbb{R}^d$ is

$$
\cos(\theta_{ij}) = \frac{\mathbf{v}_i \cdot \mathbf{v}_j}{\|\mathbf{v}_i\|\,\|\mathbf{v}_j\|}\,,
$$

and the associated cosine distance is $1 - \cos(\theta_{ij})$.

**Word2Vec default configuration in FAALProt Part B** (for each grid combination unless overridden):

- `sg = 1` (skip-gram)
- `hs = 0` (hierarchical softmax off)
- `negative = 5` (negative sampling)
- `window = 5`
- `min_count = 1`
- `workers = 48`
- `vector_size` in the user-specified list (defaults include 100 and 390)
- `epochs` in the user-specified list (e.g. 200, 500, 1500, 2500)
- fixed random seed for reproducibility.

### 5.5 Identity Bins and Pair Sampling

Pairwise identities $\mathrm{Id}(i,j)$ are grouped into **identity bins** over the range 10–100%, typically in steps of 5% or 10% (e.g. 10–20, 20–30, …, 90–100), depending on the figure.

To limit the number of **plotted** pairs while keeping coverage of the full identity space, we deterministically subsample the upper triangle of all pairs $(i,j)$, $i<j$, using a **stride**.

Let $N$ be the number of sequences, and

$$
P = \frac{N(N-1)}{2}
$$

be the total number of unordered pairs. Given a desired number of sampled pairs $S$ (e.g. $600{,}000$ or $1{,}000{,}000$), we define

$$
\text{stride} = \left\lfloor \frac{P}{S} \right\rfloor.
$$

Then we enumerate all pairs in a fixed deterministic order and keep every `stride`-th pair. This yields a deterministic, evenly spaced subsample without bias toward any region of the triangle.

### 5.6 Correlation: Identity vs Cosine Similarity

To quantify how well the embedding space reflects alignment identity, we compute the **Pearson correlation** between:

- $x_n = \mathrm{Id}(i_n, j_n)$: identity for the $n$-th sampled pair (either in % or normalized to $[0,1]$),
- $y_n = \cos(\theta_{i_n j_n})$: cosine similarity of the W2V mean embeddings for the same pair,

for sampled pairs $(i_n, j_n)$, $n = 1,\dots,N$.

Let the sample means be

$$
\bar{x} = \frac{1}{N}\sum_{n=1}^N x_n,
\qquad
\bar{y} = \frac{1}{N}\sum_{n=1}^N y_n.
$$

The Pearson correlation coefficient is

$$
r_{xy} =
\frac{\sum_{n=1}^N (x_n - \bar{x})(y_n - \bar{y})}
     {\sqrt{\sum_{n=1}^N (x_n - \bar{x})^2}\,
      \sqrt{\sum_{n=1}^N (y_n - \bar{y})^2}}\,.
$$

For each combination of subset, embedding dimension, and epochs, FAALProt computes $r_{xy}$ and reports the values in a **summary table**, which can be used to compare which configuration best captures the alignment signal in the embedding space.

---

## 6. Pipeline Stages

### 6.1 Multiple Sequence Alignment (MAFFT)

- Aligns the selected subset with **MAFFT** (default `--auto`).
- Alignment is performed separately for each subset fraction (25%, 50%, 75%, 100%) if requested.
- The aligned sequences (with gaps) are used for identity computation (Section 5.2).  
- For Word2Vec, gaps are removed before k-mer tokenization (Section 5.4).

### 6.2 Pairwise Identities & Dissimilarity Matrix

- Computes pairwise identities for all pairs in the aligned subset.
- Derives a **dissimilarity matrix** $D = [d_{ij}]$ using $d_{ij} = 1 - \mathrm{Id}(i,j)/100$.
- This dissimilarity is used for:
  - **UMAP** (`metric="precomputed"`),
  - **hierarchical clustering** (linkage on $d_{ij}$),
  - several identity-based diagnostics (barplots, bin counts).

### 6.3 Subsampling of Sequences

To handle large datasets, FAALProt can operate on subsets of the full FASTA:

- fixed **fractions**: 25%, 50%, 75%, 100% of the sequences,
- or user-defined subset sizes (depending on CLI options).

For each subset, the full pipeline (alignment, identities, MMseqs2, W2V Part B) is run, and figures are saved into separate subdirectories.

### 6.4 MMseqs2 Clustering Across Identity Cutoffs

- Clusters the subset using **MMseqs2** at identity cutoffs 10%, 20%, …, 90%.
- For each cutoff, FAALProt records cluster memberships and generates aggregations:
  - **intra-cluster** pair counts per identity bin,
  - **cluster size** distributions (mean, median, P95),
  - **largest cluster size** vs cutoff.

Diagnostic figures include:

- `06a_counts_intra_identity_bins_by_cutoff_lines_REALX` (lines across cutoffs),
- `06b_counts_intra_identity_bins_by_cutoff_heatmap` (heatmap of counts),
- `04_cluster_size_vs_identity_cutoff`,
- `05_largest_cluster_by_identity_cutoff`.

### 6.5 Word2Vec Grid Search (Part B)

When `--run-part-b` is enabled, FAALProt performs a **grid search** over Word2Vec hyperparameters:

- embedding dimensions: `--w2v-dims`, e.g. `100,390` (default list includes 100 and 390),
- training epochs: `--w2v-epochs-list`, e.g. `200,500,1500,2500`.

For each **subset fraction** (25%, 50%, 75%, 100%) and each `(dimension, epochs)` pair, FAALProt:

1. tokenizes aligned, gap-free sequences into k-mers,
2. trains a Word2Vec model with the consistent configuration described in Section 5.4,
3. builds mean sequence embeddings,
4. computes:
   - cosine similarities for sampled pairs,
   - identity vs cosine scatter plots (intra/inter-phylum),
   - an embedding-based dendrogram,
   - UMAP/t-SNE on embeddings,
   - and the Pearson correlation $r_{xy}$ (Section 5.6).

All W2V-dependent figures are generated for **every** combination and appropriately renamed to encode subset, dimension, and epoch in the filenames. A **summary table** aggregates the correlation values and other statistics to facilitate model comparison.

### 6.6 Visualizations and Figures

Key figure types (filenames refer to the *base* name; extensions include `.png`, `.tiff`, `.svg`):

- `01_bar_identity_10_100`  
  Barplot of identity distribution over bins from 10–100% (5% or 10% width).

- `03_dendrogram_branches_by_phylum`  
  Hierarchical clustering tree colored by Phylum on branches.

- `04_umap_dissimilarity_by_phylum`  
  UMAP of dissimilarity, points colored by Phylum.

- `06a_counts_intra_identity_bins_by_cutoff_lines_REALX`  
  Line plot of intra-cluster pair counts across identity bins and MMseqs2 cutoffs.

- `06b_counts_intra_identity_bins_by_cutoff_heatmap`  
  Heatmap of intra-cluster pair counts (identity bins × cutoffs).

- `07_cosine_vs_identity_intra`  
  Cosine similarity vs identity, **intra-phylum** pairs only (requires Part B).

- `08_cosine_vs_identity_inter`  
  Cosine similarity vs identity, **inter-phylum** pairs (requires Part B).

- `09_tsne_w2v_mean_by_phylum`  
  t-SNE of W2V mean embeddings, colored by Phylum (Part B).

- `10_umap_w2v_mean_by_phylum`  
  UMAP of W2V mean embeddings, colored by Phylum (Part B).

- `13_dendrogram_w2v_mean_by_phylum`  
  Dendrogram built from W2V cosine distances, branch-colored by Phylum (Part B).

For each subset and each `(dimension, epochs)` configuration, figure filenames may be prefixed or suffixed with abbreviations encoding those settings to avoid collisions.

All plots:

- have **no title** in the saved figure,
- use legends placed **below** the main panel,
- are exported in **PNG (900 dpi)**, **TIFF (900 dpi)**, and **SVG**.

---

## 7. Color Policy (Distinct Phylum Colors)

- We construct a **global palette** across all phyla using `tab20`, `tab20b`, `tab20c` as base, with HSV-based fallbacks when needed.
- **Normalization**: `Actinobacteria` → `Actinomycetota`.
- **Excluded from color plots**: labels containing `Unknown`, `uncultured`, or `uncultivated`.
- The same palette is reused across all plots (UMAP, dendrograms, identity plots) for visual consistency.
- Mixed-phylum clades in dendrograms are colored in neutral gray to highlight heterogeneity.

---

## 8. Outputs & Directory Layout

For an output directory `results_faal`, and a subset label (e.g., `subset_frac_0.25_dim_100_epochs_500`):

```text
results_faal/
  subset_0.25/
    alignment.fasta
    pairs.csv                      # all pairwise identities for this subset
    clusters_multi_threshold.csv   # MMseqs2 clusters for all cutoffs
    subset_ids.csv
    subset_phylum_map.csv
    w2v_dim_100_epochs_500/
      01_bar_identity_10_100.(png|tiff|svg)
      03_dendrogram_branches_by_phylum.(png|tiff|svg)
      04_umap_dissimilarity_by_phylum.(png|tiff|svg)
      06a_counts_intra_identity_bins_by_cutoff_lines_REALX.(png|tiff|svg)
      06b_counts_intra_identity_bins_by_cutoff_heatmap.(png|tiff|svg)
      07_cosine_vs_identity_intra.(png|tiff|svg)
      08_cosine_vs_identity_inter.(png|tiff|svg)
      09_tsne_w2v_mean_by_phylum.(png|tiff|svg)
      10_umap_w2v_mean_by_phylum.(png|tiff|svg)
      13_dendrogram_w2v_mean_by_phylum.(png|tiff|svg)
      correlation_summary.tsv      # identity vs cosine correlation for this configuration
  subset_0.50/
    ...
  subset_0.75/
    ...
  subset_1.00/
    ...
  FIGURE_MANIFEST.txt
```

Exact names may differ slightly depending on CLI options, but all plots are systematically organized by subset, dimension, and epochs.

---

## 9. Reproducibility, Performance & Practical Notes

- A global **random seed** (42) is used for:
  - subset selection,
  - UMAP/t-SNE,
  - Word2Vec initialization (via `seed`),
  - and any sampling logic.

- Pairwise identity computation streams the upper triangle and allows **deterministic subsampling** (Section 5.5) for plotting (e.g. 600k and 1M pairs).

- MSA is typically the **dominant cost** for large $N$. Using subsets (25–50%) can drastically reduce runtime while still preserving diversity.

- MMseqs2 clustering is run on the **same subset** as the MSA, so identity distributions and cluster diagnostics are directly comparable.

- If metadata is missing, Phylum is labeled `Unknown` and filtered out from color plots but **not** from numerical computations.

- Use `--mafft-bin` and `--mmseqs-bin` to force binary paths when Conda/`PATH` resolution differs across systems.

---

## 10. How to Run — CLI Version

The main script is `faalprot_heterogeneity_cli_v2.py`. Use `--help` for the full list of options.

### 10.1 Minimal Run (Alignment + Identity + MMseqs2, No W2V)

```bash
python faalprot_heterogeneity_cli_v2.py \
  --fasta FAAL_NR_MIBIG_LEGE_FFT_NS_2_raw_data.fasta \
  --table Table_S2.tsv \
  --outdir results_faal
```

### 10.2 Enable Word2Vec Part B with Default Grid

By default, the grid includes embedding dimensions 100 and 390, and a list of epoch values (e.g. 200, 500, 1500, 2500) defined in the script.

```bash
python faalprot_heterogeneity_cli_v2.py \
  --fasta FAAL_NR_MIBIG_LEGE_FFT_NS_2_raw_data.fasta \
  --table Table_S2.tsv \
  --outdir results_faal_w2v \
  --run-part-b
```

### 10.3 Custom Word2Vec Grid (Dimensions and Epochs)

```bash
python faalprot_heterogeneity_cli_v2.py \
  --fasta my_sequences.fasta \
  --table table_S2.tsv \
  --outdir resultados_grid \
  --run-part-b \
  --w2v-dims 100,390 \
  --w2v-epochs-list 200,500,1500,2500
```

### 10.4 Restrict to a Single Subset Fraction

If you only want to run, for example, the 50% subset, use the appropriate subset option (exact flag name depends on the implementation; e.g.):

```bash
python faalprot_heterogeneity_cli_v2.py \
  --fasta my_sequences.fasta \
  --table table_S2.tsv \
  --outdir results_50 \
  --subset-fractions 0.50 \
  --run-part-b \
  --w2v-dims 100,390 \
  --w2v-epochs-list 500,1500
```

### 10.5 Explicit Binary Paths

```bash
python faalprot_heterogeneity_cli_v2.py \
  --fasta FAAL_NR_MIBIG_LEGE_FFT_NS_2_raw_data.fasta \
  --table Table_S2.tsv \
  --outdir results_bins \
  --mafft-bin /usr/bin/mafft \
  --mmseqs-bin /home/USER/anaconda3/envs/faal/bin/mmseqs \
  --run-part-b
```

---

## 11. FAQ

**Q1. Why sample pairs instead of plotting all of them?**  
The number of pairs grows as $N(N-1)/2$. Even when we **compute** all identities, we only **plot** a deterministic subsample (stride-based) to keep plots interpretable and file sizes manageable.

**Q2. Why exclude “Unknown/uncultured” from color plots?**  
Color is used to encode discrete taxa (Phylum). Including unknown labels would add many arbitrary colors with limited biological meaning and clutter legends.

**Q3. Why apply Word2Vec to gap-free sequences but compute identity on the gapped alignment?**  
Word2Vec models sequence **content/context** (k-mer tokens) and benefits from natural, gap-free sequences. Alignment identity, on the other hand, must respect positional correspondence enforced by the MSA; thus it uses the **gapped** alignment.

**Q4. How are dendrogram colors assigned?**  
Branches are colored by the **unique** phylum if all descendants share it; otherwise the branch is shown in neutral gray, highlighting mixed clades and heterogeneity.

**Q5. How do I interpret the identity vs cosine correlation table?**  
For each `(subset, dimension, epochs)` configuration, FAALProt reports the Pearson correlation $r_{xy}$ between identity and cosine similarity. Higher absolute values of $r_{xy}$ indicate that the embedding space more faithfully reflects alignment similarity patterns. This can guide the choice of embedding dimension and training epochs.

---
## Reference

**If you use the FAALPred heterogeneity repository, please cite the following article and this repository:**

Diversity of FAAL enzymes and prediction of their substrate specificity using FAALPred
DOI: 10.1002/pro.70468. Authors: Leandro de Mattos Pereira†, Anne Liong†, and Pedro Leão. †These authors contributed equally to this work. Protein Science, 2026 (Article in production)

###

## 12. References

- Katoh, K. & Standley, D.M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Mol. Biol. Evol.*  
- Steinegger, M. & Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching. *Nat. Biotechnol.*  
- McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv:1802.03426.  
- van der Maaten, L. & Hinton, G. (2008). Visualizing data using t-SNE. *JMLR*.  
- Mikolov, T., Chen, K., Corrado, G., & Dean, J. (2013). Distributed representations of words and phrases and their compositionality. *NIPS* (Word2Vec).  
- Virtanen, P. *et al.* (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nat. Methods*.  
- Pedregosa, F. *et al.* (2011). Scikit-learn: Machine Learning in Python. *JMLR*.  
- Hunter, J.D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*.  
- Cock, P.J.A. *et al.* (2009). Biopython: freely available Python tools for computational biology and bioinformatics. *Bioinformatics*.
