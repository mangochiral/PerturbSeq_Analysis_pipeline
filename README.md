# CRISPRa Analysis Pipeline

## Overview

This repository provides tools and scripts for preprocessing, analyzing, and visualizing data from CRISPR experiments.

## Directory Structure

```
src/
├── basic_processing.py
│   └── preprocess_adata.py          # Data preprocessing utilities
├── guide_assignment_parallel.py     # Assign guides in parallel
├── qc_stats.py                     # Quality control statistics (standard)
├── qc_stats_heavy_load.py          # QC for large datasets (summary, t-test)
├── guide_efficiency_qc_stats.ipynb # Guide efficiency QC with t-test
├── qc_plots/
│   ├── qc_guide_type_distribution.ipynb
│   └── qc_plot_cumulative_perturbation_fraction_of_genes.ipynb
├── pseudobulk_by_lane.py           # Pseudobulk analysis per lane
│   └── prep_DE_merge_pseudobulk.ipynb # Prepare DE and merge pseudobulk results
├── Deseq2_pseudobulk.py            # Run DESeq2 analysis on pseudobulked data
```

## Getting Started

1. **Preprocess your data** using `basic_processing.py` and `preprocess.py`.
2. **Assign guides** with `guide_assignment_parallel.py`.
3. **Quality control**: Use `qc_stats.py`, `qc_stats_heavy_load.py`, or the notebooks in `qc_plots/` for QC and visualization.
4. **Pseudobulk analysis**: See `pseudobulk_by_lane.py` and related notebooks.
5. **Differential expression analysis**: Use `Deseq2_pseudobulk.py`.

## Notes

- Jupyter notebooks (`.ipynb`) provide interactive analysis and plotting.
- Each script is modular—run the steps relevant for your experiment.
