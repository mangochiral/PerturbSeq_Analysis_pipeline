# CRISPRa sgRNA re-annotation & de-duplication

Re-annotates and cleans the genome-wide CRISPRa sgRNA library (hCRISPRa-v2 / Calabrese Set A + Set B):
reconciles outdated gene names to current Ensembl IDs, de-duplicates frameshift-twin guides, aligns every
guide to hg38, and assigns/validates the target gene with full neighborhood + off-target annotation.

## Environment

Run under the `gwt-env` conda environment (has bowtie2, pysam, biopython, pyfaidx, gffutils, pandas):

```
conda activate gwt-env      # or: /Users/rzhu/miniconda3/envs/gwt-env
jupyter nbconvert --to notebook --execute --inplace <notebook>.ipynb
```

`environment.yml` describes an equivalent standalone env if `gwt-env` is unavailable.

## Run order

1. **`dedup_frameshift_twins.ipynb`** — build the unified library, reconcile gene symbols → Ensembl IDs
   (GENCODE v48 + HGNC), detect & collapse frameshift twins. *(base conda is sufficient; ~1 min)*
2. **`build_gene_models.ipynb`** — build `genome/genes_df_subset_for_sgRNA_annotation.parquet`
   (primary assembly + alt/patch contigs). *(needs internet on first run)*
3. **`align_annotate_sgRNA.ipynb`** — align to hg38 (bowtie2), PAM check, locus-based gene assignment,
   frameshift-twin revisit, neighborhood/off-target columns, and the final annotated table. *(~1 min)*

Each notebook is idempotent: input files are downloaded/copied only if missing, and the bowtie2 alignment
is skipped if its SAM already exists (delete it to force re-alignment).

## Inputs

**Auto-acquired by the notebooks** (into `genome/`, on first run):
`hgnc_complete_set.txt` (HGNC), `hg38.chromAlias.txt` (UCSC), `gencode.v48.chr_patch_hapl_scaff` gtf,
and a copy of `gene_annotations_hg38.parquet`.

**Must pre-exist** (referenced by absolute path from the sibling `GRNPerturbSeq` project — large upstream
files, not downloaded here):
- bowtie2 hg38 index + `hg38.fa` (`GRNPerturbSeq/2_files/hg38_genome/`)
- parsed GENCODE v48 gene/transcript/CDS parquets (`GRNPerturbSeq/.../5_sgRNA_annotation/genome/`)

**Library source:** `1_files/CRISPRa_library_construction/CRISPRa_pool{1,2}.csv`.

## Key outputs (`results/`)

- `CRISPRa_targeting_sgRNA_annotated.{parquet,csv}` — **the master annotation**: all 78,370 sgRNA
  (targeting + non-targeting), 57 columns, with a plain-English `note` per guide.
- `CRISPRa_targeting_sgRNA_final_deduplicated.{parquet,csv}` — the cleaned, collapsed library (`keep=True`).
- `CRISPRa_symbol_to_geneid_reconciliation`, `CRISPRa_gene_assignment_corrections`,
  frameshift-twin tables, `CRISPRa_unaligned_guides_notes`, `CRISPRa_adjacent_offset1_different_gene`, etc.

`sgRNAalign_util.py` holds the shared helper functions.
