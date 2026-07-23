"""
Utility functions for CRISPRa sgRNA re-annotation and de-duplication.

Adapted from the GRNPerturbSeq hCRISPRi-v2 annotation pipeline
(GWT_perturbseq_analysis/src/5_sgRNA_annotation/sgRNAalign_util.py).

Two groups of helpers:
  1. Library-intrinsic (Phase 1-2): FASTA export + frameshift-twin detection.
     These are pure-python / pandas and run under base conda.
  2. Alignment / annotation (Phase 3-6): SAM parsing, PAM checks, TSS/gene
     proximity. These require pysam / Bio / pyfaidx and are used later.

hCRISPRa-v2 note: the discriminating protospacer unit is 19 bp (the 5' base is
an appended/unreliable base), so alignment and twin detection operate on the
last 19 bp of the 20 bp barcode sequence.
"""

import csv
from collections import defaultdict

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Phase 1.5 : gene-name / Ensembl-ID reconciliation (base conda only)
# ---------------------------------------------------------------------------

def reconcile_symbols_to_geneid(symbols, gencode_gene_df, hgnc_df):
    """
    Map (possibly outdated) gene symbols to a current Ensembl gene_id + symbol.

    The library carries only legacy hCRISPRa-v2 / Calabrese symbols and no IDs.
    Resolution order per symbol S:
        1. ``v48_direct``   -- S is a unique gene_name in GENCODE (use its id).
        2. ``ambiguous_v48``-- S maps to >1 GENCODE gene_id (e.g. PAR X/Y genes) -> flag.
        3. ``hgnc_approved``-- S is a current HGNC approved symbol.
        4. ``hgnc_prev``    -- S is a unique HGNC *previous* symbol (a rename).
        5. ``ambiguous_prev``/``hgnc_alias``/``ambiguous_alias`` -- alias routes.
        6. ``unresolved``   -- none of the above.
    Ambiguous / unresolved rows get gene_id = NaN (caller keeps the old label and
    defers to genomic alignment in Phase 3).

    Args:
        symbols (iterable[str]): gene symbols to resolve (deduplicated internally).
        gencode_gene_df (pd.DataFrame): GENCODE gene table with 'gene_name','gene_id'
            (gene_id may carry a version suffix; it is stripped).
        hgnc_df (pd.DataFrame): HGNC complete set with 'symbol','prev_symbol',
            'alias_symbol','ensembl_gene_id' (pipe-delimited multi-values).

    Returns:
        pd.DataFrame: library_symbol, gene_id, current_symbol, resolution, n_candidates.
    """
    # GENCODE gene_name -> set(gene_id)
    g = gencode_gene_df.copy()
    g['gene_id'] = g['gene_id'].str.split('.').str[0]
    v48 = defaultdict(set)
    for name, gid in zip(g['gene_name'], g['gene_id']):
        v48[name].add(gid)
    v48_unique = {n: next(iter(ids)) for n, ids in v48.items() if len(ids) == 1}
    v48_multi = {n for n, ids in v48.items() if len(ids) > 1}

    # HGNC maps
    hg = hgnc_df[hgnc_df['ensembl_gene_id'].notna()]
    approved = dict(zip(hg['symbol'], hg['ensembl_gene_id']))
    approved_sym = dict(zip(hg['ensembl_gene_id'], hg['symbol']))
    prev, alias = defaultdict(set), defaultdict(set)
    for sym, eid, pv, al in zip(hg['symbol'], hg['ensembl_gene_id'],
                                hg['prev_symbol'], hg['alias_symbol']):
        if isinstance(pv, str):
            for x in pv.split('|'):
                prev[x].add(eid)
        if isinstance(al, str):
            for x in al.split('|'):
                alias[x].add(eid)

    def resolve(s):
        if s in v48_unique:
            return v48_unique[s], s, 'v48_direct', 1
        if s in v48_multi:
            return np.nan, s, 'ambiguous_v48', len(v48[s])
        if s in approved:
            return approved[s], s, 'hgnc_approved', 1
        if s in prev:
            c = prev[s]
            if len(c) == 1:
                e = next(iter(c))
                return e, approved_sym.get(e, s), 'hgnc_prev', 1
            return np.nan, s, 'ambiguous_prev', len(c)
        if s in alias:
            c = alias[s]
            if len(c) == 1:
                e = next(iter(c))
                return e, approved_sym.get(e, s), 'hgnc_alias', 1
            return np.nan, s, 'ambiguous_alias', len(c)
        return np.nan, s, 'unresolved', 0

    rows = [(s, *resolve(s)) for s in sorted(set(symbols))]
    return pd.DataFrame(rows, columns=['library_symbol', 'gene_id',
                                       'current_symbol', 'resolution', 'n_candidates'])


# ---------------------------------------------------------------------------
# Phase 1-2 : library-intrinsic helpers (base conda only)
# ---------------------------------------------------------------------------

def convert_csv_to_fasta(input_csv_file, output_fasta_file,
                         name_col='name', seq_col='sgRNA', drop_first_base=True):
    """
    Write a FASTA from a CSV of guides.

    By default drops the first base of each protospacer (``seq[1:]``) so the
    aligned unit is the 19 bp discriminating window used by hCRISPRa-v2.

    Args:
        input_csv_file (str): input CSV path (must contain ``name_col``/``seq_col``).
        output_fasta_file (str): output FASTA path.
        name_col (str): column holding the FASTA record id.
        seq_col (str): column holding the protospacer sequence.
        drop_first_base (bool): if True, emit ``seq[1:]`` (last 19 bp).
    """
    print(f"Reading from '{input_csv_file}'...")
    with open(input_csv_file, mode='r', encoding='utf-8') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        with open(output_fasta_file, mode='w', encoding='utf-8') as fasta_file:
            record_count = 0
            for row in csv_reader:
                name = row.get(name_col, '').strip()
                seq = row.get(seq_col, '').strip()
                if drop_first_base:
                    seq = seq[1:]
                if name and seq:
                    fasta_file.write(f">{name}\n{seq}\n")
                    record_count += 1
    print(f"Successfully converted {record_count} records -> '{output_fasta_file}'")


def frameshift_twin_pairs(df, id_col='guide_id', seq_col='protospacer',
                          gene_col='target_gene_symbol', same_gene_only=True):
    """
    Detect frameshift-twin guide pairs from sequence alone.

    A frameshift twin is a pair of guides whose 20 bp protospacers share a
    19 bp window at a 1 bp offset (same orientation). Concretely, guides A and
    B are twins when any of these hold:

        A[:19] == B[:19]   (identical 19-mer, "same-window" near-duplicate)
        A[1:]  == B[1:]     (identical 19-mer, "same-window" near-duplicate)
        A[1:]  == B[:19]   (A is B shifted +1 bp; classic frameshift twin)
        A[:19] == B[1:]     (A is B shifted -1 bp; classic frameshift twin)

    These are exactly the pairs whose RHS capture probes cross-hybridize
    (probe = revcomp of the 20-mer, complementary over 19 bp), and the pairs
    that an exact-``last19bp`` dedup misses because their last-19 differ.

    Args:
        df (pd.DataFrame): one row per guide.
        id_col, seq_col, gene_col (str): column names.
        same_gene_only (bool): keep only pairs targeting the same gene (the
            redundant-guide case the library curation should collapse). When
            False, cross-gene shared-window pairs are also returned (probe
            cross-hyb candidates).

    Returns:
        pd.DataFrame with one row per unordered twin pair:
            guide_a, guide_b, gene_a, gene_b, same_gene,
            shared_19mer, relation
        where ``relation`` is one of {'identical_head','identical_tail',
        'shift+1','shift-1'}.
    """
    # Index every guide by both of its 19 bp windows.
    #   window 'head' = protospacer[:19]
    #   window 'tail' = protospacer[1:]
    head = {}  # 19-mer -> list of guide indices whose head == 19-mer
    tail = {}
    seqs = df[seq_col].str.upper().values
    ids = df[id_col].values
    genes = df[gene_col].values

    for i, s in enumerate(seqs):
        if not isinstance(s, str) or len(s) < 20:
            continue
        head.setdefault(s[:19], []).append(i)
        tail.setdefault(s[1:], []).append(i)

    pairs = {}  # (min_i, max_i) -> dict of pair attributes

    def _add(i, j, relation, kmer):
        if i == j:
            return
        key = (i, j) if i < j else (j, i)
        # First relation found wins; identical windows take precedence label.
        if key not in pairs:
            a, b = key
            pairs[key] = dict(
                guide_a=ids[a], guide_b=ids[b],
                gene_a=genes[a], gene_b=genes[b],
                same_gene=bool(genes[a] == genes[b]),
                shared_19mer=kmer, relation=relation,
            )

    # identical 19-mer at the same window position
    for kmer, idxs in head.items():
        if len(idxs) > 1:
            for a in range(len(idxs)):
                for b in range(a + 1, len(idxs)):
                    _add(idxs[a], idxs[b], 'identical_head', kmer)
    for kmer, idxs in tail.items():
        if len(idxs) > 1:
            for a in range(len(idxs)):
                for b in range(a + 1, len(idxs)):
                    _add(idxs[a], idxs[b], 'identical_tail', kmer)

    # 1 bp frameshift: one guide's tail (proto[1:]) == other's head (proto[:19])
    for kmer, idxs_tail in tail.items():
        idxs_head = head.get(kmer)
        if not idxs_head:
            continue
        for i in idxs_tail:      # i shifted +1 relative to j
            for j in idxs_head:
                if i != j:
                    _add(i, j, 'shift+1', kmer)

    out = pd.DataFrame(list(pairs.values()))
    if out.empty:
        return pd.DataFrame(columns=['guide_a', 'guide_b', 'gene_a', 'gene_b',
                                     'same_gene', 'shared_19mer', 'relation'])
    if same_gene_only:
        out = out[out['same_gene']].reset_index(drop=True)
    return out.reset_index(drop=True)


def collapse_twin_groups(pairs_df, rank_series):
    """
    Turn a table of twin pairs into connected-component groups and pick one
    representative per group (the lowest guide_rank magnitude = highest-ranked
    guide is kept; the rest are marked collapsed).

    Args:
        pairs_df (pd.DataFrame): output of ``frameshift_twin_pairs`` (guide_a/guide_b).
        rank_series (pd.Series): index = guide_id, value = guide_rank
            (negative ints; closer to 0 == better rank).

    Returns:
        pd.DataFrame: one row per guide in any twin group with columns
            guide_id, twin_group_id, is_representative, collapsed_into.
    """
    # union-find over guide ids
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        parent[find(x)] = find(y)

    for a, b in zip(pairs_df['guide_a'], pairs_df['guide_b']):
        union(a, b)

    groups = {}
    for g in list(parent):
        groups.setdefault(find(g), []).append(g)

    rows = []
    for gi, (root, members) in enumerate(sorted(groups.items())):
        # representative = best (largest, i.e. closest to 0) guide_rank
        ranks = {m: rank_series.get(m, -np.inf) for m in members}
        rep = max(members, key=lambda m: ranks[m])
        for m in sorted(members):
            rows.append(dict(
                guide_id=m,
                twin_group_id=gi,
                is_representative=(m == rep),
                collapsed_into=(np.nan if m == rep else rep),
            ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Phase 3-6 : alignment / annotation helpers (need pysam / Bio / pyfaidx)
# ---------------------------------------------------------------------------

def sam_to_dataframe(sam_file_path):
    """Parse a SAM file into a DataFrame of alignments (requires pysam/Bio).

    One row per alignment record (a multi-mapping read yields several rows).
    Columns: sgRNA, chromosome, pos (0-based), seq, strand, nm (#mismatches).
    """
    import pysam
    from Bio.Seq import Seq

    alignments = []
    with pysam.AlignmentFile(sam_file_path, "r") as samfile:
        for read in samfile:
            if read.is_unmapped:
                alignments.append({'sgRNA': read.query_name, 'chromosome': None,
                                   'pos': np.nan, 'seq': read.query_sequence,
                                   'strand': None, 'nm': np.nan})
                continue
            alignments.append({
                'sgRNA': read.query_name,
                'chromosome': read.reference_name,
                'pos': read.reference_start,  # 0-based
                'seq': (str(Seq(read.query_sequence).reverse_complement())
                        if read.is_reverse else read.query_sequence),
                'strand': '-' if read.is_reverse else '+',
                'nm': read.get_tag('NM') if read.has_tag('NM') else np.nan,
            })
    return pd.DataFrame(alignments)


def build_tss_index(genes_df, chrom_col='chromosome', tss_col='tss'):
    """
    Build per-chromosome sorted TSS arrays for O(log n) nearest-TSS lookup.

    Returns dict: chrom (no 'chr' prefix) -> (sorted_pos_array, gene_ids, gene_names).
    """
    chrom_tss = defaultdict(list)
    for gid, gn, ch, tss in zip(genes_df['gene_id'], genes_df['gene_name'],
                                genes_df[chrom_col].str.replace('chr', '', regex=False),
                                genes_df[tss_col]):
        for t in tss:
            chrom_tss[ch].append((int(t), gid, gn))
    out = {}
    for ch, lst in chrom_tss.items():
        lst.sort()
        out[ch] = (np.array([x[0] for x in lst]),
                   [x[1] for x in lst], [x[2] for x in lst])
    return out


def nearest_tss(tss_index, chrom, pos):
    """Nearest TSS to (chrom, pos). Returns (gene_id, gene_name, distance) or NaNs."""
    chrom = str(chrom).replace('chr', '')
    if chrom not in tss_index:
        return (np.nan, np.nan, np.nan)
    posarr, gids, gns = tss_index[chrom]
    i = int(np.searchsorted(posarr, pos))
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(posarr):
            d = abs(int(posarr[j]) - pos)
            if best is None or d < best[0]:
                best = (d, gids[j], gns[j])
    return (best[1], best[2], best[0]) if best else (np.nan, np.nan, np.nan)


def genes_near(tss_index, chrom, pos, window):
    """
    All genes with a TSS within `window` bp of (chrom, pos).

    Returns dict gene_id -> min TSS distance (a gene with several TSS keeps its
    closest). Uses binary search on the per-chromosome sorted TSS array, so it is
    O(log n + hits) rather than scanning the whole gene table.
    """
    chrom = str(chrom).replace('chr', '')
    if chrom not in tss_index:
        return {}
    posarr, gids, gns = tss_index[chrom]
    lo = int(np.searchsorted(posarr, pos - window))
    hi = int(np.searchsorted(posarr, pos + window))
    out = {}
    for j in range(lo, hi):
        d = abs(int(posarr[j]) - pos)
        g = gids[j]
        if g not in out or d < out[g]:
            out[g] = d
    return out


def distance_to_gene_tss(pos, tss_list):
    """Min distance from pos to any TSS in tss_list (np.nan if none)."""
    if not isinstance(tss_list, (np.ndarray, list)) or len(tss_list) == 0:
        return np.nan
    return int(min(abs(int(t) - pos) for t in tss_list))


def has_ngg_pam(genome, chromosome, pos, strand, read_len=19):
    """
    True if an NGG PAM sits immediately 3' of a protospacer aligned at
    (chromosome, pos, strand). `genome` is a pyfaidx.Fasta. pos is 0-based
    start of the aligned read; read_len is the aligned length (19 bp here).
    """
    from Bio.Seq import Seq
    try:
        if strand == '+':
            s = genome[chromosome][pos + read_len: pos + read_len + 3].seq
        else:
            s = str(Seq(genome[chromosome][pos - 3: pos].seq).reverse_complement())
        return s[1:3].upper() == 'GG'
    except Exception:
        return False


def find_closest_target_info(row, distance_threshold=2000):
    """Return [gene_id, gene_name, min_dist] if designed TSS within threshold."""
    sgrna_pos = row['pos']
    tss_list = row['tss']
    if not isinstance(tss_list, (np.ndarray, list)):
        return pd.Series([np.nan, np.nan, np.nan])
    distances = [abs(sgrna_pos - tss_pos) for tss_pos in tss_list]
    min_distance = min(distances)
    if min_distance <= distance_threshold:
        return pd.Series([row['gene_id'], row['gene_name'], min_distance])
    return pd.Series([np.nan, np.nan, np.nan])


def find_nearby_genes(row, all_genes_df, distance_threshold=30000):
    """List gene_ids whose TSS is within threshold of the sgRNA on same chrom."""
    same_chrom = all_genes_df[all_genes_df['chromosome_norm'] == row['chromosome_norm']].copy()
    sgrna_pos = row['pos']
    same_chrom['distance'] = same_chrom['tss'].apply(
        lambda tss_list: min([abs(sgrna_pos - tss) for tss in tss_list]))
    return list(same_chrom[same_chrom['distance'] <= distance_threshold]['gene_id'])


def is_gene_nearby(sgrna_row, genes_df, distance_threshold=30000):
    """True if any TSS within threshold or sgRNA inside a gene body."""
    genes_on_chrom = genes_df[genes_df['chromosome'] == sgrna_row['chromosome']]
    for _, gene_row in genes_on_chrom.iterrows():
        for tss in gene_row['tss']:
            if abs(sgrna_row['pos'] - tss) <= distance_threshold:
                return True
        if gene_row['gene_start'] <= sgrna_row['pos'] <= gene_row['gene_end']:
            return True
    return False
