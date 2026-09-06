#!/usr/bin/env python3
"""Derive an amplicon panel from the cohort's own reads, mapped to a genome.

The standard PyEuk amplicon path needs a panel: a FASTA of the amplified regions
that reads are aligned to before ``define-windows``. When that panel is unknown or
unavailable, this step reconstructs it from the data instead of assuming it.

The idea is exactly how a targeted-amplicon panel is built in the first place, and
why every ``tas_panel_b66`` header reads ``source=peak66 loc=<scaffold>:<start>-<end>``:
targeted amplicon reads do not spread over the genome, they pile up on the handful
of amplified regions. Map the reads to the reference GENOME and the amplicons appear
as sharp, contiguous coverage peaks against an otherwise empty background. Each peak
is an amplicon; its genome sequence is the panel entry.

Input is genome-mapped BAMs (the mapping is done upstream, exactly as
``define-windows`` consumes BAMs) plus the reference genome FASTA the reads were
mapped to. Output is a panel FASTA in the same format the standard workflow expects,
so it is a drop-in for the panel input of that workflow, plus a BED and a QC table.

Panel-level, not per-specimen: amplicon locations are a property of the assay, shared
by every specimen, so a representative ``--sample`` of BAMs reconstructs the panel as
well as the whole cohort would -- the same reason ``define-windows`` samples BAMs.

The peak threshold is DATA-ADAPTIVE, not a fixed depth: amplicon peaks sit orders of
magnitude above background (measured 6e4-6e5x vs <1e3x on a CDC 8-marker cohort), so a
peak is a contiguous run whose depth clears ``--min-depth-frac`` of the run's own
maximum (with an absolute floor). That scales across datasets and sequencing depths
without retuning.
"""

import argparse
import sys

import numpy as np

# pysam is an OPTIONAL dependency, declared under the `amplicon` extra (see the sibling
# amplicon modules). Import it lazily so `pyeuk --help` works without it installed.
_pysam = None


def _require_pysam():
    global _pysam
    if _pysam is None:
        try:
            import pysam as _p
        except ImportError:
            sys.exit(
                "This step reads BAM/FASTA files and needs pysam, an optional dependency.\n"
                "  pip:   pip install 'pyeuk[amplicon]'\n"
                "  conda: conda install -c bioconda pysam\n"
                "The bioconda `pyeuk` package already carries pysam."
            )
        _pysam = _p
    return _pysam


def _sample_bams(bams, k):
    """Evenly-strided sample of k BAMs across the cohort, like define_windows.

    Even striding rather than the first k so the estimate is representative when the
    input happens to be ordered (by plate, by accession, by depth).
    """
    if k <= 0 or len(bams) <= k:
        return list(bams)
    stride = len(bams) / k
    return [bams[int(i * stride)] for i in range(k)]


def accumulate_coverage(bams):
    """Sum per-base genome coverage across BAMs, only over contigs that carry reads.

    Returns {contig: np.int64 array of length contig_len}. Contigs with no mapped read
    in any sampled BAM are never allocated -- a targeted panel touches a few dozen of a
    genome's hundreds/thousands of contigs, so this is what keeps a whole-genome scan cheap.
    """
    pysam = _require_pysam()
    cov = {}
    lengths = {}
    for bam in bams:
        with pysam.AlignmentFile(bam, "rb") as fh:
            lengths = lengths or dict(zip(fh.references, fh.lengths))
            # get_index_statistics needs the .bai; contigs with mapped>0 are the only ones worth scanning
            try:
                stats = fh.get_index_statistics()
                contigs = [s.contig for s in stats if s.mapped > 0]
            except (ValueError, AttributeError):
                contigs = list(fh.references)
            for contig in contigs:
                length = lengths[contig]
                # count_coverage returns four per-base arrays (A,C,G,T); sum for depth.
                # quality_threshold=0 counts every aligned base (amplicon calling is downstream).
                c = fh.count_coverage(contig, quality_threshold=0)
                depth = np.asarray(c, dtype=np.int64).sum(axis=0)
                if contig in cov:
                    cov[contig] += depth
                else:
                    cov[contig] = depth
    return cov, lengths


def call_peaks(cov, min_depth_frac, min_abs_depth, min_length, merge_gap, pad, lengths):
    """Call contiguous coverage peaks (amplicons) with a data-adaptive threshold.

    Threshold = max(min_abs_depth, min_depth_frac * global_max_depth). A peak is a run of
    positions at or above it; runs within merge_gap bp are joined (one amplicon can dip
    mid-peak); runs shorter than min_length are dropped as background specks; each surviving
    peak is padded by `pad` bp (clamped to the contig) to catch primer flanks.
    """
    gmax = max((int(a.max()) for a in cov.values() if a.size), default=0)
    thresh = max(min_abs_depth, int(round(min_depth_frac * gmax)))
    peaks = []
    for contig, depth in cov.items():
        above = depth >= thresh
        if not above.any():
            continue
        # contiguous runs of True
        idx = np.flatnonzero(np.diff(np.concatenate(([0], above.view(np.int8), [0]))))
        starts, ends = idx[0::2], idx[1::2]  # [start, end) per run
        # merge runs separated by < merge_gap
        merged = []
        for s, e in zip(starts, ends):
            if merged and s - merged[-1][1] < merge_gap:
                merged[-1][1] = e
            else:
                merged.append([int(s), int(e)])
        for s, e in merged:
            if e - s < min_length:
                continue
            ps = max(0, s - pad)
            pe = min(lengths[contig], e + pad)
            peaks.append({
                "contig": contig, "start": ps, "end": pe,
                "length": pe - ps,
                "max_depth": int(depth[s:e].max()),
                "mean_depth": int(round(float(depth[s:e].mean()))),
            })
    peaks.sort(key=lambda p: -p["max_depth"])
    return peaks, thresh, gmax


def write_panel_fasta(peaks, genome, out, prefix):
    """Write the derived panel FASTA in the same header format as a curated panel:
    ``>{prefix}_{NN}_L{len}bp source=derive-panel loc={contig}:{start1}-{end}``
    (1-based inclusive coords in the header, matching peak66; sequence from the genome)."""
    fa = _require_pysam().FastaFile(genome)
    try:
        with open(out, "w") as fh:
            for i, p in enumerate(sorted(peaks, key=lambda q: -q["length"]), 1):
                seq = fa.fetch(p["contig"], p["start"], p["end"])
                name = f"{prefix}_{i:02d}_L{p['length']}bp"
                fh.write(f">{name} source=derive-panel loc={p['contig']}:{p['start']+1}-{p['end']}\n")
                for j in range(0, len(seq), 70):
                    fh.write(seq[j:j + 70] + "\n")
    finally:
        fa.close()


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="pyeuk derive-panel",
        description="Reconstruct an amplicon panel FASTA from genome-mapped reads (coverage peaks = amplicons).")
    p.add_argument("--bams", nargs="+", required=True,
                   help="Genome-mapped, indexed BAMs (a representative sample of the cohort).")
    p.add_argument("--genome", required=True,
                   help="Reference genome FASTA the BAMs were mapped to (indexed .fai; used to extract peak sequences).")
    p.add_argument("--sample", type=int, default=20, help="How many BAMs to scan (evenly strided). 0 = all.")
    p.add_argument("--min-depth-frac", type=float, default=0.01,
                   help="A peak must reach this fraction of the run's MAX depth (data-adaptive threshold). Default 0.01.")
    p.add_argument("--min-abs-depth", type=int, default=30, help="Absolute depth floor for the threshold. Default 30.")
    p.add_argument("--min-length", type=int, default=80, help="Minimum peak length (bp). Default 80.")
    p.add_argument("--merge-gap", type=int, default=25, help="Join peaks separated by fewer than this many bp. Default 25.")
    p.add_argument("--pad", type=int, default=0, help="Pad each peak by this many bp on both sides (primer flanks). Default 0.")
    p.add_argument("--prefix", default="ampl", help="Amplicon name prefix for the panel FASTA. Default 'ampl'.")
    p.add_argument("--out", required=True, help="Output panel FASTA path.")
    p.add_argument("--out-bed", help="Optional output BED of the peaks.")
    p.add_argument("--out-qc", help="Optional output TSV: amplicon, contig, start, end, length, max_depth, mean_depth.")
    a = p.parse_args(argv)

    bams = _sample_bams(list(a.bams), a.sample)
    print(f"[derive-panel] scanning {len(bams)} of {len(a.bams)} BAM(s) for coverage peaks", flush=True)
    cov, lengths = accumulate_coverage(bams)
    if not cov:
        sys.exit("[derive-panel] no aligned reads in the sampled BAMs; nothing to derive.")
    peaks, thresh, gmax = call_peaks(cov, a.min_depth_frac, a.min_abs_depth,
                                     a.min_length, a.merge_gap, a.pad, lengths)
    print(f"[derive-panel] max depth {gmax:,}; adaptive threshold {thresh:,} "
          f"(max({a.min_abs_depth}, {a.min_depth_frac}*max))", flush=True)
    if not peaks:
        sys.exit("[derive-panel] no peak cleared the threshold; lower --min-depth-frac/--min-length.")
    print(f"[derive-panel] {len(peaks)} amplicon peak(s) over {len({q['contig'] for q in peaks})} contig(s)", flush=True)
    for q in peaks[:12]:
        print(f"[derive-panel]   {q['contig']}:{q['start']+1}-{q['end']} "
              f"({q['length']} bp, {q['max_depth']:,}x)", flush=True)
    if len(peaks) > 12:
        print(f"[derive-panel]   ... and {len(peaks) - 12} more", flush=True)

    write_panel_fasta(peaks, a.genome, a.out, a.prefix)
    print(f"[derive-panel] wrote panel FASTA -> {a.out}", flush=True)
    if a.out_bed:
        with open(a.out_bed, "w") as fh:
            for q in sorted(peaks, key=lambda x: (x["contig"], x["start"])):
                fh.write(f"{q['contig']}\t{q['start']}\t{q['end']}\t"
                         f"{a.prefix}\t{q['max_depth']}\n")
        print(f"[derive-panel] wrote BED -> {a.out_bed}", flush=True)
    if a.out_qc:
        with open(a.out_qc, "w") as fh:
            fh.write("amplicon\tcontig\tstart\tend\tlength\tmax_depth\tmean_depth\n")
            for i, q in enumerate(sorted(peaks, key=lambda x: -x["length"]), 1):
                fh.write(f"{a.prefix}_{i:02d}\t{q['contig']}\t{q['start']+1}\t{q['end']}\t"
                         f"{q['length']}\t{q['max_depth']}\t{q['mean_depth']}\n")
        print(f"[derive-panel] wrote QC -> {a.out_qc}", flush=True)


if __name__ == "__main__":
    main()
