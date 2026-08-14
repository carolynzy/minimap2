# IS-Aware Flank Rescue Plan

## Goal

Fix distorted alignments around reference IS copies by re-aligning the non-IS sequence within a bounded window using only the IS information already encoded in minimap2's prior alignment. No sample-side IS detection, no template scanning, no IS-vs-IS DP, no off-by-one boundary search.

If IS rescue options are not provided, minimap2 behaves as standard upstream.

## Architecture

One code path. The rescue runs as a **post-pass** over already-built regions in `mm_rescue_positional_is`. Chain gap-filling in `mm_align1` is plain upstream ksw2 — no IS awareness, no `MM_SEED_IS_SKIP`, no `is_rescue_trigger`. The earlier "legacy" gap-fill rescue path has been removed.

## Inputs

```bash
--is-seq FILE              IS sequence or consensus FASTA
--is-ref-bed FILE          BED of known IS locations in the reference
--is-config FILE           JSON/key=value config (may list IS FASTA files and tweak global thresholds)
--is-window INT            ref-IS overlap neighborhood used by the rescue trigger and chain-edge fallback (default 2000)
--is-flank-len INT         nominal flank window length for the flank check (default 50)
--is-flank-max-snp INT     mismatch tolerance per flank window (default 1)
--is-min-flank-id FLOAT    flank identity threshold (default 0.95)
```

At least one of `--is-seq`, `--is-ref-bed`, or `--is-config` must be present to activate rescue.

The options `--is-min-id`, `--is-min-ins-len`, and `--is-min-cov` are no longer in the rescue decision path. They survive only for reference-IS auto-detection from FASTA when no BED is provided (see `index.c`).

## Reference IS Intervals

BED schema is simple `chrom start end`. Adjacent intervals are **not** merged — even consecutive ref ISs within a few bp are processed as separate intervals so each can be excised, classified, and reported on its own.

If `--is-ref-bed` is not provided but IS FASTA templates are, reference IS copies are auto-detected from the FASTA at index-load time (`mm_idx_is_ref_load_fasta` in `index.c`) and the interval list is built in memory.

## Rescue Trigger

For each reference IS interval, define:

```text
IS neighborhood = [IS_reference_start - is_window, IS_reference_end + is_window]
```

Rescue is attempted only when:

- the alignment region (`base->rs`/`base->re`) overlaps the neighborhood;
- positional anchors can be found inside the same chain region (seed-or-edge).

`is_window` is used **only** to test overlap. It is not used as padding for anchor search.

Rescue is constrained to the alignment region (chain region) that triggered it — it never spans block boundaries.

## Positional Anchor Selection (`mm_find_positional_is_anchors`)

Anchor search is restricted to seed positions within the chain region on the same query, reference, and strand.

**Search neighborhood — no `is_window` padding:**

```text
left  anchor candidates: seed.ref_en <= ref_IS_start
right anchor candidates: seed.ref_st >= ref_IS_end
```

**Anchor requirements:**

- same read, same reference contig, same strand;
- collinear in query and reference;
- not repetitive, tandem, ignored, or high-occurrence;
- seed length sufficient to act as an anchor;
- right anchor downstream of left anchor in query;
- both anchors lie within the same alignment region.

High-occurrence is rescue-anchor-specific. During minimap2's ordinary seed hit
emission in `map.c`, a seed hit is marked with `MM_SEED_HIGH_OCC` when its
reference occurrence count is greater than `MM_IS_ANCHOR_MAX_OCC = 10`.
Ordinary chaining can still use seeds according to minimap2's normal
`mid_occ`/`max_occ` logic; this flag only prevents such seeds from becoming IS
rescue anchors.

Anchor seeds must also avoid annotated reference IS/repeat intervals with
padding equal to `is_flank_len`. This means a seed just outside the focal IS is
still rejected if it overlaps a nearby annotated element or its padded flank.
For this to catch adjacent mobile elements such as IS1547, those intervals must
be present in `--is-ref-bed`; an IS6110-only BED only masks IS6110 copies.

Finally, a single nearest seed is not sufficient. Each candidate anchor must be
supported by a small local collinear seed cluster:

```text
MM_IS_ANCHOR_CLUSTER_MIN_SEEDS = 3
MM_IS_ANCHOR_CLUSTER_RADIUS    = 100 bp in reference and query
MM_IS_ANCHOR_CLUSTER_MAX_DIAG  = 25 bp
```

The 100 bp cluster radius is not the rescue window size. It is only the local
support radius around a candidate seed. The rescue window is still the interval
between the accepted left and right effective anchors.

Picks: greatest `ref_en` ≤ ref-IS-start for the left, smallest `ref_st` ≥ ref-IS-end (and downstream of left in query) for the right.

**Chain-edge fallback:** if a seed anchor is unavailable on either side, the alignment region's source/reference edge (`base->rs` / `base->re`) becomes the effective anchor on that side. A chain-edge fallback anchor skips the fixed-flank check.

**Single-anchor cases:** if only one qualified seed anchor is found, it is paired with the chain-region edge on the other side.

**Rescue-window 30 kb gate (`MM_IS_RESCUE_WIN_MAX = 30000`):** applied only when chain-edge fallback is used on at least one side. If the rescue window exceeds 30 kb in either source or reference coordinates with fallback in play, the rescue is skipped. When both anchors are real seed anchors, the 30 kb threshold is ignored.

## Flank Check (`mm_check_rescue_flanks`)

Per-flank tri-state PASS / SHORT / HARD with the floor at `is_flank_len / 2`. Flank windows are clipped at chain edges. Decision:

- both SHORT → reject
- any HARD → reject
- otherwise → accept

A chain-edge fallback anchor produces a legitimately SHORT flank that is paired with the other side's PASS.

## Prior-CIGAR Pair Detection (`mm_pairs_in_window`, `mm_locate_paired_sample_chunk`)

For each ref-IS interval that overlaps the rescue window:

1. Walk the prior alignment CIGAR over those reference columns.
2. Collect the sample positions touched by match ops (`M`, `=`, `X`).
3. If at least one M op covers any of those columns, the ref IS is **paired** with the sample chunk `[lo, hi]` — the minimum and maximum sample positions touched. The pair is recorded in parallel arrays.
4. If no M op covers any of the columns, the ref IS is **unpaired** (entirely a D op in the prior CIGAR).

There is no sample-side IS template scan anywhere in the rescue path.

## DP Input Construction

The rescue DP aligns a cleaned sample window against a cleaned reference window:

- `NR` = the rescue window's reference, with **every ref IS overlapping the window excised** (both paired and unpaired).
- `NQ` = the rescue window's sample, with **only paired sample chunks excised**. Unpaired sample IS bodies — novel sample insertions identified as large I ops in the prior CIGAR — remain in `NQ` and pass through the DP as I ops in the output.

`Bq` and `Br` boundary arrays record the cumulative offsets of each excision so the post-DP walk knows where to re-emit each excised IS.

## Rescue DP Parameters

The cleaned-window DP uses the **same scoring as the surrounding whole-genome alignment** — asm5 by default, or whichever preset / overrides the user supplied on the command line. No rescue-specific scoring overrides.

Concretely, `mm_build_is_rescue_cigar` calls `mm_align_pair` with the caller's `opt` and the caller's pre-built `mat`, with bandwidth `opt->bw_long * 1.5 + 1` (matching the gap-filling convention) and `opt->zdrop`.

Rationale: keeping rescue scoring identical to the host alignment avoids two parameter sets disagreeing about what counts as a match or what gap is too long. If a particular dataset needs tuning, change the global preset (or pass `-A/-B/-O/-E/-z`) and both the main alignment and the rescue see the change consistently.

## Per-Pair Classification (`mm_flush_excised`)

The post-DP walk uses `cq` (cursor in NQ) and `cr` (cursor in NR) and the `Bq` / `Br` boundary arrays. At each boundary:

- **Both `cq == Bq[jq+1]` and `cr == Br[kr+1]`** — a sample chunk and a ref IS arrive at the same column position after the DP. Run the orientation check:
  - **Same orientation** (forward identity ≥ reverse identity in a base-by-base comparison) → emit M of `min(ql, rl)` length, plus a length residual as a trailing small I or D for any chunk-length mismatch.
  - **Opposite orientation** → emit canonical `<ql>I` + `<rl>D` at that column. Any small indels minimap2 had encoded inside the bodies are not re-derived.
- **Only `cq` boundary hit** — paired sample chunk landed at a column with no ref-IS boundary → emit the sample chunk as one I.
- **Only `cr` boundary hit** — ref IS landed at a column with no sample-chunk boundary, or it's an unpaired ref IS → emit the ref IS as one D.

The orientation check applies whenever both boundaries coincide structurally, *regardless of whether the prior CIGAR had originally paired these specific bodies*. The DP's column placement is the arbiter.

The orientation comparison is implemented in `mm_pair_is_forward` — one pass through the sample chunk against the ref IS, counting matches forward and matches against the reverse complement, returning 1 if forward ≥ reverse.

If a same-site, same-orientation sample/ref IS pair has a deletion on one flank
(for example `D + sample_IS-ref_IS_alignment + M`), the rescue keeps the IS body
as aligned sequence (`M` plus any length residual). It does not reinterpret the
sample IS as a new insertion at the left edge of the flanking deletion. This is
intentional: the aligned IS body is treated as shared sample/reference sequence,
and the flanking deletion remains a separate flank difference.

## CIGAR Splice (`mm_splice_rescue_cigar`)

The mid-CIGAR built above is spliced into the prior region's CIGAR between the two effective anchors, preserving everything outside the rescue window untouched.

## Multi-IS Within One Rescue Window

Every ref IS that overlaps the window is in the excised set. Paired ones come out with their sample partners on both sides of the DP; unpaired ones come out only on the reference side. The DP input is IS-free on the reference side regardless of how many IS bodies were in the window. Unpaired sample IS bodies stay in `NQ` and pass through the DP as I ops.

After the DP, each excised pair is classified by the rule above, independently. Adjacent ref ISs are treated as separate intervals throughout — no joint classification.

## What's Gone Versus Earlier Designs

- The legacy gap-fill rescue (`is_rescue_legacy_accept`) and its `mm_is_rescue_prepass` / `MM_SEED_IS_SKIP` machinery.
- The template scanner for sample IS detection (`mm_query_find_is_segs`) and its supporting helpers (`mm_query_has_is_like_seq`, `mm_query_find_is_like_seq`, `mm_project_query_is_search_window`).
- The IS-vs-IS DP and `mm_cigar_identity` identity check that decided same-site vs different-site by sequence comparison.
- The `is_min_id` decision dependency (the option still exists for ref-IS auto-detection in `index.c` only).
- The 20-bp adjacent ref-IS merge.
- The old `del_report_mode == 1` CIGAR reorder logic.
- The adjacent-deletion expansion/reorder path for rescued IS bodies.
- The off-by-one boundary bug — no phase-scanning sample sequence against an IS template anywhere in the rescue path.
- `mm_ref_is_in_window`, `mm_last_cigar_op`, `mm_pop_last_cigar_if_op`, `mm_last_cigar_len` — all now dead, removed.

## Flanking Gap Behavior

Flanking insertions and deletions are not allowed to shift the rescued IS
classification or event position. `mm_flush_excised` first asks whether the
sample IS and ref IS arrive at the same cleaned-window DP column. If yes, it
checks orientation: same orientation emits the aligned ref IS as `M` plus any
length residual; opposite orientation emits canonical `I + D`. If no, the sample
IS and ref IS are emitted independently at their different columns.

Surrounding I/D ops remain separate flank differences. For example, if a forward
sample IS aligns to a ref IS and has large insertions at both ends, the ref IS is
still reported as the same insertion site, and the two flanking insertions remain
two separate insertion events. The same rule applies to flanking deletions.

## Debug Output

At `-v 3`, rescue logging includes:

- `[M::is_rescue_pos_skip]` with reason `no_positional_anchor` or `window_over_threshold(30000)`;
- `[M::is_rescue_pos_reject]` with reason `flank_mismatch` (and the per-flank mismatch counts) or `cigar_splice_failed`;
- `[M::is_rescue_pos_replace]` with the ref-IS interval, left/right seed indices, block ref/query ranges, window ref/query ranges, flank length and per-flank mismatch counts, score, and cigar op count.

The log makes clear which path each ref-IS attempt took.

## Most Important Rule

Do not perform free local realignment over the whole read.

Rescue must be:

```text
triggered by an alignment region overlapping a ref-IS neighborhood
bounded by fixed effective anchors (seed or chain-edge) inside one chain region
middle-only
flank-preserving (with chain-edge clipping when flanks would cross a chain boundary)
DP-aligned on IS-excised input with a dedicated parameter set
classified by structural column coincidence after the DP, plus a single
forward-vs-reverse orientation comparison for same-column pairs
gated at 30 kb only when chain-edge fallback is used on at least one side
```

See `IS_RESCUE_EDGE_CASES.md` for worked examples and per-rule illustrations.
