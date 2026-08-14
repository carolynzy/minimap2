# IS-aware rescue & liftover — edge cases and treatment

Companion to `IS_AWARE_RESCUE_PLAN.md`. Covers the edge cases agreed on afteriterating against real samples (`anhui_2013_124182`, `jiangsu_2013_103168`,`hebei_2013_031030`, and the 1355 bp family at H37Rv IS-6110 loci).

The work is split between two layers:

- **minimap2 rescue** (`align.c`): decides whether to substitute an IS-aware
  middle alignment for the anchor-bounded interior of one chain region.
- **liftover** (`liftover_coords.sh`): turns each BED record into a single
  reference coordinate using the chain that minimap2 produced.

Each layer has its own set of edge cases below, then a worked-example
section at the end.

---

## 1. Rescue alignment (minimap2 `align.c`)

### 1.1 Anchor search neighbourhood — no `is_window` padding

An anchor candidate seed must lie **outside the IS body** on **both**
coordinate systems. `is_window` is reserved for the nearby-IS detection pass
and is NOT used to pad the anchor exclusion zone.

```
                ref IS interval
                ┌──────────────┐
ref:  ─────────[is_st        is_en]───────────────
              ▲                  ▲
              │                  │
   left anchor:                  right anchor:
   seed.ref_en ≤ is_st           seed.ref_st ≥ is_en

                query IS body
                ┌──────────────┐
query: ────────[q_is_st     q_is_en]──────────────
              ▲                  ▲
              │                  │
   left anchor:                  right anchor:
   seed.q_en ≤ q_is_st           seed.q_st ≥ q_is_en
```

Both per-coordinate checks must hold. (The query-side check is enforced by
the order of operations — the query IS search runs after anchors are picked,
inside the anchor-bounded window, so the query IS is by construction
between the anchors.)

### 1.2 Chain-edge fallback

When a seed anchor is missing on one or both flanks (no seed satisfies the
neighbourhood rule within the current chain region), the missing side falls
back to the chain region's own edge.

```
   chain region edges (one mm_reg1_t):
   r0/q0 ─────────────────────────────── r1/q1
        ↑                               ↑
        │                               │
  if no left seed:                if no right seed:
   left anchor = (r0, q0)         right anchor = (r1, q1)
```

`mm_cigar_locate` is skipped for chain-edge anchors — the chain edge lies on
the chain's CIGAR by definition.

### 1.3 Rescue-window threshold (30 kb)

The threshold gates only the cases where at least one anchor is a chain-edge
fallback. With two real seed anchors, the rescue runs regardless of window
length.

```
   if (left_is_edge || right_is_edge):
       require max(re - rs, qe - qs) ≤ MM_IS_RESCUE_WIN_MAX (30 000 bp)
       — else log [is_rescue_pos_skip … reason=window_over_threshold]
```

A real-seed pair is trusted because each seed has passed `mm_cigar_locate`
on the base region's CIGAR and is therefore a genuine aligned pair in this
region. A chain-edge anchor is just the region's boundary point, so we
limit how much sequence the rescue is allowed to re-cover with it.

### 1.4 Flank check — clipping at chain edge, tri-state per flank

`mm_check_rescue_flanks` clips its left/right flank windows at the chain
region's source and reference boundaries (and also the sequence boundaries)
— beyond the chain edge there is no alignment evidence, so we don't
compare against it.

```
   nominal flank window length:    is_flank_len
   clipping floor (rescue-able):   is_flank_len / 2

   per flank, the state is:

     PASS   clipped length ≥ floor  AND  mismatches ≤ allowed_max
     HARD   clipped length ≥ floor  AND  mismatches > allowed_max
     SHORT  clipped length < floor   (anchor sits at the chain edge)

   accept rescue if:
     left=PASS  AND right=PASS
     left=PASS  AND right=SHORT
     left=SHORT AND right=PASS

   reject rescue if:
     either side is HARD          (real mismatch — alignment unreliable)
     both sides SHORT             (no flank evidence on either side)
```

A SHORT flank is acceptable because the 30 kb rescue-window gate already
guarantees the corresponding anchor is within 30 kb of its chain edge — the
flank is short precisely because the anchor IS the chain edge.

### 1.5 Single chain region — no cross-block rescue

Anchor search and flank check run inside one `mm_reg1_t` only. The rescue
never merges adjacent chain regions, never invents an in-between
alignment, never reads seeds from a different region.

```
   chain region A             chain region B
   ┌───────────────┐  (chain   ┌───────────────┐
   │  …flanks…  IS │   split)  │ flank…        │
   └───────────────┘    ✗      └───────────────┘
                  ▲                       ▲
                  │                       │
       rescue running on A cannot pick a seed
       from B (and vice-versa).
```

The chainer's decision to split is treated as evidence the in-between
sequence is genuinely uncertain. Cases where the chain splits at the IS
boundary (`hebei_031030`, the 1355 bp family) accept this — rescue logs
`no_positional_anchor` and steps aside, and the downstream layer handles
the residue.

### 1.6 Same-site vs different-site / opposite-orientation classification

When the rescue middle reaches a column where both a query-IS boundary
(`atQ`) and a reference-IS boundary (`atR`) coincide, `mm_flush_excised`
runs a forward DP between the two IS bodies and decides by identity:

```
  forward DP identity ≥ is_min_id:
      same orientation, same-site IS event
      → emit the body alignment (M/X/I/D) verbatim

  forward DP identity < is_min_id:
      different IS family,  OR  opposite orientation,  OR  novel insertion
      → emit  I[query_IS_len] + D[ref_IS_len]
        (deletion of the reference IS  +  insertion of the sample IS)
```

The identity threshold naturally separates same-orientation from
opposite-orientation: a reverse-complement IS does not forward-align with
high identity, so it falls into the `I + D` branch. No explicit orientation
flag is needed.

### 1.7 Two consecutive sample IS at the same ref-IS locus, both opposite

Treated as three independent events:

```
   reference:  ───[flank]──[ref IS, fwd]──[flank]───

   sample:     ───[flank]──[IS rev #1][IS rev #2]──[flank]───

   reported as:
     · deletion of the reference IS
     · insertion of sample IS #1
     · insertion of sample IS #2
```

`mm_flush_excised` falls to `I + D` for each opposite-orientation pair
because the forward DP fails the identity threshold.

---

## 2. Liftover (`liftover_coords.sh`)

### 2.1 No deletion-aware position shifting

After rescue, liftover should use the CIGAR structure produced by the flank
alignment directly. A flanking insertion or deletion does not move the IS event
to the edge of that flank gap.

Same-column sample/ref IS classification is decided inside `align.c`:

```
cleaned-window DP places sample_IS and ref_IS at same column:
    forward orientation -> report the ref IS / same-site event
    reverse orientation -> report sample I + ref D

cleaned-window DP places them at different columns:
    report the sample insertion and ref deletion at their separate columns
```

Large flanking insertions or deletions remain separate events. For example, a
forward sample IS aligned to a ref IS and flanked by large insertions at both
ends reports the same IS site plus two independent insertion events.

### 2.2 Records without a rescued same-column IS

If rescue does not produce a same-column sample/ref IS pair, liftover should
fall back to the ordinary flank/block position logic for that record. It should
not snap the site to a neighboring deletion just because an `I + D` pattern is
nearby.

---

## 3. Worked examples

### 3.1 `anhui_2013_124182` — rescue produces I+D, liftover lands at `qe`

Rescue ran successfully on the single chain region. The resulting CIGAR
around the IS event:

```
… M91859 [r3028099-3119958]  I1353  D564 [r3119958-3120522]
   M3 [r3120522-3120525]      D1353 [r3120525-3121878]
   M206 [r3121878-3122084] …
```

Chain block view (after transanno):

```
   block A                block B               block C
   ref […-3119958]+       ref [3120522-3120525]+ ref [3121878-…]+
                  ↑                      ↑
            qe_A=3119958          qe_B=3120525
                  ^
            same rescue column
```

Liftover output uses the rescued same-column `I + D` position directly;
nearby target-side deletions do not move the insertion site.
### 3.2 `jiangsu_2013_103168` consec — opposite-orientation IS

Pre-rescue chain:

```
   chain3 (+)             chain16 (−)            chain13 (+)
   src […─3705298]        src [3710089─3711444]  src [3711674─…]
   ref […─3710381]        ref [3710381─3711736]  ref [3711736─…]
```

`chain16` is the small (~1355 bp) reverse mini-alignment of the sample IS
body to the reference IS-6110 copy. Under the new rescue:

- Operating on `chain16` as the base region.
- Anchor search inside `chain16`'s seeds with no `is_window` padding:
  both flanks miss (chain16 covers the IS entirely).
- Chain-edge fallback on both sides (chain length = 1355 bp ≤ 30 kb).
- `mm_flush_excised`: forward DP of sample IS body vs reference IS body
  has LOW identity (reverse-complement orientations) → emits `I + D`.
- `chain16`'s CIGAR is rewritten as `I[1355] + D[1355]`.

Post-rescue, the chain delivered to the liftover effectively has:

```
   chain3 ── [I_sample] [D_ref_IS] ── chain13
                         ^
                   same rescue column
```

Liftover output uses the rescued `I + D` column directly. Consecutive
insertions at the same site may still be separated by a stacking layer,
but they are not shifted by the neighboring deletion.

### 3.3 `hebei_2013_031030` — same-orientation, rescue keeps `M`

Chain (one region, contains the entire ref IS-6110 copy inside a long M):

```
   one mm_reg1_t region (~485 kb):
   … M36423 [r3083203-3119626]
     D74    [r3119626-3119700]
     M2178  [r3119700-3121878]   ← ref IS-6110 at [3120525-3121878] sits here
     (chain split immediately after at r=3121878 ↔ r=3122084)
```

Rescue acts:

- Anchor search finds a left seed (nearest seed before ref 3120525, inside
  early `M2178`). No right seed exists past ref 3121878 inside this
  region.
- Right anchor falls back to the chain region's right edge at ref 3121878.
- Rescue window ≈ 1.35 kb, well under 30 kb.
- `mm_flush_excised` runs forward DP of sample IS body vs ref IS body —
  identity HIGH (genuinely same orientation, same sequence) → emit the
  M-style aligned CIGAR. No `I + D` produced.

The chain handed to liftover is essentially unchanged: still
`... M D[74] M2178 ...`. The same-column rescue decision keeps the
forward sample/ref IS as aligned reference IS sequence; the neighboring
`D74` remains a separate deletion and does not shift the IS site.

This removes the old trade-off where same-orientation cases hidden inside
a long M run could be pushed to a chain block edge.
than at the IS-copy boundary inside the block.

### 3.4 1355 bp family at H37Rv IS-6110 loci — all collapse to one rule

Eight samples (`anhui_122124`, `anhui_122154`, `anhui_124055`,
`fujian_132155`, `fujian_132314`, `jiangsu_103168`, `jilin_071174`, …)
all show the pre-rescue chain pattern of §3.2 around ref 3710381↔3711736.
Under the new rescue + liftover rules they all converge to the same
result: `qe_chain3 = 3710381`.

---

## 4. Constants and gates (quick reference)

| name | value | meaning |
|---|---|---|
| `MM_IS_RESCUE_WIN_MAX` | 30 000 bp | max rescue-window length when ≥ 1 anchor is chain-edge |
| `is_flank_len` | option (default 50) | nominal flank-check window length |
| `is_flank_len / 2` | derived | clipping floor for SHORT-flank acceptance |
| `is_window` | option (default 2 000) | nearby-IS detection padding; NOT used for anchor search |
| `is_min_id` | option (default 0.85) | identity threshold for same-site vs `I + D` |
| `is_min_flank_id` / `is_flank_max_snp` | options | strict flank-identity bounds in `mm_check_rescue_flanks` |
