# kmerstash — laptop-portable DNA k-mer presence screening

A small, self-contained tool that answers one question fast and locally:

> **Is sequence X — a pathogen marker, an antimicrobial-resistance (AMR) gene, a
> contaminant — present in this sample of sequencing reads?**

No cluster, no GPU, no cloud. Pure-Rust, std-only, one binary. It is the
bioinformatics sibling of the COBOLMM trigram search engines: the same idea
(index the rare k-grams, screen in linear time) carried into a 4-letter alphabet.

```
  reference                   k-mers     found  containment   verdict
  ------------------------------------------------------------------------
  marker00                      1480      1474       99.6% ####################  PRESENT  ⚠
  marker04                      1480         0        0.0% ····················  absent
  ...
  ⚠ ALERT: 4 reference(s) present in sample (containment ≥ 50%)
```

## Why this works — and where it sits in the pipeline

Real aligners (minimap2, BWA, Kraken, mash) are **seed-and-extend**: a k-mer index
first finds *candidate* matches in time proportional to the read length, then a
slow exact stage (alignment / classification) runs only on those candidates.
`kmerstash` is exactly that **seed/screen stage** — and for the *presence/absence*
question, screening alone is the whole answer; you don't need the extend step.

The trick that makes it linear-time is the same as trigram search: a query of
length L has ~L k-mers, each a single resident lookup, so screening is O(L)
regardless of how big the reference panel or sample is.

## The engine

The trigram engine folds 3-grams into 37³ buckets and verifies candidates. For DNA
two things change, and that is the whole engine:

1. **4-letter alphabet → exact integer keys.** A k-mer of length `k ≤ 32` packs
   into a `u64` at **2 bits/base**. No hashing, no strings — the k-mer *is* its
   key. Membership in a sorted, deduped `Vec<u64>` (8 bytes/k-mer, binary search)
   is **exact**, so unlike folded trigram buckets there are no false positives and
   no separate verify stage.
2. **Double strand → canonical k-mers.** A read can come off either strand, so the
   same physical k-mer appears as a sequence or its reverse-complement. We collapse
   both to `min(forward, revcomp)`. (The 2-bit encoding is chosen so a base's
   complement is just `3 - code`.) Skip this and you miss ~half your matches.

Non-ACGT bases (`N`, gaps, ambiguity codes) reset the rolling window, so noise can
never fabricate a key.

## Robustness — the "one clean hit is enough" principle, measured

Sequencing reads have errors, and OCR taught us the relevant lesson: a screen
doesn't need to be perfect, it needs *enough* clean signal to clear a threshold.
A 150 bp read has ~130 k-mers; an error damages only the `k` k-mers spanning it.
Across typical coverage, plenty of clean k-mers survive:

| sequencing error | whole-sample containment of a present marker |
|---|---|
| 1% | 99.9% |
| 5% | 98.4% |
| 10% | 96.7% |
| 20% (brutal) | 54% — still clearly present |

So the *presence* verdict is stable far past the error rate of real instruments.

## Speed (single core, laptop)

~**35 Mbp/s** to build the sample sketch; screening a reference against it takes
**microseconds**. A few-million-read sample sketches in well under a second.

## Usage

```bash

# real use:
kmerstash gen    --out data            # or bring your own FASTA
kmerstash screen --ref panel.fa --sample reads.fa --k 21 --threshold 0.5
kmerstash screen --ref panel.fa --sample reads.fa --csv     # machine-readable

```

`panel.fa` is your reference panel (one record per marker / gene / organism).
`reads.fa` is the sample (FASTA; reads or contigs). Containment ≥ `--threshold`
fires an ALERT.

## Scope / boundaries

- **Presence/absence and classification** are the sweet spot. For *where exactly*
  a read aligns and *with what mutations*, you'd add the extend (alignment) stage —
  this tool is the fast screen in front of it, not a replacement for it.
- The demo data is **synthetic** (a seeded generator, so every figure is
  reproducible and the tool runs fully offline). The mechanism is real; to screen
  real samples, point `--ref` at a curated panel (e.g. NCBI RefSeq markers, or the
  CARD AMR-gene database) and `--sample` at your reads. Nothing else changes.
- This is a **lab / reference implementation**, in the spirit of the `tiny_pointers`
  and COBOLMM search labs — built to demonstrate and benchmark the engine, with a
  notebook (`kmerstash_demo.ipynb`) that explains and plots it end to end.

