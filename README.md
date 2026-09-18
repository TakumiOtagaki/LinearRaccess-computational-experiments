# LinearRaccess computational experiments

Benchmark protocols, scripts, provenance manifests, and compact result
bundles supporting the tables and figures in:

> Otagaki, T. and Sato, K. "LinearRaccess: linear-time approximation of RNA
> accessibility without base-pair span limits." Submitted to *IPSJ
> Transactions on Bioinformatics*.

This repository is also used as a Git submodule (mounted at `benchmark/`) of
the manuscript repository,
[`TakumiOtagaki/LinearRaccess-paper`](https://github.com/TakumiOtagaki/LinearRaccess-paper),
which is where the pinned LinearRaccess software revision
(`software.lock.json`) and the manuscript sources live.

## Layout

- `scripts/`: benchmark drivers, UGE job scripts, and figure-generation code.
- `manifests/`: experiment manifests (settings, sequence lists).
- `data/`: compact input data (sequence lists, small FASTA files).
- `results/`: result bundles. Each bundle contains at minimum
  `provenance.json` and `summary.json`; raw per-sequence measurements are
  kept out of Git where they are too large, and are regenerable from the
  provenance record.
- `tests/`: unit tests for the analysis scripts.

## Reproducing a result

Each `results/<experiment>/provenance.json` records the exact command line,
software revision (matching `software.lock.json` in the paper repository),
input file hashes, and settings used to produce that bundle's
`summary.json`. Re-running the recorded command against a checkout of
LinearRaccess at the pinned revision should reproduce the same summary
statistics.

## Relationship to LinearRaccess and Raccess

- The LinearRaccess software itself (implementation, CLI, tests) lives in
  [`TakumiOtagaki/LinearRaccess`](https://github.com/TakumiOtagaki/LinearRaccess),
  referenced here only by its pinned commit/tag, never vendored.
- Raccess is an independently installed external comparator/reference
  program; its source tree or binaries are never committed here.

---

Historical scope note (from this workspace's original benchmark protocol):
paper evidence is divided into two required experiment families --
(1) correctness and approximation quality against an explicitly identified
reference calculation, and (2) runtime and memory scaling -- with sampling
comparisons and biological case studies treated as optional additions
following the same provenance rules. Exploratory results are not paper
evidence: before a result enters a figure, its manifest, software revision,
input hashes, command-line options, hardware, compiler, repetitions, and
output hashes must be frozen.
