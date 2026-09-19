# DEPRECATED — do not use as evidence for Fig. 3 / Results S3.3

This bundle's `elapsed_sec` timings are not used in the current manuscript.

- Timed with Python-side output validation and SHA-256 hashing included in
  `elapsed_sec` (`benchmark/scripts/run_bprna_chunk.py`), not the executable
  alone.
- Compared against `raccess_w100/` (K=1 access_len), while this bundle uses
  K=4, so the two were not directly comparable for a constant-factor claim.

An independent review flagged the resulting "roughly 6-fold faster" claim as
confounded. `benchmark/results/bprna_matched_runtime/8714351/` reruns both
tools on the same 26,063 sequences with matched conditions (same K=4,
executable-only `/usr/bin/time -v` timing, sequential execution). Use that
bundle instead for runtime figures/claims.

The accessibility *values* this bundle's LinRacc runs produced are not
affected (only the timing methodology is); this bundle is not known to be
numerically wrong, just timed inconsistently with its comparator.
