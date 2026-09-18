# DEPRECATED — do not use as evidence

This bundle (production job `8661385`) is **not** used by the manuscript.

- `validation_passed=false`: 7,746/26,063 sequences have 43,775 probability
  values outside `[0, 1]` (up to `1.0009`). See `docs/bprna_audit_8661385.md`.
- Output metadata records `source_revision=unknown`, so the executing
  LinearRaccess commit cannot be confirmed.
- It predates the structural-profile normalization fix described in
  `docs/bprna_maintenance_20260914.md`.

The current, validated, current-revision (`892c93c8...`) full bpRNA run is
`benchmark/results/bprna_runs/8666393/` (26,063/26,063 succeeded, per-batch
manifests confirm the revision; see its `provenance.json`). Use that instead.
