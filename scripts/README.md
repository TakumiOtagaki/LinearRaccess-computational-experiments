# Benchmark scripts

## bpRNA on TSUBAME (UGE)

Use `uge/bprna_chunk_job.sh` for the pre-split dataset: 26,063 sequences
in `benchmark/data/bprna_chunks/chunk_01.fa` through `chunk_27.fa`.
The script requires the executable `vendor/LinearRaccess/build-raccess/LinRacc`
built from a clean revision with `--build-info` support. The default repository
path is shown below; `BPRNA_REPO_ROOT` can override it.

```bash
cd /gs/bs/tga-satolab-gtex/otgk/LinearRaccess-paper/benchmark/scripts/uge
mkdir -p logs
qsub -g tga-satolab-gtex bprna_chunk_job.sh
```

Use your allocated accounting group if it differs from `tga-satolab-gtex`.
For an initial single-chunk run, add `-t 1-1` to the `qsub` command.
The full array has 27 tasks (966 sequences each, 947 in the last),
at most 10 concurrent tasks, and 10 sequence processes per task.
Each task requests `cpu_16=1` and a four-hour wall time.
Settings are `engine=lincapr`, Turner 2004, beam 100, loop caps 30,
and accessibility lengths 1, 3, 7, and 10.

New results go to `benchmark/results/bprna_runs/<JOB_ID>/batch_<task>/`.
Each batch contains `summary.tsv`, detailed `run_summary.json`, `manifest.json`,
sequence outputs and a `logs/` directory. `.inputs/` retains the actual input
chunk and executable snapshot. Existing batch directories are rejected.
The historical `bprna_final/` results are retained unchanged.

The job exits nonzero if any sequence fails, times out, produces no/partial
output, has nonfinite probabilities, or has mismatched identity/metadata.
The manifest records the executable SHA-256 and its embedded build revision,
input SHA-256, settings, concurrency and timeout. Unknown/dirty build revisions
are rejected. Finite approximation errors outside [0,1] are counted in
`run_summary.json`, with a diagnostic tolerance of `1e-12`, but are not clipped
or treated as execution failures.

Overrides can be passed with `qsub -v`: `BPRNA_CHUNK_DIR`, `BPRNA_OUTDIR`,
`BPRNA_LINRACC`, `BPRNA_PARALLEL` (default 10), `BPRNA_TIMEOUT` (per-sequence
seconds, default 3600), and `BPRNA_NO_FAST_LOGSUMEXP=1` for the standard path.
`run_bprna_chunk.py` is the task runner; the shell wrapper propagates its exit
status with `exec`. Elapsed times now use a high-resolution wall clock, but
include output validation and concurrent execution; they are not standalone
performance benchmarks.

The `bprna_array_job.sh` (one sequence per task), `bprna_batch_job.sh`,
and `bprna_parallel*.sh` scripts are historical alternative implementations.
The failure/provenance changes above apply to `bprna_chunk_job.sh`.

### Audit a completed bpRNA run

From the repository root, run the standard-library-only audit:

```bash
python3 benchmark/scripts/audit_bprna.py --job-id 8661385 \
  --out benchmark/results/bprna_audit_8661385
```

The output directory must not exist. For a repeat audit, choose a new directory;
keep frozen audit artifacts intact. The command requires the downloaded source
FASTA, filtered FASTA, all chunk outputs and sequence logs, the current local
binary/source checkout, and `qacct` records for the specified job. It audits
`benchmark/results/bprna_final`; `--job-id` selects accounting records, not a
different result directory. Archive the original outputs before rerunning jobs.

The audit checks preprocessing and chunk reconstruction, summary ID coverage,
output sequence identity, complete zero-based positions and valid windows,
finite probabilities in [0, 1], and expected output metadata. It records
structural-profile `prob_sum` warnings separately from probability errors.
`summary.json` contains distributions and validation findings;
`per_sequence.tsv` contains lengths, integer-resolution elapsed times and warning
counts; `issues.tsv` lists per-sequence failures. Global failures are recorded
in `summary.json`. A failed audit returns nonzero. `qacct.txt` retains accounting
evidence, and `manifest.json` records input hashes and retrospective provenance.
If accounting must be collected separately, save `qacct -j 8661385` output and
pass `--accounting-file /path/to/qacct.txt`; the audit verifies job and task IDs.
The recorded production run has out-of-range probabilities, so its audit is
expected to return nonzero. See [the audit report](../../docs/bprna_audit_8661385.md).
Current binary/source hashes do not establish the executable used historically
when output metadata says `source_revision=unknown`.

Elapsed times measure wall time under concurrent execution. Zero seconds means
below timestamp resolution, and summing these values does not give job wall time
or CPU time. Per-sequence peak RSS was not instrumented. Task-level `maxvmem`
from UGE is virtual memory and must not be presented as per-sequence peak RSS.

Regression checks for the output parser:

```bash
python3 -m unittest discover -s benchmark/tests -p test_bprna_audit.py -v
```

For the follow-up fast/exact log-sum-exp comparison on deterministically selected
bpRNA examples, run `python3 benchmark/scripts/diagnose_bprna_logsumexp.py --out
/tmp/bprna-logsumexp-repeat` from the repository root, using a new output directory.
See [the diagnostic report](../../docs/bprna_logsumexp_diagnostic.md) for the
observations and limits of this small-sample comparison.

For the stratified paired speed/accuracy experiment (24 random sequences plus
3 diagnostic examples), see [the experiment report](../../docs/bprna_logsumexp_tradeoff.md).
`run_bprna_logsumexp_tradeoff.py prepare --out NEW_BUNDLE` freezes inputs and a
binary copy; `run ... --out NEW_BUNDLE` executes four pinned workers on allocated
compute resources. The prepared 2026-09-14 bundle is launched from the repository
root with `qsub -g tga-satolab-gtex benchmark/scripts/uge/bprna_logsumexp_tradeoff.sh`.
Completed bundles must not be rerun or overwritten.

## Independent Raccess comparison

`run_raccess_compare.py` refuses an existing output bundle. It writes the
aggregate summary to `summary.json` at the bundle root and retains the same
summary under `analysis/summary.json` for compatibility. A missing comparison
cell or nonfinite opening energy causes a nonzero exit after diagnostics and
provenance have been written. `comparison_complete` indicates finite coverage,
not agreement within a scientific error tolerance. Zero probabilities require
separate interpretation because their opening energies are infinite.

Before using a bundle as E1 evidence, verify the model settings, source and
binary provenance, external-reference modifications, and numerical tolerances.
The runner does not yet enforce all of these publication requirements.

Scripts in this directory run independently of the historical
`LinearRaccessCapR` tree. They receive the LinearRaccess executable and checkout
as explicit arguments and consume versioned manifests.

## Scaling runner

After building the pinned LinearRaccess checkout, run the synthetic smoke tier:

```bash
python3 benchmark/scripts/run_scaling.py \
  --manifest benchmark/manifests/scaling_smoke.json \
  --tier smoke \
  --software-dir /path/to/LinearRaccess \
  --linracc /path/to/LinearRaccess/build/make/LinRacc \
  --output-dir /tmp/linearraccess-scaling-smoke
```

The runner refuses dirty paper/software repositories and an existing output
directory. It checks the software commit against `software.lock.json`, verifies
both raw-file and normalized-sequence hashes, checks metadata emitted by the
binary, and creates an immutable-style result bundle containing:

- `manifest.json`;
- `provenance.json`;
- `raw.tsv`;
- `summary.tsv`;
- `summary.json`.

`scaling_smoke.json` and its synthetic sequence test plumbing only; they are
not paper evidence. Frozen paper/HPC manifests must use reviewed inputs,
production parameters, and at least three repetitions after one warm-up.

## Deterministic work counters

Use the same manifest to audit inside-pass recurrence work independently of
wall-clock timing:

```bash
python3 benchmark/scripts/run_work_scaling.py \
  --manifest benchmark/manifests/scaling_smoke.json \
  --tier smoke \
  --software-dir /path/to/LinearRaccess \
  --linracc /path/to/LinearRaccess/build/make/LinRacc \
  --output-dir /tmp/linearraccess-work-smoke
```

This runner checks every post-prune DP table against the beam width and checks
the observed `MB <- M1 + M2` candidate count against `N * beam^2`. The counters
explain finite-length timing behavior; fitted counter slopes are descriptive
and are not themselves the complexity proof.

## Algorithmic approximation envelope

The dependency-free algorithm-validation runner compares controlled
approximations with the same bundled energy model in unpruned exact mode:

```bash
python3 benchmark/scripts/run_algorithm_validation.py \
  --manifest benchmark/manifests/algorithm_validation_smoke_v2.json \
  --software-dir /path/to/LinearRaccess \
  --linracc /path/to/LinearRaccess/build/LinRacc \
  --output-dir /new/path/algorithm-validation-smoke-v2-fixed
```

The current smoke suite varies beam width, each loop cap separately, the
log-sum-exp implementation, and their production-default combination. It
reports absolute and thresholded relative errors, rank correlation,
probability-bound violations, and nested-window monotonicity violations.
Per-position reference and observed probabilities are retained in
`comparisons.tsv`; `raw.tsv` contains per-sequence/configuration diagnostics,
and `summary.tsv`/`summary.json` contain aggregate statistics.

This is an internal same-model test on synthetic RNA. It verifies algorithmic
behavior but is neither an independent Raccess-oracle comparison nor evidence
about a biological population. A paper evidence bundle must use a separately
reviewed manifest and retain this distinction.

`algorithm_validation_smoke.json` is retained as the manifest used to discover
the `c_hairpin` posterior inconsistency in revision `3a0d68a`.
`algorithm_validation_smoke_v2.json` pins the corrected experiment identity
used with revision `51ed3c3` or a reviewed descendant.
