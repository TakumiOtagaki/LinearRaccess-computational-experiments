# DEPRECATED — superseded by the full-dataset matched run

This 408-sequence subsample first identified and quantified the K=4-vs-K=1,
unequal-timing-boundary confound in the original Fig. 3 (~4.1-fold vs. the
confounded ~6-fold). It is no longer cited by the manuscript.

`benchmark/results/bprna_matched_runtime/8714351/` applies the same matched
protocol to the full 26,063-sequence dataset and is now the source for
Fig. 3 and Results S3.3 (fold ratio: 3.96). This bundle's own reproduction
script (`run_matched_runtime_comparison.py`) also hardcoded a scratch
workspace path rather than generating its own sample; `run_bprna_chunk_matched.py`
(used for the superseding bundle) does not have that problem.
