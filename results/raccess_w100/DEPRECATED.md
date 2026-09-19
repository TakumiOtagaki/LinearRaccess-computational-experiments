# DEPRECATED — superseded for Fig. 3

This bundle timed Raccess Window at K=1 (access_len=7 only), while
LinearRaccess's comparator (`bprna_runs/8666393/`) used K=4, making the two
not directly comparable for a runtime constant-factor claim -- flagged by an
independent review.

`benchmark/results/bprna_matched_runtime/8714351/` reruns Raccess Window at
the matched K=4, executable-only timing, on the same 26,063 sequences. Use
that bundle instead for Fig. 3 / Results S3.3 runtime and memory claims.
