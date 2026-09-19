# DEPRECATED — superseded for Fig. 3

This bundle timed Raccess Exact at K=1 (access_len=7 only), while its
Fig. 3 comparators (LinearRaccess, Raccess Window) use K=4, making all
three not directly comparable under one K.

`benchmark/results/raccess_exact_matched/8715122/` reruns Raccess Exact at
the matched K=4 on the same 25,300 sequences. Use that bundle instead for
Fig. 3's Raccess Exact curve. This bundle is still used elsewhere in the
manuscript wherever Raccess Exact serves as the ground-truth reference for
*accuracy* comparisons (Table 1, Fig. 2, Supplementary S1/S3), which is a
different quantity (probabilities) than this bundle's own runtime/memory
measurements and is unaffected by K.
