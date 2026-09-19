# DEPRECATED — superseded for Fig. 3B

This bundle's memory measurements were used for the original Fig. 3B.
`benchmark/results/bprna_matched_runtime/8714351/` reruns the same
26,063-sequence LinearRaccess memory measurement under the matched-condition
protocol (see its own `provenance.json`) and additionally covers Raccess
Window's memory at the same K=4, which this bundle does not.

The memory values themselves are consistent between the two bundles (peak
RSS is not sensitive to the K/timing-boundary confound that affected
runtime); this bundle is not numerically wrong, just superseded as the
canonical source.
