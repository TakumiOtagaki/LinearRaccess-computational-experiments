# Benchmark manifests

One JSON manifest should describe each frozen experiment family. Each manifest
must contain an integer `schema_version` and should record:

- a stable experiment identifier and purpose;
- accession/version or a content hash for every input;
- sequence-selection and preprocessing rules;
- LinearRaccess parameters and the pinned software revision;
- comparator name, version, parameters, and installation source;
- repetitions, warm-up policy, timeout, and resource limits;
- planned result bundle and archival-data identifiers.

Do not edit a manifest after a paper result has been produced. Create a new
manifest/revision and retain the old provenance instead.

The retrospective manifest for bpRNA job `8661385` is stored with its audit
artifacts at `benchmark/results/bprna_audit_8661385/manifest.json`. It records
observed file hashes and limitations; it is not an execution-time frozen
software revision. The output metadata reports `source_revision=unknown`.

The fast/standard log-sum-exp performance experiment is frozen in
`benchmark/results/bprna_logsumexp_tradeoff_20260914/manifest.json` before execution.
It includes stratified sample selection, timing protocol and executable/script
hashes; the job runs a copied executable after checking those hashes.
