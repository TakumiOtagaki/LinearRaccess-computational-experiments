# Benchmark inputs

Only compact, redistributable inputs belong here. Record the upstream database,
accession and version, retrieval date, original-file hash, preprocessing command,
and derived-file hash in the corresponding manifest.

For large or redistribution-restricted inputs, store only a retrieval recipe and
expected hash.

## Algorithm-validation suite

`algorithm_validation_v1.fa` contains six fixed synthetic RNA sequences at
lengths 60, 90, and 120 nt. They were generated with Python 3's
`random.Random(seed).choice("ACGU")`, using the seed encoded in each FASTA
identifier. The file is intended for dependency-free algorithmic regression
and approximation-envelope experiments. It is not a biological sample and
must not be used to support a biological claim.

The explicit sequences are committed so that reproduction does not depend on
the behavior of a particular Python version's pseudo-random generator.

## Long-tier input: NC_045512.2

`NC_045512.2.fa` is the NCBI Nucleotide (nuccore) record for SARS-CoV-2 isolate
Wuhan-Hu-1, complete genome, retrieved 2026-09-06 via E-utilities `efetch`.

- retrieval URL:
  `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text`
- original-file SHA-256:
  `0891c00cc2d503f58e0c600ce5ee64e9c64b939525c70d5486bc46cbc7cb073e`
- normalization command: strip FASTA header/whitespace, uppercase, replace `T`
  by `U` (identical to `benchmark_common.read_single_fasta`)
- normalized length: 29,903 nt
- normalized-sequence SHA-256:
  `d221b74def202fe4a77e14ce8695bfa5cbc7d5056ce4c8a98dafee8c2e0ab344`

This is a biological input used only for the fixed-parameter scaling experiment;
it is not evidence for a functional biological claim.

## Sampling-comparison suite

Four small accessioned RNA sequences, retrieved 2026-09-06 via E-utilities
`efetch` (NCBI `nuccore` unless noted as `structure`). They are used for the
pre-registered sampling-variance comparison (`docs/sampling_compare_plan.md`)
and are not evidence for a functional biological claim.

| accession | file | length (nt) | original-file SHA-256 | normalized-sequence SHA-256 |
|---|---|---|---|---|
| M10263.1 | `M10263_trna_phe.fa` | 76 | `9a1a8bee051fbe9816c3a7dff9f8cf9bb505784b6736eb94c44f633e3bc08841` | `28890de0dd2033379fabebc7165a219d6092ad9bbadc77f4c7528be78963ca0b` |
| 4KQY_A | `4KQY_A_sam1_riboswitch.fa` | 119 | `2f49524d78e747696bd6aa6adb3010ea199b21737ac9b4d5dc945ab499a65da4` | `6166f3c7310e9634a7daf92e8447d17de512c55e24005a52b020926eb404fa19` |
| HG323429.1 | `HG323429_srp_rna.fa` | 270 | `68e9fecc3ee7a7a928cadacf27bd9ccbabc672def02bf5eae0a55ab711458ff1` | `ccd6577793590e0e38a64f04be1d1e94c377804e724111f115648fa044981b13` |
| 7UO1_B | `7UO1_B_rnase_p_rna.fa` | 373 | `43dd5c0dcf95213dad23eb186a68c85760c90dafefa154583ce88eafafae3589` | `7c5ee6af2f8d3996e548889990796b82b56be931d619e0dd199a82ed5a2d4509` |

Retrieval note: the two nucleotide accessions use
`efetch.fcgi?db=nuccore`, the two structure accession IDs use
`efetch.fcgi?db=structure` returning the chain FASTA. E-utilities rate limit
(~3 requests/min without an API key) applies.
