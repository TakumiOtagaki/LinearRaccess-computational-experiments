# Result bundles

Each subdirectory is an immutable result bundle and must contain:

- `provenance.json`: software/input hashes, invocation, environment, hardware,
  timestamps, and links to externally archived raw data;
- `summary.json`: machine-readable values consumed by analysis and figures.

Additional compact tables and logs may live beside these files. Do not commit
large transient output merely to make Git a data archive; deposit it in Zenodo
or another durable repository and record its DOI and cryptographic hashes.
