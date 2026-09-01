# Contributing

Contributions are welcome through ordinary GitHub issues and pull requests once this repository is published.

## Development checks

Use Python 3.10 or newer and run these commands from the repository root:

```bash
python3 scripts/public_release_audit.py
python3 -m unittest discover -s plugins/seq2music/tests -v
```

Keep the implementation offline and standard-library-only unless a dependency change is discussed first. Preserve deterministic output, bounded parsing, explicit exact/lossy statuses, accessible non-audio artifacts, and the public codec/mapping contracts. Include tests for behavioral changes.

Do not commit generated audio bundles, caches, credentials, private datasets, or machine-specific paths.
