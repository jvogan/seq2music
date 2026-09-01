# Seq2Music

Seq2Music is a local-first Codex plugin for converting biological sequences and structure-derived features into audible, visible, and machine-readable music. It can also decode its self-describing MIDI, PCM WAV, and MusicXML carriers back to normalized protein, DNA, or RNA sequences, and transcribe arbitrary supported audio or scores into a selected sequence alphabet.

![Seq2Music protein-note logo](plugins/seq2music/assets/seq2music-protein-note-logo.png)

## Highlights

- Lossless normalized-sequence round trips through self-describing MIDI, PCM WAV, and MusicXML.
- Explicitly lossy transcription from metadata-free PCM WAV and plain MusicXML.
- Protein, DNA, and RNA alphabets.
- Audible reference/variant diffs, alignment summaries, and residue-order structure traces.
- Accessible SVG notation, HTML playback, CSV event ledgers, and SHA-256 manifests.
- Offline implementation using Python's standard library.
- Streamed WAV rendering with adjustable resource budgets and no preview-duration ceiling.

## Install in Codex

Clone or download this repository, open a terminal in its root, and run:

```bash
codex plugin marketplace add .
codex plugin add seq2music@seq2music
```

Start a new Codex task after installation so Codex discovers the plugin skill. No account, API key, or network service is required by Seq2Music itself.

To remove it:

```bash
codex plugin remove seq2music@seq2music
codex plugin marketplace remove seq2music
```

## Run the CLI directly

Python 3.10 or newer is required. From the repository root:

```bash
python3 plugins/seq2music/scripts/seq2music.py encode \
  --input plugins/seq2music/assets/examples/demo-protein.fasta \
  --kind protein \
  --out seq2music-output
```

The command creates a deterministic bundle containing MIDI, PCM WAV, MusicXML, accessible SVG notation, HTML, a CSV event ledger, a text summary, and a JSON run manifest.

## Data handling

All processing is local. Exact MIDI remains recoverable through its ordered notes and interpretation metadata; exact WAV and MusicXML also carry embedded sequence metadata. Sharing any exact carrier can therefore share the underlying sequence. See [PRIVACY.md](PRIVACY.md) for the complete data-handling summary.

Seq2Music supports macOS, Linux, and Windows. On platforms without secure directory-descriptor operations, new output works normally while `--force` replacement is disabled; choose a new output path instead.

## Verify the repository

```bash
python3 scripts/public_release_audit.py
python3 -m unittest discover -s plugins/seq2music/tests -v
```

The release audit checks the portable marketplace layout, manifest assets, common credential patterns, private absolute paths, symlinks, generated caches, and oversized files.

## Repository layout

```text
.
├── .agents/plugins/marketplace.json
├── plugins/seq2music/
│   ├── .codex-plugin/plugin.json
│   ├── assets/
│   ├── scripts/seq2music.py
│   ├── skills/sonify-biomolecules/
│   └── tests/
└── scripts/public_release_audit.py
```

Detailed mappings, codec statuses, supported inputs, and scientific references are documented inside the plugin at [`plugins/seq2music/README.md`](plugins/seq2music/README.md) and [`plugins/seq2music/skills/sonify-biomolecules/references/`](plugins/seq2music/skills/sonify-biomolecules/references/).

Release changes are recorded in [CHANGELOG.md](CHANGELOG.md).

## License

MIT © 2026 Jacob Vogan. See [LICENSE](LICENSE).
