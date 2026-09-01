# Seq2Music

Seq2Music is an offline Codex plugin for deterministic biological data sonification and music-to-sequence transcription. Its default single-sequence MIDI, PCM WAV container, and MusicXML score are lossless, self-describing encodings of the normalized protein, DNA, or RNA sequence. It also converts arbitrary uncompressed WAV or plain MusicXML into an explicitly lossy sequence transcription, defaulting to the protein alphabet unless the user selects DNA or RNA. Every rendering includes MIDI, WAV, accessible grand-staff pitch/event-grid SVG notation, MusicXML, HTML, a CSV event ledger, text summary, and SHA-256 provenance manifest.

It requires Python 3.10 or newer, uses only the standard library, and makes no network calls. PDB is read using fixed columns. mmCIF support is a deliberately bounded `_atom_site` subset; unsupported syntax fails explicitly rather than invoking a network dependency.

## Quick start

```bash
RUN_DIR="$(python3 scripts/seq2music.py encode --input assets/examples/demo-protein.fasta --kind protein --out seq2music-output)"
python3 scripts/seq2music.py decode --midi "$RUN_DIR/demo-protein.mid" --manifest "$RUN_DIR/demo-protein.run.json" --out recovered
python3 scripts/seq2music.py audio-decode --audio "$RUN_DIR/demo-protein.wav" --out recovered-from-audio
python3 scripts/seq2music.py score-decode --score "$RUN_DIR/demo-protein.musicxml" --out recovered-from-score
python3 scripts/seq2music.py diff --reference assets/examples/demo-protein.fasta --variant assets/examples/demo-protein-variant.fasta --kind protein --out seq2music-output
python3 scripts/seq2music.py msa --input assets/examples/demo-alignment.fasta --kind protein --out seq2music-output
python3 scripts/seq2music.py structure --input assets/examples/demo-structure.pdb --model 1 --chains A --out seq2music-output
```

PowerShell uses the same printed-directory workflow: `$runDir = py -3 scripts/seq2music.py encode --input assets/examples/demo-protein.fasta --kind protein --out seq2music-output`, followed by `py -3 scripts/seq2music.py decode --midi "$runDir/demo-protein.mid" --manifest "$runDir/demo-protein.run.json" --out recovered`.

When the plugin is installed in Codex, direct generated output to the active workspace rather than the plugin installation directory.

WAV rendering has no artificial duration ceiling. Audio is synthesized in bounded chunks rather than held in memory as one giant sample array. The default resource budget is 512 MiB per rendered WAV—roughly 100 minutes at the built-in 22.05 kHz, 16-bit stereo format—and `--max-wav-mib` can raise it to 4095 MiB (about 13.5 hours at that format), the practical limit of classic RIFF/WAV. The command reports the estimated storage requirement before rendering if the selected budget is too small. `--max-events` defaults to 5,000 as an accidental-run guard but can be explicitly raised to 100,000. Diff mode separately uses a genuinely computational bounded display alignment and requires `(len(reference)+1) × (len(variant)+1) ≤ 1,000,000` cells.

`encode` and its compatibility alias `sonify` carry bounded versioned codec records in MIDI, MusicXML, and a custom RIFF/WAV chunk. The WAV contract binds the stored normalized sequence to the exact PCM bytes; MusicXML validates parsed pitches and visible residue/base lyric labels. “Exact” means the normalized sequence symbols and resolved molecule kind are identical, not that the original FASTA formatting was reconstructed. These hashes establish internal consistency, not authorship. The SVG and MusicXML are deliberately readable pitch/event grids; audio timing, stereo pan, and rich event semantics remain authoritative in the MIDI/CSV ledger.

For metadata-free input, `audio-decode` accepts bounded mono/stereo integer PCM WAV and quantizes spectral windows into the selected alphabet. It has no duration timer; `--max-input-mib` controls the transparent input-memory budget (512 MiB by default, adjustable to 4095 MiB), and analysis holds only one spectral window at a time. `score-decode` reads one-part plain MusicXML and maps its note onsets. Both default to protein, while `--kind dna` and `--kind rna` are explicit alternatives. Their `lossy-*` status and confidence values describe transcription fit under the selected mapping. Supported import formats are integer PCM WAV and plain, one-part MusicXML. When an available local converter such as FFmpeg is present, Codex can also handle common audio and video containers by extracting the selected audio stream to PCM WAV before invoking Seq2Music; video frames are ignored.

The summary music from `diff`, `msa`, and `structure` remains intentionally noninvertible, even though those modes now receive visible SVG/MusicXML notation.

## Data handling

Processing is local. Exact MIDI is recoverable from its ordered notes plus interpretation metadata; exact WAV and MusicXML also carry embedded sequence metadata. HTML, CSV, manifests, and decoded FASTA can contain source symbols or identifiers. Treat these artifacts as copies of the source data when deciding whether to share them. See the repository-level `PRIVACY.md` for details.

macOS, Linux, and Windows are supported. On platforms without secure directory-descriptor operations, create a new output path instead of using `--force`.

MSA mode accepts the extended DNA/RNA IUPAC ambiguity alphabet when `--kind dna` or `--kind rna` is supplied. Ambiguity codes remain visible in coverage fields and are excluded from consensus and entropy. Structure mode uses the first encountered model by default; `--model N` selects a specific model and records it in the CSV and run manifest.

## Repository installation

This plugin is distributed through the repository-level Codex marketplace in the public Seq2Music repository:

```bash
codex plugin marketplace add jvogan/seq2music --ref v0.3.2
codex plugin add seq2music@seq2music
```

Start a new Codex task after installation so the skill is loaded. For local development from a clone, run `codex plugin marketplace add .` from the repository root. The CLI can also be run directly from this plugin directory with Python 3.10 or newer.
