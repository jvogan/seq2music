---
name: sonify-biomolecules
description: Reversibly encode normalized protein, DNA, or RNA sequences as auditable MIDI, PCM WAV, and MusicXML with visible SVG notation; decode those carriers; or quantize arbitrary PCM WAV/MusicXML into an explicitly lossy sequence. Use for sequence-to-music, music/audio/score-to-sequence, visible notation, audible mutation comparison, alignment audio, or structure-to-music.
---

# Sonify Biomolecules with Seq2Music

Seq2Music turns local biological data into deterministic sound, visible notation, and a traceable event ledger.

## Choose a mode

- One FASTA record, default reversible path: `encode` (`sonify` is a compatibility alias and is also reversible).
- Seq2Music sequence MIDI back to canonical FASTA: `decode`.
- Seq2Music or arbitrary PCM WAV to FASTA: `audio-decode`.
- Video or unsupported audio container: use an available local converter such as FFmpeg to extract the requested audio stream as mono or stereo integer PCM WAV, then pass that WAV to `audio-decode`; ignore the video frames.
- Seq2Music or arbitrary plain MusicXML to FASTA: `score-decode`.
- Reference and variant FASTA: `diff`.
- Pre-aligned, equal-length FASTA records: `msa`.
- Local PDB or mmCIF protein structure: `structure`.
- Input validation only: `inspect`.
- Artifact integrity check: `verify`.

Before explaining a mapping, round-trip status, audio/score import, or the papers behind the plugin, read [references/mapping-contract.md](references/mapping-contract.md), [references/roundtrip-contract.md](references/roundtrip-contract.md), [references/audio-score-import.md](references/audio-score-import.md), or [references/evidence-boundaries.md](references/evidence-boundaries.md), respectively.

## Run locally

Locate the plugin root from this skill and invoke `scripts/seq2music.py` with Python 3. The script uses only the standard library and makes no network calls.

Examples:

```bash
python3 scripts/seq2music.py inspect --input protein.fasta --kind protein
run_dir="$(python3 scripts/seq2music.py encode --input protein.fasta --kind protein --out seq2music-output)"
python3 scripts/seq2music.py decode --midi "$run_dir/protein.mid" --manifest "$run_dir/protein.run.json" --out recovered
python3 scripts/seq2music.py audio-decode --audio music.wav --out recovered-from-audio
python3 scripts/seq2music.py audio-decode --audio music.wav --kind dna --out projected-dna
python3 scripts/seq2music.py score-decode --score music.musicxml --out recovered-from-score
python3 scripts/seq2music.py diff --reference ref.fasta --variant mutant.fasta --kind protein --out seq2music-output
python3 scripts/seq2music.py msa --input aligned.fasta --kind protein --out seq2music-output
python3 scripts/seq2music.py structure --input structure.pdb --model 1 --chains A --out seq2music-output
```

Use an explicit `--kind` when an all-ACGT string could be either a short protein or DNA, and for ambiguity-only nucleotide alignments whose IUPAC letters overlap protein symbols. `encode` requires exactly one FASTA record. `diff` reads the first record from each input; `msa` requires at least two records and uses all of them. MSA accepts extended DNA/RNA IUPAC ambiguity codes, but excludes them from consensus and entropy while reporting ambiguous coverage. Diff uses a bounded display alignment of at most 1,000,000 dynamic-programming cells. `--step` is measured in beats, so BPM changes both MIDI and WAV speed. There is no duration ceiling: rendering is streamed under a 512 MiB default WAV budget, adjustable with `--max-wav-mib` up to the 4095 MiB classic-RIFF limit. `--max-events` defaults to 5,000 but can be raised explicitly to 100,000. `audio-decode` likewise uses `--max-input-mib` rather than a duration timer and analyzes one window at a time. Structure mode reads observed C-alpha atoms from the first encountered model unless `--model N` is provided, records the selection, and does not rebuild missing residues. mmCIF parsing supports a bounded `_atom_site` subset and reports unsupported syntax.

## Present the result

Return links to the generated `.html`, `.score.svg`, `.musicxml`, `.wav`, `.mid`, `.events.csv`, `.summary.txt`, and `.run.json` files. For decoding, return `sequence.fasta`, `decode.report.json`, `inference.events.csv`, and `score.svg`. State the exact status verbatim: MIDI uses `exact-embedded`/`exact-manifest-matched`; WAV uses `exact-wav-embedded`; score uses `exact-score-embedded`; arbitrary imports use `lossy-audio-quantized` or `lossy-score-transcription`. Lead with the carrier, selected alphabet, exactness, and event count.

Describe what the selected mapping encodes, which carrier was decoded, and what the reported status establishes. Keep broader interpretation tied to the event ledger and the source analysis chosen by the user.

## Output behavior

- Outputs are local and offline.
- Exact MIDI is recoverable from ordered notes plus interpretation metadata; exact WAV and MusicXML also carry embedded sequence metadata. Tell the user before sharing any carrier, HTML, CSV, manifest, or decoded FASTA derived from sensitive input.
- WAV duration is not artificially capped. Render, import, and verification byte budgets are explicit and adjustable; they protect storage/memory and classic RIFF limits rather than imposing a preview-length policy.
- A run is deterministic for the same normalized input and parameters.
- New single-sequence MIDI, WAV, and MusicXML embed separate bounded contracts. WAV exactness is container-level and PCM-bound; MusicXML exactness validates ordered pitches and visible residue/base labels. SVG/MusicXML are pitch/event-grid views, while MIDI/CSV retain timing and pan semantics.
- Exact recovery covers normalized symbols, not byte-identical FASTA formatting. Metadata-stripped/arbitrary audio and scores are lossy transcriptions; MSA, diff, and structure music are noninvertible.
- Existing run directories are refused unless `--force` is explicitly requested.
- Each run manifest stores input and artifact SHA-256 hashes.
- The HTML and CSV are first-class accessible representations; never make audio the only deliverable.

For starter prompts and mode-specific examples, read [references/examples.md](references/examples.md).
