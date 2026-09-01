# Audio and score import

## Pick the carrier

- Use `audio-decode` for uncompressed integer PCM RIFF/WAV. Mono and stereo, 8/16/24/32-bit samples, and 8–192 kHz rates are supported without a duration timer. `--max-input-mib` is an explicit memory/file budget: 512 MiB by default and adjustable up to 4095 MiB, the classic-RIFF range.
- Use `score-decode` for bounded, plain MusicXML. It handles an ordered pitch/event grid with explicit part, voice, and chord policies.
- Omitted `--kind` selects protein for a carrier without a valid embedded Seq2Music kind. Use `--kind dna` or `--kind rna` for those alphabets; their shared pitch 72 becomes `T` or `U` according to the selection.

## Exact versus lossy

An intact Seq2Music WAV stores the normalized sequence in a bounded RIFF chunk bound to the PCM-data SHA-256 and format parameters. Exactness comes from this PCM-bound container record. Metadata-free recovery uses spectral quantization. Editors and format converters may remove the record; `--allow-edited` permits fallback to lossy quantization when PCM integrity changed.

An intact Seq2Music MusicXML score stores a bounded contract and validates its ordered monophonic pitches plus visible residue/base lyric labels. With `--allow-edited`, a supported edited score can be transcribed lossily. Backup/forward time shifts and multi-voice timelines are outside this ordered pitch/event-grid importer.

For arbitrary WAV, Seq2Music removes DC per analysis window, applies a Hann window, sums Goertzel power across channels for each allowed mapping pitch and low harmonics, then chooses the strongest pitch. Adjacent equal voiced windows collapse into one symbol unless `--no-collapse` is given. Only the current analysis window is converted to floating-point samples, so long inputs keep a bounded analysis footprint. `--window-seconds`, `--silence-threshold`, and the input byte budget are explicit.

For arbitrary MusicXML, one part and at most one voice are accepted by default. Chords use `--polyphony lowest`, `highest`, or `first`; `reject` is the default. Notes outside the mapping palette are assigned to the nearest allowed MIDI pitch and recorded in the receipt.

## Interpret the receipt

For lossy audio, `confidence` measures separation between the strongest and runner-up allowed spectral candidates. For lossy score import, it measures nearest-palette pitch distance. Both are deterministic fit heuristics reported with the mapping parameters. Exact container decodes report confidence as not applicable.

Return `sequence.fasta`, `decode.report.json`, `inference.events.csv`, and `score.svg`. For a lossy result, describe the selected alphabet, mapping, parameters, and transcription status.
