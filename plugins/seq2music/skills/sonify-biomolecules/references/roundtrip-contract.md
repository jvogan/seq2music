# Seq2Music round-trip contract

## Exact scope

`encode` and the compatibility command `sonify` create three self-describing exact carriers: MIDI under `seq2music-roundtrip-v1`, PCM WAV under `seq2music-wave-roundtrip-v1`, and MusicXML under `seq2music-score-roundtrip-v1`. Exact recovery means the same normalized protein, DNA, or RNA symbols and resolved kind passed the carrier-specific integrity checks.

Exact recovery covers normalized symbols and the resolved molecule kind. Decode emits a canonical 70-column FASTA; headers are retained as labels, while source comments, spacing, case, wrapping, and file encoding remain outside the codec.

## What is carried in MIDI

The musical notes use the public `seq2music-mapping-v2` pitch map. A tick-zero Standard MIDI File sequencer-specific meta event carries bounded canonical JSON declaring the codec profile, mapping algorithm, molecule kind, normalized sequence length and SHA-256, reversible channel, timing grid, and record label. Ordinary MIDI players ignore this record and play the notes normally.

Strict decode checks the codec profile, mapping, kind, PPQ, note count, monophonic timing schedule, pitch alphabet, sequence length, and sequence digest before writing output.

## Decode statuses

- `exact-embedded`: the normalized sequence matches the valid codec embedded in MIDI.
- `exact-manifest-matched`: the embedded contract matches and the MIDI bytes also match the supplied run manifest's SHA-256 record.
- `edited-or-unverified`: the user explicitly passed `--allow-edited` (or `--legacy` plus an explicit kind) to export a changed or metadata-free note stream.
- `exact-wav-embedded`: stored normalized symbols, alphabet, sequence digest, PCM digest, and WAV parameters agree.
- `exact-score-embedded`: stored normalized symbols and digest agree with the actual parsed monophonic MusicXML note stream.
- `lossy-audio-quantized`: metadata-free or explicitly edited PCM was quantized into the selected pitch alphabet.
- `lossy-score-transcription`: metadata-free or explicitly edited MusicXML was mapped under the reported note/chord policy.

A matching manifest establishes carrier consistency with that supplied manifest.

## Noninvertible outputs

Metadata-stripped WAV is acoustically quantized, and arbitrary MusicXML is transcribed under explicit policies. `diff`, `msa`, and `structure` summarize comparisons, alignment columns, or coarse coordinate-derived features; their music and notation omit source information and are noninvertible.
