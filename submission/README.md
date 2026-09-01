# OpenAI Plugins Directory submission packet

This packet contains the public, copy-ready fields for the initial Seq2Music submission. The submission type is **Skills only**. Seq2Music has no MCP server, hosted backend, account system, authentication flow, telemetry, or network dependency of its own.

Official submission portal: <https://platform.openai.com/plugins>

## Public listing

| Field | Value |
| --- | --- |
| Plugin name | Seq2Music |
| Publisher | Jacob Vogan |
| Short description | Turn biological sequences into music—and music back into sequences. |
| Category | Scientific Research |
| Website | <https://github.com/jvogan/seq2music> |
| Support URL | <https://github.com/jvogan/seq2music/issues> |
| Privacy URL | <https://github.com/jvogan/seq2music/blob/main/PRIVACY.md> |
| Terms URL | <https://github.com/jvogan/seq2music/blob/main/TERMS.md> |
| Logo | `plugins/seq2music/assets/seq2music-protein-note-logo.png` |
| Submission type | Skills only |

### Long description

Seq2Music turns protein, DNA, and RNA into audible music, visible notation, and reversible MIDI, WAV, and MusicXML. It can recover normalized sequences from exact Seq2Music carriers, transcribe ordinary PCM audio or MusicXML into a selected biological alphabet, and create coordinated audio and visual outputs for variants, alignments, and protein structures. Processing runs locally, with no account or network service required by Seq2Music.

## Starter prompts

1. Turn this protein FASTA into reversible music, show me the notation, and decode one exact carrier to verify the normalized sequence.
2. Take the audio track from this attached video and turn it into a protein sequence with a visible score and transcription report.
3. Make the differences between these two biological sequences audible and give me the synchronized visual comparison.
4. Turn this multiple-sequence alignment into coordinated music and an accessible alignment visualization.
5. Sonify the residue-order trace in this PDB or mmCIF protein structure and return the music, score, and event data.

## Positive test cases

### P1 — Direct reversible protein round trip

- **Prompt:** “Use Seq2Music to encode the attached protein FASTA as MIDI, WAV, and MusicXML. Show the score, decode the WAV, and verify that the normalized protein sequence matches.”
- **Expected activation:** Activate `sonify-biomolecules`; run `encode`, then `audio-decode` on the generated WAV.
- **Expected result shape:** Links to MIDI, WAV, MusicXML, SVG score, HTML, CSV ledger, and run manifest; decoded FASTA and decode report; exact carrier status; equality result and hashes.
- **Fixture:** [`demo-protein.fasta`](../plugins/seq2music/assets/examples/demo-protein.fasta)

### P2 — Indirect audible variant comparison

- **Prompt:** “Make this protein mutation something I can hear and see. Compare the reference and variant, and point me to the changed event.”
- **Expected activation:** Activate `sonify-biomolecules`; infer protein input and run `diff`.
- **Expected result shape:** Coordinated reference/variant audio, notation or HTML visualization, event ledger, alignment summary, and the changed residue position.
- **Fixtures:** [`demo-protein.fasta`](../plugins/seq2music/assets/examples/demo-protein.fasta) and [`demo-protein-variant.fasta`](../plugins/seq2music/assets/examples/demo-protein-variant.fasta)

### P3 — Multiple-sequence alignment sonification

- **Prompt:** “Turn this aligned FASTA into music so I can inspect conservation and variation across the records.”
- **Expected activation:** Activate `sonify-biomolecules`; inspect the input and run `msa`.
- **Expected result shape:** Alignment audio, accessible HTML/SVG visualization, per-column event data, summary, and run manifest.
- **Fixture:** [`demo-alignment.fasta`](../plugins/seq2music/assets/examples/demo-alignment.fasta)

### P4 — Protein structure trace

- **Prompt:** “Create music and a figure from the residue-order trace in this PDB structure.”
- **Expected activation:** Activate `sonify-biomolecules`; run `structure` with the observed model and chain selection recorded.
- **Expected result shape:** Structure-derived audio, visible notation or HTML, event CSV, summary, and manifest that reports model, chain, and observed residue coverage.
- **Fixture:** [`demo-structure.pdb`](../plugins/seq2music/assets/examples/demo-structure.pdb)

### P5 — Ordinary audio to protein sequence

- **Prompt:** “Turn this ordinary PCM WAV into a protein sequence. Give me FASTA, a visible score, and the exact transcription settings you used.”
- **Expected activation:** Activate `sonify-biomolecules`; run `audio-decode` with the protein alphabet.
- **Expected result shape:** Protein FASTA, `lossy-audio-quantized` status, decode report, inference event CSV, and SVG score. The response must describe the quantization settings without claiming an exact carrier round trip.
- **Fixture:** [`seq2music-plain-audio-test.wav`](https://github.com/jvogan/seq2music/releases/download/v0.3.2/seq2music-plain-audio-test.wav)

## Negative test cases

### N1 — Speech transcription

- **Prompt:** “Transcribe the spoken interview in this MP4 into a text document.”
- **Expected behavior:** Do not activate Seq2Music. Use an available speech-transcription workflow or explain what input is needed for one.
- **Why:** The requested output is human-language speech transcription, not mapping audio into a protein, DNA, or RNA alphabet.

### N2 — Ordinary media conversion

- **Prompt:** “Convert this WAV file to MP3 and keep the audio sounding the same.”
- **Expected behavior:** Do not activate Seq2Music. Use an ordinary media-conversion workflow if one is available.
- **Why:** The user wants container or codec conversion, not sequence sonification or audio-to-sequence transcription.

### N3 — Streaming music playback

- **Prompt:** “Play my favorite album from my music streaming account.”
- **Expected behavior:** Do not activate Seq2Music. Use the relevant music-service integration if one is available.
- **Why:** The request concerns catalog playback and account access, neither of which Seq2Music provides.

## Availability

Recommended selection: all countries and regions offered by the portal where the publisher is prepared to support the plugin. The publisher must confirm this selection in the portal before attesting.

## Initial release notes

Initial public submission of Seq2Music v0.3.2 as a skills-only plugin. Seq2Music provides local biological sequence-to-music workflows for protein, DNA, and RNA; exact normalized-sequence recovery from self-describing MIDI, PCM WAV, and MusicXML; explicitly lossy transcription from ordinary PCM WAV and MusicXML; accessible notation and browser visuals; and sonification workflows for variants, alignments, and protein structures. No account, authentication, hosted service, test credentials, or network access is required by the plugin.

## Reviewer setup

- Python 3.10 or newer.
- No credentials or private network access.
- The skill bundle is the `plugins/seq2music` directory from the `v0.3.2` tag or the matching release asset.
- Public fixtures are included in the plugin; the ordinary-audio fixture is attached to the GitHub release.
- The plugin writes only to the output directory selected for the run.

## Portal-only actions

These actions require the publisher and are intentionally not pre-completed:

- Select the verified **Jacob Vogan** developer identity in the publishing organization.
- Confirm the final country and region availability.
- Review and sign the policy attestations.
- Select **Submit for Review**.
- After approval, select **Publish**.

Do not add a personal email address, phone number, mailing address, local filesystem path, credential, or private account identifier to the public listing. Use the GitHub website and Issues URLs above for public contact and support.

## Optional showcase copy

**Title:** Seq2Music — reversible biological sequence sonification in Codex

**Summary:** Seq2Music turns protein, DNA, and RNA into audible music, visible notation, and inspectable event data. Exact Seq2Music carriers can recover normalized sequences, while ordinary audio and scores can be projected into selected biological alphabets. The plugin runs locally and pairs every sound with accessible visual and tabular outputs.

Recommended showcase media: the repository banner, one short exact-round-trip audio demo, the matching score, and a concise screen recording of Codex generating and reopening the artifacts.
