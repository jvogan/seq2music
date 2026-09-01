# Seq2Music mapping contract

## Shared rules

Mapping ID: `seq2music-mapping-v2`. Round-trip codec ID: `seq2music-roundtrip-v1`. Time advances by one fixed number of beats per residue or alignment column; BPM therefore changes both MIDI and WAV duration. No hidden randomness is used. MIDI, WAV, static HTML, CSV, text summary, and run manifest are generated together.

Every render also receives a deterministic grand-staff SVG and MusicXML pitch/event-grid score. These are accessible visual and portable machine-readable views of ordered events. MIDI and the CSV ledger retain the richer timing, stereo-pan, and overlapping-track semantics. Chords and silent events remain visible, and residue/base labels appear below note positions. Single-sequence `encode`/`sonify` scores carry an exact score contract that binds pitches and visible residue/base lyric labels.

Metadata-free audio import is a deterministic projection into this same pitch alphabet. It analyzes bounded PCM WAV windows with Hann-weighted Goertzel powers, sums stereo channel powers, chooses the strongest allowed pitch, and collapses consecutive identical voiced windows by default. Metadata-free MusicXML uses parsed note onsets and an explicit chord policy.

## Protein sequence

The 20 canonical residues map one-to-one to an intentionally designed C-major pitch palette in canonical one-letter order:

```text
A C D E F G H I K L M N P Q R S T V W Y
48 50 52 53 55 57 59 60 62 64 65 67 69 71 72 74 76 77 79 81
```

This identity palette is invertible when the carrier includes the Seq2Music kind and codec contract. Velocity encodes the Kyte-Doolittle hydropathy value on a MIDI contrast scale. A small deterministic waveform variation marks broad residue classes.

## DNA and RNA

`A=60`, `C=64`, `N=66`, `G=67`, `T/U=72`. DNA and RNA therefore require the embedded resolved kind to distinguish `T` from `U`; pitch alone cannot. G/C receive higher velocity than A/T/U to provide a separate GC-content cue. The IUPAC ambiguity symbol `N` represents an unknown base and receives a quiet, distinct cue. Purines and pyrimidines use small deterministic waveform differences.

## Reference/variant diff

Seq2Music performs a deterministic, linear-gap Needleman-Wunsch display alignment. Its exact safety bound is `(len(reference)+1) × (len(variant)+1) ≤ 1,000,000` dynamic-programming cells, and the rendered event count is `len(reference)+len(variant)`. Substantially diverged inputs benefit from a suitable domain alignment first. The reference is panned left and the variant right in WAV and by MIDI pan control. Identical symbols play together; substitutions, insertions, and deletions receive a velocity accent and are labeled in the event ledger.

## MSA choir

Input must contain at least two sequences and already be aligned. Consensus identity controls the root pitch. Unweighted Shannon entropy is calculated over observed canonical, non-gap symbols and normalized to its alphabet maximum. Extended nucleotide ambiguity codes `N,R,Y,S,W,K,M,B,D,H,V` are accepted only in MSA mode with an explicit nucleotide kind; they are excluded from consensus and entropy and reported through `ambiguous_coverage`. Entropy selects an ordinal roughness-cue interval palette: unison, perfect fifth, major third, minor third, tritone, then minor second. Non-gap coverage controls velocity. Zero observed entropy is labeled “no observed variation.” Results describe the supplied alignment.

## Residue-order structure geometry trace

Structure mode reads supported C-alpha records in the first encountered PDB/mmCIF model unless the user selects a positive model number with `--model`. The selected model is retained in every structure event and the run manifest. For alternate conformers, it selects highest occupancy with a deterministic blank/A/lexical tie break. Original residue names, occupancy, and alternate IDs remain in the ledger; SEC, PYL, and MSE analog mappings are explicitly marked.

- residue identity controls the root pitch;
- tied-midrank of C-alpha distance from the selected structure centroid controls octave;
- the count of nonlocal C-alpha neighbors within 8 angstroms controls velocity and whether fifth/octave chord tones are added. Local neighbors are approximated using author residue numbers within the same chain, so missing observed residues preserve source-coordinate spacing.

Centroid distance midranks and contact counts are invariant to rigid translation and rotation up to floating-point tolerance and source-coordinate precision. The inclusive 8-angstrom comparison uses a small squared-distance tolerance. Contact count is reported as the defined geometry statistic. Source B-factor values are preserved in the event ledger.

## Provenance

The run ID hashes the canonical event stream, title, input basenames, mode, semantic details (including any round-trip record), render parameters, sample rate, synthesizer ID, and mapping version. The manifest records source hashes and hashes each generated artifact. `verify` compares listed generated artifacts with the supplied manifest. The MIDI-only and manifest-matched decode guarantees are defined separately in [roundtrip-contract.md](roundtrip-contract.md).
