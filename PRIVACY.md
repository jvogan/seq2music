# Privacy

Seq2Music runs locally and does not include network requests, telemetry, accounts, or authentication. It reads only the files selected by the user and writes only to the requested output directory.

Exact Seq2Music carriers are deliberately recoverable. MIDI encodes the sequence in its ordered notes and carries the metadata needed to interpret them; it omits a duplicate plaintext sequence field. WAV and MusicXML carry embedded sequence metadata in addition to their audible or visible events. Treat every exact carrier as a copy of the source sequence when deciding whether to share it. HTML, CSV, manifests, and decoded FASTA can also disclose source symbols, record names, file basenames, or content hashes.

Lossy transcription of ordinary WAV or MusicXML uses the local media supplied by the user. Seq2Music does not upload that media.
