#!/usr/bin/env python3
"""Seq2Music: deterministic, standard-library biological data sonification."""

from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import html
import io
import json
import math
import mmap
import os
import re
import stat
import struct
import sys
import tempfile
import wave
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

VERSION = "0.3.1"
ALGORITHM_VERSION = "seq2music-mapping-v2"
ROUNDTRIP_VERSION = "seq2music-roundtrip-v1"
ROUNDTRIP_MAGIC = b"SEQ2MUSIC\x00"
WAVE_ROUNDTRIP_VERSION = "seq2music-wave-roundtrip-v1"
SCORE_ROUNDTRIP_VERSION = "seq2music-score-roundtrip-v1"
WAVE_CODEC_CHUNK = b"s2mc"
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
DNA_IUPAC = set("ACGTRYSWKMBDHVN")
RNA_IUPAC = set("ACGURYSWKMBDHVN")
AA_PITCH = dict(zip(AA_ORDER, [48, 50, 52, 53, 55, 57, 59, 60, 62, 64, 65, 67, 69, 71, 72, 74, 76, 77, 79, 81]))
HYDROPATHY = dict(zip(AA_ORDER, [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8, 1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3]))
AA_CLASS = {
    **dict.fromkeys("AVILMFWY", "hydrophobic"),
    **dict.fromkeys("STNQ", "polar"),
    **dict.fromkeys("KRH", "positive"),
    **dict.fromkeys("DE", "negative"),
    **dict.fromkeys("CGP", "special"),
}
NUC_PITCH = {"A": 60, "C": 64, "N": 66, "G": 67, "T": 72, "U": 72}
ENTROPY_INTERVALS = (0, 7, 4, 3, 6, 1)
MAX_INPUT_BYTES = 20 * 1024 * 1024
MAX_ALIGNMENT_CELLS = 1_000_000
MAX_EVENTS_HARD = 100_000
MAX_STRUCTURE_RESIDUES = 2_500
SAMPLE_RATE = 22_050
MAX_MIDI_BYTES = 20 * 1024 * 1024
MIB = 1024 * 1024
DEFAULT_MAX_WAV_MIB = 512
DEFAULT_MAX_INPUT_WAV_MIB = 512
DEFAULT_MAX_VERIFY_MIB = 512
MAX_RIFF_MIB = 4095
MAX_RIFF_BYTES = 0xFFFFFFFF + 8
MAX_SCORE_BYTES = 64 * 1024 * 1024
MAX_SCORE_ELEMENTS = 1_500_000
MAX_WAVE_CHUNKS = 1024
MAX_MIDI_TRACKS = 64
MAX_MIDI_EVENTS = 100_000
MAX_CODEC_BYTES = 1024 * 1024
MAX_MANIFEST_ARTIFACTS = 64
THREE_TO_ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "SEC":"C","PYL":"K","MSE":"M",
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def open_regular_fd(path: Path) -> tuple[int,os.stat_result]:
    preflight=os.lstat(path); reparse=getattr(stat,"FILE_ATTRIBUTE_REPARSE_POINT",0)
    if stat.S_ISLNK(preflight.st_mode) or (reparse and getattr(preflight,"st_file_attributes",0)&reparse):
        raise ValueError(f"input is a symlink, junction, or reparse point: {path}")
    flags=os.O_RDONLY|getattr(os,"O_CLOEXEC",0)|getattr(os,"O_NOFOLLOW",0)
    fd=os.open(path,flags)
    try:
        info=os.fstat(fd)
        if not stat.S_ISREG(info.st_mode): raise ValueError(f"input is not a regular file: {path}")
        return fd,info
    except Exception:
        os.close(fd); raise


def open_new_binary(path: Path):
    flags=os.O_RDWR|os.O_CREAT|os.O_EXCL|getattr(os,"O_CLOEXEC",0)|getattr(os,"O_NOFOLLOW",0)|getattr(os,"O_BINARY",0)
    return os.fdopen(os.open(path,flags,0o600),"w+b")


def open_new_text(path: Path):
    return os.fdopen(os.open(path,os.O_WRONLY|os.O_CREAT|os.O_EXCL|getattr(os,"O_CLOEXEC",0)|getattr(os,"O_NOFOLLOW",0),0o600),
                     "w",encoding="utf-8",newline="\n")


def write_new_bytes(path: Path, data: bytes) -> None:
    with open_new_binary(path) as output: output.write(data)


def write_new_text(path: Path, data: str) -> None:
    with open_new_text(path) as output: output.write(data)


def sha256_file(path: Path, max_bytes: int | None = None) -> str:
    h = hashlib.sha256()
    fd,info=open_regular_fd(path)
    try:
        if max_bytes is not None and info.st_size>max_bytes: raise ValueError(f"file exceeds hashing safety limit: {path}")
        total=0
        while True:
            chunk=os.read(fd,1<<20)
            if not chunk: break
            total+=len(chunk)
            if max_bytes is not None and total>max_bytes: raise ValueError(f"file grew beyond hashing safety limit: {path}")
            h.update(chunk)
    finally: os.close(fd)
    return h.hexdigest()


def slug(value: str, max_length: int = 96) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", value).strip("-.") or "seq2music"
    device=cleaned.split(".",1)[0].upper()
    if device in {"CON","PRN","AUX","NUL"} or re.fullmatch(r"(?:COM|LPT)[1-9]",device):
        cleaned=f"seq-{cleaned}"
    return cleaned[:max_length].rstrip("-.") or "seq2music"


def runtime_version() -> str:
    """Return the installed plugin's exact cache-busted version when available."""
    manifest = Path(__file__).resolve().parents[1] / ".codex-plugin" / "plugin.json"
    try:
        data = strict_json_loads(read_file_bytes(manifest, 64 * 1024, "plugin manifest"), "plugin manifest")
        value = data.get("version")
        if isinstance(value, str) and re.fullmatch(r"[0-9A-Za-z.+_-]{1,80}", value):
            return value
    except (OSError, ValueError):
        pass
    return VERSION


def read_file_bytes(path: Path, limit: int, label: str = "input") -> bytes:
    fd,info=open_regular_fd(path)
    try:
        if info.st_size>limit: raise ValueError(f"{label} exceeds the {limit//(1024*1024)} MiB safety limit: {path}")
        with os.fdopen(os.dup(fd),"rb") as source:
            data=source.read(limit+1)
        if len(data)>limit: raise ValueError(f"{label} grew beyond the {limit//(1024*1024)} MiB safety limit: {path}")
        return data
    finally: os.close(fd)


def read_bytes_limited(path: Path) -> bytes:
    return read_file_bytes(path,MAX_INPUT_BYTES)


def read_text_limited(path: Path) -> str:
    return read_bytes_limited(path).decode("utf-8",errors="replace")


def parse_fasta_text(text: str, path: Path) -> list[tuple[str,str]]:
    records: list[tuple[str, str]] = []
    name = None
    parts: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts).upper().replace(" ", "")))
            name, parts = (line[1:].strip() or f"record-{len(records)+1}")[:200], []
        else:
            if name is None:
                name = (path.stem or "record")[:200]
            parts.append(re.sub(r"\s+", "", line))
    if name is not None:
        records.append((name, "".join(parts).upper()))
    if not records or any(not seq for _, seq in records):
        raise ValueError("FASTA must contain at least one non-empty sequence")
    return records


def read_fasta(path: Path) -> list[tuple[str, str]]:
    return parse_fasta_text(read_text_limited(path),path)


def read_fasta_snapshot(path: Path) -> tuple[list[tuple[str,str]],str]:
    raw=read_bytes_limited(path)
    return parse_fasta_text(raw.decode("utf-8",errors="replace"),path),sha256_bytes(raw)


def infer_kind(seq: str, requested: str) -> str:
    if requested != "auto":
        return requested
    chars = set(seq)
    if chars <= set("ACGTN-"):
        return "dna"
    if chars <= set("ACGUN-"):
        return "rna"
    # Extended nucleotide ambiguity codes overlap canonical protein letters,
    # so auto mode never guesses their molecule kind; use explicit --kind.
    return "protein"


def validate_sequence(seq: str, kind: str, allow_gaps: bool = False, allow_iupac: bool = False) -> None:
    if kind == "protein":
        valid = set(AA_ORDER)
    elif allow_iupac:
        valid = set(DNA_IUPAC if kind == "dna" else RNA_IUPAC)
    else:
        valid = set("ACGTN" if kind == "dna" else "ACGUN")
    if allow_gaps:
        valid.add("-")
    invalid = sorted(set(seq) - valid)
    if invalid:
        raise ValueError(f"unsupported {kind} symbols: {''.join(invalid)}")


def velocity_for_residue(residue: str) -> int:
    return max(24, min(112, 40 + int(48 * (HYDROPATHY[residue] + 4.5) / 9.0)))


def sequence_events(seq: str, kind: str, step: float, label: str = "sequence") -> list[dict]:
    events = []
    for i, symbol in enumerate(seq):
        if kind == "protein":
            pitch = AA_PITCH[symbol]
            velocity = velocity_for_residue(symbol)
            group = AA_CLASS[symbol]
            hydropathy = HYDROPATHY[symbol]
        else:
            pitch = NUC_PITCH[symbol]
            velocity = 38 if symbol == "N" else (84 if symbol in "GC" else 64)
            group = "ambiguous" if symbol == "N" else ("purine" if symbol in "AG" else "pyrimidine")
            hydropathy = None
        events.append({
            "event": i + 1, "position": i + 1, "symbol": symbol, "label": label,
            "start": round(i * step, 6), "duration": round(step * 0.82, 6),
            "notes": [pitch], "velocity": velocity, "pan": 0.0,
            "feature": group, "hydropathy": hydropathy,
            "mapping": "identity-pitch; hydropathy-velocity" if kind == "protein" else "base-pitch; gc-velocity; N-ambiguous-cue",
        })
    return events


def needleman_wunsch(a: str, b: str) -> tuple[str, str]:
    cells=(len(a)+1)*(len(b)+1)
    if cells > MAX_ALIGNMENT_CELLS:
        raise ValueError(f"display alignment would require {cells:,} cells; safety limit is {MAX_ALIGNMENT_CELLS:,}")
    match, mismatch, gap = 2, -1, -2
    rows, cols = len(a) + 1, len(b) + 1
    score = [[0] * cols for _ in range(rows)]
    trace = [[""] * cols for _ in range(rows)]
    for i in range(1, rows): score[i][0], trace[i][0] = i * gap, "U"
    for j in range(1, cols): score[0][j], trace[0][j] = j * gap, "L"
    for i in range(1, rows):
        for j in range(1, cols):
            choices = [(score[i-1][j-1] + (match if a[i-1] == b[j-1] else mismatch), "D"),
                       (score[i-1][j] + gap, "U"), (score[i][j-1] + gap, "L")]
            score[i][j], trace[i][j] = max(choices, key=lambda x: (x[0], {"D":2,"U":1,"L":0}[x[1]]))
    x, y, i, j = [], [], len(a), len(b)
    while i or j:
        t = trace[i][j]
        if t == "D": x.append(a[i-1]); y.append(b[j-1]); i -= 1; j -= 1
        elif t == "U": x.append(a[i-1]); y.append("-"); i -= 1
        else: x.append("-"); y.append(b[j-1]); j -= 1
    return "".join(reversed(x)), "".join(reversed(y))


def diff_events(a: str, b: str, kind: str, step: float) -> tuple[list[dict], str, str]:
    aa, bb = needleman_wunsch(a, b)
    events: list[dict] = []
    for i, (x, y) in enumerate(zip(aa, bb)):
        status = "match" if x == y else ("insertion" if x == "-" else "deletion" if y == "-" else "substitution")
        for symbol, label, pan in ((x, "reference", -0.82), (y, "variant", 0.82)):
            if symbol == "-":
                continue
            pitch = AA_PITCH[symbol] if kind == "protein" else NUC_PITCH[symbol]
            velocity = (velocity_for_residue(symbol) if kind == "protein" else (38 if symbol == "N" else (84 if symbol in "GC" else 64)))
            if status != "match": velocity = min(120, velocity + 20)
            events.append({"event": len(events)+1, "position": i+1, "symbol": symbol, "label": label,
                           "start": round(i*step, 6), "duration": round(step*0.82, 6), "notes":[pitch],
                           "velocity": velocity, "pan": pan, "feature": status,
                           "mapping": "aligned stereo diff; non-match accent"})
    return events, aa, bb


def msa_events(records: list[tuple[str, str]], kind: str, step: float) -> list[dict]:
    if kind not in ("protein","dna","rna"): raise ValueError("MSA kind must be protein, dna, or rna")
    for _,sequence in records: validate_sequence(sequence,kind,allow_gaps=True,allow_iupac=True)
    lengths = {len(s) for _, s in records}
    if len(lengths) != 1:
        raise ValueError("MSA mode requires pre-aligned sequences of equal length")
    max_entropy = math.log2(20 if kind == "protein" else 4)
    events = []
    canonical = set(AA_ORDER if kind == "protein" else ("ACGT" if kind == "dna" else "ACGU"))
    for i, column in enumerate(zip(*(s for _, s in records))):
        non_gap = [x for x in column if x != "-"]
        observed = [x for x in non_gap if x in canonical]
        ambiguous = [x for x in non_gap if x not in canonical]
        counts = Counter(observed)
        total = len(observed)
        entropy = -sum((n/total)*math.log2(n/total) for n in counts.values()) if total else 0.0
        consensus = sorted(counts, key=lambda x: (-counts[x], x))[0] if counts else "-"
        if not observed:
            ambiguous_symbols="".join(sorted(set(ambiguous)))
            symbol = (ambiguous[0] if len(set(ambiguous)) == 1 else "?") if ambiguous else "-"
            events.append({"event":i+1,"position":i+1,"symbol":symbol,"label":"alignment",
                           "start":round(i*step,6),"duration":round(step*0.82,6),"notes":[],
                           "velocity":0,"pan":0.0,"feature":"no-informative-observation","entropy_bits":0.0,
                           "entropy_ratio":0.0,"entropy_interval_semitones":0,
                           "coverage":round(len(non_gap)/len(records),6),"informative_coverage":0.0,
                           "ambiguous_coverage":round(len(ambiguous)/len(records),6),
                           "ambiguous_count":len(ambiguous),
                           "ambiguous_symbols":ambiguous_symbols,
                           "mapping":"gap/ambiguous-only column; silent ledger event"})
            continue
        base = (AA_PITCH[consensus] if kind == "protein" else NUC_PITCH[consensus]) if consensus != "-" else 48
        entropy_ratio = entropy / max_entropy if max_entropy else 0.0
        palette_index = min(len(ENTROPY_INTERVALS)-1, int(math.ceil(entropy_ratio * (len(ENTROPY_INTERVALS)-1)))) if entropy > 0 else 0
        interval = ENTROPY_INTERVALS[palette_index]
        notes = [base] if interval == 0 else [base, min(108, base + interval)]
        coverage = len(non_gap) / len(records)
        informative_coverage = total / len(records)
        feature = "no-observed-variation" if entropy == 0 and informative_coverage == 1 else ("low-information" if entropy == 0 else "variable")
        events.append({"event": i+1, "position": i+1, "symbol": consensus, "label":"alignment",
                       "start": round(i*step,6), "duration": round(step*0.82,6), "notes":notes,
                       "velocity": int(36+64*coverage), "pan":0.0, "feature":feature,
                       "entropy_bits": round(entropy,6), "entropy_ratio": round(entropy_ratio,6),
                       "entropy_interval_semitones": interval, "coverage": round(coverage,6),
                       "informative_coverage": round(informative_coverage,6),
                       "ambiguous_coverage": round(len(ambiguous)/len(records),6),
                       "ambiguous_count":len(ambiguous),"ambiguous_symbols":"".join(sorted(set(ambiguous))),
                       "mapping":"consensus-pitch; entropy-roughness-palette; coverage-velocity"})
    return events


def parse_pdb(path: Path, chain_filter: set[str] | None = None, source_text: str | None = None,
              model_selection: int | None = None) -> list[dict]:
    residues: dict[tuple, dict] = {}
    lines=(source_text if source_text is not None else read_text_limited(path)).splitlines()
    saw_model_record=any(line.startswith("MODEL") for line in lines)
    model = None if saw_model_record else 1; selected_model = model_selection
    seen_models: set[int] = set() if saw_model_record else {1}
    for line in lines:
        if line.startswith("MODEL"):
            try: model = int(line[10:14].strip())
            except ValueError: raise ValueError("PDB contains an invalid MODEL serial")
            if model < 1: raise ValueError("PDB MODEL serial must be a positive integer")
            seen_models.add(model)
            if selected_model is None: selected_model = model
            continue
        if line.startswith("ENDMDL"):
            model = None
            continue
        if saw_model_record and model is None and line.startswith(("ATOM  ","HETATM")):
            raise ValueError("PDB coordinate record appears outside a MODEL/ENDMDL block")
        if selected_model is not None and model != selected_model: continue
        if not line.startswith(("ATOM  ", "HETATM")) or line[12:16].strip() != "CA":
            continue
        alt = line[16:17]
        chain = line[21:22].strip() or "_"
        if chain_filter and chain not in chain_filter: continue
        name3 = line[17:20].strip().upper()
        if name3 not in THREE_TO_ONE: continue
        try:
            resnum = int(line[22:26]); ins = line[26:27].strip()
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            occupancy = float(line[54:60].strip() or 0)
            bfactor = float(line[60:66].strip() or 0)
        except ValueError:
            continue
        if not all(math.isfinite(v) for v in (x,y,z,occupancy,bfactor)): continue
        key = (chain, resnum, ins)
        current = residues.get(key)
        conformer_score=(occupancy, alt == " ", alt == "A", "".join(chr(255-ord(c)) for c in alt))
        current_score=None if current is None else (current["occupancy"], current["alt"] == ".", current["alt"] == "A", "".join(chr(255-ord(c)) for c in current["alt"]))
        if current is None or conformer_score > current_score:
            residues[key] = {"chain":chain,"resnum":resnum,"ins":ins,"symbol":THREE_TO_ONE[name3],"model":model,
                             "source_residue":name3,"analog_mapping":THREE_TO_ONE[name3] if name3 in ("SEC","PYL","MSE") else "",
                             "alt":alt.strip() or ".","x":x,"y":y,"z":z,"occupancy":occupancy,"bfactor":bfactor}
    result = sorted(residues.values(), key=lambda r:(r["chain"],r["resnum"],r["ins"]))
    if model_selection is not None and model_selection not in seen_models:
        preview=sorted(seen_models)[:16]
        shown=", ".join(map(str,preview))+(f", ... ({len(seen_models)} total)" if len(seen_models)>16 else "")
        suffix="; file has no MODEL records (only implicit model 1)" if not saw_model_record else f"; available models: {shown}"
        raise ValueError(f"requested PDB model {model_selection} is not present{suffix}")
    if not result:
        requested = f" in requested model {model_selection}" if model_selection is not None else ""
        if model_selection not in (None, 1) and not saw_model_record:
            requested += " (the file has no MODEL records and is treated as model 1)"
        raise ValueError(f"no supported protein C-alpha atoms found in PDB{requested}")
    return result


def cif_line_tokens(raw: str) -> list[str]:
    """Tokenize one CIF line; quotes are special only at a token boundary."""
    tokens=[]; i=0
    while i < len(raw):
        while i < len(raw) and raw[i].isspace(): i+=1
        if i >= len(raw): break
        if raw[i] in ("'", '"'):
            quote=raw[i]; i+=1; start=i
            while i < len(raw) and raw[i] != quote: i+=1
            if i >= len(raw): raise ValueError("unsupported unterminated quoted mmCIF token")
            tokens.append(raw[start:i]); i+=1
        else:
            start=i
            while i < len(raw) and not raw[i].isspace(): i+=1
            tokens.append(raw[start:i])
    return tokens


def parse_mmcif(path: Path, chain_filter: set[str] | None = None, source_text: str | None = None,
                model_selection: int | None = None) -> list[dict]:
    """Parse a bounded _atom_site loop subset without external dependencies."""
    lines=(source_text if source_text is not None else read_text_limited(path)).splitlines(); residues={}; i=0
    seen_models: set[int] = set(); selected_model=model_selection; model_column_seen=False
    while i < len(lines):
        if lines[i].strip() != "loop_": i+=1; continue
        i+=1; tags=[]
        while i < len(lines) and lines[i].strip().startswith("_"):
            tags.append(lines[i].strip().lower()); i+=1
        if not tags or not all(t.startswith("_atom_site.") for t in tags):
            while i < len(lines) and lines[i].strip() not in ("loop_", "#"): i+=1
            continue
        index={t.split(".",1)[1]:n for n,t in enumerate(tags)}
        def val(row,*names,default="?"):
            for name in names:
                name=name.lower()
                if name in index and index[name] < len(row) and row[index[name]] not in (".","?"): return row[index[name]]
            return default
        body=[]
        while i < len(lines):
            raw=lines[i].strip()
            if not raw: i+=1; continue
            if raw=="#" or raw=="loop_" or raw.startswith("_"): break
            if raw.startswith(";"): raise ValueError("semicolon-delimited mmCIF values are unsupported in _atom_site")
            body.extend(cif_line_tokens(raw)); i+=1
        if len(body) % len(tags): raise ValueError("malformed or unsupported mmCIF _atom_site loop row")
        model_column_seen = model_column_seen or "pdbx_pdb_model_num" in index
        for offset in range(0,len(body),len(tags)):
            row=body[offset:offset+len(tags)]
            if "pdbx_pdb_model_num" in index:
                model_token=row[index["pdbx_pdb_model_num"]]
                if model_token in (".","?"):
                    raise ValueError("mmCIF model number is missing in a declared model-number column")
            else:
                model_token="1"
            try: model_number=float(model_token)
            except ValueError: raise ValueError(f"mmCIF model number must be a positive integer, got {model_token!r}")
            if not math.isfinite(model_number) or not model_number.is_integer() or model_number < 1:
                raise ValueError(f"mmCIF model number must be a positive integer, got {model_token!r}")
            model=int(model_number)
            seen_models.add(model)
            if selected_model is None: selected_model=model
            if model != selected_model: continue
            atom=val(row,"auth_atom_id","label_atom_id")
            if atom != "CA": continue
            chain=val(row,"auth_asym_id","label_asym_id",default="_"); chain="_" if chain in (".","?") else chain
            if chain_filter and chain not in chain_filter: continue
            name3=val(row,"auth_comp_id","label_comp_id").upper()
            if name3 not in THREE_TO_ONE: continue
            try:
                resnum=int(float(val(row,"auth_seq_id","label_seq_id")))
                x=float(val(row,"Cartn_x")); y=float(val(row,"Cartn_y")); z=float(val(row,"Cartn_z"))
                occupancy=float(val(row,"occupancy",default="0")); bfactor=float(val(row,"B_iso_or_equiv",default="0"))
            except ValueError: continue
            if not all(math.isfinite(v) for v in (x,y,z,occupancy,bfactor)): continue
            ins=val(row,"pdbx_PDB_ins_code",default=""); ins="" if ins in (".","?") else ins
            alt=val(row,"label_alt_id",default=".")
            key=(chain,resnum,ins); current=residues.get(key)
            if current is None or (occupancy, alt in (".","?"), alt == "A", alt) > (current["occupancy"], current["alt"] in (".","?"), current["alt"] == "A", current["alt"]):
                residues[key]={"chain":chain,"resnum":resnum,"ins":ins,"symbol":THREE_TO_ONE[name3],"model":model,
                               "source_residue":name3,"analog_mapping":THREE_TO_ONE[name3] if name3 in ("SEC","PYL","MSE") else "",
                               "alt":alt,"x":x,"y":y,"z":z,"occupancy":occupancy,"bfactor":bfactor}
    result=sorted(residues.values(),key=lambda r:(r["chain"],r["resnum"],r["ins"]))
    if model_selection is not None and model_selection not in seen_models:
        preview=sorted(seen_models)[:16]
        shown=", ".join(map(str,preview))+(f", ... ({len(seen_models)} total)" if len(seen_models)>16 else "")
        suffix="; file has no model-number column (only implicit model 1)" if not model_column_seen else f"; available models: {shown}"
        raise ValueError(f"requested mmCIF model {model_selection} is not present{suffix}")
    if not result:
        requested = f" in requested model {model_selection}" if model_selection is not None else ""
        raise ValueError(f"no supported protein C-alpha atoms found in mmCIF _atom_site loop{requested}")
    return result


def structure_events(residues: list[dict], step: float) -> list[dict]:
    coords = [(r["x"],r["y"],r["z"]) for r in residues]
    nres = len(residues); sx=sum(x for x,_,_ in coords); sy=sum(y for _,y,_ in coords); sz=sum(z for _,_,z in coords)
    radial = [(nres*x-sx)**2 + (nres*y-sy)**2 + (nres*z-sz)**2 for x,y,z in coords]
    order = sorted(range(len(residues)), key=lambda i:radial[i])
    ranks = [0.0]*len(residues); cursor=0
    while cursor < len(order):
        end=cursor+1
        while end < len(order) and math.isclose(radial[order[cursor]],radial[order[end]],rel_tol=1e-12,abs_tol=1e-9): end+=1
        midrank=(cursor+end-1)/2
        for idx in order[cursor:end]: ranks[idx]=midrank
        cursor=end
    contacts = [0]*len(residues)
    for i in range(len(residues)):
        for j in range(i+1, len(residues)):
            # Approximate local sequence neighbors with author residue numbers,
            # not observed-list indices; the latter hides contacts across gaps.
            same_near = residues[i]["chain"] == residues[j]["chain"] and abs(residues[i]["resnum"]-residues[j]["resnum"]) <= 2
            dx=coords[i][0]-coords[j][0]; dy=coords[i][1]-coords[j][1]; dz=coords[i][2]-coords[j][2]
            if not same_near and dx*dx+dy*dy+dz*dz <= 64.0 + 1e-9:
                contacts[i] += 1; contacts[j] += 1
    events = []
    n = max(1, len(residues)-1)
    for i, r in enumerate(residues):
        q = ranks[i]/n
        shift = -12 if q < .25 else (12 if q >= .75 else 0)
        root = max(24,min(96,AA_PITCH[r["symbol"]]+shift))
        notes = [root]
        if contacts[i] >= 2: notes.append(min(108,root+7))
        if contacts[i] >= 5: notes.append(min(108,root+12))
        events.append({"event":i+1,"position":i+1,"source_position":f"{r['chain']}:{r['resnum']}{r['ins']}",
                       "symbol":r["symbol"],"model":r.get("model",1),"label":f"chain-{r['chain']}","start":round(i*step,6),
                       "duration":round(step*0.82,6),"notes":notes,"velocity":min(118,44+contacts[i]*7),
                       "pan":0.0,"feature":"contact-chord","radial_rank":ranks[i],"contact_count":contacts[i],
                       "bfactor":r["bfactor"],"occupancy":r.get("occupancy"),"alt":r.get("alt"),
                       "source_residue":r.get("source_residue",r["symbol"]),"analog_mapping":r.get("analog_mapping",""),
                       "mapping":"identity-pitch; radial-midrank-octave; 8A-contact-chord/velocity"})
    return events


def varlen(value: int) -> bytes:
    buf = [value & 0x7F]
    value >>= 7
    while value:
        buf.append((value & 0x7F) | 0x80); value >>= 7
    return bytes(reversed(buf))


def write_midi(path: Path, events: list[dict], bpm: int, metadata: dict | None = None) -> None:
    ppq = 480
    timeline = [
        (0, 0, bytes([0xFF,0x51,0x03]) + int(60_000_000/bpm).to_bytes(3,"big")),
        (0, 0, bytes([0xB0, 0x0A, 64])),
        (0, 0, bytes([0xB1, 0x0A, 12])),
        (0, 0, bytes([0xB2, 0x0A, 115])),
    ]
    if metadata:
        payload=ROUNDTRIP_MAGIC+json.dumps(metadata,sort_keys=True,separators=(",",":"),ensure_ascii=True).encode("utf-8")
        if len(payload)>MAX_CODEC_BYTES: raise ValueError("Seq2Music roundtrip metadata exceeds the codec safety limit")
        # 0x7F is the Standard MIDI File sequencer-specific meta event. It is
        # ignored by playback while keeping the musical note stream portable.
        timeline.append((0,0,bytes([0xFF,0x7F])+varlen(len(payload))+payload))
    codec_ticks=bool(metadata and metadata.get("profile")==ROUNDTRIP_VERSION and isinstance(metadata.get("step_ticks"),int) and isinstance(metadata.get("gate_ticks"),int))
    for event_index,event in enumerate(events):
        channel = 2 if event.get("pan",0) > .3 else (1 if event.get("pan",0) < -.3 else 0)
        start = event_index*metadata["step_ticks"] if codec_ticks else int(round(event["start"] * bpm * ppq / 60))
        end = start+metadata["gate_ticks"] if codec_ticks else int(round((event["start"]+event["duration"]) * bpm * ppq / 60))
        for note in event["notes"]:
            timeline.append((start, 1, bytes([0x90|channel,note,event["velocity"]])))
            timeline.append((end, 0, bytes([0x80|channel,note,0])))
    final_tick=(len(events)-1)*metadata["step_ticks"]+metadata["gate_ticks"] if codec_ticks and events else max((int(round((e["start"]+e["duration"]) * bpm * ppq / 60)) for e in events),default=0)
    timeline.append((final_tick,2,b"\xff\x2f\x00")); timeline.sort(key=lambda x:(x[0],x[1]))
    track = bytearray(); last = 0
    for tick, _, msg in timeline:
        track.extend(varlen(tick-last)); track.extend(msg); last=tick
    write_new_bytes(path,b"MThd"+struct.pack(">IHHH",6,0,1,ppq)+b"MTrk"+struct.pack(">I",len(track))+bytes(track))


def read_varlen(data: bytes, pos: int, end: int) -> tuple[int,int]:
    value=0
    for _ in range(4):
        if pos >= end: raise ValueError("truncated MIDI variable-length value")
        byte=data[pos]; pos+=1; value=(value<<7)|(byte&0x7F)
        if not byte&0x80: return value,pos
    raise ValueError("MIDI variable-length value exceeds four bytes")


def strict_json_loads(payload: bytes, label: str = "Seq2Music metadata") -> dict:
    def reject_duplicates(pairs):
        result={}
        for key,value in pairs:
            if key in result: raise ValueError(f"duplicate JSON key in {label}: {key}")
            result[key]=value
        return result
    def reject_constant(value): raise ValueError(f"invalid JSON constant in {label}: {value}")
    try:
        value=json.loads(payload.decode("utf-8"),object_pairs_hook=reject_duplicates,parse_constant=reject_constant)
    except (UnicodeDecodeError,json.JSONDecodeError,RecursionError) as exc:
        raise ValueError(f"invalid {label}: {exc}")
    if not isinstance(value,dict): raise ValueError(f"{label} must be an object")
    return value


def parse_midi(path: Path) -> dict:
    data=read_file_bytes(path,MAX_MIDI_BYTES,"MIDI")
    if len(data)<14 or data[:4]!=b"MThd": raise ValueError("not a Standard MIDI File")
    header_len=struct.unpack(">I",data[4:8])[0]
    if header_len<6 or 8+header_len>len(data): raise ValueError("invalid MIDI header length")
    fmt,ntracks,division=struct.unpack(">HHH",data[8:14])
    if fmt not in (0,1): raise ValueError("only MIDI format 0 or 1 is supported")
    if not 1<=ntracks<=MAX_MIDI_TRACKS: raise ValueError(f"MIDI track count must be 1..{MAX_MIDI_TRACKS}")
    if division&0x8000 or division==0: raise ValueError("SMPTE-time MIDI is not supported")
    pos=8+header_len; notes=[]; embedded=[]; serial=0; event_count=0
    for track_index in range(ntracks):
        if pos+8>len(data) or data[pos:pos+4]!=b"MTrk": raise ValueError("missing MIDI track chunk")
        size=struct.unpack(">I",data[pos+4:pos+8])[0]; pos+=8; end=pos+size
        if end>len(data): raise ValueError("truncated MIDI track chunk")
        tick=0; running=None; saw_eot=False
        while pos<end:
            event_count+=1
            if event_count>MAX_MIDI_EVENTS: raise ValueError(f"MIDI exceeds the {MAX_MIDI_EVENTS} event safety limit")
            delta,pos=read_varlen(data,pos,end); tick+=delta
            if pos>=end: raise ValueError("truncated MIDI event")
            lead=data[pos]
            if lead<0x80:
                if running is None: raise ValueError("MIDI running status without prior channel status")
                status=running
            else:
                status=lead; pos+=1
                if 0x80<=status<=0xEF: running=status
                else: running=None
            if status==0xFF:
                if pos>=end: raise ValueError("truncated MIDI meta event")
                meta_type=data[pos]; pos+=1; length,pos=read_varlen(data,pos,end)
                if pos+length>end: raise ValueError("truncated MIDI meta payload")
                payload=data[pos:pos+length]; pos+=length
                if len(payload)>MAX_CODEC_BYTES and payload.startswith(ROUNDTRIP_MAGIC):
                    raise ValueError("Seq2Music roundtrip metadata exceeds the codec safety limit")
                if meta_type != 0x7F and payload.startswith(ROUNDTRIP_MAGIC):
                    raise ValueError("Seq2Music roundtrip metadata must use MIDI sequencer-specific meta type 0x7F")
                if meta_type == 0x7F and payload.startswith(ROUNDTRIP_MAGIC):
                    value=strict_json_loads(payload[len(ROUNDTRIP_MAGIC):], "Seq2Music MIDI metadata")
                    embedded.append(value)
                if meta_type==0x2F:
                    if length!=0 or pos!=end: raise ValueError("MIDI end-of-track event must be empty and final")
                    saw_eot=True
                continue
            if status in (0xF0,0xF7):
                length,pos=read_varlen(data,pos,end)
                if pos+length>end: raise ValueError("truncated MIDI system-exclusive payload")
                pos+=length; continue
            if not 0x80<=status<=0xEF: raise ValueError(f"unsupported MIDI status 0x{status:02x}")
            high=status&0xF0; needed=1 if high in (0xC0,0xD0) else 2
            if pos+needed>end: raise ValueError("truncated MIDI channel event")
            first=data[pos]; second=data[pos+1] if needed==2 else 0; pos+=needed
            if first&0x80 or (needed==2 and second&0x80): raise ValueError("MIDI channel data bytes must be below 0x80")
            if high==0x90 and second>0:
                notes.append({"tick":tick,"track":track_index,"order":serial,"channel":status&0x0F,"pitch":first,"velocity":second}); serial+=1
        if not saw_eot: raise ValueError("MIDI track is missing its end-of-track event")
        pos=end
    if pos!=len(data): raise ValueError("unexpected bytes after declared MIDI tracks")
    if len(embedded)>1: raise ValueError("MIDI must contain exactly one Seq2Music metadata event")
    return {"format":fmt,"tracks":ntracks,"ppq":division,"notes":sorted(notes,key=lambda n:(n["tick"],n["track"],n["order"])),"metadata":embedded[0] if embedded else None,"sha256":sha256_bytes(data)}


def inverse_pitch_map(kind: str) -> dict[int,str]:
    if kind=="protein": return {pitch:symbol for symbol,pitch in AA_PITCH.items()}
    terminal="T" if kind=="dna" else "U"
    return {NUC_PITCH["A"]:"A",NUC_PITCH["C"]:"C",NUC_PITCH["N"]:"N",NUC_PITCH["G"]:"G",NUC_PITCH[terminal]:terminal}


def require_roundtrip_metadata(meta: dict) -> None:
    required={"profile":str,"class":str,"mode":str,"kind":str,"sequence_length":int,
              "normalized_sequence_sha256":str,"midi_channel":int,"algorithm":str,
              "ppq":int,"step_beats":(int,float),"step_ticks":int,"gate_ticks":int,"event_count":int}
    for key,expected in required.items():
        if key not in meta or not isinstance(meta[key],expected) or isinstance(meta[key],bool):
            raise ValueError(f"roundtrip metadata field {key!r} is missing or invalid")
    if not re.fullmatch(r"[0-9a-f]{64}",meta["normalized_sequence_sha256"]):
        raise ValueError("roundtrip metadata contains an invalid sequence SHA-256")
    if (meta["sequence_length"]<1 or meta["event_count"]<1 or meta["sequence_length"]!=meta["event_count"]
            or meta["event_count"]>MAX_EVENTS_HARD or not 0<=meta["midi_channel"]<=15
            or not 1<=meta["ppq"]<=0x7FFF or not math.isfinite(float(meta["step_beats"]))
            or not 0<float(meta["step_beats"])<=4.0 or not 1<=meta["step_ticks"]<=meta["ppq"]*4
            or not 1<=meta["gate_ticks"]<=meta["step_ticks"]):
        raise ValueError("roundtrip metadata contains invalid counts or MIDI channel")
    name=meta.get("record_name")
    if name is not None and (not isinstance(name,str) or len(name)>200 or len(name.encode("utf-8"))>800):
        raise ValueError("roundtrip metadata record_name must be a bounded string")


def decode_midi_sequence(midi_path: Path, manifest_path: Path | None = None, requested_kind: str | None = None, allow_edited: bool = False) -> dict:
    parsed=parse_midi(midi_path); embedded=parsed["metadata"]; manifest=None; manifest_meta=None; manifest_midi_match=None
    if manifest_path:
        if manifest_path.is_symlink(): raise ValueError("manifest must be a regular, non-symlink file")
        manifest=strict_json_loads(read_file_bytes(manifest_path,MAX_INPUT_BYTES,"manifest"),"manifest JSON"); manifest_meta=manifest.get("roundtrip")
        if manifest_meta is None: raise ValueError("manifest does not contain a Seq2Music roundtrip record")
        if manifest_meta is not None and not isinstance(manifest_meta,dict): raise ValueError("manifest roundtrip metadata must be an object")
        artifacts=manifest.get("artifacts")
        if not isinstance(artifacts,list) or not 1<=len(artifacts)<=MAX_MANIFEST_ARTIFACTS:
            raise ValueError(f"manifest artifacts must be a list of 1..{MAX_MANIFEST_ARTIFACTS} records")
        for artifact in artifacts:
            if (not isinstance(artifact,dict) or not isinstance(artifact.get("path"),str)
                    or not isinstance(artifact.get("sha256"),str)
                    or not re.fullmatch(r"[0-9a-f]{64}",artifact["sha256"])):
                raise ValueError("manifest contains an invalid artifact record")
        matches=[a for a in artifacts if isinstance(a,dict) and a.get("path")==midi_path.name]
        if len(matches)!=1: raise ValueError("manifest must contain exactly one artifact record for this MIDI filename")
        manifest_midi_match=parsed["sha256"]==matches[0]["sha256"]
        if manifest_midi_match is False and not allow_edited: raise ValueError("MIDI hash differs from manifest; pass --allow-edited to decode intentional note edits")
    if embedded and manifest_meta and embedded!=manifest_meta: raise ValueError("embedded MIDI and manifest roundtrip metadata conflict")
    meta=embedded or manifest_meta
    if not meta and not (requested_kind and allow_edited):
        raise ValueError("MIDI lacks Seq2Music roundtrip metadata; legacy recovery requires --kind and --allow-edited")
    if meta: require_roundtrip_metadata(meta)
    if meta and meta.get("profile")!=ROUNDTRIP_VERSION: raise ValueError(f"unsupported roundtrip profile: {meta.get('profile')}")
    if meta and meta.get("mode") not in ("encode","sonify"): raise ValueError(f"mode {meta.get('mode')} is not exactly reversible")
    if meta and meta.get("algorithm")!=ALGORITHM_VERSION: raise ValueError(f"unsupported mapping algorithm: {meta.get('algorithm')}")
    if meta and meta.get("class")!="exact-normalized-sequence": raise ValueError("MIDI does not declare an exact normalized-sequence codec")
    layout_match=True
    if meta and parsed["ppq"]!=meta["ppq"]:
        layout_match=False
        if not allow_edited: raise ValueError("MIDI PPQ differs from roundtrip metadata")
    embedded_kind=meta.get("kind") if meta else None
    if requested_kind and embedded_kind and requested_kind!=embedded_kind: raise ValueError("requested kind conflicts with roundtrip metadata")
    kind=requested_kind or embedded_kind
    if kind not in ("protein","dna","rna"): raise ValueError("decode requires protein, dna, or rna kind")
    channel=int(meta.get("midi_channel",0)) if meta else 0
    selected=[n for n in parsed["notes"] if n["channel"]==channel]
    if not selected: raise ValueError(f"no note-on events found on reversible channel {channel}")
    ticks=Counter(n["tick"] for n in selected)
    collisions=sorted(tick for tick,count in ticks.items() if count!=1)
    if collisions: raise ValueError(f"reversible channel is polyphonic at {len(collisions)} tick positions")
    if meta:
        if len(selected)!=meta["event_count"]:
            layout_match=False
            if not allow_edited: raise ValueError("reversible note count differs from roundtrip metadata")
        expected_ticks=[i*meta["step_ticks"] for i in range(len(selected))]
        if [n["tick"] for n in selected]!=expected_ticks:
            layout_match=False
            if not allow_edited: raise ValueError("reversible note timing differs from roundtrip metadata")
    mapping=inverse_pitch_map(kind); symbols=[]
    for note in selected:
        if note["pitch"] not in mapping: raise ValueError(f"MIDI pitch {note['pitch']} is not valid for {kind} roundtrip mapping")
        symbols.append(mapping[note["pitch"]])
    sequence="".join(symbols); digest=sha256_bytes(sequence.encode()); expected_hash=meta.get("normalized_sequence_sha256") if meta else None
    expected_length=meta.get("sequence_length") if meta else None
    exact=bool(expected_hash) and digest==expected_hash and (expected_length is None or len(sequence)==expected_length) and manifest_midi_match is not False and layout_match
    if meta and not exact and not allow_edited:
        raise ValueError("decoded sequence does not match embedded length/hash; pass --allow-edited to export the edited sequence")
    status=("exact-manifest-matched" if exact and manifest_midi_match else ("exact-embedded" if exact else "edited-or-unverified"))
    return {"sequence":sequence,"kind":kind,"record_name":meta.get("record_name") if meta else None,
            "sequence_length":len(sequence),"normalized_sequence_sha256":digest,"exact":exact,
            "status":status,"profile":meta.get("profile") if meta else "legacy-unverified",
            "metadata_source":"embedded-midi" if embedded else ("manifest" if manifest_meta else "explicit-kind"),
            "manifest_midi_match":manifest_midi_match,"midi":{"path":midi_path.name,"sha256":parsed["sha256"],"format":parsed["format"],"tracks":parsed["tracks"],"ppq":parsed["ppq"],"channel":channel}}


def ensure_real_directory(path: Path) -> None:
    """Create an output parent, but never traverse the requested leaf as a symlink."""
    if path.is_symlink() or (path.exists() and not path.is_dir()):
        raise ValueError(f"output parent must be a real directory, not a symlink or file: {path}")
    path.mkdir(parents=True, exist_ok=True)
    if path.is_symlink() or not path.is_dir():
        raise ValueError(f"output parent changed during creation: {path}")
    ensure_protected_output_ancestry(path)


def ensure_protected_output_ancestry(path: Path) -> None:
    """Reject POSIX output paths another account could rename during staging."""
    if not hasattr(os,"geteuid"):
        return
    owner=os.geteuid(); child=path.resolve(strict=True); child_info=os.stat(child,follow_symlinks=False)
    if not stat.S_ISDIR(child_info.st_mode) or child_info.st_uid!=owner:
        raise ValueError(f"output parent must be a directory owned by the current user: {path}")
    while child.parent!=child:
        parent=child.parent; parent_info=os.stat(parent,follow_symlinks=False)
        shared_writable=parent_info.st_mode & (stat.S_IWGRP|stat.S_IWOTH)
        sticky=parent_info.st_mode & stat.S_ISVTX
        if shared_writable and (not sticky or child_info.st_uid!=owner):
            raise ValueError(f"output ancestry is writable by another account: {parent}")
        child,child_info=parent,parent_info


HAS_ANCHORED_OUTPUT_TRANSACTION = (
    all(function in os.supports_dir_fd for function in (os.open,os.stat,os.rename,os.unlink,os.rmdir))
    and os.listdir in os.supports_fd and os.stat in os.supports_follow_symlinks
)


def supports_anchored_output_transaction() -> bool:
    return HAS_ANCHORED_OUTPUT_TRANSACTION


def output_entry_info(path: Path) -> os.stat_result:
    info=os.stat(path,follow_symlinks=False)
    reparse=getattr(stat,"FILE_ATTRIBUTE_REPARSE_POINT",0)
    if stat.S_ISLNK(info.st_mode) or (reparse and getattr(info,"st_file_attributes",0)&reparse):
        raise ValueError(f"output contains a symlink, junction, or reparse point: {path}")
    return info


def validate_generated_directory(path: Path, expected_names: set[str]) -> os.stat_result:
    info=output_entry_info(path)
    if not stat.S_ISDIR(info.st_mode): raise ValueError(f"output must be a real directory: {path}")
    names={entry.name for entry in path.iterdir()}
    if names!=expected_names: raise ValueError("output contains unknown or unsafe entries")
    for name in names:
        if not stat.S_ISREG(output_entry_info(path/name).st_mode):
            raise ValueError("output contains a non-regular file")
    return info


def cleanup_staged_directory(path: Path, expected_names: set[str]) -> None:
    """Remove only recognized regular staging files; leave anything surprising intact."""
    try:
        info=output_entry_info(path)
        if not stat.S_ISDIR(info.st_mode): return
        entries=list(path.iterdir())
        if any(entry.name not in expected_names for entry in entries): return
        for entry in entries:
            entry_info=output_entry_info(entry)
            if not stat.S_ISREG(entry_info.st_mode): return
        for entry in entries: entry.unlink()
        path.rmdir()
    except (FileNotFoundError,OSError,ValueError):
        return


def commit_generated_directory_portable(staged: Path, final: Path, expected_names: set[str], force: bool) -> None:
    """Conservative new-output fallback for platforms without directory descriptors."""
    parent_info=output_entry_info(final.parent)
    if not stat.S_ISDIR(parent_info.st_mode): raise ValueError(f"output parent must be a real directory: {final.parent}")
    stage_info=validate_generated_directory(staged,expected_names)
    try: final_info=output_entry_info(final)
    except FileNotFoundError: final_info=None
    if final_info is None:
        os.rename(staged,final)
        published=validate_generated_directory(final,expected_names)
        if not stage_info.st_ino or not published.st_ino or (published.st_dev,published.st_ino)!=(stage_info.st_dev,stage_info.st_ino):
            raise ValueError("published output directory identity changed")
        return
    if not stat.S_ISDIR(final_info.st_mode): raise ValueError(f"output must be a real directory: {final}")
    if not force: raise FileExistsError(f"output exists: {final}; pass --force to replace it")
    raise ValueError("--force replacement requires secure directory-descriptor support; choose a new output path on this platform")


def commit_generated_directory(staged: Path, final: Path, expected_names: set[str], force: bool) -> None:
    """Publish a complete staged directory; forced replacement rolls back on failure."""
    if staged.parent.resolve() != final.parent.resolve():
        raise ValueError("staged and final output directories must share a parent")
    ensure_protected_output_ancestry(final.parent)
    if not supports_anchored_output_transaction():
        commit_generated_directory_portable(staged,final,expected_names,force)
        return
    parent_flags=os.O_RDONLY|getattr(os,"O_DIRECTORY",0)|getattr(os,"O_CLOEXEC",0)|getattr(os,"O_NOFOLLOW",0)
    parent_fd=os.open(final.parent,parent_flags)
    stage_fd=None; old_fd=None
    def quarantine_final() -> Path | None:
        try: os.stat(final.name,dir_fd=parent_fd,follow_symlinks=False)
        except FileNotFoundError: return None
        quarantine=Path(tempfile.mkdtemp(prefix=f".{slug(final.name,48)}.quarantine-",dir=final.parent)); quarantine.rmdir()
        os.rename(final.name,quarantine.name,src_dir_fd=parent_fd,dst_dir_fd=parent_fd)
        return quarantine
    try:
        stage_fd=os.open(staged.name,parent_flags,dir_fd=parent_fd)
        stage_info=os.fstat(stage_fd)
        staged_names=set(os.listdir(stage_fd))
        if staged_names != expected_names:
            raise ValueError("staged output does not contain exactly the expected files")
        for name in staged_names:
            info=os.stat(name,dir_fd=stage_fd,follow_symlinks=False)
            if not stat.S_ISREG(info.st_mode): raise ValueError("staged output contains a non-regular file")
        try: final_info=os.stat(final.name,dir_fd=parent_fd,follow_symlinks=False)
        except FileNotFoundError: final_info=None
        if final_info is None:
            os.rename(staged.name,final.name,src_dir_fd=parent_fd,dst_dir_fd=parent_fd)
            published=os.stat(final.name,dir_fd=parent_fd,follow_symlinks=False)
            if (published.st_dev,published.st_ino)!=(stage_info.st_dev,stage_info.st_ino):
                quarantine=quarantine_final()
                raise ValueError(f"published output directory identity changed; unexpected directory quarantined as {quarantine.name}")
            return
        if not stat.S_ISDIR(final_info.st_mode):
            raise ValueError(f"output must be a real directory, not a symlink or file: {final}")
        if not force: raise FileExistsError(f"output exists: {final}; pass --force to replace it")
        old_fd=os.open(final.name,parent_flags,dir_fd=parent_fd)
        old_info=os.fstat(old_fd)
        backup=Path(tempfile.mkdtemp(prefix=f".{slug(final.name,48)}.backup-",dir=final.parent)); backup.rmdir()
        os.rename(final.name,backup.name,src_dir_fd=parent_fd,dst_dir_fd=parent_fd)
        moved=os.stat(backup.name,dir_fd=parent_fd,follow_symlinks=False)
        if (moved.st_dev,moved.st_ino)!=(old_info.st_dev,old_info.st_ino):
            os.rename(backup.name,final.name,src_dir_fd=parent_fd,dst_dir_fd=parent_fd)
            raise ValueError("output directory identity changed during replacement")
        published_ok=False
        try:
            unknown=[]
            for name in os.listdir(old_fd):
                try: info=os.stat(name,dir_fd=old_fd,follow_symlinks=False)
                except FileNotFoundError:
                    unknown.append(name); continue
                if name not in expected_names or not stat.S_ISREG(info.st_mode): unknown.append(name)
            if unknown:
                raise ValueError(f"refusing --force because output contains unknown or unsafe entries: {', '.join(sorted(unknown))}")
            os.rename(staged.name,final.name,src_dir_fd=parent_fd,dst_dir_fd=parent_fd)
            published=os.stat(final.name,dir_fd=parent_fd,follow_symlinks=False)
            if (published.st_dev,published.st_ino)!=(stage_info.st_dev,stage_info.st_ino):
                raise ValueError("published output directory identity changed")
            published_ok=True
        except Exception:
            quarantine=quarantine_final()
            os.rename(backup.name,final.name,src_dir_fd=parent_fd,dst_dir_fd=parent_fd)
            if quarantine is not None:
                raise ValueError(f"output publication failed; unexpected directory quarantined as {quarantine.name}")
            raise
        try:
            remaining=os.listdir(old_fd)
            unsafe=[]
            for name in remaining:
                info=os.stat(name,dir_fd=old_fd,follow_symlinks=False)
                if name not in expected_names or not stat.S_ISREG(info.st_mode): unsafe.append(name)
            current_backup=os.stat(backup.name,dir_fd=parent_fd,follow_symlinks=False)
            if unsafe or (current_backup.st_dev,current_backup.st_ino)!=(old_info.st_dev,old_info.st_ino):
                raise OSError("backup identity or contents changed; leaving it intact")
            for name in remaining: os.unlink(name,dir_fd=old_fd)
            os.close(old_fd); old_fd=None
            os.rmdir(backup.name,dir_fd=parent_fd)
        except OSError as exc:
            print(f"seq2music: warning: new output is live, but backup cleanup failed: {exc}",file=sys.stderr)
    finally:
        if old_fd is not None: os.close(old_fd)
        if stage_fd is not None: os.close(stage_fd)
        os.close(parent_fd)


def cmd_decode(args):
    midi=Path(args.midi); manifest=Path(args.manifest) if args.manifest else None
    result=decode_midi_sequence(midi,manifest,args.kind,args.allow_edited)
    out=Path(args.out) if args.out else midi.with_name(f"{slug(midi.stem)}-decoded")
    ensure_real_directory(out.parent)
    header=re.sub(r"[\x00-\x1f\x7f]+"," ",result.get("record_name") or "decoded-sequence").strip()[:200] or "decoded-sequence"
    sequence=result.pop("sequence")
    fasta=(f">{header}\n"+"\n".join(sequence[i:i+70] for i in range(0,len(sequence),70))+"\n").encode("utf-8")
    expected={"sequence.fasta","decode.report.json"}
    work=Path(tempfile.mkdtemp(prefix=f".{slug(out.name,48)}.tmp-",dir=out.parent))
    fasta_path=work/"sequence.fasta"; report_path=work/"decode.report.json"
    report={**result,"output":{"path":fasta_path.name,"sha256":sha256_bytes(fasta),"format":"canonical-fasta-70-column"},
            "integrity_scope":"carrier consistency with the supplied manifest"}
    report_bytes=(json.dumps(report,indent=2,sort_keys=True)+"\n").encode("utf-8")
    try:
        write_new_bytes(fasta_path,fasta); write_new_bytes(report_path,report_bytes)
        commit_generated_directory(work,out,expected,args.force)
    except Exception:
        cleanup_staged_directory(work,expected)
        raise
    print(json.dumps({"output":str(out.resolve()),"sequence":str((out/'sequence.fasta').resolve()),"report":str((out/'decode.report.json').resolve()),"status":result["status"],"exact":result["exact"]},indent=2))


def osc(phase: float, feature: str) -> float:
    s = math.sin(phase)
    if feature in ("positive","variable","substitution"): return 0.72*s + 0.22*math.sin(2*phase)
    if feature in ("negative","insertion","deletion"): return 0.68*s + 0.18*math.sin(3*phase)
    if feature == "special": return 0.75*s + 0.12*math.sin(4*phase)
    if feature == "contact-chord": return 0.75*s + 0.12*math.sin(2*phase)
    if feature == "purine": return 0.76*s + 0.16*math.sin(2*phase)
    if feature == "pyrimidine": return 0.76*s + 0.13*math.sin(3*phase)
    if feature == "ambiguous": return 0.55*s + 0.10*math.sin(5*phase)
    return s


PITCH_SPELLING=("C","C#","D","D#","E","F","F#","G","G#","A","A#","B")
LETTER_INDEX={"C":0,"D":1,"E":2,"F":3,"G":4,"A":5,"B":6}


def midi_components(note: int) -> tuple[str,int,int,str]:
    if not 0<=note<=127: raise ValueError(f"MIDI note out of range: {note}")
    spelling=PITCH_SPELLING[note%12]; step=spelling[0]; alter=1 if len(spelling)>1 else 0; octave=note//12-1
    return step,alter,octave,f"{spelling}{octave}"


def write_score_svg(path: Path, title: str, events: list[dict]) -> None:
    """Write a complete, accessible grand-staff view with residue/base labels."""
    per_row=24; width=1200; row_height=220; rows=max(1,math.ceil(len(events)/per_row)); height=70+rows*row_height
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',
           f'<title>{html.escape(title)} — musical score</title>',
           f'<desc>Grand-staff notation for {len(events)} Seq2Music events. Residue or base symbols and event numbers appear below each note position. The event ledger remains authoritative.</desc>',
           '<rect width="100%" height="100%" fill="#fff"/><style>text{font-family:system-ui,sans-serif;fill:#172033}.staff{stroke:#172033;stroke-width:1}.note{fill:#4f46e5;stroke:#172033;stroke-width:.8}.ledger{stroke:#172033;stroke-width:1.2}.muted{fill:#64748b}</style>',
           f'<text x="40" y="34" font-size="24" font-weight="700">{html.escape(title)}</text>',
           '<text x="40" y="56" font-size="13" class="muted">Pitch/event-grid notation; symbols and ledger are authoritative. Timing and stereo placement remain in the event ledger.</text>']
    for row in range(rows):
        top=70+row*row_height; treble_bottom=top+70; bass_bottom=top+150
        for bottom in (treble_bottom,bass_bottom):
            for line in range(5):
                y=bottom-line*10; parts.append(f'<line class="staff" x1="55" y1="{y}" x2="1175" y2="{y}"/>')
        parts.append(f'<text x="20" y="{treble_bottom-8}" font-size="34">𝄞</text><text x="20" y="{bass_bottom-8}" font-size="34">𝄢</text>')
        row_events=events[row*per_row:(row+1)*per_row]
        for boundary in range(0,len(row_events)+1,4):
            x=59+boundary*46
            parts.append(f'<line class="staff" x1="{x}" y1="{treble_bottom-40}" x2="{x}" y2="{bass_bottom}" stroke-width="1.4"/>')
        for col,event in enumerate(row_events):
            x=82+col*46; notes=event.get("notes",[])
            if not notes:
                parts.append(f'<rect x="{x-5}" y="{treble_bottom-26}" width="10" height="4" fill="#64748b"><title>Rest, event {event.get("event")}</title></rect>')
            for note_index,note in enumerate(notes):
                step,alter,octave,note_name=midi_components(int(note)); diatonic=octave*7+LETTER_INDEX[step]
                if note>=60: bottom=treble_bottom; base=30
                else: bottom=bass_bottom; base=18
                y=bottom-(diatonic-base)*5; nx=x+note_index*5
                if y>bottom:
                    ledger=bottom+10
                    while ledger<=y+1: parts.append(f'<line class="ledger" x1="{nx-8}" y1="{ledger}" x2="{nx+8}" y2="{ledger}"/>'); ledger+=10
                elif y<bottom-40:
                    ledger=bottom-50
                    while ledger>=y-1: parts.append(f'<line class="ledger" x1="{nx-8}" y1="{ledger}" x2="{nx+8}" y2="{ledger}"/>'); ledger-=10
                label=f'Event {event.get("event")}, symbol {event.get("symbol")}, note {note_name}'
                if alter: parts.append(f'<text x="{nx-12}" y="{y+4}" font-size="12">♯</text>')
                parts.append(f'<ellipse class="note" cx="{nx}" cy="{y}" rx="6" ry="4" transform="rotate(-18 {nx} {y})"><title>{html.escape(label)}</title></ellipse><line class="staff" x1="{nx+5}" y1="{y}" x2="{nx+5}" y2="{y-27}"/>')
            parts.append(f'<text x="{x}" y="{top+185}" font-size="13" text-anchor="middle" font-weight="700">{html.escape(str(event.get("symbol","-")))}</text><text x="{x}" y="{top+202}" font-size="10" text-anchor="middle" class="muted">{event.get("event",row*per_row+col+1)}</text>')
    parts.append('</svg>')
    write_new_text(path,"".join(parts))


def score_roundtrip_metadata(roundtrip: dict | None, events: list[dict]) -> dict | None:
    if not roundtrip: return None
    sequence=roundtrip.get("normalized_sequence")
    if not isinstance(sequence,str): return None
    stream=[{"notes":[int(note) for note in event.get("notes",[])],"symbol":str(event.get("symbol",""))} for event in events]
    return {"profile":SCORE_ROUNDTRIP_VERSION,"class":"exact-normalized-sequence","kind":roundtrip.get("kind"),
            "record_name":roundtrip.get("record_name"),"sequence":sequence,"sequence_length":len(sequence),
            "normalized_sequence_sha256":sha256_bytes(sequence.encode()),"algorithm":ALGORITHM_VERSION,
            "event_count":len(sequence),"note_stream_sha256":sha256_bytes(json.dumps(stream,separators=(",",":")).encode()),
            "scope":"normalized symbols only; embedded sequence bound to parsed pitches and visible lyric labels"}


def write_musicxml(path: Path, title: str, events: list[dict], roundtrip: dict | None = None) -> None:
    root=ET.Element("score-partwise",version="4.0")
    work=ET.SubElement(root,"work"); ET.SubElement(work,"work-title").text=title
    identification=ET.SubElement(root,"identification")
    meta=score_roundtrip_metadata(roundtrip,events)
    if meta:
        misc=ET.SubElement(identification,"miscellaneous")
        field=ET.SubElement(misc,"miscellaneous-field",name=SCORE_ROUNDTRIP_VERSION)
        field.text=base64.b64encode(json.dumps(meta,sort_keys=True,separators=(",",":"),ensure_ascii=True).encode()).decode("ascii")
    part_list=ET.SubElement(root,"part-list"); score_part=ET.SubElement(part_list,"score-part",id="P1"); ET.SubElement(score_part,"part-name").text="Seq2Music"
    part=ET.SubElement(root,"part",id="P1")
    for measure_index,start in enumerate(range(0,max(1,len(events)),4),1):
        measure=ET.SubElement(part,"measure",number=str(measure_index))
        if measure_index==1:
            attrs=ET.SubElement(measure,"attributes"); ET.SubElement(attrs,"divisions").text="1"
            key=ET.SubElement(attrs,"key"); ET.SubElement(key,"fifths").text="0"
            time=ET.SubElement(attrs,"time"); ET.SubElement(time,"beats").text="4"; ET.SubElement(time,"beat-type").text="4"
            ET.SubElement(attrs,"staves").text="2"
            for number,sign,line in (("1","G","2"),("2","F","4")):
                clef=ET.SubElement(attrs,"clef",number=number); ET.SubElement(clef,"sign").text=sign; ET.SubElement(clef,"line").text=line
        for event in events[start:start+4]:
            notes=event.get("notes",[])
            if not notes:
                note_el=ET.SubElement(measure,"note"); ET.SubElement(note_el,"rest"); ET.SubElement(note_el,"duration").text="1"; ET.SubElement(note_el,"type").text="quarter"
                continue
            for note_index,midi_note in enumerate(notes):
                note_el=ET.SubElement(measure,"note")
                if note_index: ET.SubElement(note_el,"chord")
                step,alter,octave,_=midi_components(int(midi_note)); pitch=ET.SubElement(note_el,"pitch")
                ET.SubElement(pitch,"step").text=step
                if alter: ET.SubElement(pitch,"alter").text=str(alter)
                ET.SubElement(pitch,"octave").text=str(octave); ET.SubElement(note_el,"duration").text="1"; ET.SubElement(note_el,"type").text="quarter"
                ET.SubElement(note_el,"staff").text="1" if midi_note>=60 else "2"
                if note_index==0:
                    lyric=ET.SubElement(note_el,"lyric"); ET.SubElement(lyric,"text").text=str(event.get("symbol",""))
    ET.indent(root,space="  ")
    write_new_bytes(path,b'<?xml version="1.0" encoding="UTF-8"?>\n'+ET.tostring(root,encoding="utf-8")+b"\n")


def xml_local(element: ET.Element) -> str:
    return element.tag.rsplit("}",1)[-1]


def bounded_codec_name(value, label: str) -> None:
    if value is not None and (not isinstance(value,str) or len(value)>200 or len(value.encode("utf-8"))>800):
        raise ValueError(f"{label} record_name must be a bounded string")


def exact_int(value, minimum: int, maximum: int) -> bool:
    return isinstance(value,int) and not isinstance(value,bool) and minimum<=value<=maximum


def valid_sha256(value) -> bool:
    return isinstance(value,str) and re.fullmatch(r"[0-9a-f]{64}",value) is not None


def parse_musicxml(path: Path) -> dict:
    raw=read_file_bytes(path,MAX_SCORE_BYTES,"MusicXML")
    if raw.startswith((b"\xff\xfe",b"\xfe\xff")) or b"\x00" in raw[:256]:
        raise ValueError("MusicXML must use UTF-8 or ASCII encoding; UTF-16 is unsupported")
    try: raw.decode("utf-8")
    except UnicodeDecodeError as exc: raise ValueError(f"MusicXML must be valid UTF-8: {exc}")
    if re.search(br"(?i)<!DOCTYPE|<!ENTITY",raw): raise ValueError("MusicXML DTD/entity declarations are not supported")
    element_count=depth=0; root_kind=None; stack=[]; element_stack=[]; current_note=None
    metadata=[]; onsets=[]; onset_symbols=[]; current_onset=None; note_count=0; voices=set()
    partwise_parts=0; timewise_part_ids=set(); timewise_part_total=0; timewise_missing_id=False
    has_time_shift=False; unsupported_note_semantics=set()
    try:
        for event,element in ET.iterparse(io.BytesIO(raw),events=("start","end")):
            if event=="start":
                local=xml_local(element); element_count+=1; depth+=1; stack.append(local); element_stack.append(element)
                if element_count>MAX_SCORE_ELEMENTS: raise ValueError("MusicXML contains too many XML elements")
                if depth>256: raise ValueError("MusicXML nesting exceeds the safety limit")
                if root_kind is None:
                    root_kind=local
                    if root_kind not in ("score-partwise","score-timewise"):
                        raise ValueError("MusicXML root must be score-partwise or score-timewise")
                if local in ("backup","forward"): has_time_shift=True
                if local in ("tie","tied","grace","cue"): unsupported_note_semantics.add(local)
                if local=="part":
                    if root_kind=="score-partwise" and depth==2: partwise_parts+=1
                    elif root_kind=="score-timewise":
                        timewise_part_total+=1; part_id=element.attrib.get("id")
                        if part_id: timewise_part_ids.add(part_id)
                        else: timewise_missing_id=True
                if local=="note":
                    note_count+=1
                    if note_count>MAX_MIDI_EVENTS: raise ValueError("MusicXML contains too many notes")
                    current_note={"rest":False,"chord":False,"step":"","alter":"0","octave":None,"voice":"","lyric":""}
                elif current_note is not None and local=="rest": current_note["rest"]=True
                elif current_note is not None and local=="chord": current_note["chord"]=True
            else:
                local=xml_local(element); value=(element.text or "").strip()
                if current_note is not None:
                    parent=stack[-2] if len(stack)>1 else ""
                    if local=="step" and parent=="pitch": current_note["step"]=value.upper()
                    elif local=="alter" and parent=="pitch": current_note["alter"]=value or "0"
                    elif local=="octave" and parent=="pitch": current_note["octave"]=value
                    elif local=="voice": current_note["voice"]=value
                    elif local=="text" and parent=="lyric": current_note["lyric"]=element.text or ""
                    elif local=="note":
                        if current_note["voice"]: voices.add(current_note["voice"])
                        if current_note["rest"]:
                            current_onset=None
                        elif current_note["octave"] is not None:
                            step=current_note["step"]
                            if step not in LETTER_INDEX: raise ValueError(f"unsupported MusicXML pitch step: {step!r}")
                            try: alter=int(current_note["alter"]); octave=int(current_note["octave"])
                            except (ValueError,TypeError): raise ValueError("MusicXML pitch requires integer alter and octave")
                            if not -2<=alter<=2 or not -1<=octave<=9: raise ValueError("MusicXML pitch alter or octave is outside the supported range")
                            semitones={"C":0,"D":2,"E":4,"F":5,"G":7,"A":9,"B":11}
                            midi=(octave+1)*12+semitones[step]+alter
                            if not 0<=midi<=127: raise ValueError("MusicXML pitch is outside MIDI range")
                            if current_note["chord"]:
                                if current_onset is None: raise ValueError("MusicXML chord note has no preceding onset")
                                current_onset.append(midi)
                            else:
                                current_onset=[midi]; onsets.append(current_onset); onset_symbols.append(current_note["lyric"])
                        current_note=None
                if local=="miscellaneous-field" and element.attrib.get("name")==SCORE_ROUNDTRIP_VERSION:
                    try: payload=base64.b64decode((element.text or "").encode("ascii"),validate=True)
                    except (ValueError,UnicodeEncodeError) as exc: raise ValueError(f"invalid Seq2Music score metadata encoding: {exc}")
                    if len(payload)>MAX_CODEC_BYTES: raise ValueError("Seq2Music score metadata exceeds the codec safety limit")
                    metadata.append(strict_json_loads(payload,"Seq2Music score metadata"))
                if len(element_stack)>1:
                    try: element_stack[-2].remove(element)
                    except ValueError: pass
                element.clear(); element_stack.pop(); stack.pop(); depth-=1
    except ET.ParseError as exc: raise ValueError(f"invalid MusicXML: {exc}")
    part_count=(partwise_parts if root_kind=="score-partwise" else
                (len(timewise_part_ids) if timewise_part_ids and not timewise_missing_id else timewise_part_total))
    if len(metadata)>1: raise ValueError("MusicXML must contain at most one Seq2Music metadata record")
    if not onsets: raise ValueError("MusicXML contains no pitched note onsets")
    return {"onsets":onsets,"metadata":metadata[0] if metadata else None,"sha256":sha256_bytes(raw),"note_count":note_count,
            "part_count":part_count,"voices":sorted(voices),"onset_symbols":onset_symbols,"has_time_shift":has_time_shift,
            "unsupported_note_semantics":sorted(unsupported_note_semantics)}


def decode_musicxml_sequence(score_path: Path, requested_kind: str | None = None, allow_edited: bool = False,
                             polyphony: str = "reject") -> dict:
    parsed=parse_musicxml(score_path); meta=parsed["metadata"]
    if len(parsed["onsets"])>MAX_EVENTS_HARD: raise ValueError(f"score decode exceeds the {MAX_EVENTS_HARD} event safety limit")
    if parsed["has_time_shift"]: raise ValueError("MusicXML backup/forward timing is unsupported")
    if parsed["unsupported_note_semantics"]: raise ValueError(f"MusicXML note semantics are unsupported: {', '.join(parsed['unsupported_note_semantics'])}")
    if len(parsed["voices"])>1: raise ValueError("MusicXML with multiple voices is unsupported")
    if meta and meta.get("profile")!=SCORE_ROUNDTRIP_VERSION: raise ValueError("unsupported Seq2Music score profile")
    if meta and (not isinstance(meta.get("kind"),str) or meta.get("kind") not in ("protein","dna","rna")):
        raise ValueError("score roundtrip metadata contains an invalid or missing molecule kind")
    if meta and meta.get("algorithm")!=ALGORITHM_VERSION: raise ValueError("unsupported Seq2Music score mapping algorithm")
    embedded_kind=meta.get("kind") if meta else None
    if meta and requested_kind and requested_kind!=embedded_kind: raise ValueError("requested kind conflicts with score metadata")
    kind=embedded_kind or requested_kind or "protein"
    if kind not in ("protein","dna","rna"): raise ValueError("score decode kind must be protein, dna, or rna")
    if parsed["part_count"]!=1: raise ValueError("MusicXML must contain exactly one logical part")
    if polyphony not in ("reject","lowest","highest","first"): raise ValueError("invalid score polyphony policy")
    inverse=inverse_pitch_map(kind); allowed=sorted(inverse); symbols=[]; events=[]; polyphonic=0
    for index,onset in enumerate(parsed["onsets"],1):
        if len(onset)>1:
            polyphonic+=1
        source=(max(onset) if polyphony=="highest" else onset[0] if polyphony=="first" else min(onset))
        nearest=min(allowed,key=lambda pitch:(abs(pitch-source),pitch)); distance=abs(nearest-source)
        symbol=inverse[nearest]; symbols.append(symbol)
        events.append({"event":index,"position":index,"symbol":symbol,"source_midi":source,"notes":[nearest],
                       "source_notes":";".join(map(str,onset)),"selected_source_midi":source,"mapped_midi":nearest,
                       "pitch_distance_semitones":distance,
                       "confidence":round(max(0.0,1.0-distance/6.0),6),"exact_pitch":distance==0,
                       "confidence_basis":"linear heuristic 1 - semitone_distance/6, clipped to 0..1; not calibrated",
                       "start":round((index-1)*.2,6),"duration":.16,"velocity":80,"label":"score-import","feature":"score-inference"})
    inferred="".join(symbols); exact=False
    if meta:
        required=(meta.get("class")=="exact-normalized-sequence" and isinstance(meta.get("sequence"),str) and bool(meta.get("sequence"))
                  and exact_int(meta.get("sequence_length"),1,MAX_EVENTS_HARD) and exact_int(meta.get("event_count"),1,MAX_EVENTS_HARD)
                  and valid_sha256(meta.get("normalized_sequence_sha256")) and valid_sha256(meta.get("note_stream_sha256")))
        if not required: raise ValueError("score roundtrip metadata is missing required fields")
        bounded_codec_name(meta.get("record_name"),"score roundtrip metadata")
        declared=meta["sequence"]; validate_sequence(declared,kind)
        parsed_stream=[{"notes":notes,"symbol":symbol} for notes,symbol in zip(parsed["onsets"],parsed["onset_symbols"])]
        stream_hash=sha256_bytes(json.dumps(parsed_stream,separators=(",",":")).encode())
        exact=(len(declared)==meta["sequence_length"]==meta.get("event_count") and sha256_bytes(declared.encode())==meta["normalized_sequence_sha256"]
               and stream_hash==meta["note_stream_sha256"]
               and inferred==declared and parsed["onset_symbols"]==list(declared) and not parsed["voices"]
               and all(len(onset)==1 for onset in parsed["onsets"]))
        if not exact and not allow_edited: raise ValueError("score notes do not match embedded sequence; pass --allow-edited for lossy recovery")
        sequence=declared if exact else inferred
    else: sequence=inferred
    if polyphonic and not exact and polyphony=="reject":
        raise ValueError("non-exact MusicXML contains chords; select --polyphony lowest, highest, or first")
    status="exact-score-embedded" if exact else "lossy-score-transcription"
    if exact:
        for event in events:
            event.pop("confidence",None); event.pop("confidence_basis",None)
            event["mapping_integrity"]="validated against embedded score contract"
    confidence=({"not_applicable":True,"meaning":"exact container integrity; no transcription confidence computed"} if exact else
                {"minimum":round(min(event["confidence"] for event in events),6),"mean":round(sum(event["confidence"] for event in events)/len(events),6),
                 "meaning":"heuristic nearest-palette pitch fit under the selected mapping"})
    return {"sequence":sequence,"kind":kind,"record_name":meta.get("record_name") if meta else "score-import",
            "sequence_length":len(sequence),"normalized_sequence_sha256":sha256_bytes(sequence.encode()),"exact":exact,"status":status,
            "kind_source":"embedded" if meta else ("explicit" if requested_kind else "default-protein"),"confidence":confidence,
            "source":{"path":score_path.name,"sha256":parsed["sha256"],"format":"MusicXML","note_onsets":len(parsed["onsets"]),
                      "polyphonic_onsets":polyphonic,"parts":parsed["part_count"],"voices":parsed["voices"],"has_time_shift":parsed["has_time_shift"],
                      "polyphony_policy":polyphony,"notation_semantics":"ordered pitch/event grid; not an audio-timing or stereo transcript"},"events":events}


def write_wav(path: Path, events: list[dict], sample_rate: int = SAMPLE_RATE, roundtrip: dict | None = None,
              max_bytes: int = DEFAULT_MAX_WAV_MIB*MIB) -> None:
    """Stream a deterministic PCM render; duration is limited by a transparent byte budget, not a preview timer."""
    duration=max((e["start"]+e["duration"] for e in events), default=.1)+.08
    length=int(duration*sample_rate); pcm_bytes=length*4
    estimated_bytes=44+pcm_bytes+(MAX_CODEC_BYTES+9 if roundtrip else 0)
    effective_limit=min(max_bytes,MAX_RIFF_BYTES)
    if estimated_bytes>effective_limit:
        needed=math.ceil(estimated_bytes/MIB)
        raise ValueError(f"rendered WAV needs about {needed} MiB for {duration:.2f}s; current --max-wav-mib is {max_bytes//MIB} (raise it up to {MAX_RIFF_MIB}, or change BPM/step)")
    prepared=[]
    for event in events:
        begin=int(event["start"]*sample_rate); count=max(1,int(event["duration"]*sample_rate)); end=min(length,begin+count)
        pan=float(event.get("pan",0)); notes=[440.0*2**((int(note)-69)/12) for note in event.get("notes",[])]
        prepared.append({"begin":begin,"count":count,"end":end,"lg":math.sqrt((1-pan)/2),"rg":math.sqrt((1+pan)/2),
                         "amp":0.14*(event["velocity"]/127)/max(1,math.sqrt(len(notes))),"frequencies":notes,
                         "feature":event.get("feature","")})
    prepared.sort(key=lambda event:event["begin"])

    def rendered_chunks(chunk_frames: int = 8192):
        active=[]; next_event=0
        for chunk_begin in range(0,length,chunk_frames):
            chunk_end=min(length,chunk_begin+chunk_frames); count=chunk_end-chunk_begin
            while next_event<len(prepared) and prepared[next_event]["begin"]<chunk_end:
                active.append(prepared[next_event]); next_event+=1
            active=[event for event in active if event["end"]>chunk_begin]
            left=[0.0]*count; right=[0.0]*count
            for event in active:
                first=max(chunk_begin,event["begin"]); last=min(chunk_end,event["end"])
                for absolute_index in range(first,last):
                    k=absolute_index-event["begin"]
                    env=min(1.0,k/max(1,.02*sample_rate))*min(1.0,(event["count"]-k)/max(1,.035*sample_rate))
                    value=0.0
                    for frequency in event["frequencies"]:
                        value+=event["amp"]*env*osc(2*math.pi*frequency*k/sample_rate,event["feature"])
                    local=absolute_index-chunk_begin; left[local]+=value*event["lg"]; right[local]+=value*event["rg"]
            yield left,right

    peak=1.0
    for left,right in rendered_chunks():
        peak=max(peak,max(map(abs,left),default=0.0),max(map(abs,right),default=0.0))
    scale=.82/peak; pcm_hasher=hashlib.sha256()
    with open_new_binary(path) as output:
        with wave.open(output,"wb") as out:
            out.setnchannels(2); out.setsampwidth(2); out.setframerate(sample_rate)
            for left,right in rendered_chunks():
                frames=bytearray(len(left)*4)
                for index,(left_value,right_value) in enumerate(zip(left,right)):
                    struct.pack_into("<hh",frames,index*4,int(max(-1,min(1,left_value*scale))*32767),int(max(-1,min(1,right_value*scale))*32767))
                pcm_hasher.update(frames); out.writeframesraw(frames)
        sequence=roundtrip.get("normalized_sequence") if roundtrip else None
        if isinstance(sequence,str):
            metadata={"profile":WAVE_ROUNDTRIP_VERSION,"class":"exact-normalized-sequence","kind":roundtrip.get("kind"),
                      "record_name":roundtrip.get("record_name"),"sequence":sequence,"sequence_length":len(sequence),
                      "normalized_sequence_sha256":sha256_bytes(sequence.encode()),"algorithm":ALGORITHM_VERSION,
                      "pcm_sha256":pcm_hasher.hexdigest(),"sample_rate":sample_rate,"channels":2,"sample_width":2,
                      "frame_count":length,"event_count":len(sequence),
                      "scope":"normalized symbols in a PCM-bound RIFF metadata chunk; not acoustic inference"}
            payload=ROUNDTRIP_MAGIC+json.dumps(metadata,sort_keys=True,separators=(",",":"),ensure_ascii=True).encode()
            if len(payload)>MAX_CODEC_BYTES: raise ValueError("Seq2Music WAV metadata exceeds the codec safety limit")
            chunk=WAVE_CODEC_CHUNK+struct.pack("<I",len(payload))+payload+(b"\x00" if len(payload)%2 else b"")
            output.seek(0,os.SEEK_END); final_size=output.tell()+len(chunk)
            if final_size>effective_limit or final_size>MAX_RIFF_BYTES: raise ValueError("rendered WAV exceeds its declared byte budget")
            output.write(chunk); output.seek(4); output.write(struct.pack("<I",final_size-8))


def parse_wave(path: Path, max_bytes: int = DEFAULT_MAX_INPUT_WAV_MIB*MIB) -> dict:
    limit=min(max_bytes,MAX_RIFF_BYTES); fd,info=open_regular_fd(path)
    try:
        if info.st_size>limit: raise ValueError(f"WAV exceeds the {limit//MIB} MiB safety limit: {path}")
        if info.st_size<12: raise ValueError("audio input must be a RIFF/WAVE file")
        raw=mmap.mmap(fd,0,access=mmap.ACCESS_READ)
    finally: os.close(fd)
    if len(raw)<12 or raw[:4]!=b"RIFF" or raw[8:12]!=b"WAVE": raise ValueError("audio input must be a RIFF/WAVE file")
    declared=struct.unpack("<I",raw[4:8])[0]+8
    if declared!=len(raw): raise ValueError("WAV RIFF size does not match file length")
    view=memoryview(raw); pos=12; fmt_chunks=[]; data_chunks=[]; metadata=[]; chunk_count=0
    while pos<len(raw):
        chunk_count+=1
        if chunk_count>MAX_WAVE_CHUNKS: raise ValueError(f"WAV exceeds the {MAX_WAVE_CHUNKS} chunk safety limit")
        if pos+8>len(raw): raise ValueError("truncated WAV chunk header")
        chunk_id=raw[pos:pos+4]; size=struct.unpack("<I",raw[pos+4:pos+8])[0]; pos+=8
        if pos+size>len(raw): raise ValueError("truncated WAV chunk payload")
        payload=view[pos:pos+size]; pos+=size
        if size%2:
            if pos>=len(raw): raise ValueError("missing WAV chunk padding byte")
            pos+=1
        if chunk_id==b"fmt ": fmt_chunks.append(payload)
        elif chunk_id==b"data": data_chunks.append(payload)
        elif chunk_id==WAVE_CODEC_CHUNK:
            if len(payload)>MAX_CODEC_BYTES: raise ValueError("Seq2Music WAV metadata exceeds the codec safety limit")
            if payload[:len(ROUNDTRIP_MAGIC)].tobytes()!=ROUNDTRIP_MAGIC: raise ValueError("invalid Seq2Music WAV metadata magic")
            metadata.append(strict_json_loads(payload[len(ROUNDTRIP_MAGIC):].tobytes(),"Seq2Music WAV metadata"))
    if len(fmt_chunks)!=1 or len(data_chunks)!=1: raise ValueError("WAV must contain exactly one fmt and one data chunk")
    if len(metadata)>1: raise ValueError("WAV must contain at most one Seq2Music metadata chunk")
    fmt=fmt_chunks[0]
    if len(fmt)<16: raise ValueError("WAV fmt chunk is too short")
    audio_format,channels,sample_rate,byte_rate,block_align,bits=struct.unpack("<HHIIHH",fmt[:16])
    if audio_format!=1: raise ValueError("only uncompressed integer PCM WAV is supported")
    if channels not in (1,2): raise ValueError("WAV must be mono or stereo")
    if not 8000<=sample_rate<=192000: raise ValueError("WAV sample rate must be 8000..192000 Hz")
    if bits not in (8,16,24,32): raise ValueError("WAV sample width must be 8, 16, 24, or 32 bits")
    width=bits//8
    if block_align!=channels*width or byte_rate!=sample_rate*block_align: raise ValueError("WAV fmt byte rate or block alignment is inconsistent")
    frames=data_chunks[0]
    if not frames or len(frames)%block_align: raise ValueError("WAV data length is empty or misaligned")
    frame_count=len(frames)//block_align; duration=frame_count/sample_rate
    return {"raw_sha256":sha256_bytes(raw),"frames":frames,"pcm_sha256":sha256_bytes(frames),"channels":channels,
            "sample_rate":sample_rate,"sample_width":width,"frame_count":frame_count,"duration_seconds":duration,
            "metadata":metadata[0] if metadata else None,"_backing":raw}


def close_wave_snapshot(parsed: dict) -> None:
    frames=parsed.pop("frames",None)
    if isinstance(frames,memoryview): frames.release()
    backing=parsed.pop("_backing",None)
    if backing is not None: backing.close()


def decode_pcm_sample(frames: memoryview, offset: int, width: int) -> float:
    scale=float(1<<(width*8-1))
    raw=frames[offset:offset+width].tobytes()
    if width==1: value=raw[0]-128; denominator=128.0
    elif width==2: value=int.from_bytes(raw,"little",signed=True); denominator=32768.0
    elif width==3:
        value=int.from_bytes(raw+(b"\xff" if raw[2]&0x80 else b"\x00"),"little",signed=True); denominator=8388608.0
    else: value=int.from_bytes(raw,"little",signed=True); denominator=scale
    return value/denominator


def goertzel_power(samples: list[float], sample_rate: float, frequency: float) -> float:
    coefficient=2.0*math.cos(2.0*math.pi*frequency/sample_rate); previous=previous2=0.0
    for sample in samples:
        current=sample+coefficient*previous-previous2; previous2=previous; previous=current
    return max(0.0,previous2*previous2+previous*previous-coefficient*previous*previous2)


def infer_audio_events(parsed: dict, kind: str, window_seconds: float, silence_threshold: float, collapse: bool) -> list[dict]:
    inverse=inverse_pitch_map(kind); allowed=sorted(inverse); stride=max(1,math.ceil(parsed["sample_rate"]/12000))
    rate=parsed["sample_rate"]/stride; frames=parsed["frames"]; width=parsed["sample_width"]; channel_count=parsed["channels"]
    block_align=width*channel_count; decimated_count=math.ceil(parsed["frame_count"]/stride)
    window=max(64,int(round(window_seconds*rate))); events=[]; active=None
    for window_index,start in enumerate(range(0,decimated_count,window),1):
        sample_count=min(window,decimated_count-start)
        if sample_count<max(32,window//3): break
        segments=[[] for _ in range(channel_count)]
        for decimated_index in range(start,start+sample_count):
            frame_index=decimated_index*stride
            if frame_index>=parsed["frame_count"]: break
            frame_offset=frame_index*block_align
            for channel in range(channel_count):
                segments[channel].append(decode_pcm_sample(frames,frame_offset+channel*width,width))
        centered_channels=[]
        for segment in segments:
            mean=sum(segment)/len(segment); centered_channels.append([value-mean for value in segment])
        rms=math.sqrt(sum(value*value for channel in centered_channels for value in channel)/(len(centered_channels)*len(centered_channels[0])))
        if rms<=silence_threshold:
            active=None; continue
        weighted_channels=[]
        for centered in centered_channels:
            weighted_channels.append([value*(.5-.5*math.cos(2*math.pi*i/(len(centered)-1))) for i,value in enumerate(centered)] if len(centered)>1 else centered)
        scores=[]
        for note in allowed:
            frequency=440.0*2**((note-69)/12); score=0.0
            for weighted in weighted_channels:
                score+=goertzel_power(weighted,rate,frequency)
                if frequency*2<rate/2: score+=.30*goertzel_power(weighted,rate,frequency*2)
                if frequency*3<rate/2: score+=.12*goertzel_power(weighted,rate,frequency*3)
            scores.append((score,note))
        scores.sort(reverse=True); best,note=scores[0]; second=scores[1][0] if len(scores)>1 else 0.0
        confidence=max(0.0,min(1.0,(best-second)/best)) if best>0 else 0.0; symbol=inverse[note]
        event_start=start/rate; duration=len(segment)/rate
        if collapse and active is not None and active["symbol"]==symbol:
            active["duration"]=round(event_start+duration-active["start"],6)
            active["confidence"]=round(min(active["confidence"],confidence),6); active["windows"]+=1
            continue
        event={"event":len(events)+1,"position":len(events)+1,"symbol":symbol,"source_midi":note,"notes":[note],
               "confidence":round(confidence,6),"rms":round(rms,6),"windows":1,"start":round(event_start,6),
               "confidence_basis":"relative Goertzel-energy margin between the two strongest allowed pitches; heuristic, not calibrated",
               "duration":round(duration,6),"velocity":80,"label":"audio-import","feature":"audio-inference"}
        events.append(event); active=event
        if len(events)>MAX_EVENTS_HARD: raise ValueError(f"audio inference exceeds the {MAX_EVENTS_HARD} event safety limit")
    if not events: raise ValueError("audio inference found no voiced windows above the silence threshold")
    return events


def decode_wave_snapshot(parsed: dict, audio_path: Path, requested_kind: str | None = None, allow_edited: bool = False,
                         window_seconds: float = .08, silence_threshold: float = .015, collapse: bool = True,
                         max_input_bytes: int = DEFAULT_MAX_INPUT_WAV_MIB*MIB) -> dict:
    meta=parsed["metadata"]; exact=False
    if meta and (not isinstance(meta.get("kind"),str) or meta.get("kind") not in ("protein","dna","rna")):
        raise ValueError("WAV roundtrip metadata contains an invalid or missing molecule kind")
    if meta and meta.get("algorithm")!=ALGORITHM_VERSION: raise ValueError("unsupported Seq2Music WAV mapping algorithm")
    embedded_kind=meta.get("kind") if meta else None
    if requested_kind and embedded_kind and requested_kind!=embedded_kind: raise ValueError("requested kind conflicts with WAV metadata")
    kind=requested_kind or embedded_kind or "protein"
    if kind not in ("protein","dna","rna"): raise ValueError("audio decode kind must be protein, dna, or rna")
    if meta:
        required=(meta.get("profile")==WAVE_ROUNDTRIP_VERSION and meta.get("class")=="exact-normalized-sequence"
                  and isinstance(meta.get("sequence"),str) and bool(meta.get("sequence"))
                  and exact_int(meta.get("sequence_length"),1,MAX_EVENTS_HARD) and exact_int(meta.get("event_count"),1,MAX_EVENTS_HARD)
                  and valid_sha256(meta.get("normalized_sequence_sha256")) and valid_sha256(meta.get("pcm_sha256"))
                  and exact_int(meta.get("sample_rate"),8000,192000) and exact_int(meta.get("channels"),1,2)
                  and exact_int(meta.get("sample_width"),1,4) and exact_int(meta.get("frame_count"),1,MAX_RIFF_BYTES))
        if not required: raise ValueError("WAV roundtrip metadata is missing required fields")
        bounded_codec_name(meta.get("record_name"),"WAV roundtrip metadata")
        sequence=meta["sequence"]; validate_sequence(sequence,kind)
        if len(sequence)>MAX_EVENTS_HARD: raise ValueError(f"WAV embedded sequence exceeds the {MAX_EVENTS_HARD} event safety limit")
        exact=(len(sequence)==meta["sequence_length"]==meta.get("event_count") and sha256_bytes(sequence.encode())==meta["normalized_sequence_sha256"]
               and parsed["pcm_sha256"]==meta["pcm_sha256"] and parsed["sample_rate"]==meta.get("sample_rate")
               and parsed["channels"]==meta.get("channels") and parsed["sample_width"]==meta.get("sample_width")
               and parsed["frame_count"]==meta.get("frame_count"))
        if not exact and not allow_edited: raise ValueError("WAV PCM does not match embedded sequence container; pass --allow-edited for lossy inference")
    if exact:
        events=sequence_events(sequence,kind,.2,meta.get("record_name") or "audio-import"); status="exact-wav-embedded"
    else:
        events=infer_audio_events(parsed,kind,window_seconds,silence_threshold,collapse); sequence="".join(event["symbol"] for event in events); status="lossy-audio-quantized"
    confidences=[event.get("confidence",1.0) for event in events]
    confidence=({"not_applicable":True,"meaning":"exact PCM-bound container integrity; no acoustic confidence computed"} if exact else
                {"minimum":round(min(confidences),6),"mean":round(sum(confidences)/len(confidences),6),
                 "meaning":"heuristic spectral-pitch margin under the selected mapping"})
    return {"sequence":sequence,"kind":kind,"record_name":meta.get("record_name") if meta else "audio-import",
            "sequence_length":len(sequence),"normalized_sequence_sha256":sha256_bytes(sequence.encode()),"exact":exact,"status":status,
            "kind_source":"embedded" if meta else ("explicit" if requested_kind else "default-protein"),"confidence":confidence,
            "source":{"path":audio_path.name,"sha256":parsed["raw_sha256"],"format":"PCM WAV","sample_rate":parsed["sample_rate"],
                      "channels":parsed["channels"],"sample_width":parsed["sample_width"],"duration_seconds":round(parsed["duration_seconds"],6),
                      "input_budget_mib":max_input_bytes//MIB,
                      "analysis":None if exact else {"window_seconds":window_seconds,"silence_threshold":silence_threshold,
                                                     "collapse_adjacent_symbols":collapse,"channel_combination":"sum per-channel spectral powers"}},"events":events}


def decode_wave_sequence(audio_path: Path, requested_kind: str | None = None, allow_edited: bool = False,
                         window_seconds: float = .08, silence_threshold: float = .015, collapse: bool = True,
                         max_input_bytes: int = DEFAULT_MAX_INPUT_WAV_MIB*MIB) -> dict:
    parsed=parse_wave(audio_path,max_input_bytes)
    try:
        return decode_wave_snapshot(parsed,audio_path,requested_kind,allow_edited,window_seconds,silence_threshold,collapse,max_input_bytes)
    finally:
        close_wave_snapshot(parsed)


def write_csv(path: Path, events: list[dict]) -> None:
    fields = sorted({k for e in events for k in e if k != "notes"}) + ["notes"]
    with open_new_text(path) as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n"); w.writeheader()
        for e in events:
            row={**e,"notes":";".join(map(str,e["notes"]))}
            for key,value in row.items():
                if isinstance(value,str) and value.startswith(("=","+","-","@")): row[key]="'"+value
            w.writerow(row)


def write_html(path: Path, title: str, events: list[dict], wav_name: str, score_name: str, summary: str) -> None:
    width=max(800,min(3200,len(events)*7)); height=240
    voiced=[note for e in events for note in e["notes"]]
    max_note=max(voiced,default=81); min_note=min(voiced,default=48)
    max_position=max((int(e["position"]) for e in events),default=1)
    rects=[]
    for e in events:
        x=20+(int(e["position"])-1)*(width-40)/max(1,max_position-1)
        if e.get("label") == "reference": x -= 2.5
        elif e.get("label") == "variant": x += 2.5
        top_note=max(e["notes"],default=min_note)
        y=20+(max_note-top_note)*(height-60)/max(1,max_note-min_note)
        color={"match":"#64748b","substitution":"#ef4444","insertion":"#f59e0b","deletion":"#8b5cf6","no-observed-variation":"#14b8a6","low-information":"#a16207","variable":"#f97316","contact-chord":"#6366f1","no-informative-observation":"#ffffff"}.get(e.get("feature"),"#6d5ef6")
        label=f'Event {e["event"]}, position {e["position"]}, {e.get("label","")}, symbol {e["symbol"]}, feature {e.get("feature","")}'
        rects.append(f'<rect tabindex="0" aria-label="{html.escape(label)}" x="{x:.1f}" y="{y:.1f}" width="5" height="{height-y-25:.1f}" fill="{color}" stroke="#172033" stroke-width="0.35"><title>{html.escape(label)}</title></rect>')
    rows="".join(f"<tr><td>{e['event']}</td><td>{e['position']}</td><td>{html.escape(str(e.get('label','')))}</td><td>{html.escape(str(e['symbol']))}</td><td>{html.escape(str(e.get('feature','')))}</td><td>{','.join(map(str,e['notes'])) or 'silent'}</td><td>{e['velocity']}</td></tr>" for e in events[:500])
    csv_name=path.stem+".events.csv"
    doc=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>{html.escape(title)}</title><style>body{{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#172033}}h1{{color:#4f46e5}}svg{{width:100%;height:auto;border:1px solid #cbd5e1}}rect:focus{{stroke:#000;stroke-width:2}}table{{border-collapse:collapse;width:100%}}th,td{{padding:.45rem;border-bottom:1px solid #ddd;text-align:left}}caption{{text-align:left;font-weight:700;padding:.5rem 0}}code{{background:#eef2ff;padding:.1rem .3rem}}object{{width:100%;min-height:480px;border:1px solid #cbd5e1}}</style></head><body><main><h1>{html.escape(title)}</h1><p>{html.escape(summary)}</p><audio controls preload="metadata" src="{html.escape(wav_name)}">Audio preview unavailable; use the event table and CSV.</audio><h2>Musical notation</h2><p><a href="{html.escape(score_name)}">Open the standalone accessible SVG score</a>. Residue/base labels and the event ledger are authoritative.</p><object data="{html.escape(score_name)}" type="image/svg+xml" aria-label="Seq2Music grand-staff notation">Your browser cannot display the score; use the SVG link.</object><h2>Static event strip</h2><p><strong>Legend:</strong> gray match; red substitution; orange insertion/variable; purple deletion/contact chord; teal no observed variation; white no informative observation. Text labels and the ledger are authoritative.</p><svg viewBox="0 0 {width} {height}" role="img"><title>Seq2Music static event feature strip</title><desc>Each outlined bar is one aligned position or residue; horizontal position follows source position and vertical position follows the highest MIDI pitch. Silent events use an outlined white bar. Focus a bar for its text label.</desc>{''.join(rects)}</svg><h2>Event ledger</h2><p>The table is a preview. <a href="{html.escape(csv_name)}">Download the complete event ledger as CSV</a>.</p><table><caption>First {min(500,len(events))} of {len(events)} events</caption><thead><tr><th>Event</th><th>Position</th><th>Label</th><th>Symbol</th><th>Feature</th><th>MIDI</th><th>Velocity</th></tr></thead><tbody>{rows}</tbody></table></main></body></html>'''
    write_new_text(path,doc)


def make_bundle(input_paths: list[Path], events: list[dict], out: Path, mode: str, title: str, bpm: int, force: bool, details: dict | None = None, parameters: dict | None = None, input_hashes: dict[str,str] | None = None) -> Path:
    parameters={"bpm":bpm, **(parameters or {})}
    d=details or {}
    roundtrip=d.get("roundtrip") if isinstance(d.get("roundtrip"),dict) else None
    public_roundtrip={key:value for key,value in roundtrip.items() if key!="normalized_sequence"} if roundtrip else None
    public_details={**d,"roundtrip":public_roundtrip} if roundtrip else d
    duration=max((e["start"]+e["duration"] for e in events),default=.1)+.08
    if len(events)>MAX_EVENTS_HARD: raise ValueError(f"render has {len(events)} events; hard safety limit is {MAX_EVENTS_HARD}")
    captured_hashes={str(p):((input_hashes or {}).get(str(p)) or sha256_file(p)) for p in input_paths}
    canonical=json.dumps({"mode":mode,"title":title,"inputs":[{"name":p.name,"sha256":captured_hashes[str(p)]} for p in input_paths],"parameters":parameters,"details":details or {},"events":events,"algorithm":ALGORITHM_VERSION,"sample_rate":SAMPLE_RATE,"synth":"seq2music-additive-v2"},sort_keys=True,separators=(",",":"),ensure_ascii=True).encode()
    run_id=sha256_bytes(canonical)[:12]
    bundle=out/f"{slug(input_paths[0].stem)}-{mode}-{run_id}"
    ensure_real_directory(out)
    base=slug(input_paths[0].stem)
    expected=[f"{base}.mid",f"{base}.wav",f"{base}.musicxml",f"{base}.score.svg",f"{base}.events.csv",f"{base}.html",f"{base}.summary.txt",f"{base}.run.json"]
    for target in (bundle/name for name in expected):
        if target.is_symlink() or (target.exists() and not target.is_file()):
            raise ValueError(f"refusing unsafe generated-artifact path: {target}")
    work=Path(tempfile.mkdtemp(prefix=f".{bundle.name}.tmp-",dir=out))
    try:
        midi=work/expected[0]; wav=work/expected[1]; musicxml=work/expected[2]; score_svg=work/expected[3]
        csvp=work/expected[4]; hp=work/expected[5]; sp=work/expected[6]
        write_midi(midi,events,bpm,public_roundtrip)
        write_wav(wav,events,roundtrip=roundtrip,max_bytes=int(parameters.get("max_wav_mib",DEFAULT_MAX_WAV_MIB))*MIB)
        write_musicxml(musicxml,title,events,roundtrip); write_score_svg(score_svg,title,events); write_csv(csvp,events)
        facts=[]
        for key in ("kind","length","changes","alignment_columns","sequences","columns","residues","chains","format"):
            if key in d: facts.append(f"{key.replace('_',' ')}: {d[key]}")
        if mode == "structure" and events:
            values=[e.get("contact_count",0) for e in events]; facts.append(f"nonlocal 8-angstrom contact count range: {min(values)}-{max(values)}")
        audible=sum(bool(e["notes"]) and e["velocity"] > 0 for e in events)
        summary=f"{mode} rendering with {len(events)} ledger events ({audible} audible) over {duration-.08:.2f} seconds under {ALGORITHM_VERSION}."
        if facts: summary += " " + "; ".join(map(str,facts)) + "."
        if roundtrip: summary += " The MIDI, PCM WAV container, and MusicXML score losslessly decode to the normalized sequence under their embedded Seq2Music contracts."
        summary += " Every event is traceable in the complete CSV ledger."
        write_new_text(sp,summary+"\n")
        write_html(hp,title,events,wav.name,score_svg.name,summary)
        manifest={"schema_version":"1.0","run_id":run_id,"plugin":{"name":"seq2music","version":runtime_version()},"algorithm":ALGORITHM_VERSION,
                  "renderer":{"sample_rate_hz":SAMPLE_RATE,"synth":"seq2music-additive-v2"},"mode":mode,
                  "inputs":[{"path":p.name,"sha256":captured_hashes[str(p)]} for p in input_paths],"parameters":parameters,
                  "details":public_details,
                  "verification_scope":"listed generated artifacts compared with the supplied manifest"}
        if public_roundtrip: manifest["roundtrip"]=public_roundtrip
        manifest["artifacts"]=[{"path":p.name,"sha256":sha256_file(p),"bytes":p.stat().st_size} for p in (midi,wav,musicxml,score_svg,csvp,hp,sp)]
        write_new_text(work/expected[7],json.dumps(manifest,indent=2,sort_keys=True)+"\n")
        commit_generated_directory(work,bundle,set(expected),force)
    except Exception:
        cleanup_staged_directory(work,set(expected))
        raise
    return bundle


def cmd_inspect(args):
    p=Path(args.input); records,_=read_fasta_snapshot(p)
    result=[]
    for name,seq in records:
        kind=infer_kind(seq,args.kind); validate_sequence(seq,kind,allow_gaps=args.aligned,allow_iupac=args.aligned)
        result.append({"name":name,"length":len(seq),"kind":kind,"sha256":sha256_bytes(seq.encode())})
    print(json.dumps({"records":result},indent=2))


def cmd_sonify(args):
    p=Path(args.input); records,source_hash=read_fasta_snapshot(p)
    if len(records) != 1: raise ValueError("reversible sequence encoding requires exactly one FASTA record")
    name,seq=records[0]; kind=infer_kind(seq,args.kind); validate_sequence(seq,kind)
    if len(seq)>args.max_events: raise ValueError(f"sequence has {len(seq)} events; raise --max-events explicitly to continue")
    seconds=args.step*60/args.bpm; events=sequence_events(seq,kind,seconds,name)
    mode="encode" if args.command == "encode" else "sonify"
    roundtrip={"profile":ROUNDTRIP_VERSION,"class":"exact-normalized-sequence","mode":mode,"kind":kind,"kind_source":"explicit" if args.kind!="auto" else "inferred",
               "record_name":name,"sequence_length":len(seq),"normalized_sequence":seq,"normalized_sequence_sha256":sha256_bytes(seq.encode()),
               "midi_channel":0,"ordering":"ascending-note-on-tick","scope":"normalized symbols only; not byte-identical FASTA formatting",
               "algorithm":ALGORITHM_VERSION,"ppq":480,"step_beats":args.step,"step_ticks":max(1,int(round(args.step*480))),
               "gate_ticks":max(1,int(round(max(1,int(round(args.step*480)))*.82))),"event_count":len(events)}
    bundle=make_bundle([p],events,Path(args.out),mode,f"Seq2Music reversible sequence: {name}",args.bpm,args.force,
                       {"kind":kind,"length":len(seq),"roundtrip":roundtrip},
                       {"step_beats":args.step,"step_seconds":seconds,"max_events":args.max_events,"max_wav_mib":args.max_wav_mib}, {str(p):source_hash})
    print(bundle.resolve())


def cmd_diff(args):
    pa,pb=Path(args.reference),Path(args.variant); ra,ha=read_fasta_snapshot(pa); rb,hb=read_fasta_snapshot(pb); na,a=ra[0]; nb,b=rb[0]
    kind=infer_kind(a,args.kind)
    if infer_kind(b,args.kind)!=kind: raise ValueError("reference and variant sequence kinds differ")
    validate_sequence(a,kind); validate_sequence(b,kind)
    event_count=len(a)+len(b)
    if event_count>args.max_events: raise ValueError(f"diff would contain {event_count} events; raise --max-events explicitly to continue")
    cells=(len(a)+1)*(len(b)+1)
    if cells>MAX_ALIGNMENT_CELLS: raise ValueError(f"display alignment would require {cells:,} cells; safety limit is {MAX_ALIGNMENT_CELLS:,}")
    seconds=args.step*60/args.bpm; events,aa,bb=diff_events(a,b,kind,seconds)
    if len(events)>args.max_events: raise ValueError(f"aligned diff expands to {len(events)} events; raise --max-events explicitly to continue")
    changes=sum(x!=y for x,y in zip(aa,bb))
    bundle=make_bundle([pa,pb],events,Path(args.out),"diff",f"Seq2Music diff: {na} vs {nb}",args.bpm,args.force,{"kind":kind,"alignment_columns":len(aa),"changes":changes,"reference_alignment":aa,"variant_alignment":bb}, {"step_beats":args.step,"step_seconds":seconds,"max_events":args.max_events,"max_wav_mib":args.max_wav_mib,"alignment":"needleman-wunsch-linear-gap-display-v1"}, {str(pa):ha,str(pb):hb})
    print(bundle.resolve())


def cmd_msa(args):
    p=Path(args.input); records,source_hash=read_fasta_snapshot(p)
    if len(records)<2: raise ValueError("MSA mode requires at least two aligned sequences")
    union="".join(seq.replace("-","") for _,seq in records); kind=infer_kind(union,args.kind)
    for _,seq in records: validate_sequence(seq,kind,allow_gaps=True,allow_iupac=True)
    if len(records[0][1])>args.max_events: raise ValueError(f"alignment has {len(records[0][1])} columns; raise --max-events explicitly to continue")
    seconds=args.step*60/args.bpm; events=msa_events(records,kind,seconds)
    bundle=make_bundle([p],events,Path(args.out),"msa",f"Seq2Music MSA choir: {p.stem}",args.bpm,args.force,{"kind":kind,"sequences":len(records),"columns":len(events)}, {"step_beats":args.step,"step_seconds":seconds,"max_events":args.max_events,"max_wav_mib":args.max_wav_mib,"input_alignment":"pre-aligned","entropy":"unweighted canonical non-gap symbols; all noncanonical IUPAC nucleotide codes excluded as ambiguous"}, {str(p):source_hash})
    print(bundle.resolve())


def cmd_structure(args):
    p=Path(args.input); chains={c.strip() for c in args.chains.split(",") if c.strip()} if args.chains else None
    raw=read_bytes_limited(p); source_hash=sha256_bytes(raw); source_text=raw.decode("utf-8",errors="replace")
    is_cif=p.suffix.lower() in (".cif",".mmcif"); residues=(parse_mmcif(p,chains,source_text,args.model) if is_cif else parse_pdb(p,chains,source_text,args.model))
    if len(residues)>MAX_STRUCTURE_RESIDUES: raise ValueError(f"structure has {len(residues)} residues; hard O(n^2) safety limit is {MAX_STRUCTURE_RESIDUES}")
    if len(residues)>args.max_events: raise ValueError(f"structure has {len(residues)} residues; raise --max-events explicitly to continue")
    seconds=args.step*60/args.bpm; events=structure_events(residues,seconds)
    selected_chains=sorted({r['chain'] for r in residues})
    selected_model=residues[0].get("model",1)
    bundle=make_bundle([p],events,Path(args.out),"structure",f"Seq2Music residue-order geometry trace: {p.stem}",args.bpm,args.force,{"format":"mmcif-subset" if is_cif else "pdb","residues":len(residues),"chains":selected_chains,"model":selected_model}, {"step_beats":args.step,"step_seconds":seconds,"max_events":args.max_events,"max_wav_mib":args.max_wav_mib,"requested_chains":sorted(chains) if chains else None,"selected_chains":selected_chains,"requested_model":args.model,"selected_model":selected_model,"model_policy":"explicit" if args.model is not None else "first encountered","contact_cutoff_angstrom":8.0}, {str(p):source_hash})
    print(bundle.resolve())


def write_import_decode_output(result: dict, out: Path, force: bool, carrier: str) -> None:
    ensure_real_directory(out.parent)
    payload=dict(result); sequence=payload.pop("sequence"); events=payload.pop("events")
    header=re.sub(r"[\x00-\x1f\x7f]+"," ",payload.get("record_name") or f"{carrier}-import").strip()[:200] or f"{carrier}-import"
    fasta=(f">{header}\n"+"\n".join(sequence[i:i+70] for i in range(0,len(sequence),70))+"\n").encode()
    expected={"sequence.fasta","decode.report.json","inference.events.csv","score.svg"}
    work=Path(tempfile.mkdtemp(prefix=f".{slug(out.name,48)}.tmp-",dir=out.parent))
    try:
        fasta_path=work/"sequence.fasta"; report_path=work/"decode.report.json"; csv_path=work/"inference.events.csv"; score_path=work/"score.svg"
        write_new_bytes(fasta_path,fasta); write_csv(csv_path,events); write_score_svg(score_path,f"Seq2Music {carrier} to {payload['kind']} sequence",events)
        report={**payload,"carrier":carrier,"output":{"path":fasta_path.name,"sha256":sha256_bytes(fasta),"format":"canonical-fasta-70-column"},
                "artifacts":{"events":csv_path.name,"score":score_path.name},
                "transcription_scope":"deterministic transcription under the reported mapping and parameters"}
        write_new_text(report_path,json.dumps(report,indent=2,sort_keys=True)+"\n")
        commit_generated_directory(work,out,expected,force)
    except Exception:
        cleanup_staged_directory(work,expected); raise
    print(json.dumps({"output":str(out.resolve()),"sequence":str((out/'sequence.fasta').resolve()),"report":str((out/'decode.report.json').resolve()),
                      "score":str((out/'score.svg').resolve()),"status":result["status"],"exact":result["exact"],"kind":result["kind"]},indent=2))


def cmd_audio_decode(args):
    audio=Path(args.audio)
    result=decode_wave_sequence(audio,args.kind,args.allow_edited,args.window_seconds,args.silence_threshold,args.collapse,args.max_input_mib*MIB)
    out=Path(args.out) if args.out else audio.with_name(f"{slug(audio.stem)}-audio-decoded")
    write_import_decode_output(result,out,args.force,"audio")


def cmd_score_decode(args):
    score=Path(args.score); result=decode_musicxml_sequence(score,args.kind,args.allow_edited,args.polyphony)
    out=Path(args.out) if args.out else score.with_name(f"{slug(score.stem)}-score-decoded")
    write_import_decode_output(result,out,args.force,"score")


def cmd_verify(args):
    p=Path(args.manifest)
    if not p.is_file() or p.is_symlink(): raise ValueError("manifest must be a regular, non-symlink file")
    data=strict_json_loads(read_file_bytes(p,MAX_INPUT_BYTES,"manifest"),"manifest JSON"); failures=[]
    if data.get("schema_version") != "1.0" or not data.get("run_id") or not isinstance(data.get("artifacts"),list) or not data["artifacts"]:
        raise ValueError("manifest is missing required schema_version, run_id, or artifacts")
    if len(data["artifacts"])>MAX_MANIFEST_ARTIFACTS: raise ValueError("manifest contains too many artifacts")
    root=p.parent.resolve(); max_artifact_bytes=getattr(args,"max_artifact_mib",DEFAULT_MAX_VERIFY_MIB)*MIB
    seen=set()
    for a in data["artifacts"]:
        if not isinstance(a,dict) or not isinstance(a.get("path"),str) or not isinstance(a.get("sha256"),str):
            raise ValueError("manifest contains an invalid artifact record")
        if not re.fullmatch(r"[0-9a-f]{64}",a["sha256"]): raise ValueError("manifest contains an invalid artifact SHA-256")
        rel=Path(a["path"])
        if rel.is_absolute() or ".." in rel.parts: raise ValueError(f"unsafe artifact path in manifest: {rel}")
        if str(rel) in seen: raise ValueError(f"duplicate artifact path in manifest: {rel}")
        seen.add(str(rel)); candidate=root/rel
        if candidate.is_symlink(): failures.append(a["path"]); continue
        target=candidate.resolve()
        if target.parent != root or not target.is_file() or target.stat().st_size>max_artifact_bytes:
            failures.append(a["path"]); continue
        declared_bytes=a.get("bytes")
        if declared_bytes is not None and (not isinstance(declared_bytes,int) or isinstance(declared_bytes,bool) or declared_bytes!=target.stat().st_size):
            failures.append(a["path"]); continue
        if sha256_file(target,max_artifact_bytes)!=a["sha256"]: failures.append(a["path"])
    print(json.dumps({"ok":not failures,"failures":failures,"scope":"generated artifacts only"},indent=2))
    if failures: raise SystemExit(1)


def parser() -> argparse.ArgumentParser:
    ap=argparse.ArgumentParser(description="Auditable biological data sonification")
    ap.add_argument("--version",action="version",version=runtime_version())
    sub=ap.add_subparsers(dest="command",required=True)
    i=sub.add_parser("inspect"); i.add_argument("--input",required=True); i.add_argument("--kind",choices=["auto","protein","dna","rna"],default="auto"); i.add_argument("--aligned",action="store_true"); i.set_defaults(func=cmd_inspect)
    def common(p):
        p.add_argument("--out",default="seq2music-output"); p.add_argument("--bpm",type=int,default=100); p.add_argument("--step",type=float,default=.30,help="beats per residue/alignment column (default: 0.30)"); p.add_argument("--max-events",type=int,default=5000)
        p.add_argument("--max-wav-mib",type=int,default=DEFAULT_MAX_WAV_MIB,help=f"rendered WAV byte budget (default: {DEFAULT_MAX_WAV_MIB} MiB; maximum: {MAX_RIFF_MIB})"); p.add_argument("--force",action="store_true")
    for command,help_text in (("encode","reversible sequence-to-music encoding"),("sonify","reversible sequence sonification (compatibility alias)")):
        s=sub.add_parser(command,help=help_text); s.add_argument("--input",required=True); s.add_argument("--kind",choices=["auto","protein","dna","rna"],default="auto"); common(s); s.set_defaults(func=cmd_sonify)
    dec=sub.add_parser("decode",help="recover a normalized sequence from Seq2Music MIDI")
    dec.add_argument("--midi",required=True); dec.add_argument("--manifest")
    dec.add_argument("--kind",choices=["protein","dna","rna"])
    dec.add_argument("--allow-edited",action="store_true",help="export an intentionally edited or metadata-free legacy note stream without claiming exact recovery")
    dec.add_argument("--legacy",dest="allow_edited",action="store_true",help=argparse.SUPPRESS)
    dec.add_argument("--out",help="decoded output directory (default: <midi-stem>-decoded)"); dec.add_argument("--force",action="store_true")
    dec.set_defaults(func=cmd_decode)
    adec=sub.add_parser("audio-decode",help="convert PCM WAV to a sequence; exact for intact Seq2Music WAV, otherwise lossy")
    adec.add_argument("--audio",required=True,help="uncompressed mono/stereo PCM WAV")
    adec.add_argument("--kind",choices=["protein","dna","rna"],help="metadata-free default: protein")
    adec.add_argument("--window-seconds",type=float,default=.08,help="lossy analysis window (default: 0.08 seconds)")
    adec.add_argument("--silence-threshold",type=float,default=.015,help="normalized RMS threshold (default: 0.015)")
    adec.add_argument("--max-input-mib",type=int,default=DEFAULT_MAX_INPUT_WAV_MIB,help=f"WAV input byte budget (default: {DEFAULT_MAX_INPUT_WAV_MIB} MiB; maximum: {MAX_RIFF_MIB})")
    adec.add_argument("--no-collapse",dest="collapse",action="store_false",help="keep repeated adjacent analysis windows as separate symbols")
    adec.add_argument("--allow-edited",action="store_true",help="fall back to lossy inference when embedded Seq2Music PCM integrity fails")
    adec.add_argument("--out"); adec.add_argument("--force",action="store_true"); adec.set_defaults(func=cmd_audio_decode,collapse=True)
    sdec=sub.add_parser("score-decode",help="convert MusicXML to a sequence; exact for intact Seq2Music score, otherwise lossy")
    sdec.add_argument("--score",required=True,help="plain MusicXML score-partwise/timewise file")
    sdec.add_argument("--kind",choices=["protein","dna","rna"],help="metadata-free default: protein")
    sdec.add_argument("--polyphony",choices=["reject","lowest","highest","first"],default="reject")
    sdec.add_argument("--allow-edited",action="store_true",help="fall back to lossy transcription when embedded score integrity fails")
    sdec.add_argument("--out"); sdec.add_argument("--force",action="store_true"); sdec.set_defaults(func=cmd_score_decode)
    d=sub.add_parser("diff"); d.add_argument("--reference",required=True); d.add_argument("--variant",required=True); d.add_argument("--kind",choices=["auto","protein","dna","rna"],default="auto"); common(d); d.set_defaults(func=cmd_diff)
    m=sub.add_parser("msa"); m.add_argument("--input",required=True); m.add_argument("--kind",choices=["auto","protein","dna","rna"],default="auto"); common(m); m.set_defaults(func=cmd_msa)
    st=sub.add_parser("structure"); st.add_argument("--input",required=True,help="local PDB or mmCIF file"); st.add_argument("--chains",help="comma-separated chain IDs"); st.add_argument("--model",type=int,help="positive PDB/mmCIF model number (default: first encountered)"); common(st); st.set_defaults(func=cmd_structure)
    v=sub.add_parser("verify"); v.add_argument("--manifest",required=True); v.add_argument("--max-artifact-mib",type=int,default=DEFAULT_MAX_VERIFY_MIB,help=f"per-artifact hashing budget (default: {DEFAULT_MAX_VERIFY_MIB} MiB; maximum: {MAX_RIFF_MIB})"); v.set_defaults(func=cmd_verify)
    return ap


def main(argv=None):
    try:
        args=parser().parse_args(argv)
        if hasattr(args,"bpm") and not 20 <= args.bpm <= 300: raise ValueError("--bpm must be 20..300")
        if hasattr(args,"step") and not .05 <= args.step <= 4.0: raise ValueError("--step must be 0.05..4.0 beats")
        if hasattr(args,"max_events") and not 1 <= args.max_events <= MAX_EVENTS_HARD: raise ValueError(f"--max-events must be 1..{MAX_EVENTS_HARD}")
        if hasattr(args,"max_wav_mib") and not 1<=args.max_wav_mib<=MAX_RIFF_MIB: raise ValueError(f"--max-wav-mib must be 1..{MAX_RIFF_MIB}")
        if hasattr(args,"max_input_mib") and not 1<=args.max_input_mib<=MAX_RIFF_MIB: raise ValueError(f"--max-input-mib must be 1..{MAX_RIFF_MIB}")
        if hasattr(args,"max_artifact_mib") and not 1<=args.max_artifact_mib<=MAX_RIFF_MIB: raise ValueError(f"--max-artifact-mib must be 1..{MAX_RIFF_MIB}")
        if hasattr(args,"model") and args.model is not None and args.model < 1: raise ValueError("--model must be a positive integer")
        if hasattr(args,"window_seconds") and not .02<=args.window_seconds<=2.0: raise ValueError("--window-seconds must be 0.02..2.0")
        if hasattr(args,"silence_threshold") and not 0<=args.silence_threshold<=.5: raise ValueError("--silence-threshold must be 0..0.5")
        args.func(args)
    except (ValueError,FileNotFoundError,FileExistsError,OSError,json.JSONDecodeError) as exc:
        print(f"seq2music: {exc}",file=sys.stderr); return 2
    return 0


if __name__ == "__main__": raise SystemExit(main())
