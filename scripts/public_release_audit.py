#!/usr/bin/env python3
"""Fail closed on common public-release leaks and packaging mistakes."""

from __future__ import annotations

import json
import os
import re
import stat
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PLUGIN = ROOT / "plugins" / "seq2music"
MARKETPLACE = ROOT / ".agents" / "plugins" / "marketplace.json"
MANIFEST = PLUGIN / ".codex-plugin" / "plugin.json"
MAX_PUBLIC_FILE_BYTES = 5 * 1024 * 1024

REQUIRED_FILES = (
    ROOT / "README.md",
    ROOT / "CHANGELOG.md",
    ROOT / "LICENSE",
    ROOT / "PRIVACY.md",
    ROOT / "SECURITY.md",
    ROOT / "TERMS.md",
    ROOT / "submission" / "README.md",
    MARKETPLACE,
    MANIFEST,
    PLUGIN / "skills" / "sonify-biomolecules" / "SKILL.md",
    PLUGIN / "scripts" / "seq2music.py",
)

FORBIDDEN_NAMES = {
    ".DS_Store",
    ".env",
    ".env.local",
    ".npmrc",
    ".pypirc",
    "credentials.json",
    "id_rsa",
    "id_ed25519",
}
FORBIDDEN_PARTS = {"__pycache__", ".pytest_cache", ".venv", "venv", "output", "seq2music-output"}
GENERATED_SUFFIXES = (".mid", ".wav", ".musicxml", ".score.svg", ".events.csv", ".run.json", ".summary.txt", ".html")

PRIVATE_PATH_PATTERNS = (
    re.compile("/" + "Users" + r"/[^/\s]+/"),
    re.compile("/" + "home" + r"/[^/\s]+/"),
    re.compile(r"[A-Za-z]:\\Users\\[^\\\s]+\\"),
)
SECRET_PATTERNS = (
    ("private key", re.compile("BEGIN " + r"(?:RSA |EC |OPENSSH )?PRIVATE KEY")),
    ("AWS access key", re.compile(r"\b" + "AKIA" + r"[0-9A-Z]{16}\b")),
    ("GitHub token", re.compile(r"\bgh[pousr]_[A-Za-z0-9_]{30,}\b")),
    ("OpenAI-style key", re.compile(r"\bsk-[A-Za-z0-9_-]{20,}\b")),
    (
        "credential assignment",
        re.compile(
            r"(?i)\b(?:api[_-]?key|access[_-]?token|client[_-]?secret|password)\b"
            r"\s*[:=]\s*['\"][^'\"]{8,}['\"]"
        ),
    ),
)


def load_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid JSON at {path.relative_to(ROOT)}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"JSON root must be an object: {path.relative_to(ROOT)}")
    return value


def scan_tree() -> list[str]:
    failures: list[str] = []
    for path in sorted(ROOT.rglob("*")):
        relative = path.relative_to(ROOT)
        if ".git" in relative.parts:
            continue
        info=os.lstat(path); reparse=getattr(stat,"FILE_ATTRIBUTE_REPARSE_POINT",0)
        if stat.S_ISLNK(info.st_mode) or (reparse and getattr(info,"st_file_attributes",0)&reparse):
            failures.append(f"symlink, junction, or reparse point is not permitted: {relative}")
            continue
        if any(part in FORBIDDEN_PARTS for part in relative.parts):
            failures.append(f"generated or environment directory is present: {relative}")
            continue
        if path.name in FORBIDDEN_NAMES:
            failures.append(f"sensitive or generated filename is present: {relative}")
        if path.name.endswith(GENERATED_SUFFIXES):
            failures.append(f"generated artifact is present: {relative}")
        if not path.is_file():
            continue
        size = path.stat().st_size
        if size > MAX_PUBLIC_FILE_BYTES:
            failures.append(f"file exceeds 5 MiB public limit: {relative} ({size} bytes)")
            continue
        raw = path.read_bytes()
        for marker in (b"/"+b"Users"+b"/",b"/"+b"home"+b"/",b":\\"+b"Users"+b"\\"):
            if marker in raw:
                failures.append(f"private absolute path bytes in {relative}")
        if b"\x00" in raw:
            continue
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            continue
        for pattern in PRIVATE_PATH_PATTERNS:
            match = pattern.search(text)
            if match:
                failures.append(f"private absolute path in {relative}: {match.group(0)!r}")
        for label, pattern in SECRET_PATTERNS:
            if pattern.search(text):
                failures.append(f"possible {label} in {relative}")
    return failures


def validate_layout() -> list[str]:
    failures: list[str] = []
    for path in REQUIRED_FILES:
        if not path.is_file():
            failures.append(f"required file is missing: {path.relative_to(ROOT)}")
    if failures:
        return failures

    try:
        marketplace = load_json(MARKETPLACE)
        manifest = load_json(MANIFEST)
    except ValueError as exc:
        return [str(exc)]

    if marketplace.get("name") != "seq2music":
        failures.append("marketplace name must be seq2music")
    plugins = marketplace.get("plugins")
    if not isinstance(plugins, list) or len(plugins) != 1:
        failures.append("marketplace must contain exactly one plugin")
    else:
        entry = plugins[0]
        source = entry.get("source") if isinstance(entry, dict) else None
        policy = entry.get("policy") if isinstance(entry, dict) else None
        if entry.get("name") != "seq2music":
            failures.append("marketplace plugin name must be seq2music")
        if source != {"source": "local", "path": "./plugins/seq2music"}:
            failures.append("marketplace source must use the portable relative plugin path")
        if policy != {"installation": "AVAILABLE", "authentication": "ON_INSTALL"}:
            failures.append("marketplace policy must be explicit and portable")
        if entry.get("category") != "Scientific Research":
            failures.append("marketplace category must be Scientific Research")

    version = manifest.get("version")
    if not isinstance(version, str) or not re.fullmatch(r"\d+\.\d+\.\d+(?:-[0-9A-Za-z.-]+)?", version):
        failures.append("public plugin version must be plain semantic versioning without a local cachebuster")
    if manifest.get("name") != "seq2music":
        failures.append("plugin manifest name must be seq2music")
    if manifest.get("license") != "MIT":
        failures.append("plugin manifest license must be MIT")

    interface = manifest.get("interface")
    if not isinstance(interface, dict):
        failures.append("plugin interface metadata is missing")
    else:
        if interface.get("category") != "Scientific Research":
            failures.append("plugin category must be Scientific Research")
        if interface.get("capabilities") != ["Interactive", "Read", "Write"]:
            failures.append("plugin capabilities must be Interactive, Read, and Write")
        for key in ("composerIcon", "logo", "logoDark"):
            value = interface.get(key)
            if not isinstance(value, str) or not value.startswith("./"):
                failures.append(f"interface.{key} must be a relative plugin asset path")
                continue
            candidate = (PLUGIN / value).resolve()
            try:
                candidate.relative_to(PLUGIN.resolve())
            except ValueError:
                failures.append(f"interface.{key} escapes the plugin root")
                continue
            if not candidate.is_file():
                failures.append(f"interface.{key} is missing: {value}")
    return failures


def main() -> int:
    failures = validate_layout() + scan_tree()
    if failures:
        for failure in sorted(set(failures)):
            print(f"FAIL: {failure}", file=sys.stderr)
        return 1
    print("Public release audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
