#!/usr/bin/env python3
"""Verify original EOG Git blobs preserved without modification in KNK."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "migration/source_inventory_20261122.json"

def git_blob_sha(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()

def main() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if len(manifest["files"]) != 83:
        raise SystemExit("unexpected migration manifest file count")
    for rec in manifest["files"]:
        p = ROOT / rec["destination"]
        if not p.is_file():
            raise SystemExit(f"missing migrated path: {rec['destination']}")
        data = p.read_bytes()
        if len(data) != rec["size_bytes"] or git_blob_sha(data) != rec["git_blob_sha"]:
            raise SystemExit(f"mismatched original bytes: {rec['destination']}")
    print(f"PASS: {len(manifest['files'])} original EOG files have identical Git blob hashes.")

if __name__ == "__main__":
    main()
