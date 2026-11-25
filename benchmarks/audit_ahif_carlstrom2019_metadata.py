"""Metadata-only audit; no abundance or effect inference."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from collections import Counter

EXPECTED = {
    "metadata53.csv": "a143701e52ac3324e5e1765d4252d2e00f52caf1",
    "metadata61.csv": "1056e7b5148b3d3825483178a29eb73f0ced4cc3",
    "metadata62.csv": "83cd360eb0c8428fad8c3586c74ec7d02c6bc840",
    "analysis_publication.r": "1bf958f12baacc790ac29c04120b266b61270c8f",
}

def audit(root):
    retained = {}
    summary = {}
    for name, expected in EXPECTED.items():
        raw = (root / name).read_bytes()
        sha = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        if sha != expected:
            raise ValueError(f"Source blob mismatch: {name} {sha}")
        if name.endswith(".r"):
            assert "meta62$initial <- droplevels(meta61$initial)" in raw.decode()
            continue
        i = int(name[8:10])
        rows = list(csv.DictReader(raw.decode("utf-8-sig").splitlines(), delimiter=";"))
        assert len({r["Sample"] for r in rows}) == len(rows)
        if i == 53:
            f = [r for r in rows if r["Treatment"] != "Ax" and r["Spray"] != "Inoc"]
            group = lambda r: "|".join(r[k] for k in ("Treatment", "Spray", "Time"))
        else:
            ax = "Axenic" if i == 61 else "Ax"
            f = [r for r in rows if r["Treatment"] != ax and r["Time"] != "t0"]
            group = lambda r: "|".join(r[k] for k in ("Treatment", "Time"))
        retained[i] = f
        summary[i] = {"raw": len(rows), "retained": len(f), "conditions": dict(Counter(map(group, f)))}
    assert [(summary[i]["raw"], summary[i]["retained"]) for i in (53,61,62)] == [(342,299),(373,311),(373,311)]
    mismatches = [
        {"filtered_position_1based": j, "R1": a, "R2": b}
        for j, (a,b) in enumerate(zip(retained[61],retained[62]),1)
        if a["Treatment"] != b["Treatment"]
    ]
    assert len(mismatches) == 1
    m = mismatches[0]
    assert (m["filtered_position_1based"],m["R1"]["Sample"],m["R2"]["Sample"]) == (72,"81","400")
    return {"status":"metadata_only_no_abundance_outcomes","source_blob_sha":EXPECTED,"experiments":summary,"cross_batch_label_mismatches":mismatches}

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--source-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit(args.source_dir), indent=2)+"\n")
