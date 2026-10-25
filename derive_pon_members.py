#!/usr/bin/env python3
"""Re-derive the static PoN member table from the 261011 evaluation sample manifest, and assert it equals the pinned
member table used by build_atlas_cnvkit_static_pon.py.

This documents (and reproduces) how the 124 members were chosen. Both files hold sample identifiers and live on
DNAnexus only; this script downloads them by pinned file ID and checks their sha256.

Membership rules (agreed 2026-10-08, 261011 evaluation, "membership D"):
  1. Eligible = manifest pon_exclusion empty. The manifest's exclusions are re-checked here:
       control       - commercial references with engineered CNVs (sample_type starts with "control")
       non_FFPE      - fresh-frozen WGS-case tumours, their dilutions, blood germline normals
       qc_depth_lt100 - UMI-consensus median depth < 100x (the manifest median_depth column)
     and every eligible sequencing must be "FFPE tumour (assumed)".
  2. One sequencing per specimen: the highest UMI-consensus median_depth; on a tie the first in manifest order
     (strict ">" as in the evaluation's plan_arms.py; there is one tie).
  3. The list is frozen from this UMI-consensus selection. It is NOT re-derived from non-UMI depths (that would swap
     the chosen replicate for 7 specimens), so non-UMI and UMI-consensus results compare pair by pair.
  4. Each member's coverage file is its non-UMI .targetcoverage.cnn (DI-4114 re-processing), cloned with the same file
     ID into 004 /coverage_nonumi/<run>/<sample_id>.targetcoverage.cnn; its sha256 is recorded.

  python3 derive_pon_members.py [--workdir DIR] [--out members.tsv]
Exit 0 only if the derived table is byte-identical to the pinned one (EXPECTED_MEMBERS_SHA256).
"""
import argparse
import os
import sys

import build_atlas_cnvkit_static_pon as B

MANIFEST_FILE = f"{B.PROJECT}:file-JFP5bPj4YggYYbjvJbgg5f40"   # atlas_cnvkit_static_pon_sample_manifest_261011.tsv
EXPECTED_MANIFEST_SHA256 = "d40c834cf0156287f0883ab9103ff63b0ce40cd9c1a0b45251fc3f746106fb3f"
COVERAGE_FOLDER = "/coverage_nonumi"
N_MANIFEST, N_ELIGIBLE = 226, 199
ELIGIBLE_TYPE = "FFPE tumour (assumed)"


def parse_manifest(text):
    lines = [l for l in text.splitlines() if l and not l.startswith("#")]
    hdr = lines[0].split("\t")
    return [dict(zip(hdr, l.split("\t"))) for l in lines[1:]]


def expected_exclusion(r):
    """the exclusion reason the rules give for a manifest row ('' = eligible)."""
    st = r["sample_type"]
    if st.startswith("control"):
        return "control"
    if st != ELIGIBLE_TYPE:
        return "non_FFPE"
    if float(r["median_depth"]) < 100:
        return "qc_depth_lt100"
    return ""


def select_members(rows):
    """returns (members in manifest order, list of rule violations)."""
    errs = [f"{r['run']} {r['sample_id']}: pon_exclusion {r['pon_exclusion']!r} != rule {expected_exclusion(r)!r}"
            for r in rows if r["pon_exclusion"] != expected_exclusion(r)]
    best = {}
    for r in rows:
        if r["pon_exclusion"] == "" and (r["specimen"] not in best
                                         or float(r["median_depth"]) > float(best[r["specimen"]]["median_depth"])):
            best[r["specimen"]] = r
    chosen = {id(r) for r in best.values()}
    return [r for r in rows if id(r) in chosen], errs


def members_table(members, coverage_ids, sha):
    """coverage_ids: (run, sample_id) -> 'project:file'; sha: (run, sample_id) -> sha256."""
    out = ["\t".join(B.MEMBER_COLUMNS)]
    for r in sorted(members, key=lambda r: (r["run"], r["sample_id"])):
        k = (r["run"], r["sample_id"])
        out.append("\t".join([r["run"], r["sample_id"], r["specimen"], coverage_ids[k], sha[k]]))
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", default="/tmp/derive_atlas_cnvkit_static_pon_members")
    ap.add_argument("--out", default=None, help="write the derived table here (default <workdir>/members.tsv)")
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    mp = os.path.join(a.workdir, "manifest.tsv")
    B.dx_download(MANIFEST_FILE, mp)
    if B.sha256(mp) != EXPECTED_MANIFEST_SHA256:
        sys.exit(f"manifest sha256 {B.sha256(mp)} != {EXPECTED_MANIFEST_SHA256}")
    rows = parse_manifest(open(mp).read())
    if len(rows) != N_MANIFEST:
        sys.exit(f"{len(rows)} manifest rows, expected {N_MANIFEST}")
    members, errs = select_members(rows)
    n_elig = sum(r["pon_exclusion"] == "" for r in rows)
    if errs or n_elig != N_ELIGIBLE or len(members) != B.N_MEMBERS:
        sys.exit(f"eligible {n_elig} (expected {N_ELIGIBLE}), members {len(members)} (expected {B.N_MEMBERS})\n  "
                 + "\n  ".join(errs))
    print(f"manifest ok: {len(rows)} sequencings, {n_elig} eligible, {len(members)} members")

    coverage_ids = {}
    for run in sorted({r["run"] for r in members}):
        for d in B.dx_json("find", "data", "--project", B.PROJECT, "--folder", f"{COVERAGE_FOLDER}/{run}",
                           "--name", "*.targetcoverage.cnn", "--class", "file"):
            sid = d["describe"]["name"][:-len(".targetcoverage.cnn")]
            coverage_ids[(run, sid)] = f"{B.PROJECT}:{d['id']}"
    missing = [(r["run"], r["sample_id"]) for r in members if (r["run"], r["sample_id"]) not in coverage_ids]
    if missing:
        sys.exit(f"no non-UMI coverage file for {len(missing)} members: {missing[:3]}")
    sha = {}
    for r in members:
        k = (r["run"], r["sample_id"])
        p = os.path.join(a.workdir, "cnn", r["run"], f"{r['sample_id']}.targetcoverage.cnn")
        os.makedirs(os.path.dirname(p), exist_ok=True)
        B.dx_download(coverage_ids[k], p)
        sha[k] = B.sha256(p)
    out = a.out or os.path.join(a.workdir, "members.tsv")
    open(out, "w").write(members_table(members, coverage_ids, sha))
    got = B.sha256(out)
    if got != B.EXPECTED_MEMBERS_SHA256:
        sys.exit(f"derived member table sha256 {got} != pinned {B.EXPECTED_MEMBERS_SHA256}")
    print(f"REPRODUCED: member table sha256 {got} == EXPECTED_MEMBERS_SHA256 ({out})")


if __name__ == "__main__":
    main()
