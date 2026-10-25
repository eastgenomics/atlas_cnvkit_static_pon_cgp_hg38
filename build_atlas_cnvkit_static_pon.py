#!/usr/bin/env python3
"""Build (and so reproduce) atlas_cnvkit_static_pon_cgp_hg38_v1.0.0.cnn, the static CNVkit panel of normals of the
Atlas CGP pipeline, from pinned inputs, and assert the result against a pinned checksum.

This script is the provenance of the file: every input is pinned by DNAnexus file ID and checksum, the tool by app ID
or image ID, and the output by md5. A clean exit (status 0) means the file was rebuilt from scratch and is
byte-identical to the controlled file. Any mismatch exits non-zero.

  python3 build_atlas_cnvkit_static_pon.py                 # build of record: runs the released DNAnexus app
  python3 build_atlas_cnvkit_static_pon.py --mode local    # independent rebuild in the pinned Docker image

Steps
  1. Download the member table (DNAnexus, pinned by file ID + sha256) and validate it.
  2. Download the 124 member coverage files and check each sha256 against the member table.
  3. Upstream provenance: check that every coverage file was made by the pinned coverage app from the pinned bins
     BED, on a BAM made by the pinned chr-prefix app from a BAM made by the pinned aligner with the pinned settings
     (non-UMI: duplicates marked), all within the expected run project. This is what makes the PoN specific to the
     non-UMI demultiplexing/alignment; any upstream change fails the build.
  4. Build: --mode app runs eggd_cgp-cnvkit-pon with the pinned inputs into a fresh DNAnexus folder and waits;
     --mode local runs `cnvkit.py reference` + gc clipping + `cnvkit.py sex` in the pinned image exactly as the app does.
  5. Check the output: md5 == EXPECTED_PON_MD5; 14,393 bins; every gc in [0.3, 0.7]; 124 samples; member sex calls
     disagree with the LIMS letter for exactly the expected, reviewed members (pinned by count and digest).
  6. Write <workdir>/build_report.json.

Identifiers (sample IDs, specimen numbers) are never stored in this repository: they live in the DNAnexus member
table and only appear in the work directory at run time.

Requires: `dx` (authenticated DNAnexus CLI, access to 004 project-JF3p1f04Yggp0vf3PYjy7Q0F); for --mode local,
`docker` with cgp-cnvkit:1.0.0 loaded from the production image tar (no registry). Python 3 standard library only.
"""
import argparse
import concurrent.futures as cf
import datetime as dt
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys

# ---------------------------------------------------------------- pinned inputs
PROJECT = "project-JF3p1f04Yggp0vf3PYjy7Q0F"   # 004_261008_cgp_cnvkit_pon_comparison (tag never-archive)
MEMBERS_FILE = f"{PROJECT}:file-JFP5bP84YggXfZvJ67xFk4yv"   # atlas_cnvkit_static_pon_members_v1.0.0.tsv
EXPECTED_MEMBERS_SHA256 = "ac14432d50701c1caeee3564bf7d8ddbf7e205f151fe015d3f9c7cb77f7beec8"
MEMBER_COLUMNS = ["run", "sample_id", "specimen", "coverage_file", "coverage_sha256"]
N_MEMBERS = 124
MEMBERS_PER_RUN = {"0734_Val1_50kb": 26, "0786_V1": 30, "0797_V2": 25, "0811_V3": 43}

# the DI-4114 non-UMI re-processing: run -> DNAnexus project the coverage jobs ran in
RUN_PROJECTS = {
    "0734_Val1_50kb": "project-JFGyVy04q4YPkj39QvFF3BBj",
    "0786_V1": "project-JFGyXP04J4q6PYkqFJZJZFY2",
    "0797_V2": "project-JFGyXv84pVZ8YZ55gpJ2J85B",
    "0811_V3": "project-JFGyY9Q4v2gzkj39QvFF3FV8",
}
COVERAGE_APP = "app-J9J2yQj441QYqP572g95kFF1"     # eggd_cgp-cnvkit-coverage 2.0.2
BINS_BED = "file-JBxvjg8433GxJf9Xypj5qyXf"        # atlas_cnvkit_bins_for_cgp_hg38_v1.1.bed (001_Reference)
CHR_PREFIX_APP = "app-J96Bgk84y9xJJ5k2YFbgfbf0"   # eggd_chr_prefix 1.1.0
ALIGNER_APP = "app-GkKxy5j0qjZ42zxJ7Y07Y246"      # sentieon-bwa 5.1.0 (non-UMI)
ALIGNER_INPUTS = {                                # inputs that define the read structure / duplicate handling
    "mark_or_remove_duplicate": "mark_duplicate",
    "genome_fastagz": "file-Gb757784XGyY3FPvkPQ74K9z",
    "genomebwaindex_targz": "file-Gb76f204XGybZ3J6F731xkBp",
    "extra_bwa_options": "-K 10000000",
}

# ---------------------------------------------------------------- pinned build
PON_APP = "app-JFB83yQ40gQb2GPyfGB39PFF"          # eggd_cgp-cnvkit-pon 2.2.0
PON_NAME = "cgp_cnvkit_static_pon_di4093_nonumi"  # app output stem (the controlled name is applied on registration)
FASTA = {  # GIABv3, chr-prefixed, bgzipped (001_Reference /genomes/b38/)
    "fasta": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-Gb73P8Q4XGyX9QZ0g9qp61xV",
    "fasta_fai": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-Gb735V84XGyYkXYPqGqgz14j",
    "fasta_gzi": "project-Fkb6Gkj433GVVvj73J7x8KbV:file-Gb735Vj4XGyj2ZPbzjfYxZxX",
}
FASTA_NAME = "GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz"
DIPLOID_PARX_GENOME = "grch38"
GC_CLIP_MIN, GC_CLIP_MAX = 0.3, 0.7
# cgp-cnvkit:1.0.0 (CNVkit git fc65941d) by immutable local image ID; the image is loaded from the production tar
# (001_Reference project-Fkb6Gkj433GVVvj73J7x8KbV:file-J8j7Vyj45FG1BbK26JQgQY6q), there is no registry to pin.
CNVKIT_IMAGE = "sha256:7999995d0d5270ae44e8b130e9cf259e3c60f77831a8bdb46cf2376db24b8aab"

# ---------------------------------------------------------------- pinned output
EXPECTED_PON_MD5 = "b4a373729d9282ce0a0ed80fd29f6d0a"
EXPECTED_N_BINS = 14393
PON_COLUMNS = ["chromosome", "start", "end", "gene", "log2", "depth", "gc", "spread"]
# members whose `cnvkit.py sex` call disagrees with their LIMS sex letter, reviewed and retained in the 261011
# evaluation (a male with a tumour chrX gain; a probable 47,XXY). Pinned by count and by the sha256 of the sorted
# specimen numbers joined with ";" (the identifiers themselves stay on DNAnexus).
EXPECTED_N_SEX_DISAGREEMENTS = 2
EXPECTED_SEX_DISAGREEMENTS_SHA256 = "6f7ee63765c2d280fc01cc0b603cf88f2094a0ca73ba1a9fbac16fe662ab0a2f"


# ================================================================ helpers
def sh(cmd, capture=False):
    r = subprocess.run(cmd, shell=isinstance(cmd, str), check=True, text=True,
                       stdout=subprocess.PIPE if capture else None)
    return r.stdout if capture else None


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def dx_download(project_file, out_path):
    sh(["dx", "download", "--no-progress", "-f", "-o", out_path, project_file])


def dx_json(*args):
    return json.loads(sh(["dx", *args, "--json"], capture=True))


def link_id(v):
    """file ID from a $dnanexus_link value (string or {id, project})."""
    v = v["$dnanexus_link"] if isinstance(v, dict) and "$dnanexus_link" in v else v
    return v if isinstance(v, str) else v["id"]


def lims_sex(sample_id):
    """LIMS sex letter = field 5 of the '-'-split sample ID (M, F or U)."""
    parts = sample_id.split("-")
    return parts[4] if len(parts) > 4 else "U"


# ================================================================ member table
def parse_members(text):
    lines = text.splitlines()
    hdr = lines[0].split("\t")
    if hdr != MEMBER_COLUMNS:
        raise ValueError(f"member table header {hdr} != {MEMBER_COLUMNS}")
    return [dict(zip(hdr, l.split("\t"))) for l in lines[1:] if l]


def validate_members(rows, n=N_MEMBERS, per_run=None):
    per_run = MEMBERS_PER_RUN if per_run is None else per_run
    errs = []
    if len(rows) != n:
        errs.append(f"{len(rows)} members, expected {n}")
    for col in ("specimen", "coverage_file"):
        vals = [r[col] for r in rows]
        if len(set(vals)) != len(vals):
            errs.append(f"duplicate {col}")
    names = [f"{r['sample_id']}.targetcoverage.cnn" for r in rows]
    if len(set(names)) != len(names):
        errs.append("duplicate coverage file names (would collide in cnn_files/)")
    counts = {}
    for r in rows:
        counts[r["run"]] = counts.get(r["run"], 0) + 1
        if not r["coverage_file"].startswith(PROJECT + ":file-"):
            errs.append(f"{r['sample_id']}: coverage_file not in {PROJECT}")
        if len(r["coverage_sha256"]) != 64:
            errs.append(f"{r['sample_id']}: bad sha256")
        if r["specimen"] not in r["sample_id"]:
            errs.append(f"{r['sample_id']}: specimen {r['specimen']} not in sample ID")
    if counts != per_run:
        errs.append(f"members per run {counts} != {per_run}")
    return errs


# ================================================================ upstream provenance
def check_provenance(row, describe_file, describe_job):
    """Return a list of problems with the chain coverage file <- coverage job <- chr-prefix job <- aligner job.
    describe_file(file_id, project) and describe_job(job_id) return the DNAnexus describe dicts (injected for tests)."""
    errs = []
    sid = row["sample_id"]
    fproj, fid = row["coverage_file"].split(":")
    f = describe_file(fid, fproj)
    cb = f.get("createdBy", {})
    if cb.get("executable") != COVERAGE_APP:
        return [f"{sid}: coverage file made by {cb.get('executable')}, expected {COVERAGE_APP}"]
    cj = describe_job(cb["job"])
    if cj.get("state") != "done":
        errs.append(f"{sid}: coverage job {cb['job']} state {cj.get('state')}")
    if cj.get("project") != RUN_PROJECTS.get(row["run"]):
        errs.append(f"{sid}: coverage job ran in {cj.get('project')}, expected {RUN_PROJECTS.get(row['run'])}")
    if link_id(cj["input"].get("baits", "")) != BINS_BED:
        errs.append(f"{sid}: coverage baits {cj['input'].get('baits')} != {BINS_BED}")
    if cj["input"].get("sample_id") not in (None, sid):
        errs.append(f"{sid}: coverage job sample_id {cj['input'].get('sample_id')}")
    bam = describe_file(link_id(cj["input"]["tumour_bam"]), cj["project"])
    pj = describe_job(bam["createdBy"]["job"]) if bam.get("createdBy", {}).get("executable") == CHR_PREFIX_APP else None
    if pj is None:
        return errs + [f"{sid}: tumour_bam not made by {CHR_PREFIX_APP}"]
    abam = describe_file(link_id(pj["input"]["input_bam"]), pj["project"])
    if abam.get("createdBy", {}).get("executable") != ALIGNER_APP:
        return errs + [f"{sid}: aligned BAM made by {abam.get('createdBy', {}).get('executable')}, expected {ALIGNER_APP}"]
    aj = describe_job(abam["createdBy"]["job"])
    for k, v in ALIGNER_INPUTS.items():
        got = aj["input"].get(k)
        got = link_id(got) if isinstance(got, dict) else got
        if got != v:
            errs.append(f"{sid}: aligner input {k}={got!r}, expected {v!r}")
    if not (cj["project"] == pj["project"] == aj["project"]):
        errs.append(f"{sid}: coverage/chr-prefix/aligner jobs not in one project")
    return errs


def provenance_all(rows, jobs=16):
    cache = {}

    def df(fid, proj):
        return dx_json("describe", f"{proj}:{fid}")

    def dj(jid):
        if jid not in cache:
            cache[jid] = dx_json("describe", jid)
        return cache[jid]

    with cf.ThreadPoolExecutor(jobs) as ex:
        res = list(ex.map(lambda r: check_provenance(r, df, dj), rows))
    return [e for r in res for e in r]


# ================================================================ build
def dx_run_args(rows, destination):
    args = ["dx", "run", PON_APP]
    for r in rows:
        args.append(f"-icoverage_files={r['coverage_file']}")
    args += [f"-ipon_name={PON_NAME}"] + [f"-i{k}={v}" for k, v in FASTA.items()]
    args += [f"-idiploid_parx_genome={DIPLOID_PARX_GENOME}", f"-igc_clip_min={GC_CLIP_MIN}", f"-igc_clip_max={GC_CLIP_MAX}",
             "--destination", destination, "--name", "atlas_cnvkit_static_pon_v1.0.0_build",
             "--property", "build=build_atlas_cnvkit_static_pon.py", "--brief", "-y", "--wait"]
    return args


def build_app(rows, workdir, destination):
    sh(["dx", "mkdir", "-p", destination])
    job = sh(dx_run_args(rows, destination), capture=True).strip().splitlines()[-1]
    d = dx_json("describe", job)
    if d["state"] != "done":
        raise RuntimeError(f"{job} state {d['state']}: {d.get('failureMessage')}")
    out = {}
    for key, name in (("cn_reference", f"{PON_NAME}.cnn"), ("reference_stats", "reference_stats.tsv"),
                      ("pon_member_sex", "pon_member_sex.tsv")):
        fid = link_id(d["output"][key])
        p = os.path.join(workdir, "out", name)
        dx_download(f"{destination.split(':')[0]}:{fid}", p)
        out[key] = {"file": f"{destination.split(':')[0]}:{fid}", "path": p}
    return job, out


def clip_gc_text(text, lo=GC_CLIP_MIN, hi=GC_CLIP_MAX):
    """The eggd_cgp-cnvkit-pon 2.2.0 clipping, byte for byte: values outside [lo, hi] set to the bound (%.10g)."""
    lines = text.splitlines()
    gi = lines[0].split("\t").index("gc")
    out = [lines[0]]
    for l in lines[1:]:
        c = l.split("\t")
        g = float(c[gi])
        if g < lo or g > hi:
            c[gi] = f"{min(max(g, lo), hi):.10g}"
        out.append("\t".join(c))
    return "\n".join(out) + "\n"


def build_local(rows, workdir, fasta_dir):
    cnn = os.path.join(workdir, "cnn")
    work = os.path.join(workdir, "local")
    os.makedirs(os.path.join(work, "cnn_files"), exist_ok=True)
    for r in rows:
        shutil.copy(os.path.join(cnn, r["run"], f"{r['sample_id']}.targetcoverage.cnn"), os.path.join(work, "cnn_files"))
    fa = os.path.abspath(os.path.join(fasta_dir, FASTA_NAME))
    for ext in ("", ".fai", ".gzi"):
        if not os.path.exists(fa + ext):
            raise FileNotFoundError(fa + ext)
    uid = f"{os.getuid()}:{os.getgid()}"
    base = ["docker", "run", "--rm", "--user", uid, "-v", f"{os.path.abspath(work)}:/work", "-w", "/work"]
    sh(" ".join(base + ["-v", f"{os.path.dirname(fa)}:/ref:ro", CNVKIT_IMAGE, "sh", "-c",
                        f"'cnvkit.py reference cnn_files/*.cnn --output ref.cnn --fasta /ref/{FASTA_NAME} "
                        f"--diploid-parx-genome {DIPLOID_PARX_GENOME}'"]))
    os.makedirs(os.path.join(workdir, "out"), exist_ok=True)
    pon = os.path.join(workdir, "out", f"{PON_NAME}.cnn")
    with open(os.path.join(work, "ref.cnn")) as f:
        open(pon, "w").write(clip_gc_text(f.read()))
    sh(" ".join(base + [CNVKIT_IMAGE, "sh", "-c",
                        f"'cd cnn_files && cnvkit.py sex ./*.cnn --diploid-parx-genome {DIPLOID_PARX_GENOME} -o ../pon_member_sex.tsv'"]))
    sex = os.path.join(workdir, "out", "pon_member_sex.tsv")
    shutil.copy(os.path.join(work, "pon_member_sex.tsv"), sex)
    return {"cn_reference": {"path": pon}, "pon_member_sex": {"path": sex}}


# ================================================================ output checks
def check_pon_text(text, n_bins=EXPECTED_N_BINS, lo=GC_CLIP_MIN, hi=GC_CLIP_MAX):
    errs = []
    lines = text.splitlines()
    hdr = lines[0].split("\t")
    if hdr != PON_COLUMNS:
        return [f"PoN columns {hdr} != {PON_COLUMNS}"]
    rows = [l.split("\t") for l in lines[1:] if l]
    if len(rows) != n_bins:
        errs.append(f"{len(rows)} bins, expected {n_bins}")
    gi = hdr.index("gc")
    bad = sum(1 for r in rows if not (lo <= float(r[gi]) <= hi))
    if bad:
        errs.append(f"{bad} gc values outside [{lo}, {hi}]")
    return errs


def specimens_digest(specimens):
    return hashlib.sha256(";".join(sorted(specimens)).encode()).hexdigest()


def sex_disagreements(sex_text):
    """specimens whose cnvkit.py sex call disagrees with a LIMS M/F letter (U is not compared)."""
    lines = sex_text.splitlines()
    hdr = lines[0].split("\t")
    si, xi = hdr.index("sample"), hdr.index("sex")
    out = set()
    for l in lines[1:]:
        c = l.split("\t")
        sid = os.path.basename(c[si]).replace(".targetcoverage.cnn", "")
        lims = lims_sex(sid)
        if lims in ("M", "F") and (c[xi] == "Male") != (lims == "M"):
            out.add(sid.split("-")[1])
    return out


def check_reference_stats(text, n_samples=N_MEMBERS, n_bins=EXPECTED_N_BINS):
    kv = dict(l.split("\t") for l in text.splitlines()[1:] if l)
    errs = []
    if kv.get("n_samples") != str(n_samples):
        errs.append(f"reference_stats n_samples {kv.get('n_samples')}")
    if kv.get("n_intervals") != str(n_bins):
        errs.append(f"reference_stats n_intervals {kv.get('n_intervals')}")
    if kv.get("gc_correction_applied") != "True":
        errs.append(f"reference_stats gc_correction_applied {kv.get('gc_correction_applied')}")
    return errs


# ================================================================ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["app", "local"], default="app")
    ap.add_argument("--workdir", default="/tmp/build_atlas_cnvkit_static_pon")
    ap.add_argument("--destination", help="app mode: DNAnexus folder for the job outputs "
                    f"(default {PROJECT}:/pon_build/runs/<UTC timestamp>)")
    ap.add_argument("--fasta-dir", help="local mode: directory holding the GIABv3 FASTA and its .fai/.gzi "
                    "(default: downloaded into <workdir>/ref)")
    ap.add_argument("--jobs", type=int, default=16, help="parallel downloads / describes")
    a = ap.parse_args()

    if os.path.exists(a.workdir) and os.listdir(a.workdir):
        sys.exit(f"--workdir {a.workdir} is not empty: use a fresh directory (no cached inputs)")
    os.makedirs(a.workdir, exist_ok=True)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    report = {"started_utc": stamp, "mode": a.mode, "script": os.path.basename(__file__),
              "dx_whoami": sh(["dx", "whoami"], capture=True).strip(), "checks": {}}

    # 1. member table
    mp = os.path.join(a.workdir, "members.tsv")
    dx_download(MEMBERS_FILE, mp)
    got = sha256(mp)
    if got != EXPECTED_MEMBERS_SHA256:
        sys.exit(f"member table sha256 {got} != {EXPECTED_MEMBERS_SHA256}")
    rows = parse_members(open(mp).read())
    errs = validate_members(rows)
    if errs:
        sys.exit("member table invalid:\n  " + "\n  ".join(errs))
    report["checks"]["members"] = f"ok: {len(rows)} members, sha256 {got}"
    print(report["checks"]["members"])

    # 2. coverage files
    def fetch(r):
        p = os.path.join(a.workdir, "cnn", r["run"], f"{r['sample_id']}.targetcoverage.cnn")
        os.makedirs(os.path.dirname(p), exist_ok=True)
        dx_download(r["coverage_file"], p)
        return None if sha256(p) == r["coverage_sha256"] else r["sample_id"]
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        bad = [x for x in ex.map(fetch, rows) if x]
    if bad:
        sys.exit(f"{len(bad)} coverage files fail sha256: {bad[:5]}")
    report["checks"]["coverage_sha256"] = f"ok: {len(rows)}/{len(rows)}"
    print("coverage files:", report["checks"]["coverage_sha256"])

    # 3. upstream provenance
    errs = provenance_all(rows, a.jobs)
    if errs:
        sys.exit(f"upstream provenance: {len(errs)} problems:\n  " + "\n  ".join(errs[:20]))
    report["checks"]["upstream_provenance"] = (f"ok: {len(rows)} chains {ALIGNER_APP} -> {CHR_PREFIX_APP} -> "
                                               f"{COVERAGE_APP} (baits {BINS_BED}); aligner {ALIGNER_INPUTS}")
    print("upstream provenance:", report["checks"]["upstream_provenance"])

    # 4. build
    if a.mode == "app":
        dest = a.destination or f"{PROJECT}:/pon_build/runs/{stamp}"
        job, out = build_app(rows, a.workdir, dest)
        report.update(job=job, destination=dest, app=PON_APP)
        print("built by", job)
    else:
        fdir = a.fasta_dir
        if not fdir:
            fdir = os.path.join(a.workdir, "ref")
            os.makedirs(fdir, exist_ok=True)
            for k, ext in (("fasta", ""), ("fasta_fai", ".fai"), ("fasta_gzi", ".gzi")):
                dx_download(FASTA[k], os.path.join(fdir, FASTA_NAME + ext))
        img = sh(["docker", "image", "inspect", "--format", "{{.Id}}", CNVKIT_IMAGE], capture=True).strip()
        if img != CNVKIT_IMAGE:
            sys.exit(f"docker image {img} != {CNVKIT_IMAGE}")
        out = build_local(rows, a.workdir, fdir)
        report.update(image=CNVKIT_IMAGE)
    report["outputs"] = out

    # 5. output checks
    pon = out["cn_reference"]["path"]
    got = md5(pon)
    errs = [] if got == EXPECTED_PON_MD5 else [f"PoN md5 {got} != {EXPECTED_PON_MD5}"]
    errs += check_pon_text(open(pon).read())
    if "reference_stats" in out:
        errs += check_reference_stats(open(out["reference_stats"]["path"]).read())
    dis = sex_disagreements(open(out["pon_member_sex"]["path"]).read())
    if len(dis) != EXPECTED_N_SEX_DISAGREEMENTS or specimens_digest(dis) != EXPECTED_SEX_DISAGREEMENTS_SHA256:
        errs.append(f"member sex disagreements: {len(dis)} (expected {EXPECTED_N_SEX_DISAGREEMENTS}), digest "
                    f"{specimens_digest(dis)} (expected {EXPECTED_SEX_DISAGREEMENTS_SHA256})")
    report["checks"]["output"] = errs or (f"ok: md5 {got}; {EXPECTED_N_BINS} bins; gc in [{GC_CLIP_MIN}, {GC_CLIP_MAX}]; "
                                          f"{len(dis)} member sex disagreements, the expected reviewed set")
    report["finished_utc"] = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    json.dump(report, open(os.path.join(a.workdir, "build_report.json"), "w"), indent=1)
    if errs:
        sys.exit("output checks failed:\n  " + "\n  ".join(errs))
    print("output:", report["checks"]["output"])
    print(f"REPRODUCED: {PON_NAME}.cnn md5 {got} == EXPECTED_PON_MD5")


if __name__ == "__main__":
    main()
