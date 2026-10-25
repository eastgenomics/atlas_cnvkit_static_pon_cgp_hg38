"""Unit tests for the pure logic of build_atlas_cnvkit_static_pon.py and derive_pon_members.py.

Synthetic fixtures only: no DNAnexus, Docker or network access, and no real sample identifiers. Integration
(dx downloads, the app job, the Docker build, describe calls) is exercised by the from-scratch reproducibility run.
"""
import copy

import pytest

import build_atlas_cnvkit_static_pon as B
import derive_pon_members as D

SHA = "a" * 64


def member(run, i, specimen=None, sex="F"):
    spec = specimen or f"99{i:03d}S{i:04d}"
    return {"run": run, "sample_id": f"1000{i:05d}-{spec}-26CGPV9-0000-{sex}-0000", "specimen": spec,
            "coverage_file": f"{B.PROJECT}:file-{i:024d}", "coverage_sha256": SHA}


def members(per_run):
    out, i = [], 0
    for run, n in per_run.items():
        for _ in range(n):
            out.append(member(run, i))
            i += 1
    return out


# ---------------------------------------------------------------- member table
def test_parse_members_roundtrip():
    rows = members({"0786_V1": 2})
    text = "\t".join(B.MEMBER_COLUMNS) + "\n" + "".join("\t".join(r[c] for c in B.MEMBER_COLUMNS) + "\n" for r in rows)
    assert B.parse_members(text) == rows


def test_parse_members_rejects_header():
    with pytest.raises(ValueError):
        B.parse_members("run\tsample\n")


def test_validate_members_ok():
    per = {"0786_V1": 2, "0797_V2": 1}
    assert B.validate_members(members(per), n=3, per_run=per) == []


def test_validate_members_pinned_counts_sum():
    assert sum(B.MEMBERS_PER_RUN.values()) == B.N_MEMBERS


@pytest.mark.parametrize("mutate, expect", [
    (lambda rs: rs.append(copy.deepcopy(rs[0])), "duplicate specimen"),
    (lambda rs: rs[1].update(coverage_file=rs[0]["coverage_file"]), "duplicate coverage_file"),
    (lambda rs: rs[0].update(coverage_file="project-Other:file-1"), "not in project"),
    (lambda rs: rs[0].update(coverage_sha256="abc"), "bad sha256"),
    (lambda rs: rs[0].update(specimen="00000S0000"), "not in sample ID"),
    (lambda rs: rs[0].update(run="0811_V3"), "members per run"),
])
def test_validate_members_errors(mutate, expect):
    per = {"0786_V1": 2, "0797_V2": 1}
    rs = members(per)
    mutate(rs)
    errs = B.validate_members(rs, n=3, per_run=per)
    assert any(expect in e for e in errs), errs


def test_validate_members_basename_collision():
    per = {"0786_V1": 1, "0811_V3": 1}
    rs = members(per)
    rs[1]["sample_id"] = rs[0]["sample_id"]           # same sequencing name in two runs
    rs[1]["specimen"] = rs[0]["specimen"] + "x"        # (specimens differ only to isolate the name check)
    assert any("collide" in e for e in B.validate_members(rs, n=2, per_run=per))


def test_lims_sex():
    assert B.lims_sex("1-99999S0001-26CGPV9-0000-M-1") == "M"
    assert B.lims_sex("short-id") == "U"


# ---------------------------------------------------------------- provenance chain
def chain(row, **over):
    """fake describe functions for a valid chain; `over` replaces fields to break it."""
    proj = B.RUN_PROJECTS[row["run"]]
    fid = row["coverage_file"].split(":")[1]
    files = {
        fid: {"createdBy": {"executable": over.get("cov_app", B.COVERAGE_APP), "job": "job-cov"}},
        "file-bam": {"createdBy": {"executable": over.get("chr_app", B.CHR_PREFIX_APP), "job": "job-chr"}},
        "file-abam": {"createdBy": {"executable": over.get("aln_app", B.ALIGNER_APP), "job": "job-aln"}},
    }
    aln_in = dict(B.ALIGNER_INPUTS)
    aln_in["genome_fastagz"] = {"$dnanexus_link": aln_in["genome_fastagz"]}
    aln_in.update(over.get("aln_inputs", {}))
    jobs = {
        "job-cov": {"state": over.get("state", "done"), "project": over.get("cov_project", proj),
                    "input": {"baits": {"$dnanexus_link": over.get("baits", B.BINS_BED)},
                              "sample_id": row["sample_id"], "tumour_bam": {"$dnanexus_link": "file-bam"}}},
        "job-chr": {"project": proj, "input": {"input_bam": {"$dnanexus_link": "file-abam"}, "mode": "add"}},
        "job-aln": {"project": over.get("aln_project", proj), "input": aln_in},
    }
    return (lambda f, p: files[f]), (lambda j: jobs[j])


def test_provenance_ok():
    r = member("0786_V1", 1)
    assert B.check_provenance(r, *chain(r)) == []


@pytest.mark.parametrize("over, expect", [
    ({"cov_app": "app-other"}, "coverage file made by"),
    ({"state": "failed"}, "state failed"),
    ({"cov_project": "project-other"}, "ran in"),
    ({"baits": "file-otherbed"}, "baits"),
    ({"chr_app": "app-other"}, "tumour_bam not made by"),
    ({"aln_app": "app-sentieon-umi"}, "aligned BAM made by"),
    ({"aln_inputs": {"mark_or_remove_duplicate": "remove_duplicate"}}, "mark_or_remove_duplicate"),
    ({"aln_inputs": {"extra_bwa_options": ""}}, "extra_bwa_options"),
    ({"aln_project": "project-other"}, "not in one project"),
])
def test_provenance_breaks(over, expect):
    r = member("0786_V1", 1)
    errs = B.check_provenance(r, *chain(r, **over))
    assert any(expect in e for e in errs), errs


def test_link_id_forms():
    assert B.link_id({"$dnanexus_link": "file-1"}) == "file-1"
    assert B.link_id({"$dnanexus_link": {"id": "file-2", "project": "project-x"}}) == "file-2"
    assert B.link_id("file-3") == "file-3"


# ---------------------------------------------------------------- build
def test_dx_run_args():
    rs = members({"0786_V1": 2})
    a = B.dx_run_args(rs, f"{B.PROJECT}:/pon_build/runs/x")
    assert a[:3] == ["dx", "run", B.PON_APP]
    assert sum(x.startswith("-icoverage_files=") for x in a) == 2
    assert f"-ifasta={B.FASTA['fasta']}" in a and f"-idiploid_parx_genome={B.DIPLOID_PARX_GENOME}" in a
    assert "-igc_clip_min=0.3" in a and "-igc_clip_max=0.7" in a
    assert a[a.index("--destination") + 1] == f"{B.PROJECT}:/pon_build/runs/x" and "--wait" in a


def test_clip_gc_text_matches_app_rule():
    text = ("chromosome\tstart\tend\tgene\tlog2\tdepth\tgc\tspread\n"
            "chr1\t1\t2\tA\t0.1\t10\t0.25\t0.1\n"
            "chr1\t3\t4\tB\t0.1\t10\t0.5\t0.1\n"
            "chr1\t5\t6\tC\t0.1\t10\t0.812345678912\t0.1\n"
            "chr1\t7\t8\tD\t0.1\t10\t0.3\t0.1\n")
    out = B.clip_gc_text(text).splitlines()
    assert [l.split("\t")[6] for l in out[1:]] == ["0.3", "0.5", "0.7", "0.3"]
    assert out[2] == text.splitlines()[2]               # values inside the bounds are untouched (string-identical)


def pon_text(n=3, gcs=None):
    gcs = gcs or [0.5] * n
    return "\t".join(B.PON_COLUMNS) + "\n" + "".join(f"chr1\t{i}\t{i + 1}\tG\t0\t1\t{g}\t0.1\n" for i, g in enumerate(gcs))


def test_check_pon_text():
    assert B.check_pon_text(pon_text(3), n_bins=3) == []
    assert any("bins" in e for e in B.check_pon_text(pon_text(2), n_bins=3))
    assert any("outside" in e for e in B.check_pon_text(pon_text(3, [0.5, 0.71, 0.5]), n_bins=3))
    assert any("columns" in e for e in B.check_pon_text("a\tb\n1\t2\n"))


def test_sex_disagreements():
    text = ("sample\tsex\tX_logratio\tY_logratio\n"
            "1-11111S0001-x-0-F-1.targetcoverage.cnn\tFemale\t0\t-29\n"      # agrees
            "1-11111S0002-x-0-M-1.targetcoverage.cnn\tFemale\t0\t-29\n"      # disagrees
            "1-11111S0003-x-0-F-1.targetcoverage.cnn\tMale\t-0.8\t-1\n"      # disagrees
            "1-11111S0004-x-0-U-1.targetcoverage.cnn\tMale\t-0.8\t-1\n")     # U not compared
    assert B.sex_disagreements(text) == {"11111S0002", "11111S0003"}


def test_specimens_digest_order_independent():
    assert B.specimens_digest({"b", "a"}) == B.specimens_digest(["a", "b"])
    assert len(B.specimens_digest(set())) == 64


def test_check_reference_stats():
    ok = "metric\tvalue\nn_intervals\t14393\nn_samples\t124\ngc_correction_applied\tTrue\n"
    assert B.check_reference_stats(ok) == []
    assert len(B.check_reference_stats(ok.replace("124", "123").replace("True", "False"))) == 2


# ---------------------------------------------------------------- derivation
def mrow(run, sid, spec, stype=D.ELIGIBLE_TYPE, depth="300", excl=""):
    return {"run": run, "sample_id": sid, "specimen": spec, "sample_type": stype, "median_depth": depth,
            "pon_exclusion": excl}


def test_expected_exclusion():
    assert D.expected_exclusion(mrow("r", "s", "p")) == ""
    assert D.expected_exclusion(mrow("r", "s", "p", stype="control (commercial)")) == "control"
    assert D.expected_exclusion(mrow("r", "s", "p", stype="blood germline (x)")) == "non_FFPE"
    assert D.expected_exclusion(mrow("r", "s", "p", depth="99.9")) == "qc_depth_lt100"
    assert D.expected_exclusion(mrow("r", "s", "p", stype="fresh-frozen tumour", depth="50")) == "non_FFPE"


def test_select_members_highest_depth_and_first_on_tie():
    rows = [mrow("V1", "a", "S1", depth="200"), mrow("V3", "a3", "S1", depth="300"),   # V3 deeper
            mrow("V2", "b", "S2", depth="400"), mrow("V3", "b3", "S2", depth="400"),   # tie: first (V2) kept
            mrow("V1", "c", "S3", depth="50", excl="qc_depth_lt100")]                   # excluded
    m, errs = D.select_members(rows)
    assert errs == [] and [r["sample_id"] for r in m] == ["a3", "b"]


def test_select_members_flags_inconsistent_exclusion():
    rows = [mrow("V1", "a", "S1", depth="50")]               # depth < 100 but not excluded
    _, errs = D.select_members(rows)
    assert errs and "qc_depth_lt100" in errs[0]


def test_members_table_sorted_and_columns():
    m = [mrow("V2", "b", "S2"), mrow("V1", "a", "S1")]
    ids = {("V1", "a"): "project-x:file-1", ("V2", "b"): "project-x:file-2"}
    sha = {("V1", "a"): SHA, ("V2", "b"): SHA}
    lines = D.members_table(m, ids, sha).splitlines()
    assert lines[0].split("\t") == B.MEMBER_COLUMNS
    assert [l.split("\t")[1] for l in lines[1:]] == ["a", "b"]


def test_build_app_creates_out_dir(tmp_path, monkeypatch):
    """the app outputs are downloaded into <workdir>/out, which must be created first (first proof run failed here)."""
    calls = []
    monkeypatch.setattr(B, "sh", lambda cmd, capture=False: "job-1\n" if capture else calls.append(cmd))
    monkeypatch.setattr(B, "dx_json", lambda *a: {"state": "done", "output": {
        "cn_reference": {"$dnanexus_link": "file-a"}, "reference_stats": {"$dnanexus_link": "file-b"},
        "pon_member_sex": {"$dnanexus_link": "file-c"}}})

    def fake_download(pf, path):
        assert (tmp_path / "out").is_dir()
        open(path, "w").write("x")
    monkeypatch.setattr(B, "dx_download", fake_download)
    job, out = B.build_app(members({"0786_V1": 1}), str(tmp_path), f"{B.PROJECT}:/pon_build/runs/t")
    assert job == "job-1" and set(out) == {"cn_reference", "reference_stats", "pon_member_sex"}
    assert out["cn_reference"]["file"] == f"{B.PROJECT}:file-a"
