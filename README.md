# atlas_cnvkit_static_pon_cgp_hg38

Builds — and so reproduces — **`atlas_cnvkit_static_pon_cgp_hg38_v1.0.0.cnn`**: the static CNVkit panel of normals
(PoN, `cn_reference`) of the Atlas CGP pipeline's CNVkit arm (Twist CGP panel + 50 kb SNP backbone, SHv1, GRCh38),
for samples processed with the **non-UMI** demultiplexing and alignment.

The build script *is* the provenance of the file: every input is pinned by DNAnexus file ID and checksum, the tool by
app ID or Docker image ID, and the output by md5 (`EXPECTED_PON_MD5 = b4a373729d9282ce0a0ed80fd29f6d0a`). A clean exit
means the file was rebuilt from scratch and is byte-identical to the controlled file.

Controlled document (Confluence, Documentation Vault):
[atlas_cnvkit_static_pon_cgp_hg38_v1.0.0.cnn](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4854120471).
Jira: DI-4093 (PoN), DI-4118 (validation on the non-UMI data), DI-4114 (non-UMI re-processing).

## What this builds

| | |
|---|---|
| Tool | `eggd_cgp-cnvkit-pon` 2.2.0 (`app-JFB83yQ40gQb2GPyfGB39PFF`) = `cnvkit.py reference` in `cgp-cnvkit:1.0.0` (CNVkit git fc65941d), then `gc` clipped to [0.3, 0.7] |
| Members | 124 FFPE specimens, one sequencing each: runs 0734 Val1 50 kb (26), 0786 V1 (30), 0797 V2 (25), 0811 V3 (43) |
| Member coverage | non-UMI `.targetcoverage.cnn` from `eggd_cgp-cnvkit-coverage` 2.0.2 with bins BED `atlas_cnvkit_bins_for_cgp_hg38_v1.1.bed` |
| Options | `--fasta` GIABv3 (chr-prefixed, bgzipped), `--diploid-parx-genome grch38`, `gc_clip_min 0.3`, `gc_clip_max 0.7` |
| Output | 14,393 bins, columns `chromosome start end gene log2 depth gc spread`, md5 `b4a373729d9282ce0a0ed80fd29f6d0a` |

**Identifiers are not in this repository.** The member table (sample IDs, specimen numbers, coverage file IDs and
sha256) and the evaluation sample manifest are on DNAnexus in 004 `project-JF3p1f04Yggp0vf3PYjy7Q0F`
(`/pon_build/inputs/`, tag never-archive) and are pinned here by file ID and sha256 (`MEMBERS_FILE`,
`MANIFEST_FILE`). The coverage inputs are pinned to their 004 `/coverage_nonumi/` copies.

### Decisions encoded

| # | Decision | Where |
|---|---|---|
| 1 | Static PoN of FFPE specimens (not per run) | 261011 evaluation, static vs per-run page |
| 2 | Members: FFPE only; controls, non-FFPE and UMI-consensus median depth < 100× excluded; one sequencing per specimen, the highest depth (first in manifest order on a tie) | `derive_pon_members.py` (`expected_exclusion`, `select_members`) |
| 3 | Member list frozen from that UMI-consensus selection, not re-derived from non-UMI depths (would swap 7 replicates), so non-UMI and UMI-consensus results compare pair by pair | `derive_pon_members.py` docstring; pinned member table |
| 4 | GIABv3 GC correction with the `gc` column clipped to [0.3, 0.7] | `FASTA`, `GC_CLIP_MIN/MAX`; 261011 GC page |
| 5 | `--diploid-parx-genome grch38` | `DIPLOID_PARX_GENOME`; 261011 edge and sex page |
| 6 | Inputs must come from the non-UMI pipeline: `sentieon-bwa` 5.1.0 with duplicates **marked**, `eggd_chr_prefix` 1.1.0, `eggd_cgp-cnvkit-coverage` 2.0.2 with bins BED v1.1, in the four DI-4114 run projects | `check_provenance` (`ALIGNER_APP`, `ALIGNER_INPUTS`, `CHR_PREFIX_APP`, `COVERAGE_APP`, `BINS_BED`, `RUN_PROJECTS`) |
| 7 | Two members disagree with their LIMS sex letter and are retained (a male with a tumour chrX gain; a probable 47,XXY) | `EXPECTED_N_SEX_DISAGREEMENTS`, `EXPECTED_SEX_DISAGREEMENTS_SHA256` (digest of the specimen numbers, which stay on DNAnexus) |

**Rebuild triggers.** Any change to the bins BED, FASTA, CNVkit image or app, the coverage app, the alignment or
demultiplexing (UMI handling, duplicate marking, read structure), the backbone/chemistry or the members requires a new
version of the file and a re-check; the pinned constants make the build fail on such a change.

## Usage — recreating the file

```bash
python3 build_atlas_cnvkit_static_pon.py
```

Build of record (`--mode app`, the default): runs the released `eggd_cgp-cnvkit-pon` 2.2.0 app on DNAnexus into a
fresh folder `project-JF3p1f04Yggp0vf3PYjy7Q0F:/pon_build/runs/<UTC timestamp>/` and waits (about 10 minutes,
about $0.04). Requires `dx` (authenticated DNAnexus CLI with access to the 004 project and billing to `org-emee_1`)
and Python 3 standard library only.

Independent rebuild without DNAnexus compute (`--mode local`): the same `cnvkit.py reference`, gc clipping and
`cnvkit.py sex` in the pinned image. Requires `docker` with `cgp-cnvkit:1.0.0` loaded from the production image tar
(001_Reference `file-J8j7Vyj45FG1BbK26JQgQY6q`; no registry), and the GIABv3 FASTA (`--fasta-dir`, else downloaded,
about 3 GB).

```bash
python3 build_atlas_cnvkit_static_pon.py --mode local --fasta-dir /path/to/GIABv3
```

Run in a fresh, empty `--workdir` (default `/tmp/build_atlas_cnvkit_static_pon`; the script refuses a non-empty one),
so every input is downloaded from scratch. The script then:

1. downloads the member table and checks its sha256 and contents (124 members, unique specimens and files, members
   per run);
2. downloads the 124 coverage files and checks each sha256;
3. checks the upstream provenance of every coverage file (decision 6);
4. builds the PoN (app job or local Docker);
5. checks the output: md5 `== EXPECTED_PON_MD5`, 14,393 bins, every `gc` in [0.3, 0.7], `reference_stats` (app mode:
   124 samples, GC applied), the member sex disagreements are the expected reviewed set (count and digest);
6. writes `<workdir>/build_report.json` and exits 0, printing `REPRODUCED: ...`. Any mismatch exits non-zero.

There is no separate test mode: a clean exit of this command is the reproducibility proof.

### Re-deriving the member table

```bash
python3 derive_pon_members.py
```

Downloads the evaluation sample manifest (pinned), re-checks every exclusion against the rules, selects the members,
looks up each member's non-UMI coverage file in 004 `/coverage_nonumi/<run>/`, records its sha256, and asserts the
result is byte-identical to the pinned member table (`EXPECTED_MEMBERS_SHA256`).

## Tests

```bash
pip install pytest
pytest test_build_atlas_cnvkit_static_pon.py -v
```

Unit tests cover the pure logic with synthetic fixtures (no real identifiers): member-table parsing and validation,
every link of the provenance chain and each way it can break, the `dx run` arguments, the gc clipping (byte-level
rule of the app), the output checks, the sex-disagreement rule, and the member-selection rules (exclusions, highest
depth, tie-break). DNAnexus, Docker and network steps are exercised by the from-scratch runs, not mocked.

### CI

`.github/workflows/ci.yml` runs on every push and pull request: a syntax check (`py_compile`) of both scripts and the
unit tests. It does not run a build (no `dx`/Docker in CI), so it cannot check the pinned checksum end to end; that is
the from-scratch run above, recorded on the Confluence page.

## Related

- Evaluation that chose the configuration: [261011 CNVkit PoN configuration for Atlas — evaluation design](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4836950081) and child pages.
- Validation of this PoN on the non-UMI data (DI-4118): logged jobs in 004 `/final_config_nonumi/`; driver and
  results in the evaluation's analysis repository (`analysis/pon_comparison/nonumi/`), linked from the Confluence page.
- Bins BED: [atlas_cnvkit_bins_for_cgp_hg38_v1.1.bed](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4805623944)
  ([eastgenomics/atlas_cnvkit_bins_for_cgp_hg38.bed](https://github.com/eastgenomics/atlas_cnvkit_bins_for_cgp_hg38.bed)).

## Files

| File | Purpose |
|---|---|
| `build_atlas_cnvkit_static_pon.py` | Build / reproduce the PoN; all pins |
| `derive_pon_members.py` | Re-derive the member table from the evaluation manifest and check it equals the pinned one |
| `test_build_atlas_cnvkit_static_pon.py` | Unit tests |
| `.github/workflows/ci.yml` | CI: syntax check and unit tests |
