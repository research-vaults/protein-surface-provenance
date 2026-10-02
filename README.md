# Protein-surface provenance diagnostics

Code and evidence for **“What Does a Protein Surface Add? Ablations Cannot
Tell You.”** The paper shows that an ordinary surface-versus-no-surface
ablation conflates four questions: how the channel was constructed, whether
fixed model weights use it, whether it adds conditional information beyond the
backbone, and whether it helps a specified finite learner.

**Paper:** [ICBINB-BIO workshop submission, 29 August
2026](paper/what-does-a-protein-surface-add-icbinb-bio-2026.pdf)  
**Code artifact:** 2 October 2026; BSD-3-Clause for project-authored code and
documentation.

The linked PDF is the exact submitted workshop version (SHA-256
`93b0c6624df2a5672c831447cfda0e0c1c32b8ea8c166bd481f69ec04277b08b`).
It is distinct from later internal research revisions. This public repository
is a research artifact, not an anonymous reviewer snapshot.

## Quick start

The safe first check uses only Python 3.8+ and the standard library. It makes
no network calls and normally finishes in seconds.

```bash
shasum -a 256 -c MANIFEST.sha256
python3 tests/test_provenance_probe.py
python3 tests/test_cache_paths.py
python3 results/identification_example/run.py
```

Expected output ends with `ALL 27 TESTS PASSED`, `ALL 9 CACHE-PATH TESTS
PASSED`, and a three-state example in which standalone accuracy is equal but
joint accuracy and conditional information differ. The manifest check covers
every released payload except the manifest itself.

## What is included

- `provenance_probe.py`: dependency-free command-line audit for controlled
  label substitution and an explicitly non-identifying observational screen.
- `examples/`: synthetic target-dependent and design-valid fixtures.
- `tests/`: offline CLI, failure-mode, and cache-path tests.
- `scripts/`: templates for fixed-coordinate interventions and
  backbone-derived controls using inputs supplied by the user.
- `reproduction/`: self-fetching scripts for the optional channel-capacity and
  frozen-model contrasts.
- `results/`: protocols and aggregate outputs for the exact identification,
  complementarity, and counterfactual-specificity analyses.
- `PROTOCOL.md`, `CHECKLIST.md`, and `REPORT_TEMPLATE.md`: reusable guidance
  for specifying, auditing, and reporting derived-input diagnostics.

Raw SurfPro/SurfDesign-derived arrays, checkpoints, and upstream repositories
are not redistributed because their repository licences did not clearly
authorize that release. See `THIRD_PARTY_NOTICES.md` and
`REPRODUCTION_LIMITS.md`.

## Minimal diagnostic examples

Target-dependent construction under a labels-only intervention:

```bash
python3 provenance_probe.py intervene \
  --original examples/target_original.jsonl \
  --substituted examples/target_substituted.jsonl
```

Known-clean construction:

```bash
python3 provenance_probe.py intervene \
  --original examples/clean_original.jsonl \
  --substituted examples/clean_substituted.jsonl
```

The first reports `TARGET_DEPENDENCE_DETECTED`; the second reports
`NO_DEPENDENCE_DETECTED`. The latter is intentionally not a certificate that
every possible substitution would leave the channel unchanged.

The observational screen is triage only:

```bash
python3 provenance_probe.py screen \
  --train examples/screen_demo.jsonl \
  --test examples/screen_demo.jsonl
```

It reports association and a train-fitted lookup score but issues no provenance
verdict.

## Result-to-command map

| Paper/evidence role | Command | What it establishes |
|---|---|---|
| Controlled provenance logic | `python3 tests/test_provenance_probe.py` | 27 synthetic assertions for intervention semantics and failure handling |
| Cache contract | `python3 tests/test_cache_paths.py` | Default, environment, and explicit-CLI path precedence without downloads |
| Finite identification example | `python3 results/identification_example/run.py` | Exact counterexample separating standalone from joint information |
| Hydropathy-channel capacity | `python3 reproduction/reproduce_channel_capacity.py` | Recomputes three held-out capacity quantities from public CATH data |
| Frozen-model contrast | `python3 reproduction/reproduce_frozen_contrast.py --smoke-test` | Small execution check; not the paper result |
| Full frozen-model contrast | `python3 reproduction/reproduce_frozen_contrast.py` | Recomputes the held-out ProteinMPNN contrast and bootstrap interval |

The complementarity and counterfactual-specificity directories contain
aggregate records and protocols, not redistributable raw arrays.

## Optional large reproductions

Install the optional dependencies in a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

Both reproduction scripts use a shared cache root. By default it is
`reproduction/.cache/`. Set `SURFACE_PROVENANCE_CACHE` to move it:

```bash
export SURFACE_PROVENANCE_CACHE=/path/with/adequate/space
python3 reproduction/reproduce_channel_capacity.py
python3 reproduction/reproduce_frozen_contrast.py --smoke-test
```

The capacity script also accepts `--data-dir /path/to/data`, which overrides
`SURFACE_PROVENANCE_CACHE/data` for that command. The frozen-contrast script
accepts `--cache-dir /path/to/cache`, which overrides the environment variable.
Inspect resolution without creating directories or downloading anything:

```bash
python3 reproduction/reproduce_channel_capacity.py --show-paths
python3 reproduction/reproduce_frozen_contrast.py --show-paths
```

Resource expectations:

| Command | Download/cache | Compute and status |
|---|---|---|
| Channel capacity | CATH JSONL plus split file, about 494 MiB total | Standard library; streaming CPU calculation. Full fresh replay: **NOT_RUN** for this release because the release host lacked disk space. |
| Frozen contrast, smoke | Same CATH data plus a pinned ProteinMPNN checkout and weights, roughly 180 MiB more | CPU; 4 training and 2 test chains; checks execution only. **NOT_RUN** for this release. |
| Frozen contrast, full | Reuses the same cache | Approximately 15–25 minutes on CPU in the original workflow; hardware-dependent. **NOT_RUN** for this release. |

Allow at least 1 GiB free for the downloaded cache, plus separate space for the
Python environment. Peak memory was not benchmarked for this release. Downloads
are checksum-verified and written through temporary files so an interrupted
transfer does not become a trusted input.

## Inputs and outputs

`intervene` expects paired JSONL files with identical chain ordering:

```json
{"chain": "id", "channel": [0.1, 0.2], "geom": ["g1", "g2"]}
```

The substituted file must recompute `channel` after changing labels while
holding preprocessing inputs, geometry, correspondence, and random state fixed.
The command prints movement rates and can write a machine-readable report with
`--json output.json`.

## Validation status and limits

Freshly verified on 2 October 2026:

- manifest integrity;
- all 27 offline provenance assertions;
- all 9 offline cache-path assertions;
- exact finite identification example;
- Python compilation and clean-clone execution.

These checks establish released software behavior, not full empirical
regeneration. The optional large replays were not completed on the release
host. Aggregate empirical records are accompanied by protocols and input
hashes, but this repository does not reproduce every training run in the paper.
See `REPRODUCTION_LIMITS.md` for the exact boundary.

## Citation

Author-bearing citation metadata has not yet been added to this public artifact.
Until a public workshop record is available, cite the linked version by its
title and version: *What Does a Protein Surface Add? Ablations Cannot Tell
You*, ICBINB: Failure Modes of AI in Biology workshop submission at NeurIPS
2026, submitted 29 August 2026.

## Licence

Project-authored code and documentation are BSD-3-Clause; see `LICENSE`.
The paper and third-party materials are not relicensed by that code licence.
