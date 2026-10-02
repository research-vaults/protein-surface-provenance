# Protein-surface provenance diagnostics

Public reproducibility release for **“What Does a Protein Surface Add? Ablations Cannot Tell You.”**

The paper separates four questions that a surface-versus-no-surface ablation
combines: how a channel was constructed, whether fixed model weights use it,
whether it adds conditional information beyond the backbone, and whether it
helps a specified finite learner. This repository provides the lightweight
source-provenance diagnostic, synthetic examples, aggregate result records,
and self-fetching reproduction scripts that can be distributed safely.

## What is included

- `provenance_probe.py`: dependency-free command-line audit for controlled
  label substitution and an explicitly non-identifying observational screen.
- `examples/`: synthetic target-dependent and design-valid fixtures. These
  contain no protein sequences, structures, personal data, or third-party
  model outputs.
- `tests/`: 27 offline assertions covering the documented CLI and failure
  modes.
- `scripts/`: templates for constructing fixed-coordinate interventions and
  backbone-derived controls from inputs supplied by the user.
- `reproduction/`: scripts that fetch and verify public inputs at runtime for
  the channel-capacity and frozen-model contrasts.
- `results/`: protocols and aggregate outputs for the exact identification
  example, complementarity analysis, and counterfactual-specificity analysis.
- `PROTOCOL.md`, `CHECKLIST.md`, and `REPORT_TEMPLATE.md`: reporting guidance.

## Deliberately excluded

The public release does not contain author identities or affiliations, internal
reviews, planning journals, roadmaps, acceptance correspondence, private data,
machine-specific paths, checkpoints, third-party repositories, venue-specific
drafts, or raw SurfPro/SurfDesign-derived fixtures and model-output arrays.
SurfPro and SurfDesign did not expose a detected repository licence when this
release was prepared, so their source, weights, data, and derived raw outputs
are not redistributed here. See `THIRD_PARTY_NOTICES.md` and
`REPRODUCTION_LIMITS.md`.

The manuscript is also deliberately omitted. The code release is public and
linked to the paper title; it is **not** an anonymous conference-review
artifact. Authoritative conference and archival manuscripts remain separate.

## Requirements

The diagnostic and tests require Python 3.8+ and the standard library only.
Optional reproduction scripts use the packages pinned in `requirements.txt`.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt  # optional, for full reproductions
```

## Quick verification

From the repository root:

```bash
shasum -a 256 -c MANIFEST.sha256
python3 tests/test_provenance_probe.py
python3 results/identification_example/run.py
```

Expected: every manifest entry is `OK`, all 27 assertions pass, and the exact
three-state example reports equal standalone accuracy but different joint
accuracy and conditional information.

## Minimal examples

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
`NO_DEPENDENCE_DETECTED`, which is intentionally not a certificate that every
possible substitution would leave the channel unchanged.

The observational screen is triage only:

```bash
python3 provenance_probe.py screen \
  --train examples/screen_demo.jsonl \
  --test examples/screen_demo.jsonl
```

It reports association and a train-fitted lookup score but issues no provenance
verdict.

## Reproducing larger quantities

Channel-only capacity:

```bash
python3 reproduction/reproduce_channel_capacity.py
```

Frozen ProteinMPNN contrast:

```bash
python3 reproduction/reproduce_frozen_contrast.py --smoke-test
python3 reproduction/reproduce_frozen_contrast.py
```

These commands download public inputs into an ignored cache, verify pinned
hashes or commits, and do not require files from the private research archive.
Set `SURFACE_PROVENANCE_CACHE=/path/to/cache` to choose the cache location.
Users remain responsible for the upstream data and model terms.

## Inputs and outputs

`intervene` expects paired JSONL files with identical chain ordering:

```json
{"chain": "id", "channel": [0.1, 0.2], "geom": ["g1", "g2"]}
```

The substituted file must recompute `channel` after changing labels while
holding every other preprocessing input, geometry, correspondence, and random
state fixed. The command prints movement rates and can write a machine-readable
report with `--json output.json`.

## Reproducibility boundary

Offline tests establish software behavior on synthetic cases. The exact finite
identification example is fully replayable. Aggregate empirical records are
published with protocols and input hashes, but raw arrays that could not be
confidently redistributed are not public. The self-fetching scripts cover two
larger quantities; they do not reproduce every training run in the paper.

## Licence

Project-authored code and documentation are BSD-3-Clause; see `LICENSE`.
Third-party materials are not relicensed by this repository.
