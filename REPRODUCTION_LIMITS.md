# Reproduction limits

This public repository is an allowlisted export, not the complete private
research archive.

## Fully replayable without network access

- the 27 diagnostic regression assertions;
- the 9 cache-path precedence assertions;
- the synthetic intervention examples; and
- the exact three-state identification example.

## Replayable with public downloads

- the channel-capacity calculation; and
- the frozen ProteinMPNN contrast.

The scripts pin or hash-check their public inputs. Their successful execution
still depends on upstream availability, platform-compatible dependencies, and
the upstream terms of use.

The shared cache root is `reproduction/.cache/` unless
`SURFACE_PROVENANCE_CACHE` is set. For channel capacity, `--data-dir` overrides
the resulting `data/` directory. For the frozen contrast, `--cache-dir`
overrides the cache root. Both scripts support `--show-paths`, which resolves
this precedence without creating directories or starting downloads.

The public CATH files occupy about 494 MiB. The frozen contrast additionally
uses a pinned ProteinMPNN checkout and weights of roughly 180 MiB. Allow at
least 1 GiB for the cache plus separate space for the Python environment. The
full frozen CPU workflow was previously estimated at 15--25 minutes; runtime
is hardware-dependent, and peak memory was not benchmarked for this release.

These two large replays were **not run to completion while preparing the public
release**. An attempted CATH download stopped when the release host ran out of
disk space. The incomplete ignored file was removed and never committed. The
offline tests and exact finite example passed; those results must not be read
as a full empirical replay.

## Aggregate records only

The complementarity and counterfactual-specificity directories contain their
protocols, input hashes, and aggregate machine-readable results. Raw retained
arrays are excluded because they derive from upstream SurfPro/SurfDesign assets
whose repositories did not expose a detected licence. The published hashes
bind the aggregate analyses to the retained private inputs but do not let a
third party reconstruct those inputs.

## Not covered

This release does not replay historic inverse-folding training, every model and
dataset cell, or the complete manuscript build. It does not claim exhaustive
privacy, security, portability, or scientific replication beyond the checks
listed above.
