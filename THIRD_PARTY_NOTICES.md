# Third-party notices and data policy

This repository contains project-authored diagnostic and reproduction code. It
does not bundle third-party repositories, model weights, benchmark archives, or
raw outputs derived from repositories without a detected licence.

## ProteinMPNN

`reproduction/reproduce_frozen_contrast.py` can clone the public ProteinMPNN
repository at a pinned revision. ProteinMPNN is not vendored here. Its upstream
repository includes an MIT licence:

- <https://github.com/dauparas/ProteinMPNN>

## SurfPro and SurfDesign

The research evaluates the public SurfPro and SurfDesign implementations:

- <https://github.com/JocelynSong/SurfPro>
- <https://github.com/smiles724/SurfDesign>

No machine-readable repository licence was detected for either repository when
this release was prepared on 2 October 2026. Their code, weights, data, and raw
derived fixtures or output arrays are therefore not redistributed here. The
project-authored scripts describe the transformations needed after users obtain
upstream inputs under applicable terms.

## CATH and benchmark inputs

Benchmark data are not bundled. Reproduction scripts fetch public copies and
verify expected hashes. Users must comply with the upstream source's terms and
should not treat this project's BSD licence as applying to downloaded data.

## Scientific attribution

Method and dataset citations appear in the associated paper and relevant
script documentation. Removing third-party files from this release does not
remove or replace those scientific attributions.
