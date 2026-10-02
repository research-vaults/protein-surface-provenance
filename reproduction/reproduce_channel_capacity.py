#!/usr/bin/env python3
"""Reproduce the target capacity of the released hydropathy channel.

The script downloads the public CATH 4.2 inverse-folding benchmark files when
needed, verifies their SHA-256 hashes, and computes three held-out quantities:
Bayes-optimal recovery from the channel alone, the uniquely resolved residue
share, and the amino-acid frequency baseline.
"""
import argparse
import collections
import hashlib
import json
import os
import urllib.request

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
HYDROPATHY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
    "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9,
    "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5,
    "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
EXPECTED = {
    "chain_set.jsonl": "1944d21b975a11543b1b92fe906a0bb4a9a276a52ca3d493e5aadafd09f464c1",
    "chain_set_splits.json": "8e9a587a50c7f6c026e4ed00f6c1c30b106100f36f7a01de47542bdfc060adc2",
}
BASE_URL = "https://people.csail.mit.edu/ingraham/graph-protein-design/data/cath/"


def ensure_data(directory):
    os.makedirs(directory, exist_ok=True)
    for name, expected in EXPECTED.items():
        path = os.path.join(directory, name)
        if not os.path.exists(path):
            print(f"fetching {name}")
            partial = path + ".partial"
            try:
                urllib.request.urlretrieve(BASE_URL + name, partial)
                os.replace(partial, path)
            finally:
                if os.path.exists(partial):
                    os.remove(partial)
        actual = hashlib.sha256(open(path, "rb").read()).hexdigest()
        if actual != expected:
            raise SystemExit(f"hash mismatch for {name}: {actual}")
        print(f"verified {name}")


def main():
    parser = argparse.ArgumentParser()
    default = os.path.join(os.path.dirname(__file__), ".cache", "data")
    parser.add_argument("--data-dir", default=default)
    args = parser.parse_args()
    ensure_data(args.data_dir)

    split = json.load(open(os.path.join(args.data_dir, "chain_set_splits.json")))
    test = set(split["test"])
    counts = collections.Counter()
    with open(os.path.join(args.data_dir, "chain_set.jsonl")) as handle:
        for line in handle:
            row = json.loads(line)
            if row["name"] in test and len(row["seq"]) <= 500:
                counts.update(a for a in row["seq"] if a in ALPHABET)

    groups = collections.defaultdict(list)
    for residue in ALPHABET:
        groups[HYDROPATHY[residue]].append(residue)
    total = sum(counts.values())
    bayes = sum(max(counts[a] for a in group) for group in groups.values()) / total
    unique = sum(counts[group[0]] for group in groups.values() if len(group) == 1) / total
    floor = max(counts[a] for a in ALPHABET) / total
    print(f"residues: {total:,}")
    print(f"channel values: {len(groups)}")
    print(f"Bayes-optimal channel-only recovery: {100 * bayes:.2f}%")
    print(f"uniquely resolved residue share: {100 * unique:.2f}%")
    print(f"frequency baseline: {100 * floor:.2f}%")


if __name__ == "__main__":
    main()
