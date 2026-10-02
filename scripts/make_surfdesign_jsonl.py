#!/usr/bin/env python3
"""Build provenance_probe input for SurfDesign's surface-chemistry channel.

SurfDesign assigns each surface vertex the Kyte-Doolittle hydropathy of the
NEAREST RESIDUE. That is a 20-row lookup collapsing to 17 classes, since D, E, N
and Q all carry -3.5. We emit one (value, residue index) pair per residue, i.e.
the full-coverage idealisation used by this project's EXP1 Bayes bound, so the
probe's check 2 is directly comparable to that independently computed number.

usage: make_surfdesign_jsonl.py chain_set.jsonl chain_set_splits.json SPLIT OUT [LIMIT]
"""
import json, sys
CS, SPL, SPLIT, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
LIMIT = int(sys.argv[5]) if len(sys.argv) > 5 else 10**9
KD = {"I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,"M":1.9,"A":1.8,"W":-0.9,"G":-0.4,"T":-0.7,
      "S":-0.8,"Y":-1.3,"P":-1.6,"H":-3.2,"N":-3.5,"D":-3.5,"Q":-3.5,"E":-3.5,"K":-3.9,"R":-4.5}
want = set(json.load(open(SPL))[SPLIT])
n = 0
with open(CS) as f, open(OUT, "w") as fo:
    for line in f:
        if n >= LIMIT: break
        r = json.loads(line)
        if r.get("name") not in want: continue
        seq = r["seq"]
        chan, rid = [], []
        for i, a in enumerate(seq):
            if a in KD:
                chan.append(KD[a]); rid.append(i)
        if not chan: continue
        fo.write(json.dumps({"seq": seq, "channel": chan, "resid": rid}) + "\n"); n += 1
print(f"wrote {n} chains to {OUT}")
