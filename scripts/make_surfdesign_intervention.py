#!/usr/bin/env python3
"""Exact do(Y) fixture for SurfDesign's surface-chemistry channel.

SurfDesign assigns each surface vertex the Kyte-Doolittle hydropathy of the
NEAREST RESIDUE, a 20-row lookup collapsing to 17 classes ({D,E,N,Q} share -3.5).

We permute residue identities among residue positions -- preserving each chain's
amino-acid composition exactly -- and hold the BACKBONE COORDINATES byte-identical.
`geom` carries the residue's N/CA/C/O coordinates, so geometric invariance is
measured rather than asserted.

usage: make_surfdesign_intervention.py chain_set.jsonl splits.json SPLIT OO OS [SEED] [LIMIT]
"""
import json, random, sys
CS, SPL, SPLIT, OO, OS = sys.argv[1:6]
SEED = int(sys.argv[6]) if len(sys.argv) > 6 else 0
LIMIT = int(sys.argv[7]) if len(sys.argv) > 7 else 10**9
KD = {"I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,"M":1.9,"A":1.8,"W":-0.9,"G":-0.4,"T":-0.7,
      "S":-0.8,"Y":-1.3,"P":-1.6,"H":-3.2,"N":-3.5,"D":-3.5,"Q":-3.5,"E":-3.5,"K":-3.9,"R":-4.5}
want = set(json.load(open(SPL))[SPLIT])
rng = random.Random(SEED)
n = nv = 0
with open(CS) as f, open(OO,"w") as fo, open(OS,"w") as fs:
    for line in f:
        if n >= LIMIT: break
        r = json.loads(line)
        if r.get("name") not in want: continue
        seq = r["seq"]; co = r.get("coords") or {}
        idx = [i for i,a in enumerate(seq) if a in KD]
        if not idx: continue
        perm = [seq[i] for i in idx]; rng.shuffle(perm)
        chan_o = [KD[seq[i]] for i in idx]
        chan_s = [KD[a] for a in perm]
        geom = []
        for i in idx:
            key = []
            for at in ("N","CA","C","O"):
                v = co.get(at)
                if v is not None and i < len(v) and v[i] is not None:
                    key.append(",".join(f"{float(x):.3f}" for x in v[i]))
            geom.append("|".join(key) if key else f"res{i}")
        cid = r["name"]
        fo.write(json.dumps({"chain":cid,"channel":chan_o,"geom":geom})+"\n")
        fs.write(json.dumps({"chain":cid,"channel":chan_s,"geom":geom})+"\n")
        n += 1; nv += len(idx)
print(f"wrote {n} chains, {nv} residues")
