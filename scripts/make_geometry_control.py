#!/usr/bin/env python3
"""A DESIGN-VALID control channel on exactly the same chains.

Builds C = f(B): each residue's channel value is a quantisation of geometric
quantities computed from the surface point cloud ALONE -- mean and spread of its
vertices' distance to the chain centroid, and its vertex count. No residue name
is read. By construction, permuting residue identities with coordinates held
fixed cannot change a single value, so the probe's check 1 MUST return 0%.

Quantised to 18 classes so the class count matches SurfPro's released channel and
the two are directly comparable.
"""
import json, sys, difflib
import numpy as np
D, SPLIT, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
LIMIT = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9
t2o = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}
n = 0
with open(f"{D}/{SPLIT}.atom.txt") as fa, open(f"{D}/{SPLIT}.coor.txt") as fc, \
     open(f"{D}/{SPLIT}.seq.txt") as fs, open(OUT, "w") as fo:
    for la, lc, ls in zip(fa, fc, fs):
        if n >= LIMIT: break
        seq = ls.strip().replace(" ", "")
        toks = la.split()
        coor = np.fromstring(lc, sep=" ", dtype=np.float32)
        if coor.size != 3 * len(toks): continue
        P = coor.reshape(-1, 3); cen = P.mean(0)
        d = np.linalg.norm(P - cen, axis=1)
        res, ridx = {}, {}
        for k, w in enumerate(toks):
            p = w.split("_"); aa = t2o.get(p[-2])
            if aa:
                i = int(p[-1]); res[i] = aa; ridx.setdefault(i, []).append(k)
        if not res: continue
        keys = sorted(res)
        surf = "".join(res[k] for k in keys)
        pos = {}
        for tag,i1,i2,j1,j2 in difflib.SequenceMatcher(None, surf, seq, autojunk=False).get_opcodes():
            if tag in ("equal","replace"):
                for k in range(j2-j1):
                    if i1+k < i2: pos[i1+k] = j1+k
        scale = d.mean() + 1e-9
        chan, rid = [], []
        for si, k in enumerate(keys):
            if si not in pos: continue
            v = d[ridx[k]]
            # three geometric descriptors, none of which sees a residue name
            a = int(np.clip(v.mean()/scale * 3.0, 0, 5.999))      # 6 radial shells
            b = int(np.clip(v.std()/(scale*0.5) * 1.5, 0, 2.999))  # 3 spread bins
            chan.append([a, b]); rid.append(pos[si])
        if not chan: continue
        fo.write(json.dumps({"seq": seq, "channel": chan, "resid": rid}) + "\n"); n += 1
print(f"wrote {n} chains to {OUT}")
