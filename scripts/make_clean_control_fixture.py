#!/usr/bin/env python3
"""A DESIGN-VALID but TARGET-CORRELATED channel, and its intervention fixture.

The hard case for any provenance test is a channel that is genuinely computed
from geometry yet correlates with residue identity -- burial correlates with
hydrophobicity, so a burial-derived channel carries real information about the
sequence without ever reading it.

We build exactly that: per vertex, a quantisation of its distance to the chain
centroid and the local vertex density. No residue name is read. Under do(Y) the
channel cannot move, so this fixture's `substituted` file is identical to its
`original` by construction -- which is the point: an honest tool must return
0.00% here even though the channel is informative about the sequence.

Also reports the measured mutual information between channel class and residue,
so "target-correlated" is a number rather than a claim.
"""
import json, math, sys, collections
import numpy as np
t2o = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}
D, SPLIT, OO, OS = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
joint = collections.Counter()
n = nv = 0
with open(f"{D}/{SPLIT}.atom.txt") as fa, open(f"{D}/{SPLIT}.coor.txt") as fc, \
     open(OO, "w") as fo, open(OS, "w") as fs:
    for ci, (la, lc) in enumerate(zip(fa, fc)):
        toks = la.split()
        c = np.fromstring(lc, sep=" ", dtype=np.float64)
        if c.size != 3 * len(toks): continue
        P = c.reshape(-1, 3); cen = P.mean(0)
        d = np.linalg.norm(P - cen, axis=1)
        scale = d.mean() + 1e-9
        # local density: vertices within 0.15*scale, computed on a coarse grid
        g = np.floor(P / (0.15 * scale)).astype(np.int64)
        cnt = collections.Counter(map(tuple, g))
        dens = np.array([cnt[tuple(x)] for x in g], dtype=np.float64)
        shell = np.clip((d / scale * 4.0).astype(int), 0, 7)          # 8 radial bins
        dbin = np.clip((np.log1p(dens) / 1.2).astype(int), 0, 2)      # 3 density bins
        chan, geom = [], []
        for k, w in enumerate(toks):
            aa = t2o.get(w.split("_")[-2])
            if aa is None: continue
            cls = int(shell[k]) * 3 + int(dbin[k])
            chan.append(cls)
            geom.append(f"{P[k,0]:.4f},{P[k,1]:.4f},{P[k,2]:.4f}")
            joint[(cls, aa)] += 1
        if not chan: continue
        cid = f"{SPLIT}:{ci}"
        rec = json.dumps({"chain": cid, "channel": chan, "geom": geom})
        fo.write(rec + "\n"); fs.write(rec + "\n")   # identical: do(Y) cannot move it
        n += 1; nv += len(chan)
tot = sum(joint.values())
pc = collections.Counter(); pa = collections.Counter()
for (c_, a), k in joint.items(): pc[c_] += k; pa[a] += k
H = lambda dd: -sum((v/tot)*math.log2(v/tot) for v in dd.values() if v)
I = H(pc) + H(pa) - (-sum((k/tot)*math.log2(k/tot) for k in joint.values() if k))
print(f"wrote {n} chains, {nv} vertices")
print(f"channel classes: {len(pc)}")
print(f"I(channel ; residue) = {I:.4f} bits   H(residue) = {H(pa):.4f} bits"
      f"   -> {100*I/H(pa):.2f}% of residue entropy")
