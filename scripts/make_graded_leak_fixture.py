#!/usr/bin/env python3
"""Instrument calibration: how much of a PARTIAL leak does `intervene` detect?

A real featuriser is rarely all-or-nothing. This builds a family of channels in
which a fixed fraction `p` of surface vertices take the target-derived lookup
value and the remaining `1-p` take a value computed from geometry alone. The
lookup/geometry assignment is drawn ONCE per vertex and is identical in the
original and substituted files, so the only thing differing between them is
residue identity.

Expected: the measured substitution rate scales linearly in `p`, reaching the
full-leak rate at p = 1 and exactly zero at p = 0.

usage: make_graded_leak_fixture.py DATA_DIR SPLIT P OUT_ORIG OUT_SUB [SEED] [LIMIT]
"""
import json, random, sys
import numpy as np
t2o = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}
HY = {"I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,"M":1.9,"A":1.8,"W":-0.9,"G":-0.4,"T":-0.7,"S":-0.8,
      "Y":-1.3,"P":-1.6,"H":-3.2,"N":-3.5,"D":-3.5,"Q":-3.5,"E":-3.5,"K":-3.9,"R":-4.5}
CH = {**{'R':1,'K':1,'D':-1,'E':-1,'H':0.1}, **{x:0 for x in 'ACFGILMNPQSTVWY'}}
D, SPLIT, P, OO, OS = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4], sys.argv[5]
SEED = int(sys.argv[6]) if len(sys.argv) > 6 else 0
LIMIT = int(sys.argv[7]) if len(sys.argv) > 7 else 100
rng = random.Random(SEED); nprng = np.random.default_rng(SEED)
n = nv = 0
with open(f"{D}/{SPLIT}.atom.txt") as fa, open(f"{D}/{SPLIT}.coor.txt") as fc, \
     open(OO, "w") as fo, open(OS, "w") as fs:
    for ci, (la, lc) in enumerate(zip(fa, fc)):
        if n >= LIMIT: break
        toks = la.split()
        c = np.fromstring(lc, sep=" ", dtype=np.float64)
        if c.size != 3 * len(toks): continue
        Pt = c.reshape(-1, 3); cen = Pt.mean(0)
        d = np.linalg.norm(Pt - cen, axis=1); scale = d.mean() + 1e-9
        geo = np.clip((d / scale * 6.0).astype(int), 0, 17)
        res, idx = [], []
        for k, w in enumerate(toks):
            p_ = w.split("_"); res.append(t2o.get(p_[-2])); idx.append(int(p_[-1]))
        uniq = sorted({i for i, a in zip(idx, res) if a})
        if not uniq: continue
        ident = {i: a for i, a in zip(idx, res) if a}
        vals = [ident[i] for i in uniq]; rng.shuffle(vals)
        perm = dict(zip(uniq, vals))
        use_lookup = nprng.random(len(toks)) < P
        co, cs_, gm = [], [], []
        for k, (a, i) in enumerate(zip(res, idx)):
            if a is None: continue
            if use_lookup[k]:
                co.append([round(HY[a]/5.0, 6), CH[a]])
                cs_.append([round(HY[perm[i]]/5.0, 6), CH[perm[i]]])
            else:
                g = int(geo[k]); co.append(g); cs_.append(g)
            gm.append(f"{Pt[k,0]:.4f},{Pt[k,1]:.4f},{Pt[k,2]:.4f}")
        if not co: continue
        cid = f"{SPLIT}:{ci}"
        fo.write(json.dumps({"chain": cid, "channel": co, "geom": gm}) + "\n")
        fs.write(json.dumps({"chain": cid, "channel": cs_, "geom": gm}) + "\n")
        n += 1; nv += len(co)
print(f"p={P}  {n} chains, {nv} vertices")
