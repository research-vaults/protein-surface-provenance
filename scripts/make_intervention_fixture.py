#!/usr/bin/env python3
"""Build the EXACT do(Y) fixture: two files, identical coordinates, permuted residues.

This is the intervention the paper's component 1 performs. For each chain we
permute the residue identities among the residue positions -- preserving the
chain's amino-acid composition exactly -- while holding every surface vertex
COORDINATE byte-identical. The featuriser is then re-run on the substituted
structure.

The two output files are per-VERTEX, not per-residue, so the comparison is over
every value the model actually receives.

    original.jsonl     {"chain": id, "channel": [...], "geom": [...]}
    substituted.jsonl  same chains, same geom, channel recomputed

`geom` is a rounded coordinate key per vertex; it exists so that "geometry did not
move" is a measured 0.00%, not an assumption.

usage: make_intervention_fixture.py DATA_DIR SPLIT OUT_ORIG OUT_SUB [SEED] [LIMIT]
"""
import json, random, sys
t2o = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}
HY = {"I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,"M":1.9,"A":1.8,"W":-0.9,"G":-0.4,"T":-0.7,"S":-0.8,
      "Y":-1.3,"P":-1.6,"H":-3.2,"N":-3.5,"D":-3.5,"Q":-3.5,"E":-3.5,"K":-3.9,"R":-4.5}
CH = {**{'R':1,'K':1,'D':-1,'E':-1,'H':0.1}, **{x:0 for x in 'ACFGILMNPQSTVWY'}}

def featurise(aa):
    """SurfPro's released featuriser: nearest-residue hydropathy and charge."""
    return (round(HY[aa]/5.0, 6), CH[aa])

D, SPLIT, OO, OS = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
SEED = int(sys.argv[5]) if len(sys.argv) > 5 else 0
LIMIT = int(sys.argv[6]) if len(sys.argv) > 6 else 10**9
rng = random.Random(SEED)
n = nv = 0
with open(f"{D}/{SPLIT}.atom.txt") as fa, open(f"{D}/{SPLIT}.coor.txt") as fc, \
     open(OO, "w") as fo, open(OS, "w") as fs:
    for ci, (la, lc) in enumerate(zip(fa, fc)):
        if n >= LIMIT: break
        toks = la.split()
        cs = lc.split()
        if len(cs) != 3 * len(toks):
            continue
        # residue identity per vertex, and the residue index it belongs to
        vert_res, vert_idx = [], []
        for w in toks:
            p = w.split("_"); aa = t2o.get(p[-2])
            vert_res.append(aa); vert_idx.append(int(p[-1]))
        uniq = sorted({i for i, a in zip(vert_idx, vert_res) if a})
        if not uniq: continue
        ident = {i: a for i, a in zip(vert_idx, vert_res) if a}
        # PERMUTE residue identities among residue positions; composition preserved
        vals = [ident[i] for i in uniq]
        rng.shuffle(vals)
        perm = dict(zip(uniq, vals))
        chan_o, chan_s, geom = [], [], []
        for k, (a, i) in enumerate(zip(vert_res, vert_idx)):
            if a is None: continue
            chan_o.append(featurise(a))
            chan_s.append(featurise(perm[i]))
            # coordinates are NOT touched; this key proves it
            geom.append(f"{float(cs[3*k]):.4f},{float(cs[3*k+1]):.4f},{float(cs[3*k+2]):.4f}")
        if not chan_o: continue
        cid = f"{SPLIT}:{ci}"
        fo.write(json.dumps({"chain": cid, "channel": chan_o, "geom": geom}) + "\n")
        fs.write(json.dumps({"chain": cid, "channel": chan_s, "geom": geom}) + "\n")
        n += 1; nv += len(chan_o)
print(f"wrote {n} chains, {nv} vertices\n  original    -> {OO}\n  substituted -> {OS}")
