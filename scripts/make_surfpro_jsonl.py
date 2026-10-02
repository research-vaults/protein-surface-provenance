#!/usr/bin/env python3
"""Convert SurfPro's released surface files into provenance_probe's JSONL format.

Each surface vertex names its nearest residue (e.g. `CG1_VAL_2`); SurfPro turns
that name into a 2-vector via a 20-row hydropathy/charge table. We emit, per
chain, one (channel value, residue index) pair per NAMED RESIDUE -- the vertex
redundancy is irrelevant to the probe, which works on the value/residue mapping.
"""
import json, sys, collections
D = sys.argv[1]; SPLIT = sys.argv[2]; OUT = sys.argv[3]
LIMIT = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9
t2o = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}
HY = {"I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,"M":1.9,"A":1.8,"W":-0.9,"G":-0.4,"T":-0.7,"S":-0.8,
      "Y":-1.3,"P":-1.6,"H":-3.2,"N":-3.5,"D":-3.5,"Q":-3.5,"E":-3.5,"K":-3.9,"R":-4.5}
CH = {**{'R':1,'K':1,'D':-1,'E':-1,'H':0.1}, **{x:0 for x in 'ACFGILMNPQSTVWY'}}
import difflib
n=0
with open(f"{D}/{SPLIT}.atom.txt") as fa, open(f"{D}/{SPLIT}.seq.txt") as fs, open(OUT,"w") as fo:
    for la, ls in zip(fa, fs):
        if n >= LIMIT: break
        seq = ls.strip().replace(" ", "")
        res = {}
        for w in la.split():
            p = w.split("_"); aa = t2o.get(p[-2])
            if aa: res[int(p[-1])] = aa
        if not res: continue
        keys = sorted(res)
        surf = "".join(res[k] for k in keys)
        # map surface order -> sequence position by the same alignment used throughout
        pos = {}
        for tag,i1,i2,j1,j2 in difflib.SequenceMatcher(None, surf, seq, autojunk=False).get_opcodes():
            if tag in ("equal","replace"):
                for k in range(j2-j1):
                    if i1+k < i2: pos[i1+k] = j1+k
        chan, rid = [], []
        for si,k in enumerate(keys):
            if si in pos:
                a = res[k]
                chan.append([round(HY[a]/5.0,6), CH[a]]); rid.append(pos[si])
        if not chan: continue
        fo.write(json.dumps({"seq":seq,"channel":chan,"resid":rid})+"\n"); n+=1
print(f"wrote {n} chains to {OUT}")
