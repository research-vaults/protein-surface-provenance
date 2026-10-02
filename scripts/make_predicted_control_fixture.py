#!/usr/bin/env python3
"""The HARD case: a channel that PREDICTS residue class from geometry alone.

This is the realistic clean-but-target-correlated control. Per residue we fit a
classifier on the training split to predict its 18-class chemistry code from
local surface geometry only -- sorted k-NN distances, local PCA shape, density,
distance to centroid -- and use the PREDICTION as the channel. It carries real
information about the sequence (that is the whole point of building it) and it is
`C = f(B)` by construction: no residue name is read at inference.

Under do(Y) with coordinates held fixed the prediction cannot change, so the
`substituted` file is identical to `original`. An honest provenance test must
return 0.00% here despite the channel being informative.
"""
import json, math, sys, collections
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.ensemble import HistGradientBoostingClassifier
D, OO, OS = sys.argv[1], sys.argv[2], sys.argv[3]
NTRAIN = int(sys.argv[4]) if len(sys.argv) > 4 else 600
t2o = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}
HY = {"I":4.5,"V":4.2,"L":3.8,"F":2.8,"C":2.5,"M":1.9,"A":1.8,"W":-0.9,"G":-0.4,"T":-0.7,"S":-0.8,
      "Y":-1.3,"P":-1.6,"H":-3.2,"N":-3.5,"D":-3.5,"Q":-3.5,"E":-3.5,"K":-3.9,"R":-4.5}
CH = {**{'R':1,'K':1,'D':-1,'E':-1,'H':0.1}, **{x:0 for x in 'ACFGILMNPQSTVWY'}}
A='ACDEFGHIKLMNPQRSTVWY'; CODE={a:(round(HY[a]/5,6),CH[a]) for a in A}
CLS=sorted({CODE[a] for a in A}); CI={c:i for i,c in enumerate(CLS)}
def feats(P):
    n=len(P); k=min(16,n-1)
    nn=NearestNeighbors(n_neighbors=k+1).fit(P); d,idx=nn.kneighbors(P)
    d=d[:,1:]; idx=idx[:,1:]; cen=P.mean(0)
    N=P[idx]-P[:,None,:]; C=np.einsum("nki,nkj->nij",N,N)/k
    ev=np.linalg.eigvalsh(C)[:,::-1]; s=ev.sum(1,keepdims=True)+1e-12; ev=ev/s
    rcn=np.linalg.norm(P-cen,axis=1)
    return np.concatenate([d[:,:8],d.mean(1,keepdims=True),d.std(1,keepdims=True),ev,
                           (ev[:,0]-ev[:,1])[:,None],np.log(s+1e-12),
                           (rcn/(rcn.mean()+1e-9))[:,None]],axis=1)
def chain(la,lc):
    toks=la.split(); c=np.fromstring(lc,sep=" ",dtype=np.float32)
    if c.size!=3*len(toks) or len(toks)<60: return None
    P=c.reshape(-1,3); F=feats(P)
    res={}; grp=collections.defaultdict(list)
    for k,w in enumerate(toks):
        p=w.split("_"); aa=t2o.get(p[-2])
        if aa: i=int(p[-1]); res[i]=aa; grp[i].append(k)
    if not res: return None
    keys=sorted(res); X=[]; Y=[]; VIDX=[]
    for i in keys:
        f=F[grp[i]]
        X.append(np.concatenate([f.mean(0),f.std(0),[len(grp[i])]]))
        Y.append(CI[CODE[res[i]]]); VIDX.append(grp[i])
    return np.array(X,dtype=np.float32), np.array(Y), VIDX, P
XS,YS=[],[]
with open(f"{D}/train.atom.txt") as fa, open(f"{D}/train.coor.txt") as fc:
    for i,(la,lc) in enumerate(zip(fa,fc)):
        if i>=NTRAIN: break
        r=chain(la,lc)
        if r: XS.append(r[0]); YS.append(r[1])
Xtr=np.concatenate(XS); Ytr=np.concatenate(YS)
print(f"train rows {Xtr.shape}",flush=True)
clf=HistGradientBoostingClassifier(max_iter=150,early_stopping=True,random_state=0).fit(Xtr,Ytr)
print("fitted",flush=True)
joint=collections.Counter(); n=nv=0; acc=[0,0]
with open(f"{D}/test.atom.txt") as fa, open(f"{D}/test.coor.txt") as fc, \
     open(OO,"w") as fo, open(OS,"w") as fs:
    for ci,(la,lc) in enumerate(zip(fa,fc)):
        r=chain(la,lc)
        if not r: continue
        X,Ytrue,VIDX,P=r; pred=clf.predict(X)
        acc[0]+=int((pred==Ytrue).sum()); acc[1]+=len(Ytrue)
        chan=[]; geom=[]
        for cls,vids in zip(pred,VIDX):
            for k in vids:
                chan.append(int(cls)); geom.append(f"{P[k,0]:.4f},{P[k,1]:.4f},{P[k,2]:.4f}")
        for cls,y,vids in zip(pred,Ytrue,VIDX):
            joint[(int(cls),int(y))]+=len(vids)
        cid=f"test:{ci}"; rec=json.dumps({"chain":cid,"channel":chan,"geom":geom})
        fo.write(rec+"\n"); fs.write(rec+"\n"); n+=1; nv+=len(chan)
tot=sum(joint.values()); pc=collections.Counter(); pa=collections.Counter()
for (c_,a),k in joint.items(): pc[c_]+=k; pa[a]+=k
H=lambda d:-sum((v/tot)*math.log2(v/tot) for v in d.values() if v)
I=H(pc)+H(pa)-(-sum((k/tot)*math.log2(k/tot) for k in joint.values() if k))
print(f"wrote {n} chains, {nv} vertices")
print(f"class accuracy from geometry: {100*acc[0]/acc[1]:.2f}%")
print(f"I(channel ; true class) = {I:.4f} bits  H(class) = {H(pa):.4f}  -> {100*I/H(pa):.2f}% of class entropy")
