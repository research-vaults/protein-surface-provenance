"""Exact enumeration; source interventions are software evaluations off support."""
import json, math
from fractions import Fraction
from collections import Counter, defaultdict
from pathlib import Path
P=Path(__file__).resolve().parent
rows=[]
for name, fn in [('valid',lambda b,y:b),('redundant_target',lambda b,y:int(y==0)),('complementary_target',lambda b,y:int(y==1))]:
 cells=[(y,int(y==0),fn(int(y==0),y)) for y in range(3)]
 def acc(indices):
  d=defaultdict(Counter)
  for cell in cells:d[tuple(cell[i] for i in indices)][cell[0]]+=1
  return Fraction(sum(max(v.values()) for v in d.values()),3)
 def entropy(indices):
  counts=Counter(tuple(cell[i] for i in indices) for cell in cells)
  return -sum(n/3*math.log2(n/3) for n in counts.values())
 cmi=entropy([0,1])+entropy([1,2])-entropy([1])-entropy([0,1,2])
 changed=sum(fn(b,y)!=fn(b,yp) for y,b,c in cells for yp in range(3))
 rows.append(dict(channel=name,states=cells,histogram=dict(Counter(c for y,b,c in cells)),cardinality=2,standalone_accuracy=str(acc([2])),joint_accuracy=str(acc([1,2])),conditional_mi_bits=cmi,changed_of_nine=changed))
assert [r['standalone_accuracy'] for r in rows]==['2/3']*3
assert [r['joint_accuracy'] for r in rows]==['2/3','2/3','1']
assert all(abs(a-b)<1e-12 for a,b in zip([r['conditional_mi_bits'] for r in rows],[0,0,2/3]))
assert [r['changed_of_nine'] for r in rows]==[0,4,4]
(P/'results.json').write_text(json.dumps({'population':'exact uniform three-state construction','rows':rows,'checks_passed':4},indent=2)+'\n')
print(json.dumps(rows,indent=2))
