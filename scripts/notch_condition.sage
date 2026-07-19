# notch_condition.sage - EXACT (a,b,c) condition for the even-m q0=q1 "notch" in SU(4) GHU
#   Author: Carles Marin  <karlesmarin@gmail.com>   (with Claude, Anthropic, as assistant)
# Claim under test (derived, not assumed):
#   Let D(t) = sum_weights t^(n2-n4) (-1)^(n3) = s_lambda(1, t, -1, 1/t).
#   Since {1,-t,-1,-1/t} = -{1,t,-1,1/t} as multisets, D(-t) = (-1)^(a+2b+3c) D(t).
#   Hence the even-m part of D vanishes  <=>  a+2b+3c odd   (i.e. Gate 1, the Z4 center charge),
#   unless D == 0 identically (degenerate case -- flagged separately).
# Notch := for every EVEN m, mult(m,q=0) == mult(m,q=1).
import json
A3 = WeylCharacterRing("A3", style="coroots")

def analyze(a,b,c):
    L=(a,b,c); boxes=a+2*b+3*c; h={}
    for w,mult in A3(*L).weight_multiplicities().items():
        n=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        key=(int(n[1]-n[3]), int(n[2]%2)); h[key]=h.get(key,0)+int(mult)
    ms=sorted(set(k[0] for k in h))
    sgn={m: h.get((m,0),0)-h.get((m,1),0) for m in ms}
    notch      = all(sgn[m]==0 for m in ms if m%2==0)
    notch_odd  = all(sgn[m]==0 for m in ms if m%2==1)
    Dzero      = all(sgn[m]==0 for m in ms)
    return dict(dim=int(A3(*L).degree()), boxes=boxes, notch=notch,
                notch_odd=notch_odd, Dzero=Dzero,
                worst=max([abs(sgn[m]) for m in ms if m%2==0]+[0]))

BOX=6; DIMCAP=2000
rows=[]; mism=[]
for a in range(BOX+1):
  for b in range(BOX+1):
    for c in range(BOX+1):
      if (a,b,c)==(0,0,0): continue
      d=A3(a,b,c).degree()
      if d>DIMCAP: continue
      r=analyze(a,b,c); r["abc"]=(a,b,c)
      pred = (r["boxes"]%2==1)
      r["pred"]=pred; r["match"]=(pred==r["notch"])
      rows.append(r)
      if not r["match"]: mism.append(r)

print("reps tested: %d (dim<=%d, a,b,c<=%d)" % (len(rows),DIMCAP,BOX))
print("MISMATCHES (notch != [a+2b+3c odd]): %d" % len(mism))
for r in mism[:40]: print("   ", r)
print("D==0 identically:", [r["abc"] for r in rows if r["Dzero"]])
print("odd-m also paired (beyond claim):", [r["abc"] for r in rows if r["notch_odd"]][:20])
# control: the 9 catalog reps
CAT={"35":(4,0,0),"60":(0,2,1),"84":(0,1,3),"140a":(1,1,2),"140b":(0,3,1),
     "216":(0,4,1),"224":(1,2,2),"280":(0,2,3),"360":(1,3,2)}
print("\n--- catalog control ---")
for k,v in CAT.items():
    r=analyze(*v); print("  %-5s %-9s boxes=%d parity=%s notch=%s worst_even_imbalance=%d" %
                         (k,str(v),r["boxes"],"odd" if r["boxes"]%2 else "even",r["notch"],r["worst"]))
json.dump([{kk:(list(vv) if isinstance(vv,tuple) else vv) for kk,vv in r.items()} for r in rows],
          open("notch_sweep.json","w"))
