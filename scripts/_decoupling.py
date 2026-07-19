# _decoupling.py - does the even-m notch make the FERMION contribution to the alpha2=1/2
#   curvature vanish identically, leaving only the gauge sector?
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
import json, math
H=json.load(open("hists.json"))
LAB={"35":(4,0,0),"60":(0,2,1),"84":(0,1,3),"140a":(1,1,2),"140b":(0,3,1),
     "216":(0,4,1),"224":(1,2,2),"280":(0,2,3),"360":(1,3,2)}
GAUGE=[(1,0,2),(-1,0,2),(1,1,2),(-1,1,2),(2,0,1),(0,0,1),(-2,0,1)]
def parse(h):
    d={}
    for k,v in h.items():
        m,q=map(int,k.split(',')); d[(m,q)]=d.get((m,q),0)+int(v)
    return d
def D_struct(d):
    """the higgs_vacii structure sum, per sector, in units of pi^2"""
    tot=0.0
    for (m,q),mu in d.items():
        if m%2: continue
        r=m%4
        s=+1 if ((r==2 and q==0) or (r==0 and q==1)) else -1
        tot+=s*m*m*mu
    return tot
def delta_even(d):
    ms=sorted(set(m for (m,q) in d if m%2==0))
    return {m:d.get((m,0),0)-d.get((m,1),0) for m in ms}
G=parse({"%d,%d"%(m,q):mu for (m,q,mu) in GAUGE})
print("gauge sector: D = %+.1f   delta(even m) = %s" % (D_struct(G), delta_even(G)))
print()
print("  rep   (a,b,c)   |L|par   D_fermion   max|delta(even m)|   verdict")
for k,v in LAB.items():
    a,b,c=v; d=parse(H[k]); de=delta_even(d); Df=D_struct(d)
    md=max(abs(x) for x in de.values())
    adm=(a+2*b+3*c)%2==1 and b>=1 and a+b+c>=3
    print("  %-5s %-9s %-6s %10.1f %14d      %s" %
          (k,str(v),"odd" if (a+2*b+3*c)%2 else "even",Df,md,
           "DECOUPLED (fermion contributes 0)" if Df==0 else "contributes"))
print("\nIdentity check: D = sum_{m=2 mod 4} m^2 d(m) - sum_{m=0 mod 4} m^2 d(m) ?")
for k,v in LAB.items():
    d=parse(H[k]); de=delta_even(d)
    alt=sum(m*m*de[m] for m in de if m%4==2)-sum(m*m*de[m] for m in de if m%4==0)
    print("   %-5s  D_struct=%10.1f   from-delta=%10.1f   %s" %
          (k,D_struct(d),alt,"OK" if abs(D_struct(d)-alt)<1e-9 else "MISMATCH"))
