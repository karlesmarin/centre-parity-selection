# _dichotomy.py - even-m vs odd-m pairing dichotomy, and whether AHMN's alpha2 fundamental
#   domain reduction transfers to SM-admissible reps.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
import json
H=json.load(open("hists.json"))
LAB={"35":(4,0,0),"60":(0,2,1),"84":(0,1,3),"140a":(1,1,2),"140b":(0,3,1),
     "216":(0,4,1),"224":(1,2,2),"280":(0,2,3),"360":(1,3,2)}
# gauge sector actually summed in the potential (higgs_vacii.py GAUGE list): adjoint 15 = (1,0,1)
GAUGE=[(1,0,2),(-1,0,2),(1,1,2),(-1,1,2),(2,0,1),(0,0,1),(-2,0,1)]
def parse(h):
    d={}
    for k,v in h.items():
        m,q=map(int,k.split(',')); d[(m,q)]=d.get((m,q),0)+int(v)
    return d
def pairing(d,parity):
    ms=set(m for (m,q) in d)
    bad=[(m,d.get((m,0),0),d.get((m,1),0)) for m in sorted(ms)
         if m%2==parity and d.get((m,0),0)!=d.get((m,1),0)]
    return (len(bad)==0, bad)
print("rep    (a,b,c)   |L|  par  odd-m pair (AHMN)   even-m pair (notch)  admissible")
G=parse({"%d,%d"%(m,q):mu for (m,q,mu) in GAUGE})
for k,v in LAB.items():
    a,b,c=v; L=a+2*b+3*c
    d=parse(H[k])
    po,_=pairing(d,1); pe,_=pairing(d,0)
    adm = (L%2==1) and b>=1 and a+b+c>=3
    print("  %-5s %-9s %3d %5s  %-18s %-19s %s" %
          (k,str(v),L,"odd" if L%2 else "even",po,pe,adm))
go,_=pairing(G,1); ge,_=pairing(G,0)
print("\ngauge adjoint 15 = (1,0,1), |L|=4 even : odd-m pair=%s  even-m pair=%s" % (go,ge))
print("\n--- TOTAL potential (gauge + fermion), which is what AHMN's argument needs ---")
for k,v in LAB.items():
    d=parse(H[k]); tot=dict(d)
    for (m,q,mu) in GAUGE: tot[(m,q)]=tot.get((m,q),0)+mu
    po,badO=pairing(tot,1)
    print("  %-5s total odd-m pairing (AHMN eq 4.15 precondition): %-5s %s" %
          (k,po,"" if po else "first failure (m,n_q0,n_q1)=%s"%(badO[0],)))
