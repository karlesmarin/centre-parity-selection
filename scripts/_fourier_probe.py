# _fourier_probe.py - is the notch literally an empty sublattice in Fourier space?
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
import json, numpy as np
H=json.load(open("hists.json"))
F={k:[(int(a.split(',')[0]),int(a.split(',')[1]),m) for a,m in v.items()] for k,v in H.items()}
order=["35","60","84","140a","140b","216","224","280","360"]
KM=6
def coeffs(hf):
    """amplitude at each (n1,n2)=(k1*m,k2*m), split by parity of k2"""
    A={}
    for (m,q,mu) in hf:
        for k1 in range(-KM,KM+1):
            for k2 in range(-KM,KM+1):
                if (k1,k2)==(0,0): continue
                w=mu/(k1*k1+k2*k2)**3*((-1)**((k2*q)%2))
                A[(k1*m,k2*m,k2%2)]=A.get((k1*m,k2*m,k2%2),0.0)+w
    return A
print("total |amplitude| carried by modes with k2 ODD, split by parity of m")
print("  %-6s %14s %14s   %s" % ("rep","m even","m odd","verdict"))
for r in order:
    A=coeffs(F[r]); ev=od=0.0
    for (m,q,mu) in F[r]:
        pass
    # regroup: recompute per (m,k2parity)
    ev=od=0.0
    for (m,q,mu) in F[r]:
        for k1 in range(-KM,KM+1):
            for k2 in range(-KM,KM+1):
                if (k1,k2)==(0,0) or k2%2==0: continue
                w=mu/(k1*k1+k2*k2)**3*((-1)**((k2*q)%2))
                if m%2==0: ev+=w
                else: od+=w
    print("  %-6s %14.6f %14.6f   %s" % (r,ev,od,"EMPTY at even m" if abs(ev)<1e-12 else "populated"))
