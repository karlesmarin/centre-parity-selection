# _kconv.py - is the "true" alpha2-curvature at 1/2 converged in KMAX?
#   Falsify the instrument before trusting the 4-of-8 mismatch table in NOTCH_THEOREM.md.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
import json, numpy as np
H=json.load(open("hists.json"))
def parse(h): return [(int(k.split(',')[0]),int(k.split(',')[1]),m) for k,m in h.items()]
FERM={k:parse(v) for k,v in H.items()}
GAUGE=[(1,0,2),(-1,0,2),(1,1,2),(-1,1,2),(2,0,1),(0,0,1),(-2,0,1)]
order=["35","60","84","140a","140b","216","224","280","360"]
def mk(KMAX):
    KS=[(x,y) for x in range(-KMAX,KMAX+1) for y in range(-KMAX,KMAX+1) if (x,y)!=(0,0)]
    K1=np.array([k[0] for k in KS]);K2=np.array([k[1] for k in KS])
    return K1,K2,1.0/(K1*K1+K2*K2)**3
def S(hist,a1,a2,K1,K2,WK):
    t=0.0
    for (m,q,mu) in hist: t+=mu*np.sum(np.cos(np.pi*(K1*m*a1+K2*(m*a2+q)))*WK)
    return t
def Uf(hf,a1,a2,K,wg=0.35): return S(hf,a1,a2,*K)-wg*S(GAUGE,a1,a2,*K)
def a1star(hf,K):
    best=(1e18,0)
    for a1 in np.linspace(0,1,201):
        v=Uf(hf,a1,0.5,K)
        if v<best[0]: best=(v,a1)
    return best[1]
def curv(hf,a1,K,h=0.004):
    return (Uf(hf,a1,0.5+h,K)-2*Uf(hf,a1,0.5,K)+Uf(hf,a1,0.5-h,K))/h**2
print("KMAX -> curvature d2V/da2^2 at (a1*,1/2).  Also a1* itself (it can move).")
print("  %-6s %s" % ("rep"," ".join("%11s"%("K=%d"%k) for k in [2,4,6,8,12,16])))
for r in order:
    hf=FERM[r]; row=[]
    for KM in [2,4,6,8,12,16]:
        K=mk(KM); a1=a1star(hf,K); row.append(curv(hf,a1,K))
    print("  %-6s %s" % (r," ".join("%11.1f"%v for v in row)))
