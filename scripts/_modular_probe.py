# _modular_probe.py - which modular / lattice transformations does V(alpha1,alpha2) actually respect?
#   The orbifold parity choice is asymmetric between the two cycles (U5=1, U6=diag(1,1,-1,1)),
#   so the square torus's modular symmetry cannot survive intact. Which subgroup does?
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
import json, numpy as np
H=json.load(open("hists.json"))
F={k:[(int(a.split(',')[0]),int(a.split(',')[1]),m) for a,m in v.items()] for k,v in H.items()}
GAUGE=[(1,0,2),(-1,0,2),(1,1,2),(-1,1,2),(2,0,1),(0,0,1),(-2,0,1)]
N=192
def Vgrid(hf,KMAX=8,wg=0.35):
    A=np.zeros((N,N))
    KS=[(x,y) for x in range(-KMAX,KMAX+1) for y in range(-KMAX,KMAX+1) if (x,y)!=(0,0)]
    for hist,sgn in ((hf,1.0),(GAUGE,-wg)):
        for (m,q,mu) in hist:
            for (k1,k2) in KS:
                A[(k1*m)%N,(k2*m)%N]+=sgn*mu/(k1*k1+k2*k2)**3*((-1)**((k2*q)%2))
    return np.real(np.fft.ifft2(A))*N*N
def resid(V,W): 
    s=np.max(np.abs(V-V.mean())); return np.max(np.abs(W-V))/s
# index maps on the period-2 grid (grid step 2/N; alpha = 2*i/N)
def apply(V,M):
    """(a1,a2) -> M @ (a1,a2) mod 2, as an index permutation"""
    i=np.arange(N); I,J=np.meshgrid(i,i,indexing="ij")
    I2=(M[0][0]*I+M[0][1]*J)%N; J2=(M[1][0]*I+M[1][1]*J)%N
    return V[I2,J2]
TESTS={
 "identity"                    : ((1,0),(0,1)),
 "S : (a1,a2)->(a2,-a1)"       : ((0,1),(-1,0)),
 "S^2 = -1 : (a1,a2)->(-a1,-a2)": ((-1,0),(0,-1)),
 "exchange (a1<->a2)"          : ((0,1),(1,0)),
 "T : (a1,a2)->(a1,a1+a2)"     : ((1,0),(1,1)),
 "T' : (a1,a2)->(a1+a2,a2)"    : ((1,1),(0,1)),
 "a1 -> -a1"                   : ((-1,0),(0,1)),
 "a2 -> -a2"                   : ((1,0),(0,-1)),
}
print("residual = max|V(g.alpha)-V(alpha)| / scale.   0 => the transformation is a symmetry")
print("  %-32s %s" % ("transformation"," ".join("%9s"%r for r in ["35","60","84","216","360"])))
Vs={r:Vgrid(F[r]) for r in ["35","60","84","216","360"]}
for name,M in TESTS.items():
    row=[resid(Vs[r],apply(Vs[r],M)) for r in Vs]
    print("  %-32s %s" % (name," ".join(("%9.2e"%x if x>1e-12 else "     0   ") for x in row)))
