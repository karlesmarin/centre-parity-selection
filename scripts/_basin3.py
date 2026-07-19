# _basin3.py - basin measure restricted to the CONNECTED COMPONENT of the global minimum.
#   Fixes _basin2: the torus has symmetry-equivalent minima, so an unrestricted sublevel set
#   measured the separation BETWEEN minima (M saturated, lambda scaled as eps). Validate against
#   the finite-difference Hessian where FD is trustworthy, then read 140b off the plateau.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
import json, numpy as np
from scipy import ndimage
exec(open("_basin.py").read().split("def basin")[0])
FD={"35":(482.0,689.1),"60":(391.7,1049.2),"84":(253.7,635.2),"140a":(332.1,1132.2),
    "140b":(8.3,3430.9),"216":(638.0,8941.1),"224":(411.7,5491.7),"280":(1177.9,2052.9),
    "360":(3019.1,14673.0)}
def basinC(V,eps):
    i=np.unravel_index(np.argmin(V),V.shape); Vmin=V[i]; dV=V.max()-Vmin
    # roll the global min to the centre so the component cannot wrap the boundary
    Vr=np.roll(V,(N//2-i[0],N//2-i[1]),axis=(0,1))
    S=Vr<Vmin+eps*dV
    lab,_=ndimage.label(S)
    comp=lab==lab[N//2,N//2]
    idx=np.argwhere(comp)
    if len(idx)<12: return None
    d=(idx-np.array([N//2,N//2]))*(2.0/N)
    M=(d.T@d)/len(d)
    ev=np.maximum(np.linalg.eigvalsh(M),1e-18)
    lam=eps*dV/(2*ev)
    return float(lam.min()),float(lam.max()),len(idx)
EPS=[1e-4,3e-4,1e-3,3e-3,1e-2,3e-2]
V={r:Vgrid(FERM[r],14) for r in order}
for which,lab in ((0,"LOW"),(1,"HIGH")):
    print("\n%s eigenvalue: basin(connected) vs FD.   [plateau = quadratic regime]"%lab)
    print("  %-6s %9s | %s" % ("rep","FD"," ".join("%9s"%("e=%.0e"%e) for e in EPS)))
    for r in order:
        row=[]
        for e in EPS:
            b=basinC(V[r],e); row.append("%9.1f"%b[which] if b else "       --")
        print("  %-6s %9.1f | %s" % (r,FD[r][which]," ".join(row)))
print("\nbasin point counts (eps=1e-3 / 1e-2):",
      {r:(basinC(V[r],1e-3)[2] if basinC(V[r],1e-3) else 0,
          basinC(V[r],1e-2)[2] if basinC(V[r],1e-2) else 0) for r in order})
