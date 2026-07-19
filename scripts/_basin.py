# _basin.py - measure the SPACE, not the point: basin geometry as a robust replacement for the
#   finite-difference Hessian eigenvalues at the vacuum (Carles's "y por que no del espacio").
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
#
# V is exactly a 2D cosine series: V = sum_{m,q,k} mult*W_k*(-1)^(k2 q) cos(pi(k1 m a1 + k2 m a2)),
# since k2*q is an integer (sin term dies). So we build the Fourier coefficients and evaluate the
# whole torus by one inverse FFT -- no derivatives anywhere.
#
# For a quadratic well the sublevel set S(eps)={V < Vmin + eps*dV} is an ellipse with semi-axes
# sqrt(2 eps dV / lambda_i); a uniform ellipse has second moment M_ii = a_i^2/4, hence
#     lambda_i = eps*dV / (2*M_ii)
# giving BOTH curvature eigenvalues from an integral over the basin.
import json, numpy as np
H=json.load(open("hists.json"))
def parse(h): return [(int(k.split(',')[0]),int(k.split(',')[1]),m) for k,m in h.items()]
FERM={k:parse(v) for k,v in H.items()}
GAUGE=[(1,0,2),(-1,0,2),(1,1,2),(-1,1,2),(2,0,1),(0,0,1),(-2,0,1)]
order=["35","60","84","140a","140b","216","224","280","360"]
N=512                      # grid over the period-2 torus
def Vgrid(hf,KMAX,wg=0.35):
    A=np.zeros((N,N))
    KS=[(x,y) for x in range(-KMAX,KMAX+1) for y in range(-KMAX,KMAX+1) if (x,y)!=(0,0)]
    for hist,sgn in ((hf,1.0),(GAUGE,-wg)):
        for (m,q,mu) in hist:
            for (k1,k2) in KS:
                w=sgn*mu/ (k1*k1+k2*k2)**3 * ((-1)**((k2*q)%2))
                A[(k1*m)%N,(k2*m)%N]+=w
    return np.real(np.fft.ifft2(A))*N*N
def basin(V,eps):
    """second-moment tensor of the sublevel set around the global min; returns (lo,hi) curvatures"""
    i=np.unravel_index(np.argmin(V),V.shape); Vmin=V[i]; dV=V.max()-Vmin
    S=V < Vmin+eps*dV
    idx=np.argwhere(S)
    # minimal-image displacement on the period-2 torus (grid step = 2/N)
    d=(idx-np.array(i)+N//2)%N-N//2
    d=d*(2.0/N)
    if len(d)<12: return None
    M=(d.T@d)/len(d)
    ev=np.linalg.eigvalsh(M)          # ev = a_i^2/4 for a uniform ellipse
    ev=np.maximum(ev,1e-15)
    lam=eps*dV/(2*ev)                 # -> curvature eigenvalues, LARGEST ev = flattest dir = lo lam
    return float(lam.min()),float(lam.max()),len(d)/V.size
print("BASIN GEOMETRY vs FINITE-DIFFERENCE HESSIAN   (grid %d^2, torus period 2)"%N)
print("  %-6s %-13s %-24s %-24s" % ("rep","FD eig (K=14)","basin eps=0.05 (lo,hi)","basin eps=0.15 (lo,hi)"))
FD={"35":(482.0,689.1),"60":(391.7,1049.2),"84":(253.7,635.2),"140a":(332.1,1132.2),
    "140b":(8.3,3430.9),"216":(638.0,8941.1),"224":(411.7,5491.7),"280":(1177.9,2052.9),
    "360":(3019.1,14673.0)}
res={}
for r in order:
    V=Vgrid(FERM[r],14)
    b5=basin(V,0.05); b15=basin(V,0.15)
    res[r]=(b5,b15)
    print("  %-6s (%7.1f,%8.1f) (%8.1f,%9.1f) f=%.3f (%8.1f,%9.1f) f=%.3f" %
          (r,FD[r][0],FD[r][1],b5[0],b5[1],b5[2],b15[0],b15[1],b15[2]))
print("\nKMAX STABILITY of the basin measure (eps=0.15):")
print("  %-6s %-22s %-22s  drift" % ("rep","KMAX=4","KMAX=14"))
for r in order:
    a=basin(Vgrid(FERM[r],4),0.15); b=basin(Vgrid(FERM[r],14),0.15)
    dl=abs(a[0]-b[0])/max(abs(b[0]),1e-9)*100; dh=abs(a[1]-b[1])/max(abs(b[1]),1e-9)*100
    print("  %-6s (%8.1f,%9.1f) (%8.1f,%9.1f)  lo %5.1f%%  hi %5.1f%%"%(r,a[0],a[1],b[0],b[1],dl,dh))
