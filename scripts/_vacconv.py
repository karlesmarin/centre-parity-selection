# _vacconv.py - does the TRUE full-torus vacuum move with KMAX? (engine default is 4)
#   The naturalness ranking is computed from the vacuum, so this is the decisive convergence test.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
import json, numpy as np
H=json.load(open("hists.json"))
def parse(h): return [(int(k.split(',')[0]),int(k.split(',')[1]),m) for k,m in h.items()]
FERM={k:parse(v) for k,v in H.items()}
GAUGE=[(1,0,2),(-1,0,2),(1,1,2),(-1,1,2),(2,0,1),(0,0,1),(-2,0,1)]
order=["35","60","84","140a","140b","216","224","280","360"]
def mk(K):
    KS=[(x,y) for x in range(-K,K+1) for y in range(-K,K+1) if (x,y)!=(0,0)]
    return (np.array([k[0] for k in KS])[None,:], np.array([k[1] for k in KS])[None,:],
            (1.0/(np.array([k[0] for k in KS])**2+np.array([k[1] for k in KS])**2)**3)[None,:])
def Ugrid(hist,A1,A2,K):
    K1,K2,WK=K; out=np.zeros(A1.shape[0])
    for (m,q,mu) in hist:
        out+=mu*np.sum(np.cos(np.pi*(K1*m*A1[:,None]+K2*(m*A2[:,None]+q)))*WK,axis=1)
    return out
def U(hf,A1,A2,K,wg=0.35): return Ugrid(hf,A1,A2,K)-wg*Ugrid(GAUGE,A1,A2,K)
def vac(hf,K,n=96):
    g=np.linspace(0,1,n+1); X,Y=np.meshgrid(g,g,indexing='ij')
    A1=X.ravel(); A2=Y.ravel(); V=U(hf,A1,A2,K)
    i=int(np.argmin(V)); a1,a2=A1[i],A2[i]; d=1.0/n
    for _ in range(5):
        d*=0.3; o=np.linspace(-4,4,9)*d
        P1,P2=np.meshgrid((a1+o)%1.0,(a2+o)%1.0,indexing='ij')
        p1=P1.ravel();p2=P2.ravel(); VV=U(hf,p1,p2,K); j=int(np.argmin(VV)); a1,a2=p1[j],p2[j]
    return a1,a2
print("full-torus vacuum (a1,a2), folded a2 -> [0,1/2].  engine default is KMAX=4")
print("  %-6s %-18s %-18s %-18s" % ("rep","KMAX=4","KMAX=8","KMAX=14"))
for r in order:
    hf=FERM[r]; cells=[]
    for KM in [4,8,14]:
        a1,a2=vac(hf,mk(KM)); a2f=min(a2%1.0,1.0-(a2%1.0))
        cells.append("(%.4f, %.4f)"%(a1,a2f))
    print("  %-6s %-18s %-18s %-18s" % (r,*cells))

# --- does the RANKING input (Hessian eigenvalues at the vacuum) converge? ---
def hess(hf,a1,a2,K,h=0.004):
    f=lambda x,y: U(hf,np.array([x]),np.array([y]),K)[0]
    U11=(f(a1+h,a2)-2*f(a1,a2)+f(a1-h,a2))/h**2
    U22=(f(a1,a2+h)-2*f(a1,a2)+f(a1,a2-h))/h**2
    U12=(f(a1+h,a2+h)-f(a1+h,a2-h)-f(a1-h,a2+h)+f(a1-h,a2-h))/(4*h*h)
    tr=U11+U22; disc=np.sqrt((U11-U22)**2+4*U12*U12)
    return 0.5*(tr-disc),0.5*(tr+disc)
print("\nHessian eigenvalues at the vacuum -- THE ranking input (rho = eig/|alpha|^2)")
print("  %-6s %-22s %-22s %-22s  %s" % ("rep","KMAX=4 (lo,hi)","KMAX=8","KMAX=14","drift 4->14"))
for r in order:
    hf=FERM[r]; row=[]; 
    for KM in [4,8,14]:
        K=mk(KM); a1,a2=vac(hf,K); row.append(hess(hf,a1,a2,K))
    d_lo=abs(row[2][0]-row[0][0])/max(abs(row[2][0]),1e-9)*100
    d_hi=abs(row[2][1]-row[0][1])/max(abs(row[2][1]),1e-9)*100
    print("  %-6s %-22s %-22s %-22s  lo %5.1f%%  hi %5.1f%%" %
          (r,"(%.1f, %.1f)"%row[0],"(%.1f, %.1f)"%row[1],"(%.1f, %.1f)"%row[2],d_lo,d_hi))
