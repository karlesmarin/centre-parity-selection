# _basin4.py - KMAX stability of the VALIDATED basin measure (connected component, eps in the
#   quadratic plateau). This is the number that decides whether the basin replaces the Hessian.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
import numpy as np
exec(open("_basin3.py").read().split("EPS=[")[0])
print("KMAX stability at eps=5e-4 (quadratic plateau), connected component")
print("  %-6s %-20s %-20s %-14s %s" % ("rep","KMAX=4 (lo,hi)","KMAX=14 (lo,hi)","drift lo/hi","FD-Hessian drift"))
FDD={"35":"0.1%","60":"0.1%","84":"0.7%","140a":"0.4%","140b":"68.9%","216":"0.9%",
     "224":"4.5%","280":"1.2%","360":"0.1%"}
for r in order:
    a=basinC(Vgrid(FERM[r],4),5e-4); b=basinC(Vgrid(FERM[r],14),5e-4)
    dl=abs(a[0]-b[0])/max(abs(b[0]),1e-9)*100; dh=abs(a[1]-b[1])/max(abs(b[1]),1e-9)*100
    print("  %-6s (%7.1f,%8.1f) (%7.1f,%8.1f) %5.1f%% /%5.1f%%   %s" %
          (r,a[0],a[1],b[0],b[1],dl,dh,FDD[r]))
