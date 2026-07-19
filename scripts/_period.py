# _period.py - SOCRATIC ROOT CHECK. We verified AHMN's HYPOTHESIS (odd-m pairing) fails for
#   admissible reps and concluded their CONCLUSION (period 1 in alpha2) fails. That is denying
#   the antecedent. Test the conclusion DIRECTLY on the potential itself.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
import numpy as np
exec(open("_basin.py").read().split("def basin")[0])   # Vgrid, FERM, GAUGE, order, N
LAB={"35":(4,0,0),"60":(0,2,1),"84":(0,1,3),"140a":(1,1,2),"140b":(0,3,1),
     "216":(0,4,1),"224":(1,2,2),"280":(0,2,3),"360":(1,3,2)}
# grid spans alpha in [0,2) with N points -> shift by 1 in alpha2 == roll by N//2 along axis 1
print("Is V actually periodic with period 1 in alpha2?  (AHMN's CONCLUSION, tested directly)")
print("  %-6s %-9s %-6s %12s %12s %12s   %s" %
      ("rep","(a,b,c)","|L|","max|V(a2+1)-V|","scale max|V|","rel","period-1?"))
for r in order:
    a,b,c=LAB[r]; V=Vgrid(FERM[r],10)
    d=np.max(np.abs(np.roll(V,N//2,axis=1)-V)); s=np.max(np.abs(V-V.mean()))
    per = d/s < 1e-10
    print("  %-6s %-9s %-6s %12.4e %12.4e %12.2e   %s" %
          (r,str((a,b,c)),"odd" if (a+2*b+3*c)%2 else "even",d,s,d/s,"YES" if per else "NO"))
print("\nAlso test the other candidate symmetries on the SAME grid:")
print("  %-6s %-14s %-14s %-14s" % ("rep","a2 -> -a2","a1 -> a1+1","(a1,a2)->(a1+1,a2+1)"))
for r in order:
    V=Vgrid(FERM[r],10); s=np.max(np.abs(V-V.mean()))
    f=lambda W: "%.1e"%(np.max(np.abs(W-V))/s)
    print("  %-6s %-14s %-14s %-14s" %
          (r, f(V[:,::-1]), f(np.roll(V,N//2,axis=0)), f(np.roll(np.roll(V,N//2,axis=0),N//2,axis=1))))
