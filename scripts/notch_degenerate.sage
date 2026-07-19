# notch_degenerate.sage - characterize the degenerate class D(t) == 0 for the SU(4) notch
#   Author: Carles Marin  <karlesmarin@gmail.com>   (with Claude, Anthropic, as assistant)
# D(t) = s_lambda(1, t, -1, 1/t), lambda = (a+b+c, b+c, c, 0).  Bialternant, exact in QQ(t).
# Notch(even-m) holds  <=>  D odd in t  <=>  (a+2b+3c odd)  OR  D == 0.  Here we fit D==0.
R.<t> = LaurentPolynomialRing(QQ)
X = [R(1), t, R(-1), t**-1]
def D(a,b,c):
    lam=[a+b+c, b+c, c, 0]
    M=matrix(R,4,4,lambda i,j: X[i]**(lam[j]+3-j))
    V=matrix(R,4,4,lambda i,j: X[i]**(3-j))
    q=M.det()/V.det()
    return R(q) if V.det().is_unit() or True else q
BOX=10
zero=[]; nz=0
for a in range(BOX+1):
  for b in range(BOX+1):
    for c in range(BOX+1):
      if (a,b,c)==(0,0,0): continue
      num=matrix(R,4,4,lambda i,j: X[i]**([a+b+c,b+c,c,0][j]+3-j)).det()
      if num==0: zero.append((a,b,c))
      else: nz+=1
print("box<=%d : %d reps, D==0 for %d" % (BOX,nz+len(zero),len(zero)))
print("D==0 set:", zero)
# candidate fits
def fits(pred): return all((( (a,b,c) in set(zero) )==pred(a,b,c))
        for a in range(BOX+1) for b in range(BOX+1) for c in range(BOX+1) if (a,b,c)!=(0,0,0))
cands={
 "b odd and a==c"                : lambda a,b,c: b%2==1 and a==c,
 "b odd and a,c both odd"        : lambda a,b,c: b%2==1 and a%2==1 and c%2==1,
 "b odd and (a==c or (a odd and c odd))": lambda a,b,c: b%2==1 and (a==c or (a%2==1 and c%2==1)),
 "b odd and a+c even and |a-c| in {0} or a,c odd": lambda a,b,c: b%2==1 and ((a%2==1 and c%2==1) or a==c),
 "b odd and a==c mod 2 and a==c when even": lambda a,b,c: b%2==1 and ((a%2)==(c%2)) and (a%2==1 or a==c),
}
for name,f in cands.items(): print("  fit? %-45s %s" % (name, fits(f)))
print("counterexamples to 'b odd and a==c(mod2)':",
      [(a,b,c) for a in range(BOX+1) for b in range(BOX+1) for c in range(BOX+1)
       if (a,b,c)!=(0,0,0) and (b%2==1 and (a%2)==(c%2)) != ((a,b,c) in set(zero))][:20])
