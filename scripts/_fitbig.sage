R.<t> = LaurentPolynomialRing(QQ)
X=[R(1),t,R(-1),t**-1]
BOX=14; zero=set(); n=0
for a in range(BOX+1):
 for b in range(BOX+1):
  for c in range(BOX+1):
   if (a,b,c)==(0,0,0): continue
   n+=1
   lam=[a+b+c,b+c,c,0]
   if matrix(R,4,4,lambda i,j: X[i]**(lam[j]+3-j)).det()==0: zero.add((a,b,c))
pred=lambda a,b,c: b%2==1 and ((a%2==1 and c%2==1) or a==c)
bad=[(a,b,c) for a in range(BOX+1) for b in range(BOX+1) for c in range(BOX+1)
     if (a,b,c)!=(0,0,0) and pred(a,b,c)!=((a,b,c) in zero)]
print("BOX=%d  reps=%d  |D==0|=%d  mismatches to fit: %d" % (BOX,n,len(zero),len(bad)))
print("  all degenerate have EVEN n-ality:", all((a+2*b+3*c)%2==0 for (a,b,c) in zero))
print("  any degenerate ALSO admissible (odd & b>=1 & a+b+c>=3):",
      [z for z in zero if (z[0]+2*z[1]+3*z[2])%2==1])
