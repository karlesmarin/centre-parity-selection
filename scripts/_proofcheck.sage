R.<t> = LaurentPolynomialRing(QQ)
X=[R(1),t,R(-1),t**-1]   # order: x1=1, x2=t, x3=-1, x4=1/t  -> reorder to (1,-1,t,1/t)
Y=[R(1),R(-1),t,t**-1]
def N(a,b,c):
    lam=[a+b+c,b+c,c,0]
    return matrix(R,4,4,lambda i,j: Y[i]**(lam[j]+3-j)).det()
def br(x): return t**x - t**(-x)
def PQR(a,b,c): return (a+b+c+3, b+c+2, c+1)
bad=[]
for a in range(9):
 for b in range(9):
  for c in range(9):
   if (a,b,c)==(0,0,0): continue
   P,Q,Rr=PQR(a,b,c); n=N(a,b,c)
   al,be,ga=a%2,b%2,c%2
   if be==1 and ga==1 and al==1: pred=R(0)                      # |E|=4 : all mu even
   elif be==1 and ga==0 and al==0: pred=2*(br(Q)-br(Rr)+br(P-Q)-br(P-Rr))   # E={1,4}
   elif be==0 and ga==1 and al==1: pred=2*(-br(Q)+br(Q-Rr)+br(P)-br(P-Rr))  # E={3,4}
   elif be==0 and ga==0 and al==0: pred=2*(br(Rr)+br(Q-Rr)-br(P)+br(P-Q))   # E={2,4}
   else: pred=None                                              # |E| in {1,3}: nonvanishing arg
   if pred is not None and n!=pred: bad.append((a,b,c,"formula"))
   if pred is None and n==0: bad.append((a,b,c,"|E|=1,3 vanished!"))
print("checked a,b,c<=8. failures:",len(bad)); print(bad[:15])
# the E={1,4} vanishing criterion: {Q,P-Q}=={R,P-R}  <=> a==c
print("E={1,4} branch: N==0 <=> a==c ?",
  all(((N(a,b,c)==0)==(a==c)) for a in range(0,9,2) for b in range(1,9,2) for c in range(0,9,2)))
