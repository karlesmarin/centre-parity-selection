# _derive_mq.sage - DERIVE the (m,q) projection from the T^2/Z2 boundary conditions instead of
#   importing it. m = weight under the diagonalized Wilson-line generator T = diag(0,1,0,-1);
#   q = phase under the translation matrix U6 = P2 P0^-1 = diag(1,1,-1,1), i.e. (-1)^n3.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
A3 = WeylCharacterRing("A3", style="coroots")
P0=vector([1,1,1,-1]); P2=vector([1,1,-1,-1])
U5=vector([P0[i]*P0[i] for i in range(4)])          # P1 P0^-1 with P1 = P0
U6=vector([P2[i]*P0[i] for i in range(4)])          # P2 P0^-1
T =vector([0,1,0,-1])                               # diagonalized Wilson-line generator
print("U5 =",list(U5),"  (identity -> no discrete twist along x5, so no q in eq 3.3)")
print("U6 =",list(U6),"  (-> phase (-1)^n3 along x6, so q = n3 mod 2 in eq 3.4)")
print("T  =",list(T), "  (-> m = n2 - n4)")
def content(w,boxes): return [ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
def hist(L):
    boxes=sum((i+1)*L[i] for i in range(3)); h={}
    for w,mult in A3(*L).weight_multiplicities().items():
        n=content(w,boxes)
        m=sum(T[i]*n[i] for i in range(4))                     # DERIVED m
        u6=prod(U6[i]**n[i] for i in range(4))                 # DERIVED U6 phase
        q=0 if u6==1 else 1
        assert q==(n[2]%2), "q != n3 mod 2"
        assert m==(n[1]-n[3]), "m != n2-n4"
        h[(int(m),int(q))]=h.get((int(m),int(q)),0)+int(mult)
    return h
h35=hist((4,0,0))
print("\n35 = (4,0,0), derived (m,q) multiplicities, m>=1 (AHMN Table 1 range):")
for m in range(1,5):
    print("   (%d,0): %2d   (%d,1): %2d" % (m,h35.get((m,0),0),m,h35.get((m,1),0)))
print("\nAHMN Table 1, 35 fermion column, decoded to multiplicities:")
print("   (1,0): 2(--)x2 + 4(--)      = 2*2+4 =  8      (1,1): 2(+-)x2 + 4(+-)      = 8")
print("   (2,0): 3(++)x2 + 5(++)      = 3*2+5 = 11      (2,1): 3(-+)               = 3")
print("   (3,0): 4(--)                =         4      (3,1): 4(+-)               = 4")
print("   (4,0): 5(++)                =         5      (4,1): -                   = 0")
tab={(1,0):8,(1,1):8,(2,0):11,(2,1):3,(3,0):4,(3,1):4,(4,0):5,(4,1):0}
ok=all(h35.get(k,0)==v for k,v in tab.items())
print("\nDERIVED == AHMN Table 1 ?", ok)
# and the gauge adjoint, whose Table 1 column is the (c>0) one
h15=hist((1,0,1))
print("adjoint 15 derived, m>=1:", {k:v for k,v in sorted(h15.items()) if k[0]>=1})
