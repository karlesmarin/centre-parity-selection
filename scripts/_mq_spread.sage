# _mq_spread.sage - why does #weights per (m,q) equal AHMN's #SU(2)_L multiplets per (m,q)?
#   Hypothesis: T=diag(0,1,0,-1) does not commute with SU(2)_L (which acts on indices 1,2), so a
#   multiplet of size s spreads over s consecutive m, contributing exactly ONE weight to each.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
A3 = WeylCharacterRing("A3", style="coroots")
T=[0,1,0,-1]
def data(L):
    boxes=sum((i+1)*L[i] for i in range(3)); out=[]
    for w,mult in A3(*L).weight_multiplicities().items():
        n=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        # SU(2)_L acts on indices 1,2: multiplet is labelled by (n1+n2, n3, n4); T3 = (n1-n2)/2
        out.append((int(n[0]+n[1]),int(n[2]),int(n[3]),int(n[0]-n[1]),int(n[1]-n[3]),int(n[2]%2),int(mult)))
    return out
d=data((4,0,0))
from collections import defaultdict
mult_of=defaultdict(list)          # SU(2)_L multiplet key -> list of (T3, m, q)
for (s,n3,n4,t3,m,q,mu) in d:
    for _ in range(mu): mult_of[(s,n3,n4)].append((t3,m,q))
print("35 = (4,0,0): each SU(2)_L multiplet, its size, and the m values its members occupy")
print("  %-14s %-5s %-6s %s" % ("(n1+n2,n3,n4)","size","q","m values"))
bad=0
for k in sorted(mult_of):
    v=sorted(mult_of[k]); ms=[x[1] for x in v]; qs=set(x[2] for x in v)
    onecons = (len(set(ms))==len(ms)) and (max(ms)-min(ms)==len(ms)-1)
    if not onecons: bad+=1
    print("  %-14s %-5d %-6s %s%s" % (str(k),len(v),str(sorted(qs)),ms,"" if onecons else "   <-- NOT one-per-consecutive-m"))
print("\nmultiplets violating 'one weight per consecutive m':",bad)
# consequence: count of weights at (m,q) == count of multiplets spanning that m with that q
from collections import Counter
wc=Counter(); mc=Counter()
for (s,n3,n4,t3,m,q,mu) in d: wc[(m,q)]+=mu
for k in mult_of:
    v=mult_of[k]; q=v[0][2]
    for m in set(x[1] for x in v): mc[(m,q)]+=1
print("weights per (m,q) == multiplets spanning (m,q) ?", wc==mc)
print("m>=1 slice:", {k:wc[k] for k in sorted(wc) if k[0]>=1})
