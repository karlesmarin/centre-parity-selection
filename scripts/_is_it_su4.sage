# _is_it_su4.sage - is the centre-charge selection rule an accident of SU(4)?
#   Derived condition: the alphabet A = {eps_i t^{c_i}} (eps from the Z2 twist U6, c from the
#   Wilson-line charge) must satisfy  -A = A|_{t->-t}  as multisets, which is what makes the
#   shift equal a central element up to Weyl. Entries with c = +-1 pair automatically; entries
#   with c = 0 need their twist signs balanced. Brute-force every alphabet for SU(3)..SU(8).
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
from itertools import product as iproduct
R.<t> = LaurentPolynomialRing(QQ)

def works(eps, c):
    """does -A equal A at t -> -t, as multisets?"""
    A  = [eps[i]*t**c[i] for i in range(len(eps))]
    nA = sorted([-a for a in A], key=str)
    sA = sorted([eps[i]*(-1)**c[i]*t**c[i] for i in range(len(eps))], key=str)
    return nA == sA

print("Sweeping every (twist signs, Wilson charges) alphabet for SU(N), N = 3..8.")
print("Wilson charges drawn from {0,+1,-1} with equal numbers of +1 and -1 (a Cartan direction),")
print("twist signs from {+1,-1}. Reporting which N admit ANY alphabet with the property.\n")
print("  %-4s %-10s %-12s %-10s %s" % ("N","#alphabets","#with rule","neutral slots","verdict"))
summary = {}
for N in range(3, 9):
    total = good = 0
    good_examples = []
    for c in iproduct([0, 1, -1], repeat=N):
        if sum(1 for x in c if x == 1) != sum(1 for x in c if x == -1):  continue
        if all(x == 0 for x in c):                                       continue
        for eps in iproduct([1, -1], repeat=N):
            total += 1
            if works(eps, c):
                good += 1
                if len(good_examples) < 2: good_examples.append((eps, c))
    # how many neutral slots do the successful ones have?
    ns = sorted({sum(1 for x in c if x == 0) for _, c in good_examples}) if good_examples else []
    summary[N] = good
    print("  %-4d %-10d %-12d %-10s %s" % (N, total, good, ns if ns else "-",
          "rule EXISTS" if good else "no such alphabet"))

print("\n  N with a selection rule :", [n for n, g in summary.items() if g])
print("  N without               :", [n for n, g in summary.items() if not g])

print("\nAnd the SU(4) case we actually have:")
eps = (1, 1, -1, 1); c = (0, 1, 0, -1)     # U6 = diag(1,1,-1,1), Wilson along diag(0,1,0,-1)
print("   twist signs", eps, " Wilson charges", c, " -> rule holds:", works(eps, c))
print("   neutral slots:", [i+1 for i in range(4) if c[i] == 0],
      "with twist signs", [eps[i] for i in range(4) if c[i] == 0], "-- balanced, as required")
