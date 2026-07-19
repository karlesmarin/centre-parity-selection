# _central_element2.sage - state the principle PRECISELY. The naive claim "alpha_2 -> alpha_2+1
#   acts as -1 in SU(4)" is not literally true as matrices; it is true up to a Weyl element,
#   which is invisible to a character. Pin down exactly which statement holds.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
R.<t> = LaurentPolynomialRing(QQ)

# the group element the orbifold hands us, slot by slot:
#   slot 1 : 1          (free index)
#   slot 2 : t          )  Wilson line along T = diag(0,1,0,-1)
#   slot 3 : -1         (the U6 twist, diag(1,1,-1,1))
#   slot 4 : 1/t        )
g      = diagonal_matrix(R, [1,  t, -1,  t**-1])
g_shift= diagonal_matrix(R, [1, -t, -1, -t**-1])   # the same at alpha_2 -> alpha_2 + 1  (t -> -t)
minus1 = -identity_matrix(R, 4)                    # the central element -1 of SU(4)
prod   = minus1 * g                                # what a central transformation would give

print("g(alpha2)        diag:", g.diagonal())
print("g(alpha2 + 1)    diag:", g_shift.diagonal())
print("(-1) * g(alpha2) diag:", prod.diagonal())
print()
print("equal as MATRICES  :", g_shift == prod)
print("equal as MULTISETS :", sorted(g_shift.diagonal(), key=str) == sorted(prod.diagonal(), key=str))

# find the permutation relating them
from itertools import permutations
d1 = g_shift.diagonal(); d2 = prod.diagonal()
perms = [p for p in permutations(range(4)) if [d2[p[i]] for i in range(4)] == d1]
print("\npermutations P with  (P . (-1)g . P^-1) = g(alpha2+1):", perms)
print("  -> a Weyl element of SU(4) (S4 acts by permuting the 4 eigenvalues)")
print()
print("PRECISE STATEMENT:")
print("  g(alpha_2 + 1)  is CONJUGATE to  (-1) * g(alpha_2)  by a Weyl element.")
print("  Characters are class functions, so for every irrep lambda:")
print("      chi_lambda(g(alpha_2+1)) = chi_lambda(-1 . g(alpha_2)) = (-1)^|lambda| chi_lambda(g(alpha_2))")
print("  which is exactly the parity lemma. The Weyl conjugation is invisible to chi.")
print()
# and check the naive matrix-level claim really does fail, so we do not overstate it
print("SANITY: is the shift itself central?  exp-style element diag(1,-1,1,-1) central?",
      diagonal_matrix(QQ,[1,-1,1,-1]) in [z*identity_matrix(QQ,4) for z in (1,-1)])
print("  -> No. Only the COMBINATION with a Weyl reflection is the central action.")
