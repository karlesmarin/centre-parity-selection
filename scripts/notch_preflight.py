#!/usr/bin/env python3
"""notch_preflight.py - is your fundamental-domain reduction legal for YOUR representation?

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

WHAT IT IS
  A one-shot check to run BEFORE minimizing a Wilson-line effective potential on T^2/Z2.
  Gauge-Higgs constructions routinely halve the search region in alpha_2 using an argument
  that holds only if the spectrum pairs its twist multiplicities at ODD m. That hypothesis
  is not a property of the model: it is the parity of the representation's centre charge,
  and it INVERTS on every representation able to host a Standard-Model quark generation.
  A minimization run on the halved region for such a representation can return a boundary
  point that is not a vacuum.

WHAT IT ANSWERS
  Given SU(4) Dynkin labels (a,b,c) -- or an explicit (m,q) multiplicity table for any
  spectrum -- it reports which pairing holds, hence which fundamental region is legal, and
  whether the representation can host a quark cell at all.

USAGE
  python notch_preflight.py 0 2 1              # by Dynkin labels (needs SageMath for weights)
  python notch_preflight.py --hist hists.json --key 60      # from a precomputed (m,q) table
  python notch_preflight.py --scan 6           # sweep all labels up to a+b+c <= 6

THE RULE IT IMPLEMENTS (proved; see NOTCH_THEOREM.md)
  D(t) = sum_weights t^(n2-n4) (-1)^n3 = s_lambda(1,-1,t,1/t)  obeys  D(-t) = (-1)^|L| D(t),
  |L| = a+2b+3c.  Hence
      |L| odd   -> pairs at EVEN m -> odd-m hypothesis FAILS -> full torus required
      |L| even  -> pairs at ODD  m -> odd-m hypothesis HOLDS -> half domain legal
  with both holding on the degenerate class Z = { b odd and ((a,c both odd) or a=c even) }.
  Admissibility (Part II): |L| odd and b>=1 and a+b+c>=3 -- so ADMISSIBLE ALWAYS IMPLIES
  the half domain is illegal.
"""
import sys, json, argparse

def nality(a, b, c):            return (a + 2*b + 3*c) % 2
def degenerate(a, b, c):        return b % 2 == 1 and ((a % 2 == 1 and c % 2 == 1) or a == c)
def admissible(a, b, c):        return nality(a,b,c) == 1 and b >= 1 and a + b + c >= 3

def verdict(a, b, c):
    odd = nality(a,b,c) == 1
    deg = degenerate(a,b,c)
    return {
        "labels": (a,b,c), "centre_charge_parity": "odd" if odd else "even",
        "pairs_at": "both (degenerate)" if deg else ("even m" if odd else "odd m"),
        "odd_m_hypothesis": "HOLDS" if (not odd or deg) else "FAILS",
        "legal_domain": "alpha_2 in [0,1/2]" if (not odd or deg) else "alpha_2 in [0,1]  (FULL TORUS)",
        "hosts_quark_cell": admissible(a,b,c),
        "degenerate": deg,
    }

def from_hist(hist):
    """hist: {(m,q): multiplicity}. Returns which pairings hold, measured not assumed."""
    d = {}
    for k, v in hist.items():
        m, q = (k if isinstance(k, tuple) else tuple(map(int, k.split(','))))
        d[(int(m), int(q))] = d.get((int(m), int(q)), 0) + int(v)
    ms = sorted({m for (m, _) in d})
    bad = lambda par: [(m, d.get((m,0),0), d.get((m,1),0)) for m in ms
                       if m % 2 == par and d.get((m,0),0) != d.get((m,1),0)]
    return {"even_m_notch": not bad(0), "odd_m_pairing": not bad(1),
            "first_odd_failure": (bad(1)[0] if bad(1) else None),
            "first_even_failure": (bad(0)[0] if bad(0) else None)}

def report(v):
    print(f"  labels (a,b,c)        : {v['labels']}")
    print(f"  centre charge a+2b+3c : {v['centre_charge_parity']}")
    print(f"  tower pairs at        : {v['pairs_at']}")
    print(f"  odd-m hypothesis      : {v['odd_m_hypothesis']}")
    print(f"  LEGAL SEARCH REGION   : {v['legal_domain']}")
    print(f"  hosts a quark cell    : {v['hosts_quark_cell']}")
    if v['hosts_quark_cell']:
        print("  >> WARNING: this representation can host the Standard Model, therefore its")
        print("     centre charge is odd, therefore the half-domain reduction is NOT valid.")
        print("     Scan the full torus, and check the Hessian is positive semi-definite.")

def main():
    p = argparse.ArgumentParser(description="fundamental-domain preflight for T^2/Z2 gauge-Higgs")
    p.add_argument("labels", nargs="*", type=int, help="Dynkin labels a b c")
    p.add_argument("--hist", help="JSON file of (m,q) multiplicity tables")
    p.add_argument("--key", help="which entry of --hist to test")
    p.add_argument("--scan", type=int, help="sweep all labels with a+b+c <= N")
    args = p.parse_args()
    if args.hist:
        H = json.load(open(args.hist))
        keys = [args.key] if args.key else list(H)
        for k in keys:
            print(f"\n[{k}]  measured directly from the (m,q) table:")
            for kk, vv in from_hist(H[k]).items(): print(f"  {kk:20s}: {vv}")
        return
    if args.scan:
        N = args.scan
        print("  (a,b,c)      centre  pairs at            legal region             quark cell")
        for a in range(N+1):
            for b in range(N+1):
                for c in range(N+1):
                    if not (1 <= a+b+c <= N): continue
                    v = verdict(a,b,c)
                    print("  %-12s %-7s %-19s %-24s %s" % (str((a,b,c)),
                          v['centre_charge_parity'], v['pairs_at'], v['legal_domain'],
                          "YES" if v['hosts_quark_cell'] else ""))
        return
    if len(args.labels) != 3:
        p.print_help(); sys.exit(1)
    print()
    report(verdict(*args.labels))

if __name__ == "__main__":
    main()
