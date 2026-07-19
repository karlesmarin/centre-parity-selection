# make_release_data.py - emit the tabular data behind every table and figure of Part III as CSV,
#   so the numbers are readable without running anything.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
import json, csv, math, pathlib, numpy as np
OUT = pathlib.Path("release_iii/data"); OUT.mkdir(parents=True, exist_ok=True)
H = json.load(open("hists.json"))
F = {k: [(int(a.split(',')[0]), int(a.split(',')[1]), m) for a, m in v.items()] for k, v in H.items()}
GAUGE = [(1,0,2),(-1,0,2),(1,1,2),(-1,1,2),(2,0,1),(0,0,1),(-2,0,1)]
LAB = {"35":(4,0,0),"60":(0,2,1),"84":(0,1,3),"140a":(1,1,2),"140b":(0,3,1),
       "216":(0,4,1),"224":(1,2,2),"280":(0,2,3),"360":(1,3,2)}
ORDER = ["35","60","84","140a","140b","216","224","280","360"]
def parse(h):
    d = {}
    for k, v in h.items():
        m, q = map(int, k.split(',')); d[(m,q)] = d.get((m,q),0)+int(v)
    return d
def pairs_at(d, par):
    ms = {m for (m,_) in d}
    return all(d.get((m,0),0) == d.get((m,1),0) for m in ms if m % 2 == par)
def Dstruct(d):
    t = 0.0
    for (m,q),mu in d.items():
        if m % 2: continue
        r = m % 4
        t += (1 if ((r==2 and q==0) or (r==0 and q==1)) else -1)*m*m*mu
    return t*math.pi**2
N = 384
def Vgrid(hf, KMAX=10, wg=0.35):
    A = np.zeros((N,N))
    KS = [(x,y) for x in range(-KMAX,KMAX+1) for y in range(-KMAX,KMAX+1) if (x,y)!=(0,0)]
    for hist,sgn in ((hf,1.0),(GAUGE,-wg)):
        for (m,q,mu) in hist:
            for (k1,k2) in KS:
                A[(k1*m)%N,(k2*m)%N] += sgn*mu/(k1*k1+k2*k2)**3*((-1)**((k2*q)%2))
    return np.real(np.fft.ifft2(A))*N*N

def write(name, header, rows, note):
    with open(OUT/name, "w", newline="", encoding="utf8") as f:
        w = csv.writer(f); w.writerow(["# " + note]); w.writerow(header); w.writerows(rows)
    print("  %-34s %d rows" % (name, len(rows)))

G = parse({"%d,%d"%(m,q): mu for (m,q,mu) in GAUGE})
# --- Table 2 (the inversion) + the per-rep dichotomy -------------------------------
rows = [["gauge adjoint 15","(1,0,1)",4,"even",pairs_at(G,1),pairs_at(G,0),False]]
for r in ORDER:
    a,b,c = LAB[r]; d = parse(H[r]); L = a+2*b+3*c
    adm = L%2==1 and b>=1 and a+b+c>=3
    rows.append([r,"(%d,%d,%d)"%(a,b,c),L,"odd" if L%2 else "even",
                 pairs_at(d,1),pairs_at(d,0),adm])
write("table2_dichotomy.csv",
      ["sector","abc","boxes_a+2b+3c","parity","odd_m_pairing","even_m_notch","hosts_quark_cell"],
      rows, "Table 2: AHMN's hypothesis holds for the gauge adjoint and the 35, and fails for every admissible rep")

# --- section 5: direct periodicity test --------------------------------------------
rows = []
for r in ORDER:
    V = Vgrid(F[r]); sc = np.max(np.abs(V-V.mean()))
    rel = float(np.max(np.abs(np.roll(V,N//2,axis=1)-V))/sc)
    rows.append([r,str(LAB[r]).replace(" ",""),round(rel,5),"yes" if rel<1e-12 else "no"])
write("sec5_periodicity_residuals.csv",
      ["rep","abc","max|V(a2+1)-V|/scale","period_1_holds"], rows,
      "Section 5: measured on the potential itself, not inferred from the hypothesis failing")

# --- section 5: which lattice transformations survive --------------------------------
TR = {"identity":((1,0),(0,1)),"S":((0,1),(-1,0)),"S^2=-1":((-1,0),(0,-1)),
      "exchange a1<->a2":((0,1),(1,0)),"T: a2->a1+a2":((1,0),(1,1)),
      "T': a1->a1+a2":((1,1),(0,1)),"a1->-a1":((-1,0),(0,1)),"a2->-a2":((1,0),(0,-1))}
M = 192
def Vg2(hf, KMAX=8, wg=0.35):
    A = np.zeros((M,M))
    KS = [(x,y) for x in range(-KMAX,KMAX+1) for y in range(-KMAX,KMAX+1) if (x,y)!=(0,0)]
    for hist,sgn in ((hf,1.0),(GAUGE,-wg)):
        for (m,q,mu) in hist:
            for (k1,k2) in KS:
                A[(k1*m)%M,(k2*m)%M] += sgn*mu/(k1*k1+k2*k2)**3*((-1)**((k2*q)%2))
    return np.real(np.fft.ifft2(A))*M*M
reps5 = ["35","60","84","216","360"]
Vs = {r: Vg2(F[r]) for r in reps5}
i = np.arange(M); I,J = np.meshgrid(i,i,indexing="ij")
rows = []
for name,Mx in TR.items():
    row = [name]
    for r in reps5:
        V = Vs[r]; s = np.max(np.abs(V-V.mean()))
        W = V[(Mx[0][0]*I+Mx[0][1]*J)%M, (Mx[1][0]*I+Mx[1][1]*J)%M]
        row.append(round(float(np.max(np.abs(W-V))/s), 5))
    rows.append(row)
write("sec5_modular_residuals.csv", ["transformation"]+reps5, rows,
      "Section 5: residual 0 means the transformation is a symmetry of V. Only the reflections and S^2 survive")

# --- section 6: the decoupling ------------------------------------------------------
Dg = Dstruct(G); rows = []
for r in ORDER:
    d = parse(H[r]); Df = Dstruct(d)
    rows.append([r,str(LAB[r]).replace(" ",""),round(Df,4),round(Df-0.35*Dg,4)])
write("sec6_decoupling.csv", ["rep","abc","D_fermion","D_full_gauge_plus_fermion"], rows,
      "Section 6: the fermion contribution is identically zero for every admissible rep; D_full is gauge-only")

# --- the even-m imbalance behind figure 3 -------------------------------------------
rows = []
for r in ORDER:
    d = parse(H[r])
    for m in sorted({m for (m,_) in d}):
        rows.append([r,m,d.get((m,0),0),d.get((m,1),0),d.get((m,0),0)-d.get((m,1),0)])
write("fig3_delta_m.csv", ["rep","m","mult_q0","mult_q1","delta_m"], rows,
      "Figure 3: the signed imbalance delta(m) whose vanishing pattern is the notch")

# --- the classification landscape behind figure 4 ------------------------------------
inZ = lambda a,b,c: b%2==1 and ((a%2==1 and c%2==1) or a==c)
adm = lambda a,b,c: (a+2*b+3*c)%2==1 and b>=1 and a+b+c>=3
rows = []
for a in range(13):
    for b in range(13):
        for c in range(13):
            if a+b+c == 0: continue
            rows.append([a,b,c,a+2*b+3*c,"odd" if (a+2*b+3*c)%2 else "even",
                         inZ(a,b,c), adm(a,b,c),
                         "even m" if (a+2*b+3*c)%2 else "odd m"])
write("fig4_classification.csv",
      ["a","b","c","boxes","parity","in_Z_degenerate","hosts_quark_cell","pairs_at"], rows,
      "Figure 4: the full classification to a,b,c <= 12. Z never intersects the quark-hosting set")
print("\ndata/ written")
