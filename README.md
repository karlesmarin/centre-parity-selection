# 🚪 A Centre-Charge Selection Rule for the Wilson-Line Potential

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21438227-1B6F8C?logo=doi&logoColor=white)](https://doi.org/10.5281/zenodo.21438227)
[![License](https://img.shields.io/badge/License-Apache_2.0-B5530F)](LICENSE)
[![Lean](https://img.shields.io/badge/Lean_4-sorry--free-2C2C2C?logo=lean)](lean/NotchCentreCharge.lean)
[![Verified](https://img.shields.io/badge/SageMath-exact_to_10648_reps-2C2C2C)](https://www.sagemath.org/)
[![Language](https://img.shields.io/badge/paper-EN_%2B_ES-1B6F8C)](.)

**The fundamental domain of gauge–Higgs unification is representation-dependent.**

**📄 Paper (EN + ES), Lean certificate and every verification script on Zenodo → https://doi.org/10.5281/zenodo.21438227**

> ### 📚 Part **III** of a series
> - **Part I — *Anomaly- and Tadpole-Compatible Fermion Completion of 6D SU(4) GHU***
>   → [github.com/karlesmarin/ghu-su4-completion](https://github.com/karlesmarin/ghu-su4-completion) · [Zenodo 10.5281/zenodo.21432625](https://doi.org/10.5281/zenodo.21432625)
> - **Part II — *Three Gates to a Quark Generation***
>   → [github.com/karlesmarin/su4-sm-cell-criterion](https://github.com/karlesmarin/su4-sm-cell-criterion) · [Zenodo 10.5281/zenodo.21432627](https://doi.org/10.5281/zenodo.21432627)
> - **Part III — *A Centre-Charge Selection Rule for the Wilson-Line Potential*** (this repo)

## 🎯 The principle

> **Centre-Parity Selection.** The $\mathbb{Z}_2$ reduction of the $SU(4)$ centre charge
> determines which Fourier sector of the Wilson-line potential is allowed to be non-zero.
>
> Advancing the Wilson line by one period is, **up to a Weyl reflection**, multiplication by
> the central element $-\mathbf{1}\in Z(SU(4))$; a representation answers that element with
> the scalar $(-1)^{a+2b+3c}$; and the sector of the potential carrying the opposite sign is
> therefore *identically empty* — not suppressed, absent.

Matter admissibility, Kaluza–Klein twist pairing, the Fourier support of the potential and
the geometry of its vacuum are **four readings of that one bit**.

For fifty years the centre of the gauge group has been read as a restriction on *matter* —
which representations may carry which charges. Here it also restricts the *potential*.

## ⚠️ The consequence, and why it matters if you are building one of these models

The Higgs in these models is a Wilson line. AHMN halve the search region to
$\alpha_2\in[0,\tfrac12]$ using a hypothesis they verify by inspecting their Table 1: equal
multiplicities of $(m,0)$ and $(m,1)$ at **odd** $m$. That hypothesis is not a property of
the model — it is decided by the centre charge, and **it inverts on exactly the class of
representations able to host the Standard Model** (Part II forces their centre charge odd,
while the gauge adjoint is unavoidably even).

Verified directly on the potential: the $\mathbf{35}$ is period-1 in $\alpha_2$ *exactly*,
while every admissible representation violates it by **2–30 % of its own scale**.

**The half-domain reduction does not transfer. The full torus must be scanned.**

This is not an erratum for AHMN, whose field content satisfies its own hypothesis. It is a
warning for everything downstream, and the failure is silent: a minimizer restricted to the
half-domain returns a boundary point that looks like a vacuum. We paid this ourselves.

## 🔧 The tool — run it before you minimize

```bash
python scripts/notch_preflight.py 0 2 1                        # by Dynkin labels
python scripts/notch_preflight.py --hist hists.json --key 60   # measured from an (m,q) table
```

```
labels (a,b,c)        : (0, 2, 1)
centre charge a+2b+3c : odd
tower pairs at        : even m
odd-m hypothesis      : FAILS
LEGAL SEARCH REGION   : alpha_2 in [0,1]  (FULL TORUS)
hosts a quark cell    : True
>> WARNING: this representation can host the Standard Model, therefore its
   centre charge is odd, therefore the half-domain reduction is NOT valid.
```

It costs one parity bit.

## 🌍 Is it an accident of SU(4)?

No, and not universal either. The shift equals a central element only when the
Wilson-line-neutral slots carry balanced twist signs, which for a single Wilson-line
direction requires **$N$ even**. Brute-forced over every alphabet of this shape:

| $SU(N)$ | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|
| alphabets admitting the rule | 0 | **192** | 0 | **4880** | 0 | **134400** |
| verdict | none | **exists** | none | **exists** | none | **exists** |

**$SU(4)$ is the smallest group that can carry this selection rule.**

## ✅ Machine-checked core

`lean/NotchCentreCharge.lean` — sorry-free; axioms `propext`, `Classical.choice`, `Quot.sound` only.

| theorem | content |
|---|---|
| `notch_degenerate_iff` | $N\,a\,b\,c=0 \leftrightarrow (\mathrm{Odd}\,b \wedge ((\mathrm{Odd}\,a\wedge\mathrm{Odd}\,c)\vee a=c))$ |
| `notch_parity` | the parity lemma, division-free on the bialternant |
| `degenerate_centre_charge_even` | every member of $\mathcal{Z}$ has $a+2b+3c$ even |
| `degenerate_disjoint_admissible` | $\lvert\lambda\rvert$ odd $\Rightarrow N\neq0$ |

No `decide` on the determinant: Laplace along the two constant rows as a ring identity, then
Laurent coefficient extraction at the strictly largest bracket argument.

## 📂 Contents

```
paper/     ghu_notch.tex/.pdf (EN), ghu_notch_es.tex/.pdf (ES), and the four data figures
lean/      NotchCentreCharge.lean
data/      the numbers behind every table and figure, as CSV
scripts/   every script that regenerates a number quoted in the paper
```

| script | regenerates |
|---|---|
| `_central_element.sage`, `_central_element2.sage` | the principle: the shift as a central element, up to Weyl |
| `_is_it_su4.sage` | the $SU(N)$ sweep above |
| `_derive_mq.sage`, `_mq_spread.sage` | §2 the $(m,q)$ projection derived from the boundary conditions |
| `notch_condition.sage`, `_fourier_probe.py` | §3 the notch predicate, and the empty Fourier sublattice |
| `notch_degenerate.sage`, `_fitbig.sage`, `_proofcheck.sage` | §4 the class $\mathcal{Z}$ and the proof's own identities |
| `_dichotomy.py`, `_period.py`, `_modular_probe.py` | §5 the inversion, the periodicity test, the surviving symmetry group |
| `_decoupling.py`, `_kconv.py` | §6 the exact decoupling and its convergence |
| `_vacconv.py`, `_basin*.py` | §7 the basin measure replacing the finite-difference Hessian |
| `notch_preflight.py` | §8 the tool |
| `make_notch_figs.py`, `make_fig_spectrum.py`, `make_release_data.py` | the figures and the CSVs |
| `hists.json` | the $(m,q)$ weight data every script reads |

## 📖 Citing

```bibtex
@misc{Marin2026CentreParity,
  author = {Mar\'in, Carles},
  title  = {A Centre-Charge Selection Rule for the Wilson-Line Potential:
            The Fundamental Domain of Gauge--Higgs Unification Is Representation-Dependent},
  year   = {2026},
  doi    = {10.5281/zenodo.21438227},
  note   = {Part III of a series}
}
```

## ⚖️ Honesty ledger

The centre-charge congruence is classical and we claim none of it. The mechanism behind the
vanishing theorem is classical too — $\{1,-1,t,t^{-1}\}$ has determinant $-1$, so it
parametrizes the improper component of $O(4)$, where an irrep and its associate
$\mu\otimes\det$ have opposite characters — and the determinant technique is described by
Ayyer and Behrend as routine. **What is ours:** the principle and its four readings, the
dichotomy, the explicit classification of $\mathcal{Z}$ with its Lean certificate, the
inversion of the fundamental domain with direct verification, the decoupling corollary, the
derivation of the $(m,q)$ projection from the boundary conditions, and the $SU(N)$ sweep.

Two classical sources remain unread and are flagged as such in the paper: King's 1971
modification-rules paper in its own exposition, and Littlewood's *Theory of Group Characters*
Ch. XI. Neither is load-bearing for the physics.

---

Carles Marín · `karlesmarin@gmail.com` · Claude (Anthropic) as AI research assistant
