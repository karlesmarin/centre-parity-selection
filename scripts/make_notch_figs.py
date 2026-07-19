# make_notch_figs.py - publication figures for Part III.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
# Palette: validated categorical/sequential/diverging set (all six checks PASS, light surface).
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["hatch.linewidth"]=0.5
matplotlib.rcParams["hatch.color"]="#ffffff"
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Patch

SURF="#fcfcfb"; INK="#0b0b0b"; SEC="#52514e"; MUT="#898781"; GRID="#e1e0d9"; BASE="#c3c2b7"
BLUE="#2a78d6"; ORANGE="#eb6834"; VIOLET="#4a3aa7"; RED="#e34948"; GREEN="#008300"
SEQ=["#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7","#3987e5","#2a78d6",
     "#256abf","#1c5cab","#184f95","#104281","#0d366b"]
cmap_seq=LinearSegmentedColormap.from_list("seqblue",SEQ)
cmap_div=LinearSegmentedColormap.from_list("div",[ "#0d366b","#2a78d6","#9ec5f4","#f0efec","#f3a9a8","#e34948","#8f2020"])
plt.rcParams.update({"font.family":"sans-serif","font.size":8.5,"axes.edgecolor":BASE,
    "axes.labelcolor":SEC,"xtick.color":MUT,"ytick.color":MUT,"text.color":INK,
    "axes.facecolor":SURF,"figure.facecolor":SURF,"savefig.facecolor":SURF,
    "axes.linewidth":0.7,"xtick.major.width":0.7,"ytick.major.width":0.7})

H=json.load(open("hists.json"))
def parse(h):
    d={}
    for k,v in h.items():
        m,q=map(int,k.split(',')); d[(m,q)]=d.get((m,q),0)+int(v)
    return d
GAUGE=[(1,0,2),(-1,0,2),(1,1,2),(-1,1,2),(2,0,1),(0,0,1),(-2,0,1)]
def delta(d):
    ms=sorted(set(m for (m,q) in d))
    return ms,[d.get((m,0),0)-d.get((m,1),0) for m in ms]

# ---------------- FIG 1: the complementary combs ----------------
fig,axes=plt.subplots(1,2,figsize=(7.0,2.55),sharey=False)
for ax,(name,lab,sub) in zip(axes,[("35",r"$\mathbf{35}=(4,0,0)$",r"$|\lambda|=4$ even $\Rightarrow$ pairs at odd $m$"),
                                   ("60",r"$\mathbf{60}=(0,2,1)$",r"$|\lambda|=7$ odd $\Rightarrow$ pairs at even $m$")]):
    ms,dl=delta(parse(H[name]))
    cols=[BLUE if m%2==0 else ORANGE for m in ms]
    ax.bar(ms,dl,width=0.62,color=cols,zorder=3)
    ax.axhline(0,color=BASE,lw=0.9,zorder=2)
    ax.set_xlabel(r"Kaluza--Klein twist level $m$")
    ax.set_title(lab+"\n"+sub,fontsize=8.5,color=INK,pad=7,linespacing=1.5)
    ax.grid(axis="y",color=GRID,lw=0.6,zorder=0); ax.set_axisbelow(True)
    ax.set_xticks(range(min(ms),max(ms)+1,2))
    for s in ("top","right"): ax.spines[s].set_visible(False)
    # annotate the zero comb
    zer=[m for m,v in zip(ms,dl) if v==0 and (m%2==0)==(name!="35")]
    if zer:
        ax.plot(zer,[0]*len(zer),marker="o",ms=3.4,ls="none",
                mfc=SURF,mec=BLUE if name!="35" else ORANGE,mew=1.1,zorder=5)
axes[0].set_ylabel(r"$\delta(m)=\mathrm{mult}(m,0)-\mathrm{mult}(m,1)$")
leg=[Patch(facecolor=BLUE,label="even $m$"),Patch(facecolor=ORANGE,label="odd $m$")]
axes[1].legend(handles=leg,frameon=False,fontsize=8,loc="upper right",labelcolor=SEC)
fig.tight_layout(); fig.savefig("paper/fig_combs.pdf",bbox_inches="tight"); plt.close(fig)
print("fig_combs.pdf")

# ---------------- FIG 2: the dichotomy landscape ----------------
inZ =lambda a,b,c: b%2==1 and ((a%2==1 and c%2==1) or a==c)
adm =lambda a,b,c: (a+2*b+3*c)%2==1 and b>=1 and a+b+c>=3
MAXB,MAXAC=9,11
fig,ax=plt.subplots(figsize=(6.0,3.4))
for b in range(0,MAXB+1):
    for ac in range(0,MAXAC+1):
        # representative (a,c) with a+c = ac; colour by class (parity is a function of a+c)
        reps=[(a,ac-a) for a in range(ac+1)]
        deg=any(inZ(a,b,c) for a,c in reps); par=(ac)%2
        col = VIOLET if deg else (BLUE if par==1 else ORANGE)
        hat = "xxx" if deg else ("///" if par==1 else None)   # redundant channel: texture
        ax.add_patch(plt.Rectangle((ac-0.42,b-0.42),0.84,0.84,facecolor=col,
                                   edgecolor=SURF,lw=1.2,zorder=2,alpha=0.92,
                                   hatch=hat))
        if any(adm(a,b,c) for a,c in reps):
            ax.plot(ac,b,marker="o",ms=5.0,mfc=SURF,mec=INK,mew=1.4,zorder=4)
ax.set_xlim(-0.8,MAXAC+0.8); ax.set_ylim(-0.8,MAXB+0.8)
ax.set_xlabel(r"$a+c$"); ax.set_ylabel(r"$b$")
ax.set_xticks(range(0,MAXAC+1)); ax.set_yticks(range(0,MAXB+1))
for s in ("top","right"): ax.spines[s].set_visible(False)
ax.set_aspect("equal")
leg=[Patch(facecolor=BLUE,hatch="///",label=r"even-$m$ notch  ($|\lambda|$ odd)"),
     Patch(facecolor=ORANGE,label=r"odd-$m$ pairing  ($|\lambda|$ even)"),
     Patch(facecolor=VIOLET,hatch="xxx",label=r"degenerate $\mathcal{Z}$  (both)"),
     plt.Line2D([],[],marker="o",ls="none",ms=5,mfc=SURF,mec=INK,mew=1.3,
                label="hosts a quark cell")]
ax.legend(handles=leg,frameon=False,fontsize=8,loc="upper left",
          bbox_to_anchor=(1.01,1.0),labelcolor=SEC,handlelength=1.2)
fig.tight_layout(); fig.savefig("paper/fig_landscape_notch.pdf",bbox_inches="tight"); plt.close(fig)
print("fig_landscape_notch.pdf")

# ---------------- FIG 3: the fundamental domain (the central claim) ----------------
N=384
def Vgrid(hf,KMAX=10,wg=0.35):
    A=np.zeros((N,N))
    KS=[(x,y) for x in range(-KMAX,KMAX+1) for y in range(-KMAX,KMAX+1) if (x,y)!=(0,0)]
    for hist,sgn in ((hf,1.0),(GAUGE,-wg)):
        for (m,q,mu) in hist:
            for (k1,k2) in KS:
                A[(k1*m)%N,(k2*m)%N]+=sgn*mu/(k1*k1+k2*k2)**3*((-1)**((k2*q)%2))
    return np.real(np.fft.ifft2(A))*N*N
FERM={k:[(int(a.split(',')[0]),int(a.split(',')[1]),m) for a,m in v.items()] for k,v in H.items()}
fig,axes=plt.subplots(2,2,figsize=(7.0,6.1),gridspec_kw={"height_ratios":[1,1],"hspace":0.42,"wspace":0.18})
ext=[0,2,0,2]
for j,(name,lab) in enumerate([("35",r"$\mathbf{35}=(4,0,0)$  —  AHMN's fermion, $|\lambda|$ even"),
                               ("60",r"$\mathbf{60}=(0,2,1)$  —  admissible, $|\lambda|$ odd")]):
    V=Vgrid(FERM[name]); Vn=(V-V.min())/(V.max()-V.min())
    ax=axes[0,j]
    ax.imshow(Vn.T,origin="lower",extent=ext,cmap=cmap_seq,aspect="auto",interpolation="bilinear")
    ax.contour(np.linspace(0,2,N),np.linspace(0,2,N),Vn.T,levels=7,colors=SURF,linewidths=0.45,alpha=0.55)
    ax.axhline(0.5,color=RED,lw=1.4,ls="--",zorder=5)
    ax.axhline(1.0,color=INK,lw=1.0,ls=":",zorder=5)
    ax.text(0.04,0.545,r"AHMN cutoff $\alpha_2=1/2$",color=RED,fontsize=7.0,ha="left",va="bottom")
    ax.text(0.04,1.045,r"$\alpha_2=1$",color=INK,fontsize=7.0,ha="left",va="bottom")
    ax.set_title(lab,fontsize=8.5,color=INK,pad=6)
    ax.set_xticklabels([])
    if j==0: ax.set_ylabel(r"$\alpha_2$")
    # difference panel
    D=V-np.roll(V,N//2,axis=1); s=np.max(np.abs(V-V.mean()))
    ax2=axes[1,j]
    lim=max(np.abs(D).max(),1e-12)
    im=ax2.imshow((D/s).T,origin="lower",extent=ext,cmap=cmap_div,aspect="auto",
                  vmin=-0.32,vmax=0.32,interpolation="bilinear")
    ax2.set_xlabel(r"$\alpha_1$")
    if j==0: ax2.set_ylabel(r"$\alpha_2$")
    rel=np.abs(D).max()/s
    ax2.set_title((r"$V(\alpha_2{+}1)-V(\alpha_2)$   " +
                   (r"$\equiv 0$  (period 1 holds)" if rel<1e-10 else
                    r"$\max = %.0f\%%$ of scale"%(100*rel))),
                  fontsize=8.5,color=(GREEN if rel<1e-10 else RED),pad=6)
    for a in (axes[0,j],axes[1,j]):
        for sp in a.spines.values(): sp.set_color(BASE)
cb=fig.colorbar(im,ax=axes[1,:].tolist(),fraction=0.026,pad=0.025)
cb.set_label("difference / potential scale",color=SEC,fontsize=7.5)
cb.outline.set_edgecolor(BASE); cb.ax.tick_params(labelsize=7,color=MUT)
fig.savefig("paper/fig_domain.pdf",bbox_inches="tight"); plt.close(fig)
print("fig_domain.pdf")
