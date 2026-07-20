# make_fig_spectrum.py - the calling-card figure: the notch as an empty sublattice in Fourier space.
#   Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)
# Every mode the potential COULD carry is drawn as a ghost seat; the ones it actually carries are
# filled. Eight representations leave the same seats empty. One does not.
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
SURF="#fcfcfb"; INK="#0b0b0b"; SEC="#52514e"; MUT="#898781"; GRID="#e1e0d9"; BASE="#c3c2b7"
BLUE="#2a78d6"; ORANGE="#eb6834"; RED="#e34948"; GHOST="#d8d6cf"
plt.rcParams.update({"font.family":"sans-serif","font.size":8,"axes.edgecolor":BASE,
 "axes.labelcolor":SEC,"xtick.color":MUT,"ytick.color":MUT,"text.color":INK,
 "axes.facecolor":SURF,"figure.facecolor":SURF,"savefig.facecolor":SURF,"axes.linewidth":0.6})
H=json.load(open("hists.json"))
F={k:[(int(a.split(',')[0]),int(a.split(',')[1]),m) for a,m in v.items()] for k,v in H.items()}
LAB={"35":"(4,0,0)", "60":"(0,2,1)", "84":"(0,1,3)", "140a":"(1,1,2)", "140b":"(0,3,1)", "224":"(0,2,3)", "280":"(0,4,1)", "360":"(1,2,2)", "756":"(1,3,2)"}
order = ["35", "60", "84", "140a", "140b", "224", "280", "360", "756"]
KM=5; LIM=13
def modes(hf):
    """amplitude at (n1,n2) coming from k2-ODD shells only -- the sector the notch controls."""
    A={}
    for (m,q,mu) in hf:
        for k1 in range(-KM,KM+1):
            for k2 in range(-KM,KM+1):
                if (k1,k2)==(0,0) or k2%2==0: continue
                n1,n2=k1*m,k2*m
                if abs(n1)>LIM or abs(n2)>LIM: continue
                A[(n1,n2,m%2)]=A.get((n1,n2,m%2),0.0)+mu/(k1*k1+k2*k2)**3*((-1)**((k2*q)%2))
    return A
fig,axes=plt.subplots(3,3,figsize=(7.2,7.4))
for ax,r in zip(axes.ravel(),order):
    A=modes(F[r]); adm = r!="35"
    seats=set((n1,n2) for (n1,n2,p) in A if p==0)          # m even -> notch-protected seats
    filled=[(n1,n2,abs(v)) for (n1,n2,p),v in A.items() if p==0 and abs(v)>1e-12]
    other =[(n1,n2,abs(v)) for (n1,n2,p),v in A.items() if p==1 and abs(v)>1e-12]
    mx=max([v for *_,v in other+filled]+[1e-9])
    # ghost seats: positions a protected mode could occupy
    if seats:
        gx=[s[0] for s in seats]; gy=[s[1] for s in seats]
        ax.scatter(gx,gy,s=26,facecolors="none",edgecolors=GHOST,linewidths=0.9,zorder=2)
    if other:
        ax.scatter([o[0] for o in other],[o[1] for o in other],
                   s=[16+150*(o[2]/mx) for o in other],color=BLUE,alpha=0.85,
                   edgecolors=SURF,linewidths=0.6,zorder=3)
    if filled:
        ax.scatter([f[0] for f in filled],[f[1] for f in filled],
                   s=[16+150*(f[2]/mx) for f in filled],color=RED,alpha=0.92,
                   edgecolors=SURF,linewidths=0.6,zorder=4)
    ax.set_xlim(-LIM-1,LIM+1); ax.set_ylim(-LIM-1,LIM+1); ax.set_aspect("equal")
    ax.axhline(0,color=GRID,lw=0.6,zorder=1); ax.axvline(0,color=GRID,lw=0.6,zorder=1)
    ax.set_xticks([-10,0,10]); ax.set_yticks([-10,0,10])
    for s in ("top","right"): ax.spines[s].set_visible(False)
    n_occ=len(filled)
    ax.set_title(r"$\mathbf{%s}$  %s"%(r.replace("a","a").replace("b","b"),LAB[r]),
                 fontsize=8.5,color=INK,pad=3)
    ax.text(0.5,-0.30,("seats empty" if n_occ==0 else "%d seats FILLED"%n_occ),
            transform=ax.transAxes,ha="center",fontsize=8,
            color=(SEC if n_occ==0 else RED),
            fontweight=("normal" if n_occ==0 else "bold"))
for ax in axes[:, 0]: ax.set_ylabel(r"$n_2$")
for ax in axes[-1, :]: ax.set_xlabel(r"$n_1$")
leg=[plt.Line2D([],[],marker="o",ls="none",ms=6,color=BLUE,label=r"modes present ($m$ odd)"),
     plt.Line2D([],[],marker="o",ls="none",ms=6,color=RED,label=r"notch-protected modes, present ($m$ even)"),
     plt.Line2D([],[],marker="o",ls="none",ms=6,mfc="none",mec=GHOST,mew=1.2,
                label="seat a protected mode could occupy")]
fig.legend(handles=leg,frameon=False,fontsize=8.2,loc="lower center",ncol=3,
           labelcolor=SEC,bbox_to_anchor=(0.5,-0.005))
fig.tight_layout(rect=[0,0.045,1,1]); fig.subplots_adjust(hspace=0.52)
fig.savefig("paper/fig_spectrum.pdf",bbox_inches="tight"); plt.close(fig)
print("fig_spectrum.pdf")
