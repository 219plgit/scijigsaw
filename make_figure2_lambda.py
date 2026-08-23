"""Figure 2 (lambda sweep): two panels, direct labels, plain-language axes.
A: concentration (N_eff vs lambda). B: robustness of the biological conclusion."""
import math, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
exec(open('run_eif3_centrepiece.py').read()
     .split('# ------------------------------------------------------------------ report')[0]
     .split('print("LEMMA')[0])
from sj_prob import TypedModel, merger_F

LAMS = [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4]
SETS = [("declared (literature)", DECLARED_C, "-",  "o"),
        ("AlphaFold 10/10",       AF3_10,     "--", "s"),
        (r"AlphaFold $\geq$8/10", AF3_8,      "-.", "^")]
H = frozenset({"Tif35", "Tif34"})

fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.6, 3.9))

for label, C, ls, mk in SETS:
    m = TypedModel(V=V, C=C, P=PREREQ)
    n0 = merger_F(m)(frozenset(V))
    rows = {l: weighted_trees(m, l) for l in LAMS}
    ys = [h2_neff([p for _, _, p in rows[l]])[1] for l in LAMS]
    ax.plot(LAMS, ys, ls, marker=mk, color="black", ms=4, mfc="white", lw=1.2)
    ss = [support_weighted(rows[l], H) for l in LAMS]
    bx.plot(LAMS, ss, ls, marker=mk, color="black", ms=4, mfc="white", lw=1.2)

# Panel A dressing
ax.set_xlabel(r"weight on structural evidence, $\lambda$")
ax.set_ylabel("effective number of\ncompeting mechanisms, $N_{\\mathrm{eff}}$")
ax.set_ylim(0, 33); ax.set_xlim(-0.15, 4.2)
ax.annotate("flat: acyclic graph (provable)",
            xy=(2.6, 14.0), xytext=(1.9, 7.5), fontsize=8.5, style="italic",
            arrowprops=dict(arrowstyle="-", lw=0.7))
ax.text(0.55, 12.6, "declared (literature), 14 trees", fontsize=8.5)
ax.text(1.62, 24.3, "AlphaFold 10/10, 27 trees", fontsize=8.5, rotation=-13)
ax.text(1.32, 29.2, r"AlphaFold $\geq$8/10, 30 trees", fontsize=8.5, rotation=-11)
ax.set_title("A   Narrowing the field", loc="left", fontsize=10)
ax.grid(alpha=0.25, lw=0.5)

# Panel B dressing
bx.axhspan(0.357, 0.500, color="0.88", zorder=0)
bx.set_xlabel(r"weight on structural evidence, $\lambda$")
bx.set_ylabel("probability that Tif35--Tif34\nforms as a module")
bx.set_ylim(0, 0.65); bx.set_xlim(-0.15, 4.2)
bx.annotate("band: 0.36–0.50",
            xy=(2.1, 0.43), xytext=(1.6, 0.12), fontsize=8.5, style="italic",
            arrowprops=dict(arrowstyle="-", lw=0.7))
bx.set_title("B   Robustness of the conclusion", loc="left", fontsize=10)
bx.grid(alpha=0.25, lw=0.5)

fig.tight_layout()
fig.savefig("figure2_lambda_sweep.pdf")
fig.savefig("figure2_lambda_sweep.png", dpi=200)
print("saved")
