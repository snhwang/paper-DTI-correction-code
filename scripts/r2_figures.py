"""
R2 figures 2 and 3, from the R2 number files.

Figure 2  head pose and age (a, b), and the division of the age decline into
          radial anisotropy and axis error (c), from r2_decomposition.csv
Figure 3  sensitivity to genuine repositioning between visits, from
          r2_repositioning_pairs.csv and r2_repositioning.csv

Palette: the dataviz reference categorical slots 1 to 4, validated together
(light surface, all checks pass). Slots 3 and 4 fall below 3:1 contrast against
the page, so they appear only where the index is also named in text on the axis.

    python r2_figures.py

Writes fig_r2_posture.{png,pdf} and fig_r2_repositioning.{png,pdf} beside the
manuscript.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
OUT = HERE.parent

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e2e1dc"
SERIES = {"Classic": "#2a78d6", "Refined": "#eb6834",
          "Measured axis": "#1baf7a", "lambda2/lambda3": "#eda100"}
OBLIQUE, ALIGNED = "#3b3b3b", "#b4b3ad"
TISSUE, AXIS = "#a9a8a2", "#4a3aa7"

plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 9,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.edgecolor": "#8a8984", "axes.linewidth": 0.7, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "figure.dpi": 300, "lines.linewidth": 1.6, "axes.spines.top": False,
    "axes.spines.right": False,
})


def key(d):
    d = d.copy()
    d["Subject_ID"] = d.Subject_ID.astype(str)
    d["Visit"] = d.Visit.astype(str)
    return d


def first_with_pose(mf, hf):
    h = key(pd.read_csv(HERE / hf))
    h = h[h.total <= 45.0]  # failed registrations (Methods 3.3)
    d = key(pd.read_csv(HERE / mf)).merge(h, on=["Subject_ID", "Visit"])
    d = d.sort_values(["Subject_ID", "Visit"]).groupby("Subject_ID").first().reset_index()
    return d.dropna(subset=["Age", "pitch", "classic"])


def grid(ax, axis="both"):
    ax.grid(axis=axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def figure2():
    D = first_with_pose("measured_pvs_axis_dlbs_all.csv", "head_rotation_dlbs.csv")
    H = first_with_pose("measured_pvs_axis_hcpa_b1500_all.csv", "head_rotation_hcpa.csv")
    dec = pd.read_csv(HERE / "r2_decomposition.csv")

    fig = plt.figure(figsize=(7.2, 2.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 0.7, 1.45], wspace=0.62,
                          left=0.065, right=0.99, bottom=0.17, top=0.86)

    ax = fig.add_subplot(gs[0, 0])
    for d, c, lab, y in ((H, ALIGNED, "HCP-A, aligned", 0.86), (D, OBLIQUE, "DLBS, oblique", 0.96)):
        x, v = d.Age.to_numpy(float), d.pitch.abs().to_numpy(float)
        ax.scatter(x, v, s=5, color=c, alpha=0.45, linewidths=0)
        lr = stats.linregress(x, v)
        xx = np.linspace(x.min(), x.max(), 50)
        ax.plot(xx, lr.intercept + lr.slope * xx, color=c, linewidth=2)
        ax.text(0.03, y, f"{lab}  r = {np.corrcoef(x, v)[0, 1]:+.2f}", transform=ax.transAxes,
                va="top", fontsize=7.5, color=INK if c == OBLIQUE else INK2)
    ax.set_xlabel("Age (years)")
    ax.set_ylabel("Absolute head pitch (degrees)")
    ax.set_ylim(0, 27)
    ax.set_title("(a) Pitch rises with age", loc="left")
    grid(ax)

    ax = fig.add_subplot(gs[0, 1])
    rs = [np.corrcoef(D.Age, D[c].abs())[0, 1] for c in ("pitch", "roll", "yaw")]
    ax.bar(range(3), rs, 0.62, color=[SERIES["Classic"], ALIGNED, ALIGNED])
    for i, r in enumerate(rs):
        ax.text(i, r + 0.012, f"{r:+.2f}", ha="center", va="bottom", fontsize=7.5, color=INK)
    ax.axhline(0, color="#8a8984", linewidth=0.7)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["Pitch", "Roll", "Yaw"])
    ax.set_ylabel("Correlation with age, DLBS")
    ax.set_ylim(-0.05, 0.42)
    ax.set_title("(b) By axis", loc="left")
    grid(ax, "y")

    ax = fig.add_subplot(gs[0, 2])
    rows = [("DLBS", "classic", "DLBS classic"), ("DLBS", "cross", "DLBS refined"),
            ("HCP-A", "classic", "HCP-A classic"), ("HCP-A", "cross", "HCP-A refined")]
    q = lambda c, i, k: float(dec[(dec.cohort == c) & (dec["index"] == i) & (dec.quantity == k)].value.iloc[0])
    gap = 0.03
    for y, (c, i, lab) in enumerate(rows[::-1]):
        r, a = -q(c, i, "ratio"), -q(c, i, "attenuation")
        share = q(c, i, "axis_share_pct")
        ax.barh(y, r, 0.6, color=TISSUE)
        if a > 0:
            ax.barh(y, a - gap, 0.6, left=r + gap, color=AXIS)
        ax.text(r + max(a, 0) + 0.08, y, f"{share:.0f}%", va="center", fontsize=7.5, color=INK)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[2] for r in rows[::-1]])
    ax.set_xlabel("Age decline of the index (% per decade)")
    longest = max(-q(c, i, "ratio") + max(-q(c, i, "attenuation"), 0) for c, i, _ in rows)
    ax.set_xlim(0, longest + 0.75)
    ax.tick_params(axis="y", pad=6)
    ax.set_title("(c) Age decline, tissue and axis error", loc="left")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=TISSUE, label=r"$\lambda_2/\lambda_3$ term"),
                       Patch(color=AXIS, label="Axis-error term, share labeled")],
              loc="upper center", bbox_to_anchor=(0.42, -0.2), ncol=2, frameon=False,
              handlelength=1.2, columnspacing=1.2)
    grid(ax, "x")

    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig_r2_posture.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)


def figure3():
    p = pd.read_csv(HERE / "r2_repositioning_pairs.csv")
    s = pd.read_csv(HERE / "r2_repositioning.csv").set_index("index")
    fig = plt.figure(figsize=(7.2, 2.7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.5,
                          left=0.07, right=0.98, bottom=0.17, top=0.88)

    ax = fig.add_subplot(gs[0, 0])
    xx = np.linspace(0, p.d_theta_SCR.quantile(0.98), 50)
    for name in ("Classic", "Refined"):
        c = SERIES[name]
        ax.scatter(p.d_theta_SCR, p[name], s=6, color=c, alpha=0.35, linewidths=0)
        b, a0 = np.polyfit(p.d_theta_SCR, p[name], 1)
        ax.plot(xx, a0 + b * xx, color=c, linewidth=2,
                label=f"{name}, {b:+.2f}% per degree")
    ax.set_xlim(0, xx.max())
    ax.set_ylim(0, p[["Classic", "Refined"]].quantile(0.98).max())
    ax.set_xlabel(r"Change in projection-tract angle between visits (degrees)")
    ax.set_ylabel("Change in index between visits (%)")
    ax.set_title("(a) Classic and refined, every visit pair", loc="left")
    ax.legend(loc="upper left", frameon=False)
    grid(ax)

    ax = fig.add_subplot(gs[0, 1])
    names = ["Classic", "Refined", "Measured axis", "lambda2/lambda3"]
    labels = ["Classic", "Refined", "Measured axis", r"$\lambda_2/\lambda_3$"]
    for y, name in enumerate(names[::-1]):
        r = s.loc[name]
        ax.plot([r.ci_lo, r.ci_hi], [y, y], color=SERIES[name], linewidth=2, solid_capstyle="round")
        ax.plot(r.slope_pct_per_deg, y, "o", ms=6, color=SERIES[name],
                markeredgecolor="white", markeredgewidth=1.2)
    ax.axvline(0, color="#8a8984", linewidth=0.8, linestyle="--")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(labels[::-1])
    ax.set_xlabel("Change per degree (%)")
    ax.set_title("(b) Slope for each index", loc="left")
    grid(ax, "x")

    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"fig_r2_repositioning.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    figure2()
    figure3()
    print(f"wrote {OUT / 'fig_r2_posture.png'} and {OUT / 'fig_r2_repositioning.png'}")
