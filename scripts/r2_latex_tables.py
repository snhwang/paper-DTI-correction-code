"""
R2: LaTeX for Tables 2 to 5, generated from the R2 number files.

Reads r2_numbers.csv (r2_tables.py) and r2_decomposition.csv
(r2_decomposition.py), so no table value is typed by hand.

    python r2_latex_tables.py

Writes r2_tables.tex, which is pasted into mri_revision.tex.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
N = pd.read_csv(HERE / "r2_numbers.csv")
DEC = pd.read_csv(HERE / "r2_decomposition.csv")
BS = "\\"

NAMES = {"classic": "Classic", "cross": "Refined", "v2_slab": "Measured axis",
         "pv_perp": BS + r"(\lambda_2/\lambda_3\)", "ALPS-PAS": "ALPS-PAS " + BS + r"cite{ref3}",
         "LD-ALPS": "LD-ALPS " + BS + r"cite{ref4}",
         "reoriented": "Template reorientation " + BS + r"cite{ref2}"}
ORDER = ["classic", "cross", "v2_slab", "pv_perp", "ALPS-PAS", "LD-ALPS", "reoriented"]


def val(q, cohort, index):
    s = N[(N.quantity == q) & (N.cohort == cohort) & (N["index"] == index)]
    return None if s.empty or pd.isna(s.value.iloc[0]) else float(s.value.iloc[0])


def f(v, d=3, sign=True):
    if v is None:
        return "--"
    s = f"{v:+.{d}f}" if sign else f"{v:.{d}f}"
    return "$" + s.replace("-", "-") + "$"


def p_fmt(p):
    if p is None:
        return "--"
    if p < 0.001:
        return r"$<0.001$"
    return f"${p:.3f}$" if p < 0.01 else f"${p:.2f}$"


def star(r, p):
    return f(r) + ("*" if p is not None and p < 0.05 else "")


def npose(cohort):
    s = N[(N.quantity == "median_total_rotation") & (N.cohort == cohort)]
    return int(s.n.iloc[0])


def table2():
    rows = [("Median total rotation", "median_total_rotation", None),
            ("Participants beyond 10" + BS + r"textdegree{} (DLBS) or 5" + BS + r"textdegree{} (HCP-A)", None, None)]
    out = [BS + r"begin{table}[pos=tb]",
           BS + r"caption{Head pose and its association with age, first session of each participant. DLBS keeps the head position of each session and HCP-A has it removed during preprocessing. Each axis enters as its absolute value.}",
           BS + r"label{tbl:pose-age}",
           BS + r"begin{tabular*}{" + BS + r"tblwidth}{@{}LCC@{}}",
           BS + "toprule",
           BS + r"textbf{Measure} & " + BS + f"textbf{{DLBS ($n={npose('DLBS')}$)}} & " + BS + f"textbf{{HCP-A ($n={npose('HCP-A')}$)}} " + BS + BS,
           BS + "midrule"]
    deg = "^{" + BS + "circ}"
    d_med, h_med = val("median_total_rotation", "DLBS", "pose"), val("median_total_rotation", "HCP-A", "pose")
    out.append(f"Median total rotation & ${d_med:.1f}{deg}$ & ${h_med:.1f}{deg}$ " + BS + BS)
    d_b, h_b = val("pct_beyond_10deg", "DLBS", "pose"), val("pct_beyond_5deg", "HCP-A", "pose")
    out.append(f"Total rotation beyond $10{deg}$ (DLBS) or $5{deg}$ (HCP-A) & ${d_b:.1f}" + BS
               + f"%$ & ${h_b:.1f}" + BS + "%$ " + BS + BS)

    def pp(p):
        return "<0.001" if p < 0.001 else (f"={p:.3f}" if p < 0.01 else f"={p:.2f}")

    for a, lab in (("pitch", "Pitch"), ("yaw", "Yaw"), ("roll", "Roll"), ("total", "Total rotation")):
        cells = [f"$r={val(f'age_r_{a}', c, 'pose'):+.3f}$, $p{pp(val(f'age_p_{a}', c, 'pose'))}$"
                 for c in ("DLBS", "HCP-A")]
        out.append(f"{lab} against age & {cells[0]} & {cells[1]} " + BS + BS)
    out += [BS + "bottomrule", BS + r"end{tabular*}", BS + r"end{table}"]
    return "\n".join(out)


def table3():
    q = lambda c, i, k: DEC[(DEC.cohort == c) & (DEC["index"] == i) & (DEC.quantity == k)].iloc[0]
    out = [BS + r"begin{table*}[pos=tb]",
           BS + r"caption{Division of each index's age decline into radial anisotropy and axis error, using Eq.~" + BS + r"eqref{eq:factor}. Slopes of the logarithm on age, in percent per decade, first session of each participant. The total is the sum of the $" + BS + r"lambda_2/" + BS + r"lambda_3$ term and the axis-error term. The pose-adjusted column adds absolute pitch and total rotation to the axis-error model. The share is the axis-error term divided by the total. Brackets are $95" + BS + r"%$ intervals from $2000$ participant bootstrap resamples.}",
           BS + r"label{tbl:decomposition}",
           BS + r"begin{tabular*}{" + BS + r"tblwidth}{@{}LCCCCC@{}}",
           BS + "toprule",
           BS + r"textbf{Index} & " + BS + r"textbf{Total} & " + BS + r"textbf{$" + BS + r"lambda_2/" + BS + r"lambda_3$ term} & "
           + BS + r"textbf{Axis-error term} & " + BS + r"textbf{Axis error, pose-adjusted} & " + BS + r"textbf{Axis-error share (" + BS + r"%)} " + BS + BS,
           BS + "midrule"]
    for c in ("DLBS", "HCP-A"):
        n = int(DEC[DEC.cohort == c].n.iloc[0])
        out.append(BS + r"multicolumn{6}{@{}l}{" + BS + r"textit{" + f"{c}, $n={n}$" + "}} " + BS + BS)
        for i in ("classic", "cross", "v2_slab"):
            t, r = q(c, i, "total"), q(c, i, "ratio")
            a, ap, s = q(c, i, "attenuation"), q(c, i, "attenuation_pose"), q(c, i, "axis_share_pct")
            ci = lambda x, d=2: f"${x.value:+.{d}f}$ [${x.lo:+.{d}f}$, ${x.hi:+.{d}f}$]"
            share = f"${s.value:.0f}$ [${s.lo:.0f}$, ${s.hi:.0f}$]"
            out.append(f"{NAMES[i]} & ${t.value:+.2f}$ & ${r.value:+.2f}$ & {ci(a)} & {ci(ap)} & {share} " + BS + BS)
    out += [BS + "bottomrule", BS + r"end{tabular*}", BS + r"end{table*}"]
    return "\n".join(out)


def table4():
    out = [BS + r"begin{table*}[pos=tb]",
           BS + r"caption{Each index against the bound $" + BS + r"bar" + BS + r"rho$ ($" + BS + r"lambda_2/" + BS + r"lambda_3$), first session of each participant. Attained is the median ratio of the index to $" + BS + r"lambda_2/" + BS + r"lambda_3$. The age columns give the Pearson correlation with age before and after $" + BS + r"lambda_2/" + BS + r"lambda_3$ is partialled out, with an asterisk marking $p<0.05$. LD-ALPS places its own regions, so its ratio to $" + BS + r"lambda_2/" + BS + r"lambda_3$ compares different voxels. Template reorientation applies only to DLBS, since HCP-A was already aligned.}",
           BS + r"label{tbl:bound}",
           BS + r"begin{tabular*}{" + BS + r"tblwidth}{@{}LCCCCCCCC@{}}",
           BS + "toprule",
           r" & " + BS + r"multicolumn{2}{c}{" + BS + r"textbf{Attained}} & " + BS + r"multicolumn{2}{c}{" + BS + r"textbf{$r$ with $" + BS + r"lambda_2/" + BS + r"lambda_3$}} & "
           + BS + r"multicolumn{2}{c}{" + BS + r"textbf{Age, raw}} & " + BS + r"multicolumn{2}{c}{" + BS + r"textbf{Age, partialled}} " + BS + BS,
           BS + r"textbf{Index} & HCP-A & DLBS & HCP-A & DLBS & HCP-A & DLBS & HCP-A & DLBS " + BS + BS,
           BS + "midrule"]
    for i in ORDER:
        cells = []
        for q in ("attainment", "r_with_ratio"):
            for c in ("HCP-A", "DLBS"):
                v = val(q, c, i)
                cells.append(f(v, 2, sign=False) if i != "pv_perp" else ("$1$" if v is not None else "--"))
        for c in ("HCP-A", "DLBS"):
            cells.append(f(val("age_r_raw", c, i)))
        for c in ("HCP-A", "DLBS"):
            cells.append("--" if i == "pv_perp" else star(val("age_r_partial_ratio", c, i),
                                                          val("age_p_partial_ratio", c, i)))
        if i == "ALPS-PAS":
            out.append(BS + "midrule")
        out.append(f"{NAMES[i]} & " + " & ".join(cells) + " " + BS + BS)
    out += [BS + "bottomrule", BS + r"end{tabular*}", BS + r"end{table*}"]
    return "\n".join(out)


def table5():
    out = [BS + r"begin{table}[pos=tb]",
           BS + r"caption{Reliability between visits, ICC(1,1), over participants with repeat visits ($628$ in HCP-A and $156$ in DLBS).}",
           BS + r"label{tbl:reliability}",
           BS + r"begin{tabular*}{" + BS + r"tblwidth}{@{}LCC@{}}",
           BS + "toprule",
           BS + r"textbf{Index} & " + BS + r"textbf{HCP-A} & " + BS + r"textbf{DLBS} " + BS + BS,
           BS + "midrule"]
    for i in ORDER:
        if i == "ALPS-PAS":
            out.append(BS + "midrule")
        out.append(f"{NAMES[i]} & {f(val('icc', 'HCP-A', i), 3, False)} & {f(val('icc', 'DLBS', i), 3, False)} " + BS + BS)
    out += [BS + "bottomrule", BS + r"end{tabular*}", BS + r"end{table}"]
    return "\n".join(out)


if __name__ == "__main__":
    tex = "\n\n".join([table2(), table3(), table4(), table5()]) + "\n"
    (HERE / "r2_tables.tex").write_text(tex, encoding="utf-8", newline="\n")
    print(tex)
