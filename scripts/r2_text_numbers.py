"""
R2: the numbers quoted in the text that are not in a table.

Same conventions as r2_tables.py: warped atlas regions, first session of each
participant for cross-sectional values. The R1 body-mass split averaged pose
over each participant's visits. Here it uses the first visit, like everything
else.

    slice prescription   slab angulation from the raw header, against age
    body mass split      head pitch in DLBS participants above and below median BMI
    Williams' tests      classic against each corrected index, age association
    hemispheres          left against right, and classic against refined per side
    FA floor             largest change in each quantity across floors 0.15 to 0.25
    tract departure      median angle of each tract from its assumed scanner axis

    python r2_text_numbers.py

Writes r2_text_numbers.csv (quantity, cohort, index, value, n).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import atomic_io  # noqa: F401  writes become atomic on import

HERE = Path(__file__).resolve().parent
DIFF = HERE.parent.parent / "diffusion"

# The DLBS main sample is every participant (284), matching HCP-A. Cross-sectional
# values use each participant's first session. Reliability and repositioning use the
# 156 with repeat visits, selected inside each analysis. R2_DLBS=repeat restores the
# R1 restriction to those 156 for every analysis, and suffixes the outputs.
import os as _os
DLBS_MAIN = ("measured_pvs_axis_dlbs.csv" if _os.environ.get("R2_DLBS") == "repeat"
             else "measured_pvs_axis_dlbs_all.csv")
SUFFIX = "_dlbs_repeat" if _os.environ.get("R2_DLBS") == "repeat" else ""
MAX_ROTATION = 45.0   # above this a registration is treated as failed (Methods 3.3)
rows: list[dict] = []


def put(quantity, cohort, index, value, n=None):
    rows.append({"quantity": quantity, "cohort": cohort, "index": index,
                 "value": float(value), "n": n})
    print(f"  {cohort:6s} {index:14s} {quantity:34s} {value:+.4g}  (n {n})")


def key(d):
    d = d.copy()
    d["Subject_ID"] = d.Subject_ID.astype(str)
    d["Visit"] = d.Visit.astype(str)
    return d


def first(d):
    return d.sort_values(["Subject_ID", "Visit"]).groupby("Subject_ID").first().reset_index()


def williams(r12, r13, r23, n):
    """Williams' test for r12 against r13, which share variable 1. Two-sided p."""
    det = 1 - r12**2 - r13**2 - r23**2 + 2 * r12 * r13 * r23
    rbar = (r12 + r13) / 2
    t = (r12 - r13) * np.sqrt((n - 1) * (1 + r23)
                              / (2 * (n - 1) / (n - 3) * det + rbar**2 * (1 - r23) ** 3))
    return float(2 * stats.t.sf(abs(t), n - 3))


def main() -> None:
    main_files = {"HCP-A": "measured_pvs_axis_hcpa_b1500_all.csv",
                  "DLBS": DLBS_MAIN}
    D = {c: first(key(pd.read_csv(HERE / f))) for c, f in main_files.items()}

    print("slice prescription, DLBS")
    s = key(pd.read_csv(HERE / "slab_prescription_dlbs.csv"))
    s = first(s.merge(D["DLBS"][["Subject_ID", "Visit", "Age"]], on=["Subject_ID", "Visit"]))
    put("median_slab_angulation_deg", "DLBS", "prescription", s.slab_tilt.median(), len(s))
    r, p = stats.pearsonr(s.Age, s.slab_pitch)
    put("age_r_slab_pitch", "DLBS", "prescription", r, len(s))
    put("age_p_slab_pitch", "DLBS", "prescription", p, len(s))

    print("body mass split, DLBS")
    t = pd.read_csv(DIFF / "DLBS" / "ds004856_participants.tsv", sep="\t", low_memory=False)
    t["Subject_ID"] = t.participant_id.astype(str)
    t["BMI_W1"] = pd.to_numeric(t.BMI_W1, errors="coerce")
    h = key(pd.read_csv(HERE / "head_rotation_dlbs.csv"))
    h = h[h.total <= MAX_ROTATION]
    b = (D["DLBS"][["Subject_ID", "Visit", "Age"]]
         .merge(h, on=["Subject_ID", "Visit"]).merge(t[["Subject_ID", "BMI_W1"]], on="Subject_ID")
         .dropna(subset=["BMI_W1", "pitch"]))
    hi = b.BMI_W1 > b.BMI_W1.median()
    g1, g0 = b.pitch.abs()[hi], b.pitch.abs()[~hi]
    put("median_abs_pitch_high_bmi", "DLBS", "pose", g1.median(), int(hi.sum()))
    put("median_abs_pitch_low_bmi", "DLBS", "pose", g0.median(), int((~hi).sum()))
    put("welch_p_pitch_by_bmi", "DLBS", "pose", stats.ttest_ind(g1, g0, equal_var=False).pvalue, len(b))
    put("age_r_bmi", "DLBS", "pose", stats.pearsonr(b.Age, b.BMI_W1)[0], len(b))

    print("lambda2/lambda3 and the classic attenuation against head pitch, first session")
    for c, hf in (("DLBS", "head_rotation_dlbs.csv"), ("HCP-A", "head_rotation_hcpa.csv")):
        hh = key(pd.read_csv(HERE / hf))
        d = D[c].merge(hh[hh.total <= MAX_ROTATION], on=["Subject_ID", "Visit"]).dropna(
            subset=["Age", "pv_perp", "classic", "pitch"])
        ap = d.pitch.abs().to_numpy(float)
        r, p = stats.pearsonr(d.pv_perp, ap)
        put("r_ratio_abs_pitch", c, "pv_perp", r, len(d))
        put("p_ratio_abs_pitch", c, "pv_perp", p, len(d))
        X = np.column_stack([np.ones(len(d)), d.Age.to_numpy(float)])
        res = lambda v: v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
        rp = float(np.corrcoef(res(d.pv_perp.to_numpy(float)), res(ap))[0, 1])
        put("r_ratio_abs_pitch_given_age", c, "pv_perp", rp, len(d))
        r, p = stats.pearsonr(d.classic / d.pv_perp, ap)
        put("r_attenuation_abs_pitch", c, "classic", r, len(d))
        put("p_attenuation_abs_pitch", c, "classic", p, len(d))

    print("Williams' tests, age association, first session")
    for c, d in D.items():
        for other in ("cross", "v2_slab", "pv_perp"):
            s = d.dropna(subset=["Age", "classic", other])
            r12 = np.corrcoef(s.Age, s.classic)[0, 1]
            r13 = np.corrcoef(s.Age, s[other])[0, 1]
            r23 = np.corrcoef(s.classic, s[other])[0, 1]
            put("williams_p_vs_classic", c, other, williams(r12, r13, r23, len(s)), len(s))

    print("hemispheres, first session")
    for c, f in (("HCP-A", "hemisphere_age_hcpa_b1500.csv"), ("DLBS", "hemisphere_age_dlbs.csv")):
        # restricted to the main sample, so DLBS uses the same 156 participants as every
        # other cross-sectional value rather than all 284 in the hemisphere file
        d = key(pd.read_csv(HERE / f))
        d = first(d[d.Subject_ID.isin(D[c].Subject_ID)]).dropna()
        for v in ("classic", "refined"):
            put("left_minus_right_pct", c, v, 100 * ((d[f"{v}_L"] - d[f"{v}_R"]) / d[f"{v}_R"]).mean(), len(d))
            put("paired_p_left_right", c, v, stats.ttest_rel(d[f"{v}_L"], d[f"{v}_R"]).pvalue, len(d))
        for side in ("L", "R"):
            r12 = np.corrcoef(d.Age, d[f"classic_{side}"])[0, 1]
            r13 = np.corrcoef(d.Age, d[f"refined_{side}"])[0, 1]
            r23 = np.corrcoef(d[f"classic_{side}"], d[f"refined_{side}"])[0, 1]
            put(f"age_r_classic_{side}", c, "classic", r12, len(d))
            put(f"age_r_refined_{side}", c, "refined", r13, len(d))
            put(f"williams_p_{side}", c, "classic_vs_refined", williams(r12, r13, r23, len(d)), len(d))

    print("FA floor, largest change across 0.15 to 0.25")
    f = pd.read_csv(HERE / "fa_threshold_sweep.csv")
    for c in ("hcpa", "dlbs"):
        for m in ("attainment", "icc", "age", "beyond"):
            s = f[(f.cohort == c) & (f.measure == m)]
            span = s.groupby("variant").value.agg(lambda v: v.max() - v.min()).max()
            put(f"fa_floor_max_change_{m}", c.upper() if c == "dlbs" else "HCP-A", "all", span,
                int(s.n.max()))

    print("tract departure from assumed scanner axis, all sessions")
    for c, f in (("HCP-A", "roi_placement_quality_hcpa_b1500.csv"),
                 ("DLBS", "roi_placement_quality_dlbs_all.csv")):
        d = pd.read_csv(HERE / f)
        for a in ("theta_scr", "theta_slf", "theta_pvs", "theta_interfiber"):
            put(f"median_{a}", c, "tract", d[a].median(), len(d))
    v = pd.read_csv(HERE / "tract_direction_variation.csv")
    for _, r in v.iterrows():
        put(f"icc_{r.angle}", r.cohort, "tract", r.icc, int(r.n_repeat))

    pd.DataFrame(rows).to_csv(HERE / f"r2_text_numbers{SUFFIX}.csv", index=False)
    print(f"\nwrote {HERE / f'r2_text_numbers{SUFFIX}.csv'}")


if __name__ == "__main__":
    main()
