"""
R2: every number in the five R2 tables, recomputed from the current data.

Checking the R1 numbers for reuse found two kinds of drift. Table 6 (r of each
variant with lambda2/lambda3) had been pinned as text by the verifier rather
than recomputed, and its DLBS values no longer match the data. Supplement A.4
quoted reorientation reliability from the redrawn 5 mm spheres while every
table used the warped atlas regions. R2 therefore takes no number from the R1
text. Each table is regenerated here, on one region definition (the warped
atlas regions, measured_pvs_axis_*.csv), with one session per participant for
cross-sectional values and the longitudinal subset for reliability.

Indices carried in R2:
    classic   fixed scanner axes
    cross     refined, perpendicular to both band-estimated tract directions
    v2_slab   measured axis, pooled v2 over the tract band
    pv_perp   lambda2/lambda3 as the ratio of regional sums (Theory)
plus the published comparators ALPS-PAS, LD-ALPS and template reorientation.

    python r2_tables.py

Writes r2_numbers.csv (quantity, cohort, index, value, n), which the R2
manuscript and verifier read.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import atomic_io  # noqa: F401  writes become atomic on import

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent

# The DLBS main sample is every participant (284), matching HCP-A. Cross-sectional
# values use each participant's first session. Reliability and repositioning use the
# 156 with repeat visits, selected inside each analysis. R2_DLBS=repeat restores the
# R1 restriction to those 156 for every analysis, and suffixes the outputs.
import os as _os
DLBS_MAIN = ("measured_pvs_axis_dlbs.csv" if _os.environ.get("R2_DLBS") == "repeat"
             else "measured_pvs_axis_dlbs_all.csv")
SUFFIX = "_dlbs_repeat" if _os.environ.get("R2_DLBS") == "repeat" else ""
MAX_ROTATION = 45.0   # above this a registration is treated as failed (Methods 3.3)
sys.path.insert(0, str(HERE))
from estimator_variants import variance_components  # noqa: E402

N_PERM = 2000
OURS = ["classic", "cross", "v2_slab", "pv_perp"]
FILES = {"HCP-A": ("measured_pvs_axis_hcpa_b1500_all.csv", "head_rotation_hcpa.csv",
                   "comparators_hcpa.csv", "ld_alps_hcpa.csv"),
         "DLBS": (DLBS_MAIN, "head_rotation_dlbs.csv",
                  "comparators_dlbs.csv", "ld_alps_dlbs.csv")}

rows: list[dict] = []


def put(quantity, cohort, index, value, n=None):
    rows.append({"quantity": quantity, "cohort": cohort, "index": index,
                 "value": None if value is None else float(value), "n": n})


def key(d):
    d = d.copy()
    d["Subject_ID"] = d.Subject_ID.astype(str)
    d["Visit"] = d.Visit.astype(str)
    return d


def one_per(d):
    return d.sort_values(["Subject_ID", "Visit"]).groupby("Subject_ID").first().reset_index()


def z(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / v.std(ddof=1)


def beta_age(d, col, extra=()):
    X = np.column_stack([np.ones(len(d)), z(d.Age)] + [z(e) for e in extra])
    return float(np.linalg.lstsq(X, z(d[col]), rcond=None)[0][1])


def partial_r(d, col, ctrl):
    """Correlation of age with col after regressing ctrl out of both."""
    X = np.column_stack([np.ones(len(d)), d[ctrl]])
    ra = d.Age - X @ np.linalg.lstsq(X, d.Age, rcond=None)[0]
    rc = d[col] - X @ np.linalg.lstsq(X, d[col], rcond=None)[0]
    r = float(np.corrcoef(ra, rc)[0, 1])
    df = len(d) - 3
    p = float(2 * stats.t.sf(abs(r) * np.sqrt(df / (1 - r * r)), df))
    return r, p


def load(cohort):
    mf, hf, cf, lf = FILES[cohort]
    d = key(pd.read_csv(HERE / mf))
    c = key(pd.read_csv(HERE / cf))[["Subject_ID", "Visit", "ALPS-PAS"]]
    d = d.merge(c, on=["Subject_ID", "Visit"], how="left")
    if (HERE / lf).exists():
        L = key(pd.read_csv(HERE / lf))[["Subject_ID", "Visit", "ALPS_overall"]]
        d = d.merge(L.rename(columns={"ALPS_overall": "LD-ALPS"}),
                    on=["Subject_ID", "Visit"], how="left")
    if cohort == "DLBS":
        v = pd.read_csv(HERE / "vecreg_comparison.csv").rename(columns={"Session": "Visit"})
        v = key(v)[["Subject_ID", "Visit", "reoriented"]]
        d = d.merge(v, on=["Subject_ID", "Visit"], how="left")
    h = key(pd.read_csv(HERE / hf))
    h = h[h.total <= MAX_ROTATION]
    return d, h


def main() -> None:
    for cohort in ("HCP-A", "DLBS"):
        d, h = load(cohort)
        idx = OURS + [c for c in ("ALPS-PAS", "LD-ALPS", "reoriented") if c in d]
        one = one_per(d)
        print(f"\n===== {cohort}: {len(d)} sessions, {len(one)} participants")

        # Table 1: attainment of the bound and correlation with it
        print("\nTable 1  median fraction of lambda2/lambda3, and r with it (one per participant)")
        for c in idx:
            s = one.dropna(subset=[c, "pv_perp"])
            att = float((s[c] / s.pv_perp).median())
            r = float(np.corrcoef(s[c], s.pv_perp)[0, 1])
            above = float(100 * (d[c] > d.pv_perp).mean()) if c != "pv_perp" else np.nan
            put("attainment", cohort, c, att, len(s))
            put("r_with_ratio", cohort, c, r, len(s))
            put("pct_sessions_above_bound", cohort, c, above, int(d[c].notna().sum()))
            print(f"   {c:11s} attains {att:.3f}   r {r:.3f}   above bound {above:5.2f}% of sessions")

        # Table 2: head pose against age
        m = one_per(h.merge(d[["Subject_ID", "Visit", "Age"] + OURS], on=["Subject_ID", "Visit"]))
        m = m.dropna(subset=["Age", "pitch", "total", "classic"])
        thr = 10 if cohort == "DLBS" else 5
        put("median_total_rotation", cohort, "pose", m.total.median(), len(m))
        put(f"pct_beyond_{thr}deg", cohort, "pose", 100 * (m.total > thr).mean(), len(m))
        print(f"\nTable 2  pose, {len(m)} participants: median total {m.total.median():.1f} deg, "
              f"{100 * (m.total > thr).mean():.1f}% beyond {thr} deg")
        for a in ("pitch", "total", "yaw", "roll"):
            r, p = stats.pearsonr(m.Age, m[a].abs())
            put(f"age_r_{a}", cohort, "pose", r, len(m))
            put(f"age_p_{a}", cohort, "pose", p, len(m))
            print(f"   |{a}| vs age r {r:+.3f}  p {p:.1e}")

        # Table 3: age coefficient absorbed by pose adjustment
        print("\nTable 3  standardized age coefficient, raw -> pose-adjusted (|pitch|, total)")
        for c in OURS:
            b0 = beta_age(m, c)
            b1 = beta_age(m, c, (m.pitch.abs(), m.total))
            pct = 100 * (1 - abs(b1) / abs(b0))
            put("age_beta_raw", cohort, c, b0, len(m))
            put("age_beta_pose", cohort, c, b1, len(m))
            put("pct_absorbed_pose", cohort, c, pct, len(m))
            print(f"   {c:9s} {b0:+.3f} -> {b1:+.3f}   {pct:5.1f}% absorbed")
        rng = np.random.default_rng(20260811)
        obs = 100 * (1 - abs(beta_age(m, "classic", (m.pitch.abs(), m.total)))
                     / abs(beta_age(m, "classic")))
        null = []
        for _ in range(N_PERM):
            pp = m.pitch.abs().to_numpy()[rng.permutation(len(m))]
            tt = m.total.to_numpy()[rng.permutation(len(m))]
            null.append(100 * (1 - abs(beta_age(m, "classic", (pp, tt)))
                               / abs(beta_age(m, "classic"))))
        null = np.array(null)
        pperm = float((null >= obs).mean())
        for q, v in (("perm_mean", null.mean()), ("perm_p95", np.percentile(null, 95)),
                     ("perm_max", null.max()), ("perm_p", pperm)):
            put(q, cohort, "classic", v, len(m))
        print(f"   permutation: mean {null.mean():+.2f}%  95th {np.percentile(null, 95):.1f}%  "
              f"max {null.max():.1f}%  p {pperm:.4f}")

        # Table 4: age association after partialling lambda2/lambda3
        print("\nTable 4  age r raw | partial on lambda2/lambda3 (one per participant)")
        for c in idx:
            s = one.dropna(subset=[c, "pv_perp", "Age"])
            raw = float(np.corrcoef(s.Age, s[c])[0, 1])
            put("age_r_raw", cohort, c, raw, len(s))
            if c == "pv_perp":
                print(f"   {c:11s} {raw:+.3f} |   --")
                continue
            pr, pp = partial_r(s, c, "pv_perp")
            put("age_r_partial_ratio", cohort, c, pr, len(s))
            put("age_p_partial_ratio", cohort, c, pp, len(s))
            print(f"   {c:11s} {raw:+.3f} | {pr:+.3f}  p {pp:.3f}")

        # Table 5: reliability and age over all sessions
        lon = d[d.Subject_ID.isin(d.Subject_ID.value_counts()[lambda s: s >= 2].index)]
        print(f"\nTable 5  ICC(1,1) over {lon.Subject_ID.nunique()} repeat participants, "
              f"age r over all sessions")
        for c in idx:
            s = lon.dropna(subset=[c])
            icc = variance_components(s, c)["icc"] if len(s) > 20 else np.nan
            a = d.dropna(subset=[c, "Age"])
            r = float(np.corrcoef(a.Age, a[c])[0, 1])
            put("icc", cohort, c, icc, int(s.Subject_ID.nunique()))
            put("age_r_all_sessions", cohort, c, r, len(a))
            print(f"   {c:11s} ICC {icc:.3f}   age r {r:+.3f}  (n {len(a)})")

    # Table 3, second sample: DLBS first sessions of all 284 participants
    q = key(pd.read_csv(HERE / "roi_placement_quality_dlbs_all.csv"))
    h = key(pd.read_csv(HERE / "head_rotation_dlbs.csv"))
    # a total rotation above 45 degrees is a failed registration, not a posture
    # (one session, sub-4794 wave 1, 61 degrees of mostly yaw); see r2_composite.py
    h = h[h.total <= 45.0]
    m = one_per(q.merge(h[["Subject_ID", "Visit", "pitch", "total"]], on=["Subject_ID", "Visit"]))
    m = m.dropna(subset=["Age", "classic", "refined_slab", "pitch", "total"])
    print(f"\nTable 3, full DLBS sample ({len(m)} participants)")
    for c, lab in (("classic", "classic"), ("refined_slab", "cross")):
        b0, b1 = beta_age(m, c), beta_age(m, c, (m.pitch.abs(), m.total))
        pct = 100 * (1 - abs(b1) / abs(b0))
        put("pct_absorbed_pose_full", "DLBS", lab, pct, len(m))
        print(f"   {lab:9s} {b0:+.3f} -> {b1:+.3f}   {pct:5.1f}% absorbed")

    pd.DataFrame(rows).to_csv(HERE / f"r2_numbers{SUFFIX}.csv", index=False)
    print(f"\nwrote {HERE / f'r2_numbers{SUFFIX}.csv'}")


if __name__ == "__main__":
    main()
