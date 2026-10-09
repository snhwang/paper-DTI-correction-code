"""
R2 headline: how much of each index's age decline is axis error.

Theory gives index = rho_bar * a, where rho_bar is lambda2/lambda3 as the ratio
of regional sums (pv_perp) and a = index / rho_bar is the attenuation set by
the angle between the measurement axis and v2. On a log scale the age slopes
add exactly,

    slope(log index) = slope(log rho_bar) + slope(log a),

so the age decline splits into a tissue term, which no rotation can change,
and an axis-error term. Adjusting each term for head pose (|pitch| and total
rotation, as everywhere else) then shows how much of the axis-error term is
head position.

This replaces the standardized-coefficient attenuation (45%) as the headline,
because in DLBS head pose also tracks rho_bar itself (r = -0.32), so the 45%
mixes axis error with an association between pose and tissue that no choice
of axes can reach.

One session per participant, slopes in percent per decade, 95% intervals from
2000 participant bootstrap resamples.

    python r2_decomposition.py

Writes r2_decomposition.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import atomic_io  # noqa: F401  writes become atomic on import

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
BOOT = 2000
COHORTS = {"DLBS": (DLBS_MAIN, "head_rotation_dlbs.csv"),
           "HCP-A": ("measured_pvs_axis_hcpa_b1500_all.csv", "head_rotation_hcpa.csv")}
INDICES = {"classic": "Classic", "cross": "Refined", "v2_slab": "Measured axis"}


def key(d):
    d = d.copy()
    d["Subject_ID"] = d.Subject_ID.astype(str)
    d["Visit"] = d.Visit.astype(str)
    return d


def load(cohort):
    mf, hf = COHORTS[cohort]
    h = key(pd.read_csv(HERE / hf))
    h = h[h.total <= MAX_ROTATION]
    d = key(pd.read_csv(HERE / mf)).merge(h, on=["Subject_ID", "Visit"])
    d = d.sort_values(["Subject_ID", "Visit"]).groupby("Subject_ID").first().reset_index()
    return d.dropna(subset=["Age", "pv_perp", "pitch", "total", *INDICES]).reset_index(drop=True)


def slope(y, X):
    """Age coefficient, in percent per decade of a log quantity."""
    A = np.column_stack([np.ones(len(y))] + X)
    return 1000 * np.linalg.lstsq(A, y, rcond=None)[0][1]


def terms(d, c):
    age = d.Age.to_numpy(float)
    pose = [age, np.abs(d.pitch.to_numpy(float)), d.total.to_numpy(float)]
    li = np.log(d[c].to_numpy(float))
    lr = np.log(d.pv_perp.to_numpy(float))
    la = li - lr
    t = {"total": slope(li, [age]), "ratio": slope(lr, [age]), "attenuation": slope(la, [age]),
         "total_pose": slope(li, pose), "ratio_pose": slope(lr, pose),
         "attenuation_pose": slope(la, pose)}
    t["axis_share_pct"] = 100 * t["attenuation"] / t["total"]
    t["pose_share_of_axis_pct"] = 100 * (1 - t["attenuation_pose"] / t["attenuation"])
    t["pose_share_of_total_via_axis_pct"] = 100 * (t["attenuation"] - t["attenuation_pose"]) / t["total"]
    return t


def main() -> None:
    rows = []
    rng = np.random.default_rng(20261006)
    for cohort in COHORTS:
        d = load(cohort)
        n = len(d)
        print(f"\n{cohort}, {n} participants (percent per decade, 95% bootstrap interval)")
        r_rp = stats.pearsonr(d.pv_perp, d.pitch.abs())
        print(f"   lambda2/lambda3 vs |pitch|  r {r_rp[0]:+.3f}  p {r_rp[1]:.1e}")
        rows.append({"cohort": cohort, "index": "ratio", "quantity": "r_with_abs_pitch",
                     "value": r_rp[0], "lo": np.nan, "hi": np.nan, "n": n})
        boots = [d.iloc[rng.integers(0, n, n)].reset_index(drop=True) for _ in range(BOOT)]
        for c, name in INDICES.items():
            t = terms(d, c)
            bt = pd.DataFrame([terms(b, c) for b in boots])
            ra = stats.pearsonr(d[c] / d.pv_perp, d.pitch.abs())
            print(f"   {name}")
            for q, v in t.items():
                lo, hi = np.percentile(bt[q], [2.5, 97.5])
                rows.append({"cohort": cohort, "index": c, "quantity": q, "value": v,
                             "lo": lo, "hi": hi, "n": n})
                print(f"      {q:34s} {v:+7.2f}   [{lo:+.2f}, {hi:+.2f}]")
            rows.append({"cohort": cohort, "index": c, "quantity": "attenuation_r_with_abs_pitch",
                         "value": ra[0], "lo": np.nan, "hi": np.nan, "n": n})
            print(f"      attenuation vs |pitch| r {ra[0]:+.3f}  p {ra[1]:.1e}")
    pd.DataFrame(rows).to_csv(HERE / f"r2_decomposition{SUFFIX}.csv", index=False)
    print(f"\nwrote {HERE / f'r2_decomposition{SUFFIX}.csv'}")


if __name__ == "__main__":
    main()
