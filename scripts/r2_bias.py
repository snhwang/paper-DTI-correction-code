"""
R2 headline: how much does head position inflate the age decline of the classic index?

The age effect an index reports is the slope of log(index) on age, in percent per
decade, first session of each participant. The classic index is compared with
indices that carry no head-position error: template reorientation (DLBS only,
since HCP-A is already aligned), the refined index, and lambda2/lambda3. The
measured axis is reported but not used as a reference, because its own
attenuation drifts with age (Table 3).

    overstatement = slope(classic) / slope(reference) - 1

with 2000 paired participant bootstrap resamples, so the classic index and the
reference are always resampled together. Template reorientation reaches the same
conclusion as the refined index by an independent route (registration rather
than measured tract directions), which is what shows that the bias does not
depend on the correction used to remove it.

    python r2_bias.py

Writes r2_bias.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import atomic_io  # noqa: F401  writes become atomic on import

HERE = Path(__file__).resolve().parent
BOOT = 2000


def key(d):
    d = d.copy()
    d["Subject_ID"] = d.Subject_ID.astype(str)
    d["Visit"] = d.Visit.astype(str)
    return d


def first(d):
    return d.sort_values(["Subject_ID", "Visit"]).groupby("Subject_ID").first().reset_index()


def load(cohort):
    if cohort == "DLBS":
        m = key(pd.read_csv(HERE / "measured_pvs_axis_dlbs_all.csv"))
        v = key(pd.read_csv(HERE / "vecreg_comparison.csv").rename(columns={"Session": "Visit"}))
        m = m.merge(v[["Subject_ID", "Visit", "reoriented"]], on=["Subject_ID", "Visit"], how="left")
        return first(m), ("reoriented", "cross", "pv_perp", "v2_slab")
    m = key(pd.read_csv(HERE / "measured_pvs_axis_hcpa_b1500_all.csv"))
    return first(m), ("cross", "pv_perp", "v2_slab")


def slope(age, x):
    return 1000 * np.polyfit(age, np.log(x), 1)[0]


def main() -> None:
    rows = []
    rng = np.random.default_rng(20261007)
    for cohort in ("DLBS", "HCP-A"):
        d, refs = load(cohort)
        d = d.dropna(subset=["Age", "classic", *refs]).reset_index(drop=True)
        n = len(d)
        a = d.Age.to_numpy(float)
        X = {c: d[c].to_numpy(float) for c in ("classic", *refs)}
        boots = [rng.integers(0, n, n) for _ in range(BOOT)]
        print(f"\n{cohort}, {n} participants, age effect in % per decade")
        for c in ("classic", *refs):
            s = slope(a, X[c])
            b = [slope(a[i], X[c][i]) for i in boots]
            lo, hi = np.percentile(b, [2.5, 97.5])
            rows.append({"cohort": cohort, "index": c, "quantity": "age_effect_pct_per_decade",
                         "value": s, "lo": lo, "hi": hi, "n": n})
            print(f"   {c:10s} {s:+.2f} [{lo:+.2f}, {hi:+.2f}]")
        for c in refs:
            over = 100 * (slope(a, X["classic"]) / slope(a, X[c]) - 1)
            b = [100 * (slope(a[i], X["classic"][i]) / slope(a[i], X[c][i]) - 1) for i in boots]
            lo, hi = np.percentile(b, [2.5, 97.5])
            rows.append({"cohort": cohort, "index": c, "quantity": "classic_overstates_pct",
                         "value": over, "lo": lo, "hi": hi, "n": n})
            print(f"   classic overstates vs {c:10s} {over:5.0f}% [{lo:.0f}, {hi:.0f}]")
        if "reoriented" in refs:
            dif = slope(a, X["reoriented"]) - slope(a, X["cross"])
            b = [slope(a[i], X["reoriented"][i]) - slope(a[i], X["cross"][i]) for i in boots]
            lo, hi = np.percentile(b, [2.5, 97.5])
            rows.append({"cohort": cohort, "index": "reoriented-cross", "quantity": "difference_pct_per_decade",
                         "value": dif, "lo": lo, "hi": hi, "n": n})
            print(f"   reorientation minus refined {dif:+.2f} [{lo:+.2f}, {hi:+.2f}]")
    pd.DataFrame(rows).to_csv(HERE / "r2_bias.csv", index=False)
    print(f"\nwrote {HERE / 'r2_bias.csv'}")


if __name__ == "__main__":
    main()
