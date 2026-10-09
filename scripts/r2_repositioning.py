"""
R2: sensitivity to genuine repositioning, on the same pipeline as every table.

The R1 figure (Fig S1) took its slopes from the older sphere pipeline
(longitudinal_reliability.py --cohort spheres), whose refined index estimates
tract directions from the measurement spheres. Every R2 table uses the
band-estimated directions instead, so the figure is recomputed here on that
basis, with the measured axis and the eigenvalue ratio added.

Design as before. For every pair of visits of the same participant, the
relative change in each index, 100 |b - a| / mean(a, b), is regressed on the
absolute change in theta_SCR, the angle between the projection-tract direction
and scanner z. Within a participant the tract anatomy is fixed, so that change
is head repositioning. Session level, hemispheres averaged. Confidence
intervals from a participant-clustered bootstrap.

    python r2_repositioning.py

Writes r2_repositioning.csv (slopes) and r2_repositioning_pairs.csv (the pairs,
for the figure).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import atomic_io  # noqa: F401  writes become atomic on import

HERE = Path(__file__).resolve().parent
BOOT = 2000
INDICES = {"Classic": "classic", "Refined": "cross", "Measured axis": "v2_slab",
           "lambda2/lambda3": "pv_perp"}


def load() -> pd.DataFrame:
    q = pd.read_csv(HERE / "roi_placement_quality_dlbs_all.csv")
    m = pd.read_csv(HERE / "measured_pvs_axis_dlbs_all.csv")
    for d in (q, m):
        d["Subject_ID"] = d.Subject_ID.astype(str)
        d["Visit"] = d.Visit.astype(str)
    d = q[["Subject_ID", "Visit", "theta_scr", "classic", "refined_slab"]].merge(
        m[["Subject_ID", "Visit", "classic", "cross", "v2_slab", "pv_perp"]],
        on=["Subject_ID", "Visit"], suffixes=("_q", ""))
    # the two files are the same pipeline; confirm rather than assume
    assert np.allclose(d.classic_q, d.classic) and np.allclose(d.refined_slab, d.cross)
    d["wave"] = d.Visit.str.extract(r"(\d+)").astype(int)
    return d


def pairs(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for s, g in d.sort_values("wave").groupby("Subject_ID"):
        g = g.reset_index(drop=True)
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                a, b = g.iloc[i], g.iloc[j]
                r = {"Subject_ID": s, "wave_a": a.wave, "wave_b": b.wave,
                     "d_theta_SCR": abs(b.theta_scr - a.theta_scr)}
                for name, c in INDICES.items():
                    r[name] = 100 * abs(b[c] - a[c]) / ((a[c] + b[c]) / 2)
                rows.append(r)
    return pd.DataFrame(rows)


def main() -> None:
    p = pairs(load())
    subs = p.Subject_ID.unique()
    groups = {s: p[p.Subject_ID == s] for s in subs}
    rng = np.random.default_rng(0)
    print(f"{len(p)} visit pairs from {len(subs)} participants, "
          f"median |d theta_SCR| {p.d_theta_SCR.median():.2f} deg\n")
    out = []
    for name in INDICES:
        slope = float(np.polyfit(p.d_theta_SCR, p[name], 1)[0])
        boot = []
        for _ in range(BOOT):
            r = pd.concat([groups[x] for x in rng.choice(subs, len(subs), replace=True)])
            boot.append(float(np.polyfit(r.d_theta_SCR, r[name], 1)[0]))
        lo, hi = np.percentile(boot, [2.5, 97.5])
        out.append({"index": name, "slope_pct_per_deg": round(slope, 4),
                    "ci_lo": round(float(lo), 4), "ci_hi": round(float(hi), 4),
                    "n_pairs": len(p), "n_participants": len(subs)})
        print(f"  {name:16s} {slope:+.3f}% per degree  CI [{lo:+.2f}, {hi:+.2f}]")
    pd.DataFrame(out).to_csv(HERE / "r2_repositioning.csv", index=False)
    p.to_csv(HERE / "r2_repositioning_pairs.csv", index=False)
    print(f"\nwrote {HERE / 'r2_repositioning.csv'}")


if __name__ == "__main__":
    main()
