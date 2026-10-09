"""
R2: where each index's between-visit reliability comes from.

log index = log(lambda2/lambda3) + log(attenuation), so for both the between- and
the within-participant variance component

    V(index) = V(ratio) + V(attenuation) + 2 C(ratio, attenuation).

Components from the one-way random-effects ANOVA over participants with repeat
visits (156 in DLBS, 628 in HCP-A), on the log scale. The covariance component
is obtained exactly as [V(ratio + attenuation) - V(ratio) - V(attenuation)] / 2.

The question is whether the classic index's higher reliability in DLBS reflects
a better measurement of tissue or an error that repeats. Two checks:
  - within-participant noise of the attenuation, classic against refined
  - whether each participant's mean attenuation follows their mean |pitch|
    (habitual posture), with the between-visit ICC of |pitch| itself

    python r2_reliability_parts.py

Writes r2_reliability_parts.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import atomic_io  # noqa: F401  writes become atomic on import

HERE = Path(__file__).resolve().parent
MAX_ROTATION = 45.0
COHORTS = {"DLBS": ("measured_pvs_axis_dlbs_all.csv", "head_rotation_dlbs.csv"),
           "HCP-A": ("measured_pvs_axis_hcpa_b1500_all.csv", "head_rotation_hcpa.csv")}
INDICES = {"classic": "classic", "cross": "refined", "v2_slab": "measured axis"}


def vc(d, col):
    """Between and within variance components, unbalanced one-way ANOVA (unclipped)."""
    d = d[["Subject_ID", col]].dropna()
    k = d.Subject_ID.nunique()
    n_i = d.groupby("Subject_ID")[col].size().to_numpy(float)
    m_i = d.groupby("Subject_ID")[col].mean().to_numpy(float)
    N = n_i.sum()
    yb = (n_i * m_i).sum() / N
    msb = (n_i * (m_i - yb) ** 2).sum() / (k - 1)
    msw = d.groupby("Subject_ID")[col].transform(lambda s: s - s.mean()).pow(2).sum() / (N - k)
    n0 = (N - (n_i ** 2).sum() / N) / (k - 1)
    return (msb - msw) / n0, msw


def main() -> None:
    rows = []
    for cohort, (mf, hf) in COHORTS.items():
        d = pd.read_csv(HERE / mf)
        d["Subject_ID"] = d.Subject_ID.astype(str)
        d["Visit"] = d.Visit.astype(str)
        d = d[d.Subject_ID.isin(d.Subject_ID.value_counts()[lambda s: s >= 2].index)].copy()
        k = d.Subject_ID.nunique()
        d["r"] = np.log(d.pv_perp)
        print(f"\n{cohort}: {k} participants, {len(d)} sessions (variances x1e3, log scale)")
        br, wr = vc(d, "r")
        rows.append({"cohort": cohort, "index": "ratio", "between": br, "within": wr,
                     "icc": br / (br + wr), "n_participants": k})
        print(f"  ratio       between {1e3*br:6.2f}  within {1e3*wr:6.2f}  ICC {br/(br+wr):.3f}")
        for c, name in INDICES.items():
            d["a"] = np.log(d[c] / d.pv_perp)
            d["s"] = d.r + d.a
            ba, wa = vc(d, "a")
            bs, ws = vc(d, "s")
            cb, cw = (bs - br - ba) / 2, (ws - wr - wa) / 2
            rows.append({"cohort": cohort, "index": c, "between": bs, "within": ws,
                         "icc": bs / (bs + ws), "atten_between": ba, "atten_within": wa,
                         "atten_icc": ba / (ba + wa), "cov_between": cb, "cov_within": cw,
                         "n_participants": k})
            print(f"  {name:13s} between {1e3*bs:6.2f} = {1e3*br:5.2f} + atten {1e3*ba:5.2f} + 2cov {2e3*cb:+6.2f}   "
                  f"within {1e3*ws:6.2f} = {1e3*wr:5.2f} + atten {1e3*wa:5.2f} + 2cov {2e3*cw:+6.2f}   "
                  f"ICC {bs/(bs+ws):.3f}")
        h = pd.read_csv(HERE / hf)
        h["Subject_ID"] = h.Subject_ID.astype(str)
        h["Visit"] = h.Visit.astype(str)
        m = d.merge(h[h.total <= MAX_ROTATION], on=["Subject_ID", "Visit"])
        m["ap"] = m.pitch.abs()
        bp, wp = vc(m, "ap")
        rows.append({"cohort": cohort, "index": "abs_pitch", "between": bp, "within": wp,
                     "icc": bp / (bp + wp), "n_participants": m.Subject_ID.nunique()})
        print(f"  |pitch| ICC between visits {bp/(bp+wp):.3f}")
        for c, name in (("classic", "classic"), ("cross", "refined")):
            m["a"] = np.log(m[c] / m.pv_perp)
            g = m.groupby("Subject_ID")[["ap", "a"]].mean()
            r = float(np.corrcoef(g.ap, g.a)[0, 1])
            rows.append({"cohort": cohort, "index": c, "habitual_pitch_r": r,
                         "n_participants": len(g)})
            print(f"  participant means: {name} attenuation against |pitch| r {r:+.3f} ({len(g)} participants)")
    pd.DataFrame(rows).to_csv(HERE / "r2_reliability_parts.csv", index=False)
    print(f"\nwrote {HERE / 'r2_reliability_parts.csv'}")


if __name__ == "__main__":
    main()
