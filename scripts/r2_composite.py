"""
NOT USED IN THE MANUSCRIPT (decided 6 October). Kept as a record.

Why it was dropped. Any regression of age on these indices cannot cleanly test
whether the classic index is a composite of tissue and pose. The model
classic + pose counts pose twice, because the classic index already contains
it. The model corrected + pose against classic alone compares explicit pose
terms with pose that enters only through the attenuation, so it is close to
true by construction. The near-identical R2 of classic + pose and refined + pose
shows a tie, not a mechanism. The composite question is answered exactly, with
no regression, by the log decomposition in r2_decomposition.py (Table 3).

R2: how much of the classic index's stronger age association is the head pose it contains?

The classic index equals lambda2/lambda3 times an attenuation that follows head
pose, and head pose rises with age in DLBS. Its raw advantage over a corrected
index may therefore be pose. The test compares the variance in age that each
index explains, without and with head pose in both models:

    raw advantage        R2(age ~ classic) - R2(age ~ corrected)
    given pose           R2(age ~ classic + pose) - R2(age ~ corrected + pose)
    reduction            raw advantage - advantage given pose

Pose is absolute pitch and total rotation, as everywhere else. All three use
the same 2000 participant bootstrap resamples, so the reduction has a paired
interval. Run on the 156-participant DLBS sample (refined and lambda2/lambda3),
on the first sessions of all DLBS participants (refined only, from the
placement-quality table), and on HCP-A as the control with no pose to remove.

A stress test (6 October) showed that the earlier statement "head pose alone
predicts age better than the classic index" did not hold. It compared two pose
terms with one index, and with one pose term, out of sample, on ranks or in the
larger sample the two were equal or the index was stronger. It is not reported.
The reduction is reported because it held under every pose specification, on
ranks, out of sample and in both DLBS samples.

Registrations with a total rotation above 45 degrees are treated as failures
(one DLBS session, sub-4794 wave 1, 61 degrees of mostly yaw).

    python r2_composite.py

Writes r2_composite.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import atomic_io  # noqa: F401  writes become atomic on import

HERE = Path(__file__).resolve().parent
BOOT = 2000
MAX_ROTATION = 45.0


def key(d):
    d = d.copy()
    d["Subject_ID"] = d.Subject_ID.astype(str)
    d["Visit"] = d.Visit.astype(str)
    return d


def first(d):
    return d.sort_values(["Subject_ID", "Visit"]).groupby("Subject_ID").first().reset_index()


def r2(y, cols):
    X = np.column_stack([np.ones(len(y))] + cols)
    res = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    return 1 - res @ res / ((y - y.mean()) @ (y - y.mean()))


def williams(r12, r13, r23, n):
    det = 1 - r12**2 - r13**2 - r23**2 + 2 * r12 * r13 * r23
    rbar = (r12 + r13) / 2
    t = (r12 - r13) * np.sqrt((n - 1) * (1 + r23)
                              / (2 * (n - 1) / (n - 3) * det + rbar**2 * (1 - r23) ** 3))
    return float(2 * stats.t.sf(abs(t), n - 3))


def samples():
    h = key(pd.read_csv(HERE / "head_rotation_dlbs.csv"))
    h = h[h.total <= MAX_ROTATION]
    m = key(pd.read_csv(HERE / "measured_pvs_axis_dlbs.csv"))
    yield "DLBS", "main", first(m.merge(h, on=["Subject_ID", "Visit"])), ("cross", "pv_perp")
    q = key(pd.read_csv(HERE / "roi_placement_quality_dlbs_all.csv")).rename(columns={"refined_slab": "cross"})
    yield "DLBS", "all participants", first(q.merge(h, on=["Subject_ID", "Visit"])), ("cross",)
    hh = key(pd.read_csv(HERE / "head_rotation_hcpa.csv"))
    mh = key(pd.read_csv(HERE / "measured_pvs_axis_hcpa_b1500_all.csv"))
    yield "HCP-A", "main", first(mh.merge(hh, on=["Subject_ID", "Visit"])), ("cross", "v2_slab", "pv_perp")


def main() -> None:
    rows = []
    rng = np.random.default_rng(20261006)
    for cohort, sample, d, corrected in samples():
        d = d.dropna(subset=["Age", "classic", *corrected, "pitch", "total"]).reset_index(drop=True)
        n = len(d)
        y = d.Age.to_numpy(float)
        V = {c: d[c].to_numpy(float) for c in ("classic", *corrected)}
        P = [np.abs(d.pitch.to_numpy(float)), d.total.to_numpy(float)]
        boots = [rng.integers(0, n, n) for _ in range(BOOT)]

        def put(q, idx, v, lo=np.nan, hi=np.nan, p=np.nan):
            rows.append({"cohort": cohort, "sample": sample, "index": idx, "quantity": q,
                         "value": v, "lo": lo, "hi": hi, "p": p, "n": n})

        put("r2_pose", "pose", r2(y, P))
        put("r2", "classic", r2(y, [V["classic"]]))
        print(f"\n{cohort}, {sample}, {n} participants: R2 pose {r2(y, P):.3f}, classic {r2(y, [V['classic']]):.3f}")
        for c in corrected:
            def raw(i):
                return r2(y[i], [V["classic"][i]]) - r2(y[i], [V[c][i]])

            def cond(i):
                Pi = [p[i] for p in P]
                return r2(y[i], Pi + [V["classic"][i]]) - r2(y[i], Pi + [V[c][i]])

            a = np.arange(n)
            br = np.array([raw(i) for i in boots])
            bc = np.array([cond(i) for i in boots])
            red = br - bc
            pw = williams(np.corrcoef(y, V["classic"])[0, 1], np.corrcoef(y, V[c])[0, 1],
                          np.corrcoef(V["classic"], V[c])[0, 1], n)
            put("r2", c, r2(y, [V[c]]))
            put("raw_advantage", c, raw(a), *np.percentile(br, [2.5, 97.5]), p=pw)
            put("advantage_given_pose", c, cond(a), *np.percentile(bc, [2.5, 97.5]))
            put("reduction", c, raw(a) - cond(a), *np.percentile(red, [2.5, 97.5]))
            share = 100 * (1 - cond(a) / raw(a)) if raw(a) > 0 else np.nan
            put("share_removed_pct", c, share)
            print(f"  vs {c:8s} R2 {r2(y, [V[c]]):.3f}  raw advantage {raw(a):+.4f} "
                  f"[{np.percentile(br, 2.5):+.4f}, {np.percentile(br, 97.5):+.4f}] (Williams p {pw:.3f})  "
                  f"given pose {cond(a):+.4f} [{np.percentile(bc, 2.5):+.4f}, {np.percentile(bc, 97.5):+.4f}]  "
                  f"reduction {raw(a) - cond(a):+.4f} [{np.percentile(red, 2.5):+.4f}, {np.percentile(red, 97.5):+.4f}]  "
                  f"share {share:.0f}%")
    pd.DataFrame(rows).to_csv(HERE / "r2_composite.csv", index=False)
    print(f"\nwrote {HERE / 'r2_composite.csv'}")


if __name__ == "__main__":
    main()
