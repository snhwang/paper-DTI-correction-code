"""
R2: template reorientation on the same sample and regions as the comparison table.

The R1 text (Supplement A.4) quoted ICC 0.607 classic, 0.604 reorientation and
0.486 refined. Those came from vecreg_comparison_sphere5.csv, the redrawn 5 mm
spheres, while the comparison table (Table 9, build_variant_table.py) used the
warped atlas regions. R2 uses the warped regions throughout, so the
reorientation row is taken from vecreg_comparison.csv, which evaluates it in
the identical regions and voxels as classic.

Reorientation here is the native-space form described in vecreg_comparison.py:
the fixed axes are rotated by the transpose of the rotation of each session's
subject-to-template affine, which is equivalent to reorienting the tensors and
involves no resampling. The fully resampled FSL vecreg version is in
vecreg_summary.py.

Endpoints match build_variant_table.py: ICC(1,1) over participants with repeat
visits, and a Pearson age correlation over every session.

    python r2_reorientation.py

Writes r2_reorientation.csv.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import atomic_io  # noqa: F401  writes become atomic on import

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from estimator_variants import variance_components  # noqa: E402


def main() -> None:
    m = pd.read_csv(HERE / "measured_pvs_axis_dlbs.csv")
    v = pd.read_csv(HERE / "vecreg_comparison.csv").rename(columns={"Session": "Visit"})
    for d in (m, v):
        d["Subject_ID"] = d.Subject_ID.astype(str)
        d["Visit"] = d.Visit.astype(str)
    d = m[["Subject_ID", "Visit", "Age", "classic", "cross", "pv_perp"]].merge(
        v[["Subject_ID", "Visit", "classic", "reoriented"]],
        on=["Subject_ID", "Visit"], suffixes=("", "_v"))
    assert np.allclose(d.classic, d.classic_v), "classic differs between files"
    lon = d[d.Subject_ID.isin(d.Subject_ID.value_counts()[lambda s: s >= 2].index)]
    print(f"{len(d)} sessions, {d.Subject_ID.nunique()} participants, "
          f"{lon.Subject_ID.nunique()} with repeat visits\n")
    base = variance_components(lon, "classic")["var_within"] / d.classic.mean() ** 2
    rows = []
    for name, c in (("Classic", "classic"), ("Template reorientation", "reoriented"),
                    ("Refined", "cross")):
        vc = variance_components(lon, c)
        rel = vc["var_within"] / d[c].mean() ** 2
        r_age = float(np.corrcoef(d.Age, d[c])[0, 1])
        r_ratio = float(np.corrcoef(d.pv_perp, d[c])[0, 1])
        rows.append({"index": name, "icc": round(vc["icc"], 3),
                     "rel_within_var_vs_classic": round(rel / base, 3),
                     "age_r": round(r_age, 3), "r_with_ratio": round(r_ratio, 3),
                     "n_sessions": len(d), "n_participants": d.Subject_ID.nunique()})
        print(f"  {name:24s} ICC {vc['icc']:.3f}  within-var x{rel / base:.2f}  "
              f"age r {r_age:+.3f}  r with ratio {r_ratio:.3f}")
    pd.DataFrame(rows).to_csv(HERE / "r2_reorientation.csv", index=False)
    print(f"\nwrote {HERE / 'r2_reorientation.csv'}")


if __name__ == "__main__":
    main()
