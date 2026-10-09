"""
R2 verifier: every number the R2 manuscript states, checked against the R2 outputs.

verify_manuscript.py pins R1 sentences as text and no longer applies to the
restructured paper. This replaces it for R2. Each check takes a value from an
r2_*.csv file, formats it the way the manuscript prints it, and requires that
string to appear within a window after an anchor phrase, so that a number is
checked in the sentence that states it and not anywhere in the file. The four
generated tables are compared verbatim with r2_tables.tex.

Run the r2_*.py scripts first (r2_tables, r2_decomposition, r2_text_numbers,
r2_repositioning, r2_bias, r2_reliability_parts, r2_latex_tables).

    python verify_r2.py
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
# The public code repository ships this script and the r2_*.csv files but not
# the manuscript, which is not ours to publish before acceptance. Without it the
# text checks cannot run, so they are skipped loudly rather than crashing, and
# the table regeneration still runs. In the manuscript repository the file is
# always there and nothing is skipped.
_tex = HERE.parent / "mri_revision.tex"
if not _tex.exists():
    _tex = HERE / "mri_revision.tex"
HAVE_TEX = _tex.exists()
TEX = _tex.read_text(encoding="utf-8") if HAVE_TEX else ""
FLAT = re.sub(r"\s+", " ", TEX)
N = pd.read_csv(HERE / "r2_numbers.csv")
DEC = pd.read_csv(HERE / "r2_decomposition.csv")
TXT = pd.read_csv(HERE / "r2_text_numbers.csv")
BIAS = pd.read_csv(HERE / "r2_bias.csv")
REP = pd.read_csv(HERE / "r2_repositioning.csv").set_index("index")
REL = pd.read_csv(HERE / "r2_reliability_parts.csv")

fails: list[str] = []
count = 0


def check(desc, anchor, expected, window=1200):
    global count
    # Every check here reads the manuscript. Without it there is nothing to
    # read, so they are not counted as passing and not reported as failing.
    if not HAVE_TEX:
        return
    count += 1
    i = FLAT.find(anchor)
    if i < 0:
        fails.append(f"{desc}: anchor not found: {anchor!r}")
        return
    if expected not in FLAT[i: i + len(anchor) + window]:
        fails.append(f"{desc}: expected {expected!r} after {anchor!r}")


def v(df, **kw):
    s = df
    for k, val in kw.items():
        s = s[s[k] == val]
    if len(s) != 1:
        raise SystemExit(f"lookup {kw} matched {len(s)} rows")
    return s.iloc[0]


def num(df, quantity, cohort, index):
    return float(v(df, quantity=quantity, cohort=cohort, index=index).value)


def r(x, d=3, sign=True):
    s = f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    return s


def pexp(p):
    """p as the manuscript prints a small value, e.g. 2\\times10^{-10}."""
    e = math.floor(math.log10(p))
    m = round(p / 10 ** e)
    if m == 10:
        m, e = 1, e + 1
    return f"{m}\\times10^{{{e}}}"


def pct(x, d=0):
    return f"{x:.{d}f}\\%"


# ---- 4.1 head pose and age, DLBS ----
a = "In DLBS, head pitch increased with age ($r=+"
check("pitch r", a, f"r={r(num(N, 'age_r_pitch', 'DLBS', 'pose'))}")
check("pitch p", a, pexp(num(N, "age_p_pitch", "DLBS", "pose")))
check("pitch n", a, f"${int(v(N, quantity='age_r_pitch', cohort='DLBS', index='pose').n)}$ participants")
check("yaw r", a, f"r={r(num(N, 'age_r_yaw', 'DLBS', 'pose'))}")
check("yaw p", a, pexp(num(N, "age_p_yaw", "DLBS", "pose")))
check("total r", a, f"r={r(num(N, 'age_r_total', 'DLBS', 'pose'))}")
check("roll r", a, f"r={r(num(N, 'age_r_roll', 'DLBS', 'pose'))}")
check("roll p", a, f"p={num(N, 'age_p_roll', 'DLBS', 'pose'):.2f}")
check("median rotation", "The median total rotation from the template was",
      f"{num(N, 'median_total_rotation', 'DLBS', 'pose'):.1f}^{{\\circ}}")
check("beyond 10", "The median total rotation from the template was",
      pct(num(N, "pct_beyond_10deg", "DLBS", "pose")))
check("HCP-A median rotation", "In HCP-A, where preprocessing had removed head position, the median was",
      f"{num(N, 'median_total_rotation', 'HCP-A', 'pose'):.1f}^{{\\circ}}")
hcp_r = max(abs(num(N, f"age_r_{q}", "HCP-A", "pose")) for q in ("pitch", "yaw", "roll", "total"))
check("HCP-A no rotation with age", "no rotation was associated with age", f"|r|\\le{math.ceil(hcp_r * 100) / 100:.2f}")
abs_r = abs(num(N, "age_r_pitch", "DLBS", "pose"))
check("abstract pitch r", "In DLBS, head pitch increased with age ($r=0.", f"r={abs_r:.2f}")
check("abstract pitch p", "In DLBS, head pitch increased with age ($r=0.", pexp(num(N, "age_p_pitch", "DLBS", "pose")))
a = "In DLBS the median slab angulation was"
check("slab median", a, f"{num(TXT, 'median_slab_angulation_deg', 'DLBS', 'prescription'):.1f}^{{\\circ}}")
check("slab r", a, f"r={r(num(TXT, 'age_r_slab_pitch', 'DLBS', 'prescription'))}")
check("slab p", a, pexp(num(TXT, "age_p_slab_pitch", "DLBS", "prescription")))
a = "Splitting DLBS at the median body mass index"
check("BMI high", a, f"{num(TXT, 'median_abs_pitch_high_bmi', 'DLBS', 'pose'):.1f}^{{\\circ}}")
check("BMI low", a, f"{num(TXT, 'median_abs_pitch_low_bmi', 'DLBS', 'pose'):.1f}^{{\\circ}}")
check("BMI p", a, f"p={num(TXT, 'welch_p_pitch_by_bmi', 'DLBS', 'pose'):.3f}")
check("BMI n", a, f"${int(v(TXT, quantity='welch_p_pitch_by_bmi', cohort='DLBS', index='pose').n)}$ participants")
check("BMI age r", a, f"r={num(TXT, 'age_r_bmi', 'DLBS', 'pose'):.2f}")

# ---- 4.2 bias ----
b = lambda c, i, q: v(BIAS, cohort=c, index=i, quantity=q)
a = "In DLBS the classic index declined by"
cl = b("DLBS", "classic", "age_effect_pct_per_decade")
check("classic decline", a, f"{-cl.value:.1f}\\%")
check("classic decline CI", a, f"{-cl.hi:.1f}$ to ${-cl.lo:.1f}")
check("reorientation decline", a, f"declined by ${-b('DLBS', 'reoriented', 'age_effect_pct_per_decade').value:.1f}\\%")
check("refined decline", a, f"declined by ${-b('DLBS', 'cross', 'age_effect_pct_per_decade').value:.1f}\\%")
for i in ("reoriented", "cross"):
    o = b("DLBS", i, "classic_overstates_pct")
    check(f"overstates vs {i}", a, f"${o.value:.0f}\\%$ (${o.lo:.0f}$ to ${o.hi:.0f}\\%$)")
o = b("HCP-A", "cross", "classic_overstates_pct")
check("HCP-A overstates", a, f"${o.value:.0f}\\%$ (${o.lo:.0f}$ to ${o.hi:.0f}\\%$)", window=1200)
o = b("DLBS", "reoriented-cross", "difference_pct_per_decade")
check("reorientation minus refined", a, f"{-o.value:.1f}\\%$ per decade more, ${-o.hi:.2f}$ to ${-o.lo:.2f}")
for i, lab in (("reoriented", "against template reorientation"), ("cross", "against a simple correction")):
    o = b("DLBS", i, "classic_overstates_pct")
    check(f"abstract overstates {i}", "In DLBS the classic index overstated the age decline by",
          f"${o.value:.0f}\\%$", window=300)

# ---- 4.2 decomposition ----
d = lambda c, i, q: v(DEC, cohort=c, index=i, quantity=q)
a = "Table~\\ref{tbl:decomposition} shows where the excess comes from."
check("ratio decline", a, f"${-d('DLBS', 'classic', 'ratio').value:.1f}\\%$ per decade of the decline")
att = d("DLBS", "classic", "attenuation")
check("attenuation decline", a, f"${-att.value:.1f}\\%$ per decade ($95\\%$ interval ${-att.hi:.1f}$ to ${-att.lo:.1f}$)")
sh = d("DLBS", "classic", "axis_share_pct")
check("axis share", a, f"${sh.value:.0f}\\%$ (${sh.lo:.0f}$ to ${sh.hi:.0f}\\%$)")
check("attenuation vs pitch", a, f"r={d('DLBS', 'classic', 'attenuation_r_with_abs_pitch').value:+.2f}".replace("+", "") if False else
      f"r=-{abs(d('DLBS', 'classic', 'attenuation_r_with_abs_pitch').value):.2f}")
check("attenuation vs pitch p", a, pexp(num(TXT, "p_attenuation_abs_pitch", "DLBS", "classic")))
ap_ = d("DLBS", "classic", "attenuation_pose")
check("attenuation after pose", a, f"${-ap_.value:.1f}\\%$ per decade (${-ap_.hi:.1f}$ to ${-ap_.lo:.1f}$)")
ps = d("DLBS", "classic", "pose_share_of_axis_pct")
check("pose share of axis", a, f"${ps.value:.0f}\\%$ (${ps.lo:.0f}$ to ${ps.hi:.0f}\\%$)")
sh = d("HCP-A", "classic", "axis_share_pct")
check("HCP-A axis share", a, f"${sh.value:.0f}\\%$ (${round(sh.lo):.0f}$ to ${sh.hi:.0f}\\%$)", window=1400)
check("abstract axis share is in the Discussion", "The excess was attenuation by axis error that followed head pitch",
      f"${d('DLBS', 'classic', 'axis_share_pct').value:.0f}\\%$")
a = "Adjustment for head pose can also be expressed"
for i, lab in (("classic", "classic"), ("cross", "refined"), ("v2_slab", "measured"), ("pv_perp", "ratio")):
    check(f"absorbed {lab}", a, f"{num(N, 'pct_absorbed_pose', 'DLBS', i):.1f}\\%")
check("absorbed beta raw", a, f"{num(N, 'age_beta_raw', 'DLBS', 'classic'):.3f}")
check("absorbed beta pose", a, f"{num(N, 'age_beta_pose', 'DLBS', 'classic'):.3f}")
check("perm p95", a, f"{num(N, 'perm_p95', 'DLBS', 'classic'):.1f}\\%")
check("perm max", a, f"{num(N, 'perm_max', 'DLBS', 'classic'):.1f}\\%")
check("HCP-A perm p", a, f"p={num(N, 'perm_p', 'HCP-A', 'classic'):.2f}")
check("ratio vs pitch", a, f"r=-{abs(num(TXT, 'r_ratio_abs_pitch', 'DLBS', 'pv_perp')):.2f}")
check("ratio vs pitch given age", a, f"r=-{abs(num(TXT, 'r_ratio_abs_pitch_given_age', 'DLBS', 'pv_perp')):.2f}")

# ---- 4.3 repositioning ----
a = "Within a participant the anatomy is fixed"
for idx, lab in (("Classic", "classic"), ("Refined", "refined"), ("Measured axis", "measured"), ("lambda2/lambda3", "ratio")):
    s = REP.loc[idx]
    check(f"slope {lab}", a, f"{s.slope_pct_per_deg:.2f}\\%")
    check(f"slope CI {lab}", a, f"{s.ci_lo:.2f}$ to ${s.ci_hi:.2f}")
check("pairs", a, f"${int(REP.loc['Classic'].n_pairs)}$ visit pairs from ${int(REP.loc['Classic'].n_participants)}$")
check("theta ICC HCP-A", a, f"{num(TXT, 'icc_theta_scr', 'HCP-A', 'tract'):.2f}")
check("theta ICC DLBS", a, f"{num(TXT, 'icc_theta_scr', 'DLBS', 'tract'):.2f}")

# ---- 4.4 bound ----
a = "gives the fraction of $\\lambda_2/\\lambda_3$ that each index attained"
for i in ("classic", "cross", "v2_slab", "ALPS-PAS", "LD-ALPS", "reoriented"):
    check(f"attainment {i}", a, f"{num(N, 'attainment', 'DLBS', i):.2f}")
a = "Across participants the measured axis tracked"
check("r ratio v2 HCP-A", a, f"{num(N, 'r_with_ratio', 'HCP-A', 'v2_slab'):.2f}")
check("r ratio v2 DLBS", a, f"{num(N, 'r_with_ratio', 'DLBS', 'v2_slab'):.2f}")
for i in ("classic", "cross"):
    for c in ("HCP-A", "DLBS"):
        check(f"r ratio {i} {c}", a, f"{num(N, 'r_with_ratio', c, i):.2f}")
a = "In DLBS the classic index never exceeded it"
check("exceed refined", a, f"{num(N, 'pct_sessions_above_bound', 'DLBS', 'cross'):.1f}\\%")
check("exceed v2", a, f"{num(N, 'pct_sessions_above_bound', 'DLBS', 'v2_slab'):.1f}\\%")
check("exceed HCP-A classic", a, f"{num(N, 'pct_sessions_above_bound', 'HCP-A', 'classic'):.1f}\\%")
check("exceed HCP-A v2", a, f"{num(N, 'pct_sessions_above_bound', 'HCP-A', 'v2_slab'):.1f}\\%")
check("tract departure HCP-A scr", "In HCP-A the projection and association directions departed",
      f"{num(TXT, 'median_theta_scr', 'HCP-A', 'tract'):.1f}^{{\\circ}}")
check("interfiber HCP-A", "The two tracts were", f"{num(TXT, 'median_theta_interfiber', 'HCP-A', 'tract'):.1f}^{{\\circ}}")
check("interfiber DLBS", "The two tracts were", f"{num(TXT, 'median_theta_interfiber', 'DLBS', 'tract'):.1f}^{{\\circ}}")

# ---- 4.5 partialling ----
a = "In DLBS only the measured axis lost its association"
check("v2 raw DLBS", a, f"{num(N, 'age_r_raw', 'DLBS', 'v2_slab'):.3f}")
check("v2 partial DLBS", a, f"{num(N, 'age_r_partial_ratio', 'DLBS', 'v2_slab'):.3f}")
for i in ("classic", "cross", "ALPS-PAS", "LD-ALPS", "reoriented"):
    check(f"partial DLBS {i}", a, f"{num(N, 'age_r_partial_ratio', 'DLBS', i):.3f}")
a = "In HCP-A the associations of the refined index"
for i in ("cross", "v2_slab", "ALPS-PAS", "classic", "LD-ALPS"):
    check(f"partial HCP-A {i}", a, f"{num(N, 'age_r_partial_ratio', 'HCP-A', i):+.3f}".replace("+-", "-"))
a = "In DLBS the refined index was less strongly associated"
check("williams refined", a, f"p={num(TXT, 'williams_p_vs_classic', 'DLBS', 'cross'):.3f}")
lo = min(num(TXT, "williams_p_vs_classic", "DLBS", i) for i in ("v2_slab", "pv_perp"))
check("williams v2/ratio", a, f"p\\ge{math.floor(lo * 100) / 100:.2f}")

# ---- 4.6 reliability ----
a = "In DLBS reliability was lower for every index"
for i in ("classic", "reoriented", "pv_perp", "v2_slab", "ALPS-PAS", "cross"):
    check(f"ICC DLBS {i}", a, f"{num(N, 'icc', 'DLBS', i):.2f}")
rp = lambda c, i: v(REL[REL["between"].notna()], cohort=c, index=i)
a = "The higher reliability of the classic index in DLBS came from axis error that repeats"
check("atten within classic", a, f"{1e3 * rp('DLBS', 'classic').atten_within:.1f}")
check("atten within refined", a, f"{1e3 * rp('DLBS', 'cross').atten_within:.1f}")
check("atten between classic", a, f"{1e3 * rp('DLBS', 'classic').atten_between:.1f}")
check("atten between refined", a, f"{1e3 * rp('DLBS', 'cross').atten_between:.1f}")
hp = REL[REL.habitual_pitch_r.notna() & (REL.cohort == "DLBS")].set_index("index").habitual_pitch_r
check("habitual classic", a, f"r=-{abs(hp['classic']):.2f}")
check("habitual refined", a, f"r=-{abs(hp['cross']):.2f}")
check("pitch ICC", a, f"{rp('DLBS', 'abs_pitch').icc:.2f}")
check("ratio ICC in Discussion", "Its reliability was", f"${num(N, 'icc', 'DLBS', 'pv_perp'):.2f}$ in DLBS against ${num(N, 'icc', 'HCP-A', 'pv_perp'):.2f}$", window=80)

# ---- 4.7 sensitivity ----
a = "The left hemisphere gave the higher index in both cohorts"
for c in ("HCP-A", "DLBS"):
    for i in ("classic", "refined"):
        check(f"hemisphere {c} {i}", a, f"{num(TXT, 'left_minus_right_pct', c, i):.1f}\\%")
check("hemi williams L HCP-A", a, f"{num(TXT, 'williams_p_L', 'HCP-A', 'classic_vs_refined'):.2f}")
check("hemi williams R HCP-A", a, f"{num(TXT, 'williams_p_R', 'HCP-A', 'classic_vs_refined'):.2f}")
check("hemi classic R DLBS", a, f"{num(TXT, 'age_r_classic_R', 'DLBS', 'classic'):.3f}")
check("hemi refined R DLBS", a, f"{num(TXT, 'age_r_refined_R', 'DLBS', 'refined'):.3f}")
check("hemi williams L DLBS", a, f"p={num(TXT, 'williams_p_L', 'DLBS', 'classic_vs_refined'):.2f}")
fa = {c: {m: num(TXT, f"fa_floor_max_change_{m}", c, "all") for m in ("attainment", "icc", "age", "beyond")}
      for c in ("HCP-A", "DLBS")}
a = "Repeating the analyses with fractional anisotropy floors"
for m in ("attainment", "icc", "age", "beyond"):
    check(f"FA {m}", a, f"{max(fa['HCP-A'][m], fa['DLBS'][m]):.3f}")

# ---- generated tables, verbatim ----
# With the manuscript here, each generated table must appear in it word for
# word. Without it, what can still be checked is that the generator
# reproduces the shipped tables from the shipped numbers, which is the claim
# the code repository makes.
gen_path = HERE / "r2_tables.tex"
gen = gen_path.read_text(encoding="utf-8").strip().split("\n\n")
if HAVE_TEX:
    for block in gen:
        label = re.search(r"\\label\{([^}]+)\}", block).group(1)
        count += 1
        if block.strip() not in TEX:
            fails.append(f"table {label} differs from r2_tables.tex")
else:
    import shutil as _sh, subprocess as _sp, tempfile as _tf
    _keep = gen_path.read_text(encoding="utf-8")
    _bak = Path(_tf.mkdtemp()) / "r2_tables.tex"
    _sh.copy2(gen_path, _bak)
    _r = _sp.run([sys.executable, str(HERE / "r2_latex_tables.py")],
                 capture_output=True, cwd=str(HERE))
    count += 1
    if _r.returncode:
        fails.append("r2_latex_tables.py failed: " + _r.stderr.decode()[-200:])
        _sh.copy2(_bak, gen_path)
    elif gen_path.read_text(encoding="utf-8") != _keep:
        fails.append("r2_latex_tables.py does not reproduce r2_tables.tex")
        _sh.copy2(_bak, gen_path)
    else:
        print(f"{len(gen)} tables regenerate from the shipped numbers, "
              "byte for byte")

# ---- author names, against the ORCID records they are submitted with ----
# The manuscript, the title page and the cover letter each carry the author
# list separately, and nothing kept them in step. "Jonathon" was spelled
# "Jonathan" in some of them. A co-author's name is not something to get
# wrong, so the spelling is checked wherever a live document states it.
ORCID_NAMES = {"0000-0003-2401-0345": "Jonathon K. Maffie",
               "0009-0009-9493-0957": "Rafael Kuc",
               "0000-0002-6496-2087": "Scott N. Hwang",
               "0000-0002-4192-7414": "Sangam G. Kanekar"}
LIVE = ("mri_revision.tex", "mri_title_page.tex", "revision_cover_letter.tex",
        "revision2_cover_letter.tex", "response_to_reviewers_r2.tex")
for _f in LIVE:
    _p = HERE.parent / _f
    if not _p.exists():
        continue
    _t = " ".join(_p.read_text(encoding="utf-8", errors="ignore").split())
    for _orcid, _name in ORCID_NAMES.items():
        _sur = _name.rsplit(" ", 1)[1]
        if _sur not in _t:
            continue
        count += 1
        if _name not in _t:
            fails.append(f"{_f}: {_sur} is not spelled '{_name}' (ORCID {_orcid})")
    # and the misspelling specifically, in case both forms appear
    count += 1
    if "Jonathan" in _t:
        fails.append(f"{_f}: spells Maffie's first name Jonathan, not Jonathon")

if not HAVE_TEX:
    print("mri_revision.tex is not here, so the checks that read the "
          "manuscript text were skipped.")
    print("Those run in the manuscript repository. The table "
          "regeneration below does not need it.")
print(f"{count} checks, {len(fails)} failed")
for f in fails:
    print("  FAIL", f)
sys.exit(1 if fails else 0)
