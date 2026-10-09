r"""Every second-round reviewer point, and whether the paper carries the change.

The reply says 14 things were done. A reply is easy to write and a change is
easy to forget, so each claim is checked against mri_revision.tex rather than
taken on the reply's word. Where a claim is a judgement that no string can
settle, it is listed as needing an eye rather than reported as passing.

    python check_reviewer_coverage_r2.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEX = (ROOT / "mri_revision.tex").read_text(encoding="utf-8")
FLAT = " ".join(TEX.split())
AUX = (ROOT / "mri_revision.aux").read_text(encoding="utf-8", errors="ignore")
B = chr(92)

rows: list[tuple[str, str, bool | None, str]] = []


def check(who: str, point: str, ok: bool | None, what: str) -> None:
    rows.append((who, point, ok, what))


def section(title: str) -> str:
    """The source of one numbered unit, by its title."""
    m = re.search(re.escape(B) + r"(sub)?section\{" + re.escape(title) + r"\}", TEX)
    if not m:
        return ""
    nxt = re.search(re.escape(B) + r"(sub)?section\{", TEX[m.end():])
    return TEX[m.end():m.end() + (nxt.start() if nxt else len(TEX))]


def label(name: str) -> str:
    m = re.search(re.escape(B) + r"newlabel\{" + re.escape(name) + r"\}\{\{([^}]*)\}", AUX)
    return m.group(1) if m else ""


# ---------------------------------------------------------------- Reviewer 1
# No further comments. The reply lists where the first-round analyses now sit,
# and a wrong pointer there is the kind of thing a reviewer checks first.
for what, title, needle in (
        ("reliability in Section 4.6", "Reliability", "ICC"),
        ("hand-drawn comparison in Section 3.2", "Measurement regions", "hand"),
        ("failure rates in Section 3.1", "Cohorts and quality control",
         "completed processing"),
        ("inter-tract angle in Section 4.4",
         "Each index approaches the radial anisotropy as its axis error falls",
         "The two tracts were"),
        ("FA threshold and covariates in Section 4.7", "Sensitivity analyses",
         "fractional anisotropy floors")):
    body = section(title)
    check("R1", "pointers", bool(body) and needle.lower() in body.lower(), what)

# ---------------------------------------------------------------- Reviewer 3
check("R3", "0", B + "appendix" not in TEX and "Supplementary material" not in TEX,
      "the supplement is gone")

sec2 = section("Theory")
check("R3", "1", bool(re.search(re.escape(B) + r"subsection\{Notation\}", TEX)),
      "Section 2 opens with the notation")
check("R3", "1", "begin{proposition}" in TEX and "begin{proof}" in TEX,
      "the bound is a proposition with a proof")
# cas-sc numbers \begin{equation}; an unnumbered display is \[ ... \]
theory = TEX[TEX.index(B + "section{Theory}"):TEX.index(B + "section{Methods}")]
check("R3", "1", theory.count(B + "[") == 0, "no unnumbered display in Section 2")

notation = section("Notation")
check("R3", "1a", "rangle" in notation and "langle" in notation,
      "the region mean is defined in Section 2.1")
check("R3", "1a", label("eq:rhobar") == "7",
      "the pooling bound is Eq. 7")
check("R3", "1b", bool(re.search(r"v\}_2.{0,200}eigenvector|eigenvector.{0,200}v\}_2",
                                 notation, re.S)),
      "v2 and v3 are defined as eigenvectors in Section 2.1")
check("R3", "1c", label("eq:axes") == "3",
      "the axis error is defined at Eq. 3")
check("R3", "1c", "We call $1-A/" + B + "rho$ the shortfall" in FLAT,
      "the shortfall is defined")
check("R3", "1c", label("eq:expansion") == "6",
      "the expansion is Eq. 6")
check("R3", "1c", B + "frac{" + B + "rho^2-1}{" + B + "rho}" in FLAT,
      "the second-order term is stated")

intro = section("Introduction")
check("R3", "2", "correct" in intro.lower(),
      "the Introduction introduces the corrected indices")
s34 = section("Tract directions and the corrected indices")
check("R3", "2", "deliberately simple" in s34, "Section 3.4 says the correction is simple")

abstract = TEX[TEX.index(B + "begin{abstract}"):TEX.index(B + "end{abstract}")]
check("R3", "3", "widely used" in abstract or "used for" in abstract,
      "the abstract says what the index is used for")
check("R3", "3", "fixed in the scanner" in abstract or "axes fixed" in abstract,
      "and why fixed axes make it depend on head position")

for frag in ("where preprocessing has removed head position, it does not",
             "variants are indistinguishable"):
    check("R3", "4", frag not in FLAT, f"removed: {frag[:44]}")

check("R3", "5", "report white-matter geometry more than fluid clearance" not in FLAT,
      "the fluid-clearance phrase is gone")
disc = section("What DTI-ALPS measures")
check("R3", "5", "diffusivit" in disc.lower(),
      "the Discussion states it in terms of diffusivities")

check("R3", "6", "Wilcoxon" not in FLAT or "signed-rank" in FLAT,
      "no bare Wilcoxon reference")
check("R3", "6", "per-voxel" not in FLAT.lower(),
      "the per-voxel comparison is gone")
s35 = section("Statistical analysis")
check("R3", "6", bool(s35.strip()), "Section 3.5 exists and names the tests")

check("R3", "7", "hyperkyphosis" in intro.lower(),
      "hyperkyphosis is in the Introduction")
results = TEX[TEX.index(B + "section{Results}"):TEX.index(B + "section{Discussion}")]
check("R3", "7", "hyperkyphosis" not in results.lower(),
      "and no longer in Results")

check("R3", "8", None, "interpretation moved out of Results (needs an eye)")

# ---------------------------------------------------------------- Reviewer 4
check("R4", "1", label("prop:bound") == "1", "Proposition 1 is numbered 1")
check("R4", "1", label("eq:expansion") == "6", "the expansion is Eq. 6")
for num in ("32.9", "-0.450", "-0.302", "12.8", "+0.369", "-0.543", "-0.430"):
    check("R4", "1", num.lstrip("+") in FLAT, f"the paper states {num}")


def main() -> None:
    fail = [r for r in rows if r[2] is False]
    eye = [r for r in rows if r[2] is None]
    for who, pt, ok, what in rows:
        mark = "ok  " if ok else ("EYE " if ok is None else "FAIL")
        print(f"  {mark} {who}.{pt:<3s} {what}")
    print(f"\n  {len(rows) - len(fail) - len(eye)}/{len(rows)} verified, "
          f"{len(eye)} need an eye, {len(fail)} failed")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
