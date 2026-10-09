# paper-DTI-correction-code

Every script that produced a number, table or figure in *Head Position
Confounds DTI-ALPS, While Its Correction Approaches Radial Anisotropy*
(Magnetic Resonance Imaging, MAGRESIMAGING-D-26-00371), together with the
derived values that can be redistributed.

These are the working scripts, not a packaged tool. They are published so that
any result in the paper can be traced to the code that made it. Expect to edit
paths before anything that reads images will run.

## What runs in a fresh clone

Every table in the paper, and the check on every number in it:

    uv venv
    uv pip install -r requirements.txt
    cd scripts
    python r2_latex_tables.py     # regenerates the four tables as r2_tables.tex
    python verify_r2.py           # confirms they regenerate byte for byte

`verify_r2.py` reads the manuscript when it is beside it and then checks every
number the text states against the CSV it came from. The manuscript is not
published here, so in this repository it checks the tables instead and says
which checks it skipped.

`scripts/sorting_bias_floor.py` also needs nothing: it simulates the
acquisition and reproduces the eigenvalue-sorting noise floor quoted in the
paper.

## What needs the images

There is a chain: source images, then processing, then per-session CSVs, then
the analyses that read those CSVs. The per-session CSVs for DLBS are here. The
ones for HCP-A are not, because they carry participant identifiers and the AABC
Data Use Terms do not permit redistributing them. `scripts/WITHHELD_DATA.txt`
names the six files and what regenerates them.

So of the second-round analyses, `r2_latex_tables.py`, `r2_reorientation.py`
and `r2_repositioning.py` run here, and the other seven need your own HCP-A
access. Scripts that read images need the cohorts below, FSL, and the processed
sessions on disk.

## Data

- **HCP-A** — the BALSA repository, <https://balsa.wustl.edu>, to qualified
  registered users under the AABC Consortium Data Use Terms. Identifier-linked
  derived values may not be redistributed under those Terms, which is why none
  are here.
- **DLBS** — OpenNeuro, <https://openneuro.org/datasets/ds004856>, version 1.2.0.
- **Trigeminal neuralgia** — OpenNeuro,
  <https://openneuro.org/datasets/ds005713>, version 2.0.2, CC0. Used in
  earlier work, not in this revision. Its scripts are kept for the record.

## Install

    uv venv
    uv pip install -r requirements.txt

Install **FSL** separately, <https://fsl.fmrib.ox.ac.uk>. It is not a Python
package. Several scripts shell out to it for registration and the tract atlases.

## Data directories must be edited

Paths are hardcoded to the machine the analyses ran on: `Q:/dti_output` for the
processed sessions, `M:/` for the trigeminal release, `F:/` and `P:/` for the
HCP-A structural packages, `C:/Users/.../fsl` for FSL.

`scripts/data_paths.py` translates those between Windows and WSL, but the drive
letters are literals. Change them to your own locations first.

    grep -rn "winpath\|:/" scripts/*.py

## Scripts that take a cohort

Some scripts analyse more than one region set and take a flag, because running
them without one silently gives you whichever cohort the default names. Outputs
are suffixed by the choice so two runs cannot overwrite each other.

    python scripts/longitudinal_reliability.py --cohort manual   # hand-drawn regions
    python scripts/longitudinal_reliability.py --cohort auto     # atlas placement
    python scripts/longitudinal_reliability.py --cohort spheres  # 5 mm spheres

This one produces the between-visit reliability and the repositioning
sensitivity, including the change in the index per degree of change in the
scanner-to-anatomy angle. `--cohort manual` writes the unsuffixed files;
the others append `_auto` or `_spheres`.

`scripts/reliability_analysis.py` takes `--input` and `--label` instead, and
`scripts/manual_vs_atlas_icc.py` computes both region sets on one estimator so
the two are directly comparable.

## Three things that fail quietly

Each of these produced a wrong answer that looked like a clean result, so they
are worth knowing before you trust anything you run.

**FSL through WSL has no `FSLDIR`.** On Windows the FSL commands run via
`wsl bash -lc`, and that non-interactive login shell does not export `FSLDIR`.
A `flirt` given an unresolvable `-ref` **exits zero and writes no matrix**, so
the failure surfaces later as an empty result rather than an error. Resolve
template paths at runtime and raise when they are missing;
`scripts/ds001907_t1_pose.py` shows the pattern.

**HCP-A is two shells and the paper uses one.** It carries b=1500 and b=3000,
the tensors are fitted on b=1500 and written with a suffix, and the published
sample is the 1706 sessions that have them, not the 2742 that passed region
placement. Select the shell with the environment variable, and the sample from
the published table:

    ALPS_TENSOR_SUFFIX=_b1500 python scripts/measured_pvs_axis.py --cohort hcpa

`alps_variants` in `scripts/tn_alps.py` honours the same variable. Getting
either wrong is invisible in the output, which is why every script that
recomputes a published quantity also recomputes the classic index in the same
pass and checks it reproduces its printed value before the new numbers are
used.

**`atomic_io` breaks `to_csv(mode="a")`.** Importing it replaces `to_csv` with a
write-to-temp-and-rename, which is right in general and silently ignores append
mode: each "append" rewrites the file with only the new rows and drops the
header. Accumulate by read, concatenate, write instead.

**Partialling a variable out of itself returns an arbitrary number.** `pv_perp`
is the eigenvalue ratio, so in any analysis that adjusts for the ratio its
residual is floating-point residue rather than signal. Correlating that residue
against age returns whatever the rounding happens to give, and at n=809 it can
clear p<0.05. It produced +0.081 in one script and +0.081 with −0.087 flagged
significant in another, either of which reads as a variant retaining something
beyond the ratio. Guarding on an exactly zero residual does not catch it,
because the residue is never exactly zero. Compare the residual standard
deviation against the variable's own instead, and report a collapsed residual
as undefined. `scripts/beyond_ratio_adjusted.py` shows the guard.

## Region placement

The regions are 5 mm spheres drawn in native space at the centre of the warped
template mask. This is the primary analysis and the default, so no environment
variable is needed to reproduce the manuscript.

The first submission measured inside the warped mask itself, which arrives
distorted in both size and shape. Region size varied 7.8-fold across HCP-A and
3.1-fold across DLBS, which is why the age models carried a region-volume
covariate. Redrawing the sphere holds size to within 6% and 15% respectively,
so the covariate is no longer load-bearing. Registration still decides where
the region sits, which is what it is good at, and no longer decides how big it
is or what shape.

Set `ALPS_SPHERE_MM=0` to restore the warped masks. Those runs write to
`_warpedmask` filenames, so neither placement can overwrite the other:

    ALPS_SPHERE_MM=0 ALPS_TENSOR_SUFFIX=_b1500 \
        python scripts/measured_pvs_axis.py --cohort hcpa --all-sessions

`scripts/resphere_impact.py` reports both placements side by side for every
quantity the manuscript states. `alps_variants` in `scripts/tn_alps.py` reads
the same variable, so the comparators cannot drift to a different placement
from the variants they are tabulated beside.

## Running

Each script is self-contained and documented at the top of the file, which is the
best starting point for tracing a single number.

`scripts/regenerate_and_diff.py` runs the chain and diffs each output against the
committed version, reporting any column that moved. It is 41 steps and roughly
three hours, starting from the two index tables that everything else reads:

    python scripts/regenerate_and_diff.py --dry-run    the plan and its timings
    python scripts/regenerate_and_diff.py              run all of it
    python scripts/regenerate_and_diff.py --only ratio one group
    python scripts/regenerate_and_diff.py --skip-slow  omit the long steps

The groups run in dependency order, so `index` first and the rest after it.
Roughly 90 minutes of the total is the four `index` steps, which read images.
Everything downstream of them reads CSVs and takes about a minute each.

The repository holds more scripts than the chain calls. The extra ones are
figure builders, exploratory analyses, and work that did not reach the final
manuscript. Every script behind a number in the paper is in the chain, and each
one is documented at the top of its own file, which is the best place to start
when tracing a single value.

`scripts/verify_r2.py` checks every number the manuscript states against those
outputs, in the sentence that states it. It needs the manuscript source, which is
not published here. Without it, it checks only that the tables regenerate. The
first revision's verifier, and the two tools that maintained it, were removed
because they apply only to the first revision's text.

## License

MIT. See `LICENSE`.
