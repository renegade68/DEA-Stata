# DEA-Stata

Public regression tests work from a Git clone or source ZIP. From this directory,
run `python3 -B tests/run.py --suite smoke`; see [test setup](tests/README.md) for
Python dependencies and Stata discovery. Personal settings and generated test
results are excluded from the public inventory. [Publication checks](tests/PUBLICATION.txt)
cover paths, credentials, staged files and commit email metadata.

`dea_crosseff` 1.10.4 provides failure-only certified recovery in original data
units for explicit simplex. Accepted evaluator weights and self scores remain
unchanged. Recovered rows use directly exported raw weights and bounded exact
peer arithmetic, including input VRS Lim-Zhu. Unverified primary solutions stay
unavailable; a failed secondary retains its certified primary weights.
Version 1.10.4 adds exact cancellation and common-factor division for failed
raw secondary certificates. It discovers a positive integer objective factor
in Mata and keeps the original objective units with the certificate. The extra
attempt shares prior work and pivot limits. Recovery path 5 identifies this
attempt; its outcome distinguishes a finite optimum, a verified improving ray,
and an unresolved result. Existing accepted solutions are preserved.

Version 1.10.3 allows one bounded exact primary fit after both numerical
proposals fail, using a valid lossless representation. Existing accepted weights
retain their selection. Finite input data can still require original-unit weights
outside Stata's numeric range; unavailable weights do not establish mathematical
infeasibility. Scaled-weight output is a separate research topic.

Version 1.10.2 factors common powers of two in selected exact arithmetic
checks while preserving the existing width, work and tolerance limits. Additional
secondary optima and improving rays can be certified; an exactly unbounded
secondary retains certified primary weights. Other limited attempts remain
unknown. Run `python3 -B tests/run.py --suite crosseff_dyadic --suite crosseff_raw`
to check the portable arithmetic and recovery references.

Version 1.10.1 extends finite coefficient decoding only for failed raw
primaries and their qualifying secondaries, with the same arithmetic limits.
Matrix diagnostics distinguish recovery, denominator policy and arithmetic
limits. Native IPM remains the default. See
[the help](ado/dea_crosseff.sthlp) and [installation](INSTALL_CROSSEFF.txt).

`dea_crosseff` 1.9.0 adds explicit `cerange`: bounded exact peer CE ranges
on the original primary/secondary optimal faces. It includes `optimadiag`
weight/denominator diagnostics and retains existing scores, weights and engine.
Finite attained bounds, limit-only bounds, unbounded ratios, zero denominators,
unknown results and stored-weight linkage have separate states. Mean/rank ranges
and a weight-selection policy remain separate work. See
[the help](ado/dea_crosseff.sthlp), [installation](INSTALL_CROSSEFF.txt), and
[CE-range verification](tests/CROSSEFF_CERANGE.txt).

Cross-efficiency update (`dea_crosseff` 1.9.0): explicit simplex now uses
bounded original-data optimal-face certificates, failure-only exact recovery,
and per-stage recovery/cost diagnostics. Version1.6.1 checks both native
stages for nonfinite prepared arrays before optimizer setup and reports
unavailable stages/partial peer means with warnings and observation-level
diagnostics. A skipped secondary retains primary weights; a skipped primary
stays unavailable. Native IPM remains the default.

Version1.7.0 builds secondary objectives from peer-only sums. Explicit
simplex certificates check the exact intended original-data sum, including
coefficient contributions lost by floating aggregation. Native objectives
remain numerical-only; rounding and fallback are reported separately.
The primary calculation and CRS/VRS peer conventions are unchanged. See
[objective regression](tests/CROSSEFF_OBJECTIVE.txt).

Crossgame update (`dea_crossgame` 2.4.0): failed secondary/game solves can
use bounded original-data certificates. Diagnostics distinguish native failures,
certified recovery, unresolved cases and complete Nash iteration acceptance.
Original self scores and accepted native solutions are preserved. Failed secondary
solves use the certified primary optimal face; game recovery retains the same theta.
Automatic helper loading and temporary outputs preserve existing result columns
on compute errors. See [installation and verification](INSTALL_CROSSGAME.txt) and
[the recovery test contract](tests/CROSSGAME_RECOVERY.txt).

Install the complete matched 92-file `ado/` directory; see
[the cross-efficiency installation guide](INSTALL_CROSSEFF.txt).

Package update (`dea` 2.0.73): the explicit Cooper simplex path now uses
bounded exact certificates and optional reference selection. Native IPM shares
the objective-weight preflight. See the Cooper section and INSTALL_COOPER.txt.
The simplex path returns `r(integration_status)="integrated"`; read the
certificate and generated status fields for each DMU's availability.

Development update (`dea` 2.0.72): optional `recoverylog` reports common
simplex calls, recovery attempts and final model availability for ordinary TE,
SBM, super-SBM and joint weak-disposal DDF. Recording is off by default.
See `help dea` and [the installation and verification guide](INSTALL_RECOVERYLOG.txt).

Ordinary simplex TE now checks the radial factor passed from Stage 1 to
Stage 2. A bounded retry is permitted only after arithmetic error intervals
prove a negative component in the original-row basis solution. A corrected
factor must pass the change budget and original primal/dual quality checks.
The retry uses min(tolerance(),1e-10) for that DMU; the public default tolerance
is unchanged. Stage 2 rebuilds its LP from the same original-data-derived,
column-normalized inputs and outputs, with a fresh basis and the accepted
factor. A failed corrected Stage 2 retains the score and withholds references.
Stage(1) uses the same radial decision and skips slack optimization.

r(stage_linkage) records the before/after factors, reason, retry tolerance and
Stage-2 quality for the main ordinary simplex TE calculation. It is available
without recoverylog or reference; VRS auxiliary calls do not replace its rows.
r(numerical_diagnostics)'s radial_initial is the accepted factor used to build
Stage 2, so use stage_linkage's radial_before for the value before this retry.
Factors are theta (input) and phi (output); output efficiency is 1/phi, subject
to the existing near-one presentation cleanup. The matrix is not added to
saving() files. The IPM projection-recovery contract remains separate.

The basis check is bounded by m<=32,n+2m<=4096,m*n<=16384 and arithmetic-range
guards. Here m and n are the original Stage-1 LP row and real-variable counts.
It requires a complete original-row basis; unsupported/uncertain cases retain
the prior result. The dual fallback is not covered by this correction. Shared radial TE
subproblems also use the guard, but only ordinary model(te) exposes this
public diagnostic matrix.
One fresh retry is allowed; existing solver-internal recovery caps still apply.
This does not guarantee a feasible exact fixed target for every dataset, nor
remove all tolerance sensitivity or nonuniqueness of reference weights.

See [the stage-link regression note](tests/STAGE_LINKAGE.txt).

Successful simplex solves now receive an additional basis-quality check.
Already-certified results are preserved; a verified defect may be corrected
only after the refined point passes original primal and primal-dual checks.
The public LP adapter uses the original data units, bounds and objective
constant, including when `scaling(row)` is requested. If certification is
unavailable or refinement cannot be verified, the original success is retained.
Thus, a successful return is not a universal certificate of exact optimality.
A conservative fast screen falls back to compensated checks when needed.
No new public option is required. See [the quality-check regression note](tests/SUCCESS_QUALITY.txt).
The failure-only recovery descriptions below retain their ordering; their
successful outputs also undergo this additional check.

The simplex engine also supports bounded exact selection of independent
original rows, including the right-hand side, and analytic checks for a trivial
empty basis. After every earlier recovery path fails, it may solve a reduced
row system once; acceptance still requires checks against all original rows,
columns and bounds. Near dependence alone never authorizes row removal.
The selector is limited to256 rows/columns,rank32,32 expansion components and
one million work units, with additional arithmetic-range guards. Unsupported
cases retain the existing outcome. No new public option or default strict
policy is introduced. These finite checks do not guarantee exact optimality
for every dataset. See tests/ROW_RECOVERY.txt for scope and reproducible tests.

The simplex engine (2.0.25) now checks both basis-inverse identities with
compensated residuals, using a conservative error bound for fast acceptance
when safe. The internal inverse threshold and user tolerances are unchanged.
Difficult cases can change availability or select alternative optimal reference
weights; the checks do not guarantee universal numerical reliability.
See [the numerical verification note](tests/INVERSE_VALIDATION.txt).

The common simplex solver now retries a numerical failure with a factorized
basis and checks the result against the original LP. Initially accepted results
proceed to the successful-basis quality check described above. If these paths still fail, an artificial-free retry can refine the
original basis solution and check its primal-dual quality with arithmetic error
bounds. Uncertain candidates stay unavailable; user tolerances are unchanged.
This does not promise relative objective accuracy for every accepted result.
See [the precision-recovery verification note](tests/PRECISION_RECOVERY.txt).

If all earlier paths fail after a blocked Phase-1 pivot, a final numerical retry
can reconstruct a basis from original columns. Original rows and the final
primal/dual quality checks remain required. Earlier accepted solutions and
reference weights reach the successful-basis quality check first. Failed retries can add iterations; every
attempt is included in the reported total. See the
[reconstruction regression note](tests/RECONSTRUCTION.txt) for exact-witness
controls and dimension limits.


If every preceding numerical path fails, an additional bounded retry can remove
exact positive power-of-two duplicate rows or columns. Rows include their RHS;
columns include their cost and every original constraint, including bound rows.
The result is restored to the original variables and checked against all original
constraints and primal/dual quality requirements. Existing accepted solutions
and reference weights reach the successful-basis quality check first. This does not merge integer variables;
MILP branching and integrality checks use the restored original variables.
Unsupported scales, uncertain certificates and degenerate supports can still
return numerical failure. Tolerances are unchanged. See
[the duplicate-recovery regression note](tests/DUPLICATE_RECOVERY.txt).
This priority applies within the common simplex solver. With `engine(auto)` or
`lpengine(auto)`, a newly recovered simplex result can replace a previously
successful native fallback; scores, reference weights, iteration totals and the
reported backend can therefore change across versions.



Repeated Phase-2 Bland states can now end a numerical path early. Existing
restarts keep priority, preserving their choice among alternative optimal
weights. If all preceding paths fail after an observed repeat, an additional
retry uses reduced-cost intervals that include basis-solve error. Uncertain
prices remain unresolved. This does not change user tolerances or certify the
accuracy of all previously accepted results. See
[the certified-pricing regression note](tests/PRECISION_PRICING.txt).


Failure-only transitions now guard trial basis exchanges and can reconstruct a
real-column basis while retaining every original row. Small QR selection is
preserved; wider problems use thin Householder column selection. The separate
precision screen still has the standard-form limit `n+2m<=4096`.

The optional standalone LP/MILP `scaling(row)` / `lpscaling(row)` policy now uses
power-of-two row and objective multipliers with exact coefficient roundtrip.
A transformation that loses a stored coefficient is skipped. The default
`none` is unchanged. Changing this initial representation can change previously
accepted numerical results; it is distinct from failure-only recovery and is
not a guarantee that every result becomes more accurate. See
[the wide-recovery and scaling tests](tests/WIDE_SCALING.txt).
A separate conservative dual-bound recovery is limited to the
VRS/input-direction weak-disposal model. Previously missing results may become
available; these mechanisms do not guarantee that every numerical failure can
be recovered, and they do not change Stata's native IPM implementation.

Ordinary TE validates nonnegative quantities and at least one positive input
and output in every selected DMU before replacing results or saving. Individual
zero components and zero columns are allowed under that per-DMU condition.
Selection and missing required quantities determine the checked sample.
Earlier ordinary TE results from negative or all-zero input/output vectors
require review. Legacy IDEA transformations retain their own internal domain.
VRS super-efficiency can still have no mathematical solution.

Stata software for the article:

> Lee, C. Forthcoming. *dea: An Update and Extension to Data Envelopment Analysis*.
> Accepted for publication in the [Stata Journal](https://www.stata-journal.com/).

This repository maintains the article's DEA command family, subsequent updates,
and independent LP/MILP tools supporting the earlier optimization workflow.
The current development version may differ from the software submitted with the
article. The acceptance notice refers to the article; it does not extend the
article's validation claims to later code changes.

Requires Stata 16 or newer. The command index below defines the software scope.
The game-theoretic commands are exploratory diagnostics, as described in the
article, and are not fully specified structural game estimators. Command-specific
assumptions and limitations are documented in the linked help files.

`dea_crossgame` 2.4.0 validates and scales CRS input LPs, uses nonnegative
weights, and implements `matrix`/`detail` output with stage and iteration
results. Nash requires complete initialization and complete pairwise iterations;
a failed update retains the last complete state and reports failure, while an
iteration limit is distinguished from convergence. Aggressive/benevolent
strategies retain the documented weighted-output surrogates. Their scores are
not guaranteed to bracket Nash scores. See help for averaging and fallback rules.

`dea_crosseff` 1.10.4 keeps native IPM as the default. Explicit
`engine(simplex)` obtains the numerical primary formulation, then uses
original quantities for bounded exact primary/secondary certificates. The
secondary face uses positive exact primary dual constraints; a rounded
numerical self score is not inserted as a new equality RHS. The numerical
self score is preserved and checked against the certificate and exported
weights using the requested tolerance.

The exact path is limited to 16 variables after free-variable splitting and
32 selected observations. Equality-first and auxiliary-feasibility recovery
share 16 exact exchanges and finite arithmetic budgets per LP. A failed
secondary retains primary weights. After bounded recovery fails, an exact
feasible point and improving recession direction can prove unboundedness of
the original-data LP. An unverified proposal or missing numerical primary
does not establish infeasibility.

`matrix` returns weights, stage/solver/link diagnostics, normalization-loss
locations, and primary/secondary recovery journals. `detail` prints stage
outcomes and costs; summary counters are available without `matrix`. Internal
exact certificates, rounded exports and numerical self are distinct.
The separate failure-only raw path can recover unavailable numerical primaries;
normalization-loss diagnostics alone do not repair accepted evaluator rows.
Acceleration experiments allocating additional exchanges are excluded.

The legacy VRS peer ratio omits the intercept and remains the default.
`rts(vrs) vrscross(limzhu)` implements the Lim-Zhu input peer formula;
translated diagonal appraisal is distinct from ordinary self efficiency.
Undefined ratios remain missing. Existing secondary surrogates,
off-diagonal averaging and self/peer maverick ratio are unchanged. These
surrogates do not optimize corrected peer ratios. See help and
INSTALL_CROSSEFF.txt for mathematical scope and return definitions.

`dea_malmquist` 1.3.7 and `dea_window` 1.3.8 now apply `tolerance()` to native
LP solves and primal checks, and reject invalid tolerances or quantity domains
before replacing results. Malmquist resets summaries for period pairs with no
common DMUs and returns `r(pair_summary)` with component failure counts.
Window analysis reports attempted, successful and unavailable evaluations.
Available components and scores remain usable when other fits fail. The help
files distinguish the selected RTS index, pooled best-practice gap change,
natural cross-period infeasibility, and the weighting of window summaries.

In `dea` 2.0.58, SBM reduces the averaging divisor when omitting zero terms,
and SBM/super-SBM rescale columns internally to improve unit invariance.
`model(sbm)` permits nonnegative data with a positive input and output per DMU;
zero outputs use the documented structural-zero interpretation.
`model(supersbm)` requires strictly positive data and now applies explicit
input/output orientation. Its default remains non-oriented. Earlier versions
ignored the super-SBM orientation option; those results should be recomputed.
Oriented super-SBM can be mathematically infeasible under non-CRS technologies;
unavailable scores remain missing. See `help dea` for the precise conventions.

In `dea` 2.0.57, `model(allocative)` validates the required price list and
its dimensions before creating or replacing results. Input orientation requires
one `w()` variable per input; output orientation requires one `p()` variable per
output. Every supplied price list must match the corresponding quantities in
length and order. Earlier short or mismatched lists could cause computation
errors or missing scores without rejecting the price specification.
Valid specifications retain their calculation rules.

In `dea` 2.0.56, numeric interval IDEA uses a consistent realization of the
evaluated DMU in both its target and reference row. Both pessimistic and
optimistic results from earlier versions must be recomputed. With positive
ordered bounds, the mathematical scores lie between 0 and 1; a singleton has
both scores equal to 1. Point intervals reduce to ordinary DEA. These are
worst/best relative efficiencies over independently varying data intervals,
not confidence intervals or scores against one fixed frontier. The legacy
bounded/ordinal transformation is unchanged; its internal DEA calls now also
receive the requested `tolerance()` instead of silently using the default. See the
[interval formulation tests](tests/README.md#numeric-interval-idea-formulation).

In `dea` 2.0.55, core non-TE LP models forward `tolerance()` and check
finite solutions against their submitted constraints and bounds. Native
convergence keeps the stricter of the request and the previous model setting;
simplex retains its pivot safeguards. A converged native candidate that fails
the primal checks may get one stricter restart of the identical LP; an accepted
initial solution is retained. A rejected result remains missing.
This checks primal residuals, not a guarantee of objective accuracy or an
infeasibility certificate. SBM with undesirable outputs preserves every
positive denominator and normalizes columns with their slack units. Small
positive multiplier weights are retained; failed native weight solves return
missing weights, while a separately available radial score may remain usable.
Legacy Form A `tol1()` remains a deprecated standalone alias; combining it
with `tolerance()`, or supplying obsolete `tol2()`, now returns `r(198)`.

In `dea` 2.0.54, ordinary TE returns `r(numerical_diagnostics)` with the
initial and projection radial factors, original-unit and normalized residuals,
and a separate check against the requested tolerance. Both engines expose the
same diagnostics, including after automatic fallback or explicit radial recovery.
Numerical acceptance does not certify exact LP feasibility or optimality.
Accepted projections are not reoptimized or moved to another radial target.
Small positive simplex scores
are preserved instead of being rounded to zero by the old 1e-8 cleanup rule.
Reference/slack cleanup is checked against the original LP; if rounding breaks
the requested residual criterion, the original solver-accepted values are retained.

In `dea` 2.0.53, congestion LP failures remain missing instead of becoming
an efficiency of 1. Congestion and SBM-undesirable exclude missing values
from efficiency counts and publish actual failure counts. Congestion also
reports failures by component, retaining a usable main measure after a
CRS-only failure. Undefined ratios remain missing. Core `saving()` includes
all five congestion auxiliary columns or the SBM model's bad-output variables.
See the command help for the distinction between total and main-score failures.

## Start using the commands

Download or clone this repository. Keep **all files in `ado/` together**, including
the three implementation `.do` files. In Stata, add that folder for this session:

```stata
adopath ++ "/your/path/DEA-Stata/ado"
which dea
help dea
```

Replace the example path with the location on your computer. Run `which` to check
which copy Stata uses if you already have another DEA installation. This working
copy uses a local-folder installation; it does not yet provide `net install` metadata.
After replacing solver files, start a fresh Stata session so previously compiled
Mata routines do not retain references to the old solver.

Keep the complete `ado/` directory together when updating, including `milp_bb_impl.do`,
`_dea_output_check.ado`, `_dea_parse_io.ado`, `_dea_bootstrap_core.ado`,
`_dea_sw2007_seq.ado`, `_dea_recovery.ado`, and `_dea_weak_bound.ado`. All public commands reject existing result-variable
names by default. Use `replacevars` to explicitly refresh existing results;
model variables and identifiers remain protected. For `dea`, `replace` still
controls only the dataset written by `saving()`. Repeated VRS calculations also
need `replacevars` when their fixed auxiliary columns already exist. Reuse an
existing automatic identifier explicitly with `dmu(_dmu)`.

For radial `dea, bootstrap`, `reps1()` now drives actual pseudo-frontier DEA
replications and both bias and confidence intervals. It requires at least 50
replications and creates `_bc`, `_bias`, `_se`, `_lo`, and `_hi` suffixes. These
columns are included in `saving()`. The shared smoother now applies variance
correction before boundary reflection; earlier bootstrap inference results
should be recomputed. Adding `envars()` selects sequential SW2007 Algorithm 2:
`reps1()` controls score correction and `reps2()` controls regression inference.
This path creates only the extra `_bc` score column; regression results are in `e()`.
The dedicated `dea_sw2007` defaults to `procedure(sequential)` and retains the
previous shared-loop estimator as the explicit secondary `procedure(joint)`.
Its `backend(simarwilson)` calls the published implementation. See its help for
score conventions, boundary differences, and the Mac plugin architecture.
From `dea_sw2007` 1.2.1, `saving()` writes selected observation results for
both internal procedures and the external backend. Internal files contain
`prefix_eff`/`prefix_eff_bc`; external files contain the published inverse-score
`prefix_delta`/`prefix_delta_bc` (Algorithm 1 has no corrected-score column).
The file records its backend and score scale. `replacevars` does not authorize
file replacement; `replace` permits it and retains the legacy variable-replacement
alias. Estimation results remain in `e()`; use `estimates save` to retain them.
The dedicated `dea_bootstrap` shares the numerical routines. See the help for
score scales, engine behavior, and inference limitations.

Bootstrap inference needs assumptions beyond a successfully solved DEA model.
The radial smoother uses a common efficiency distribution and independent draws;
SW2007 additionally uses a separable common technology and a constant-variance,
truncated-normal regression model. The panel extension still resamples
observations independently. Small multidimensional samples can have poor
interval coverage, and bias correction can increase RMSE for some DMUs.
Increasing bootstrap counts does not resolve these issues. See the
[inference and convention checks](tests/README.md#sequential-sw2007-and-preserved-joint-procedure)
and each command's help before interpreting a nominal confidence level.

In `dea_sw2007` 1.2.2, a failed short regression-bootstrap fit gets one retry
on the same pseudo-data with 200 iterations before that sample is discarded.
Already valid fits are retained. `e(boot_retried)`, `e(boot_recovered)`, and
`e(boot_fit_attempts)` expose this work; `e(boot_failed)` counts rejected samples.
The original/corrected scores and point regression are unchanged. Recompute
earlier sequential covariance matrices and intervals if regression failures
were reported. The joint procedure and the external backend keep their own rules.

In `dea_crossgame`, aggressive, benevolent and all-strategy calls require
at least two observations after `if`/`in` and missing-value exclusions. A
singleton returns `r(2001)` before result replacement; an empty sample retains
`r(2000)`. Explicit `strategy(nash)` remains available for a single DMU as a
self-appraisal case. The numerical routines for accepted samples are unchanged.

Fuzzy DEA version 1.0.4 corrects the previously reversed `nirs` and `ndrs`
aliases. They now select `drs` (sum of reference weights <= 1) and `irs`
(sum >= 1), respectively. Recompute earlier fuzzy results that used either
alias. Version 1.0.5 additionally applies `tolerance()` (default `1e-10`),
scales LP columns, and checks feasibility and interval ordering. It retries
numerical failures before publishing results and protects existing result
columns if validation still fails. These numerical changes can alter the last
digits for every RTS specification. See the help for diagnostics and limits.
Version 1.0.7 supports close alpha levels such as `alpha(.001 .004)` by
adding `_1`, `_2`, etc. within each colliding rounded suffix group. Existing
names are retained when the rounded suffix is unique. The solver uses each
parsed alpha without rounding it to a percentage. Use `r(alpha)`,
`r(pessvars)`, and `r(optvars)` together to identify the results; the lists
are sorted by alpha. Repeated parsed levels are rejected before changing data.

In `dea` 2.0.41 with simplex 2.0.8, `stage(1)` skips slack maximization and
returns the radial LP's reference weights. Earlier simplex versions also ran
the second LP, so stage-1 references may change even when scores agree.
`stage(2)` remains the default. A dual-only recovered score has missing
reference weights until a primal solution is recovered. See the current
status contract below and the help for fallback behavior and slack interpretation.

Simplex 2.0.10 checks both basis-inverse identities, rejects singular
pseudo-inverses, and propagates failed refactorizations. It removes dependent
Phase-I rows instead of inserting arbitrary basis columns. Bland selection
uses minimum-ratio ties; the Harris step is bounded by basic feasibility.
After a numerical failure, one fresh-basis retry uses row, column and objective
equilibration and Bland pricing. This adopts the recovery/scaling approach
observed in Benchmarking/lpSolveAPI with an independent Mata implementation;
R is not a runtime dependency. Recovered variables are mapped back and checked
against the original constraints, and the objective is recomputed there.
These changes can alter scores or reference weights previously affected by an
invalid basis. They do not make floating-point acceptance an exact certificate
or change the fixed radial target. See [stability tests](tests/README.md).

In `dea` 2.0.51, `r(slack_in)` and `r(slack_out)` use the expanded input and
output variable names in the specified order. The common return path applies
to simplex, IPM and auto, including after VRS analysis, bootstrap and saving.
This changes matrix column names only; stage-1 calls still return no public
slack matrices. Use variable names instead of the former temporary names.

In `dea` 2.0.46 with IPM 2.0.8, TE stage 2 uses the correct CRS/VRS/NIRS/NDRS
technology and returns its projection weights and original-unit slack matrices.
Earlier IPM versions restricted the lambda sum incorrectly outside VRS, returned
stage-1 weights, and omitted the slack returns. Recompute affected references.
In `dea` 2.0.52, both TE engines report `r(stage_status)` and
`r(solution_source)`: radial/stage-2 outcomes and score/reference provenance.
A failed requested stage 2 retains the radial score but returns missing
references and slacks; simplex no longer substitutes stage-1 projections.
`r(n_failed)` counts each incomplete requested result once, and separate
radial/stage-2 counts identify the cause. `saving()` includes four diagnostic
columns without adding them to the active dataset. The existing `auto` engine
selection policy is unchanged; stage-2 failure alone does not switch engines.

`dea` 2.0.49 / IPM 2.0.11 adds explicit `radialtolerance(#)` recovery for
`engine(ipm) model(te) stage(2)`. The default zero preserves previous behavior.
A positive finite allowance below 1 permits one conservative radial target
after fixed-target stage 2 fails; it does not change original efficiency scores
or already accepted projections. `r(projection)` and `r(projection_trace)`
separate original results from the changed LP. `_dea_radial`, `_dea_radial2`
and `_dea_radial_adj` retain the original factor, projection factor and change,
including in `saving()` files. Output factors use phi, while scores use 1/phi.
This initial implementation requires explicit IPM and does not apply to
super-efficiency, simplex/auto or other models. Exact feasibility is not
certified; an accepted approximate projection is not retrospectively rejected.
See [dea help](ado/dea.sthlp) for the endpoint formula, counters and limitations.

`dea_malmquist` 1.3.5 and `dea_window` 1.3.6 reject repeated DMU-period
keys in the selected complete-case sample before replacing results (`r(459)`).
Duplicates outside the sample and unbalanced panels remain supported.
Malmquist 1.3.6 and window 1.3.7 number only selected observations when `dmu()`
is omitted, in the incoming within-period row order. Excluded rows no longer
shift IDs. Positional matching still assumes the same selected position denotes
the same real DMU in every period; supply an explicit, stable `dmu()` when order
or membership varies. `r(dmu_method)` reports `variable` or `position`.
Recompute earlier results that omitted `dmu()`.

Core `dea` 2.0.42 excludes missing `bad()` values in both undesirable-output
models and missing explicit `gx()`/`gy()` values in `model(ddf)` from evaluation
and reference samples. Numeric IDEA endpoints must be finite and ordered in
both interval-only and mixed calls; invalid selected endpoints return `r(459)`
before results or files are replaced. With both endpoints supplied, an unused
nominal value may be missing. Recompute earlier DDF, SBM-undesirable or interval
results that used these incomplete or invalid inputs.

Core `dea` 2.0.43 preserves user scalars and matrices during `reference`
output. Data columns named `_dmu_k`, `_dmu_j` or `_lam_val` no longer distort
the displayed reference lists. Earlier calls could overwrite same-named user
objects; the fix preserves computed scores, returned weights and saved data.

Core `dea` 2.0.44 requires the effective `tolerance()` to be positive and
finite across all core models and both public syntaxes. Invalid values return
`r(198)` before sample evaluation, result replacement or file writing. The
internal LP/TE contract commands also reject missing tolerances. The default
`1e-7` default is retained. Version 2.0.55 revises the legacy Form A options;
see the help for migration from `tol1()`/`tol2()` and `r(tolerance)`.

The Examples section in each help file generates its own teaching data, or uses
Stata's bundled `auto` data. Copy its command lines, omitting leading dots, into
the Do-file Editor and run the complete block from `preserve` through `restore`.
No separate dataset or example script is required. Small simulation samples and
replication counts demonstrate syntax, not an empirical research design.

```stata
preserve
clear
input id x y
1 1 1
2 2 3
3 4 4
4 5 3
end
dea x = y, dmu(id) rts(vrs) orientation(input) generate(eff)
list id eff
restore
```

## Choose a command

Use the help link to inspect syntax and examples on GitHub, or type
`help command_name` in Stata for the formatted version. The
[game-model overview](ado/dea_game.sthlp) compares the game-related commands.

| Family | Command | Purpose | Files |
| --- | --- | --- | --- |
| Compatibility | `dea_lp` | DEA with historical VRS/IPM defaults | [code](ado/dea_lp.ado) · [help](ado/dea_lp.sthlp) |
| Independent optimization | `lp` | Dataset-table LP and MILP | [code](ado/lp.ado) · [help](ado/lp.sthlp) |
| Independent optimization | `milp_bb` | Matrix-based MILP branch and bound | [code](ado/milp_bb.ado) · [help](ado/milp_bb.sthlp) |
| Core DEA | `dea` | Radial, SBM, economic, directional and super-efficiency models | [code](ado/dea.ado) · [help](ado/dea.sthlp) |
| Core DEA | `dea_congestion` | Cooper decomposition / legacy equality ratio | [code](ado/dea_congestion.ado) · [help](ado/dea_congestion.sthlp) |
| Core DEA | `dea_sbmundesirable` | SBM with undesirable outputs | [code](ado/dea_sbmundesirable.ado) · [help](ado/dea_sbmundesirable.sthlp) |
| Time and productivity | `dea_malmquist` | Malmquist productivity | [code](ado/dea_malmquist.ado) · [help](ado/dea_malmquist.sthlp) |
| Time and productivity | `dea_window` | Window analysis | [code](ado/dea_window.ado) · [help](ado/dea_window.sthlp) |
| Cross-efficiency and networks | `dea_crosseff` | Cross-efficiency | [code](ado/dea_crosseff.ado) · [help](ado/dea_crosseff.sthlp) |
| Cross-efficiency and networks | `dea_network` | Two-stage network DEA | [code](ado/dea_network.ado) · [help](ado/dea_network.sthlp) |
| Inference and uncertainty | `dea_bootstrap` | Bootstrap bias correction | [code](ado/dea_bootstrap.ado) · [help](ado/dea_bootstrap.sthlp) |
| Inference and uncertainty | `dea_sw2007` | Two-stage bootstrap inference | [code](ado/dea_sw2007.ado) · [help](ado/dea_sw2007.sthlp) |
| Inference and uncertainty | `dea_fuzzy` | Fuzzy DEA | [code](ado/dea_fuzzy.ado) · [help](ado/dea_fuzzy.sthlp) |
| Exploratory game diagnostics | `dea_crossgame` | Cross-efficiency strategy comparison | [code](ado/dea_crossgame.ado) · [help](ado/dea_crossgame.sthlp) |
| Exploratory game diagnostics | `dea_nash` | Historical capped allocation heuristic | [code](ado/dea_nash.ado) · [help](ado/dea_nash.sthlp) |
| Exploratory game diagnostics | `dea_cooperative` | Exact additive-game Shapley values | [code](ado/dea_cooperative.ado) · [help](ado/dea_cooperative.sthlp) |
| Cross-efficiency and networks | `dea_stackelberg` | Independent stages or CRS leader-follower DEA | [code](ado/dea_stackelberg.ado) · [help](ado/dea_stackelberg.sthlp) |

## Independent optimization and legacy DEA entry point

`dea_lp` now loads in a fresh session and retains the former VRS/IPM defaults.
Explicit `rts()`, `engine()` and `method()` choices take precedence; the rest of
its behavior follows current `dea`, including `replacevars`. Ordinary `dea`
still defaults to CRS/simplex. This restores an entry point, not every old option
or every historical numerical result. See [dea_lp help](ado/dea_lp.sthlp).

`lp` accepts the historical dataset-table layout: the first selected row is the
objective, and later rows are constraints. `intvars()` selects mixed-integer
optimization. It uses the shared continuous-LP adapter and `milp_bb`, with
simplex, IPM or automatic fallback. `milp_bb` is also available directly for
matrix input. Keep its loader, implementation file and the other `ado/` files
together. The older compiled `ldea.mlib` is not required.

These standalone optimization tools were previously preserved locally outside
the journal package allowlist. They are now included as maintained extensions;
this does not extend the article's validation claims. Read [lp help](ado/lp.sthlp)
for tolerance migration, stored-result differences and the solution-file format.
Always inspect `r(status)`; a normal command return can report infeasibility,
unboundedness, a solver failure or a node limit instead of an optimal solution.
The shared LP adapter now makes at most three native retries after an initial
failure or rejected candidate: equivalent equality inequalities, row/objective
scaling, and a tighter convergence setting. Each candidate must satisfy the
original constraints at the requested tolerance. Successful initial solutions
are returned immediately. `scaling(none)` controls the initial representation;
retries may still scale it. Explicit `engine(ipm)` never switches to simplex.

Continuous LP returns `r(attempts)`, `r(initial_status)`, `r(native_retries)`,
`r(status_conflict)` and `r(certificate)`. Exact sign/zero checks resolve some
infeasible or unbounded cases, including the independently reproduced Stata
16.1 false-infeasible report for max x, x>=0, -x<=0. Unresolved contradictions
return status 3 with no LP solution; a conflicted MILP node stops with status 5.
MILP and the internal radial translator count native retries and conflicts.
See the help for the narrow certificate scope and failure semantics.

This policy covers the shared LP adapter and its callers. The main `dea` TE
and auxiliary VRS paths now use native retries with separate per-DMU histories.
General optimality/status certification and other IPM numerical issues remain open;
the default table-LP engine is simplex.

## Versions and feedback

Core `dea` 2.0.40 connects numeric bound lists with `ordinal()`, `xord()`,
`yord()` or `bounded()` through the legacy DMU-wise transformation. Numeric
and string encodings of the same bounds return the same single transformed
score; `_opt` and `_mid` remain copies. `r(idea_method)` is `legacy`, and
`r(idea_numeric_vars)` identifies dimensions supplied by numeric bound lists.
The selected engine and automatic retry apply to this path.

Bound lists still follow the complete input/output order. At ordinal or
`bounded()` positions, repeat that variable's own name as a placeholder;
other endpoint variables in those positions are rejected to avoid ambiguous
specifications. String `bounded()` placeholders may coexist with numeric
bounds on the same side. An omitted endpoint uses the nominal model variable.
Selected numeric endpoints must be finite and ordered; validation precedes
result replacement and saving. Nominal values supplied by both endpoints
may be missing, while missing required endpoints cause r(459). Ordinary and
ordinal missing values are excluded before this validation. Numeric bound
lists alone use the separate interval calculation described above.
Earlier integrated versions could silently ignore ordinal options in numeric
mixed calls; versions 2.0.37–2.0.39 rejected that combination. Recompute earlier
mixed results using the now-explicit legacy score definition.
See the [mixed bounds tests](tests/README.md#numeric-bounds-and-legacy-mixed-data).

Core `dea` 2.0.39 fixes observation-order dependence in legacy ordinal output
transformations. Each evaluated DMU's original comparison value is held in a
double-precision scalar before modifying the data. **Recompute results from
earlier versions that used `yord()` or `ordinal()` on output variables**,
including calls combining them with `bounded()`: reordering the same data
could previously change scores. Input/output comparison rules and ties retain
their existing definitions. See the
[ordinal transformation and permutation tests](tests/README.md#ordinal-transformation-and-observation-order).

Core `dea` 2.0.38 restores string intervals in `bounded()` for both inputs
and outputs, using either `inputs = outputs` or `outputs, inputs(...)`.
Only variables explicitly listed in `bounded()` may contain string intervals;
other model variables must remain numeric. Selected bounds require finite,
ordered endpoints in `[lower,upper]` form (decimal/scientific notation and
surrounding whitespace are accepted). Invalid selected bounds stop with
r(459) before automatic IDs, result replacement or saving. Rows excluded by
`if`/`in` or missing ordinary model variables do not undergo this validation.
See the [bounded syntax tests](tests/README.md#bounded-string-syntax-and-validation).

Core `dea` 2.0.37 activates IDEA for `xord()` and `yord()` alone, with the same
legacy transformation used by `ordinal()`. **Recompute earlier ordinal-only
calls that omitted `model(imprecise)` and other IDEA activation options**:
they could silently produce ordinary TE scores. Ordinal names repeated within
or across options are transformed once. `xord()` must name input variables,
`yord()` output variables, and every bounded/ordinal variable must have one
model role. Distinct bounded and ordinal variables can still be combined.

Conflicting model selections, invalid variable roles and ambiguous
bounded/ordinal assignments are rejected before result replacement.
See the [IDEA option tests](tests/README.md#idea-activation-and-option-combinations).

Core `dea` 2.0.36 supplies FDH/IDEA summary statistics and failure counts,
reports the effective dispatched model in `r(model)`, and retains the normalized
model option in `r(model_requested)`. FDH spellings now share the VRS label and
reject unsupported NIRS/NDRS requests. IDEA no longer triggers the ordinary
radial VRS decomposition. Existing raw score names and solver calculations
are retained. Interval IDEA generates the primary pessimistic score and
`_opt`; legacy bounded/ordinal IDEA also generates `_mid`, but its `_opt` and
`_mid` are copies of the transformed primary score. `r(idea_method)` identifies
the path and `r(idea_vars)` lists its columns, all now included in `saving()`.
Missing primary or `_opt` results count as failed rows once; separate counts
and `_opt` summaries are available. See the
[FDH/IDEA contract tests](tests/README.md#fdh-and-idea-results).
Do not interpret the legacy copies as independently estimated uncertainty
bounds. Ordinal activation and option validation are corrected in 2.0.37,
as described above.
Version 2.0.56 corrects the separate numeric interval formulation's own-row
inconsistency, which previously caused optimistic scores above 1 or infeasible
optimistic LPs and could alter pessimistic scores. Failed solves remain missing;
the correction does not replace failures with an efficiency of 1.

`dea` 2.0.29 fixes persistent `saveref` variables, refreshes returned results
after automatic IPM retries, and saves only the final results once. Automatic
DMU identifiers now preserve an existing `_dmu` by requiring an explicit
`dmu()` choice. See the [reproduction tests](tests/README.md).
The core `dea` parser accepts the help's documented minimum abbreviations and
support legacy Form A `using filename` as an alternative to `saving()`.

The eleven extensions with Form A input/output syntax now share
`_dea_parse_io.ado`. Their `if` expressions can contain equality comparisons,
quoted strings and functions with commas in both Form A and Form B. Network,
Malmquist and window commands also retain positional output varlists; omitting
that varlist requires `outputs()`. See the syntax regression tests for sample
and historical-result comparisons.
Extensions now accept the minimum option abbreviations advertised in their
help, including the corrected `int()`, `r()`, and `or()` entries. The option
tests cover all allowed prefix lengths and preserve existing full names and
aliases. Specify only one of `ort()` and `orientation()`; fuzzy version 1.0.6
also rejects this duplicate choice before creating or replacing results.

Core `dea` 2.0.33 corrects the common display's DDF counts and ranks and the
output-oriented super-efficiency counts. DDF/additive scores rank from low to
high; both super-efficiency orientations rank from high to low (output radial
scores are stored as `1/phi`). `r(N_eff)` includes super-model ties at 1;
`r(N_super)` separately counts scores strictly above `1+tolerance`.
Frontier ties affect only the displayed ranking key, never the stored scores.
See [score classification tests](tests/README.md#score-counts-and-ranks).

Core `dea` 2.0.34 and engine modules 2.0.7 correct the good-output beta
coefficient in `model(undesirable)`. **Recompute earlier `yb`/`xyb` results**:
identical positive CRS DMUs could previously receive beta=1 instead of 0.
`yb` expands good outputs and contracts bad outputs; `xyb` also contracts
inputs. Column scaling, tighter internal precision and primal-solution checks
improve numerical stability. Missing `bad()` values now leave the estimation
sample; negative selected data, zero observed directions and invalid tolerances
are rejected before result replacement. See the
[independent undesirable-DDF tests](tests/README.md#undesirable-output-ddf).
This correction is separate from `model(sbmundesirable)` and does not change
the specified RTS restrictions or settle alternative technology formulations.

Core `dea` 2.0.35 preserves one VRS-decomposition row per selected observation,
including repeated or missing DMU labels. It no longer averages by label or
leaves unused zero rows. `r(vrsfrontier)` follows `r(efficiency)` order and uses
row names `r1` through `rN`; the new `r(vrsobs)` maps these rows to source
observation numbers at the time of estimation. Scripts matching old VRS row
names to DMU labels should use the mapping instead. See the
[VRS observation tests](tests/README.md#vrs-observation-mapping).
Core `dea` 2.0.50 displays string DMU labels literally in score and reference
tables, including quotation marks, macro punctuation and formatting directives.
FDH and numeric interval IDEA tables use the same display rule. Earlier versions
could stop with `r(198)` after creating scores when a label contained a quote.
Original identifiers and numerical results are preserved.
Saved reference names also replace backticks; the full peer identifier remains
available in each reference variable's `dea_peer_id` characteristic.

This is a development working copy. Individual command versions are recorded in
their `.ado` headers; help revisions have a separate dated marker. A help-only
revision does not change the computational command version. A future tagged
release will identify a tested combination of files.

For a problem report, include the command name, Git commit or release tag,
`which command_name` output, Stata version, a small reproducible example, and
expected versus actual results. Use [GitHub Issues](https://github.com/renegade68/DEA-Stata/issues).
Create synthetic example data instead of attaching private research data.

## Attribution and licensing

Preserve the author and reference information in the code and help files.
Cite the article above for the software and each command's references for the
method. The journal citation and article link will be updated when a DOI and
complete publication details become available.

The software is licensed under the [Academic Free License 3.0](LICENSE),
continuing the DEAS license used on SourceForge. See [the software notice](LICENSE_CODE.txt)
for scope and attribution. Article manuscripts and publisher PDFs are not included
and are not covered by this software license.

`dea_network` and `dea_stackelberg` validate nonnegative quantities with a positive
input, intermediate, and output in each selected observation, pass positive
finite `tolerance()` to checked/scaled LPs, and expose component failure counts.
`matrix detail` returns weights or reference intensities and solve diagnostics.
Network's default CRS/input multiplicative model shares intermediate weights;
its secondary decomposition maximizes Stage 1 at the system optimum. Other
multiplicative combinations explicitly report the historical independent-stage
product. Additive is an equal average of independent stages.

`dea_stackelberg, method(noncooperative) rts(crs)` holds the selected leader's
independent optimal score fixed and optimizes the follower using shared
intermediate weights. Both leader choices are supported. A single active
intermediate uses the equivalent independent-stage reduction, with shared
weights checked against the full constraints; `r(follower_route)` records this. The default
`method(independent)` preserves historical stage ordering and supports CRS/VRS.
A failed joint solve retains the leader/benchmarks and leaves follower/system
scores missing. Check `r(formulation)`, failure counts and `r(follower_status)`;
VRS noncooperative optimization is not yet supported.


`dea_cooperative` and `dea_nash` (1.3.0) share input-oriented CRS/VRS DEA
scores with bounded exact certification against the original binary64 inputs.
Successful original-data certificates replace the numerical score; native LP
status is retained separately. Certification currently covers n<=32 and
inputs+outputs+2*(VRS)<=16. Outside this scope, an available native estimate is
labelled unverified. Read `r(score_diagnostics)` and the quality/reason fields;
a returned number alone is not evidence of certification.

`dea_cooperative` retains the fixed-frontier characteristic function
`v(S)=sum(e_i)`: Shapley values equal individual efficiencies and contributions
are zero. Exact score fractions determine ranks and totals, even when different
scores round to the same exported double. `detail` skips enumeration for n>15.

`dea_nash` retains the historical scenario transform and capped proportional
allocation. It is not a Nash-product optimizer or a centralized DEA model.
Certified fractions feed the scenario, gain and allocation formulas directly.
In particular, gain is exported from its exact expression: it can be positive
even when exported bargain and disagreement scores round to the same double.
Per-context arithmetic limits and export precision checks may withhold a derived
field without discarding independently certified scores or other valid fields.
A missing input score prevents complete allocation. Read `r(result_quality)`,
`r(result_reason)` and `r(aggregate_diagnostics)` with the result matrices.
The logarithmic product is a numerical diagnostic, not a rational certificate.

Install the complete matched folder, including `_dea_games_load.ado`,
`_dea_games_exact.ado` and `_dea_games_quality.ado`. See
[game installation and verification](INSTALL_GAMES.txt),
[portable test contracts](tests/GAMES_QUALITY.txt), and both help files for
formulas, arithmetic/export limits and precise quality codes. A substantive
pooling game or formal feasible-set Nash model requires a separate method choice.


`dea_fuzzy` 1.2.0 preserves the triangular Kao-Liu best/worst alpha-cut model.
The joint cuts are Cartesian products of marginal intervals; the finite grid
samples a membership function and does not supply probabilities, defuzzification
or alpha-weighted rankings. At the smallest requested cut, each row must have
at least one positive lower-cut input and output. Other zero entries remain
allowed. Domain failures are rejected before result replacement.

The summary now uses the same complete-grid DMU cohort at every alpha and
returns per-cut availability counts in `r(summary)`. These are descriptive
marginal-bound averages; different DMUs' extrema need not be simultaneously
attainable. `matrix` returns both endpoint grids, complete-row flags and original
observation indices. Nesting checks bridge missing cuts, alpha=1 uses the same
center matrix, and quoted/numeric identifiers are displayed safely. See
[dea_fuzzy help](ado/dea_fuzzy.sthlp) for the formulas and missing-score contract.

Requested alpha=1 centers also receive bounded original-data simplex and exact
primal/dual checks, for both orientations and all four RTS choices. A verified
score replaces both center endpoints while native scores and failure reports
remain in `r(center_diagnostics)`. Counts distinguish certified, unverified and
unavailable DMUs; unverified or unavailable centers produce warnings. The
current scope is at most 32 selected DMUs and 16 multiplier variables including
RTS extras, with additional arithmetic/work/export limits. Certification uses
the stored binary data and does not cover alpha<1 cuts or summary means.
See [fuzzy installation and verification](INSTALL_FUZZY.txt) and
[portable center regressions](tests/FUZZY_QUALITY.txt).


Core `dea` 2.0.59 makes three environmental-model conventions explicit.
SBM-undesirable uses separable reduction slacks and strictly positive denominator
data; NIRS/NDRS now enforce the transformed RTS restrictions, so earlier results
with those options must be recomputed. Its objective is non-oriented.
The observed-direction bad-equality DDF is retained; its VRS/NDRS technology is
not a general joint weak-disposal model with activity-specific abatement.

`dea_congestion` 2.1.0 retains the historical ratio of ordinary VRS efficiency to
an input-and-output-equality radial score. The ratio is neither sufficient nor
necessary to identify congestion, and it is not the Cooper slack decomposition.
The help gives false-positive and false-negative examples. Column scaling and
nonnegative/positive-row validation improve its numerical contract without
changing that historical model. Specialized backend/formulation returns identify
what was calculated. See the three command help files before interpreting results.


### Cooper congestion decomposition (explicit option)

`dea` 2.0.73 / `dea_congestion` 2.3.0 retain `congmethod(legacy)` by default.
Explicit `congmethod(cooper)` decomposes strictly positive-data VRS input slacks
into technical inefficiency and congestion. The primary result is mean relative
congestion, with zero meaning none; it is not the legacy equality ratio.

```stata
dea_congestion x1 x2 = y1 y2, dmu(id) congmethod(cooper) congweight(normalized) engine(ipm) congmatrix generate(cp)
tempname A M
matrix `A'=r(lambda_stage1)
matrix `M'=r(lambda_stage2)
dea y1 y2, inputs(x1 x2) model(congestion) congmethod(cooper) congweight(normalized) engine(simplex) conganchor(`A') congpolicy(stage2) congreference(`M') congranges dmu(id) generate(cpn)
```

Simplex now means package revised-simplex proposals plus bounded integer
primal/dual certificates, with exact optimal-face linkage. Its domain requires
`2*N + inputs + outputs + 3 <= 64`; the optional stage2 selector defaults to16 exact basis exchanges and can
be explicitly extended to128. Main-stage repair retains16. This is separate from numerical simplex
iterations. `congpolicy(stage1)` keeps the baseline decomposition; optional
stage2 selects a nearest reference on the original technical-optimal face.
References are positional N-by-N matrices matching the selected sample/order.
Missing reference rows or uncertified selection return the certified baseline
with explicit status. No automatic IPM prepass/fallback or budget escalation is
used. Cooper `engine(auto)` is rejected. Existing defaults differ: `dea` defaults
simplex and the wrapper defaults ipm, so specify the engine for comparisons.

IPM retains native `LinearProgram()` and numerical retries. Both engines reject
column/objective scaling that underflows to zero or overflows. Their success
criteria differ: simplex certifies the internal rational point of the normalized
binary64 LP; native availability is numerical. `r(actual_engine)` and
`r(actual_backend)` report the actual path. Returned simplex components and
ranges are exact fractions rescaled to original units then rounded once;
returned binary64 coordinates themselves are not exact certificates.

`congmatrix` returns components, targets, peers and statuses. Simplex adds
`r(certificate_status)` and explicit selection/conversion counts. Its six-variable
point bundle is atomic, and a range export failure does not erase a usable point.
`congranges` uses the original slack and conditional technical faces, not a
narrowed reference-selection face. These are not global congestion ranges over
all stage1 projections. Components need not attain their extrema jointly.
For simplex, `tolerance()` classifies summaries rather than relaxing exact LP
acceptance. See the [full contract](ado/dea_congestion.sthlp) and
[matched-package installation guide](INSTALL_COOPER.txt).


`dea` 2.0.61 adds explicit `model(undesirable) disposability(weak)` for a
convex joint weak-disposal DDF. The VRS construction separates active and idle
intensities following Kuosmanen (2005); NIRS/NDRS restrict their total intensity.
The default `disposability(equality)` keeps the existing technology and results.

```stata
dea x1 x2 = y1 y2, model(undesirable) bad(b1 b2) rts(vrs) disposability(weak) weakmatrix dmu(id) generate(wd)
```

The score is directional distance beta (zero is efficient). `weakmatrix`
returns active/idle/total intensities, directional targets, remaining slacks and
primal diagnostics. This is one LP, with no secondary slack optimization; both
stage labels give the same result. The new option honors `engine(simplex)` and
`engine(ipm)`; `engine(auto)` tries native only after a simplex failure.
Per-DMU diagnostics distinguish the actual solver and equivalent retries.
Unavailable scores stay missing with an explicit `_status` field. Nonnegative
data require a positive input and good output per selected DMU; zero-bad data
are permitted but do not imply null-jointness. The model imposes no materials
balance. See help for direction, RTS, output protection and numerical limits.

### Two-stage numerical quality

Network and Stackelberg commands retain native diagnostics while verifying
small models against the original stored inputs. Independent products and means
are composed exactly before export. Unrepresentable positive aggregate values
are withheld with a warning; certified stage scores remain available.
Use `matrix` to inspect `r(score_quality)`, `r(score_reason)` and native provenance.
See [two-stage installation and verification](INSTALL_TWOSTAGE.txt).
