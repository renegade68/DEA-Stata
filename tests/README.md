# Running the public regression tests

Use Python 3.10+ and a licensed Stata 16+ batch executable. Run the following
from the package root; no source-code edits or personal research folders are needed.
Current-package suites work in both a Git clone and an extracted source ZIP.
The batch launch protocol is intended for macOS/Linux; this revision was
executed on macOS. Windows execution has not been validated.

```sh
python3 -m pip install -r tests/requirements.txt
python3 -B tests/run.py --list
python3 -B tests/run.py --suite smoke
python3 -B tests/run.py --suite crosseff_raw
python3 -B tests/run.py --suite crosseff_decoder
python3 -B tests/run.py --suite crosseff_dyadic
python3 -B tests/check_publication.py
```

Install dependencies into your own virtual environment. The runner never installs
software, downloads packages, or modifies Stata settings. It finds Stata on PATH
or in a standard macOS application location. If necessary, set DEA_STATA once to
your installed executable, or pass --stata. DEA_RSCRIPT/--rscript is the equivalent
R setting; R comparisons run only when --with-r or --rscript is requested.
R packages required by a suite are documented below.

The default output is a unique OS temporary directory outside the package.
Use --output with a new directory outside the package for durable local evidence.
Suites run sequentially. A suite failure stops the runner and returns a nonzero
exit code. Smoke selects three help-example cases and the full fixed raw-recovery
suite; it is not the full regression collection. --list identifies suites needing
a separate historical --baseline-ado. Missing historical inputs are errors, not
passed or silently skipped tests. Do not use the current ado directory as a
historical baseline. Historical Git comparisons require the referenced commit
history; they cannot be reconstructed from a source ZIP alone.

Only shareable-summary.json from the common runner is designed for sharing:
it contains suite names and pass/fail status, not local paths, commands or logs.
Raw console files, Stata logs, generated do-files and per-suite manifests can
contain local paths and should remain local unless reviewed and redacted.
See PUBLICATION.txt for the publication audit and commit-metadata checks.

The decoder suite verifies exact finite binary64 decomposition, the unchanged
default coefficient domain, and invalid-input/arithmetic-budget boundaries.
The raw suite preserves fixed legacy selections and independently certified
references for the additional recoveries in version1.10.1.

The dyadic suite generates signed synthetic integers to check exact products,
differences, integer quotients, restored values and residual signs. It verifies
zero and negative denominators, noninteger quotients, width/budget failures and
preservation of inputs and ordinary-path integer work counts. Its fixed seed
does not require captured research files. The raw references additionally check
finite secondary certification versus exact unboundedness in version1.10.2;
the selected weights and efficiency values are unchanged.

The detailed suite-specific commands below remain available for advanced use.
They assume DEA_STATA (and, if needed, DEA_RSCRIPT) is already set and a new
output directory is chosen. Prefer the common runner for automatic discovery.
Optional historical snapshot paths are examples to replace with a reviewed
baseline, not required private files. Requirements and limitations differ by suite.

## Detailed suite reference

## Economic efficiency and price specifications

```sh
python3 tests/test_economic_models.py \
  --stata "$DEA_STATA" \
  --rscript "$DEA_RSCRIPT" \
  --output ../dea-test-results/dea-economic-models \
  --baseline-ado /path/to/pre-R11A/ado
```

Requires NumPy/SciPy; optional R checks require Benchmarking. Independently
formulated lambda-only LPs check cost, revenue, restricted profit, and both
allocative decompositions. The grid covers unequal dimensions, individual prices,
four RTS restrictions and aliases, both syntaxes, three engine selectors, unit
and price changes, ordering, duplicates, zeros and incomplete/subset samples.
Small analytic fixtures check the score conventions. Price-list failures must
preserve caller data, existing scores, labels and saved files. The optional
baseline compares valid results exactly with the pre-fix source.

R cost/revenue and radial scores are compared after converting output expansion
to the stored reciprocal convention. The documented profit model restricts
improvements to observed input/output bounds, so R's unrestricted `profit.opt`
is a different optimization problem. This suite does not establish all R11
models or universal numerical reliability.

For `test_output_protection.py`, use `--baseline-ado /path/to/pre-change/ado`
when comparing against a reviewed working snapshot. The default Git baseline
predates the R08/R09 numerical and interval-formulation repairs; equality with
those superseded results is not a regression requirement. The chosen snapshot
and its hashes are recorded without relaxing the comparison assertions.

## Checked inverse identities and cutoff boundaries

```sh
python3 tests/test_inverse_validation.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-inverse-validation
```

Requires NumPy. The included hexadecimal fixture preserves exact input doubles
without depending on private research folders. Fraction arithmetic checks 120
matrix decisions, 350 input pairs, eight invalid-input controls, and 504
extreme-scale cutoff neighbors. Cases cover cancellation, both inverse identities,
fast-path boundaries, exact binary scaling, and conservative unavailable results.
See [INVERSE_VALIDATION.txt](INVERSE_VALIDATION.txt) for the error-bound assumptions,
implementation limits and sources. The full maximum supported matrix dimension
is not exercised, and these finite tests do not guarantee every numerical case.

The cached-engine checks in `test_sbm_models.py` and
`test_nonte_tolerance.py` verify the readiness
marker declared by the tested package, so they work with subsequent compatible
engine revisions as well as the original SBM repair.

## Failure-only original-basis precision recovery

```sh
python3 tests/test_precision_recovery.py --stata "$DEA_STATA" --output ../dea-test-results/dea-precision-lp
python3 tests/test_precision_arithmetic.py --stata "$DEA_STATA" --output ../dea-test-results/dea-precision-arithmetic
```

Requires NumPy/SciPy and Stata 16+. The LP fixture stores exact binary64 coefficients
as hexadecimal strings and checks known primal/dual optima independently.
The arithmetic suite checks compensated products, outward bounds, cutoff
neighbors, signed/cancelling sums, near-zero objectives and finite-range rejection.
These tests use only the packaged sources and fixtures; no private research
folders, downloads or R packages are required. A recorded preexisting accepted
objective-boundary case is reported explicitly, not counted as a new recovery.
See [PRECISION_RECOVERY.txt](PRECISION_RECOVERY.txt) for the contracts and limits.

## Wide transition recovery and exact initial row scaling

```sh
python3 tests/test_wide_scaling.py --stata "$DEA_STATA" --output ../dea-test-results/dea-wide-scaling
```

Requires NumPy/SciPy and Stata 16+. The 269 checks cover thin column selection,
binary multiplier boundaries, 21 exact-witness wide LPs and 124 public LP
commands. All data are packaged; rational witnesses check original units.
The separate precision fixture now requires 279 exact recoveries (114 original,
153 transition, one certified-price and11 reconstruction recoveries) and396
available results among396 conditions.
Known failures and the preexisting accepted objective boundary remain explicit.
See [WIDE_SCALING.txt](WIDE_SCALING.txt) for scope and supported limits.

## Certified prices and repeated pivot states

```sh
python3 tests/test_precision_pricing.py --stata "$DEA_STATA" --output ../dea-test-results/dea-precision-pricing
```

Requires NumPy/SciPy and Stata16+. The packaged hexadecimal fixture contains
1,588 pricing helper controls and18 positive LP conditions. Exact rational
arithmetic checks price enclosures and every certified entering/optimal choice.
The LP controls check new exact recoveries, old weights on duplicate columns,
skipping certified prices when unnecessary, and recoverylog OFF/ON parity.
Synthetic pricing helpers include negative costs; they are not positive-data
DEA examples. Instrumentation is confined to a private copy, whose numeric
outputs must match the tested package. See [PRECISION_PRICING.txt](PRECISION_PRICING.txt).

## Original-column reconstruction after failed pivots

```sh
python3 tests/test_reconstruction.py --stata "$DEA_STATA" --output ../dea-test-results/dea-reconstruction
```

The packaged exact-witness controls include96 reconstruction recoveries,15 later
exact-duplicate recoveries and3 earlier accepted results. The private
observational copy checks actual iteration sums and logging OFF/ON parity.
Fault controls and wide boundary checks reject false recovery. See
[RECONSTRUCTION.txt](RECONSTRUCTION.txt). These finite conditions overlap other
suites and are not independent empirical datasets.

## Simplex basis stability and equilibrated recovery

```sh
python3 tests/test_simplex_stability.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-simplex-stability
```

Requires NumPy. Tests singular and ill-conditioned basis inverses, failed
refactorization/pivots, exact Bland ratio ties, bounded Harris steps, tiny eta
updates, dependent/zero rows, infeasibility, unboundedness, Beale's cycling
example, objective units and tiny positive pivots. Another 160 LPs have
independently checked exact binary-rational primal/dual certificates, including
degeneracy, duplicate constraints and row/variable unit changes. The resulting
solutions are checked in original coordinates against certified objectives.
`--source /path/to/old/ado` compares the LP cases with an older package, omitting
new helper assertions. It is expected to fail if the old defects reproduce.
This suite is a finite regression set, not a proof for every possible LP.

The output-protection and sample-selection suites allow `1e-12` differences
from committed radial scores in the documented VRS/SW2007/bootstrap fixtures,
with separate analytic frontier checks. Current clean fits, protected
replacement, and qualified versus restricted samples still compare exactly.
These historical tolerances accommodate final-bit changes from refactorization.

## Reference weights and saved results

Run with Python 3, Git, and a licensed Stata 16+ batch executable. No additional
Python packages or network access are needed. Choose a new output directory:

```sh
python3 tests/test_saveref.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-saveref-validation \
  --baseline-ref 41ea1dccda72b4a49c0964e358a0f76196709994
```

The optional baseline argument reproduces the old defects and compares 32
solver/specification combinations against that committed version. Omit it to
test the current working files alone. Each case runs in a fresh Stata batch
session with isolated ado paths. The runner returns a nonzero status if any
case fails; it checks a completion marker because Stata can exit with status
zero even when a do-file fails.

Coverage includes:

- Simplex and IPM; CRS, VRS, NIRS and NDRS; input/output orientation; stages 1/2.
- Reference variables equal to the returned matrix, primal feasibility,
  convexity restrictions, and unchanged baseline scores. Stage-2 weights also
  match; stage-1 simplex weights may change after the F12 routing fix.
- Sample alignment, missing values, names, metadata, collisions and automatic IDs.
- New and existing files, write failures, internal-variable exclusion and
  preservation of unrelated variables in memory.
- Automatic retries, multiplier weights, radial super-efficiency, partially
  infeasible VRS super-efficiency, selected other saving paths, and 120-DMU samples.

Retry tests use a separate generated copy that invalidates the initial simplex
attempt before the display/return step. The real IPM solver remains unchanged.
These tests establish control-flow behavior, not the frequency of numerical
failures in empirical data. VRS super-efficiency additionally exercises natural
infeasibility; an infeasible DMU remaining missing is an expected model outcome.

Generated source copies, `manifest.json`, `results.json`, and each case's
`test.do`, `test.log`, and datasets remain in the specified output directory.
Instrumented copies must not be distributed as the package. Successful old-defect
tests mean the defect was reproduced, not that the old version was correct.

These checks do not establish correctness of every DEA model or performance on
very large samples. Cross-engine weights need not be identical when optima are
nonunique; baseline comparisons use the same engine and the same data.

## Reference display and caller objects

```sh
python3 tests/test_reference_state.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-reference-state \
  --baseline-ado /path/to/pre-F15/ado
```

Uses Python's standard library and Stata. The optional baseline is core 2.0.42.
The grid checks scalar/matrix inventories, values, string types and matrix
labels before and after `reference`, including absent objects, numeric/string
scalars, matrices and mixed types. It covers both syntaxes, all three engine
selectors, four RTS restrictions, both orientations and stages, if/in, numeric,
string, Unicode, repeated and missing identifiers, automatic IDs and identifiers
using the former scratch names. Hostile data columns with the same names must
not suppress or mislabel peers. Protected object names can also be score outputs.

A Python renderer compares reference text with every returned lambda row;
returned scores, scalars, saved weights and original data must match controls.
A fabricated matrix checks the six-peer display limit, strict tolerance cutoff,
negative/zero/missing weights, and rounded printed values independently of LP
solutions. Additional cases cover natural super-efficiency infeasibility,
private forced retry, output collisions, file errors and successful quoted-label
display. Caller objects remain intact on both successful and failing calls.
Historical controls compare scores, reference weights and reference text exactly.
The help-example harness also protects numeric/string scalars and a matrix
while running each literal Examples block. Production code is never instrumented.

## DMU label display

```sh
python3 tests/test_id_display.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-id-display \
  --baseline-ado /path/to/pre-F27/ado
```

The optional baseline enables exact historical comparisons for ordinary labels.
Current cases compare special-character labels with numeric-ID controls across
core TE, FDH, numeric interval IDEA and legacy ordinal calls. They check scores,
returned matrices/scalars, source data, saved datasets and caller objects.
Quotes, compound quotes, macro punctuation, literal SMCL, Unicode, long strL,
blank/duplicate labels, fractional/large numbers and extended missings are covered.
Both syntax forms, if/in, replacement, engines, orientations and stages are used.
The legacy path has aggregate-only output; it is checked for data/results rather
than a per-DMU table. Fabricated references check quoted failed-score rows and
the six-peer limit. A private FDH failure injection covers its FAILED row.
All 95 printable ASCII characters are checked in saved reference names and
verbatim peer-ID metadata. No computational engine is modified by these tests.

## Slack matrix column identities

```sh
python3 tests/test_slack_names.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-slack-names \
  --baseline-ado /path/to/pre-R02/ado
```

The optional baseline compares every returned numeric matrix and scalar as
binary doubles against the pre-R02 code. Current checks require original input
and output column names, their specified order, and observation-aligned rows.
Both syntax forms, all RTS families, both orientations, simplex/IPM/auto,
reordered/asymmetric/32-character varlists, wildcard/range expansion and
`dea_lp` forwarding are covered. Repeated calls allocate different temporary
variables; names and slack values must remain stable. Projection equations are
checked against the named original-unit data when references are available.
The wildcard/range fixture contains linearly dependent columns and explicitly
requests `tolerance(1e-9)` for that check. At the default tolerance, native IPM
gave the same normalized residual of about `1.039e-6` before and after R02;
this existing numerical behavior is recorded separately from the name fix.

Private copies exercise forced auto fallback and preserved extended missing
slacks. Additional cases cover saving, replacement, ordinary/environmental
bootstrap and absence of stale slack matrices after stage 1 or FDH. The option
abbreviation suite also compares the complete stage-2 `dea_lp` returns,
including matrix stripes, across full and abbreviated spellings.

## Simplex stages and reference solutions

```sh
python3 tests/test_simplex_stages.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-simplex-stages \
  --baseline-ado /path/to/pre-F12/ado
```

Requires NumPy and SciPy. The optional baseline is the pre-F12 package
(`dea` 2.0.40), including all earlier local fixes. Test-only counters measure
actual radial, dual-fallback and slack-maximization LP calls. Public stage-1
calls must never run the second LP, including auxiliary VRS fits, legacy
ordinal/mixed routes, and multiplier score paths. Existing optional Mata
signatures retain their stage-2 defaults.

Independent HiGHS LPs verify radial scores, reference feasibility and the
stage-2 maximum normalized slack sum. The grid includes all eight RTS
spellings, both orientations, simplex/auto/IPM, both public syntaxes, selected
samples, repeated observations and heterogeneous units. Stage-1 references
are allowed to differ from historical stage-2 solutions; successful stage-2
scores and weights are compared exactly against the snapshot. The deliberately
nonunique fixture demonstrates equal radial scores with distinct reference
solutions and a larger normalized slack sum after Stage 2.

Fault injection in private engine copies tests primal failure with dual score
recovery, both solvers failing, second-stage failure, automatic IPM retry,
missing unrecovered references/slacks, and equality of saved reference columns
to returned matrices. It establishes control flow, not empirical failure rates.
The production source is never instrumented. Actual help examples also check
equal stage-1/stage-2 scores and absence of public stage-1 slack matrices.

## Panel DMU-period keys

```sh
python3 tests/test_panel_keys.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-panel-keys \
  --baseline-ado /path/to/pre-F13/ado
```

Requires NumPy and SciPy. The optional baseline is the pre-F13 package
(Malmquist 1.3.4 and window 1.3.5). Tests reject repeated selected keys before
result replacement, including identical observations, different repeated values,
multiple groups, both public syntaxes, numeric/fractional/large IDs, Unicode,
quoted strings, whitespace and strL labels longer than 2,045 characters.
Rejection checks preserve complete data values, types, labels, variable order
and metadata. Duplicates excluded by if/in or missing model/key values are
compared with physically restricted samples. String time variables fail early.

Normal adjacent/global Malmquist and window results are compared exactly with
the snapshot and with separately constructed HiGHS LPs, including both
orientations, all supported RTS values, one/two dimensions, unbalanced panels,
and gaps in numeric time. Row permutations check score-to-observation mapping.
Historical duplicate acceptance and incomplete per-row summary storage are
reproduced separately. `_dea_panel_check.ado` must be distributed with the two
commands. The shared checker does not sort or modify the Stata dataset.

A successful duplicate-key check cannot establish real-world identity for
automatically numbered observations. Use an explicit stable panel identifier
when order or membership varies; the selected-position checks below cover
omitted `dmu()` calls.

## Selected positional panel IDs

```sh
python3 tests/test_panel_autoid.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-panel-autoid \
  --baseline-ado /path/to/pre-F31/ado
```

Requires Python's standard library and Stata. The optional baseline is the
pre-F31 package (Malmquist 1.3.5 and window 1.3.6). It reproduces the excluded-row
ID shift; the current package must match explicitly supplied selected-position
IDs exactly. Python independently computes expected IDs in original row order.
Test-only ado copies expose the actual selected sample and generated IDs for
comparison; production files are never instrumented.

The grid covers adjacent/global Malmquist and window analysis, both syntaxes,
all supported RTS values and orientations, one/two inputs and outputs, if/in,
missing values, unbalanced selection, four row orders, fractional time gaps,
and two Stata sort seeds. It checks returned scalars, every result column,
original data and row order, metadata, sorted-by state, and temporary cleanup.
Additional cases cover empty/short samples with existing results, replacement
when the qualifier refers to an old result, and numeric/string explicit IDs
under row reversal. Literal help examples compare positional and explicit IDs.
The unbalanced fixture verifies positional numbering, not real-world identity:
removing different DMUs can still invalidate the positional assumption.

## Help examples and public syntax

```sh
python3 tests/test_help_syntax.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-help-syntax-validation
```

This runs the marked example blocks directly from all packaged help files,
including their original replication counts and `preserve`/`restore` sequences.
It decodes SMCL brace/indentation escapes and records the help source line for
each command in `help_line_map.json`; it does not silently insert `dmu()` or
modify example commands to make them pass. The fuzzy help's additional alpha=1
check in Remarks is also executed.

The additional cases compare both syntax forms, `if`/`in` with numeric and
string equality, solver choices, reference/saving combinations, documented
abbreviations and aliases, all 16 core models, bootstrap, wildcard/range
varlists, legacy `using`, and expected rejection of conflicting options.
Use `--only REGEX` to rerun selected cases in a new output directory.

## Output-variable protection

```sh
python3 tests/test_output_protection.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-output-protection-validation
```

This uses commit `97cbb4c6c35be93a425934c05649f5682a99bfdd` as the default
pre-change numerical control. It covers all 14 public commands, both core
syntax forms and engines, all 16 core models' input-collision rejection,
primary and auxiliary output collisions, full and partial explicit replacement,
protected model/identifier variables, automatic IDs, exact-name versus
abbreviated-name collisions, and the distinction between file and variable
replacement. Rejection checks compare all data, variable order, types, formats,
labels, value-label assignments, and custom variable characteristics.

Normal outputs and explicitly replaced outputs are compared with the committed
control on the same fixture and seed. These checks validate data protection and
unchanged results; they do not validate the statistical theory of every model.
Known review issues outside output protection remain separate work items.

Sequential help examples and existing numerical regression cases now spell out
`replacevars` where an earlier call created output columns to be refreshed.
The runner still executes the literal help commands. Source snapshots include
new, non-ignored package files so an uncommitted helper is tested before staging.

## Positive finite numerical tolerance

The canonical LP loader/implementation is now 0.2.2; malformed matrix shapes
also return an invalid-input status before entering typed diagnostic functions.

```sh
python3 tests/test_tolerance_validation.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-tolerance \
  --baseline-ado /path/to/pre-F16/ado
```

Uses Python's standard library and Stata. The optional baseline is core 2.0.43.
Invalid-value cases cover both public syntaxes, all engine requests, 25 model
profiles (including FDH, numeric/string/ordinal/mixed IDEA and both bootstrap
paths), zero/negative/extended missing values, option abbreviations, active
legacy `tol1()`, empty and random if/in expressions, and forwarding wrappers.
They require `r(198)` and the tolerance-specific diagnostic. Data, types,
metadata, caller objects, RNG state and an existing output file are preserved.
The internal LP and TE contracts must also reject missing tolerances.

Positive parser checks include very small and very large finite values without
requiring numerical optimization at those extremes. Historical controls compare
normal scores for engine/RTS/orientation/stage/syntax/tolerance combinations and
all model profiles, using the same bootstrap seeds. Simple LP and radial DEA
cases have independent analytic answers. A pre-existing internal TE-contract
IPM stage-2 failure at `1e-10` is reproduced and compared with the baseline:
its reported missing score/status is retained, and its finite scores are checked
against the analytic answer. This is a separate numerical limitation, not a
successful solve. Legacy `tol1()` is tested as a standalone alias; conflicting
`tolerance()`/`tol1()` and obsolete `tol2()` are rejected with preservation checks. Literal help examples exercise invalid
values with existing results and return the effective valid tolerance.

## Model-specific missing values and IDEA endpoints

```sh
python3 tests/test_model_missing.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-model-missing \
  --baseline-ado /path/to/pre-F14/ado
```

Requires NumPy and SciPy. The optional baseline is core 2.0.41. Dedicated
DDF/SBM-undesirable cases compare complete-case calls with physically restricted
data and independently constructed HiGHS LPs. They cover both syntaxes, all
three engine requests, applicable RTS values, one/two dimensions, individual
and combined custom directions, ordinary/extended missing codes, if/in,
`ddfbeta()`, replacement, saved samples, returned scalars and efficiency arrays.
The standalone SBM-undesirable wrapper must agree with the core command.

Interval cases cover eight endpoint combinations, both orientations and syntax
forms, RTS variants, unused nominal missing values and nominal endpoint fallback.
Missing/reversed selected endpoints must fail before changing existing results
or files; invalid values outside the selected complete sample are ignored.
Shared nominal variables with an ordinary input/output role remain required.
Empty samples preserve results; unrelated bad/direction options do not alter
ordinary TE samples. Healthy historical scores must remain exactly equal.
Historical failures and silent acceptance are reproduced as separate cases.
Literal help examples also exercise missing directions, bad outputs, unused
nominal values and endpoint rejection. The separate mixed/bounded/ordinal suites
check that existing legacy IDEA behavior remains intact.

## Sample selection in extension commands

```sh
python3 tests/test_sample_selection.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-sample-selection-validation
```

This checks `dea_bootstrap`, `dea_cooperative`, `dea_crossgame`, and `dea_fuzzy`.
Qualified calls are compared with calls on physically restricted data, using
the same order, options, and bootstrap seed. Every generated result column and
returned scalar/local is compared, excluded observations must have missing
results, and original data/order must be preserved. An excluded, highly
productive DMU makes reference-set leakage detectable. Additional single-input
CRS cases compare scores with an analytic frontier calculation.

Coverage includes both available syntax forms, numeric/string conditions,
equality and inequality operators, functions containing commas, combined
`if`/`in`, negative observation ranges, varlist wildcards/ranges, numeric/string
and automatic identifiers, missing model/bound values, unrelated missing data,
all supported RTS/orientation combinations, game strategies, and explicit
replacement of a result variable also used in the condition. Empty samples
must return `r(2000)` without changing data or existing results.

The runner also checks unchanged full-sample outputs and reproduces the old
sample-selection defects. By default it uses Git HEAD as the pre-fix control;
after committing the fix, pass `--baseline-ado /path/to/pre-fix/ado` to retain
those old-defect checks. The source snapshot, hashes, do-files, and logs are
retained in the output directory. `--only REGEX` selects cases. These checks
validate sample handling; they do not establish every estimator's theoretical
correctness or fix the separate bootstrap option inside the core `dea` command.
The historical control still reproduces crossgame's missing singleton peer
averages. Current aggressive/benevolent/all calls instead require `r(2001)`
and unchanged data when only one DMU remains. Singleton Nash remains supported.

## Crossgame minimum sample by strategy

```sh
python3 tests/test_crossgame_sample_size.py --stata "$DEA_STATA" \
  --output ../dea-test-results/crossgame-sample-size --baseline-ado /path/to/pre-R04/ado
```

The optional baseline checks binary equality of accepted singleton-Nash and
multi-DMU scores and returned counts. Cases cover both syntaxes, all strategies
and the omitted/default strategy, numeric/string/automatic identifiers,
whole/qualified samples, ordinary and extended missing inputs/outputs,
replacement and qualifiers using previous result columns. Empty samples retain
`r(2000)`; peer-excluding singleton calls require `r(2001)` before data changes.
Private entry counters verify that rejected calls never start the computation.
Two-DMU single-input/output fixtures have an independent ratio-frontier oracle;
duplicate IDs and identical observations still count as two observations.
`matrix`, `detail` and `debug` are exercised in combination for sample handling
and score preservation; this suite does not certify their display behavior or
the general game-theoretic model.

## Core bootstrap and independent LP audit

```sh
python3 tests/test_core_bootstrap.py \
  --stata "$DEA_STATA" \
  --baseline-ado /path/to/pre-F03/ado \
  --output ../dea-test-results/dea-core-bootstrap-validation
```

Use a frozen source directory from before the core bootstrap repair for the
required historical controls. This suite compares core and dedicated bootstrap
calculations across both syntax forms, all RTS/orientation choices, and initial
simplex/IPM/auto engines. It checks the original score/reference contract,
physical-subset equivalence, bias identities, actual confidence intervals and
standard errors, saved data, replacement, unsupported-option rejection, and
the smoother's support near the efficiency boundary.

For a separate numerical audit, use a Python environment with NumPy and SciPy:

```sh
python3 tests/verify_bootstrap_lp.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-bootstrap-lp-audit
```

The audit instruments a private source copy to record actual pseudo-data draws.
SciPy HiGHS independently solves every initial and pseudo-frontier LP for one-
and two-input/output fixtures. NumPy independently recomputes bias, correction,
standard errors, and inverse-CDF pivot quantiles. Instrumented copies must never
be distributed as the package. Passing these checks does not establish nominal
confidence-interval coverage for every data-generating process.

F03 intentionally changes bootstrap inference: the old core approximation is
replaced, and the dedicated smoother's boundary treatment and CI quantiles are
repaired. Historical controls therefore compare original scores for these two
profiles. Output-protection tests separately compare every inference column
against a current, unreplaced control; they do not require old incorrect
corrections to remain unchanged. The SW2007 profile likewise compares original scores with historical controls
and compares all new corrected outputs with a current control. Its separate
legacy-procedure tests compare the complete old numerical results.


## Sequential SW2007 and preserved joint procedure

Observation-file saving has a dedicated regression suite:

```sh
python3 tests/test_sw2007_saving.py --stata "$DEA_STATA" \
  --output ../dea-test-results/sw2007-saving --baseline-ado /path/to/pre-R03/ado \
  --external --intel
```

Omit `--external --intel` to run without the separately installed public backend.
`--intel` applies only to external cases. The optional baseline reproduces the
formerly ignored `saving()` request and compares every numeric estimation
matrix/scalar and internal output column as binary doubles. Current saved versus
unsaved fits must preserve values, equation stripes, macros, `e(sample)` and RNG
state (temporary identifiers are normalized for metadata comparisons).
Cases span both procedures, all supported RTS restrictions, orientations, syntax
forms, supplied scores, panel periods, `nobcorrected`, selected samples and paths
with spaces, punctuation and Korean characters. External score columns are also
compared directly to the unmodified published command. File guards, replacement,
protected/invalid names and private final-write faults check caller-data recovery.
The files hold observation results; inference matrices remain in `e()`.

```sh
python3 tests/test_sw2007.py --stata "$DEA_STATA" --output ../dea-test-results/sw2007-current
python3 tests/test_sw2007_failures.py --stata "$DEA_STATA" --output ../dea-test-results/sw2007-failures
python3 tests/verify_sw2007_lp.py --stata "$DEA_STATA" --output ../dea-test-results/sw2007-lps
```

The last command requires NumPy and SciPy. It independently solves 26,880
initial/pseudo-frontier LPs across both orientations, four RTS restrictions and
one/two-input/output datasets. It also recomputes inverse-score correction,
bootstrap covariance, empirical percentile intervals and pivot intervals.

`test_sw2007.py` covers the corrected-score refit, independent B1/B2 counts,
level(), tau(), if/in and equality syntax, panel periods, supplied scores,
output suppression, core integration, reference saving/replacement, default
1,000/2,000 replications, invalid inputs, and extreme truncated-normal draws.
Use `--baseline-ado /path/to/pre-integration/ado` to compare the complete
`procedure(joint)` outputs with the historical implementation. Use `--only`
to select case names. The fault tests verify whole-replication rejection,
retry counts, and preservation of existing columns after persistent failure.
For a newer frozen baseline that already defaults to sequential estimation,
add `--baseline-procedure joint`. Use a baseline with the intended solver
version: an exact comparison against a much older package also measures all
intervening LP repairs, not just changes to the bootstrap procedure.

For the separately installed published implementation:

```sh
python3 tests/test_sw2007.py --stata "$DEA_STATA" --output ../dea-test-results/sw2007-public --external --intel
python3 tests/test_sw2007.py --stata "$DEA_STATA" --output ../dea-test-results/sw2007-aligned --external --intel --align-boundary
```

`--intel` is for the x86_64 Mac plugin under Rosetta; omit it when native Stata
and the installed plugin already match. Reference version/hash, logs and the
boundary-alignment patch are retained. No installed reference file is edited.
The unmodified VRS-input fixture deliberately confirms a boundary classification
difference (LP scores differ by about 1e-15) rather than asserting false equality
of the downstream estimates. The second run applies one documented rounding
line to a private reference copy and tests numerical agreement on that common
sample. Other original-reference fixtures test score/coefficient/covariance and
bootstrap-draw agreement directly. Wrapper versus direct-call comparisons
always use the same reference variant. Public interpolated percentiles differ
from this package's empirical inverse-CDF convention.

These checks establish finite-case numerical and control-flow behavior. They
are not a general coverage guarantee or validation of a clustered panel model.

Boundary selection and confidence-interval conventions have a separate audit:

```sh
python3 tests/test_sw2007_conventions.py --stata "$DEA_STATA" \
  --output ../dea-test-results/sw2007-conventions \
  --reference-ado /path/to/installed/simarwilson.ado
```

Omit `--reference-ado` to test only the internal convention. With the public
2.22 source, this checks 77 combinations of bootstrap counts (2 through 399)
and confidence levels (0.1 through 99.9), including tied and nonlinear draws.
Python independently calculates both empirical inverse-CDF endpoints and the
public command's weighted adjacent-order-statistic rule. The original public
percentile program is extracted verbatim; no installed code is edited and its
DEA plugin is unnecessary. Different valid conventions are recorded, not forced
to agree. The full reference source and SHA-256 are preserved with the results.

Eight further checks run the real internal score handling and initial truncated
regression with exact boundary scores and scores at half/twice the requested
tolerance, both orientations, and `tau(1)`/`tau(2)`. A test-only interception
stops before bootstrapping and compares the selected sample and initial fit
with an independently constructed regression. This interception tests boundary
selection; the Algorithm 2 and LP suites above test subsequent calculations.

The same-data regression retry has a focused fault-injection suite:

```sh
python3 tests/test_sw2007_retries.py --stata "$DEA_STATA" \
  --output ../dea-test-results/sw2007-retries
```

It forces one short-fit failure after an otherwise valid fit, verifies recovery
without redrawing, and compares score columns, 15 estimation matrices and RNG
state exactly against an unmodified control. Persistent failure must terminate
within the existing pseudo-sample budget and preserve caller columns. It also
verifies that an already larger iteration limit is not replaced by a smaller
retry limit. The older
failure suite also checks invalid pseudo-data and LP failures. No installed
files are instrumented.

Statistical validation must report independently generated datasets as its
Monte Carlo units, rather than treating DMUs sharing a frontier as independent
replications. Report bias/RMSE, empirical dispersion versus bootstrap standard
errors, interval coverage with Monte Carlo uncertainty, and estimation failures.
Keep failure-inclusive coverage separate from coverage conditional on success.
Use a known true frontier and coefficients; compare unmodified public results
and boundary-aligned diagnostic copies separately. Numerical tests alone cannot
validate heteroskedastic or cluster-dependent inference.


## Extension input/output syntax

```sh
python3 tests/test_io_syntax.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-io-syntax \
  --baseline-ado /path/to/pre-F06/ado
```

Requires only the Python standard library, Git and Stata 16+. The optional
baseline must be a snapshot from before the seven remaining extensions adopted
`_dea_parse_io`. The recorded run has 502 cases, including 25 historical
controls/reproductions; without a baseline, it runs 477 current-code cases.

The seven profiles are `dea_crosseff`, `dea_nash`, `dea_network`,
`dea_malmquist`, `dea_window`, `dea_congestion`, and `dea_sbmundesirable`.
Checks cover Form A/B and the three supported positional-output forms;
equality/inequality, string and function conditions; `if` plus `in`; relative
observation ranges; input/output wildcards and ranges; missing data; string
IDs; output replacement and protection; invalid/duplicate options; and 92
model-option combinations. Panels include three periods with an excluded
productive DMU, so reference-set leakage can change the answer.

Each qualified run is compared with a physically restricted control dataset.
Scores, returned scalars/macros/matrices, matrix labels, and original data
must match exactly. The delegated core's `r(cmdline)` is excluded from string
equality because its temporary sample-marker name legitimately differs.
An entirely missing ordinary result is not accepted as a successful match.
Historical qualifier-free results also require exact equality. This is a
syntax/sample regression audit, not independent certification of the model
objectives, solvers, or game-theoretic interpretation. This suite uses full
option names; abbreviation coverage is in the separate suite below.

## Extension option abbreviations

The suite also covers the `dea_lp` forwarding wrapper: its required positional
specification, named option prefixes, retained `*` options, and duplicate
rejection after delegation. Its radial-only complete-command fixture is
complemented by stage-2 checks in the LP compatibility and stage suites.

```sh
python3 tests/test_option_abbreviations.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-option-abbreviations
```

Requires only the Python standard library, Git and Stata 16+. Optional
`--baseline-ado /path/to/pre-F24/ado` reproduces the historical parser fixes;
its exact numerical comparisons also require matching numerical kernels.
The suite reads packaged extension and
compatibility-wrapper syntax blocks and their help's marked minimum spellings. It compiles the exact
blocks into isolated Stata parser probes, checks every accepted prefix length
and its resolved value, checks rejection or passthrough of duplicate/short
options according to the parser contract, and verifies that the old accepted prefix sets are retained.

Separate complete-command cases compare full names and documented abbreviations
in Form A/B (where supported), including equality/string/function conditions,
`in()` with `int()`, `or()` with the `ort()` alias, returned scalars/macros/matrices,
original-data preservation, and protected replacement on errors. Fixed seeds
make bootstrap result comparisons exact. `r(cmdline)`/`e(cmdline)` is excluded
from text equality because temporary sample-marker names can differ. Historical
full-name scores must be identical. The 14 former help/parser mismatches and
fuzzy's former silent `ort()` override are reproduced on the baseline.

Parser probes test option routing, not model calculations. Complete-command
checks test equivalence for the supplied fixtures; they do not certify every
model's numerical or statistical properties. The dedicated SW2007 `saving()`
path is implemented and has its own file-content and failure-path tests above.

## Fuzzy returns-to-scale aliases

```sh
python3 tests/test_fuzzy_rts.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-fuzzy-rts \
  --baseline-ado /path/to/pre-F05/ado
```

This script requires NumPy and SciPy. It independently constructs the alpha-cut
pessimistic/optimistic LPs and solves them using HiGHS. It covers all eight RTS
spellings, both orientations, one/two inputs and outputs, six bound specifications,
five alpha levels, and full/if/if+in samples. It checks aliases and returned
metadata, uppercase spellings, alpha ordering, interval nesting, crisp limits
against both core engines, replacement, and protected variables. The optional
pre-F05 baseline snapshot checks that six already-correct RTS spellings remain
within 1e-6 of their prior results while the two reversed aliases change.
Version 1.0.5 deliberately changes precision and scaling, so historical
bit equality is no longer required; the independent oracle threshold remains
1e-6. Alias/canonical comparisons and replacement/current-control comparisons
still require exact equality.

`nirs` means sum(lambda)<=1 and maps to `drs`; `ndrs` means sum(lambda)>=1
and maps to `irs`. Results obtained using the two aliases before fuzzy version
1.0.4 must be recomputed. This audit checks the documented interval LPs; it does
not establish equivalence to every fuzzy-DEA formulation in the literature.

The F05 audit exposed two pre-existing precision failures:
`1d_input_outputs_ifin` and `2d_output_inputs_if`. Version 1.0.5 resolves them
with the original 1e-6 oracle threshold unchanged. The full 77-case suite now
passes, with 18,240 independently constructed LPs and maximum score difference
1.02e-8 in the recorded Stata 16.1 run.

## Fuzzy precision and recovery

```sh
python3 tests/test_fuzzy_precision.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-fuzzy-precision
```

Requires NumPy and SciPy. Tests cover actual tolerance delivery to all four
LP branches, explicit-tolerance counterexamples, extreme changes of units,
random positive triangular fixtures, invalid tolerances and bounds,
single-DMU and zero-column cases. The independent LP checks use the original
unscaled data. Fault injection in private ado copies tests primal/objective
validation, nonconvergence, missing parameters, consistent/inconsistent
failure statuses, whole-grid retries and terminal-error data preservation.
Injected failures are control-flow tests, not evidence of real solver failure
rates. No installed or production ado is modified by the test runner.
The precision suite also checks interval ordering and nesting at 1e-8.
The output-protection and sample-selection suites use bounded historical
fuzzy comparisons because of the deliberate numerical change; current
replacement and qualified/physically restricted comparisons remain exact.

## Fuzzy alpha names and output mapping

```sh
python3 tests/test_fuzzy_alpha_names.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-fuzzy-alpha-names \
  --baseline-ado /path/to/pre-F11/ado
```

Requires NumPy and SciPy. The optional baseline is a package snapshot immediately
before the F11 naming fix (fuzzy 1.0.6), including its local numerical fixes.
It checks all RTS spellings, both orientations, crisp/partial/full bounds,
selected samples, reversed alpha order, decimal and numlist step syntax,
very close levels, rounding boundaries, and 1,001 levels in a single suffix
group. Independent HiGHS LPs use each actual alpha. Historical controls compare
whole supported grids and individually calculated close levels exactly.
Returned lists, labels, double storage, default/legacy names, explicit replacement,
automatic IDs, protected inputs, invalid/repeated levels, and the 32-character
name limit are checked. Rejection tests compare the entire dataset and metadata.
The help suite executes the new close-alpha example directly from the help file.
Stata's numlist parsing precision still applies: values parsed identically are
rejected before outputs change. A rounded variable suffix is not the alpha used
by the LP; read `r(alpha)`, `r(pessvars)`, and `r(optvars)` as parallel lists.

## Score counts and ranks

```sh
python3 tests/test_score_classification.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-score-classification \
  --baseline-ado /path/to/pre-F07/ado
```

Uses the Python standard library. The optional baseline must be a pre-F07
snapshot, not the older published commit. Explicit score fixtures check the
frontier, inclusive tolerance endpoints, competition ties, ordinary/extended
missings, excluded observations and all-missing results for all twelve common
display models. The real-command grid covers DDF, additive, undesirable,
radial super-efficiency and super-SBM across three engines, four RTS choices
and both orientations, with varying stages, syntax forms, qualifiers, numeric
or string IDs, saving, replacement and DDF directions. Other common models
serve as controls. Raw scores, samples, returned matrices, saved scores and
original data are checked; console ranks are parsed and verified independently.

Original counterexamples and duplicate frontier DMUs have explicit expected
counts/ranks. Private-package failure injection checks that an actual IPM retry
refreshes the final counts, including `r(N_super)`. The LP kernels are unchanged
by F07. This suite validates score classification and reporting, not optimality
of every underlying model. F25 deliberately changes the undesirable-output
`yb`/`xyb` LPs, so their historical scores are no longer equality controls;
the `x` comparison allows 1e-6 numerical variation. The independent suite below
checks the corrected LP. Nonpositive tolerance handling in other core models
remains a separate review item.

## Undesirable-output DDF

```sh
python3 tests/test_undesirable_ddf.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-undesirable-ddf \
  --baseline-ado /path/to/pre-F25/ado
```

Requires NumPy and SciPy. The optional baseline is a pre-F25 package snapshot.
The suite independently constructs the projected-input/good-output/bad-output
LPs and solves them with HiGHS. It checks three dimension combinations, all
three directions, CRS/VRS/NIRS/NDRS, all engine choices, RTS aliases and default
direction, both core syntax forms, qualifiers, stages and orientations. A
second command with private instrumented engines must reproduce each production
score exactly and exposes lambda/beta for a separate primal-residual check.
Both production engine routes use Stata's LinearProgram for this model; the
HiGHS comparison is therefore essential to an independent numerical check.

Score errors must be below 1e-6 and normalized primal residuals below 1e-7.
Unit changes, individual zeros, a redundant zero bad-output column, extended
missing values, singleton and unbounded examples, negative inputs, zero
directions and invalid tolerances are covered. Tests also check saved scores,
returned matrices, input preservation, automatic retry, and rejection of
controlled invalid solver outputs. Those injections test control flow, not
real-world solver-failure probabilities. All instrumentation stays in the
output directory; installed and production ado files are never instrumented.

The historical `yb`/`xyb` sign error is reproduced explicitly. Historical `x`
scores are compared with a 1e-6 tolerance; any larger change must be explained
by the independent oracle showing a pre-existing numerical error and a more
accurate current result. This suite checks the documented finite-sample LP and
does not establish equivalence among alternative weak-disposability models.

## VRS observation mapping

```sh
python3 tests/test_vrs_observations.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-vrs-observations \
  --baseline-ado /path/to/pre-F08/ado
```

Uses the Python standard library. The optional pre-F08 snapshot reproduces
duplicate-ID averaging, empty zero rows, and ID-sorted matrix order. The grid
checks fifteen identifier patterns across three engines, two orientations and
two stages. Numeric, fractional/large, repeated, missing, string, Unicode,
long and automatic labels are checked against otherwise identical calls with
unique row identifiers. Qualifiers, reversed data order, missing model values,
both syntax forms, replacement, saved scores, reference weights and peer
metadata are included. Every VRS row must match the generated decomposition
variables and the new source-row map exactly; the raw efficiency scores and
RTS counts must match the unique-ID control.
Twelve pre-F08/current pairs also require bit-exact preservation of the raw
score and generated decomposition variables, independent of matrix row names.

Separate cases check extended missing auxiliary results without averaging or
zero filling, real IPM retries after test-only first-attempt invalidation,
ordinary and environmental bootstrap, and removal of stale maps after CRS.
Neither the production LP kernels nor installed ado files are instrumented.
This suite validates observation identity and result assembly, not the VRS
methodology itself. The separate DMU-label suite covers literal quotes and other
special characters; those display failures were corrected in core dea 2.0.50.


## FDH and IDEA results

```sh
python3 tests/test_fdh_idea_returns.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-fdh-idea \
  --baseline-ado /path/to/pre-F09/ado
```

Uses the Python standard library and an optional pre-F09 snapshot. The grid
covers FDH spellings, engines, orientations, RTS labels, both syntax forms,
qualifiers, output replacement and saving. An exhaustive positive-data
dominance oracle checks FDH scores. Interval cases vary input/output bounds,
four RTS restrictions and three engine selectors. Legacy cases cover numeric
and string bounded data, ordinal inputs/outputs and explicit xord/yord calls.

Checks compare returned statistics to the actual score columns, verify selected
observation order and saved values, protect existing VRS variables, and reject
unsupported FDH RTS labels before changing outputs. Test-only missing results
exercise partial/all failures, extended missing values and exact tolerance
boundaries. A separate forced simplex failure verifies a real legacy IPM retry.
Pre-F09 FDH/legacy comparisons require bit-exact preservation of raw scores.
Point-interval historical comparisons allow 1e-6 for the tighter R09 convergence
setting; the separate interval oracle still requires 1e-7 accuracy. Historical
implicit VRS comparisons preload the simplex module needed by the erroneous
old decomposition, with cold-session failures tested separately.

This suite validates result contracts and the specified historical comparisons.
It does not certify interval or ordinal formulations. Ordinal activation and
combination validation are covered separately by the IDEA option tests below. Installed and
production ado files are never instrumented.

A separate one-DMU case now requires both endpoints to equal 1 (the 2.0.56
formulation correction). The optional pre-F09 control retains the old error:
CRS/input optimistic 2.25, with missing optimistic scores for CRS/output and
both VRS orientations. Interval grid cases export high-precision score CSVs.
The independent formulation checks are described below.


## IDEA activation and option combinations

```sh
python3 tests/test_idea_options.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-idea-options \
  --baseline-ado /path/to/pre-F10/ado
```

Uses the Python standard library. The optional pre-F10 snapshot reproduces
ignored `xord()`/`yord()` and repeated ordinal transformations. Valid calls
compare the actual score columns, matrices and summary returns with unique
`ordinal()` controls across engines, RTS restrictions, orientations, syntax
forms, qualifiers and both stages. Tests include wildcards, repeated ordinal
names, output replacement, saved columns, reversed data order, missing ordinary
model values, numeric/string bounded dimensions combined with distinct ordinal
variables, monotone relabelling and actual IPM retries.

Numeric lists containing point placeholders now match the legacy control.
Conflicting models must return r(198) before result replacement or automatic IDs. Tests compare the
entire dataset, existing disk file, labels and characteristics after rejection.
Variable membership, input/output roles, bounded/ordinal overlap, malformed bound
list lengths and protected result names are checked separately. Historical
comparisons retain bit-exact scores for unaffected ordinal-input and
bounded-plus-ordinal-input calls. F30 deliberately corrects ordinal-output
scores; those paths are checked by the independent oracle and permutation
tests below, rather than requiring equality with erroneous historical scores.

Legacy `bounded()` plus distinct ordinal variables and numeric mixed bounds
remain supported. Non-point numeric mixed cases and ambiguous placeholders
are checked by the mixed-bounds suite below. These tests verify dispatch and data preservation;
they do not establish the mathematical validity of the legacy transformations
or replace the independent numeric interval tests below.


## Bounded string syntax and validation

```sh
python3 tests/test_bounded_syntax.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-bounded-syntax \
  --baseline-ado /path/to/pre-F29/ado
```

Uses the Python standard library. Both public syntax forms admit string bounded
inputs, outputs or both. The optional pre-F29 snapshot reproduces the former
r(109) failures and supplies private-worker controls for combinations of
ordinal input/output roles, four RTS restrictions, two orientations and three
engine selectors. Historical controls are used only when ordinal outputs are
absent. Ordinal-output controls use the current worker after F30; the separate
ordinal-order suite independently checks their transformations. Without a
historical snapshot, all controls use the current worker. Canonical interval
scores must match their respective controls exactly.

Tests cover selected rows, saved columns and summary returns, automatic IDs,
result replacement and source-data protection. Decimal/scientific notation,
Unicode whitespace and long strL intervals must agree with canonical intervals.
Numeric point bounds agree with matching string point intervals. A separate
one-input/one-output CRS oracle computes the maximal transformed productivity
ratio directly, without the legacy worker or an LP solver. Malformed,
missing, overflowing and reversed selected bounds must fail with r(459),
preserving the entire dataset, labels, characteristics and existing disk file.
Invalid bounds outside the effective sample are ignored. Unmarked string model
variables remain prohibited in ordinary DEA/FDH/DDF/SBM calls. Core help examples
also exercise string inputs and outputs in both forms through the help suite.

Historical-score agreement alone is not evidence of invariance to reordering.
The F30 suite below tests corrected ordinal transforms and observation-order
invariance explicitly. Numeric lists combined with legacy ordinal/bounded
options are covered by the mixed-bounds suite below.


## Ordinal transformation and observation order

```sh
python3 tests/test_ordinal_order.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-ordinal-order \
  --baseline-ado /path/to/pre-F30/ado
```

Uses the Python standard library and licensed Stata 16+. The optional pre-F30
snapshot reproduces the original order-dependent scores in both orientations
and all three engine selectors. Corrected scores must match a four-DMU analytic
solution through all 24 permutations, four RTS restrictions, both orientations
and three engines. Both stages appear in the permutation blocks.

A separate grid covers single/multiple ordinal outputs, ordinal inputs and
outputs together, mixed string-bounded and ordinal dimensions, ordinal-only
outputs and input-only controls. Every combination of profile, RTS, orientation,
engine and stage is checked under original, reversed and seeded shuffled order.
Both syntax forms, aliases and duplicate ordinal names, if/in selection,
ordinary missing values, output replacement, saved data, observation-aligned
return matrices and summary statistics are included. In-based samples are
permuted within the specified positions to retain the same effective sample.

Only isolated test copies export each transformed dataset before the ordinary
DEA call. An independent Python oracle checks every transformed cell from the
original values and the evaluated DMU's fixed threshold. Fixtures include ties,
all-equal ordinal outputs, negative ranks, close double values near 1e12 and
large ranks near 1e20. Identity cases include duplicate, string and missing DMU
labels; separate cases check missing ordinal values and singleton mixed data.

Score invariance and analytic agreement use 1e-7; transformed cells must match
exactly. This verifies the comparison and mutation fix. It does not certify all
IDEA methodologies; numeric interval formulation has its own tests below.


## Numeric bounds and legacy mixed data

```sh
python3 tests/test_mixed_bounds.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-mixed-bounds \
  --baseline-ado /path/to/pre-mixed-2.0.39/ado
```

Uses Python's standard library and licensed Stata 16+. The pre-mixed snapshot
must include F30's fixed ordinal comparisons. Optional historical controls
encode numeric bounds as `[lower,upper]` strings and use `bounded()` plus
`ordinal()`; without the snapshot the current legacy worker provides controls.
Independent Python calculations verify the transformed values from the raw
fixture, in addition to exact equality between numeric and string encodings.
All exports and failure injection occur only in isolated test copies.

The grid covers eight input/output bound patterns, input/output/both ordinal
roles, four RTS restrictions, two orientations and three engines. Both syntaxes,
two stages, duplicate aliases, effective models, reversed order, if/in,
ordinal missing values, unused missing nominal values, replacement and saving
are distributed across it. Every score column, source row, return statistic,
saved dataset and numeric-dimension metadata must agree with its control.

Separate tests cover mixed numeric/string representations on the same side,
shared and cross-dimension bound sources, double endpoints applied to narrower
nominal storage types, omitted endpoints, invalid selected bounds, ambiguous
placeholders, duplicate model roles, endpoint/output name collisions, stale
metadata, permutation invariance, using()/automatic IDs, endpoint names matching
VRS result variables, and actual IPM retry. Rejection must preserve
the full dataset, labels, characteristics and existing destination file.
This suite verifies equivalence to the legacy transformation; numeric-only
interval LPs are independently checked by the interval formulation suite below.

## Independent LP/MILP and DEA compatibility entry point

```sh
python3 tests/test_lp_compatibility.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-lp-compatibility \
  --legacy-lp /path/to/old/lp-and-ldea-library-directory \
  --legacy-dea /path/to/old/dea_lp-6.1.0-directory
```

The legacy directories are optional and are read without modification. The
runner copies only the package's allowed ado/help/runtime files and uses a
fresh Stata process and isolated installation directories for every case.
It checks the standalone `dea_lp` entry point before and after loading `dea`,
both syntax forms, solver precedence, RTS/stage combinations, qualifiers,
output protection, and historical VRS/IPM defaults. Explicit equivalent calls
to current `dea` must return identical scores. Historical scores are compared
within an absolute tolerance of 1e-6.

Table `lp` and matrix `milp_bb` tests cover continuous, integer and mixed-integer
solutions with analytic optima, min/max, equality and inequality relations,
bounds, row scaling, selected objective rows, `using`, quoted filenames,
solution files, caller-data/object preservation, and input validation.
Malformed matrices must return a failure status without a Mata conformability
error. Tests distinguish infeasible/unbounded problems from node limits with
and without an incumbent. Small bounded integer optima are independently found
by exhaustive enumeration with fixed seed 1719. Solution/objective checks use
absolute tolerance 1e-6 and integer-value checks 1e-7.

Stata/MP 16.1 reports the redundant-constraint problem max x, -x<=0, x>=0 as
infeasible in a direct `LinearProgram()` call. The test retains that native
control (code 1 or 2 depending on runtime), while all wrapper engines must now
report status 2 from the exact origin/ray check. A native code 1 is not a
mathematically correct diagnosis or a successful solve. General status
certification remains part of the separate IPM investigation.
All three new help files contain literal executable examples and are included
in `test_help_syntax.py` automatically.

## IPM TE stage-2 corrections

```sh
python3 tests/test_ipm_stages.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-ipm-stages
```

Requires NumPy/SciPy. Seven fixtures, four RTS specifications and
their aliases, both orientations, both stages, two requested tolerances and
three engine selectors are exercised. Independent radial/slack LPs and original
production constraints check scores, RTS, projection identities and slack
objectives. Cases with native stage-2 failures require missing projections and
truthful status/counts; they are not counted as successful numerical solves.
`oracle_nonoptimal` separately records cases where the independent LP cannot
certify the fixed approximate radial point. Residual/objective checks for
accepted projections use normalized columns and an absolute 1e-6 threshold.

Private copies inject primal, dual and single-row stage-2 failures, verify that
stage 1 skips stage 2, and force automatic retries. Returns must survive selected
samples, missing/duplicate rows, VRS decomposition, bootstrap, `saveref`, and
file saving. Baseline comparisons require byte-for-byte equality of the exported
IEEE-double radial scores and stage-1 reference weights, in addition to the CSV
checks. The optional `--baseline-ado` mode is a strict historical diagnostic: a
pre-P03 baseline is expected to differ after the 2.0.47 numerical-policy change
and is not a passing criterion for this release. The older direct Mata
signature is also exercised.
The current runner independently checks the shared scaled projection-residual
criterion at the requested tolerance. Native convergence receives that same
tolerance. Failed projections and dual-only radial recoveries are validated as
incomplete results; they are not counted as successful primal LPs. Native retry integration has its own suite below; a general
optimality-certificate policy remains outside these tests.


## Shared LP residual acceptance (P03)

```sh
python3 tests/test_lp_residuals.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-lp-residuals
```

Requires NumPy and SciPy. The runner uses a fresh package copy and isolated
Stata paths. It checks an exactly infeasible boundary system (with a separating
row combination), an explicit feasible witness, three tolerances, three engines,
row scaling and changes of row units. NumPy independently recomputes diagnostics
from returned candidates in original coordinates. Private fault copies test
rejection after backend success, auto retry, stage-2 status 90, and checked dual
recovery. Independent HiGHS and analytic frontiers test both dual orientations;
native `getTol()` assertions test actual option delivery. Near-integer MILPs
verify that rounding cannot publish an infeasible incumbent. Missing candidate
and diagnostic values are preserved for native failures with no solution.
Zero rows, free bounds, nonfinite inputs and overflow near `8e307` are also checked.

`test_simplex_stages.py` and `test_saveref.py` distinguish checked dual scores
from unavailable primal references. Retry tests compare status matrices and
missing patterns with a direct IPM control. Historical congestion values use
1e-6 comparison because TE convergence now uses the requested tolerance;
replacement behavior is still compared exactly with a current clean run.


## Bounded native retries and status evidence (P03)

```sh
python3 tests/test_native_retries.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-native-retries
```

Requires NumPy and SciPy. This isolated suite checks an exact binary-rational
primal/dual certificate for the tight feasible counterexample, independent
HiGHS optima, and bit-for-bit preservation of accepted initial solution vectors.
It exercises structural proofs and their guards, including tiny nonzero
coefficients, both scaling options, all engines, bounded retry counts, original
residual rejection, public history after saving, MILP counters and conflicts.
Fault injection is confined to copied sources. Artificial solver reports test
status policy only; they are not evidence about those problems' mathematical
status. `test_lp_residuals.py` separately checks recovered internal DEA
projections in original coordinates and retains unresolved independent LP
statuses in its evidence. The subsequent core integration is covered by `test_core_retries.py`.


## Core TE/VRS native retry integration (P03)

```sh
python3 tests/test_core_retries.py \
  --stata "$DEA_STATA" \
  --baseline-ado /path/to/pre-core-retry/ado \
  --output ../dea-test-results/dea-core-retries
```

Use the immediately preceding reviewed snapshot (`dea` 2.0.47 / IPM 2.0.9)
as baseline. The suite checks both stages and syntaxes, four technologies,
two orientations, two tolerances, scaling/duplicates and selected observations.
Initially accepted primal scores, stage-1 weights and stage-2 projections must
remain bit-for-bit identical. Recovered rows are not required to match a failed
historical result. Per-phase histories are checked against returned counters,
stage status and the corresponding VRS primal. Private initial-only faults
exercise native recovery in each path, including output dual reciprocity;
conflicts and unresolved unbounded reports must not be hidden by dual fallback.
Saving, bootstrap, contractcheck, auto fallback, repeated calls and a small
user matsize exercise preservation of returned histories. Existing persistent
fault tests disable/corrupt retries as well as the initial solve. The ordinary
`test_ipm_stages.py` suite independently verifies recovered projections using
HiGHS and original data, and records unresolved independent LP statuses.


## Explicit IPM projection recovery (P03)

```sh
python3 tests/test_projection_recovery.py \
  --stata "$DEA_STATA" \
  --baseline-ado /path/to/dea-2.0.48/ado \
  --output ../dea-test-results/dea-projection-recovery
```

The baseline is the immediately preceding reviewed snapshot (`dea` 2.0.48 /
IPM 2.0.10). Omitted and zero allowances must preserve all original binary
scores, references, slacks and histories, including failed rows. Positive
allowances must preserve original scores and accepted projections while
reporting the exact factor used for recovered reference/slack identities.
Independent residual and HiGHS objective checks use that declared factor.
The suite covers both syntax forms, technologies, orientations, tolerances,
selected observations, domain/scope rejection before output replacement,
output-name collisions, saved metadata, bootstrap, shadow checks, a previously
loaded engine, a band too small to change a double, and a shifted LP that
still fails. Native fault injection affects private source copies only.


## TE stage status and solution provenance (R05)

```sh
python3 tests/test_stage_status.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-stage-status \
  --baseline-ado /path/to/dea-2.0.51/ado
```

Requires NumPy/SciPy and the immediately preceding complete ado snapshot.
The tests reproduce the old stage-2 failure fallback in an instrumented copy,
compare unchanged radial scores and successful projections against that snapshot,
and check both syntax forms, RTS aliases, both stages/orientations, selected
rows, string duplicate IDs, optional reference output and saving.
The saved diagnostic columns must equal `r(stage_status)` and
`r(solution_source)` in observation order, without leaking into active data.
Collisions with retained scientific columns are rejected before data changes.
Mixed row-level faults verify that overlapping failures are counted once, and
that an automatic IPM retry replaces all initial simplex metadata. Bootstrap
retains the main fit's failure metadata even when its score is valid.
Non-TE models must not inherit either metadata matrix from a preceding TE call.

Run `test_simplex_stages.py` and `test_ipm_stages.py` alongside this suite.
Their private fault injections verify primal/dual recovery, unavailable
stage-2 references/slacks, source codes, failure counts, saved metadata and
legacy optional Mata signatures. Independent HiGHS LPs check radial scores,
projection feasibility and slack objectives. Natural boundary cases and forced
failures are reported separately; a truthful incomplete result passes the
reporting contract without being counted as a successful projection solve.
`test_projection_recovery.py` also checks the explicit radial adjustment path.

The pre-F12 baseline option of `test_simplex_stages.py` still refers to its
older stage-dispatch change, not to the R05 snapshot above. Omit that option
when running current-only regressions for R05.


## Congestion and SBM-undesirable failures (R06)

```sh
python3 tests/test_specialized_failures.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-specialized-failures \
  --baseline-ado /path/to/dea-2.0.52/ado
```

Requires NumPy/SciPy and the complete immediately preceding ado snapshot.
Each case uses a fresh Stata process with isolated ado paths. Positive fixtures
cover regular, repeated, linearly dependent and changed-unit data; both public
syntaxes, core/wrapper entry points, selected samples, CRS/VRS and accepted
engine selectors are exercised. Independent HiGHS programs verify the actual
congestion/SBM equations and CRS/VRS radial components. An explicit self-DMU
solution verifies feasibility of each original LP. Test-only copies append
`0=1` to produce truly infeasible LPs; these are distinguished from forced
numerical-failure branches and from failures on unmodified data.

The suite checks partial/all failures, ordinary/extended missing values,
undefined denominators, overlapping component failures, finite efficiency
counts, means/extrema, wrapper aliases and saved data. A CRS-only failure must
preserve the main congestion measure. Congestion `n_failed` is a row union
across requested components/ratios; `n_failed_primary` isolates the main score.
Core saving retains all congestion auxiliary columns and SBM bad-output data.
Historical score comparisons exclude no finite current values: every available
value must exactly match the pre-change value; newly unavailable values are
counted explicitly. Failure reporting does not certify solver accuracy or the
methodological interpretation of the specialized model.

`test_model_missing.py` also retains bad-output columns in the expected saved
SBM dataset. Run the shared help, IO, output-protection and TE stage-status
regressions to check that the common return and saving paths remain compatible.


## R07 TE numerical diagnostics and small radial scores

```sh
python3 tests/test_numerical_diagnostics.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-numerical-diagnostics \
  --baseline-ado /path/to/pre-R07/ado
```

Requires NumPy/SciPy and the complete R06 package (dea 2.0.53). Runs both
stages, both orientations, four RTS technologies and aliases, simplex/IPM/auto,
two tolerances, duplicates, selected samples, heterogeneous units, dependent
columns and boundary fixtures. Compares unchanged score/reference/slack/status
files with the baseline and independently recomputes the returned residuals.
Tests explicit radial recovery, missing/dual-only results, stale metadata,
old-engine cache reload, direct malformed audit candidates and small positive
score preservation. Tiny-weight/slack fixtures also verify that cleanup cannot
break the original LP residual criterion: the accepted raw vector is retained
when necessary, without re-solving the LP. Faults exist only in test copies. The small-score cases
separate extraction rounding from pre-existing extreme-scale solver error;
they do not claim an exact optimum for every computed score. Numerical residual
checks likewise do not certify exact feasibility. R07's separate historical
66-case replay retains rational infeasibility certificates in the review archive.


## Core non-TE numerical tolerance

```sh
python3 tests/test_nonte_tolerance.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-nonte-tolerance
```

Requires NumPy/SciPy. Isolated instrumentation captures original submitted
LPs and requested/convergence tolerances across economic, allocative, DDF,
additive, SBM, super-efficiency, multiplier, interval, congestion, and
SBM-undesirable paths. HiGHS independently solves these LPs after positive
row rescaling. Each accepted candidate is also checked against original rows
and bounds in Python, separately from the objective comparison (relative
threshold 1e-6, not a tolerance-based optimality claim).
Both syntax forms, engines, RTS choices, stricter tolerances, duplicates,
selected samples and unit changes are covered. A tiny positive SBM-undesirable
fixture is checked against an independently constructed model. Test-only
corruption checks reject bound/equality/nonfinite candidates as missing in
both solver paths. Passive recording is checked against uninstrumented scores;
restart traces verify the two-attempt limit and the stricter convergence value.
Small multiplier weights are checked through normalization and supporting
inequalities; failed native weight LPs must return missing weights while
preserving a usable radial score. Cached older engine signatures must reload. The suite reports feasible LPs
that our solvers still reject separately from genuinely infeasible problems;
it does not certify all model formulas or all datasets.

## Numeric interval IDEA formulation

```sh
python3 tests/test_interval_formulation.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-interval-formulation \
  --rscript "$DEA_RSCRIPT"
```

Requires NumPy/SciPy and Stata 16+. The optional R comparison requires installed
`deaR`; it sets `RGL_USE_NULL=TRUE` and reads raw efficiencies without the
five-decimal rounding in `efficiencies()`. The optional `--baseline-ado` accepts
the pre-R09 2.0.55 package for singleton error reproduction and cached-code
reload checks. Production and installed package files are not instrumented.

The oracle constructs a complete realized dataset for each endpoint and solves
ordinary radial DEA with HiGHS. It does not reuse production matrices. Checks
cover both orientations, all eight RTS names, both command forms, engine
selectors, stages, full/partial/point bounds, if/in selection, permutation,
duplicate rows, positive and admissible zero components, unit changes, random
data and unused missing nominal values. Point bounds are compared with ordinary
DEA under all engines, full bounds with `dea_fuzzy` alpha(0), and optimistic
scores with the legacy bounded-only transformation under all engines.

A two-DMU experiment enumerates every data corner and interior realizations;
its minimum/maximum ordinary DEA scores must equal the reported endpoints.
Nested intervals and singletons are checked separately. R's
`modelfuzzy_kaoliu(..., kaoliu_modelname="basic", alpha=0)` is run on the same
intervals; output-oriented expansion factors are inverted before comparison.
Column unit normalization in the independent LP/R controls preserves the model.
The comparisons use 1e-7 score accuracy, distinct from the requested LP residual
tolerance. Crisp/legacy controls explicitly request tolerance(1e-9). Test-only
engine assertions verify tolerance delivery through bounded, ordinal and mixed
legacy calls under all engines. These tests cover independent box uncertainty; they do not certify
models with dependent bounds, ordinal restrictions or a fixed unified frontier.

## SBM and super-SBM formulation

```sh
python3 tests/test_sbm_models.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-sbm-models \
  --rscript "$DEA_RSCRIPT"
```

Requires NumPy/SciPy and Stata 16+. The optional R comparison uses unmodified
installed `deaR` with `RGL_USE_NULL=TRUE`, reading unrounded scores. The ordinary
oracle eliminates slacks and uses only reference intensities (plus the
Charnes-Cooper variable when needed). The super oracle uses dimensionless
projection ratios. Neither copies production LP matrices. Tests cover the
Tone (2001) example, random data, all RTS/orientation/engine combinations,
both syntax forms, stages, unit changes, duplicates, permutations, selected
samples, zero components/columns, empty super reference sets, and natural
oriented super infeasibility. Expected missing scores must agree with the
independent solver, not merely with another production engine.

Zero cases are checked against the documented structural-zero model and
analytic examples. They are excluded from deaR comparisons because that
package substitutes small positive quantities for zeros. The combined super
score retains ordinary inefficient scores; raw deaR super scores alone have a
different contract. Domain rejections must preserve existing results and avoid
creating files. Optional `--baseline-ado` checks reloading previously cached
engines; `--source-ado` runs the same checks against a frozen package.
The 3e-6 score comparison threshold is separate from LP feasibility tolerance
and does not certify optimality for arbitrary datasets.

## Malmquist and window formulation

```sh
python3 tests/test_panel_models.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-panel-models --rscript "$DEA_RSCRIPT"
python3 tests/test_panel_lp.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-panel-lp
```

Requires NumPy/SciPy and Stata 16+. `test_panel_models.py` constructs independent
radial LPs and explicitly matches observation IDs and reference periods. It
checks adjacent/global decomposition, all window widths 2–4, both orientations,
supported RTS assumptions, balanced/unbalanced and disjoint panels, decimal and
high-precision time values, unit changes, zero components, selected samples,
row permutation, generated values, summary scalars and every pair-summary row.
Optional unmodified `Benchmarking::malmq` uses all available reference units
(`SAMEREF=FALSE`). Its output-oriented index is inverted to match the command's
improvement-above-one convention. Global and window comparisons use
`Benchmarking::dea` against the explicitly specified pooled references.
`--source-ado` can run the checks against an earlier package snapshot.

`test_panel_lp.py` passively records actual submitted LPs and verifies requested
tolerance, accepted residuals and objectives with HiGHS. It separately counts
independently feasible LPs rejected by the native path. Test-only candidate
corruption and selected-distance failures check component preservation,
all-missing summaries, and failure denominators. Analytic CRS growth and code
reloads are checked separately. Source hashes and raw evidence are retained.
The output-protection suite allows window historical differences below 1e-12
from column scaling, checks the analytic pooled CRS score within 1e-8, and
requires exact equality between clean current runs and `replacevars` runs.
These comparisons do not certify arbitrary panels, native solver behavior, or
TFP interpretations for every RTS/decomposition choice.


## Cross-efficiency weights and diagnostics

```sh
python3 tests/test_crosseff_models.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-crosseff-models --rscript "$DEA_RSCRIPT"
python3 tests/test_crosseff_lp.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-crosseff-lp
```

Requires NumPy/SciPy and Stata 16+. Independent multiplier LPs verify primary
and secondary objectives, returned weights, reconstructed peer ratios,
available-peer means, and stage counters for both orientations and CRS/VRS.
The grid includes unit changes, zeros, duplicates, near-efficient observations,
sample exclusions, permutations and a singleton. Alternative optimal weights
are checked by feasibility and objective value rather than requiring identical
cross-efficiency matrices. The optional unmodified deaR comparison aligns
CRS input Method II and distinguishes its other secondary and VRS conventions.
The native LP audit records actual constraints and explicitly counts feasible
secondary problems that use the disclosed primary fallback. Separate fault
copies check unavailable primary weights, secondary fallback, and invalid
candidates. No faults are introduced into production sources.


## Input VRS cross-efficiency correction

```sh
python3 tests/test_crosseff_vrs.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-crosseff-vrs --rscript "$DEA_RSCRIPT" \
  --baseline-ado /path/to/reviewed/ado
```

Requires NumPy/SciPy and Stata 16+. The `vrscross(limzhu)` grid checks both
syntax forms, all three existing weight methods, three tolerances, mixed
measurement units, zero components/columns, duplicated or reordered rows,
selected samples, near-efficient observations and singleton samples. Independent
LPs check primary and accepted secondary objectives; returned weights reconstruct
the guarded translated ratio, its diagonal, off-diagonal means and diagnostics.
Exact rational evaluation of the published seven-DMU example verifies the
intercept sign and the distinction between translated and ordinary self appraisal.
The optional deaR comparison includes its 37-project dataset, independent primary
scores and unchanged `cross_efficiency()` called with the Stata weight solution.
The latter checks the formula, not independent weight selection. R replaces the
diagonal with ordinary efficiency; this command preserves the translated diagonal.
Alternative optimal weights and secondary conventions preclude blanket equality
of native and R cross-efficiency matrices or reproduction of published rankings.

Invalid CRS/output combinations must preserve existing result data. Test-only
fault copies check failure/fallback accounting, and passive capture checks actual
LPs. Optional `--baseline-ado` verifies bit-identical legacy scores and matrices
against a reviewed predecessor; `--source-ado` supports frozen sources. A passing
suite does not imply every secondary LP succeeds, every ratio exists, or the
paper's output-oriented extension has been implemented.


## Cross-efficiency game calculation and stopping states

```sh
python3 tests/test_crossgame_models.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-crossgame-models --rscript "$DEA_RSCRIPT"
python3 tests/test_crossgame_lp.py \
  --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-crossgame-lp
```

Requires NumPy/SciPy and Stata 16+. Independent LPs verify primary weights,
accepted secondary surrogate objectives, reconstructed peer matrices, averaging
sets and every complete Nash iteration starting from the returned initial state.
The grid covers all strategies, both public syntax forms, unit changes, zeros,
duplicates, row orders, selected samples, near-efficient units and tolerances.
The five-DMU Liang et al. example also runs with short and long iteration limits.
Optional R `lpSolve` constructs and iterates model (5) independently; it is a
separate reference implementation, not a packaged R game-DEA estimator. The
published table uses stopping epsilon 0.001, so a tighter solution is compared
as a numerical diagnostic rather than asserted to reproduce all printed digits.

The LP suite captures actual submissions and tests missing primary weights,
secondary fallback, corrupted candidates, and failed first/later Nash iterations.
A failed Nash attempt must not overwrite the last complete score/matrix state
or be marked converged. `test_crossgame_sample_size.py` verifies minimum samples;
run it without `--baseline-ado` after the numerical changes. Historical exact
weight equality is no longer promised. The output-protection suite checks the
analytic unique 1x1 CRS solution and exact clean-versus-replacement outputs,
with a 1e-6 allowance for the historical native solve in that fixed fixture.

Two-stage network and Stackelberg models (Stata 16+, NumPy/SciPy; optional R
DJL and Benchmarking):

```sh
python tests/test_twostage_models.py --stata "$DEA_STATA" --output ../dea-test-results/twostage-models --rscript "$DEA_RSCRIPT"
python tests/test_twostage_lp.py --stata "$DEA_STATA" --output ../dea-test-results/twostage-lp
```

The model suite builds independent multiplier-dual stage LPs and independent
shared-weight LPs, tests both leaders, validates decomposition faces and returned
weights/intensities, and compares original R packages. DJL centralized stage
weights can select a different optimum: only the system score is compared
unless the decomposition selection rule also agrees. Public option/domain guards,
workspace cleanup, sample selection, units, duplicates, zeros and reloading are
included. The LP suite records submitted models in isolated copies and injects
failures to verify component preservation and missing follower results. Passing
these tests does not mean that all native LP attempts succeed.


`test_cooperative_nash.py` audits the additive fixed-frontier Shapley identity
with independent coalition enumeration, and the historical capped allocation
rule with independent arithmetic. It compares CRS/VRS scores with independently
constructed SciPy/HiGHS multiplier duals and optionally R Benchmarking; captures
actual native LPs and residuals (including explicit self-reference witnesses
and a HiGHS IPM check without presolve for near-duplicate VRS cases); covers units, zero entries/columns, duplicates,
near duplicates, sample syntax, replacement, reload, large identifiers, n>20,
forced score/candidate failures, endpoint adjustments, gain caps and product
underflow. This does not validate a general pooling game or Nash optimizer.

```sh
python tests/test_cooperative_nash.py --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-games-review --rscript "$DEA_RSCRIPT" \
  --baseline-ado /path/to/reviewed-before/ado
```

NumPy/SciPy and Stata 16+ are required; the baseline and R comparison are optional.
Output directories must be new. Historical fault controls demonstrate defects
in an isolated baseline copy; they never alter installed or archived ado files.


`test_fuzzy_models.py` checks the triangular joint-cut convention with independent
multiplier duals for CRS/VRS/NIRS/NDRS and both orientations, full-box corner
extrema on a tiny fixture, alpha nesting/crisp collapse, partial bounds, units,
zero entries, sample selection, complete-grid cohort summaries, missing-cut
failures and quoted identifiers. An optional original deaR 1.5.4 comparison
uses `RGL_USE_NULL=TRUE` for headless operation; it does not modify the package.
`r(498)` numerical aborts, if encountered, are recorded separately and checked
for preservation of the dataset, rather than counted as available DEA scores.

```sh
python tests/test_fuzzy_models.py --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-fuzzy-model-review --rscript "$DEA_RSCRIPT" \
  --baseline-ado /path/to/reviewed-before/ado
```

Use Python 3.10+ with NumPy/SciPy. The output directory must be new. The historical
controls and forced failures affect isolated copies only. Existing fuzzy RTS,
alpha-name and precision suites remain separate regression checks.


`test_environment_models.py` compares SBM-undesirable to an independently derived
LP after eliminating the slacks, undesirable DDF to a multiplier dual, and the
legacy equality ratio to a dual equality LP plus independent radial scores.
It covers RTS, three observed directions, engine spellings, units, zeros,
duplicates, selection, wrapper parity and two interpretation counterexamples.
Original `deaR` checks read SBM efficiency and DDF **beta**, not DDF efficiency.

```sh
python tests/test_environment_models.py --stata "$DEA_STATA" \
  --output ../dea-test-results/dea-environment-review --baseline-ado /path/to/before/ado \
  --rscript "$DEA_RSCRIPT"
```

Use NumPy/SciPy and optional unmodified deaR. Output paths must be new. Missing
native results are reported separately; passing failure-handling checks does
not mean every numerical LP converged. The congestion examples check the
historical ratio's limitations, not an implementation of a Cooper decomposition.


Cooper congestion: `test_cooper_congestion.py` checks the author manuscript's
A--G example, two legacy-ratio counterexamples, nonradial output gains, multiple
optima, singleton samples, China industry data, random/unit-scaled/duplicate/near
samples, both public syntax forms, output protection, saving, and injected phase
failures. It uses independently eliminated-slack SciPy LPs and optional R lpSolve
(the latter is an independent transcription, not a third-party congestion package).
Range comparisons condition on the selected projection and attained optimal face;
main slack objectives are compared with the independently optimized radial value.
Native unavailability is counted separately from a test-contract failure.
The independent SciPy objective threshold is 2e-7; R lpSolve's default-precision
comparison uses 1e-6 and separately reports cells differing by more than 2e-7.
The automobile raw-objective case exposes a roughly 2.3e-7 difference in R's
near-zero technical slack; the original values are retained for review.

```sh
python tests/test_cooper_congestion.py --stata "$DEA_STATA" --output /new/results --rscript "$DEA_RSCRIPT"
```


## Convex joint weak-disposal DDF

```sh
python3 tests/test_weak_disposability.py \
  --stata "$DEA_STATA" \
  --rscript "$DEA_RSCRIPT" \
  --output ../dea-test-results/dea-weak-disposal
```

Requires NumPy/SciPy. Compares active/idle Stata LPs with an independently
formulated multiplier dual and a total/active-intensity primal. It checks
CRS/NIRS equivalence and VRS/NDRS inclusion relative to the existing equality
technology, joint-contraction witnesses, analytic examples, all four RTS and
aliases, observed directions, dimensions, zeros, duplicates, unit/order changes,
syntax forms, stages, engine labels and atomic result publication. Optional R
uses an independent transcription through lpSolve, not an external complete DEA
command. Missing native results are counted separately even when independent
LPs succeed. The existing `test_environment_models.py`, `test_undesirable_ddf.py`
and help/output suites cover the default paths. The new worker uses distinct package-simplex and Stata-native routes; auto
retains accepted simplex solutions and falls back only for unavailable DMUs.
Private copies test both solver paths, equivalent retries, fallback and failures.


Recovery-log contract: `test_recoverylog.py` exercises both command
syntaxes, all four RTS choices, TE stages, SBM orientations, weak-DDF directions,
all three engines, selection/string labels, per-LP context, error cleanup and
row caps. It compares every pre-existing returned numeric matrix/scalar and score
byte-for-byte against a supplied frozen compatible ado directory, with recording OFF/ON.
Use the frozen P3 candidate when checking installation packaging; a production
version before weak-disposal diagnostics has a different return contract.
Private fault copies distinguish core success, adapter rejection, model result,
stage-2 references and native fallback. See ../INSTALL_RECOVERYLOG.txt for command.


`test_te_domain.py` checks ordinary TE's selected-quantity domain, before-output
failure behavior, 0 cells/columns, extended missing values, excluded invalid
rows, extreme positive scales and small positive scores. It compares all existing
numeric returns against frozen P2 with recording OFF/ON and independently solves
the accepted zero-component/zero-column radial LPs. It does not impose this
ordinary-TE contract on FDH, IDEA or super-efficiency. See the installation note.

`test_te_legacy_boundary.py` compares 288 legacy bounded/ordinal/mixed
combinations with the previous candidate, then checks 576 invalid public TE
commands immediately after legacy calls. Internal transformed reference rows
retain the legacy domain; they do not disable validation in later public calls.
Run with `--stata`, `--output`, and `--baseline-ado` as for the domain test.

The synthetic `fixtures/weak_units2.csv` is a 10-DMU, two-input, two-good-output,
two-bad-output regression fixture. It preserves the normalized doubles from
the near-duplicate units2 failure (original simulation seed 202611032). Use it
with `test_recoverylog.py --fixture fixtures/weak_units2.csv` from this folder.
It is a numerical regression example, not observed firm data.


## Exact duplicate recovery after numerical failure

```sh
python3 tests/test_duplicate_recovery.py --stata "$DEA_STATA" --output ../dea-test-results/dea-duplicates
```

The packaged 1,215 conditions derive from five seeds. The 388 newly recovered
conditions require exact original constraints and objective; 827 prior accepted
outputs (including objective and audit record) retain their exact stored expectations. Private observation verifies
actual iteration totals, one aggregate audit record and no recursive reduction.
One-ULP and injected-failure controls reject unsafe recovery. See
[DUPLICATE_RECOVERY.txt](DUPLICATE_RECOVERY.txt). These are finite local tests,
not 1,215 independent empirical datasets or a universal reliability guarantee.


## Regression source selection

Some broad command regression scripts select release files with `git ls-files`.
Run those scripts in a Git checkout, or run `git init` in an extracted complete
package before using them; no staging, commit or remote account is required.
The Stata commands and the self-contained duplicate-recovery test do not require
Git. The numerical fixture tests documented above use the packaged inputs.

## Successful-basis quality regression

Run `python tests/test_success_quality.py --stata "$DEA_STATA"
--output /new/output` as one shell command from the package root. The packaged
fixtures cover original-coordinate public LP checks, lower-bound constants,
exact zero rows and conservative fast bounds. See [SUCCESS_QUALITY.txt](SUCCESS_QUALITY.txt)
for the counts, exact-arithmetic comparisons and finite limits.

The duplicate regression now requires43 previously approximate successes to
reach their exact witnesses, while retaining1172 already-exact results.
Previous failure-only preservation statements in these historical test notes
are subject to the successful-basis quality check in simplex2.0.23.

## Exact original-row recovery regression

Run `python tests/test_row_recovery.py --stata "$DEA_STATA" --output /new/output`
from the package root. This self-contained suite uses packaged fixtures and
Python Fraction arithmetic; no research directories or network are required.
See [ROW_RECOVERY.txt](ROW_RECOVERY.txt) for exact arithmetic,original-constraint,
empty-basis,recovery-depth,affine-bound and MILP checks and explicit limitations.

## DEA stage-link regression

Run `python tests/test_stage_linkage.py --stata "$DEA_STATA" --output /new/output`.
The packaged [STAGE_LINKAGE.txt](STAGE_LINKAGE.txt) describes the exact basis/target
checks, public diagnostics, fault controls, sample order and scope limits.

Cooper2.0.73 integration: `test_cooper_certified.py` tests the shared native/simplex
weight-domain rejection and the explicit16/128selection contract using portable
synthetic fixtures. `test_cooper_congestion.py` observes the simplex LPs in a
private copy and verifies each original rational primal/dual certificate and
exact-face endpoint export with `cooper_certified/exact.py`. Floating SciPy
range equality is not an acceptance oracle for exact optimal faces. Native
LP/oracle/fault coverage and public syntax/output protections remain included.


## Cross-efficiency original-data certificates and recovery diagnostics

Run `python tests/test_crosseff_recovery.py --stata "$DEA_STATA" --output /new/output`.
Portable synthetic fixtures exercise both orientations/RTS choices and all
methods at actual requested tolerances. A private observation copy records
original LPs and integer witnesses. Fraction arithmetic checks accepted
primal/dual points, optimal-face links and unbounded rays; public journal
sums and numerical-primary preservation are checked. Artifacts never
replace installed code. See INSTALL_CROSSEFF.txt for the 78-file contract.


Cross-efficiency native preflight and warning checks:
    python tests/test_crosseff_preflight.py --stata "$DEA_STATA" --output /new/output
The related synthetic fixtures distinguish nonfinite primary preparation
from secondary aggregate overflow. Both stages are exercised at three actual
requested tolerances with matrix/no-matrix calls, unchanged primary fallback,
warning counters and cold native loading. See CROSSEFF_PREFLIGHT.txt for scope.


## Exact intended peer objectives

Run `python tests/test_crosseff_objective.py --stata "$DEA_STATA" --output /new/output`.
The56related fixtures check exact peer sums, original primary primal/dual
identities, original optimal-face constraints and intended-objective
certificates independently with Fraction arithmetic. The52historical cases
include48rounded vectors;four additional cases separate objectives whose
peer sums round to equal doubles. The test enforces existing work/limb and
exchange bounds. See CROSSEFF_OBJECTIVE.txt for scope and native limitations.
Portable recovery now verifies integer intended costs for secondary witnesses.
The preflight regression checks the old aggregate-overflow fixture using the
new finite peer-only builder and retains first-stage nonfinite guards.


## Optional optimal-face diagnostics

Run `python tests/test_crosseff_optima.py --stata "$DEA_STATA" --output /new/output`.
The12sensitive command fixtures ×2engines compare option-off/on original outputs,
independent rational F1/F2 classifications, coordinate/denominator endpoints,
cutoff states and exact stored-weight linkage. Fixtures carry immutable exact
references; no research directory or external solver is needed at runtime.
See CROSSEFF_OPTIMA.txt. The previous full regression remains separate from
this scoped diagnostic test; no claim that all multiple optima are errors.

### Explicit peer CE ranges

`test_crosseff_cerange.py` verifies the `cerange` option against independent
exact references, including limit-only/unbounded/undefined/unavailable states,
summary counts and unchanged native/optima outputs. See
[CROSSEFF_CERANGE.txt](CROSSEFF_CERANGE.txt). Run with `--stata` and a new
`--output` directory; optional `--baseline-ado` compares a complete older package.

The version1.10.3 raw-recovery references include four failed-proposal primary
models certified in original units, while finite-weight normalization obstructions
remain unavailable. Secondary arithmetic-limit cases retain primary weights.
New weights are derived from exact certificates and independently rounded.
