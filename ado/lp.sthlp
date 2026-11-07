{smcl}
{* Help revision 2026-11-05: bounded exact original-row recovery and empty-basis checks.}{...}
{* Help revision 2026-11-04: original-column reconstruction after failed Phase-1 pivots.}{...}
{* Help revision 2026-11-04: cycle-gated certified simplex prices.}{...}
{* Help revision 2026-11-03: wide recovery and exact power-of-two row scaling.}{...}
{* *! version 2.0.5  05nov2026}{...}
{title:lp — independent linear and mixed-integer optimization}

{p 8 16 2}{cmd:lp} {it:coefficient_variables} {ifin} [{cmd:using} {it:filename}]
{cmd:,} {cmd:min}|{cmd:max} [{it:options}]

{pstd}
Successful simplex results receive a basis-quality check. Already-certified
values are retained. A verified defect may be corrected only when the refined
point passes the original primal and conservative primal-dual checks. If
certification is unavailable or refinement fails, the original success is
retained; this does not certify every returned optimum. The public LP adapter
checks original coefficients and bounds, including lower-bound objective
constants. Exact zero rows are handled separately. A conservative fast check
falls back to compensated arithmetic when necessary. No new public option is
required, and native IPM calculations are unchanged. Earlier-result priority
in the failure-only descriptions is subject to this additional check.

{pstd}
The simplex engine also supports bounded exact selection of independent
original rows, including the right-hand side, and analytic checks for a trivial
empty basis. After every earlier recovery path fails, it may solve a reduced
row system once; acceptance still requires checks against all original rows,
columns and bounds. Near dependence alone never authorizes row removal.
The selector is limited to256 rows/columns,rank32,32 expansion components and
one million work units, with additional arithmetic-range guards. Unsupported
cases retain the existing outcome. No new public option or default strict
policy is introduced. These finite checks do not guarantee exact optimality
for every dataset. The packaged ROW_RECOVERY.txt note describes the regression scope.

{title:Input table}
{pstd}Requires Stata 16 or later. Each numeric variable is a decision-variable
coefficient column. After {cmd:if}/{cmd:in}, the first selected row specifies
the objective; subsequent selected rows are constraints. There is no implicit
intercept. Variables are nonnegative by default. This preserves the main
dataset-table interface of the former standalone {cmd:lp} command.

{pstd}{cmd:rel(}{it:variable}{cmd:)} names the string relation column (default
{cmd:rel}). Constraints accept {cmd:=}, {cmd:<=}, and {cmd:>=}; the historical
{cmd:<}/{cmd:>} spellings mean non-strict inequalities. {cmd:rhs(}{it:variable}{cmd:)}
names the numeric right-hand-side column (default {cmd:rhs}). The objective
row's relation and RHS are ignored; RHS is represented as zero in {cmd:r(problem)}.
Selected coefficients and constraint RHS values must be finite. Missing values
are rejected, not silently removed. With {cmd:using}, variables and conditions
are evaluated in that file; the caller's dataset is restored on success or error.

{title:Options}
{phang}{cmd:min} or {cmd:max}: specify exactly one objective direction.
{phang}{cmd:intvars(}{it:varlist}{cmd:)}: decision columns constrained to integers.
The list must be a subset of the coefficient variables. Other columns remain
continuous. Without this option the problem is a continuous LP.
{phang}{cmd:engine(simplex|ipm|auto)}: LP engine, default {cmd:simplex}.
{cmd:auto} tries simplex first. Without a structural proof, any unsuccessful
simplex result is checked with IPM. MILP uses the same policy at each node.
Explicit {cmd:ipm} uses native IPM only; see the bounded retry policy below.
{phang}{cmd:scaling(none|row)}: optional objective/constraint row scaling; default
{cmd:none}. Returned solutions and objective values use the original units.
This option uses exactly reversible power-of-two objective and row multipliers;
A row and its right-hand side share one multiplier. A row or objective is left
unscaled if reversal would lose a stored coefficient or exceed the supported
numeric range. Nonzero maxima normally fall in [1,2).
Compared with earlier max-absolute scaling, this optional policy can change
previously accepted numerical results.
This option controls initial contract scaling. After a numerical failure,
the simplex core can make one fresh-basis retry with row, column and objective
equilibration, including under {cmd:scaling(none)}. The recovered candidate is
checked in original coordinates. Native-IPM retry counters do not count this
internal simplex restart; its iterations are included in the simplex total.
If the established simplex paths still fail numerically, a further artificial-free
basis retry can refine a candidate in original coordinates and check its quality
with rounding-error bounds. Within the selected initial representation, this
failure-only retry retains existing accepted results. It preserves original rows
and makes at most two real-column basis proposals. The extra
check has finite arithmetic limits and does not guarantee recovery or relative
objective accuracy for every LP; {cmd:tolerance()} retains its existing meaning.
A detected repeated Phase-2 Bland pivot state ends only that path. Existing
restarts retain priority; a further retry with basis-error-aware price intervals
is attempted only after an observed repeat and failure of the preceding paths.
Uncertain prices remain unresolved. Its iterations are included in the simplex
total; it does not change the native-IPM counters or certify every LP solution.

{phang}{cmd:lb(}{it:matrix}{cmd:)} and {cmd:ub(}{it:matrix}{cmd:)}: bounds in coefficient
variable order. Lower bounds must be finite; default zero. Upper bounds default
to missing ({cmd:.}), meaning unbounded above. Negative finite lower bounds
are allowed. A fully free variable must be represented by a difference of two
nonnegative variables.
{phang}{cmd:tolerance(}{it:#}{cmd:)}: positive finite LP tolerance, default {cmd:1e-8}.
The old spellings {cmd:tol1()} and {cmd:tol2()} are accepted as aliases for this
single modern tolerance. They no longer specify separate pivot and inverse-matrix
criteria. Conflicting explicit values are rejected; migrate such calls to one
{cmd:tolerance()} request. The solver implementation and default precision differ
from the historical library; numerical identity is not promised.
{phang}{cmd:inttol(}{it:#}{cmd:)} and {cmd:maxnodes(}{it:#}{cmd:)}: MILP integrality
tolerance (default {cmd:1e-7}) and node limit (default 10000). Both must be positive
and finite; {cmd:maxnodes()} is an integer. See {help milp_bb} for stopping rules.
{phang}{cmd:rownames(}{it:variable}{cmd:)}: observation labels for {cmd:r(problem)},
using Stata's {cmd:mkmat} row-name rules.
{phang}{cmd:trace}: show solver details.
{phang}{cmd:saving(}{it:filename}{cmd:)} [{cmd:replace}]: save a solution dataset,
with one row per decision variable and columns {cmd:parameter}, {cmd:value},
{cmd:objective}, {cmd:status}. The filename gets a {cmd:.dta} extension if absent.
Existing files require {cmd:replace}. Failed solves and node-limited runs retain
their status in the file; check it before using values as an optimal solution.

{pstd}
After all earlier simplex paths fail, an observed Phase-1 pivot blockage can
trigger original-column basis reconstruction. Original rows, bounds, and final
primal/dual checks are retained. Earlier accepted solutions keep priority.
All attempted iterations count toward the reported total, including unsuccessful
reconstruction retries. Dimension/range limitations still apply.

{pstd}
After all earlier simplex recovery paths fail, one bounded retry may remove
exact positive power-of-two duplicate rows or columns. Row comparisons include
the right-hand side; column comparisons include the objective and every original
constraint and bound row. The solution is restored to the original variables
and must pass the original primal and primal/dual quality checks. Earlier
accepted common-simplex solutions and reference weights retain priority. All attempted
iterations are included in the reported total, with one aggregate recovery
record per common simplex call.

{pstd}
This retry uses exact stored-coefficient comparisons, not tolerance-based row
deletion. Supported ratios are powers of two with exponents from -400 to 400;
finite-arithmetic checks can further restrict this range. The original
standard-form size limit is n+2m <= 4096, and duplicate comparisons have a
1,000,000-element budget. Unverified or unsupported candidates remain numerical
failures; the retry does not publish a new infeasible or unbounded certificate.
It does not guarantee recovery of every degenerate LP or improve every already
accepted result. User tolerances and native IPM calculations are unchanged.

{pstd}
With {cmd:engine(auto)}, a new simplex recovery can replace a previously
successful native fallback. Scores, reference weights, iteration totals and the
reported backend can therefore change across versions. Earlier-result priority
within the common simplex solver does not promise identical automatic-engine
results across versions.

{title:Results and compatibility}
{pstd}Data, observation order and caller matrices are preserved. No automatic
{cmd:lp.log} is written. Returns are in {cmd:r()}, not {cmd:e()}.
{cmd:r(x)} contains decision values; {cmd:r(objective)} is the objective.
{cmd:r(lprslt)} contains objective {cmd:z}, decision values and nonnegative
inequality slacks {cmd:slack1}, ... in their original constraint order.
Equality rows have no slack column. {cmd:r(problem)} is the selected coefficient
table with the final RHS column. {cmd:r(nvars)}, {cmd:r(nconstraints)}, and
{cmd:r(nslacks)} count the dimensions. {cmd:r(variables)}, {cmd:r(intvars)},
{cmd:r(problem_type)}, {cmd:r(tolerance)}, and {cmd:r(cmd)} describe the call.
There is no {cmd:r(tableau)} final basis tableau: IPM does not supply a simplex
basis. The old library's artificial-variable columns are not reproduced.

{pstd}Inspect {cmd:r(status)} and {cmd:r(status_text)} even when the Stata command
returns normally. For continuous LP: 0 optimal, 1 infeasible, 2 unbounded,
3 invalid input, solver failure, residual rejection, or conflicting status. For MILP: 0 optimal within tolerances,
1 infeasible, 2 unbounded LP relaxation, 3 node limit with incumbent,
4 node limit without incumbent, 5 invalid input or LP failure. An unbounded
relaxation does not certify an unbounded integer problem. Backend diagnostics
are also returned: LP residuals/engine, or MILP node and solver counts.

{pstd}Continuous LP candidates are checked against the original constraints and
bounds after solving, including when {cmd:scaling(row)} is used. For a row
{it:a*x=b}, the scaled violation is {it:abs(a*x-b)/D}; for {it:a*x<=b}, it is
{it:max(0,a*x-b)/D}, where {it:D=max(abs(b),sum_j abs(a_j)*max(1,abs(x_j)))}.
An all-zero row and RHS has zero violation. Each finite bound uses denominator
{it:max(1,abs(bound),abs(x_j))}. All scaled violations must be no greater than
{cmd:tolerance()}; there is no larger hidden acceptance tolerance. Absolute
violations remain in original units. This checks primal feasibility within the
requested tolerance, not exact feasibility, optimality, or infeasibility.

{pstd}For continuous LPs, {cmd:r(max_eq_violation)},
{cmd:r(max_ineq_violation)}, and {cmd:r(max_bound_violation)} are absolute
violations; the corresponding {cmd:r(max_scaled_eq_violation)},
{cmd:r(max_scaled_ineq_violation)}, and {cmd:r(max_scaled_bound_violation)}
are dimensionless. {cmd:r(feasibility_tolerance)} records the acceptance limit.
If no retry recovers an acceptable solution, residual rejection gives status 3
and {cmd:r(status_text)} {cmd:residual check failed}, with missing {cmd:r(x)}
and {cmd:r(objective)}.
{cmd:r(candidate_x)} retains a finite candidate for diagnosis when available;
{cmd:r(reported_status)} retains the selected adapter report (zero after a
successful recovery). Residual diagnostics describe the selected candidate;
{cmd:r(attempts)} retains the separate per-attempt evidence. These candidate/diagnostic returns are for
continuous LPs, not the collection of MILP node relaxations.

On unsuccessful solves values may be missing; a MILP incumbent is not necessarily
optimal. General LP/MILP tools are maintained extensions of the package; their
addition is not a claim of validation by the DEA journal article.

{title:Native retries and status evidence}

{pstd}A successful initial solve that passes the residual check is returned
immediately. Otherwise IPM makes at most three additional attempts:
(1) replace equality rows by both inequality directions;
(2) apply positive row and objective scaling to that representation;
(3) use the original representation with convergence tolerance
{it:min(requested tolerance,1e-12)}. Inapplicable or identical attempts are
skipped. Every candidate is checked against the original constraints, bounds
and requested acceptance tolerance; no RHS or bound is relaxed. The native
iteration limit is not increased. The retry policy also applies with
{cmd:scaling(none)}; that option controls the initial representation.

{pstd}Continuous LP returns include {cmd:r(initial_status)} (first raw solver
code), {cmd:r(native_retries)} (number of extra native solves),
{cmd:r(status_conflict)} (one for contradictory reports), and
{cmd:r(attempts)}. The latter has one row per actual attempt and columns
{cmd:engine variant native checked tolerance iterations residual objective}.
Engine codes are 1=simplex and 2=IPM. Variant 0 is the initial representation;
1, 2 and 3 are the retries above. The checked code is zero only for an accepted
candidate, 90 for a native success rejected by the original-coordinate check,
or a failure code. The tolerance column is the actual convergence setting;
residual is the maximum scaled original-coordinate violation. A missing
residual or objective indicates an unavailable or nonfinite diagnostic quantity.
For MILP, {cmd:r(native_retries)} aggregates extra native solves across nodes,
and {cmd:r(status_conflicts)} counts conflicted nodes. Adapter counts exclude
these additional solves.

{pstd}{cmd:r(certificate)} is {cmd:none} unless a narrow structural proof is
available: {cmd:zero_row} or {cmd:nonnegative_row} proves infeasibility from
an impossible constraint; {cmd:zero_origin_ray} proves LP unboundedness from
an exactly feasible origin and an improving coordinate direction. These tests
use exact signs and zeros of the supplied coefficients, never tolerance-based
zeroing. They are not a general infeasibility/unboundedness certificate system.

{pstd}In the tested Stata/MP 16.1 runtime, a direct native call labels max x,
x>=0, -x<=0 infeasible. The wrapper now reports status 2 from the structural
ray proof while retaining the native report in {cmd:r(attempts)}. Without such
a proof, contradictory infeasible/unbounded reports, or an unbounded report
followed by a finite candidate, give status 3 and no solution. An infeasible
report may be superseded by an accepted native solution. If checking ends in
numerical failure, the result remains status 3. Agreement between solvers is
not an independent mathematical certificate. Status 0 still means a solver
optimum that passes the primal tolerance check, not a separate optimality proof.

{pstd}These retries cover the shared LP adapter, its internal radial DEA
translator and standalone MILP relaxations. From {cmd:dea} 2.0.48, its native
TE and auxiliary VRS paths also use bounded native retries and per-DMU evidence;
see {help dea}. Other native model kernels retain their individual policies.

{title:Examples}
{pstd}Run this complete block; it restores the original data.
{* DEAS_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input double x1 x2 str2 rel double rhs}{p_end}
{phang2}{cmd:3 2 "=" 0}{p_end}
{phang2}{cmd:2 1 "<=" 4}{p_end}
{phang2}{cmd:1 2 "<=" 4}{p_end}
{phang2}{cmd:. end}{p_end}
{phang2}{cmd:. lp x1 x2, max}{p_end}
{phang2}{cmd:. assert r(status)==0 & abs(r(objective)-20/3)<1e-6}{p_end}
{phang2}{cmd:. lp x1 x2, max intvars(x1 x2) engine(auto)}{p_end}
{phang2}{cmd:. assert r(status)==0 & abs(r(objective)-6)<1e-6}{p_end}
{phang2}{cmd:. lp x1 x2, max intvars(x2) engine(ipm)}{p_end}
{phang2}{cmd:. assert r(status)==0 & abs(r(objective)-6.5)<1e-6}{p_end}
{phang2}{cmd:. matrix list r(lprslt)}{p_end}
{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}
