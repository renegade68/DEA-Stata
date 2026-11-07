{smcl}
{* Help revision 2026-11-04: bounded exact duplicate recovery after numerical failure.}{...}
{* Help revision 2026-11-04: original-column reconstruction after failed Phase-1 pivots.}{...}
{* Help revision 2026-11-04: common simplex certified-price recovery.}{...}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.2.4 04nov2026}{...}
{title:Title}

{phang}{bf:milp_bb} {hline 2} solve a small mixed-integer linear program by branch-and-bound

{title:Syntax}

{p 8 16 2}{cmd:milp_bb}, {cmd:c(}{it:C}{cmd:)}
{cmd:integer(}{it:numlist}{cmd:)}
[{cmd:sense(max|min)} {cmd:aeq(}{it:Aeq}{cmd:)} {cmd:beq(}{it:beq}{cmd:)}
{cmd:aub(}{it:Aub}{cmd:)} {cmd:bub(}{it:bub}{cmd:)}
{cmd:lb(}{it:lb}{cmd:)} {cmd:ub(}{it:ub}{cmd:)}
{cmd:maxnodes(}{it:#}{cmd:)} {cmd:inttol(}{it:#}{cmd:)}
{cmd:lptol(}{it:#}{cmd:)} {cmd:lpengine(ipm|simplex|auto)}
{cmd:lpscaling(none|row)} {cmd:trace}]

{title:Description}

{pstd}{cmd:milp_bb} optimizes c*x subject to Aeq*x=beq, Aub*x<=bub,
and lb<=x<=ub.  {cmd:integer()} contains one-based column indices.  Variables
not listed there remain continuous.  Defaults are maximization, lower bounds
zero, no finite upper bounds, 10,000 nodes, integrality tolerance 1e-7, and
LP tolerance 1e-10, {cmd:lpengine(ipm)}, and {cmd:lpscaling(none)}. Stata 16
or newer is required.

{pstd}The implementation uses deterministic depth-first branch-and-bound. Node
relaxations pass through the shared internal LP contract. {cmd:lpengine(ipm)}
uses Mata's {cmd:LinearProgram()} class;
{cmd:lpengine(simplex)} uses the package's revised-simplex core; and
{cmd:lpengine(auto)} tries simplex first and checks unsuccessful results with
IPM unless a structural proof is available. Native retries use the bounded
policy documented in {help lp}. It is intended for teaching, reproducibility, and small or
medium experiments; it is not a replacement for a production MILP solver.

{pstd}{cmd:lpscaling(row)} uses power-of-two multipliers
for each equality and inequality row and for the objective. Each transformation
must exactly recover every stored coefficient on reversal; otherwise its multiplier
is one. Nonzero maxima normally fall in [1,2).
Compared with earlier max-absolute scaling, this optional policy can change
previously accepted numerical results. The reported
solution, objective, bounds, and feasibility diagnostics remain in the
original units. This option controls the initial representation. After failure,
bounded native retries may also use scaling, including with {cmd:lpscaling(none)}.

{pstd}
Simplex node relaxations use the common LP recovery policy documented in
{help lp}. Repeated Phase-2 states can trigger a bounded price retry after the
existing restarts. Original LP checks and the MILP integrality test still apply.

{pstd}
The common simplex policy also permits original-column reconstruction after
all earlier paths fail at a blocked Phase-1 pivot. Node bounds and original
constraints remain in force, followed by the usual integrality checks.

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
With {cmd:lpengine(auto)}, a new simplex recovery can replace a previously
successful native fallback. Scores, reference weights, iteration totals and the
reported backend can therefore change across versions. Earlier-result priority
within the common simplex solver does not promise identical automatic-engine
results across versions.

{pstd}
Duplicate reduction applies only to a continuous node relaxation. Restored
original variables undergo the usual integrality, branching, bounds and node
limit checks. It does not aggregate integer variables.

{title:Stored results}

{synoptset 22 tabbed}{...}
{synopt:{cmd:r(x)}}incumbent solution row vector{p_end}
{synopt:{cmd:r(objective)}}incumbent objective{p_end}
{synopt:{cmd:r(status)}}0 optimal; 1 infeasible; 2 unbounded relaxation;
3 node limit with incumbent; 4 node limit without incumbent; 5 error{p_end}
{synopt:{cmd:r(status_text)}}status description{p_end}
{synopt:{cmd:r(nodes)}}processed nodes{p_end}
{synopt:{cmd:r(lp_solves)}}LP relaxations solved{p_end}
{synopt:{cmd:r(simplex_solves)}}simplex adapter calls{p_end}
{synopt:{cmd:r(ipm_solves)}}IPM adapter calls{p_end}
{synopt:{cmd:r(lp_fallbacks)}}simplex-to-IPM fallbacks{p_end}
{synopt:{cmd:r(native_retries)}}extra native solves across all nodes{p_end}
{synopt:{cmd:r(status_conflicts)}}nodes with contradictory solver reports{p_end}
{synopt:{cmd:r(objective_scale)}}objective divisor used by row scaling{p_end}
{synopt:{cmd:r(row_scale_min)}}smallest constraint-row multiplier{p_end}
{synopt:{cmd:r(row_scale_max)}}largest constraint-row multiplier{p_end}
{synopt:{cmd:r(rootbound)}}root LP objective{p_end}
{synopt:{cmd:r(hasinc)}}whether an incumbent exists{p_end}
{synopt:{cmd:r(integrality_gap)}}zero after an exhaustive optimal search,
missing otherwise{p_end}
{synopt:{cmd:r(cmd)}}{cmd:milp_bb}{p_end}
{synopt:{cmd:r(sense)}}requested objective direction{p_end}
{synopt:{cmd:r(engine)}}relaxation-engine family{p_end}
{synopt:{cmd:r(lpengine)}}requested LP adapter policy{p_end}
{synopt:{cmd:r(lpscaling)}}requested LP scaling policy{p_end}
{synopt:{cmd:r(integer)}}declared integer-variable columns{p_end}
{synopt:{cmd:r(maxnodes)}}configured node limit{p_end}
{synopt:{cmd:r(inttol)}}integrality tolerance{p_end}
{synopt:{cmd:r(lptol)}}LP convergence and residual acceptance tolerance{p_end}

{pstd}Input coefficients, constraint RHS values and lower bounds must be finite.
Invalid matrix dimensions or values return status 5 before solving nodes.
Upper bounds may be missing ({cmd:.}) for no upper limit. Lower bounds default
 to zero and may be negative finite values. Tolerances must be positive and
finite; the node limit must be a positive finite integer. Row or column RHS
vectors are accepted. Always check {cmd:r(status)} and {cmd:r(hasinc)}.
Status 2 describes the LP relaxation and does not certify integer unboundedness.
The objective is computed from the returned incumbent after integer rounding.
From version 1.2.3, both node LP solutions and rounded incumbents must satisfy
the scaled primal-residual criterion described in {help lp}, using {cmd:lptol()}.
If a near-integer candidate fails after rounding, branching continues on a
fractional integer variable. This prevents an infeasible rounded incumbent;
it does not supply an independent optimality or infeasibility certificate.

{pstd}For the historical coefficient-table interface and integer variables
specified by variable name, use {help lp}. These independent optimization tools
are maintained package extensions; their inclusion does not extend the journal
article's validation claims to every optimization problem.

{pstd}From version 1.2.3, the shared LP adapter applies bounded native retries
and exact sign/zero structural checks; see {help lp}. Conflicting reports stop
the search with status 5, rather than pruning a disputed node as infeasible.
{cmd:r(ipm_solves)} counts adapter calls; {cmd:r(native_retries)} separately
counts additional native solves. The known Stata/MP 16.1 misclassification of
max x, x>=0, -x<=0 is resolved by an exact LP ray check; MILP status 2 still
refers only to the continuous relaxation. General status certification and
other native-runtime numerical limitations remain open.

{title:Example}

{pstd}Copy the command lines below (omit leading dots) from {cmd:preserve}
through {cmd:restore} into the Do-file Editor and run them as one block.
Examples use synthetic
data or Stata's bundled {cmd:auto} data and require no downloaded datasets.
The block restores the data in memory. Run dependent lines in order.{p_end}
{* DEAS_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}



{phang2}{cmd:. matrix C=(3,2)}{p_end}
{phang2}{cmd:. matrix A=(2,1 \ 1,2)}{p_end}
{phang2}{cmd:. matrix b=(4 \ 4)}{p_end}
{phang2}{cmd:. milp_bb, c(C) aub(A) bub(b) integer(1 2) sense(max)}{p_end}
{phang2}{cmd:. milp_bb, c(C) aub(A) bub(b) integer(1 2) sense(max) lpengine(auto)}{p_end}
{phang2}{cmd:. assert r(status)==0 & abs(r(objective)-6)<1e-6}{p_end}
{phang2}{cmd:. matrix list r(x)}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}
