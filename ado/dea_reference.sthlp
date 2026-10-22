{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.14.0  20jul2026}{...}
{vieweralsosee "dea" "help dea"}{...}
{viewerjumpto "Syntax" "dea_reference##syntax"}{...}
{viewerjumpto "Description" "dea_reference##description"}{...}
{viewerjumpto "Options" "dea_reference##options"}{...}
{viewerjumpto "Stored results" "dea_reference##results"}{...}
{title:Title}

{p2colset 5 24 26 2}{...}
{p2col:{cmd:dea_reference} {hline 2}}Classify peer weights over the two-stage radial DEA optimal set{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:dea_reference} {ifin}{cmd:,}
{cmd:inputs(}{it:varlist}{cmd:)} {cmd:outputs(}{it:varlist}{cmd:)}
[{cmd:rts(crs|vrs|nirs|ndrs|drs|irs)} {cmd:orientation(input|output)}
{cmd:dmu(}{it:varname}{cmd:)} {cmd:evaluate(}{it:indicator}{cmd:)}
{cmd:candidates(}{it:indicator}{cmd:)}
{cmd:lexicographic(min|max)} {cmd:lexorder(}{it:varname}{cmd:)}
{cmd:policyweight(}{it:varname}{cmd:)} {cmd:policymatrix(}{it:matname}{cmd:)}
{cmd:policysense(min|max)}
{cmd:sparse} {cmd:broad} {cmd:relativeinterior} {cmd:basic} {cmd:audit} {cmd:rtsonly} {cmd:compare} {cmd:interceptbounds}
{cmd:maxnodes(}{it:#}{cmd:)} {cmd:maxauxlp(}{it:#}{cmd:)} {cmd:inttol(}{it:#}{cmd:)}
{cmd:engine(simplex|ipm|auto)} {cmd:scaling(none|row)}
{cmd:tolerance(}{it:#}{cmd:)} {cmd:peertolerance(}{it:#}{cmd:)}
{cmd:slacknormalization(raw|range)} {cmd:radialtolerance(}{it:#}{cmd:)}
{cmd:opttolerance(}{it:#}{cmd:)} {cmd:generate(}{it:prefix}{cmd:)}
{cmd:saving(}{it:filename}{cmd:)} {cmd:replace} {cmd:trace}]

{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_reference} first solves the radial technical-efficiency model with
stage 2.  For every evaluated DMU and candidate peer, it then minimizes and
maximizes that peer's intensity while preserving the radial optimum, the
maximum total slack within {cmd:opttolerance()}, the technology, and the
stage-2 feasibility equations.

{pstd}
By default, the stage-2 objective divides each input and output slack by the
largest absolute observed value of that variable.  This {cmd:range}
normalization makes the declared projection convention invariant to positive
variable-specific changes of measurement units.  {cmd:raw} retains the
unweighted sum for replication of older results.  If exact radial fixing
fails, the command retries once at the conservative endpoint of the declared
{cmd:radialtolerance()} band; the adjustment is returned and is never hidden.

{pstd}
With {cmd:rts(crs)}, the command also minimizes and maximizes
sum(lambda) over that preserved two-stage optimal set.  The returned interval
implements the alternate-optimum RTS test of Cooper, Seiford, and Tone (2007):
an interval wholly below one implies IRS, an interval wholly above one implies
DRS, and an interval containing one implies that a CRS alternate optimum
exists.  {cmd:rtsonly} performs only these two auxiliary LPs per evaluated DMU
and skips the peer-by-peer intensity bounds.

{pstd}
With {cmd:rtsonly generate(prefix)}, the command creates
{it:prefix}{cmd:_lb}, {it:prefix}{cmd:_ub}, {it:prefix}{cmd:_code},
{it:prefix}{cmd:_status}, {it:prefix}{cmd:_margin}, and
{it:prefix}{cmd:_stable}.  The margin is the distance from the relevant bound
to one; for a CRS interval it is the smaller distance from one to either
endpoint.  Stable equals one when this margin exceeds the DMU-specific local
optimality tolerance.  It diagnoses numerical separation under the declared
model and tolerance, not sampling or statistical stability.

{pstd}
With {cmd:rtsonly compare}, the command additionally estimates standard VRS
and NIRS radial scores on the unchanged technology sample. It compares the
resulting CRS--VRS--NIRS code with the optimal-set robust code and returns the
standard code, robust code, and agreement indicator. With {cmd:generate()},
these use suffixes {cmd:_standard}, {cmd:_robust}, and {cmd:_agree}. A
disagreement is reported for diagnosis; neither code is silently substituted
for the other.

{pstd}
With {cmd:rtsonly interceptbounds}, the command obtains the corresponding
two-stage BCC projection, minimizes and maximizes the Cooper-signed BCC
intercept {it:u0} over its supporting-hyperplane optimal set, and compares the
resulting sign classification with the CCR sum(lambda) interval. An unbounded
endpoint is admissible when the finite endpoint already determines the sign.
Always require {cmd:r(u0_bound_status)==0} and
{cmd:r(u0_rts_agree)==1} before reporting the cross-check.

{pstd}
A peer is {it:possible} when its maximum intensity exceeds
{cmd:peertolerance()}, {it:necessary} when its minimum intensity exceeds that
tolerance, {it:optional} when possible but not necessary, and {it:impossible}
otherwise.  These classifications describe the stated two-stage optimal set.
They are not claims that a generic IPM solution is an analytic center or that
the result is a global reference set for a different DEA model.

{pstd}
With all observations selected, the bounds require 2N-squared auxiliary LPs;
{cmd:evaluate()} and {cmd:candidates()} reduce this to 2*N_e*N_c.
Lexicographic selection adds N_e*N sequential LPs because it must order the
full peer vector. The command is therefore intended for reference-set
diagnosis rather than automatic use in every large DEA run.

{pstd}
With {cmd:lexicographic()}, the command additionally constructs one
reproducible intensity vector for each evaluated DMU. It optimizes every peer
intensity sequentially over the same two-stage optimal set and fixes each
attained value within the documented optimality tolerance before proceeding.
This is an explicit tie-breaking rule, not an attempt to reproduce an
undocumented simplex pivot path. Sequential fixing can accumulate numerical
residuals; always require {cmd:r(lex_failures)==0}. If a well-scaled problem
fails, rerun with a justified {cmd:opttolerance()} and report that value.

{pstd}
With {cmd:sparse}, the command minimizes the number of positive peer weights
over the same optimal set by calling {cmd:milp_bb}. For VRS and DRS the valid
linking bound is lambda_j<=1. For CRS and IRS the command uses each peer's
previously verified finite {cmd:r(lambda_max)}; it never substitutes an
arbitrary big-M. A missing or unbounded linking value produces a reported
failure. Always require {cmd:r(sparse_failures)==0}. Sparse selection uses the
complete peer universe and therefore may not be combined with
{cmd:candidates()}.

{pstd}
With {cmd:broad}, the command averages one admissible maximizing solution for
each possible peer. Convexity keeps this average on the same two-stage optimal
set, and every possible peer receives positive weight without introducing a
user-chosen delta or big-M. This is a reproducible broad-support construction;
its support is solver-independent under the stated tolerance, but its component
weights can vary when a peer-specific maximum has alternative solutions. It is
not claimed to be the analytic center of the full optimal face. Always
require {cmd:r(broad_failures)==0}. The construction uses the complete peer
universe and may not be combined with {cmd:candidates()}.

{pstd}
With {cmd:relativeinterior}, the command first finds a feasible point in the
tolerance-defined two-stage optimal polyhedron. It places an auxiliary finite
box one unit above that point, maximizes every nonnegative coordinate, and
maximizes every inequality slack. Averaging the base point and all positive
witnesses gives a point that is strict for every nonfixed lower bound and
nonimplicit inequality. Convexity keeps the average in the original optimal
polyhedron; the box is only a device for obtaining finite witnesses and does
not enter the returned model. This constructs a relative-interior witness, not
the analytic center. Require {cmd:r(relative_failures)==0}. The option uses the
complete peer universe and may not be combined with {cmd:candidates()}.

{pstd}
With {cmd:basic}, the command applies the marked-sample order weights
w_j=(N-j+1)/N to the peer block and solves that third objective with the
released revised-simplex adapter. The returned point is a basic feasible
optimum of the same tolerance-defined optimal polyhedron. This is simplex basic
recovery, not a crossover that inherits an IPM basis, not minimum-cardinality
selection, and not a guarantee that another solver's pivot rule will return
the same vertex. Require {cmd:r(basic_failures)==0} and report the data order.

{pstd}
With {cmd:policyweight()}, the command uses the variable's value for peer j as
the coefficient on lambda_j in a common third objective over the same optimal
set. With {cmd:policymatrix()}, a 1-by-N row supplies common coefficients or an
N-by-N matrix supplies evaluated-DMU-specific coefficients from row o.
The default {cmd:policysense(min)} selects the lowest weighted policy cost;
{cmd:max} selects the highest weighted policy preference. Values may encode
organizational incompatibility, risk, or an externally defined managerial
score; an N-by-N matrix can encode pair-specific distance. The variable must
be complete on the full marked
technology sample; missing values cause an error rather than silently changing
the frontier. If the policy objective has alternative optima, the returned
vector remains solver-selected. Require {cmd:r(policy_failures)==0}.

{pstd}
{cmd:audit} runs {cmd:basic}, {cmd:sparse}, {cmd:broad}, and
{cmd:relativeinterior} together; policy results are included when policy
coefficients are supplied. It reports support counts, L1 distances from the
original stage-2 selected decomposition, the number of possible peers omitted
by that selected decomposition, and mean/maximum lambda identification-range
width. It does not require component equality across representatives. Because
the audit diagnoses the complete peer universe, it may not be combined with
{cmd:candidates()}. Every component failure count must be zero before the audit
table is interpreted.

{pstd}
Before auxiliary optimization begins, the command reports a planned continuous
LP count. This includes peer bounds and requested lexicographic,
relative-interior, basic, and policy LPs. Broad support reuses peer-maximizing
LPs; data-dependent sparse B&B relaxations are counted separately.
{cmd:maxauxlp()} can reject an oversized continuous plan before it runs.
Automatic engine retries can make the realized continuous count exceed the
plan.

{marker options}{...}
{title:Options}

{phang}{cmd:rts()} specifies returns to scale; default {cmd:vrs}. The preferred
one-sided names are {cmd:nirs} and {cmd:ndrs}; {cmd:drs} and {cmd:irs} remain
backward-compatible aliases.{p_end}
{phang}{cmd:rtsonly} requires {cmd:rts(crs)} and computes only the robust
min/max sum(lambda) RTS diagnostic. It may be combined with {cmd:evaluate()}
but not with peer-selection, {cmd:candidates()}, or {cmd:saving()} options.{p_end}
{phang}{cmd:generate()} creates observation-level robust RTS diagnostic
variables and currently requires {cmd:rtsonly}. Existing variables are never
overwritten.{p_end}
{phang}{cmd:compare} requires {cmd:rtsonly} and compares the ordinary
CRS--VRS--NIRS score rule with the CCR optimal-set bound rule. It adds two
stage-1 LPs per technology DMU.{p_end}
{phang}{cmd:orientation()} specifies input or output orientation; default {cmd:input}.{p_end}
{phang}{cmd:engine()} specifies the auxiliary LP engine; default {cmd:auto}.{p_end}
{phang}{cmd:scaling()} specifies canonical row scaling; default {cmd:row}.{p_end}
{phang}{cmd:dmu()} specifies a numeric or string DMU identifier to preserve as
{cmd:dmu_id} and {cmd:peer_id} in the saved long-format file.{p_end}
{phang}{cmd:evaluate()} specifies a numeric nonzero indicator selecting the
DMUs for which bounds are calculated.{p_end}
{phang}{cmd:candidates()} specifies a numeric nonzero indicator selecting the
peer columns to diagnose. Neither indicator changes the frontier technology,
which continues to use the complete {cmd:if}/{cmd:in} estimation sample.{p_end}
{phang}{cmd:lexicographic(min|max)} requests a reproducible representative.
{cmd:max} gives earlier peers priority for larger intensity; {cmd:min} drives
earlier peer intensities downward before later peers are considered.{p_end}
{phang}{cmd:lexorder()} supplies a numeric priority variable. Smaller values
are processed first; equal priorities are resolved by marked-sample order.
Without this option, marked-sample order is used. The option requires
{cmd:lexicographic()}.{p_end}
{phang}{cmd:sparse} requests a minimum-cardinality optimal peer decomposition.
It may be combined with {cmd:evaluate()} but not {cmd:candidates()}.{p_end}
{phang}{cmd:broad} requests an equal average of the peer-specific maximizing
solutions, giving positive weight to every possible peer. It may be combined
with {cmd:evaluate()} but not {cmd:candidates()}.{p_end}
{phang}{cmd:relativeinterior} requests an LP witness average that is positive
for every nonfixed nonnegative coordinate and strict for every nonimplicit
inequality of the tolerance-defined optimal polyhedron. It may be combined
with {cmd:evaluate()} but not {cmd:candidates()}.{p_end}
{phang}{cmd:basic} requests a revised-simplex basic solution selected by the
documented descending marked-sample order weights. It may be combined with
{cmd:evaluate()} and {cmd:candidates()}; neither changes the full third
objective or frontier technology.{p_end}
{phang}{cmd:policyweight()} supplies a numeric peer-level coefficient for an
explicit policy third objective. It must be complete on the marked technology
sample.{p_end}
{phang}{cmd:policymatrix()} supplies either a 1-by-N common coefficient row or
an N-by-N matrix whose row o contains the coefficients for evaluated DMU o.
It may not be combined with {cmd:policyweight()}.{p_end}
{phang}{cmd:policysense(min|max)} chooses cost minimization or preference
maximization; default {cmd:min}. It requires {cmd:policyweight()} or
{cmd:policymatrix()}.{p_end}
{phang}{cmd:audit} requests the integrated basic/sparse/broad/relative-interior
comparison and automatically enables those four selection policies. It may not
be combined with {cmd:candidates()}.{p_end}
{phang}{cmd:maxnodes()} sets the B&B node limit for each sparse MILP; default
{cmd:10000}.{p_end}
{phang}{cmd:maxauxlp()} sets the maximum planned number of continuous auxiliary
LP calls. The default {cmd:-1} is unlimited. Sparse B&B relaxations remain
governed by {cmd:maxnodes()}; zero is invalid.{p_end}
{phang}{cmd:inttol()} sets the binary integrality tolerance; default
{cmd:1e-7}.{p_end}
{phang}{cmd:tolerance()} specifies the LP tolerance; default {cmd:1e-8}.{p_end}
{phang}{cmd:peertolerance()} defines positive peer support; default {cmd:1e-7}.{p_end}
{phang}{cmd:slacknormalization()} selects {cmd:range} (default) or the legacy
unweighted {cmd:raw} stage-2 slack objective.{p_end}
{phang}{cmd:radialtolerance()} sets the relative, scale-aware retry band for
radial fixing. The default is 100*tolerance; zero disables repair.{p_end}
{phang}{cmd:opttolerance()} defines the scaled band around the stage-2 optimum.
The default base is max(100*tolerance, peertolerance), multiplied for each DMU
by 1+abs(total slack).{p_end}
{phang}{cmd:saving()} saves one row per evaluated-DMU--candidate-peer pair with
indices, intensity bounds, numeric status, and status text.  Specify
{cmd:replace} to overwrite an existing file.{p_end}

{marker results}{...}
{title:Stored results}

{p2colset 9 30 32 2}{...}
{p2col:{cmd:r(lambda_selected)}}stage-2 selected intensity matrix{p_end}
{p2col:{cmd:r(lambda_min)}}minimum admissible intensities{p_end}
{p2col:{cmd:r(lambda_max)}}maximum admissible intensities{p_end}
{p2col:{cmd:r(rts_sum_bounds)}}minimum and maximum sum(lambda) on the preserved CCR optimal set{p_end}
{p2col:{cmd:r(rts_robust_code)}}robust RTS code: -1 DRS, 0 CRS, 1 IRS{p_end}
{p2col:{cmd:r(rts_bound_status)}}0 when both RTS-bound LPs succeeded{p_end}
{p2col:{cmd:r(rts_margin)}}distance from the classification boundary at one{p_end}
{p2col:{cmd:r(rts_stable)}}1 when the margin exceeds local optimality tolerance{p_end}
{p2col:{cmd:r(rts_compare)}}standard code, robust code, and agreement indicator{p_end}
{p2col:{cmd:r(n_rts_disagree)}}number of evaluated DMUs with disagreeing codes{p_end}
{p2col:{cmd:r(rts_compare_lp_planned)}}VRS/NIRS comparison LPs planned{p_end}
{p2col:{cmd:r(u0_bounds)}}minimum and maximum Cooper-signed BCC intercept{p_end}
{p2col:{cmd:r(u0_rts_code)}}RTS code implied by the intercept interval{p_end}
{p2col:{cmd:r(u0_bound_status)}}intercept-bound solver status{p_end}
{p2col:{cmd:r(u0_rts_agree)}}agreement with the CCR intensity-interval code{p_end}
{p2col:{cmd:r(radial_adjustment)}}stage-2 radial value minus the reported stage-1 value{p_end}
{p2col:{cmd:r(radial_tolerance)}}declared radial retry tolerance{p_end}
{p2col:{cmd:r(slack_normalization)}}stage-2 slack normalization rule{p_end}
{p2col:{cmd:r(lambda_lexicographic)}}lexicographically selected full peer vectors{p_end}
{p2col:{cmd:r(lambda_sparse)}}minimum-cardinality full peer vectors{p_end}
{p2col:{cmd:r(sparse_support)}}number of positive weights in each sparse vector{p_end}
{p2col:{cmd:r(sparse_status)}}MILP status for each evaluated DMU{p_end}
{p2col:{cmd:r(lambda_broad)}}broad-support full peer vectors{p_end}
{p2col:{cmd:r(broad_support)}}number of possible peers represented{p_end}
{p2col:{cmd:r(broad_status)}}0 success, 1 no possible peer, 5 auxiliary failure{p_end}
{p2col:{cmd:r(lambda_relative)}}relative-interior peer vectors{p_end}
{p2col:{cmd:r(relative_solution)}}full lambda/input-slack/output-slack witnesses{p_end}
{p2col:{cmd:r(relative_status)}}0 success, 5 auxiliary failure{p_end}
{p2col:{cmd:r(relative_witnesses)}}number of averaged witnesses by DMU{p_end}
{p2col:{cmd:r(lambda_basic)}}simplex-basic peer vectors{p_end}
{p2col:{cmd:r(basic_solution)}}full lambda/input-slack/output-slack basic solutions{p_end}
{p2col:{cmd:r(basic_status)}}simplex status for each evaluated DMU{p_end}
{p2col:{cmd:r(lambda_policy)}}policy-selected peer vectors{p_end}
{p2col:{cmd:r(policy_solution)}}full lambda/input-slack/output-slack policy solutions{p_end}
{p2col:{cmd:r(policy_status)}}policy LP status for each evaluated DMU{p_end}
{p2col:{cmd:r(policy_objective)}}attained policy objective by DMU{p_end}
{p2col:{cmd:r(audit_support)}}selected/basic/sparse/broad/relative/policy support counts{p_end}
{p2col:{cmd:r(audit_l1_selected)}}L1 distance of each alternative from selected lambda{p_end}
{p2col:{cmd:r(audit_identification)}}possible count, omitted count, mean and maximum range{p_end}
{p2col:{cmd:r(possible)}}possible-peer indicator matrix{p_end}
{p2col:{cmd:r(necessary)}}necessary-peer indicator matrix{p_end}
{p2col:{cmd:r(reference_status)}}0 impossible, 1 optional, 2 necessary{p_end}
{p2col:{cmd:r(aux_status)}}auxiliary min/max solver status pairs{p_end}
{p2col:{cmd:r(aux_failures)}}number of failed peer-bound pairs{p_end}
{p2col:{cmd:r(lex_failures)}}number of failed lexicographic selections{p_end}
{p2col:{cmd:r(sparse_failures)}}number of failed sparse selections{p_end}
{p2col:{cmd:r(sparse_nodes)}}total B&B nodes over sparse selections{p_end}
{p2col:{cmd:r(lp_solves_planned)}}continuous auxiliary LP calls planned before execution{p_end}
{p2col:{cmd:r(aux_lp_solves)}}continuous auxiliary LP calls actually dispatched{p_end}
{p2col:{cmd:r(sparse_lp_solves)}}LP relaxations solved inside sparse B&B{p_end}
{p2col:{cmd:r(total_lp_solves)}}continuous plus sparse LP calls{p_end}
{p2col:{cmd:r(maxauxlp)}}requested continuous-LP planning limit{p_end}
{p2col:{cmd:r(broad_failures)}}number of failed broad selections{p_end}
{p2col:{cmd:r(relative_failures)}}number of failed relative-interior selections{p_end}
{p2col:{cmd:r(basic_failures)}}number of failed simplex-basic selections{p_end}
{p2col:{cmd:r(policy_failures)}}number of failed policy selections{p_end}
{p2col:{cmd:r(evaluated_index)}}evaluated rows in the marked technology sample{p_end}
{p2col:{cmd:r(candidate_index)}}candidate rows in the marked technology sample{p_end}
{p2col:{cmd:r(N_evaluated)}}number of evaluated DMUs{p_end}
{p2col:{cmd:r(N_candidates)}}number of diagnosed peer candidates{p_end}
{p2colreset}{...}

{title:Example}

{pstd}Copy the command lines below (omit leading dots) from {cmd:preserve}
through {cmd:restore} into the Do-file Editor and run them as one block.
Examples use synthetic
data or Stata's bundled {cmd:auto} data and require no downloaded datasets.
The block restores the data in memory. Run dependent lines in order.{p_end}
{* DEAS_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input id x1 x2 y1 y2 region}{p_end}
{phang2}{cmd:. 1 1 2 1 1 1}{p_end}
{phang2}{cmd:. 2 2 4 1 1 1}{p_end}
{phang2}{cmd:. 3 3 6 1 1 2}{p_end}
{phang2}{cmd:. 4 4 8 1 1 2}{p_end}
{phang2}{cmd:. end}{p_end}
{phang2}{cmd:. generate priority = id}{p_end}
{phang2}{cmd:. generate distance = id}{p_end}
{phang2}{cmd:. generate cost = id}{p_end}
{phang2}{cmd:. generate byte review = id == 2}{p_end}
{phang2}{cmd:. matrix D = (1,2,3,4)}{p_end}
{phang2}{cmd:. tempfile peer_bounds peer_audit}{p_end}


{phang2}{cmd:. dea_reference, inputs(x1 x2) outputs(y1) rts(crs) rtsonly compare generate(rtsdiag)}{p_end}
{phang2}{cmd:. matrix list r(rts_sum_bounds)}{p_end}
{phang2}{cmd:. matrix list r(rts_robust_code)}{p_end}
{phang2}{cmd:. generate byte audit = region==1}{p_end}
{phang2}{cmd:. dea_reference, inputs(x1 x2) outputs(y1 y2) dmu(id) evaluate(audit) rts(vrs) saving("`peer_bounds'") replace}{p_end}
{phang2}{cmd:. dea_reference, inputs(x1 x2) outputs(y1 y2) dmu(id) evaluate(audit) rts(vrs) lexicographic(max) lexorder(priority) saving("`peer_bounds'") replace}{p_end}
{phang2}{cmd:. dea_reference, inputs(x1 x2) outputs(y1 y2) dmu(id) evaluate(audit) rts(vrs) sparse saving("`peer_bounds'") replace}{p_end}
{phang2}{cmd:. dea_reference, inputs(x1 x2) outputs(y1 y2) dmu(id) evaluate(audit) rts(vrs) broad saving("`peer_bounds'") replace}{p_end}
{phang2}{cmd:. dea_reference, inputs(x1 x2) outputs(y1 y2) dmu(id) evaluate(audit) rts(vrs) relativeinterior saving("`peer_bounds'") replace}{p_end}
{phang2}{cmd:. dea_reference, inputs(x1 x2) outputs(y1 y2) dmu(id) evaluate(audit) rts(vrs) basic saving("`peer_bounds'") replace}{p_end}
{phang2}{cmd:. dea_reference, inputs(x1 x2) outputs(y1 y2) dmu(id) evaluate(audit) rts(vrs) policyweight(distance) policysense(min) saving("`peer_bounds'") replace}{p_end}
{phang2}{cmd:. dea_reference, inputs(x1 x2) outputs(y1 y2) dmu(id) evaluate(audit) rts(vrs) policymatrix(D) policysense(min) saving("`peer_bounds'") replace}{p_end}
{phang2}{cmd:. dea_reference, inputs(x1 x2) outputs(y1 y2) dmu(id) evaluate(review) rts(vrs) audit policyweight(cost) saving("`peer_audit'") replace}{p_end}
{phang2}{cmd:. matrix list r(lambda_max)}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

