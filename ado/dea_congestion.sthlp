{smcl}
{* Help revision 2026-11-06: certified simplex engine, shared weight guard and explicit reference selection.}{...}
{* Help revision 2026-11-02: finite-score counts, component failures and saved decomposition.}{...}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* *! version 2.3.0 06nov2026}{...}
{viewerjumpto "Related article" "dea_congestion##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_sbmundesirable" "help dea_sbmundesirable"}{...}
{viewerjumpto "Syntax" "dea_congestion##syntax"}{...}
{viewerjumpto "Description" "dea_congestion##description"}{...}
{viewerjumpto "Options" "dea_congestion##options"}{...}
{viewerjumpto "Remarks" "dea_congestion##remarks"}{...}
{viewerjumpto "Examples" "dea_congestion##examples"}{...}
{viewerjumpto "Stored results" "dea_congestion##stored_results"}{...}
{viewerjumpto "References" "dea_congestion##references"}{...}
{viewerjumpto "Authors" "dea_congestion##authors"}{...}
{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{cmd:dea_congestion} {hline 2}}Cooper congestion decomposition and legacy radial diagnostic{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea_congestion}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}
[{cmd:,} {it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_congestion}
{ifin}{cmd:,}
{cmdab:in:puts(}{it:varlist}{cmd:)}
{cmdab:out:puts(}{it:varlist}{cmd:)}
[{it:options}]


{pstd}
Both forms support {cmd:if} and {cmd:in}. Equality tests such as
{cmd:if group==1}, quoted strings, and functions such as
{cmd:if inlist(group,1,2)} are parsed as sample restrictions.
The restrictions apply to both evaluated observations and the reference data.
Input and output varlists may use standard variable ranges and wildcards.

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmdab:in:puts(}{it:varlist}{cmd:)}}input variable(s){p_end}
{synopt:{cmdab:out:puts(}{it:varlist}{cmd:)}}output variable(s){p_end}

{syntab:Method}
{synopt:{cmdab:congm:ethod(}{it:string}{cmd:)}}{cmd:legacy} (default) or {cmd:cooper}{p_end}
{synopt:{cmdab:congw:eight(}{it:string}{cmd:)}}Cooper slack objective: {cmd:raw} (default) or {cmd:normalized}{p_end}
{synopt:{cmdab:congmat:rix}}return Cooper component, projection, peer and status matrices{p_end}
{synopt:{cmdab:congr:anges}}also compute optimal-face ranges; implies {cmd:congmatrix}{p_end}
{synopt:{cmdab:congp:olicy(}{it:stage1|stage2}{cmd:)}}simplex reference-selection policy; default {cmd:stage1}{p_end}
{synopt:{cmdab:conga:nchor(}{it:A}{cmd:)}}stage-1 reference lambda matrix, aligned with the selected sample{p_end}
{synopt:{cmdab:congref:erence(}{it:M}{cmd:)}}stage-2 reference mu matrix; required with {cmd:congpolicy(stage2)}{p_end}
{synopt:{cmdab:congs:electlimit(}{it:16|128}{cmd:)}}stage2-selector basis-exchange limit; default 16, explicit extension 128{p_end}


{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for generated variables; default is {cmd:_cong}{p_end}
{synopt:{cmdab:sav:ing(}{it:filename}{cmd:)}}save selected observations and generated results{p_end}
{synopt:{cmd:replace}}allow replacement of the saved file{p_end}
{synopt:{cmd:debug}}show per-observation Cooper status and adjustment counts{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}

{syntab:Technical}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable{p_end}
{synopt:{cmd:engine(}{it:string}{cmd:)}}solver engine; default is {cmd:engine(ipm)}{p_end}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:congmethod(cooper)} implements the two-stage hybrid congestion procedure
in Cooper, Deng, Gu, Li, and Thrall's author manuscript (revised November 2000),
equations (1)--(7). It first expands outputs and maximizes remaining slacks,
then separates input slacks into technical inefficiency and congestion.
The default {cmd:congmethod(legacy)} retains the historical diagnostic described
below. The two methods have different output meanings and suffixes.
See {help dea_congestion##cooper:Cooper method} for the full specification.

{pstd}{bf:Legacy method (default)}{p_end}

{pstd}
{cmd:dea_congestion} retains the package's historical radial ratio for
compatibility. It compares the standard input-oriented VRS score theta_VRS
with theta_C from the more restrictive LP
min theta subject to X'lambda=theta*x0, Y'lambda=y0, sum(lambda)=1,
lambda>=0 and 0<=theta<=1. CRS efficiency is also computed for scale efficiency.

{pstd}
The ratio theta_VRS/theta_C measures the gap caused by imposing both sets
of equalities. It is {bf:not a sufficient or necessary test of input congestion}.
It is not the Cooper et al. slack-decomposition procedure and is not the
usual comparison that changes only input disposability while retaining
output inequalities. The cited congestion literature is background.

{pstd}
A simple counterexample has equal output 1 and input vectors (1,1) and (2,3).
The second DMU has theta_VRS=.5, theta_C=1 and ratio=.5, although output
cannot increase anywhere in the observed convex hull. Conversely, inputs
(1,1), (2,1) with outputs 2,1 allow an input reduction and output increase,
but the second DMU's ratio is 1. The ratio does not identify which input
is congested or its amount; it must not be interpreted as a causal effect.

{pstd}
The primary variable named by {cmd:generate()} stores the legacy ratio.
Five auxiliary variables are also created with that prefix:

{p2colset 9 32 34 2}
{p2col:{it:prefix}{cmd:_crs}}CRS technical efficiency{p_end}
{p2col:{it:prefix}{cmd:_vrs}}VRS technical efficiency{p_end}
{p2col:{it:prefix}{cmd:_cong_eff}}congestion-model efficiency ({it:theta_C}){p_end}
{p2col:{it:prefix}{cmd:_congestion}}legacy ratio ({it:theta_VRS / theta_C}){p_end}
{p2col:{it:prefix}{cmd:_se}}scale efficiency ({it:theta_CRS / theta_VRS}){p_end}
{p2colreset}{...}

{pstd}
{cmd:dea_congestion} is a convenience wrapper around
{cmd:dea ..., model(congestion)}.


{pstd}
From core {cmd:dea} 2.0.53, an unsuccessful congestion LP leaves
{it:prefix}{cmd:_cong_eff} missing. It is never replaced by 1.
The primary measure and {cmd:_congestion} require finite VRS and congestion
scores and a congestion denominator greater than {cmd:tolerance()}.
{cmd:_se} likewise requires finite CRS/VRS scores and a VRS denominator
greater than the tolerance. Undefined ratios remain missing.
Available components are retained when another component is unavailable.

{pstd}
{cmd:r(n_failed)} counts each selected row once if any requested score or ratio
is unavailable. {cmd:r(n_failed_primary)} counts missing congestion measures;
a CRS-only failure can make {cmd:r(n_failed)}>0 while the main measure remains
available. Means use nonmissing scores; efficiency counts require a finite
main score within {cmd:tolerance()} of 1. Missing values do not indicate no congestion.

{pstd}
Use the core command with {cmd:saving()} to save the ID, inputs, outputs,
primary measure and all five auxiliary result columns. The wrapper does not
accept {cmd:saving()}. Rows with unavailable results remain in the saved sample.
Both entry points use the specialized IPM routines; accepting {cmd:engine()}
does not change that computational backend.

{marker options}{...}
{title:Options}

{dlgtab:Output}

{phang}
{cmdab:gen:erate(}{it:prefix}{cmd:)} specifies the prefix for all generated
variables.  Default is {cmd:_cong}.

{dlgtab:Technical}

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies the variable identifying each DMU.
Accepts both string and numeric variables.

{phang}
{cmd:engine(}{it:string}{cmd:)} specifies the solver engine.
{cmd:ipm}, {cmd:simplex} and {cmd:auto} are accepted by the shared parser,
but this specialized calculation always uses Stata {cmd:LinearProgram()}.
The selected engine name is retained in {cmd:r(engine)}; the actual backend
is identified by {cmd:r(lp_backend)}.

{phang}
{cmdab:tol:erance(}{it:#}{cmd:)} sets the numerical tolerance for the LP solver.
Default is {cmd:1e-7}.


{marker remarks}{...}
{title:Remarks}

{pstd}
Congestion arises when continued application of an input reduces rather
than increases output.  This can occur in situations of resource overuse:
too many workers on a fixed production line, excess fertilizer on a crop,
or overcrowding in a facility.

{pstd}
Selected quantities must be nonnegative, with at least one positive input
and output in each DMU. Individual zeros and redundant zero columns are allowed.
The equality LP scales each data column by its positive maximum, preserving
theta and lambda while reducing sensitivity to measurement units.
Self-reference provides a feasible theta=1 solution; a missing numerical
result is not a mathematical infeasibility certificate.

{pstd}
The core interface uses input orientation for this diagnostic regardless
of {cmd:orientation()} and computes its own component LPs regardless of
{cmd:stage()}. The core normalizes {cmd:rts()} to VRS for the primary ratio;
the auxiliary CRS efficiency is calculated separately. The wrapper does not
accept these model options. {cmd:r(orientation)} reports the effective input
orientation; {cmd:r(orientation_requested)} records the parsed request.

{pstd}
The model requires Stata 16 or later (for the {cmd:LinearProgram()}
class used in the interior-point solver).


{pstd}
{bf:Protecting existing data.} All output names are checked before results are
created. An existing result variable causes an error unless {cmd:replacevars}
is specified. The option applies only to variables produced by the current
command; it does not remove other variables sharing the prefix. Replacement
recreates whole result columns, including their storage types and labels;
observations not scored by the command do not retain old scores. Variables
used as model data or identifiers cannot be replaced, even with
{cmd:replacevars}. Invalid or duplicate output names are rejected before any
existing result column is changed. This name check does not promise rollback
of later computation or file-write errors.

{marker cooper}{...}
{title:Cooper method}

{pstd}
Use {cmd:congmethod(cooper)} with either wrapper syntax, or with
{cmd:dea ..., model(congestion) congmethod(cooper)}. The estimation sample is
also the reference sample, including each evaluated DMU itself. All selected
inputs and outputs must be strictly positive. Missing quantities or IDs are
excluded; one observation is allowed. With no {cmd:dmu()}, a new {cmd:_dmu}
contains the original observation numbers; an existing {cmd:_dmu} requires
explicit {cmd:dmu(_dmu)} or another identifier. Existing result variables require
{cmd:replacevars}; inputs, outputs and IDs are protected.

{pstd}
There are three sequential LPs (two substantive stages). First maximize phi
subject to X'lambda+sminus=x0, Y'lambda-splus=phi*y0,
sum(lambda)=1, lambda,sminus,splus>=0 and phi>=1. Fix the attained phi and
maximize the input/output slack objective. Let xhat=x0-sminus and
yhat=phi*y0+splus denote the selected projection. Finally maximize the
technical input slack delta subject to X'mu-delta=xhat, Y'mu=yhat,
sum(mu)=1, mu>=0 and 0<=delta<=sminus. Input congestion is
c=sminus-delta. This is an empirical frontier decomposition, not a causal estimate.
The manuscript's later managerial-efficiency extension is not implemented here.

{pstd}
{cmd:congweight(raw)} uses the unweighted sum of original-unit slacks at each
slack optimization, matching the displayed manuscript objectives.
{cmd:congweight(normalized)} is an explicit variant: phase 1b maximizes
sum(sminus_i/x0_i)+sum(splus_r/y0_r), and phase 2 maximizes
sum(delta_i/x0_i). Raw objectives may change their chosen optimum when units
change. Normalized objectives are invariant to positive column rescaling in
exact arithmetic; neither option guarantees a unique projection or decomposition.
The primary summary is always the mean relative congestion, regardless of
which objective was chosen.

{pstd}
In the November 2000 author manuscript's two 17-year industry examples,
four Table 3 components differ from the displayed VRS/raw-slack procedure:
textile 1989 capital (0 versus 0.52), textile 1993 labor (0 versus 1.7875),
automobile 1990 labor (34.34 versus about 35.08985383), and automobile 1996
labor (0 versus about 1.50372480). These are component quantities, not the
primary relative score. Exact checks of the printed-decimal, stored-binary64
and normalized-binary64 models establish unique projections and target
components for these four cases; alternative optima and these representation
differences do not explain the discrepancies. This is a bounded comparison
with the author manuscript, not a verified erratum in the final publication.
The packaged {cmd:tests/COOPER_TABLE.txt} describes the sources and reproducible
{cmd:cooper_table} regression. No result is adjusted to match the manuscript table.

{p2colset 9 34 36 2}
{p2col:{it:prefix}}mean_i(c_i/x0_i); 0 means no detected input congestion{p_end}
{p2col:{it:prefix}{cmd:_phi}}output expansion phi (at least 1){p_end}
{p2col:{it:prefix}{cmd:_eff}}1/phi; ignores nonradial slacks{p_end}
{p2col:{it:prefix}{cmd:_total}}mean_i(sminus_i/x0_i){p_end}
{p2col:{it:prefix}{cmd:_technical}}mean_i(delta_i/x0_i){p_end}
{p2col:{it:prefix}{cmd:_status}}0 complete; 1 radial LP unavailable; 2 slack LP/projection unavailable; 3 decomposition LP unavailable; 4 exact-point numeric export unavailable (simplex){p_end}
{p2colreset}{...}

{pstd}
With {cmd:engine(ipm)}, a failed phase leaves downstream results missing and
retains accepted upstream components. With {cmd:engine(simplex)}, the complete
point is returned atomically; an unavailable or unrepresentable point leaves its
six-variable numeric result bundle missing apart from the explicit status. Failure is not an infeasibility certificate. No default
value of 1 or 0 replaces an unavailable congestion result. The core's
{cmd:stage(1)} and explicit input orientation are rejected for this method;
VRS and the complete output-oriented procedure are required. FDH, bootstrap,
bad/direction/prices, imprecise data, reference/weight exports, positive
{cmd:radialtolerance()}, {cmd:contractcheck}, and {cmd:trace} do not apply.
Use {cmd:debug} for per-observation diagnostics.

{pstd}
{cmd:engine(simplex)} uses package-supplied revised-simplex proposals followed by
bounded integer primal/dual certification and, when needed, exact basis exchange.
Later stages use original constraint rows and exact optimal-face certificates;
rounded phi or projected targets are not fixed as new right-hand sides. The
certification domain is strictly positive data with
{cmd:2*N + inputs + outputs + 3 <= 64}, at most 256 constraint rows, and bounded
integer work/limb budgets. Larger models may use {cmd:engine(ipm)} explicitly.
An unavailable certificate is not an infeasibility proof. Main-stage exact repair retains16 exchanges. The optional stage2 selector
also defaults to16 and can explicitly use128. These are separate from the
numerical simplex iteration limit.

{pstd}
{cmd:engine(ipm)} retains Stata's native {cmd:LinearProgram()} procedure and
bounded algebraically equivalent native retries, original-constraint acceptance,
bound-roundoff projection and rechecked peer normalization. These are numerical
checks, not exact optimality proofs. Native stages and endpoints may remain
unavailable. The certified simplex path does not call or load the IPM engine.
{cmd:engine(auto)} and {cmd:method(auto)} are explicitly rejected for Cooper;
no implicit engine fallback or IPM reference prepass is performed. Other models'
auto behavior is unchanged. {cmd:r(actual_engine)} and {cmd:r(actual_backend)}
describe the path actually used. Existing defaults are preserved: {cmd:dea}
defaults to simplex; this {cmd:dea_congestion} wrapper defaults to ipm. Specify
the engine and slack objective when comparing results.

{pstd}
Both engines reject column-scaling underflow and objective-weight overflow or
zero underflow before result publication. Respecify numerical units or the
slack objective when this error is reported. Inputs being positive does not
guarantee that every intermediate ratio is representable.

{pstd}
For simplex, {cmd:congpolicy(stage1)} keeps the certified baseline decomposition
vertex; it does not mean the public {cmd:stage(1)} procedure. All three substantive
LP steps still run. Without {cmd:conganchor()}, the chosen exact slack vertex
can differ from a native-IPM projection in a nonunique model. Optional
{cmd:conganchor(A)} minimizes the L-infinity distance to row i of A on the exact
slack-optimal face. {cmd:congpolicy(stage2) congreference(M)} additionally minimizes
the distance to row i of M on the original technical-optimal face. Both reference
matrices must be named N-by-N Stata matrices in the selected observation order;
IDs are not automatically matched or sorted. Copy {cmd:r(lambda_stage1/2)} into
named matrices first. A missing anchor row follows the unanchored path; a missing
mu-reference row skips stage-2 selection. An uncertified selector returns the
certified baseline with an explicit status. {cmd:congselectlimit(128)} is an
explicit extended-work option; the default is never automatically increased.
Selection options require simplex. Existing {cmd:congr} means ranges; the new
reference abbreviation is {cmd:congref} or longer.

{pstd}
Simplex certificates concern internal rational points of the normalized binary64
LP, not rounded returned coordinates or exact rational rescaling of raw data.
Weights and scores are dimensionless. Components, targets and ranges multiply
the exact LP fraction by the original column maximum before one nearest/even
binary64 rounding. The symmetric export range is |value| less than 2^1023.
For simplex, {cmd:r(acceptance_tolerance)} is 0 and {cmd:r(tolerance_role)} is
{cmd:summary_classification_only}; {cmd:tolerance()} affects the efficiency,
congestion and multiplicity summaries, using min(tolerance,1e-9), not the exact
acceptance or basis budget. Native acceptance retains its numerical tolerance.

{pstd}
{cmd:congmatrix} returns original-unit input totals, technical slack and
congestion; output slack; input/output targets; peer weights; and statuses.
Rows follow the selected observations in their existing order, as recorded by
{cmd:r(obs_index)}. Peer-matrix columns follow the same order. Increasing
{cmd:set matsize} may be necessary; peer matrices are returned only when requested. Internal certified calculations
may retain dense arrays within the bounded sample domain. {cmd:saving()} saves the selected IDs, model variables and six result
variables, not these optional matrices. File failures preserve in-memory data.

{pstd}
{cmd:congranges} minimizes and maximizes each input/output target over the
phase-1 optimal face, fixing both phi and the selected weighted slack optimum.
It also minimizes/maximizes each congestion component over the phase-2 optimal
face, {bf:conditional on the selected phase-1 projection}. These second ranges
are not global congestion ranges over all possible phase-1 projections.
Component extrema need not be jointly attainable. No lexicographic tie-break
among individual inputs/outputs is imposed. For simplex, a failed endpoint certificate or conversion makes the whole range
bundle missing, with an explicit range status, while preserving a valid main
point. Native ranges retain their per-endpoint numerical availability. Skipped
ranges remain missing. These are conditional ranges, not global ranges over all
stage-1 projections.

{title:Cooper stored results}

{pstd}
With {cmd:congmethod(cooper)}, {cmd:r(N)} is the selected sample size;
{cmd:r(N_complete)} and {cmd:r(n_failed)} count complete and incomplete
complete point sequences. {cmd:r(n_congested)} counts finite primary summaries
above min(tolerance,1e-9). {cmd:r(n_efficient)} (alias {cmd:r(N_eff)}) counts
solutions with phi near 1 and relative slacks near 0 under the summary threshold.
Native counts may use an accepted upstream phase even if decomposition fails;
simplex counts require a usable complete point. Zero congestion alone does not mean DEA
efficiency. {cmd:r(mean_congestion)} and {cmd:r(mean_phi)} use their respective
available results. {cmd:r(n_endpoint_adjusted)} counts adjusted coordinates
across accepted LP candidates and projection reconstruction, not DMUs.

{pstd}
{cmd:r(formulation)} is {cmd:cooper_two_stage_vrs}; {cmd:r(orientation)} is
{cmd:output}; {cmd:r(rts)} is {cmd:vrs}; {cmd:r(congmethod)} and
{cmd:r(congweight)} identify the specification. Native projection policy is
{cmd:selected_slack_optimum}; simplex reports {cmd:exact_slack_vertex} or
{cmd:nearest_reference_stage1_on_exact_slack_face}. {cmd:r(congpolicy)} identifies
the simplex decomposition-selection policy. No legacy ratio summary aliases are supplied.

{pstd}
With {cmd:congmatrix}: {cmd:r(input_total)}, {cmd:r(input_technical)},
{cmd:r(input_congestion)}, {cmd:r(output_slack)}, {cmd:r(target_inputs)},
{cmd:r(target_outputs)}, {cmd:r(lambda_stage1)}, {cmd:r(lambda_stage2)},
{cmd:r(obs_index)} and {cmd:r(status)}. The seven status columns are radial,
slack, decomposition (0 accepted, 1 unavailable, missing if skipped),
adjusted coordinate count, and maximum absolute phase-1 projection
reconstruction change in scaled coordinates, equivalent-retry count, and
recovered-retry count. The last two count main LPs, excluding optional range LPs.
Their totals are {cmd:r(n_native_retry)} and {cmd:r(n_native_recovered)}.

{pstd}
For the integrated simplex path, {cmd:r(integration_status)} is {cmd:integrated}.
This identifies the command path; per-DMU availability is recorded in
{cmd:r(certificate_status)} and the generated status variable.

{pstd}
For simplex, {cmd:congmatrix} also returns {cmd:r(certificate_status)} with columns:
usable, base_certificate, selection_status, selection_applied,
selection_certificate, conversion_status, rounded_values, selection_maxpivots,
reference_distance, ranges_complete, range_conversion_status,
range_rounded_values, radial, slack, anchor, decomposition, projection, congestion.
Selection status is 0 not requested, 1 certified, 2 reference missing,
3 uncertified, or 4 association/base unavailable. Conversion status is
0 normal, 1 invalid exact representation, 2 resource bound, or 3 export-range
overflow. Stage entries record certificates/counts; anchor 2 means skipped.
Read usable and selection/conversion status together. Legacy seven-column
{cmd:r(status)} remains available; coordinate-adjustment and native-retry counts
are zero for the exact route. Additional scalar counts identify certified,
missing-reference, uncertified and mismatched selections and conversion failures.
Simplex matrix stripes use unique original observation numbers such as obs3;
{cmd:r(obs_index)} maps them to data and IDs. Projection columns use input/output
equation stripes, including when the same variable has both roles.

{pstd}
With {cmd:congranges}: {cmd:r(projection_min/max)} (inputs then outputs),
{cmd:r(congestion_min/max)} and {cmd:r(range_status)} (phase 1, phase 2;
0 complete, 1 at least one endpoint unavailable, missing if skipped).
{cmd:r(n_range_failed)} counts incomplete requested simplex range bundles, or
attempted native range failures;
{cmd:r(n_range_complete)} counts rows with both sets complete.
{cmd:r(n_projection_multiple)} and {cmd:r(n_decomposition_multiple)} count
complete range sets with a width, relative to the observed quantity, greater
than 10 times min(tolerance,1e-9). They are lower-bound diagnostics when
any range endpoints are unavailable. {cmd:r(range_scope)} records conditioning.

{marker examples}{...}
{title:Examples}

{pstd}Copy the command lines below (omit leading dots) from {cmd:preserve}
through {cmd:restore} into the Do-file Editor and run them as one block.
Examples use synthetic
data or Stata's bundled {cmd:auto} data and require no downloaded datasets.
The block restores the data in memory. Run dependent lines in order.{p_end}
{* DEAS_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}


{pstd}Setup

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str10 hospital beds nurses patients}{p_end}
{phang2}{cmd:"Hospital A" 120 55 5100}{p_end}
{phang2}{cmd:"Hospital B" 150 70 6200}{p_end}
{phang2}{cmd:"Hospital C" 180 95 6400}{p_end}
{phang2}{cmd:"Hospital D" 110 50 5000}{p_end}
{phang2}{cmd:"Hospital E" 210 120 6600}{p_end}
{phang2}{cmd:"Hospital F" 160 75 6100}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Historical equality-ratio diagnostic

{phang2}{cmd:. dea_congestion beds nurses = patients, dmu(hospital)}{p_end}

{pstd}Form B syntax

{phang2}{cmd:. dea_congestion, inputs(beds nurses) outputs(patients) dmu(hospital) replacevars}{p_end}

{pstd}Custom prefix for output variables

{phang2}{cmd:. dea_congestion beds nurses = patients, generate(hosp) dmu(hospital)}{p_end}
{phang2}{cmd:. list hospital hosp_vrs hosp_congestion if hosp_congestion < 0.99}{p_end}

{pstd}Equivalently using dea directly

{phang2}{cmd:. dea beds nurses = patients, model(congestion) generate(hosp) dmu(hospital) replacevars}{p_end}



{pstd}Equivalent sample restrictions in both syntax forms

{phang2}{cmd:. dea_congestion beds nurses = patients if inrange(_n,2,5), dmu(hospital) generate(qa)}{p_end}
{phang2}{cmd:. dea_congestion if inrange(_n,2,5), inputs(beds nurses) outputs(patients) dmu(hospital) generate(qb)}{p_end}
{phang2}{cmd:. assert qa==qb}{p_end}


{pstd}Explicit Cooper decomposition (separate result prefix):{p_end}
{phang2}{cmd:. dea_congestion beds nurses = patients, dmu(hospital) congmethod(cooper) congweight(normalized) engine(ipm) congmatrix generate(cp)}{p_end}
{phang2}{cmd:. tempname ca cm}{p_end}
{phang2}{cmd:. matrix `ca'=r(lambda_stage1)}{p_end}
{phang2}{cmd:. matrix `cm'=r(lambda_stage2)}{p_end}
{phang2}{cmd:. list hospital cp cp_phi cp_total cp_technical cp_status}{p_end}
{phang2}{cmd:. matrix list r(input_congestion)}{p_end}
{phang2}{cmd:. dea patients, inputs(beds nurses) model(congestion) congmethod(cooper) congweight(normalized) engine(simplex) conganchor(`ca') congpolicy(stage2) congreference(`cm') congranges dmu(hospital) generate(cpn)}{p_end}
{phang2}{cmd:. matrix list r(congestion_min)}{p_end}
{phang2}{cmd:. matrix list r(congestion_max)}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
With {cmd:congmethod(legacy)}, {cmd:dea_congestion} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:r(N)}}number of DMUs{p_end}
{synopt:{cmd:r(n_congested)}}legacy count of finite ratios below 1-tolerance; not a congestion certificate{p_end}
{synopt:{cmd:r(mean_vrs)}}mean VRS technical efficiency{p_end}
{synopt:{cmd:r(mean_efficiency)}}mean congestion measure{p_end}
{synopt:{cmd:r(n_efficient)}}finite main scores within tolerance of 1; alias of {cmd:r(N_eff)}{p_end}
{synopt:{cmd:r(N_eff)}}finite main scores within tolerance of 1{p_end}
{synopt:{cmd:r(mean_eff)}}mean finite congestion measure{p_end}
{synopt:{cmd:r(n_failed)}}rows with any unavailable requested component or ratio, counted once{p_end}
{synopt:{cmd:r(n_failed_primary)}}missing main congestion measures{p_end}
{synopt:{cmd:r(n_failed_crs)}}missing CRS scores{p_end}
{synopt:{cmd:r(n_failed_vrs)}}missing VRS scores{p_end}
{synopt:{cmd:r(n_failed_cong)}}missing congestion-model scores{p_end}
{synopt:{cmd:r(n_failed_se)}}missing scale-efficiency ratios{p_end}
{p2colreset}{...}


{pstd}
{cmd:r(formulation)} is {cmd:legacy_equalities_ratio} and
{cmd:r(lp_backend)} is {cmd:stata_linearprogram}. Historical variable names
and {cmd:r(n_congested)} are retained for compatibility with existing scripts.

{marker references}{...}
{title:References}

{phang}
Cooper, W. W., H. Deng, B. Gu, S. Li, and R. M. Thrall.
{it:Using DEA to Improve the Management of Congestion in Chinese Industries (1981-1997)}.
Author manuscript, revised November 2000, equations (1)--(7), pp. 4--14.
{browse "https://iiif.library.cmu.edu/file/Cooper_box00010_fld00034_bdl0001_doc0001/Cooper_box00010_fld00034_bdl0001_doc0001.pdf":CMU author manuscript}.


{phang}
Cooper, W. W., H. Deng, L. M. Seiford, and J. Zhu. 2004.
Congestion: Its identification and management with DEA.
In {it:Handbook on Data Envelopment Analysis},
ed. W. W. Cooper, L. M. Seiford, and J. Zhu, 177--201.
Boston: Kluwer Academic Publishers.

{phang}
Cooper, W. W., L. M. Seiford, and K. Tone. 2007.
{it:Data Envelopment Analysis: A Comprehensive Text with Models, Applications,
References and DEA-Solver Software}. 2nd ed. New York: Springer.


{marker authors}{...}
{title:Authors}

{pstd}
Choonjoo Lee (corresponding author){break}
Korea National Defense University{break}
Seoul, Korea{break}
bloom.rampike@gmail.com

{marker alsosee}{...}
{title:Also see}

{pstd}
{it:DEA package commands:}

{synoptset 28 tabbed}{...}
{synopt:{helpb dea}}Data envelopment analysis{p_end}
{synopt:{helpb dea_bootstrap}}Bootstrap bias-corrected DEA{p_end}
{synopt:{helpb dea_congestion}}Congestion DEA{p_end}
{synopt:{helpb dea_cooperative}}Cooperative game DEA{p_end}
{synopt:{helpb dea_crosseff}}Cross-efficiency DEA{p_end}
{synopt:{helpb dea_crossgame}}Cross-efficiency game DEA{p_end}
{synopt:{helpb dea_game}}Game-theoretic DEA overview{p_end}
{synopt:{helpb dea_malmquist}}Malmquist productivity index{p_end}
{synopt:{helpb dea_nash}}Nash bargaining DEA{p_end}
{synopt:{helpb dea_network}}Network DEA{p_end}
{synopt:{helpb dea_sbmundesirable}}SBM with undesirable outputs{p_end}
{synopt:{helpb dea_stackelberg}}Stackelberg game DEA{p_end}
{synopt:{helpb dea_sw2007}}Simar-Wilson two-stage DEA{p_end}
{synopt:{helpb dea_window}}Window analysis DEA{p_end}
{synoptline}

{marker publication}{...}
{title:Related article and software version}

{pstd}
Lee, C. Forthcoming. {it:dea: An Update and Extension to Data Envelopment Analysis}.
Accepted for publication in the {browse "https://www.stata-journal.com/":Stata Journal}.
{p_end}

{pstd}
This author-maintained software includes subsequent updates and may differ from
the software submitted with the article. The acceptance notice refers to the
article; it does not extend the article's validation claims to later code changes.
{p_end}
