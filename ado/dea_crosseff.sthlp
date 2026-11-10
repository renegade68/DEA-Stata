{smcl}
{* Help revision 2026-11-02: explicit Lim-Zhu input VRS formula and diagonal convention.}{...}
{* Help revision 2026-11-02: validated nonnegative weights and explicit stage diagnostics.}{...}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* *! version 1.10.5 10nov2026}{...}
{viewerjumpto "Related article" "dea_crosseff##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_bootstrap" "help dea_bootstrap"}{...}
{vieweralsosee "dea_malmquist" "help dea_malmquist"}{...}
{vieweralsosee "dea_window" "help dea_window"}{...}
{viewerjumpto "Syntax" "dea_crosseff##syntax"}{...}
{viewerjumpto "Description" "dea_crosseff##description"}{...}
{viewerjumpto "Options" "dea_crosseff##options"}{...}
{viewerjumpto "Remarks" "dea_crosseff##remarks"}{...}
{viewerjumpto "Examples" "dea_crosseff##examples"}{...}
{viewerjumpto "Stored results" "dea_crosseff##stored_results"}{...}
{viewerjumpto "References" "dea_crosseff##references"}{...}
{viewerjumpto "Authors" "dea_crosseff##authors"}{...}
{pstd}Cross-efficiency 1.10.5 uses native IPM by default. Explicit {cmd:engine(simplex)} adds bounded original-data optimal-face certification, failure-only recovery and stage diagnostics. Install the complete matched {cmd:ado/} directory; see INSTALL_CROSSEFF.txt in the package root. Secondary objectives exclude the evaluator before summation. Simplex certificates use the exact intended peer sums; numerical rounding is diagnosed separately. No acceleration is enabled.{p_end}

{pstd}Version 1.10.5 corrects internal point-export failure reasons. Rejected exports remain unavailable, and the score, weight-selection and fallback policies are unchanged. Upgrade the complete matched package.{p_end}

{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:dea_crosseff} {hline 2}}Cross-efficiency DEA{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea_crosseff}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}
[{cmd:,} {it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_crosseff}
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

{syntab:Model}
{synopt:{cmdab:r:ts(}{it:rtstype}{cmd:)}}returns to scale; default is {cmd:rts(crs)}{p_end}
{synopt:{cmdab:or:ientation(}{it:type}{cmd:)}}orientation; default is {cmd:orientation(input)}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}
{synopt:{cmd:method(}{it:methodname}{cmd:)}}secondary objective method; default is {cmd:method(arbitrary)}{p_end}

{synopt:{cmdab:vrsc:ross(}{it:rule}{cmd:)}}VRS peer formula: {cmd:legacy} (default) or {cmd:limzhu}; latter requires input VRS{p_end}

{synopt:{cmdab:optima:diag}}additional bounded exact optimal-face weight and denominator diagnostics; reported DEA results are retained{p_end}

{synopt:{cmdab:cer:ange}}additional bounded exact CE ranges; implies optimadiag; retains existing results{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:newvar}{cmd:)}}name for cross-efficiency variable; default is {cmd:_ce_eff}{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable (string or numeric){p_end}
{synopt:{cmd:matrix}}display full cross-efficiency matrix{p_end}
{synopt:{cmd:detail}}display detailed results{p_end}

{syntab:Technical}
{synopt:{cmdab:eng:ine(}{it:ipm|simplex}{cmd:)}}solver; default {cmd:ipm}; no automatic engine switch{p_end}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmd:debug}}display debugging information{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:rtstype}: {cmd:crs} | {cmd:vrs}{p_end}
{p 4 6 2}
{it:type}: {cmd:input} | {cmd:output}{p_end}
{p 4 6 2}
{it:methodname}: {cmd:arbitrary} | {cmd:benevolent} | {cmd:aggressive}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_crosseff} computes cross-efficiency scores for DEA.
In standard DEA, each DMU selects multiplier weights that maximize its own
efficiency, which may result in extreme weighting schemes and poor discrimination
among efficient DMUs. Cross-efficiency evaluation uses each DMU's optimal weights
to also evaluate all other DMUs, producing a cross-efficiency matrix.

{pstd}
The reported score is the column average of available evaluations by other
selected observations, excluding the diagonal. The matrix uses evaluator rows
and evaluated-unit columns. Alternative optimal weights and ranking ties can
remain, including after a secondary objective is applied.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{cmdab:r:ts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
{cmd:crs} assumes constant returns to scale; {cmd:vrs} assumes variable returns
to scale. Default is {cmd:rts(crs)}. The original cross-efficiency formulation
is CRS-based. VRS uses VRS self-appraisal weights and the peer formula selected
by {cmd:vrscross()}; see the interpretation limits below.

{phang}
{cmdab:or:ientation(}{it:type}{cmd:)} specifies whether to minimize inputs
({cmd:input}) or maximize outputs ({cmd:output}). Default is
{cmd:orientation(input)}.

{phang}
{cmdab:vrsc:ross(}{it:rule}{cmd:)} selects the peer formula. The default
{cmd:legacy} preserves u_i*y_j/(v_i*x_j) in either orientation.
{cmd:limzhu} requires {cmd:rts(vrs)} and input orientation and evaluates
u_i*y_j/(v_i*x_j-w_i), using the input multiplier convention
u_i*y_j+w_i<=v_i*x_j. This implements equation (4) of
{help dea_crosseff##LZ2015:Lim and Zhu (2015)}, with w_i=-xi_i in their notation.
It changes peer appraisal only: the primary LP, secondary weight selection,
self score and off-diagonal averaging remain unchanged. It is not an
implementation of the paper's output-oriented formula. Unsupported combinations
return error 198 before replacing results.

{phang}
{cmd:method(}{it:methodname}{cmd:)} specifies the secondary objective when
multiple optimal weight sets exist.
{cmd:method(arbitrary)} uses the standard DEA solution without a secondary
objective ({help dea_crosseff##SSH1986:Sexton et al. 1986}).
The two secondary methods retain the first-stage objective. For input
orientation, {cmd:method(benevolent)} minimizes the sum over peers of
v*x_j-u*y_j; {cmd:method(aggressive)} maximizes it. For output orientation,
they respectively minimize or maximize the sum of v*x_j. The VRS intercept
is excluded from these secondary objectives. These are the command's existing
linear surrogate objectives, not optimization of the sum of peer ratios.
They remain unchanged with {cmd:vrscross(limzhu)} and do not optimize the
corrected peer ratios; they are not deaR Method II/III under VRS.
The input CRS objective corresponds to the aggregate-slack formulation
(Method II) discussed by {help dea_crosseff##DG1994:Doyle and Green (1994)};
the output and VRS conventions are not interchangeable with every published
secondary formulation. Default is {cmd:method(arbitrary)}.

{dlgtab:Output}

{phang}
{cmdab:gen:erate(}{it:newvar}{cmd:)} specifies the name of the saved mean
cross-efficiency variable. The command also saves {it:newvar}{cmd:_self}
and {it:newvar}{cmd:_maverick}. Default is {cmd:_ce_eff}.

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies a variable identifying each DMU.
Accepts both string and numeric variables.

{phang}
{cmd:matrix} displays and returns the full n×n cross-efficiency matrix, along
with selected weights, stage statuses and available-peer counts. Matrix rows
follow the selected observation order; the weights are in the original data
units. Stata's matrix size limits apply.

{phang}
{cmd:detail} displays the weight-stage status and available-peer count for each
selected observation.

{phang}
{cmd:tolerance(}{it:#}{cmd:)} must be positive and finite. From version 1.4.5,
requested tolerance is applied to normalized LP constraint and bound checks.
Native convergence uses min(tolerance,1e-12) to reduce error before fixing the
first-stage objective in the second LP. The first-stage value is not rounded
to 1. Column and evaluated-row reparameterization and positive row rescaling
preserve the mathematical model. This tolerance does not bound peer-ratio
error or certify exact optimality. Earlier versions did not pass this option
to the native optimizer.


{marker remarks}{...}
{title:Remarks}

{pstd}
With {cmd:vrscross(legacy)}, peer ratios are u_i*y_j/(v_i*x_j) in both orientations. A finite nonnegative
numerator and finite positive denominator are required. CRS keeps small
positive denominators; VRS retains the legacy requirement v_i*x_j>tolerance
under the command's self-normalized weights. Ineligible ratios are missing.
The VRS cutoff limits division by numerically tiny denominators but does not
resolve the underlying interpretation or non-uniqueness. The default averaging
excludes the evaluator's self-appraisal and uses only available peers. One
selected observation can provide a self score but has no peer mean.

{pstd}
With {cmd:vrscross(legacy)} under CRS, the diagonal agrees with the self-efficiency score up to numerical
error. Under VRS, self-appraisal includes the free intercept but peer ratios
omit it. Thus VRS peer ratios may exceed 1 and their diagonal need not equal
the self score. This is an exploratory ratio convention, not the translated
coordinate/corrected VRS method of {help dea_crosseff##LZ2015:Lim and Zhu (2015)}.
Use {cmd:vrscross(limzhu)} for the input VRS correction. Legacy VRS peer means
and the maverick ratio should be interpreted in light of this distinction.

{pstd}
With {cmd:vrscross(limzhu)}, the diagonal uses the same translated appraisal
formula as every other cell; it can differ from the ordinary VRS self score
in {it:generate}{cmd:_self}. This follows the paper's distinction between
simple efficiency and translated appraisal. The command still averages only
available off-diagonal evaluations, whereas the paper also includes the
translated diagonal in its mean. R package deaR's input correction uses the
same off-diagonal formula but replaces the diagonal with simple efficiency.
Alternative optimal weights can produce different matrices in either package.

{pstd}
The corrected input formula is in [0,1] for feasible nonnegative weights when
its denominator is positive. This implementation allows zero weights (rather
than enforcing a non-Archimedean positive bound), so positive data alone do
not ensure a defined ratio. Zero divided by a positive denominator is zero;
0/0 remains missing. Let scale=max(abs(v_i*x_j),abs(w_i),u_i*y_j).
A corrected cell requires a finite nonnegative numerator and weighted input,
a finite denominator greater than zero and tolerance*scale, and a finite ratio
no greater than 1+tolerance. This relative guard rejects cancellation-sensitive
or out-of-range ratios; it does not clip them or substitute small positive
values. Ratios slightly above 1 within tolerance remain unrounded.
Missing off-diagonal and diagonal cells are counted separately. Missingness can
also result from unavailable primary weights; it is not proof of infeasibility.

{pstd}
The generated maverick variable retains the legacy definition self/peer_mean;
it is missing when either score is unavailable or the peer mean is zero.
It is not the alternative definition abs(self-peer_mean)/peer_mean used by
some packages. With alternative optimal weights, different solvers may
produce different peer matrices while satisfying the same objectives.

{pstd}
All selected quantities must be nonnegative, with at least one positive input
and one positive output per observation. Individual zero components are allowed.
Empty samples, invalid quantities and invalid tolerances are rejected before
result replacement. From version 1.4.5 weights have lower bound zero; the old
fixed 1e-6 lower bound could change the model or make it infeasible after a unit
change. Old peer scores and ranks may therefore change and should be recomputed.

{pstd}
An unavailable primary solve leaves that evaluator's weights and self score
missing; the observation can still receive scores from other evaluators.
If a requested secondary solve is unavailable, the accepted primary weights
are retained and the fallback is explicitly counted and displayed. A fallback
is not completion of the requested secondary objective. Secondary LPs may be
unbounded (particularly with zero quantities or the VRS convention), or may
fail numerically when fixing an approximate primary value. Counts do not
certify which cause applies. Available-peer counts are returned because means
based on different sets of evaluators are not directly comparable without care.


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
{phang2}{cmd:. input str6 firm capital labor output}{p_end}
{phang2}{cmd:"Firm A" 12  8  95}{p_end}
{phang2}{cmd:"Firm B" 18 12 125}{p_end}
{phang2}{cmd:"Firm C" 15 10 118}{p_end}
{phang2}{cmd:"Firm D" 22 15 132}{p_end}
{phang2}{cmd:"Firm E" 10  7  90}{p_end}
{phang2}{cmd:"Firm F" 20 13 128}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Basic cross-efficiency with arbitrary secondary objective

{phang2}{cmd:. dea_crosseff capital labor = output, dmu(firm)}{p_end}

{pstd}Benevolent cross-efficiency

{phang2}{cmd:. dea_crosseff capital labor = output, dmu(firm) method(benevolent) matrix replacevars}{p_end}

{pstd}VRS cross-efficiency

{phang2}{cmd:. dea_crosseff capital labor = output, dmu(firm) rts(vrs) replacevars}{p_end}

{pstd}Lim-Zhu input VRS peer formula, preserving the separate ordinary self score

{phang2}{cmd:. dea_crosseff capital labor = output, dmu(firm) rts(vrs) vrscross(limzhu) gen(lz) matrix detail}{p_end}
{phang2}{cmd:. dea_crosseff, in(capital labor) out(output) dmu(firm) r(v) ort(i) vrsc(limzhu) gen(lzshort)}{p_end}
{phang2}{cmd:. assert lz==lzshort}{p_end}

{pstd}String DMU identifiers

{phang2}{cmd:. dea_crosseff capital labor = output, dmu(firm) method(aggressive) replacevars}{p_end}



{pstd}Equivalent sample restrictions in both syntax forms

{phang2}{cmd:. dea_crosseff capital labor = output if inrange(_n,2,5), dmu(firm) generate(qa)}{p_end}
{phang2}{cmd:. dea_crosseff if inrange(_n,2,5), inputs(capital labor) outputs(output) dmu(firm) generate(qb)}{p_end}
{phang2}{cmd:. assert qa==qb}{p_end}


{pstd}Documented minimum option abbreviations

{phang2}{cmd:. dea_crosseff, inputs(capital labor) outputs(output) dmu(firm) rts(vrs) orientation(output) generate(fullopt)}{p_end}
{phang2}{cmd:. dea_crosseff, in(capital labor) out(output) dmu(firm) r(vrs) or(output) gen(shortopt)}{p_end}
{phang2}{cmd:. assert fullopt==shortopt}{p_end}

{pstd}Original-data certificate and numerical diagnostics

{phang2}{cmd:. dea_crosseff capital labor = output, dmu(firm) method(benevolent) engine(simplex) gen(certdemo) matrix detail}{p_end}
{phang2}{cmd:. matrix list r(computation_diagnostics)}{p_end}
{phang2}{cmd:. matrix list r(primary_cert_diagnostics)}{p_end}
{phang2}{cmd:. matrix list r(secondary_cert_diagnostics)}{p_end}
{phang2}{cmd:matrix list r(primary_recovery_diagnostics)}{p_end}
{phang2}{cmd:matrix list r(secondary_recovery_diagnostics)}{p_end}
{phang2}{cmd:. assert r(diagnostic_version)==6201}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea_crosseff} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:r(n)}}number of DMUs{p_end}
{synopt:{cmd:r(mean_ce)}}mean of available peer means{p_end}
{synopt:{cmd:r(mean_self)}}mean of available self-appraisal scores{p_end}
{synopt:{cmd:r(n_failed_primary)}}evaluators with unavailable primary weights{p_end}
{synopt:{cmd:r(n_failed_secondary)}}requested secondary solves that used primary fallback{p_end}
{synopt:{cmd:r(n_failed_ce)}}observations without an available peer mean{p_end}
{synopt:{cmd:r(n_missing_peer_pairs)}}unavailable off-diagonal ratios, out of n*(n-1){p_end}
{synopt:{cmd:r(n_missing_diagonal)}}unavailable translated/legacy diagonal cells, out of n{p_end}
{synopt:{cmd:r(tolerance)}}requested feasibility tolerance{p_end}
{syntab:Macros}
{synopt:{cmd:r(method)}}secondary objective used{p_end}
{synopt:{cmd:r(rts)}}returns-to-scale assumption{p_end}
{synopt:{cmd:r(orientation)}}model orientation{p_end}
{synopt:{cmd:r(peer_rule)}}{cmd:ratio_no_intercept} or {cmd:limzhu_input}{p_end}
{synopt:{cmd:r(vrscross)}}{cmd:legacy} or {cmd:limzhu}{p_end}
{syntab:Matrices}
{synopt:{cmd:r(ce_matrix)}}cross-efficiency matrix; stored with {cmd:matrix}{p_end}
{synopt:{cmd:r(weights_input)}}n×m selected input weights in original units; with {cmd:matrix}{p_end}
{synopt:{cmd:r(weights_output)}}n×s selected output weights in original units; with {cmd:matrix}{p_end}
{synopt:{cmd:r(intercept)}}n×1 VRS intercept, zero for CRS; with {cmd:matrix}{p_end}
{synopt:{cmd:r(stage_status)}}n×1: 0 primary unavailable, 1 primary only, 2 secondary accepted, 3 primary fallback; with {cmd:matrix}{p_end}
{synopt:{cmd:r(peer_count)}}n×1 available off-diagonal evaluations; with {cmd:matrix}{p_end}
{p2colreset}{...}


{marker warnings}{...}
{title:Warnings and unavailable results}

{pstd}An unresolved primary failure leaves that evaluator's primary weights and self
score missing. Explicit simplex first attempts the bounded raw recovery below. A failed requested secondary solve retains the accepted
primary weights. Undefined peer ratios remain missing, and peer means use
available off-diagonal evaluations. These conditions, positive values lost
during scaling, and unverified requested certificates print {cmd:Warning:}
messages. These warnings also appear under {cmd:quietly}; returned counters
remain available for scripts. Use {cmd:detail} and {cmd:matrix} to locate affected observations;
do not treat partial peer means as evaluations by every DMU.{p_end}

{pstd}Native IPM checks both stages' prepared coefficient and right-hand-side
arrays before setting up the optimizer. If intermediate arithmetic has made
an array nonfinite, that stage is skipped and a warning reports the count.
Other evaluators continue. This check does not recover the missing solution
or establish infeasibility of the original-data model. Missing bound markers
used for free or unbounded variables are intentional.{p_end}

{pstd}Always returned are {cmd:r(has_numerical_warning)} (0/1),
{cmd:r(numerical_diagnostic_version)} (6601),
{cmd:r(n_native_nonfinite_primary)}, and {cmd:r(n_native_nonfinite_secondary)}.
The last two count skipped native stages and are zero with simplex;
{cmd:r(native_numerical_scope)} is {cmd:prepared_array_preflight} or
{cmd:not_evaluated_simplex}. Zero warnings does not certify exact optimality
or guarantee representable original-unit weights for every dataset.{p_end}

{pstd}With {cmd:matrix}, {cmd:r(native_numerical_diagnostics)} has columns
{cmd:obs primary_state primary_rc secondary_state secondary_rc}.
Rows follow selected observation order; {cmd:obs} is the original observation
number, so duplicate or string DMU labels can be joined to the data.
Native state 0 means not attempted, 1 means accepted by the existing numerical
checks, 2 means a nonfinite prepared array stopped the stage before the
optimizer, 3 means the numerical solver result was not accepted, and 4 means
the primary self score failed the positive finite check. State 4 is used only
for the primary. Solver return codes are missing if the optimizer was not
called, 0 for acceptance, and otherwise the existing helper's return code.
They do not prove infeasibility. Simplex state/return-code columns are missing;
use its existing solver and certificate diagnostics. {cmd:detail} prints native
stage failures with observation numbers and DMU labels.{p_end}

{pstd}No automatic tolerance increase, engine switch, replacement by a score
of one is performed. Explicit simplex has the failure-only raw recovery described below. Scaling loss, intermediate
overflow, solver failure, and the finite range of normalized weights remain
distinct limits. Inspect the data units and diagnostics before interpreting
or comparing results.{p_end}

{title:Peer objective construction}

{pstd}Benevolent and aggressive objectives sum only the other selected DMUs.
The evaluator is excluded before addition, avoiding cancellation in
total-minus-self arithmetic. The surrogate definitions remain unchanged:
input uses aggregate peer input minus output appraisal, output uses peer
input appraisal, and aggressive reverses the objective sign.{p_end}

{pstd}With {cmd:engine(simplex)}, the numerical objective proposes a basis;
bounded exact stationarity, exchanges and recession checks use the exact
sum of the original binary64 peer quantities after the documented variable
map. A rounded numerical vector is never treated as the exact intended
objective. Objective construction is charged to the existing work limits.
Unverified or excluded stages retain primary fallback. The initially accepted primary
and its optimal face are preserved.{p_end}

{pstd}Native IPM uses a finite dyadically scaled peer-only numerical objective.
It retains its numerical acceptance policy and does not load simplex or
claim exact secondary optimality. Rounded or lost objective terms are
reported. Tiny secondary-objective differences can select different weights
and peer scores even when the objective error is below tolerance. A failed
secondary solve or unavailable objective retains primary weights. Neither
the engine nor tolerance is automatically changed.{p_end}

{pstd}Always returned: {cmd:r(objective_diagnostic_version)} (7001),
{cmd:r(n_objective_attempted)}, {cmd:r(n_objective_build_failed)},
{cmd:r(n_objective_numeric_rounded)}, {cmd:r(n_objective_numeric_lost)},
{cmd:r(n_objective_exact_verified)}, and {cmd:r(n_objective_warning)}.
The latter contributes to {cmd:r(has_numerical_warning)}. Numerical rounding
by itself does not raise that flag when an intended-objective exact
certificate was obtained.{p_end}

{pstd}With {cmd:matrix}, {cmd:r(secondary_objective_diagnostics)} has columns
{cmd:obs attempted build_state numeric_state scope build_ops max_limbs}.
Build state 0 is not attempted, 1 built, -1 unavailable, -2 outside the
bounded domain. Numerical state 0 is not attempted, 1 exactly represented,
2 rounded, 3 lost positive terms, 4 unavailable. Scope 0 is not attempted,
1 intended original objective exactly certified, 2 native numerical-only,
3 requested exact scope unverified. Exact scope describes the certificate,
not a guarantee that rounded exported weights were applied; inspect the
existing stage and link diagnostics as well. {cmd:obs} is the original
observation number. Counts are available without matrices. {cmd:detail}
prints exceptional objective rows.{p_end}

{marker diagnostics}{...}
{title:Numerical and certificate diagnostics}

{pstd}
The default {cmd:engine(ipm)} keeps the native numerical calculation.
Explicit {cmd:engine(simplex)} preserves the numerical primary score and, for
{cmd:method(benevolent)} or {cmd:method(aggressive)}, attempts the bounded
original-data optimal-face certificates described below and in INSTALL_CROSSEFF.txt.
The bounds are 16 variables after free-variable splitting and 32 observations.
Unavailable certificates retain primary weights. Arbitrary method skips
certification. A numerical status is not an exact infeasibility certificate.

{pstd}
Always returned scalars are {cmd:r(diagnostic_version)} (6201),
{cmd:r(n_normalization_lost_positive)}, {cmd:r(n_normalization_lost_inputs)},
{cmd:r(n_normalization_lost_outputs)}, and {cmd:r(n_normalization_lost_obs)}.
They count original positive quantities that become zero in a separate
observation of the existing global column-max scaling. Original zeros are
not counted. No transformed diagnostic data are passed to a solver.
Input and output roles are counted separately, including when the same dataset
variable is used on both sides.
{cmd:r(n_numeric_primary_unavailable)} equals the existing primary-failure count.
This reports numerical availability, not original-model feasibility.

{pstd}
Certificate counters are {cmd:r(n_cert_primary_attempted)},
{cmd:r(n_cert_primary_verified)}, {cmd:r(n_cert_secondary_attempted)},
{cmd:r(n_cert_secondary_verified)}, {cmd:r(n_cert_repair_eq_stops)},
{cmd:r(n_cert_repair_both_stops)}, and {cmd:r(n_cert_mapping_failed)}.
Native calls return zero certificate counters and scope
{cmd:r(certificate_diagnostic_scope)}={cmd:not_evaluated_native}.
The simplex scope is {cmd:bounded_original_data_optimal_face}; it does not
mean that every requested certificate was verified.

{pstd}
With {cmd:matrix}, {cmd:r(computation_diagnostics)} contains original
observation number, numerical primary availability (0/1), existing weight-stage
status, and number of positive quantities lost for that observation.
{cmd:r(normalization_input_loss)} and {cmd:r(normalization_output_loss)} are
0/1 flag matrices with the original variable names. All diagnostic matrices
follow the selected observation order, including after {cmd:if}/{cmd:in}.

{pstd}
{cmd:r(primary_cert_diagnostics)} and {cmd:r(secondary_cert_diagnostics)} each
have nine columns: {cmd:attempted status reason repair_stop repair_steps maxlimbs eq_violation iq_violation dual_violation}. Attempted is 0/1; unattempted fields
are missing, except status -2 for dimension exclusion and -1 for an unavailable
coordinate map. A certified status 1 refers to the internal original-model
fraction, while the returned doubles are checked separately. Status 0 is an
unverified basis proposal; -1 means unknown. Native calls return attempted 0
and missing status/reason/violation fields. No exact feasibility label is inferred
when the numerical primary is unavailable and certification was not attempted.

{pstd}
Certificate reason 0 means verified, 1 means proposal unavailable or outside
supported inputs, 2 means rank/basis selection unavailable, 3 means the proposed
point or dual was rejected; 11 and 12 indicate operation and integer-width
limits, and 13--15 indicate arithmetic, representation, or identity-check failures.
These are certificate-process reasons, not original-model solution classifications.

{pstd}
Repair stop is missing when no repair ran, 0 on repair success,
1 when the proposal is no longer eligible, 2 when the proposed basis has both
primal and dual violations, 3 when it violates an original equality,
4 when no admissible exact exchange is available, 5 on a repeated basis,
6 at the pivot cap, 7 at the repair operation cap, 8 on repair arithmetic/resource
failure, and 9 when an updated certificate is unavailable.
Known equality, inequality/bound, and dual-sign violation flags are 0/1;
missing means not established. A violation flag concerns the proposed basis.

{pstd}
{cmd:r(certification_map_reason)} is missing before a map is attempted,
0 for a successful dyadic map, 1 for the exponent-step cap, 2 for an invalid
scale, and 3 for coefficient round-trip loss. Maps are used for certification,
not for replacing the preserved numerical primary model.
Existing {cmd:r(link_diagnostics)} and {cmd:r(solver_diagnostics)} keep their
existing columns and values. Link columns are {cmd:p_cert s_cert applied reason p_ops s_ops p_solver s_solver p_reason face_rows}; link reason 0 means not
requested or numerical primary unavailable, 1 applied, 2 primary uncertified,
3 primary-score inconsistency, 4 secondary uncertified, 5 export/score rejection,
and 6 dimension exclusion. A raw primary certificate may therefore be verified
while the secondary remains unapplied. Solver diagnostics keep the primary
and secondary blocks {cmd:status solver obj iters scaled absolute}; these
numerical fields do not establish original-data feasibility.
{cmd:detail} additionally prints affected observations,
variables, and unverified-certificate repair stops. New counters are also
available without {cmd:matrix}; diagnostic matrices are returned only with it.

{marker references}{...}

{marker rawrecovery}{...}
{title:Failure-only original-unit recovery (simplex)}

{pstd}Version 1.10.0 attempts a bounded original-data primary LP only when
the initial simplex primary has no accepted weights. Accepted evaluator weights,
self scores, and existing peer rows retain their original computation.
The default native IPM engine is unchanged. Recovery is restricted to
32 selected observations and 16 variables after free-variable splitting.{p_end}

{pstd}At most two exact dyadic representations propose a feasible point.
Version 1.10.3 preserves an accepted proposal. If neither proposal is accepted,
it selects the last valid lossless representation (column-centred when valid)
and attempts the same bounded exact primary fit once. Numerical solver failure
does not by itself prevent an exact certificate. There is at most one exact
primary fit per evaluator; no fresh budget is granted to try another representation.
A merely feasible numerical point is insufficient. Coefficients must round-trip without loss.
The requested tolerance and engine are unchanged. Certified fractions are
converted directly to original-unit weights, combining the free-variable
split before rounding. Normalization, term-relative inequalities, bounds and
objective are checked again in original units. Unexportable weights leave the
primary unavailable; the command does not publish a score without weights.{p_end}

{pstd}Finite positive data do not guarantee representable original-unit
weights. A mathematically finite normalized weight may exceed Stata's finite
double range or a nonzero weight may round to zero. When no representable weights
can satisfy normalization, extra iterations or tolerance changes cannot restore
the current original-unit output. This is a representation limit, not proof that
the mathematical LP is infeasible. The export-obstruction count records failed
proposals or certified-point conversions; by itself it does not prove that every
alternative optimum is unexportable. Alternative units or scaled-weight returns
require separate validation and are not substituted automatically.{p_end}

{pstd}Version 1.10.1 can decode finite original-data coefficients outside
the earlier 2^-900 to 2^900 absolute range in this failure-only primary path,
including subnormal coefficients. A secondary receives the broader domain only
when its primary model required it. Existing accepted evaluator selections and
the default coefficient domain of other certificate/diagnostic paths are
preserved. Arithmetic limits are unchanged: a broader coefficient domain does
not guarantee a certificate or a representable result. Such attempts retain
the existing unknown/fallback states and warnings.{p_end}

{pstd}Version 1.10.2 factors common powers of two before selected exact
products, differences, quotients, and residual checks that would otherwise exceed
the arithmetic width limit. Integer quotients and restored values are checked
exactly; this does not change the data, tolerance, pivot tie rules, or the
128-limb and work limits. The bounded checks can certify additional secondary
optima or improving rays. An exactly unbounded secondary has no finite optimum:
the certified primary weights are retained and the exact-ray note is displayed.
Other resource-limited attempts remain unknown and retain their warnings.
This is not a guarantee for every numerical scale or problem dimension.{p_end}

{pstd}A requested secondary uses the primary's exact dual-defined optimal
face, not a rounded self equality. Existing certificate budgets and the
peer-sum objective remain unchanged. Unknown, unbounded or unexportable
secondary results retain certified primary weights (stage3).
Unknown does not mean that the original model is infeasible.{p_end}

{pstd}Only recovered evaluator rows use bounded exact products/sums of the
exported binary64 raw weights for peer ratios. CRS still requires a positive
denominator; legacy VRS requires denominator greater than the original
tolerance and omits the intercept. Lim-Zhu retains its positive translated
denominator, relative cutoff and [0,1+tolerance] rule. Comparisons use exact
arithmetic on the supplied weights and binary64 tolerance. The diagonal is
calculated by the selected formula and is not replaced by self.
An unrepresentable nonzero ratio stays missing. Existing available-peer means
and ranking rules are unchanged; newly available rows can change those means.
Recovery does not repair scaling losses in already accepted rows.{p_end}

{pstd}A failed raw secondary certificate with a recorded integer-width refusal may
make one additional bounded attempt after exact positive common-factor removal
from its intended objective. Original objective coefficients and their binary
exponent are retained. The positive factor changes certificate units only;
primal constraints, maximizing points and improving-ray signs are unchanged.
The common factor is discovered in Mata with exact integer division, not a
floating approximation. Existing accepted secondary solutions are retained.
The extra attempt shares prior work and the 16-exchange limit; it does not
reset resource limits. A finite certificate is used only after original-unit
weight conversion and guard checks. An exact unbounded secondary objective
retains the certified primary weights and the ray diagnostic.{p_end}

{pstd}Simplex returns {cmd:r(raw_recovery_version)}=8401,
{cmd:r(raw_recovery_scope)}={cmd:failed_primary_only_certified_raw_weights},
and {cmd:r(raw_peer_scope)}={cmd:recovered_rows_only_exported_weights}.
Scalars {cmd:r(n_raw_attempted)}, {cmd:r(n_raw_recovered)},
{cmd:r(n_raw_unavailable)}, {cmd:r(n_raw_primary_unknown)},
{cmd:r(n_raw_export_obstructed)}, {cmd:r(n_raw_secondary_unknown)} and
{cmd:r(n_raw_secondary_unbounded)} count evaluator rows.
{cmd:r(n_raw_peers_evaluated)}, {cmd:r(n_raw_peers_unavailable)},
{cmd:r(n_raw_peers_range)} and {cmd:r(n_raw_peers_arithmetic)} count cells
including the diagonal. These returns are absent with native IPM.{p_end}

{pstd}With {cmd:matrix}, {cmd:r(raw_recovery_diagnostics)} has columns:
{cmd:obs attempted attempts representation primary_status primary_reason primary_export secondary_status secondary_reason secondary_export final_stage stop eq_error iq_error obj_error scale_ops export_ops export_obstructed}.
Representation1 is evaluator-centred and2 is column-centred. Zero means
no selected representation. Exact statuses/reasons use the existing certificate
convention. Export flags mean both finite conversion and raw guard acceptance.
The obstruction flag records a numeric proposal or certified point with a
nonrepresentable raw weight, even when a later attempt has another final stop.{p_end}

{pstd}Stop0 means not attempted;1 primary recovered;2 secondary recovered;
3 secondary unknown with primary fallback;4 exact secondary improving ray;
5 coefficient construction rejected;6 no eligible proposal or representation;
7 primary certificate unknown;8 raw export unavailable;9 raw guard rejected;
10 dimension limit;11 secondary export/guard rejected with primary fallback;
12 secondary objective/map unavailable with primary fallback.
The first observation column preserves the original selected-observation index.
Build/export counts describe binary scaling and fraction mapping, exclude
internal rounding-converter work, and are not elapsed-time measures.{p_end}

{pstd}{cmd:r(raw_scaling_diagnostics)} gives original observation, objective
exponent, input/output variable exponents, VRS intercept exponent if present,
then normalization and peer-row exponents. They define
qc=L*c*D, qA=R*A*D and qb=R*b with positive powers of two.
Unselected representations have missing entries. Outside supported dimensions,
this matrix contains only the original observation column.
{cmd:r(raw_peer_status)} is an evaluator-by-peer matrix:
0 original path/not evaluated;1 finite ratio;2 denominator policy excludes it;
3 ratio export range/loss;4 arithmetic budget;5 invalid operands;
6 Lim-Zhu appraisal outside its allowed interval.{p_end}

{pstd}The existing solver matrix describes the initial normalized numerical
path; it can therefore show a failure for a subsequently recovered row.
Use the new raw matrix for that recovery's outcome. Existing certificate,
recovery and link matrices include the raw path's primary/secondary fits and
retain their schema versions. Exact LP limits are unchanged.
Raw export/guard and peer arithmetic have separate bounded working storage;
they do not enlarge the LP's 128-limb limit.{p_end}

{title:References}

{marker DG1994}{...}
{phang}
Doyle, J., and R. Green. 1994.
Efficiency and cross-efficiency in DEA: Derivations, meanings and uses.
{it:Journal of the Operational Research Society} 45(5): 567--578.
{browse "https://doi.org/10.1057/jors.1994.84":doi:10.1057/jors.1994.84}

{marker LZ2015}{...}
{phang}
Lim, S., and J. Zhu. 2015.
DEA cross-efficiency evaluation under variable returns to scale.
{it:Journal of the Operational Research Society} 66(3): 476--487.
{browse "https://doi.org/10.1057/jors.2014.13":doi:10.1057/jors.2014.13}

{marker SSH1986}{...}
{phang}
Sexton, T. R., R. H. Silkman, and A. J. Hogan. 1986.
Data envelopment analysis: Critique and extensions.
{it:New Directions for Program Evaluation} 1986(32): 73--105.
{browse "https://doi.org/10.1002/ev.1441":doi:10.1002/ev.1441}


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

{title:Recovery path diagnostics}
{pstd}{cmd:r(primary_recovery_diagnostics)} and {cmd:r(secondary_recovery_diagnostics)} are returned with {cmd:matrix}. Each has one row per selected observation, in the order mapped by {cmd:r(computation_diagnostics)}. The 31 columns are:{p_end}
{pstd}{cmd:requested outcome path stop initial_status initial_reason original_steps eq_attempted eq_steps aux_attempted aux_steps final_steps total_steps ray_verified initial_ops ray_ops eq_ops aux_ops final_ops total_ops maxlimbs returned_stop eq_status eq_reason eq_stop aux_status aux_reason aux_stop final_status final_reason final_stop}.{p_end}
{pstd}Requested means a failed initial certificate entered the recovery dispatcher, including a skipped request. Outcome is 0 when recovery was not requested, 1 for a recovered finite optimum, 2 for an exact feasible-point/ray proof, 3 for an unresolved bounded recovery, and -1 for an excluded request. Path 1 is equality-first recovery, 2 is auxiliary feasibility followed by the original objective, 4 is post-failure ray inspection, 5 is the bounded positive-objective-factor attempt, and 0 is no successful recovery path. Path 5 can end in a finite optimum, an exact ray or an unresolved result; inspect outcome as well. Stop 0 means recovery/proof returned, 1 means eligibility or prior resource exclusion, 2 means no verified recovery, and 3 means unsupported canonical form. These process codes do not establish infeasibility.{p_end}
{pstd}Stage status/reason and stop use the certificate/repair codes above. Missing stage status means that stage was not evaluated; a rank-selection failure may have status -1 and reason 2 before any repair. Returned_stop describes the retained certificate. The eq/aux/final stop columns describe the actual stage, even when its trial was not retained.{p_end}
{pstd}Total_steps equals original_steps+eq_steps+aux_steps+final_steps and never exceeds the configured 16 exchanges. Total_ops equals initial_ops+ray_ops+eq_ops+aux_ops+final_ops. Ray inspection has no pivots and runs only after an unsuccessful bounded recovery path; all work shares the existing 80-million operation allowance. Individual checks retain the 10-million operation and 128-limb settings. The counter can include the last rejected integer operation. Path 5 accumulates earlier and new work and exchanges. Initial_status, initial_reason, original_steps and initial_ops retain the first attempt. New factor discovery and its first certificate/repair are included in final_steps/final_ops; eq/aux/final status fields describe the latest reached check. Maxlimbs includes rejected widths from earlier attempts, so it can exceed 128 even when the eventual certificate stays within the configured limit.{p_end}
{pstd}Without {cmd:matrix}, summary scalars are still returned: {cmd:r(n_recovery_requested)}, {cmd:r(n_recovery_finite)}, {cmd:r(n_recovery_unknown)}, {cmd:r(n_recovery_skipped)}, {cmd:r(n_recovery_equality)}, {cmd:r(n_recovery_auxiliary)}, {cmd:r(n_exact_unbounded)}, {cmd:r(recovery_integer_ops)}, {cmd:r(recovery_steps)}, and {cmd:r(certificate_integer_ops)}. Counts include both stages; operation counts cover certificate, recovery and ray integer arithmetic; they exclude the numerical LP, ordinary matrix work and final binary64 export, and do not measure elapsed time. {cmd:detail} prints observation/stage path summaries.{p_end}
{pstd}The diagnostic version is 6201. Native calls return zero work/counters and scope {cmd:not_evaluated_native}; they do not load the simplex/certificate engine unless optimadiag or cerange is requested. Unattempted or dimension-excluded rows remain explicit. No infeasibility conclusion is inferred from a failed numerical primary or from a bad proposed basis.{p_end}


{title:Optional optimal-face diagnostics}

{pstd}{cmd:optimadiag} independently diagnoses the original-data primary optimal
face F1 and, for {cmd:method(benevolent)} or {cmd:method(aggressive)}, the secondary
optimal face F2. It retains the selected engine, efficiency scores, weights,
fallback decisions and existing numerical diagnostics. Without optimadiag or cerange no
additional diagnostic solves or optima returns are produced.{p_end}

{pstd}The supported diagnostic domain is at most 32 selected observations and
at most 16 nonnegative transformed columns (inputs + outputs + two columns for
a VRS intercept), with {cmd:vrscross(legacy)}. Other cases retain DEA results and
report unsupported diagnostics. Explicit {cmd:optimadiag} loads the simplex and
exact-certificate helpers even if the reported DEA engine is IPM. This is an
independent diagnostic calculation, not an automatic engine switch.{p_end}

{pstd}Verified positive dual multipliers identify binding original constraints
that define each entire optimal face. Rounded objective values are never used
as new equality right-hand sides. Actual input/output weights and the unsplit
VRS intercept are minimized and maximized on each face. Different artificial
split variables or bases alone do not establish multiple original weights.
Uniqueness requires all actual coordinate ranges to be verified singletons.
Verified variation establishes multiplicity even if another range is unknown.
Unverified ranges are not treated as unique.{p_end}

{pstd}Denominator states concern {it:v X_j}, including diagonal cells. State 1
means strictly positive throughout the face, 2 identically zero, 3 both zero
and positive values are possible, and -1 unverified. The separate cutoff state
uses the existing threshold (0 for CRS, {cmd:tolerance()} for legacy VRS):
1 always above, 2 always at or below, 3 can cross, -1 unverified. Bounds are
rounded for display; exact signs and comparisons determine the states.
These are denominator ranges, not CE ranges.{p_end}

{pstd}Multiple optimal weights need not imply different CE values or ranks.
The optimadiag option alone does not compute CE ranges. Add cerange for
bounded cellwise CE ranges; no option here selects new weights or computes rank ranges.
Always-zero denominators remain undefined. A mathematical exact-face diagnostic
does not automatically certify rounded native output. Stored-weight linkage
is checked separately; ordinary rounding or primary fallback can place stored
weights off F2. Linkage concerns the stored weights, not the independently
rounded self score or the final floating-point CE division.{p_end}

{pstd}Exact work is bounded. Each LP uses the existing certificate/recovery
limits, with separate diagnostic ceilings of 224 LPs and 200 million counted
integer operations per evaluator, and 4096 LPs and one billion operations per
command. A solve starts only with reserve for its existing bounded path.
Additional exact comparisons have their own bounded work. Exhaustion or an
unsupported representation is reported as unverified; no tolerance is enlarged.
These are algorithmic work bounds, not wall-clock guarantees.{p_end}

{pstd}When requested, {cmd:r(optima_diagnostic_version)} is 7301.
{cmd:r(n_optima_requested)} counts requested evaluator/face records, F1 only
for arbitrary or F1 and F2 otherwise. {cmd:r(n_optima_certified)},
{cmd:r(n_optima_unique)}, {cmd:r(n_optima_multiple)}, {cmd:r(n_optima_unknown)},
{cmd:r(n_optima_unsupported)}, and {cmd:r(n_optima_unbounded)} distinguish face
availability and weight classifications. A known unbounded objective has no
finite optimal face; it is not counted as an unknown weight classification.
{cmd:r(n_optima_linked)}, {cmd:r(n_optima_offface)},
{cmd:r(n_optima_no_weights)} and {cmd:r(n_optima_link_unknown)} describe stored
weights. {cmd:r(n_optima_den_positive)}, {cmd:r(n_optima_den_zero)},
{cmd:r(n_optima_den_possible_zero)}, {cmd:r(n_optima_den_unknown)} and
{cmd:r(n_optima_cutoff_variable)} count computed evaluator/face/peer records,
including diagonal entries. Unavailable faces have no peer records.
{cmd:r(n_optima_lp_calls)} and {cmd:r(n_optima_integer_ops)} report work.
{cmd:r(optima_ce_ranges_available)} and {cmd:r(optima_rank_ranges_available)}
are 0 with optimadiag alone (cerange sets the CE request flag to 1). Existing {cmd:r(has_numerical_warning)} retains its prior definition;
optima diagnostics and their warnings are separate.{p_end}

{pstd}With {cmd:matrix}, {cmd:r(optima_faces)} has columns obs, face,
face_state, weight_state, coords_done, den_done, weight_link, public_stage,
reason, row_lp_calls, row_integer_ops. Face is 1 or 2; face_state is 1 verified,
2 exact unbounded, -1 unverified, -2 unsupported. Weight_state is 1 unique,
2 multiple, -1 unverified, 0 not applicable. Weight_link is 1 exactly on face,
2 off face, 3 no finite stored weights, 0 not checked, -1 unverified.
Public_stage retains 0/1/2/3 from the original calculation. Work columns are
cumulative within that evaluator, not additive across its two face rows.
Reason is 0 completed, 1 domain unsupported, 2 exact mapping unavailable,
3 certificate unavailable, 4 diagnostic budget, 5 exact arithmetic/resource
unavailable, 6 prerequisite face unavailable. Completion of a face or weight
classification does not imply completion of every coordinate/denominator/link;
inspect the corresponding states and completed counts.{p_end}

{pstd}{cmd:r(optima_coordinates)} has obs, face, coordinate, min_state,
max_state, range_state, lower, upper. Coordinates follow inputs, outputs, then
the actual VRS intercept. Endpoint states are 1 finite certified, 2 certified
unbounded, -1 unverified; range_state is 1 fixed, 2 variable, -1 unverified.
{cmd:r(optima_denominators)} has obs, face, peer_obs, min_state, max_state,
denominator_state, cutoff_state, lower, upper. Observation numbers are original
dataset row numbers after if/in selection, not the values of dmu(). Missing
display bounds can mean an unbounded or unavailable endpoint; consult states.
Unsupported whole-command diagnostics return one aggregate face row with
obs missing/face 0 and all-missing coordinate/denominator placeholders.
Any other empty detail table also uses an all-missing placeholder.
Without {cmd:matrix}, counts and warnings remain available without these tables.{p_end}

{pstd}Example with multiple optimal weights and proportional quantities:{p_end}
{phang2}{cmd:clear}{p_end}
{phang2}{cmd:input double id x1 x2 y}{p_end}
{phang2}{cmd:1 1 1 1}{p_end}
{phang2}{cmd:2 2 2 2}{p_end}
{phang2}{cmd:3 3 3 3}{p_end}
{phang2}{cmd:end}{p_end}
{phang2}{cmd:dea_crosseff x1 x2 = y, dmu(id) method(benevolent) optimadiag matrix}{p_end}
{phang2}{cmd:matrix list r(optima_faces)}{p_end}
{phang2}{cmd:matrix list r(optima_denominators)}{p_end}


{title:Optional CE ranges}

{pstd}{cmd:cerange} (abbreviation {cmd:cer}) implies {cmd:optimadiag} and adds
bounded exact ranges of the legacy peer ratio {it:u Y_j / v X_j} on each
verified original-data face. F1 is the primary face; F2 is the secondary face
for benevolent/aggressive methods. An unavailable or unbounded secondary
objective produces unavailable F2 rows; F1 is never substituted. Existing
weights, scores, fallback, means and ranks are retained. A range is a property
of the entire exact face, not an error interval around the stored score.{p_end}

{pstd}Ranges use the domain {it:v X_j > 0}. Denominators that are identically
zero give undefined ratios. If zero and positive denominators both occur,
the reported bounds concern only the positive domain. The numerical cutoff
(0 for CRS, tolerance() for legacy VRS) is diagnosed separately; ranges are
not conditional on passing that cutoff. In VRS the legacy ratio omits the
intercept. {cmd:vrscross(limzhu)} is unsupported by this diagnostic.{p_end}

{pstd}At most 32 selected observations and 15 nonnegative original-face
columns (inputs + outputs + two for a VRS intercept) are supported. The
linear-fractional transformation adds one column, keeping the LP limit at 16.
Each exact LP may have at most 64 equality/inequality rows before sign bounds;
the fractional normalization and endpoint checks count toward this limit.
The complete LP has at most 256 rows. CE work has separate ceilings of 512
LPs/200 million counted integer operations per evaluator and 4096 LPs/one
billion per command. F1 and F2 share that evaluator budget. Bounds include
reserved conversion work and do not measure elapsed time. Arithmetic and
representation limits can leave results unknown even inside these dimensions.
In particular, the inherited exact LP decoder does not cover the most extreme
subnormal coefficients. Unknown does not mean zero, infeasible, or unique.
Tolerance is not enlarged, and there is no automatic weight selection.{p_end}

{pstd}Both native engines use the same independent simplex/exact diagnostic
path. Existing {cmd:r(has_numerical_warning)} is unchanged. New
{cmd:r(has_cerange_warning)} includes variable, undefined, unknown, unsupported
or unavailable ranges, limit-only/unbounded endpoints, zero-possible or
cutoff-dependent denominators, display loss, and stored/reconstructed weight
link issues. Warnings and counts remain available without {cmd:matrix},
including under {cmd:quietly}. Exact mathematical states, displayed numeric
bounds and linkage to native weights must be read separately.{p_end}

{pstd}{cmd:r(cerange_version)} is 7701; scope is
{cmd:bounded_original_exact_faces}, engine {cmd:independent_simplex_exact}.
{cmd:r(optima_ce_ranges_available)} becomes 1 when cerange is requested;
this flag does not assert that all cells were verified. Mean/rank range flags
{cmd:r(cerange_mean_ranges_available)}, {cmd:r(cerange_rank_ranges_available)}
and {cmd:r(optima_rank_ranges_available)} remain 0. Cellwise bounds cannot
simply be combined into attainable mean or rank bounds.{p_end}

{pstd}The {cmd:r(n_cerange_*)} scalars are: requested, verified, fixed,
variable, undefined, unknown, unsupported, face_unavailable, limit_only,
unbounded, den_zero_possible, cutoff_variable, display_loss, display_unknown,
offface, no_weights, link_unknown, raw_link_unknown, lp_calls, integer_ops.
Counts include diagonal cells and both faces when applicable. Verified is
fixed + variable. Requested is the sum of verified, undefined, unknown,
unsupported and face_unavailable. Auxiliary counts overlap; limit_only and
unbounded count cells with at least one such endpoint. Work is separate from
{cmd:n_optima_lp_calls}/{cmd:n_optima_integer_ops}.{p_end}

{pstd}With {cmd:matrix}, {cmd:r(ce_ranges)} has these columns in order:{p_end}
{p 8 8 2}obs face peer_obs face_state range_state denominator_state cutoff_state
lower upper lower_state upper_state lower_display upper_display
raw_lower_link raw_upper_link stored_weight_link public_stage reason
row_lp_calls row_integer_ops{p_end}

{pstd}Observation indices are original dataset rows after if/in selection,
not dmu() values. Face is 1 or 2. Face_state is 1 verified finite face,
2 exactly unbounded objective, -1 unverified, -2 unsupported. Range_state is
1 fixed, 2 variable, 0 undefined everywhere, -1 unknown, -2 unsupported,
-3 prerequisite face unavailable. Denominator/cutoff states follow the
optimadiag definitions. Whole-command exclusion returns one aggregate row
with obs/peer_obs missing, face 0, range_state -2 and reason 2; its unsupported
count is the full number of requested cells.{p_end}

{pstd}Endpoint states: 1 attained by finite positive-denominator weights,
2 finite limit only, 3 exactly unbounded, -1 unavailable. Display states:
0 not applicable, 1 exact binary64, 2 rounded, 3 underflow, 4 outside the
conservative export range (absolute value below 2^1023), 5 unverified.
Read endpoint and display states before interpreting missing or zero bounds.
Raw links are 1 for verified reconstruction of an attained endpoint in
original weights, -1 unverified, 0 not applicable. Stored_weight_link is
copied from the corresponding optima face (1 on face, 2 off face, 3 no finite
stored weights, -1 unknown, 0 not checked). This does not certify the rounded
CE division or the native mean/rank. Public_stage retains its native status.
Work columns are cumulative within an evaluator across peers and both faces;
do not add them. Reason codes: 0 completed, 2 unsupported input/dimension,
3 solver unverified, 4 budget, 5 arithmetic/representation, 6 mapping loss,
7 negative denominator, 8 prerequisite face unavailable.{p_end}

{pstd}Self-contained example: multiple optimal input weights, fixed CE ratios.{p_end}
{phang2}{cmd:clear}{p_end}
{phang2}{cmd:input double id x1 x2 y}{p_end}
{phang2}{cmd:1 1 1 1}{p_end}
{phang2}{cmd:2 2 2 2}{p_end}
{phang2}{cmd:3 3 3 3}{p_end}
{phang2}{cmd:end}{p_end}
{phang2}{cmd:dea_crosseff x1 x2=y, dmu(id) method(benevolent) cerange matrix gen(ce)}{p_end}
{phang2}{cmd:assert r(n_optima_multiple)==6}{p_end}
{phang2}{cmd:assert r(n_cerange_fixed)==18 & r(n_cerange_variable)==0}{p_end}
{phang2}{cmd:matrix list r(ce_ranges)}{p_end}
