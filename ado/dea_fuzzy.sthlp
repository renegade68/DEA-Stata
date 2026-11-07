{smcl}
{* Help revision 2026-11-01: correct and document NIRS/NDRS aliases.}{...}
{* Help revision 2026-11-01: if/in defines both evaluation and reference samples.}{...}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-11-01: explicit crisp-DEA identifiers in the sequential example.}{...}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* *! version 1.1.0  03nov2026}{...}
{viewerjumpto "Related article" "dea_fuzzy##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_bootstrap" "help dea_bootstrap"}{...}
{vieweralsosee "dea_window" "help dea_window"}{...}
{viewerjumpto "Syntax" "dea_fuzzy##syntax"}{...}
{viewerjumpto "Description" "dea_fuzzy##description"}{...}
{viewerjumpto "Options" "dea_fuzzy##options"}{...}
{viewerjumpto "Remarks" "dea_fuzzy##remarks"}{...}
{viewerjumpto "Examples" "dea_fuzzy##examples"}{...}
{viewerjumpto "Stored results" "dea_fuzzy##stored_results"}{...}
{viewerjumpto "References" "dea_fuzzy##references"}{...}
{viewerjumpto "Authors" "dea_fuzzy##authors"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:dea_fuzzy} {hline 2}}Fuzzy data envelopment analysis (α-cut approach){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:dea_fuzzy}
{ifin}{cmd:,}
{cmdab:in:puts(}{it:varlist}{cmd:)}
{cmdab:out:puts(}{it:varlist}{cmd:)}
[{it:options}]


{pstd}
{cmd:if} and {cmd:in} restrict both the DMUs evaluated and the reference
set. Observations with missing input or output values are excluded; missing
values in supplied fuzzy-bound variables also exclude observations. Unrelated
variables do not affect sample selection. Generated results are missing outside
the selected sample. An empty estimation sample returns {cmd:r(2000)}.
When {cmd:replacevars} is used, it replaces entire result columns, including
setting observations outside the selected sample to missing.{p_end}

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmdab:in:puts(}{it:varlist}{cmd:)}}center (modal) values of input variables{p_end}
{synopt:{cmdab:out:puts(}{it:varlist}{cmd:)}}center (modal) values of output variables{p_end}

{syntab:Fuzzy bounds}
{synopt:{cmd:xlo(}{it:varlist}{cmd:)}}lower-bound variables for inputs at α = 0; default uses {cmd:inputs()}{p_end}
{synopt:{cmd:xhi(}{it:varlist}{cmd:)}}upper-bound variables for inputs at α = 0; default uses {cmd:inputs()}{p_end}
{synopt:{cmd:ylo(}{it:varlist}{cmd:)}}lower-bound variables for outputs at α = 0; default uses {cmd:outputs()}{p_end}
{synopt:{cmd:yhi(}{it:varlist}{cmd:)}}upper-bound variables for outputs at α = 0; default uses {cmd:outputs()}{p_end}

{syntab:Model}
{synopt:{cmdab:alp:ha(}{it:numlist}{cmd:)}}α-cut levels ∈ [0,1]; default is {cmd:alpha(0 0.25 0.5 0.75 1)}{p_end}
{synopt:{cmdab:r:ts(}{it:rtstype}{cmd:)}}returns to scale; default is {cmd:rts(crs)}{p_end}
{synopt:{cmdab:or:ientation(}{it:type}{cmd:)}}orientation; default is {cmd:orientation(input)}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for efficiency variables; default is {cmd:fdea_eff}{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable{p_end}
{synopt:{cmd:detail}}display DMU-level efficiency table{p_end}
{synopt:{opt mat:rix}}return endpoint, complete-cohort and observation-index matrices{p_end}

{syntab:Technical}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}solver convergence tolerance; default is {cmd:1e-10}{p_end}
{synopt:{cmd:debug}}display debugging information{p_end}
{synoptline}

{p 4 6 2}
{it:rtstype}: {cmd:crs} {cmd:ccr} {cmd:vrs} {cmd:bcc} {cmd:nirs} {cmd:drs} {cmd:ndrs} {cmd:irs}{p_end}
{p 4 6 2}
{it:type}: {cmd:input} {cmd:i} {cmd:output} {cmd:o}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_fuzzy} evaluates triangular fuzzy data using a finite grid of
alpha cuts and the Kao-Liu best/worst reference convention
({help dea_fuzzy##KL2000:Kao and Liu 2000}). Its radial model corresponds to
{cmd:deaR::modelfuzzy_kaoliu(..., kaoliu_modelname="basic")} with the matching
orientation and returns to scale. The evaluated DMU uses the same realization
in its target and reference row. Other fuzzy/possibilistic DEA formulations,
including the cited {help dea_fuzzy##LFJ2003:Lertworasirikul et al. (2003)} and
{help dea_fuzzy##GT2001:Guo and Tanaka (2001)} methods, are background rather
than separate implemented model options.

{pstd}
The joint alpha-cut is the Cartesian product of the supplied marginal
intervals (the usual minimum-membership, noninteractive convention).
Cross-variable or cross-DMU dependence constraints are not supplied to the
model. Each evaluated DMU's extrema generally use different realizations of
the dataset; all DMU bounds need not be attained simultaneously.

{pstd}
A triangular fuzzy number for variable {it:x} is defined by three values
({it:x_lo}, {it:x_center}, {it:x_hi}) where {it:x_center} is the modal
(full-membership) value, {it:x_lo} is the lower support, and {it:x_hi} is the
upper support.  At α-cut level α ∈ [0,1], the fuzzy number is
represented by the crisp interval [{it:x_L}(α), {it:x_U}(α)] where

{p 8 8 2}
{it:x_L}(α) = {it:x_lo} + α × ({it:x_center} − {it:x_lo})

{p 8 8 2}
{it:x_U}(α) = {it:x_hi} − α × ({it:x_hi} − {it:x_center})

{pstd}
At α = 1 both bounds collapse to {it:x_center}, so the result is identical
to standard (crisp) DEA.  At α = 0 the full support [{it:x_lo}, {it:x_hi}]
is used, representing maximum uncertainty.

{pstd}
For each α-level, two LPs are solved per DMU:

{phang2}
{it:Pessimistic LP.}  DMU {it:k} is evaluated at its worst-case values
(high inputs, low outputs).  Peer DMUs enter the reference set at favorable
bounds (low inputs, high outputs), while DMU {it:k} enters the reference set
at the same worst-case realization.  This gives the {it:lower efficiency bound}
θ_L(α).

{phang2}
{it:Optimistic LP.}  DMU {it:k} is evaluated at its best-case values
(low inputs, high outputs).  Peer DMUs enter the reference set at unfavorable
bounds (high inputs, low outputs), while DMU {it:k} enters the reference set
at the same best-case realization.  This gives the {it:upper efficiency bound}
θ_U(α).

{pstd}
Under the stated joint-cut convention and positive-row domain, the pair
[θ_L(α), θ_U(α)] describes the efficiency range at that alpha. Requested cuts
sample the membership function. A finite grid does not determine its values
between cuts, and linear interpolation of plotted points is an approximation.
An alpha level is a membership level, not a probability or confidence level.
No defuzzified score, alpha-weighted average or fuzzy ranking is computed.

{pstd}
For output-oriented models, the raw expansion factor φ (≥ 1) is
converted to 1/φ ∈ (0,1] for consistency with the input-oriented scale,
following the {cmd:dea} package convention for imprecise DEA (IDEA).


{marker options}{...}
{title:Options}

{dlgtab:Fuzzy bounds}

{phang}
{cmd:xlo(}{it:varlist}{cmd:)} specifies the lower-bound variables for inputs
at α = 0 (the full lower support of the fuzzy number).  The number of
variables must equal the number of variables in {cmd:inputs()}, and each
corresponds to the matching input variable in order.  If not specified,
the center values in {cmd:inputs()} are used, giving zero left spread.
An input is crisp only when both bounds equal its center.

{phang}
{cmd:xhi(}{it:varlist}{cmd:)} specifies the upper-bound variables for inputs
at α = 0.  Same correspondence rules as {cmd:xlo()}.

{phang}
{cmd:ylo(}{it:varlist}{cmd:)} specifies the lower-bound variables for outputs
at α = 0.  If not specified, the center values in {cmd:outputs()} are used.

{phang}
{cmd:yhi(}{it:varlist}{cmd:)} specifies the upper-bound variables for outputs
at α = 0.  If not specified, the center values in {cmd:outputs()} are used.

{pstd}
Each supplied bound list must contain one variable for every input or output,
in the matching order. To fuzzify only some variables on that side, put the
center variables themselves in the other positions; for example,
{cmd:inputs(x1 x2) xlo(x1_lo x2) xhi(x1_hi x2)}. An omitted bound option defaults
to the center on that side of the triangle. If only one bound option is supplied,
the resulting fuzzy number may be one-sided.  If none of {cmd:xlo()},
{cmd:xhi()}, {cmd:ylo()}, {cmd:yhi()} are specified, all variables are
treated as crisp and the results are identical to standard DEA at every
α-level.

{dlgtab:Model}

{phang}
{cmd:alpha(}{it:numlist}{cmd:)} specifies the α-cut levels at which to
evaluate fuzzy efficiency.  Each value must be in [0,1].  Multiple values
are supplied as a space-separated list; the list is automatically sorted.
The default is {cmd:alpha(0 0.25 0.5 0.75 1)}.  Use {cmd:alpha(1)} to
obtain only the standard DEA result.

{phang}
{cmd:rts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
{cmd:crs} (constant, the default; synonym {cmd:ccr}) imposes no constraint
on the intensity variables λ.  {cmd:vrs} (variable; synonym {cmd:bcc})
adds sum(λ) = 1. {cmd:nirs} (nonincreasing; alias {cmd:drs}) adds
sum(λ) ≤ 1. {cmd:ndrs} (nondecreasing; alias {cmd:irs}) adds sum(λ) ≥ 1.
The legacy {cmd:drs}/{cmd:irs} spellings select these technology restrictions;
they do not classify an individual DMU's local returns to scale.
The same RTS assumption is
applied to both the pessimistic and optimistic LPs at every α-level.
The normalized specification returned in {cmd:r(rts)} is {cmd:crs}, {cmd:vrs},
{cmd:drs}, or {cmd:irs}.

{pstd}
Versions before 1.0.4 connected {cmd:nirs} and {cmd:ndrs} to the opposite
restrictions. Recompute results previously obtained with either spelling.
Calculations explicitly using {cmd:drs} or {cmd:irs} are unchanged.

{phang}
Specify only one of {cmd:orientation()} and {cmd:ort()}, including their
abbreviations; specifying both returns {cmd:r(198)}, even with identical
values. Versions before 1.0.6 silently preferred {cmd:orientation()}.

{phang}
{cmd:orientation(}{it:type}{cmd:)} or {cmd:ort(}{it:type}{cmd:)} specifies
the model orientation.  {cmd:input} (the default) minimizes the radial
contraction factor θ, giving input efficiency scores in (0,1].
{cmd:output} maximizes the radial expansion factor φ ≥ 1, stored as
1/φ ∈ (0,1].  Abbreviations {cmd:i} and {cmd:o} are accepted.

{dlgtab:Output}

{phang}
{cmd:generate(}{it:prefix}{cmd:)} specifies the variable name prefix for
the efficiency results.  The default prefix is {cmd:fdea_eff}.  For each
α-level, two variables are created:

{p2colset 9 42 44 2}{...}
{p2col :{it:prefix}{cmd:_pess_a}{it:N}}pessimistic efficiency at the specified α{p_end}
{p2col :{it:prefix}{cmd:_opt_a}{it:N}}optimistic efficiency at the specified α{p_end}
{p2colreset}{...}

{pstd}
where {it:N} = round(α × 100).  For example, {cmd:generate(eff)} with
{cmd:alpha(0 0.5 1)} creates {cmd:eff_pess_a0}, {cmd:eff_opt_a0},
{cmd:eff_pess_a50}, {cmd:eff_opt_a50}, {cmd:eff_pess_a100}, and
{cmd:eff_opt_a100}.

{pstd}
When several distinct levels have the same rounded {it:N}, all members of
that group receive {cmd:_1}, {cmd:_2}, etc., in ascending alpha order.
For example, {cmd:alpha(.004 .001)} with {cmd:generate(eff)} creates
{cmd:eff_pess_a0_1}/{cmd:eff_opt_a0_1} for .001 and
{cmd:eff_pess_a0_2}/{cmd:eff_opt_a0_2} for .004.
Other groups retain their historical names. The suffix is a name only:
the solver uses the specified alpha, not {it:N}/100. Adding or removing
levels can change these suffixes. Use the parallel lists {cmd:r(alpha)},
{cmd:r(pessvars)}, and {cmd:r(optvars)} for the exact mapping, and copy them
to local macros before running another r-class command.
All completed names must fit Stata's 32-character variable-name limit;
choose a shorter prefix if necessary.

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies a string or numeric variable
containing DMU identifiers used in output displays.  If not specified,
DMUs are labelled {cmd:dmu1}, {cmd:dmu2}, etc. using original observation
numbers, and the successful command creates {cmd:_dmu}. An existing {cmd:_dmu}
requires an explicit {cmd:dmu()} choice. Here identifiers are display labels:
missing or repeated labels do not remove otherwise eligible observations.
{cmd:r(obs_index)} identifies matrix rows without relying on unique labels.

{phang}
{cmd:detail} displays a DMU-level table of pessimistic and optimistic
efficiency scores at every α-level.

{dlgtab:Technical}

{phang}
{cmd:tolerance(}{it:#}{cmd:)} sets the initial {cmd:LinearProgram.setTol()}
convergence tolerance; the default is {cmd:1e-10}. It must be positive and
nonmissing. Versions before 1.0.5 accepted this option without applying it.
Each input/output column and its bounds are divided by their common maximum
over the selected sample before solving. This preserves the LP under changes
of measurement units. Selected data must satisfy
0 <= lower <= center <= upper; otherwise the command returns {cmd:r(459)}.
Omitted bounds equal the center. Individual zero entries and redundant
all-zero columns are allowed. At the smallest requested alpha, every selected
row must have at least one strictly positive lower-cut input and one strictly
positive lower-cut output. This guarantees those conditions throughout every
requested cut. An all-zero side returns {cmd:r(459)} before result replacement.
A support with an all-zero side at alpha=0 can still be used at strictly positive
requested levels if the lower cuts there satisfy the condition. This domain
check was added in version 1.1.0; earlier versions could return degenerate
scores or missing unbounded output scores for such rows.

{phang}
Converged solutions are checked for finite parameters, scaled primal
constraint/bound residuals, objective consistency, and efficiency in
[0, 1+1e-8]. A failed LP is retried once with
min(current tolerance/100, 1e-12). If both attempts report native status 1 or 2, its score is missing and a note
is displayed. These are solver reports, not mathematical infeasibility or
unboundedness certificates. With the positive-row domain and matching
self-reference, the mathematical radial problems have finite feasible solutions;
a native failure in that domain is retained as an unavailable numerical result.
Other persistent numerical failures return {cmd:r(498)} before any permanent
result variables are changed, including when {cmd:replacevars} is specified.

{phang}
The complete alpha grid is checked for pessimistic <= optimistic and nested
intervals, allowing absolute roundoff of {cmd:1e-8}. Comparisons bridge missing
interior cuts and use all earlier available levels. At alpha=1, both known
endpoints must coincide within that threshold. A violation triggers one
complete grid retry with min(requested tolerance/100, 1e-12); a persistent
violation returns {cmd:r(498)}. Scores are not clipped or sorted to repair
violations. In the underflow case where division by 100 becomes zero, the
positive current tolerance is retained. The convergence tolerance and the
fixed validation threshold are not guarantees of absolute score accuracy or
proofs of optimality. Larger requested tolerances can reduce accuracy.

{phang}
{cmd:debug} displays additional status messages during execution.


{marker remarks}{...}
{title:Remarks}

{pstd}
{it:Relationship to IDEA.}
{cmd:dea_fuzzy} generalizes the Imprecise DEA ({cmd:model(imprecise)}) in
{cmd:dea} to a parametric family of interval problems indexed by α.
At a single α-level the LP structure is identical to IDEA: pessimistic
uses (X_U[k], Y_L[k]) for DMU {it:k} and (X_L, Y_U) for the frontier,
while optimistic reverses the roles.  The novelty of {cmd:dea_fuzzy} is
the α-cut parameterization, which reduces the full fuzzy number to a
linear interpolation between the support bounds and the modal value.

{pstd}
{it:Interpreting efficiency intervals.}
For input-oriented models, the efficiency interval [θ_L(α), θ_U(α)] has
the following interpretation:

{p2colset 9 36 38 2}{...}
{p2col :{it:Condition}}{it:Interpretation}{p_end}
{p2col :{hline 60}}{p_end}
{p2col :θ_L(α) = 1}DMU is {it:certainly efficient} at level α{p_end}
{p2col :θ_L(α) < 1, θ_U(α) ≥ 1}DMU is {it:possibly efficient} at level α{p_end}
{p2col :θ_U(α) < 1}DMU is {it:certainly inefficient} at level α{p_end}
{p2colreset}{...}

{pstd}
The same classification applies for output-oriented models using 1/φ in
place of θ.

{pstd}
{it:Membership function.}
The endpoint values at the requested levels sample the lower and upper
membership branches. The efficiency membership function can be nonlinear even
when every data triangle is piecewise linear. At alpha=1 the LPs use exactly
the same center matrix; at alpha=0 they use the support endpoints. Neither
unrequested cuts nor a continuous interpolated function are returned.

{pstd}
{it:Summary cohort and aggregation.}
The table and {cmd:r(summary)} use the same rows at every alpha: a row is
included only if both endpoints are available at every requested level.
{cmd:r(N_complete)} gives this cohort size and {cmd:r(N_partial)} counts the
other selected rows. Known individual endpoint scores remain available even
for partial rows. Each summary row also shows N_pess, N_opt and N_pair at that
cut. If no row is complete, summary means/ranges are missing. Changing the
requested alpha list can change the complete cohort.

{pstd}
These means are descriptive averages of marginal efficiency bounds. They are
not scores of an aggregated production unit or the exact alpha-cut of the
mean efficiency under one common realization: different DMUs may require
incompatible worst/best realizations. The command does not average across alpha
levels. Earlier versions displayed means from separate available subsets for
the two endpoints, which could produce misleading comparisons after failures.

{pstd}
{opt mat:rix} returns {cmd:r(eff_pess)} and {cmd:r(eff_opt)} with selected
observations in their original order and columns alpha1, alpha2, ... matching
{cmd:r(summary)} and the returned name lists. {cmd:r(complete)} is the 0/1
complete-cohort flag; {cmd:r(obs_index)} contains original observation numbers.
Matrices retain missing scores and full stored precision and are subject to
Stata's matrix size limits. {cmd:r(summary)} is returned without this option.

{pstd}
{it:LP formulation — input-oriented.}
At each α-level, {cmd:dea_fuzzy} solves the following two LPs for each
DMU {it:k}:

{p 8 8 2}
{it:Pessimistic:} {space 2} min θ{break}
subject to  X_W(α)′ λ − θ × x_{k,U}(α) ≤ 0  {space 2} (m input constraints){break}
{space 12}  −Y_W(α)′ λ {space 12} ≤ −y_{k,L}(α)  (s output constraints){break}
{space 12}  [convexity constraint per {cmd:rts()}]{break}
{space 12}  λ ≥ 0,  θ ≥ 0

{p 8 8 2}
{it:Optimistic:} {space 3} min θ{break}
subject to  X_B(α)′ λ − θ × x_{k,L}(α) ≤ 0{break}
{space 12}  −Y_B(α)′ λ {space 12} ≤ −y_{k,U}(α){break}
{space 12}  [convexity constraint per {cmd:rts()}]{break}
{space 12}  λ ≥ 0,  θ ≥ 0

{pstd}
where X_W/Y_W use favorable bounds for peers and DMU {it:k}'s unfavorable
realization for DMU {it:k}; X_B/Y_B use unfavorable bounds for peers and
DMU {it:k}'s favorable realization for DMU {it:k}. This is the
Kao-Liu/deaR convention for best/worst alpha-cut evaluation.

{pstd}
{it:LP formulation — output-oriented.}
The output-oriented model maximizes φ (stored as 1/φ):

{p 8 8 2}
{it:Pessimistic:} {space 2} max φ  (equivalently, min −φ){break}
subject to  X_W(α)′ λ {space 14} ≤ x_{k,U}(α){break}
{space 12}  −Y_W(α)′ λ + φ × y_{k,L}(α) ≤ 0{break}
{space 12}  [convexity constraint per {cmd:rts()}]{break}
{space 12}  λ ≥ 0,  φ ≥ 1

{p 8 8 2}
{it:Optimistic:} {space 3} max φ{break}
subject to  X_B(α)′ λ {space 14} ≤ x_{k,L}(α){break}
{space 12}  −Y_B(α)′ λ + φ × y_{k,U}(α) ≤ 0{break}
{space 12}  [convexity constraint per {cmd:rts()}]{break}
{space 12}  λ ≥ 0,  φ ≥ 1

{pstd}
{it:Solver.}
All LPs are solved using Stata 16's {cmd:LinearProgram()} interior point
method.  The {cmd:dea_fuzzy} command therefore requires Stata 16.0 or later.

{pstd}
{it:Number of LPs solved.}
For {it:n} DMUs, {it:K} α-levels, and two orientations per level
(pessimistic and optimistic), the total number of LPs solved is
2 × {it:n} × {it:K}.  For example, with 50 DMUs and 5 α-levels,
500 LPs are solved.

{pstd}
{it:Verifying results at α = 1.}
Because α = 1 collapses all bounds to the center value, results at
α = 1 are identical for pessimistic and optimistic, and should
match standard DEA results from {cmd:dea} with the same {cmd:rts()} and
{cmd:orientation()}.  This provides a useful cross-check:

{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) ort(i) gen(crs_eff) dmu(dmu)}{p_end}
{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) alpha(1) gen(f) dmu(dmu)}{p_end}
{phang2}{cmd:. assert abs(crs_eff - f_pess_a100) < 1e-6}{p_end}


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

{pstd}
An existing {cmd:_dmu} is never replaced automatically, including with
{cmd:replacevars}; specify {cmd:dmu(_dmu)} to reuse it. Repeated alpha
levels are rejected with {cmd:r(198)} before changing data. This includes
values that become identical when Stata parses the numlist; the interface
retains Stata's numlist precision. Distinct parsed levels sharing a rounded
suffix are supported from version 1.0.7 onward. Their alpha/variable mapping
is displayed automatically; {cmd:detail} displays the mapping for any grid.

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


    {title:Setup}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str7 dmu i_employee i_area o_sales o_profits}{p_end}
{phang2}{cmd:. "store A" 10 20 70 6}{p_end}
{phang2}{cmd:. "store B" 15 15 100 3}{p_end}
{phang2}{cmd:. "store C" 20 30 80 5}{p_end}
{phang2}{cmd:. "store D" 25 15 100 2}{p_end}
{phang2}{cmd:. "store E" 12 9 90 8}{p_end}
{phang2}{cmd:. "store F" 16 18 75 5}{p_end}
{phang2}{cmd:. end}{p_end}

{phang2}{cmd:. // Variables: dmu, i_employee, i_area, o_sales, o_profits}{p_end}
{phang2}{cmd:. // Create ±10% triangular fuzzy spreads around observed values}{p_end}
{phang2}{cmd:. foreach v of varlist i_employee i_area o_sales o_profits {c -(}}{p_end}
{phang2}{cmd:. {space 4}gen `v'_lo = `v' * 0.90}{p_end}
{phang2}{cmd:. {space 4}gen `v'_hi = `v' * 1.10}{p_end}
{phang2}{cmd:. {c )-}}{p_end}

    {title:Basic fuzzy DEA: CRS, input-oriented, 5 α-levels}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) rts(crs) ort(input) gen(crs_f) dmu(dmu)}{p_end}

    {title:VRS fuzzy DEA with DMU labels and detailed table}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) rts(vrs) ort(input) dmu(dmu) gen(vrs_f) detail}{p_end}

    {title:Custom α-levels}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) alpha(0 0.2 0.4 0.6 0.8 1) gen(f6) dmu(dmu)}{p_end}

    {title:Close α-levels and returned variable mapping}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) alpha(.004 .001 1) gen(close) dmu(dmu) detail}{p_end}
{phang2}{cmd:. local levels "`r(alpha)'"}{p_end}
{phang2}{cmd:. local pessimistic "`r(pessvars)'"}{p_end}
{phang2}{cmd:. local optimistic "`r(optvars)'"}{p_end}
{phang2}{cmd:. assert "`pessimistic'"=="close_pess_a0_1 close_pess_a0_2 close_pess_a100"}{p_end}
{phang2}{cmd:. assert "`optimistic'"=="close_opt_a0_1 close_opt_a0_2 close_opt_a100"}{p_end}
{phang2}{cmd:. list dmu close_pess_a0_1 close_pess_a0_2, clean}{p_end}

    {title:Output-oriented (stored as 1/φ)}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) rts(crs) ort(output) gen(out_f) dmu(dmu)}{p_end}

    {title:Fuzzy inputs only (crisp outputs)}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) gen(fi) dmu(dmu)}{p_end}

    {title:Verification: α=1 should equal standard DEA}

{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) ort(i) gen(theta) dmu(dmu)}{p_end}
{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) alpha(1) gen(chk) dmu(dmu)}{p_end}
{phang2}{cmd:. list dmu theta chk_pess_a100 chk_opt_a100, clean}{p_end}
{phang2}{cmd:. assert abs(theta-chk_pess_a100) < 1e-6 & abs(theta-chk_opt_a100) < 1e-6}{p_end}

    {title:NIRS/NDRS aliases and crisp limits}

{phang2}{cmd:. foreach spec in nirs ndrs {c -(}}{p_end}
{phang2}{cmd:. {space 4}dea i_employee i_area = o_sales o_profits, rts(`spec') ort(output) engine(simplex) gen(c_`spec') dmu(dmu)}{p_end}
{phang2}{cmd:. {space 4}dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) rts(`spec') ort(output) gen(f_`spec') dmu(dmu)}{p_end}
{phang2}{cmd:. {space 4}assert abs(c_`spec'-f_`spec'_pess_a100)<1e-6 & abs(c_`spec'-f_`spec'_opt_a100)<1e-6}{p_end}
{phang2}{cmd:. {c )-}}{p_end}

    {title:Compute a finer membership grid}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) alpha(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1) gen(mf) dmu(dmu)}{p_end}
{phang2}{cmd:. list dmu mf_pess_a0 mf_opt_a0 mf_pess_a100 mf_opt_a100, clean}{p_end}



{pstd}Restrict both evaluation and reference DMUs{p_end}
{phang2}{cmd:. dea_fuzzy if i_employee>=15 in 2/5, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) alpha(0 1) dmu(dmu) gen(sub)}{p_end}
{phang2}{cmd:. assert r(N)==3}{p_end}
{phang2}{cmd:. assert missing(sub_pess_a100) if i_employee<15 | _n<2 | _n>5}{p_end}


{pstd}Documented minimum option abbreviations

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) dmu(dmu) alpha(1) rts(crs) orientation(output) generate(fullopt)}{p_end}
{phang2}{cmd:. dea_fuzzy, in(i_employee i_area) out(o_sales o_profits) dmu(dmu) alpha(1) r(crs) or(output) gen(shortopt)}{p_end}
{phang2}{cmd:. assert fullopt_pess_a100==shortopt_pess_a100}{p_end}

{pstd}Inspect full-precision endpoints and the common summary cohort{p_end}
{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) dmu(dmu) alpha(0 .37 1) gen(audit) matrix}{p_end}
{phang2}{cmd:. matrix list r(summary)}{p_end}
{phang2}{cmd:. matrix list r(eff_pess)}{p_end}
{phang2}{cmd:. assert audit_pess_a100==audit_opt_a100}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea_fuzzy} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of selected DMUs{p_end}
{synopt:{cmd:r(N_complete)}}DMUs with both endpoints at every requested cut{p_end}
{synopt:{cmd:r(N_partial)}}selected DMUs outside that complete cohort{p_end}
{synopt:{cmd:r(m_in)}}number of inputs{p_end}
{synopt:{cmd:r(m_out)}}number of outputs{p_end}
{synopt:{cmd:r(n_alpha)}}number of distinct alpha levels{p_end}
{synopt:{cmd:r(tolerance)}}requested initial convergence tolerance{p_end}
{synopt:{cmd:r(check_tolerance)}}fixed validation threshold, 1e-8{p_end}
{synopt:{cmd:r(lp_count)}}LP problems evaluated, including grid reruns{p_end}
{synopt:{cmd:r(lp_retries)}}extra single-LP attempts{p_end}
{synopt:{cmd:r(grid_retries)}}complete alpha-grid retries, 0 or 1{p_end}
{synopt:{cmd:r(lp_missing)}}missing individual scores in the final grid{p_end}
{synopt:{cmd:r(lp_infeasible)}}legacy name: repeated native status 1 reports across all passes, not a certificate{p_end}
{synopt:{cmd:r(lp_unbounded)}}legacy name: repeated native status 2 reports across all passes, not a certificate{p_end}
{synopt:{cmd:r(max_primal_violation)}}largest scaled residual among accepted solutions, including bound and objective-consistency checks{p_end}

{pstd}
The total number of solver attempts is {cmd:r(lp_count)+r(lp_retries)}.
Missing scores are retained; nesting checks compare all known levels across gaps.

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}}{cmd:dea_fuzzy}{p_end}
{synopt:{cmd:r(formulation)}}{cmd:triangular_kaoliu_basic}{p_end}
{synopt:{cmd:r(summary_sample)}}{cmd:complete_alpha_grid}{p_end}
{synopt:{cmd:r(orientation)}}orientation ({cmd:input} or {cmd:output}){p_end}
{synopt:{cmd:r(rts)}}returns to scale ({cmd:crs}, {cmd:vrs}, {cmd:drs}, or {cmd:irs}){p_end}
{synopt:{cmd:r(alpha)}}parsed α-levels used, in ascending order{p_end}
{synopt:{cmd:r(pessvars)}}pessimistic variable names, in the same order as {cmd:r(alpha)}{p_end}
{synopt:{cmd:r(optvars)}}optimistic variable names, in the same order as {cmd:r(alpha)}{p_end}
{synopt:{cmd:r(generate)}}generate prefix used{p_end}

{pstd}
Generated variables (prefix {it:P} = {cmd:r(generate)}, name component {it:N} = round(α×100)):

{synoptset 28 tabbed}{...}
{synopt:{it:P}{cmd:_pess_a}{it:N}}pessimistic efficiency at the specified α{p_end}
{synopt:{it:P}{cmd:_opt_a}{it:N}}optimistic efficiency at the specified α{p_end}
{synoptline}

{pstd}
Colliding groups append {cmd:_1}, {cmd:_2}, etc. to both names above.
Use the returned lists to identify the actual alpha corresponding to each
variable; a rounded suffix alone does not specify the calculation's alpha.

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:r(summary)}}one row per sorted alpha; columns alpha, N_pess, N_opt, N_pair, N_complete, mean_pess, mean_opt, min_pess, max_pess, min_opt, max_opt{p_end}
{synopt:{cmd:r(eff_pess)}}pessimistic scores; requires {cmd:matrix}{p_end}
{synopt:{cmd:r(eff_opt)}}optimistic scores; requires {cmd:matrix}{p_end}
{synopt:{cmd:r(complete)}}complete-cohort flag; requires {cmd:matrix}{p_end}
{synopt:{cmd:r(obs_index)}}original observation number; requires {cmd:matrix}{p_end}

{marker references}{...}
{title:References}

{pstd}
{it:Fuzzy DEA}

{marker KL2000}{...}
{phang}
Kao, C., and S.-T. Liu. 2000.
Fuzzy efficiency measures in data envelopment analysis.
{it:Fuzzy Sets and Systems} 113: 427-437.

{marker LFJ2003}{...}
{phang}
Lertworasirikul, S., S.-C. Fang, J. A. Joines, and H. L. W. Nuttle. 2003.
Fuzzy data envelopment analysis (DEA): a possibility approach.
{it:Fuzzy Sets and Systems} 139: 379-394.

{marker GT2001}{...}
{phang}
Guo, P., and H. Tanaka. 2001.
Fuzzy DEA: a perceptual evaluation method.
{it:Fuzzy Sets and Systems} 119: 149-160.

{pstd}
{it:Imprecise DEA (related)}

{marker CPY1999}{...}
{phang}
Cooper, W. W., K. S. Park, and G. Yu. 1999.
IDEA and AR-IDEA: Models for dealing with imprecise data in DEA.
{it:Management Science} 45(4): 597-607.
{browse "https://doi.org/10.1287/mnsc.45.4.597":doi:10.1287/mnsc.45.4.597}

{pstd}
{it:Foundational DEA}

{marker CCR1978}{...}
{phang}
Charnes, A., W. W. Cooper, and E. Rhodes. 1978.
Measuring the efficiency of decision making units.
{it:European Journal of Operational Research} 2: 429-444.

{marker BCC1984}{...}
{phang}
Banker, R. D., A. Charnes, and W. W. Cooper. 1984.
Some models for estimating technical and scale inefficiencies in data
envelopment analysis.
{it:Management Science} 30: 1078-1092.

{pstd}
{it:Software}

{marker JL2010}{...}
{phang}
Ji, Y.-B., and C. Lee. 2010.
Data envelopment analysis.
{it:Stata Journal} 10(2): 267-280.


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
{synopt:{helpb dea_fuzzy}}Fuzzy DEA{p_end}
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
