{smcl}
{* Help revision 2026-11-02: same-sample MLE retries and inference assumptions.}{...}
{* Help revision 2026-11-02: observation-file saving for internal and external backends.}{...}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* *! version 1.2.2  02nov2026}{...}
{viewerjumpto "Related article" "dea_sw2007##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_bootstrap" "help dea_bootstrap"}{...}
{vieweralsosee "dea_crosseff" "help dea_crosseff"}{...}
{vieweralsosee "dea_malmquist" "help dea_malmquist"}{...}
{viewerjumpto "Syntax" "dea_sw2007##syntax"}{...}
{viewerjumpto "Description" "dea_sw2007##description"}{...}
{viewerjumpto "Options" "dea_sw2007##options"}{...}
{viewerjumpto "Remarks" "dea_sw2007##remarks"}{...}
{viewerjumpto "Examples" "dea_sw2007##examples"}{...}
{viewerjumpto "Stored results" "dea_sw2007##stored_results"}{...}
{viewerjumpto "References" "dea_sw2007##references"}{...}
{viewerjumpto "Authors" "dea_sw2007##authors"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:dea_sw2007} {hline 2}}Two-stage DEA with environmental variables (Simar and Wilson 2007 Algorithm 2){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea_sw2007}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}{cmd:,}
{cmdab:env:ars(}{it:varlist}{cmd:)}
[{it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_sw2007}
{ifin}{cmd:,}
{cmdab:in:puts(}{it:varlist}{cmd:)}
{cmdab:out:puts(}{it:varlist}{cmd:)}
{cmdab:env:ars(}{it:varlist}{cmd:)}
[{it:options}]


{pstd}
Form C (pre-computed scores)

{p 8 16 2}
{cmd:dea_sw2007}
{ifin}{cmd:,}
{cmdab:sc:ores(}{it:varname}{cmd:)}
{cmdab:env:ars(}{it:varlist}{cmd:)}
[{it:options}]


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmdab:env:ars(}{it:varlist}{cmd:)}}environmental (Z) variables for second-stage regression{p_end}

{syntab:DEA model (Forms A and B)}
{synopt:{cmd:rts(}{it:rtstype}{cmd:)}}returns to scale; default is {cmd:rts(vrs)}{p_end}
{synopt:{cmdab:ori:entation(}{it:type}{cmd:)}}orientation; default is {cmd:orientation(output)}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}

{syntab:Panel}
{synopt:{cmdab:pa:nel(}{it:varname}{cmd:)}}panel identifier (cross-section if omitted){p_end}
{synopt:{cmdab:ti:me(}{it:varname}{cmd:)}}time identifier (required with {cmd:panel()}){p_end}

{syntab:Bootstrap}
{synopt:{cmd:reps(}{it:#}{cmd:)}}bootstrap replications for inference; default is {cmd:reps(2000)}{p_end}
{synopt:{cmd:bcreps(}{it:#}{cmd:)}}bootstrap replications for bias correction; default is {cmd:bcreps(1000)}{p_end}
{synopt:{cmd:level(}{it:#}{cmd:)}}confidence level; default is {cmd:level(95)}{p_end}
{synopt:{cmd:seed(}{it:#}{cmd:)}}random-number seed for reproducibility{p_end}
{synopt:{cmd:tau(}{it:#}{cmd:)}}truncation point; default is {cmd:tau(1.0)}{p_end}

{syntab:Procedure and backend}
{synopt:{cmdab:proc:edure(}{it:name}{cmd:)}}{cmd:sequential} (default, Algorithm 2) or {cmd:joint} (legacy){p_end}
{synopt:{cmdab:back:end(}{it:name}{cmd:)}}{cmd:internal} (default) or {cmd:simarwilson}{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for result variables; default is {cmd:sw07}{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable (string or numeric){p_end}
{synopt:{cmdab:sav:ing(}{it:filename}{cmd:)}}save selected observations and score results as a Stata dataset{p_end}
{synopt:{cmd:replace}}allow file replacement; also retains the legacy alias for {cmd:replacevars}{p_end}
{synopt:{cmdab:nobc:orrected}}suppress corrected-score output; see procedure-specific behavior below{p_end}

{syntab:Technical}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmd:trace}}display per-iteration progress{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:rtstype}: {cmd:crs} | {cmd:vrs} | {cmd:drs} | {cmd:irs};
{cmd:nirs} is an alias for {cmd:drs} internally{p_end}
{p 4 6 2}
{it:type}: {cmd:input} | {cmd:output}{p_end}
{p 4 6 2}
{it:name}: {cmd:internal} | {cmd:simarwilson}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_sw2007} provides sequential Simar-Wilson (2007) Algorithm 2 and a
separately selectable legacy estimator. The default {cmd:procedure(sequential)}
computes original DEA scores, fits an initial truncated regression, performs
{cmd:bcreps()} pseudo-frontier DEA replications, corrects the inverse scores,
refits the regression, and performs {cmd:reps()} parametric regression replications.
The second bootstrap does not re-solve DEA.

{pstd}
Both orientations store efficiency theta in (0,1]. The regression uses
delta = tau/theta, with left truncation at tau (default 1). Score correction is
performed as delta_bc = 2*delta - mean(delta_boot), and the generated corrected
efficiency is tau/delta_bc. This is generally different from correcting theta
directly. {cmd:e(b)} contains the regression fitted to delta_bc; {cmd:e(V)}
is the covariance of the second bootstrap coefficient and sigma estimates.

{pstd}
{cmd:backend(internal)} uses the package simplex solver. The optional
{cmd:backend(simarwilson)} calls Badunenko and Tauchmann's separately installed
command. These are comparison paths with explicit differences in boundary
handling and percentile interpolation; numerical identity is not promised.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{cmdab:env:ars(}{it:varlist}{cmd:)} lists the environmental variables Z that
enter the second-stage truncated regression.  At least one Z variable is
required.

{dlgtab:DEA model}

{phang}
{cmd:rts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption
for the stage-one DEA.  Default is {cmd:rts(vrs)}, which matches Simar and
Wilson (2007).  With {cmd:backend(simarwilson)}, only
{cmd:crs}/{cmd:vrs}/{cmd:nirs} are accepted.

{phang}
{cmdab:ori:entation(}{it:type}{cmd:)} selects input or output orientation
for the Farrell distance measure.  Default is {cmd:orientation(output)}.

{dlgtab:Panel}

{phang}
{cmdab:pa:nel(}{it:varname}{cmd:)} {cmd:time(}{it:varname}{cmd:)} declare
panel data.  When both are specified, DEA is solved year-by-year with the
within-period reference set. Each selected period must have at least 3 DMUs.
The regression pools periods and draws independently across observations; this
is a package extension, not a cluster bootstrap for repeated DMUs. Serially
correlated panel inference is not validated. Panel data are unavailable with
{cmd:backend(simarwilson)}.

{dlgtab:Bootstrap}

{phang}
{cmd:reps(}{it:#}{cmd:)} sets the number of bootstrap replications used to
construct confidence intervals (B2 in Algorithm 2).
Default is {cmd:reps(2000)}; at least 2 are required.

{phang}
{cmd:bcreps(}{it:#}{cmd:)} sets the number of bootstrap replications for
bias correction of stage-one DEA scores (B1 in Simar and Wilson 2007).
Default is {cmd:bcreps(1000)}; at least 2 are required.

{phang}
{cmd:seed(}{it:#}{cmd:)} sets the random-number seed for reproducible
results.  Setting the same seed across backends does not in general
produce bit-identical estimates because internal and simarwilson call the
RNG at different points. Matching model, score scale, boundary classification,
and estimation sample is necessary before comparing simulation results.

{dlgtab:Procedure and backend}

{phang}
{cmd:procedure(sequential)} is the default and follows Algorithm 2's separate
bootstrap loops and corrected-score regression refit. Both loops require the
requested number of successful replications; they stop with an error if this
cannot be achieved in twice that many attempts. A DEA replication is accepted
only when every selected DMU has a valid score. Failed or nonconverged regression
replications are not included. The initial regression-bootstrap iteration limit
is max(25,3*iterations of the corrected-score fit). If that limit is below 200
and the fit fails, the command retries the same pseudo-data once with 200
iterations before rejecting the replication. It draws no new errors for this
retry and retains an initially valid fit. Recovered fits must converge, contain
all selected observations and parameters, and have finite estimates.
{cmd:e(boot_retried)} and {cmd:e(boot_recovered)} count these extra fits and
recoveries; {cmd:e(boot_failed)} counts rejected pseudo-samples after this check.

{phang}
From version 1.2.2, a recovered regression fit can change the bootstrap covariance
and intervals, and the subsequent RNG state. Original and corrected DEA scores
and the corrected-score point regression are unchanged for the same data and
seed. Recompute earlier sequential inference if regression failures were
reported. The historical joint procedure and the external backend retain their
own fitting rules.

{phang}
{cmd:procedure(joint)} preserves the previous implementation for reproduction.
It uses max(reps,bcreps) common replications, re-solves DEA and regresses its
inverse scores in each replication, corrects theta afterward, and returns the
initial regression in {cmd:e(b)} and its model-based covariance in {cmd:e(V)}.
It does not refit on corrected scores before inference. It is a legacy
estimator, not Simar-Wilson Algorithm 2; its statistical validity is not
established by the numerical reproduction tests. It is available internally only.


{phang}
{cmdab:back:end(}{it:name}{cmd:)} selects the computational engine.
Default is {cmd:backend(internal)}.

{phang2}
{cmd:internal} uses the package's own Mata implementation and solves the DEA
linear programs through {cmd:_dea_simplex.ado}.  The four technology restrictions and the {cmd:nirs} alias
are supported, as is the period-specific panel extension.

{phang2}
{cmd:simarwilson} dispatches to the community-contributed {cmd:simarwilson}
command (Badunenko and Tauchmann 2019, {it:Stata Journal} 19(4)).  Requires
{cmd:simarwilson} to be installed
({stata ssc install simarwilson:ssc install simarwilson}).
Cross-sectional data only; {cmd:rts()} restricted to
{cmd:crs}/{cmd:vrs}/{cmd:nirs}.  The score-only {cmd:scores()} syntax is
also available through this backend as Algorithm 1 inference. Supply theta in
(0,1]; the wrapper converts it to 1/theta and passes {cmd:nounit}. With inputs and
outputs but no {cmd:scores()}, it calls Algorithm 2. For input orientation it
passes {cmd:invert}, so both orientations use inverse efficiency above 1.
This backend requires {cmd:tau(1)} and {cmd:procedure(sequential)}.  When dispatch succeeds, the saved results
follow simarwilson's conventions
(see {help dea_sw2007##stored_results:Stored results} below).

{dlgtab:Output}

{phang}
{cmdab:gen:erate(}{it:prefix}{cmd:)} specifies a prefix for the generated
result variables.  Default prefix is {cmd:sw07}.

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies a variable identifying each DMU.
Accepts both string and numeric variables.

{phang}
{cmdab:sav:ing(}{it:filename}{cmd:)} writes a Stata dataset after successful
estimation. Quote paths containing spaces. The extension is {cmd:.dta} if
none is supplied. The file contains observations selected by {cmd:if}/{cmd:in}
after excluding missing values in the supplied model variables and panel/time
keys. This is the wrapper's analysis sample, which can differ from the
external backend's regression {cmd:e(sample)}. Original observation order,
variable types and labels are retained. No identifier is invented if
{cmd:dmu()} is omitted.

{pstd}
Columns are the supplied {cmd:dmu()}, {cmd:panel()}, {cmd:time()}, input,
output, {cmd:scores()} and environmental variables (duplicates removed),
followed by the result columns. Unrelated variables and internal temporary
variables are excluded. With {cmd:backend(internal)}, both procedures save
{it:prefix}{cmd:_eff} and, unless {cmd:nobcorrected} is specified,
{it:prefix}{cmd:_eff_bc}, exactly as generated in memory.

{pstd}
With {cmd:backend(simarwilson)}, Algorithm 2 saves the public implementation's
native inverse-efficiency measures as {it:prefix}{cmd:_delta} and
{it:prefix}{cmd:_delta_bc}; {cmd:nobcorrected} omits the latter from the file
but does not change estimation. These use the inverse-score scale (normally
at least 1), unlike the internal efficiency columns. No reciprocation,
clipping or replacement of missing values is applied. With {cmd:scores()},
the external backend uses Algorithm 1 and saves only {it:prefix}{cmd:_delta}
equal to {cmd:1/scores()}; Algorithm 1 does not produce a corrected DEA score.
The external result columns exist only in the saved file. The requested
prefix must be valid and must not conflict with exported model variables.

{pstd}
Dataset characteristics {cmd:dea_sw2007_backend}, {cmd:dea_sw2007_procedure},
{cmd:dea_sw2007_score_scale} and {cmd:dea_sw2007_sample} describe the file.
The observation file does not store regression matrices or bootstrap draws;
use {cmd:estimates save} to retain the estimation results. Saving preserves
the estimation results, regression equation names and random-number state.

{pstd}
An existing file requires {cmd:replace}; {cmd:replacevars} alone never permits
file replacement. For compatibility, {cmd:replace} also authorizes replacing
the internal generated result variables, subject to model-data protection.
File/name preflight occurs before seeding and estimation. When {cmd:saving()}
is requested, an estimation or write error restores the caller's data;
a final write error also clears the new {cmd:e()} results. An I/O error is
reported as an error, not as a successfully saved analysis.


{marker remarks}{...}
{title:Remarks}

{pstd}
The Simar and Wilson (2007) Algorithm 2 procedure supports inference about
environmental variables under its data-generating assumptions: a common
technology separable from the environmental variables, a linear conditional
inverse-efficiency model with truncated-normal disturbances, constant sigma,
and independent disturbances across DMUs. Conventional two-step
approaches (DEA followed by Tobit or OLS) are flagged by Simar and Wilson
(2007) as producing invalid standard errors due to (a) the bounded support
of DEA scores and (b) the cross-DMU dependence inherent in DEA.

{pstd}
These assumptions distinguish dependence among estimated DEA scores, induced
by their common estimated frontier, from dependence among the underlying
disturbances. The independent parametric draws address the former; they do not
implement a cluster bootstrap or a heteroskedastic regression model. Increasing
{cmd:bcreps()} or {cmd:reps()} cannot repair such misspecification. In a small,
multidimensional sample, finite-sample bias and confidence-interval coverage
still require sensitivity analysis even when the numerical steps succeed.

{pstd}
{cmd:backend(simarwilson)} provides an independent implementation for
comparison. The published Mac plugin examined here is Intel x86_64; on Apple
Silicon it runs with the Intel version of Stata under Rosetta. An unrecognized
plugin command in a native ARM session is an execution-environment problem.

{pstd}
The internal sequential path snaps original scores within {cmd:tolerance()}
of 1 to exactly 1 before excluding efficient DMUs from the initial regression.
The published reference may classify values such as 1+2e-16 differently.
A fixed-seed VRS input example exhibited this difference despite matching LP
scores to machine precision; further examples also occurred under CRS and VRS
output orientation. Compare initial regression samples as well as
scores. The tests retain the unmodified public results and separately document
a one-line boundary-alignment change in an isolated reference copy.

{pstd}
Even after boundary alignment, tiny changes in pseudo-data can change whether
{cmd:truncreg} converges within a short iteration limit. Consequently, the two
implementations can accept different bootstrap draws. The internal same-data
retry above is another explicit difference from the public 2.22 implementation.
Compare accepted-draw counts and convergence as well as coefficient estimates.

{pstd}
Sequential percentile intervals use the empirical inverse CDF (ceil(B*p));
{cmd:e(ci_pivot)} also returns basic/pivot intervals. The public command uses
interpolated percentiles, so endpoints need not agree even for identical draws.
The joint path retains its historical pivot quantiles. Numerical LP agreement
and finite example comparisons do not establish nominal coverage across data
generating processes.

{pstd}
{cmd:nobcorrected} with the sequential procedure suppresses only the
{cmd:generate()}_eff_bc variable; score correction and the refit still run.
With {cmd:procedure(joint)}, it retains the old behavior of skipping score
correction accumulation. Sequential estimation requires strictly positive
inputs and outputs. All supplied {cmd:scores()} must be in (0,1]. Invalid or
missing computed scores cause an error instead of being recoded as efficient.
{cmd:tau()} scales the inverse response and truncation point together; it does
not move the efficiency frontier away from theta=1.

{pstd}
{bf:Protecting existing data.} All output names are checked before results are
created. An existing result variable causes an error unless {cmd:replacevars}
is specified. The option applies only to variables produced by the current
command; it does not remove other variables sharing the prefix. Replacement
recreates whole result columns, including their storage types and labels;
observations not scored by the command do not retain old scores. Variables
used as model data or identifiers cannot be replaced, even with
{cmd:replacevars}. Invalid or duplicate output names are rejected before any
existing result column is changed. The name check itself does not promise
rollback of later computation errors; calls with {cmd:saving()} additionally
preserve the caller's data until successful completion.

{pstd}
For {cmd:dea_sw2007}, the previously supported {cmd:replace} spelling remains
an alias for {cmd:replacevars}, and permits overwriting the requested
{cmd:saving()} file. The same model-variable protection applies.
The external {cmd:backend(simarwilson)} does not create or replace the internal
{cmd:generate()} result columns.

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

{pstd}Synthetic teaching data. The small replication counts keep the example
short; increase both counts for substantive inference. The internal backend
requires no community-contributed package.{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. set obs 40}{p_end}
{phang2}{cmd:. generate long country = _n}{p_end}
{phang2}{cmd:. generate double exp_k = 5 + 10*runiform()}{p_end}
{phang2}{cmd:. generate double ptr = 10 + 10*runiform()}{p_end}
{phang2}{cmd:. generate double env_exp = runiform()}{p_end}
{phang2}{cmd:. generate double pisa_math = (exp_k^.6 * ptr^.4) / (1 + .3*env_exp + abs(rnormal(0,.2)))}{p_end}
{phang2}{cmd:. generate double pisa_read = pisa_math*1.1}{p_end}
{phang2}{cmd:. generate double pisa_sci = pisa_math*.9}{p_end}
{phang2}{cmd:. dea_sw2007 exp_k ptr = pisa_math, envars(env_exp) rts(vrs) ort(out) reps(50) bcreps(50) seed(12345) generate(swint)}{p_end}
{phang2}{cmd:. expand 2}{p_end}
{phang2}{cmd:. bysort country: generate year = 2020 + _n}{p_end}
{phang2}{cmd:. replace pisa_math = pisa_math * 1.01 if year == 2022}{p_end}
{phang2}{cmd:. dea_sw2007 exp_k ptr = pisa_math, envars(env_exp) panel(country) time(year) rts(vrs) reps(50) bcreps(50) seed(12345) generate(swpan)}{p_end}

{pstd}Pre-computed scores on a fresh cross-section:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. set obs 40}{p_end}
{phang2}{cmd:. generate long country = _n}{p_end}
{phang2}{cmd:. generate double exp_k = 5 + 10*runiform()}{p_end}
{phang2}{cmd:. generate double ptr = 10 + 10*runiform()}{p_end}
{phang2}{cmd:. generate double env_exp = runiform()}{p_end}
{phang2}{cmd:. generate double pisa_math = (exp_k^.6 * ptr^.4) / (1 + .3*env_exp + abs(rnormal(0,.2)))}{p_end}
{phang2}{cmd:. generate double pisa_read = pisa_math*1.1}{p_end}
{phang2}{cmd:. generate double pisa_sci = pisa_math*.9}{p_end}
{phang2}{cmd:. dea exp_k ptr = pisa_math, rts(vrs) ort(out) engine(simplex) gen(theta)}{p_end}
{phang2}{cmd:. dea_sw2007 pisa_math, scores(theta) inputs(exp_k ptr) envars(env_exp) reps(50) bcreps(50) seed(12345) generate(swint2)}{p_end}

{pstd}Legacy reproduction is an explicit secondary option:{p_end}
{phang2}{cmd:. dea_sw2007 exp_k ptr = pisa_math, envars(env_exp) rts(vrs) reps(50) bcreps(50) seed(12345) generate(swold) procedure(joint)}{p_end}

{pstd}Save observation results while retaining the estimation results:{p_end}
{phang2}{cmd:. tempfile sw_results}{p_end}
{phang2}{cmd:. dea_sw2007 exp_k ptr = pisa_math, envars(env_exp) dmu(country) rts(vrs) reps(50) bcreps(50) seed(12345) generate(swsave) saving("`sw_results'")}{p_end}
{phang2}{cmd:. confirm file "`sw_results'"}{p_end}
{phang2}{cmd:. assert "`e(cmd)'"=="dea_sw2007"}{p_end}
{phang2}{cmd:. matrix list e(b)}{p_end}

{pstd}The optional {cmd:backend(simarwilson)} requires that separately installed package; it is not used by these examples.{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
With {cmd:backend(internal)}, {cmd:dea_sw2007} stores the following in
{cmd:e()}:

{synoptset 22 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}}number of DMUs{p_end}
{synopt:{cmd:e(N_trunc)}}DMUs entering the truncated regression{p_end}
{synopt:{cmd:e(reps)}}bootstrap replications for CI{p_end}
{synopt:{cmd:e(bcreps)}}bootstrap replications for bias correction{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}
{synopt:{cmd:e(sigma)}}sigma from the corrected-score fit (initial fit for joint){p_end}
{synopt:{cmd:e(sigma0)}}initial sigma{p_end}
{synopt:{cmd:e(N_initial)}}initial regression sample size (sequential){p_end}
{synopt:{cmd:e(bc_success)}}successful DEA replications (sequential){p_end}
{synopt:{cmd:e(boot_success)}}successful regression replications{p_end}
{synopt:{cmd:e(bc_failed)}}failed DEA attempts (sequential){p_end}
{synopt:{cmd:e(boot_failed)}}rejected regression pseudo-samples after retries (sequential){p_end}
{synopt:{cmd:e(boot_fit_attempts)}}regression fits, including same-data retries (sequential){p_end}
{synopt:{cmd:e(boot_retried)}}same-data regression retries (sequential){p_end}
{synopt:{cmd:e(boot_recovered)}}valid fits recovered by those retries (sequential){p_end}

{syntab:Macros}
{synopt:{cmd:e(cmd)}}{cmd:dea_sw2007}{p_end}
{synopt:{cmd:e(backend)}}{cmd:internal}{p_end}
{synopt:{cmd:e(rts)}}returns-to-scale specification{p_end}
{synopt:{cmd:e(orientation)}}{cmd:input} or {cmd:output}{p_end}

{syntab:Matrices}
{synopt:{cmd:e(b)}}corrected-score regression coefficients and sigma; initial fit for joint{p_end}
{synopt:{cmd:e(V)}}bootstrap covariance; initial model-based covariance for joint{p_end}
{synopt:{cmd:e(b_initial)}}initial fit (sequential){p_end}
{synopt:{cmd:e(V_initial)}}initial model-based covariance (sequential){p_end}
{synopt:{cmd:e(V_modelbased)}}corrected-score model-based covariance (sequential){p_end}
{synopt:{cmd:e(b_reps)}}all B2 coefficient/sigma draws, in draw order (sequential){p_end}
{synopt:{cmd:e(ci_percentile)}}2 x parameter-count percentile limits (sequential){p_end}
{synopt:{cmd:e(ci_pivot)}}2 x parameter-count pivot limits (sequential){p_end}
{synopt:{cmd:e(procedure)}}{cmd:sequential} or {cmd:joint}{p_end}
{synopt:{cmd:e(beta_bc)}}additional coefficient bias correction, 2*point - bootstrap mean; not e(b){p_end}
{synopt:{cmd:e(ci_lo)}}coefficient percentile lower limits (pivot for joint){p_end}
{synopt:{cmd:e(ci_hi)}}coefficient percentile upper limits (pivot for joint){p_end}
{p2colreset}{...}

{pstd}
With {cmd:backend(simarwilson)}, the saved results follow the conventions
of {cmd:simarwilson} (see {help simarwilson} after installation), with
two additional macros:

{synoptset 22 tabbed}{...}
{synopt:{cmd:e(sw07_backend)}}{cmd:simarwilson}{p_end}
{synopt:{cmd:e(sw07_wrapper)}}{cmd:dea_sw2007}{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker bt2019}{...}
{phang}
Badunenko, O., and H. Tauchmann. 2019.
Simar and Wilson two-stage efficiency analysis for Stata.
{it:Stata Journal} 19(4): 950--988.

{marker sw1998}{...}
{phang}
Simar, L., and P. W. Wilson. 1998.
Sensitivity analysis of efficiency scores: How to bootstrap in nonparametric
frontier models.
{it:Management Science} 44: 49--61.

{marker sw2007}{...}
{phang}
Simar, L., and P. W. Wilson. 2007.
Estimation and inference in two-stage, semi-parametric models of production
processes.
{it:Journal of Econometrics} 136: 31--64.

{phang}
McDonald, J. 2009.
Using least squares and tobit in second stage DEA efficiency analyses.
{it:European Journal of Operational Research} 197: 792--798.


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
