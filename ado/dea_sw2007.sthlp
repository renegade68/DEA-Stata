{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.1.0  01may2026}{...}
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

{syntab:Backend}
{synopt:{cmdab:back:end(}{it:name}{cmd:)}}{cmd:internal} (default) or {cmd:simarwilson}{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for result variables; default is {cmd:sw07}{p_end}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable (string or numeric){p_end}
{synopt:{cmdab:nobc:orrected}}suppress bias-correction display{p_end}

{syntab:Technical}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmd:trace}}display per-iteration progress{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:rtstype}: {cmd:crs} | {cmd:vrs} | {cmd:drs} | {cmd:irs};
{cmd:nirs} is available with {cmd:backend(simarwilson)} only{p_end}
{p 4 6 2}
{it:type}: {cmd:input} | {cmd:output}{p_end}
{p 4 6 2}
{it:name}: {cmd:internal} | {cmd:simarwilson}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_sw2007} implements the two-stage data envelopment analysis estimator
of {help dea_sw2007##sw2007:Simar and Wilson (2007), Algorithm 2}.  The
procedure (i) computes DEA efficiency scores in stage one, (ii) regresses
inverse-efficiency on a vector of environmental variables Z using truncated
maximum likelihood, and (iii) returns bootstrap bias-corrected coefficient
estimates and percentile confidence intervals that are consistent with the
underlying data-generating process.

{pstd}
The command provides two interchangeable backends.  The default
{cmd:backend(internal)} uses a self-contained Mata implementation that solves
the DEA linear programs through {cmd:_dea_simplex.ado} and supports panel
data and the full set of returns-to-scale options
({cmd:crs}/{cmd:vrs}/{cmd:drs}/{cmd:irs}).  The optional
{cmd:backend(simarwilson)} dispatches to the community-contributed
{help dea_sw2007##bt2019:Badunenko and Tauchmann (2019)} {cmd:simarwilson}
package and is intended primarily as a verification path: identical inputs
produce numerically equivalent estimates within Monte Carlo error.


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
within-period reference set, matching Simar and Wilson (2007, Section 3.1).
Panel data are not supported with {cmd:backend(simarwilson)}.

{dlgtab:Bootstrap}

{phang}
{cmd:reps(}{it:#}{cmd:)} sets the number of bootstrap replications used to
construct percentile confidence intervals (B in Simar and Wilson 2007).
Default is {cmd:reps(2000)}.

{phang}
{cmd:bcreps(}{it:#}{cmd:)} sets the number of bootstrap replications for
bias correction of stage-one DEA scores (B1 in Simar and Wilson 2007).
Default is {cmd:bcreps(1000)}.

{phang}
{cmd:seed(}{it:#}{cmd:)} sets the random-number seed for reproducible
results.  Setting the same seed across backends does not in general
produce bit-identical estimates because internal and simarwilson call the
RNG at different points; statistical equivalence holds.

{dlgtab:Backend}

{phang}
{cmdab:back:end(}{it:name}{cmd:)} selects the computational engine.
Default is {cmd:backend(internal)}.

{phang2}
{cmd:internal} uses the package's own Mata implementation and solves the DEA
linear programs through {cmd:_dea_simplex.ado}.  All five {cmd:rts()}
settings and panel data are supported.

{phang2}
{cmd:simarwilson} dispatches to the community-contributed {cmd:simarwilson}
command (Badunenko and Tauchmann 2019, {it:Stata Journal} 19(4)).  Requires
{cmd:simarwilson} to be installed
({stata ssc install simarwilson:ssc install simarwilson}).
Cross-sectional data only; {cmd:rts()} restricted to
{cmd:crs}/{cmd:vrs}/{cmd:nirs}.  The score-only {cmd:scores()} syntax is
also available through this backend as Algorithm 1 inference.  When dispatch succeeds, the saved results
follow simarwilson's conventions
(see {help dea_sw2007##stored_results:Stored results} below).

{dlgtab:Output}

{phang}
{cmdab:gen:erate(}{it:prefix}{cmd:)} specifies a prefix for the generated
result variables.  Default prefix is {cmd:sw07}.

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies a variable identifying each DMU.
Accepts both string and numeric variables.


{marker remarks}{...}
{title:Remarks}

{pstd}
The Simar and Wilson (2007) Algorithm 2 procedure is the appropriate
inferential tool when the goal is to test hypotheses about the effect of
environmental variables on technical efficiency.  Conventional two-step
approaches (DEA followed by Tobit or OLS) are flagged by Simar and Wilson
(2007) as producing invalid standard errors due to (a) the bounded support
of DEA scores and (b) the cross-DMU dependence inherent in DEA.

{pstd}
{cmd:backend(simarwilson)} is provided for verification rather than for
routine use.  Running the same model under both backends and comparing the
estimated coefficients, sigma, and percentile CI bounds is a useful sanity
check.  Numerical equivalence within ~1% relative tolerance is expected;
exact agreement is precluded by differences in random-number streams.

{pstd}
For panel data ({cmd:panel()} {cmd:time()}), only {cmd:backend(internal)}
is available because {cmd:simarwilson} is documented as cross-sectional
only (Badunenko and Tauchmann 2019, footnote 10).  In the panel mode,
{cmd:dea_sw2007} solves DEA separately within each period, then runs the
truncated regression on the pooled inverse-efficiency series with Z, and
finally bootstraps the entire two-stage procedure.


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
{synopt:{cmd:e(sigma)}}initial sigma-hat{p_end}

{syntab:Macros}
{synopt:{cmd:e(cmd)}}{cmd:dea_sw2007}{p_end}
{synopt:{cmd:e(backend)}}{cmd:internal}{p_end}
{synopt:{cmd:e(rts)}}returns-to-scale specification{p_end}
{synopt:{cmd:e(orientation)}}{cmd:input} or {cmd:output}{p_end}

{syntab:Matrices}
{synopt:{cmd:e(b)}}coefficient vector (initial truncated regression){p_end}
{synopt:{cmd:e(beta_bc)}}bias-corrected coefficient vector{p_end}
{synopt:{cmd:e(ci_lo)}}lower CI bounds{p_end}
{synopt:{cmd:e(ci_hi)}}upper CI bounds{p_end}
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
bloom.rampike@gmail.com; bloom.rampike@gmail.com

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
