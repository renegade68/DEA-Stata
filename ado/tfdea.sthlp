{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 2.1.1 20aug2026}{...}
{vieweralsosee "dea" "help dea"}{...}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:tfdea} {hline 2}}Technology forecasting using data envelopment analysis{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:tfdea} {it:inputvars} {cmd:=} {it:outputvars} {ifin}{cmd:,}
{cmd:tk(}{it:timevar}{cmd:)} {cmd:dmu(}{it:varname}{cmd:)}
[{cmd:rts(}{it:string}{cmd:)} {cmd:ort(}{it:string}{cmd:)}
{cmd:tf(}{it:#}{cmd:)} {cmd:time(}{it:#}{cmd:)}
{cmd:rocfrontier(final|first)} {cmd:edate(weighted|min|max)}
{cmd:continuation} {cmd:detail}]

{marker description}{...}
{title:Description}

{pstd}
{cmd:tfdea} estimates the historical annual rate of technological change from
sequential DEA frontiers.  A unit that was on the frontier when released and
is later superseded contributes an annual rate
gamma = theta^(1/(t_eff-t_k)).  The command reports the mean of valid rates.
When {cmd:tf()} leaves later observations outside the training frontier, it
also evaluates those observations and reports forecast effective dates.

{pstd}
The implementation follows the TFDEA formulation used in the third-edition
reproduction and requires Stata 16 or newer.  It calls the packaged
{helpb dea} command and Mata's {cmd:LinearProgram()} class.

{marker options}{...}
{title:Options}

{phang}{cmd:tk(}{it:timevar}{cmd:)} is required and identifies the numeric
release-period variable.

{phang}{cmd:dmu(}{it:varname}{cmd:)} is required and identifies observations.
String identifiers, including values containing spaces, are accepted.

{phang}{cmd:rts(vrs|bcc|crs|ccr)} specifies returns to scale.  The default is
{cmd:vrs}; {cmd:bcc} and {cmd:ccr} are accepted aliases.

{phang}{cmd:ort(output|out|input|in)} specifies orientation.  The default is
{cmd:output}.

{phang}{cmd:tf(}{it:#}{cmd:)} sets the final training-frontier period.  If it
is omitted, the last period in the estimation sample is used and no
out-of-training-sample forecast is produced.

{phang}{cmd:time(}{it:#}{cmd:)} records a nonnegative forecast horizon.  The
default is zero.

{phang}{cmd:rocfrontier(final|first)} selects the sequential frontier used for
ROC calculation. {cmd:final}, the default, follows the archived original R
TFDEA implementation: every DMU efficient at release but inefficient at the
final training frontier is evaluated against that final frontier. {cmd:first}
uses the first sequential frontier on which the DMU becomes inefficient and
preserves the previous Stata behavior.

{phang}{cmd:edate(weighted|min|max)} selects the effective date when alternate
optimal lambda vectors yield the same efficiency score. {cmd:weighted}, the
previous Stata behavior, uses the lambda vector returned by the primary
efficiency LP. {cmd:min}, the default, or {cmd:max} runs a
lexicographic second LP: the optimal efficiency is held fixed while the
lambda-weighted reference date is minimized or maximized. Under VRS this is
the minimum or maximum effective date directly. {cmd:min} corresponds to the
{cmd:second="min"} choice in the archived R TFDEA example.

{phang}{cmd:detail} displays intermediate frontier iterations and diagnostic
matrices in addition to the result matrices.

{phang}{cmd:continuation} enables an experimental fallback when a radial VRS
forecasting LP has no solution. It retains the VRS convexity constraint and
selects training-frontier weights that minimize the sum of normalized input
excesses and output shortfalls. If this minimum frontier gap is {it:g}, the
fallback forecasting index is {cmd:1+g} for input orientation and
{cmd:1/(1+g)} for output orientation. The ordinary TFDEA forecasting equation
then converts {cmd:ln(1+g)} to a time increment using the estimated mean ROC.
This continuation index is not a BCC super-efficiency score. Feasible BCC
results are unchanged, and the option has no fallback effect under CRS.

{marker continuation}{...}
{title:Experimental VRS continuation}

{pstd}
Xue and Harker (2002) show that infeasibility in an input-oriented VRS
super-efficiency model can carry substantive information: under their stated
assumptions, an infeasible DMU lies beyond the finite radial reach of the
reduced VRS frontier and may be assigned the highest efficiency class. Their
paper resolves the ranking interpretation; it does not construct finite
super-efficiency scores, reference weights, or TFDEA forecast dates for such
DMUs.

{pstd}
The {cmd:continuation} option is one experimental operational proposal
motivated by that interpretation. It is not a generally established
always-feasible DEA model, and it is not a model proposed by Xue and Harker.
When the ordinary radial VRS forecasting problem has no solution, the option
solves the following normalized additive continuation problem conceptually:

{p 8 12 2}
minimize  sum_m d_xm/a_xm + sum_s d_ys/a_ys

{p 8 12 2}
subject to  X lambda <= x_0 + d_x,
Y lambda >= y_0 - d_y,
sum_j lambda_j = 1,
and lambda, d_x, d_y >= 0.

{pstd}
The scale factors {it:a_xm} and {it:a_ys} are the larger of the absolute
target value, the corresponding mean absolute training value, and a small
positive numerical floor. This makes input excesses and output shortfalls
comparable across differently scaled variables. Among alternate solutions
with the same minimum gap, {cmd:edate(min|max)} applies the same secondary
effective-date rule used by the ordinary forecasting LP.

{pstd}
Let {it:g} denote the minimized normalized gap. For compatibility with the
existing TFDEA forecast-time equation, the option reports a continuation
index of 1+{it:g} under input orientation or 1/(1+{it:g}) under output
orientation. Consequently, both orientations imply the additional forecast
interval ln(1+{it:g})/ln(mean ROC). The selected convex weights determine the
reported effective frontier date and are returned in
{cmd:r(continuation_lambda)}.

{pstd}
Results produced by this fallback should therefore be reported as
{it:VRS continuation estimates}, not as BCC super-efficiency estimates.
Users should inspect {cmd:cont_used}, report the normalization and ROC
assumptions, and examine sensitivity to orientation, variable selection, and
the effective-date rule. The option is intended for experimentation and
methodological evaluation rather than as a generally accepted replacement
for BCC or CCR.

{marker remarks}{...}
{title:Remarks on comparison with R TFDEA}

{pstd}
The default combination {cmd:rocfrontier(final) edate(min)} is intended to
match the archived original R TFDEA workflow.  The former Stata calculation
can be requested explicitly with
{cmd:rocfrontier(first) edate(weighted)}.

{pstd}
Efficiency scores must be compared using the same orientation and convention.
For an output-oriented ordinary DEA model, {cmd:tfdea} reports the expansion
factor theta, normally theta >= 1.  Some R output reports its reciprocal,
normally between zero and one.  Thus the two results should be compared after
taking the reciprocal of one convention where necessary.

{pstd}
DEA linear programs can have alternate optimal intensity vectors.  Two
solvers may therefore report the same efficiency score but different
reference weights and effective dates.  The {cmd:edate(min)} and
{cmd:edate(max)} options make this choice explicit by solving a secondary
lexicographic LP.  Small numerical differences may remain because of solver
tolerances or degenerate optimal solutions.

{marker results}{...}
{title:Stored results}

{synoptset 24 tabbed}{...}
{synopt:{cmd:r(mean_roc)}}mean annual rate of change{p_end}
{synopt:{cmd:r(num_roc_dmu)}}number of DMUs contributing a valid rate{p_end}
{synopt:{cmd:r(edate)}}effective-date selection method{p_end}
{synopt:{cmd:r(rocfrontier)}}ROC evaluation-frontier method{p_end}
{synopt:{cmd:r(continuation)}}{cmd:continuation} when enabled{p_end}
{synopt:{cmd:r(rocrslt)}}DMU-level matrix: {cmd:tk theta_tk theta_tf tk_eff t_eff_roc ROC}{p_end}
{synopt:{cmd:r(theta_all)}}sequential-frontier efficiency matrix{p_end}
{synopt:{cmd:r(lambda)}}last sequential-frontier intensity matrix{p_end}
{synopt:{cmd:r(forecasting)}}forecast matrix: {cmd:tk SE_tf tf_eff tf_exp}; with
{cmd:continuation}, adds {cmd:cont_gap cont_used}{p_end}
{synopt:{cmd:r(continuation_lambda)}}forecast-row by training-row continuation
weight matrix; returned with {cmd:continuation}{p_end}

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

{pstd}Synthetic products illustrate forecasting and the continuation option;
they are not an empirical technology study.{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str4 product year price weight performance}{p_end}
{phang2}{cmd:. "A" 2000 10 20 10}{p_end}
{phang2}{cmd:. "B" 2001 9 18 11}{p_end}
{phang2}{cmd:. "C" 2002 8 16 12}{p_end}
{phang2}{cmd:. "D" 2003 6 12 13}{p_end}
{phang2}{cmd:. end}{p_end}
{phang2}{cmd:. tfdea price weight = performance, tk(year) dmu(product)}{p_end}
{phang2}{cmd:. tfdea price weight = performance, tk(year) dmu(product) tf(2002) rts(crs) ort(input)}{p_end}
{phang2}{cmd:. tfdea price = performance, tk(year) dmu(product) tf(2002) rts(vrs) ort(output) edate(min)}{p_end}
{phang2}{cmd:. tfdea price = performance, tk(year) dmu(product) tf(2002) rts(vrs) ort(output) continuation}{p_end}
{phang2}{cmd:. matrix list r(forecasting)}{p_end}
{phang2}{cmd:. matrix list r(rocrslt)}{p_end}
{phang2}{cmd:. return list}{p_end}

{pstd}Defaults are {cmd:rocfrontier(final) edate(min)}. To reproduce the former policy, specify {cmd:rocfrontier(first) edate(weighted)}. Without continuation, the held-out VRS target in this example has missing forecasts because it lies outside the training hull.{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker references}{...}
{title:References}

{phang}Inman, L., D. Anderson, and R. Harmon. 2006. Predicting U.S. jet
fighter aircraft introductions from 1944 to 1982: A dogfight between
regression and TFDEA. {it:Technological Forecasting and Social Change}
73: 1178-1187.{p_end}

{phang}Xue, M., and P. T. Harker. 2002. Note: Ranking DMUs with infeasible
super-efficiency DEA models. {it:Management Science} 48(5): 705-710.{p_end}

{title:Authors}

{pstd}Choonjoo Lee and Byung-Ki Jung.  Initial Stata implementation by
Byung-Ki Jung; release integration and maintenance by Choonjoo Lee.{p_end}
