{smcl}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* version 1.2.0  02nov2026}{...}
{viewerjumpto "Related article" "dea_nash##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_cooperative" "help dea_cooperative"}{...}
{vieweralsosee "dea_crossgame" "help dea_crossgame"}{...}
{vieweralsosee "dea_stackelberg" "help dea_stackelberg"}{...}
{title:Title}

{phang}
{bf:dea_nash} {hline 2} DEA historical bargaining-style allocation

{marker syntax}{...}
{title:Syntax}

{pstd}
Form A: direct input-output syntax

{p 8 17 2}
{cmd:dea_nash} {it:inputvars} {cmd:=} {it:outputvars} {ifin}{cmd:,}
[{it:options}]

{pstd}
Form B: option syntax

{p 8 17 2}
{cmd:dea_nash} {ifin}{cmd:,}
{cmd:inputs(}{it:varlist}{cmd:)}
{cmd:outputs(}{it:varlist}{cmd:)}
[{it:options}]

{pstd}
Both forms support {cmd:if} and {cmd:in}. Equality tests such as
{cmd:if group==1}, quoted strings, and functions such as
{cmd:if inlist(group,1,2)} are parsed as sample restrictions.
The restrictions apply to both evaluated observations and the reference data.
Input and output varlists may use standard variable ranges and wildcards.

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt inputs(varlist)}}input variables{p_end}
{synopt:{opt outputs(varlist)}}output variables{p_end}

{syntab:Optional}
{synopt:{opt dmu(varname)}}DMU identifier variable{p_end}
{synopt:{opt rts(string)}}{cmd:crs} (default) or {cmd:vrs}{p_end}
{synopt:{opt gen:erate(string)}}prefix for generated variables; default {cmd:_nash}{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}
{synopt:{opt tol:erance(#)}}positive finite LP tolerance; default {cmd:1e-7}{p_end}
{synopt:{opt mat:rix}}return result and LP diagnostic matrices{p_end}
{synopt:{opt det:ail}}display the historical allocation rule{p_end}
{synopt:{opt debug}}display debug information{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_nash} retains a historical bargaining-style heuristic. It does not
optimize the Nash product over a specified feasible utility set, solve a
centralized DEA model, or estimate gains from pooling resources. The command
name and existing output suffixes are retained for compatibility.
{cmd:r(formulation)} is {cmd:historical_capped_proportional}.

{pstd}
The exact implemented rule is as follows. Let d_i be the input-oriented DEA
score. Define the scenario c_i = d_i/(0.8+0.2*d_i), except c_i=1 when
{cmd:d_i>=1-tolerance}. The constants 0.2 and 0.01 below are historical choices,
not parameters estimated from the data. Compute:

{p 8 8 2}
G = sum_i(c_i-d_i), w_i = max(0.01,1-d_i),
q_i = G*w_i/sum_j(w_j), b_i = min(1,d_i+q_i), g_i = b_i-d_i.

{pstd}
These inefficiency-dependent weights are unequal in general. The cap may discard
part of G, including the positive share assigned to an efficient DMU by the
0.01 weight floor. The discarded gain is not redistributed: it is reported in
{cmd:r(unallocated_gain)} = G - sum_i(g_i). This preserves the historical rule
while exposing its limitation. Floating-point summation may leave a tiny signed
roundoff residual even when no cap binds. The scenario is not an individual
upper bound on the allocation. The final gain is {cmd:bargain-disagree}, not
{cmd:coop-disagree}.

{pstd}
If any DEA score is unavailable, that row's scenario is missing and all
allocation/gain values and aggregate gain/product summaries are missing. Known
DEA scores and their scenarios remain visible. The sample is not silently
reduced for allocation. {cmd:r(allocation_complete)} is then zero.

{pstd}
{cmd:r(nash_product)} is the unweighted product of the computed gains, not an
optimized objective. If any gain is exactly zero, the product is zero,
{cmd:r(n_zero_gain)}>0, and {cmd:r(log_product)} is missing (log zero).
Otherwise {cmd:r(log_product)} stores sum(log(g_i)); if exponentiation
underflows, the product is zero and {cmd:r(product_underflow)} is one.
Near-efficient scores slightly below one can produce very small positive gains
and products. A nonzero product does not establish a Nash solution.

{pstd}
{help dea_nash##NASH1950:Nash (1950)} and
{help dea_nash##DLLCZ2011:Du et al. (2011)} are methodological background,
not claims of reproducing their bargaining optimizations.
{helpb dea_cooperative} instead uses an explicit additive characteristic
function whose Shapley values equal its DEA scores.


{pstd}
Both commands use input-oriented DEA against the entire selected sample,
including each evaluated DMU itself. Quantities must be numeric and nonnegative;
each selected row must contain a positive input and a positive output (zeros in
other entries or whole columns are allowed). Negative values or an all-zero
input/output row return {cmd:r(459)} before replacement. Missing quantities or a
missing explicit numeric/string DMU identifier exclude that observation.
An empty sample returns {cmd:r(2000)}. Results are missing outside the sample.

{pstd}
{opt tol:erance(#)} must be positive and finite. Columns are scaled by their
positive maxima. Each unchanged, equivalent LP is checked for convergence,
finite values and normalized primal residuals against
{cmd:min(tolerance,1e-8)}, matching the numerical endpoint envelope. The native
convergence tolerance uses the same cap; a stricter requested tolerance is honored.
A converged but rejected candidate may receive one stricter restart of the same
LP. Rejection means unavailable, not a proof that the model is infeasible.
The self-reference implies 0 <= efficiency <= 1. The LP includes this upper
bound. An accepted score just outside this interval by at most
{cmd:min(tolerance,1e-8)} is moved to the endpoint; the raw score and this
adjustment are retained in {cmd:r(lp_diagnostics)}. Larger excursions are
rejected. Scores already inside the interval are unchanged, including values
slightly below 1. Primal checks alone do not certify optimality.

{pstd}
{opt mat:rix} returns {cmd:r(results)} (the four generated columns in the order
listed below) and {cmd:r(lp_diagnostics)} with columns {cmd:raw_score score
accepted endpoint_adjusted}. Rows follow selected observation order, not unique
DMU labels; identifiers may repeat. The raw score is available only when the LP
passed its convergence and residual checks. Matrix results are subject to
Stata's matrix size limits. {opt debug} displays per-row LP acceptance and
endpoint adjustments.

{marker options}{...}
{title:Options}

{phang}
{opt rts(string)} specifies returns to scale: {cmd:crs} (default) or {cmd:vrs}.

{phang}
{opt det:ail} displays the allocation formula. The result table and allocated/
unallocated gain summary are displayed by default.

{phang}
{opt gen:erate(string)} sets the prefix for generated variables (default {cmd:_nash}).
The following variables are created:

{p2colset 9 26 28 2}
{p2col:{it:prefix}{cmd:_disagree}}disagreement point (standard DEA efficiency){p_end}
{p2col:{it:prefix}{cmd:_coop}}historical scenario c_i{p_end}
{p2col:{it:prefix}{cmd:_bargain}}capped proportional allocation b_i{p_end}
{p2col:{it:prefix}{cmd:_gain}}allocated gain (bargain - disagree){p_end}

{marker results}{...}
{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n)}}number of DMUs{p_end}
{synopt:{cmd:r(nash_product)}}unweighted product of allocated gains{p_end}
{synopt:{cmd:r(mean_gain)}}mean allocated gain{p_end}
{synopt:{cmd:r(total_gain)}}scenario gain G{p_end}
{synopt:{cmd:r(allocated_gain)}}sum of final gains{p_end}
{synopt:{cmd:r(unallocated_gain)}}G minus allocated gains; cap shortfall plus roundoff{p_end}
{synopt:{cmd:r(log_product)}}sum of log gains if every gain is positive{p_end}
{synopt:{cmd:r(n_zero_gain)}}number of exactly zero gains{p_end}
{synopt:{cmd:r(n_capped)}}number with proposed d_i+q_i > 1{p_end}
{synopt:{cmd:r(product_underflow)}}1 if positive-gain product underflowed, otherwise 0{p_end}
{synopt:{cmd:r(allocation_complete)}}1 for a complete allocation, 0 if any score unavailable{p_end}

{synopt:{cmd:r(n_lp)}}number of initial DEA LP calls (n); a residual-check restart is not another DMU{p_end}
{synopt:{cmd:r(n_failed)}}number of unavailable individual DEA scores{p_end}
{synopt:{cmd:r(n_boundary)}}number of accepted scores moved to a numerical endpoint{p_end}
{synopt:{cmd:r(complete)}}1 if all individual DEA scores are available, otherwise 0{p_end}
{synopt:{cmd:r(tolerance)}}requested tolerance{p_end}


{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(rts)}}returns to scale used{p_end}
{synopt:{cmd:r(orientation)}}{cmd:input}{p_end}
{synopt:{cmd:r(formulation)}}implemented convention, described above{p_end}

{p2col 5 20 24 2: Matrices (with matrix)}{p_end}
{synopt:{cmd:r(results)}}four result columns, with {cmd:matrix}{p_end}
{synopt:{cmd:r(lp_diagnostics)}}raw/used score, acceptance and endpoint adjustment, with {cmd:matrix}{p_end}


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


{pstd}Historical allocation scenario with gain-shortfall diagnostics{p_end}

{pstd}Setup{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str8 agency budget staff output}{p_end}
{phang2}{cmd:"Agency A" 120 24  90}{p_end}
{phang2}{cmd:"Agency B" 160 30 118}{p_end}
{phang2}{cmd:"Agency C" 140 27 105}{p_end}
{phang2}{cmd:"Agency D" 190 36 132}{p_end}
{phang2}{cmd:"Agency E" 110 22  86}{p_end}
{phang2}{cmd:"Agency F" 175 33 126}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Option syntax{p_end}
{phang2}{cmd:. dea_nash, inputs(budget staff) outputs(output) dmu(agency)}{p_end}

{pstd}Direct input-output syntax with VRS{p_end}
{phang2}{cmd:. dea_nash budget staff = output, dmu(agency) rts(vrs) detail replacevars}{p_end}

{pstd}Compare disagreement and cooperative efficiency{p_end}
{phang2}{cmd:. dea_nash budget staff = output, dmu(agency) replacevars}{p_end}
{phang2}{cmd:. list agency _nash_disagree _nash_coop _nash_bargain _nash_gain}{p_end}
{phang2}{cmd:. sum _nash_gain}{p_end}


{pstd}Equivalent sample restrictions in both syntax forms

{phang2}{cmd:. dea_nash budget staff = output if inrange(_n,2,5), dmu(agency) generate(qa)}{p_end}
{phang2}{cmd:. dea_nash if inrange(_n,2,5), inputs(budget staff) outputs(output) dmu(agency) generate(qb)}{p_end}
{phang2}{cmd:. assert qa_bargain==qb_bargain}{p_end}

{pstd}Inspect gain allocation and LP diagnostics{p_end}
{phang2}{cmd:. dea_nash budget staff=output, dmu(agency) gen(audit) matrix}{p_end}
{phang2}{cmd:. assert audit_gain==audit_bargain-audit_disagree}{p_end}
{phang2}{cmd:. display r(total_gain) " " r(allocated_gain) " " r(unallocated_gain)}{p_end}
{phang2}{cmd:. matrix list r(lp_diagnostics)}{p_end}
{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker references}{...}
{title:References}

{marker DLLCZ2011}{...}
{phang}
Du, J., L. Liang, Y. Chen, W. D. Cook, and J. Zhu. 2011. A bargaining game model
for measuring performance of two-stage network structures.
{it:European Journal of Operational Research} 210(2): 390-397.
{browse "https://doi.org/10.1016/j.ejor.2010.08.025":doi:10.1016/j.ejor.2010.08.025}

{marker NASH1950}{...}
{phang}
Nash, J. F. 1950. The bargaining problem.
{it:Econometrica} 18(2): 155-162.
{browse "https://doi.org/10.2307/1907266":doi:10.2307/1907266}

{marker ZHU2014}{...}
{phang}
Zhu, J. 2014.
{it:Quantitative Models for Performance Evaluation and Benchmarking:
Data Envelopment Analysis with Spreadsheets}. 3rd ed. Cham: Springer.
{browse "https://doi.org/10.1007/978-3-319-06647-9":doi:10.1007/978-3-319-06647-9}

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
