{smcl}
{* Help revision 2026-11-01: if/in defines both evaluation and reference samples.}{...}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* version 1.3.0  09nov2026}{...}
{viewerjumpto "Related article" "dea_cooperative##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_nash" "help dea_nash"}{...}
{vieweralsosee "dea_crossgame" "help dea_crossgame"}{...}
{title:Title}

{phang}
{bf:dea_cooperative} {hline 2} Additive DEA game (exact Shapley allocation)

{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 17 2}
{cmd:dea_cooperative}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}
[{cmd:,} {it:options}]

{pstd}
Form B (extended syntax)

{p 8 17 2}
{cmd:dea_cooperative} {ifin}{cmd:,}
{cmd:inputs(}{it:varlist}{cmd:)}
{cmd:outputs(}{it:varlist}{cmd:)}
[{it:options}]

{pstd}
{cmd:if} and {cmd:in} restrict both the DMUs evaluated and the reference
set. Observations with missing input or output values are excluded. Unrelated
variables do not affect sample selection. Generated results are missing outside
the selected sample. An empty estimation sample returns {cmd:r(2000)}.
When {cmd:replacevars} is used, it replaces entire result columns, including
setting observations outside the selected sample to missing.{p_end}

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt inputs(varlist)}}input variables{p_end}
{synopt:{opt outputs(varlist)}}output variables{p_end}

{syntab:Optional}
{synopt:{opt dmu(varname)}}DMU identifier variable{p_end}
{synopt:{opt rts(string)}}{cmd:crs} (default) or {cmd:vrs}{p_end}
{synopt:{opt gen:erate(string)}}prefix for generated variables; default {cmd:_coop}{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}
{synopt:{opt tol:erance(#)}}positive finite LP tolerance; default {cmd:1e-7}{p_end}
{synopt:{opt mat:rix}}return result and LP diagnostic matrices{p_end}
{synopt:{opt det:ail}}list additive coalition values for n<=15{p_end}
{synopt:{opt debug}}display debug information{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_cooperative} computes the exact Shapley allocation for the additive,
fixed-frontier game {cmd:v(S)=sum(i in S) e_i}, with {cmd:v(empty)=0}.
Here e_i is input-oriented DEA efficiency against the full selected frontier,
not against coalition S alone. Therefore every marginal contribution
{cmd:v(S+i)-v(S)} equals e_i, and the Shapley value is exactly e_i.
The generated contribution {cmd:Shapley-e_i} is zero wherever e_i is available.
This is an additive efficiency summary: it does not estimate pooling gains,
synergies, an endogenous coalition frontier, or variable importance.

{pstd}
The Shapley formula ({help dea_cooperative##SHAPLEY1953:Shapley 1953}) reduces
analytically to this identity. Only n initial DEA LPs are required, instead of
re-solving the same member LPs for all coalitions. There is no n=20 score limit.
{cmd:detail} enumerates coalition masks only for n<=15 (bit i corresponds to
selected row i). For n>15 it prints a notice and skips the listing; individual scores
are still computed. {cmd:r(detail_shown)} records whether it was displayed.
This listing includes the empty coalition and reuses the computed scores.

{pstd}
If any LP score is unavailable, its individual/Shapley/contribution values are
missing. Known members retain their values. The grand-coalition total and all
ranks are missing because a complete sample total/ranking cannot be established.
A listed coalition is missing if it contains an unavailable member or its
exact sum cannot pass an arithmetic/export check; the
empty coalition is always zero. With all scores certified, ranks compare exact rational efficiencies; distinct
efficiencies can have identical displayed values. Otherwise complete scores use
numerical descending competition ranks (ties share a rank, with gaps).

{pstd}
{cmd:r(formulation)} is {cmd:additive_fixed_frontier}. The cited DEA game papers
({help dea_cooperative##NT2006:Nakabayashi and Tone 2006},
{help dea_cooperative##LOZ2012:Lozano 2012}, and
{help dea_cooperative##LL2010:Li and Liang 2010}) provide background; their
coalition definitions and empirical results are not reproduced by this command.
No general core-membership optimization is performed.


{pstd}
Both commands use input-oriented DEA against the entire selected sample,
including each evaluated DMU itself. Quantities must be numeric and nonnegative;
each selected row must contain a positive input and a positive output (zeros in
other entries or whole columns are allowed). Negative values or an all-zero
input/output row return {cmd:r(459)} before replacement. Missing quantities or a
missing explicit numeric/string DMU identifier exclude that observation.
An empty sample returns {cmd:r(2000)}. Results are missing outside the sample.


{marker quality}{...}
{title:Numerical quality and bounded exact verification}

{pstd}
The native LP path first scales columns and checks convergence and primal
residuals, using min(tolerance,1e-8) and at most one stricter restart. Independently,
a bounded rational path verifies the original stored input/output data using
primal and dual certificates. Its scope is n<=32 and
inputs+outputs+2*(rts==vrs)<=16. Scope, work, and representation limits can prevent
certification; a native failure is not a proof of model infeasibility.
No rounding of the original data or change of the economic model is authorized
by these checks. Each stored binary64 input and the requested tolerance are
interpreted exactly; this does not recover an unrecorded decimal measurement.

{pstd}
An available certified score replaces the native proposal. Otherwise an accepted
native score remains available but is explicitly unverified. Missing native
scores remain missing if certification cannot recover them. {cmd:r(complete)}
means that every individual score is available, not that all outputs are certified.
{cmd:r(n_certified)}, {cmd:r(n_unverified)}, and {cmd:r(n_unavailable)} count these
three score states. {cmd:r(n_native_failed)} and {cmd:r(n_native_boundary)} preserve
the native diagnostics; {cmd:r(n_failed)} and {cmd:r(n_boundary)} describe the
scores actually used. {cmd:r(n_lp)} counts initial native LP calls, not all exact
verification work or numerical restart attempts.

{pstd}
With {cmd:matrix}, {cmd:r(result_quality)} and {cmd:r(result_reason)} have the same
four columns as {cmd:r(results)}. Quality codes are 0=unavailable,
1=numerical/unverified, 2=exact expression with checked binary64 export,
3=withheld after an arithmetic or export check. Quality 2 does not mean that a
non-dyadic rational is stored without rounding: its relative export error is
checked against {cmd:tolerance}. Logarithms remain numerical (quality 1).
Derived reasons are 0=passed, 1001=unverified or unavailable source,
1004=log of an exactly zero product, 1005=numerical logarithm,
1101=positive exact value rounded to zero, 1102=invalid/nonfinite export,
1103=relative export error above tolerance, and 1201--1204=bounded-arithmetic
failure. Score-certification reasons are separate internal diagnostic codes.
A withheld output is missing and triggers a warning; other valid outputs remain.
Exactly zero and an unrepresentable positive value are distinguished.

{pstd}
{cmd:r(aggregate_diagnostics)} contains columns
{cmd:quality reason value relative_error export_ops max_limbs}.
The work columns describe export checking, or the failing arithmetic stage;
they are not the command's total work. Each aggregate also has scalar returns
{it:name}{cmd:_quality} and {it:name}{cmd:_reason}, available without {cmd:matrix}.
{cmd:r(n_withheld_results)} counts withheld result cells and
{cmd:r(n_withheld_aggregates)} counts withheld aggregate rows.
Arithmetic is bounded to 10 million charged operations and 512 base-2^26 limbs
per work context. Nash uses separate bounded contexts for each scenario/term row,
each allocation row, the total, the weight sum, their common ratio, cap-loss sum,
allocated total/mean, and product. Summation uses a balanced tree within its
context. Each export/log check also has its own context. No failed context is
retried with a larger budget. With n<=32, the derived Nash path uses at most
6*n+11 such contexts, including exports/logs; this bound excludes the independent
DEA-score verification work. A cooperative listing shares one context across
all coalition sums; the grand total has a separate context. Limits can still
withhold results inside the dimension limits.

{pstd}
The exact Nash calculation forms the common ratio once, and uses
allocated total = scenario gain - sum of exact cap losses. This is algebraically
the same capped proportional rule. Exact division checks the integer remainder
and reconstructs dividend = quotient * divisor; at the integer width boundary
it retains bounded aligned binary division. If a row fails, previously checked
rows remain available, while dependent full-sample aggregates and zero/cap counts
are withheld. An available scenario gain is retained. Failure of cap-loss or
mean aggregation does not erase valid rows or an independently computed product.
There is no switch to unverified rounded results after an exact-derived failure.

{pstd}
{cmd:r(lp_diagnostics)} retains columns {cmd:raw_score score accepted
endpoint_adjusted} for the score used. In a certified row both score columns
contain the certified export and endpoint_adjusted is zero; the original native
proposal is instead in {cmd:r(native_lp_diagnostics)}, with the same four columns.
The word raw in the compatibility matrix must not be read as native in a
certified row. {cmd:r(score_diagnostics)} has 36 columns: selected row index;
native failure site (0=accepted, 1=native optimizer, 2=converged but failed
numerical checks, 3=endpoint rejection);
native converged, error, value; native accepted score; final quality, reason,
score; and 27 internal certificate diagnostics. Rows follow selected observation
order, not distinct DMU labels. {cmd:r(lambda)} and {cmd:r(weights)} are returned
when at least one row is certified; uncertified rows are missing. Lambda columns
refer to selected reference rows. Weights contain input multipliers, output
multipliers, and (under VRS) the two nonnegative components of the free intercept.

{pstd}
{cmd:r(numerical_results)} and {cmd:r(numerical_aggregates)} retain the old
floating-point derived calculations on the final rounded scores for comparison;
they are not certificates or additional estimates. Matrix outputs remain subject
to Stata's matrix-size limits. Generated variables are binary64 values.


{pstd}
The single aggregate row is v_grand. {cmd:r(v_grand_quality)} and
{cmd:r(v_grand_reason)} describe its quality. {cmd:r(rank_quality)} and
{cmd:r(rank_reason)} describe the ranking. {cmd:detail} prints coalition quality
codes beside the values. The four scalars {cmd:r(n_unavailable_coalitions)},
{cmd:r(n_numerical_coalitions)}, {cmd:r(n_certified_coalitions)}, and
{cmd:r(n_withheld_coalitions)} count displayed masks only (zero if not listed).
A coalition can be certified when all its members are certified even if another
member outside it is unavailable. The empty coalition is exactly zero.
The grand total and ranks may be withheld independently of individual scores.
{cmd:r(allocation_quality)} is zero for this command; inspect the result and
aggregate quality codes instead.

{marker options}{...}
{title:Options}

{phang}
{opt rts(string)} specifies returns to scale: {cmd:crs} (default) or {cmd:vrs}.

{phang}
{opt gen:erate(string)} sets the prefix for generated variables (default {cmd:_coop}).
The following variables are created:

{p2colset 9 26 28 2}
{p2col:{it:prefix}{cmd:_indiv}}full-selected-frontier DEA efficiency e_i{p_end}
{p2col:{it:prefix}{cmd:_shapley}}Shapley value (φ_i){p_end}
{p2col:{it:prefix}{cmd:_contrib}}Shapley allocation minus individual efficiency{p_end}
{p2col:{it:prefix}{cmd:_rank}}competition rank; missing if any DEA score is unavailable{p_end}

{marker results}{...}
{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n)}}number of DMUs{p_end}
{synopt:{cmd:r(v_grand)}}sum of all scores; missing if any is unavailable{p_end}
{synopt:{cmd:r(detail_shown)}}1 if the requested coalition listing was displayed{p_end}
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
{synopt:{cmd:r(lp_diagnostics)}}used-score diagnostic matrix; see numerical quality{p_end}


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

{pstd}Additive Shapley allocation for six firms{p_end}
{phang2}{cmd:. dea_cooperative, inputs(capital labor) outputs(output) dmu(firm)}{p_end}

{pstd}VRS cooperative DEA with detailed coalition output{p_end}
{phang2}{cmd:. dea_cooperative capital labor = output, dmu(firm) rts(vrs) detail replacevars}{p_end}

{pstd}Rank DMUs by coalition value{p_end}
{phang2}{cmd:. dea_cooperative capital labor = output, dmu(firm) replacevars}{p_end}
{phang2}{cmd:. gsort _coop_rank}{p_end}
{phang2}{cmd:. list firm _coop_indiv _coop_shapley _coop_contrib _coop_rank}{p_end}

{pstd}
{ul:Note:} For datasets with more than 15 DMUs, computation time may be substantial
(2^n coalitions). Consider running on a subsample or using {helpb dea_nash}
as a computationally lighter alternative.


{pstd}Restrict both evaluation and reference DMUs{p_end}
{phang2}{cmd:. dea_cooperative capital labor = output if firm!="Firm A" & firm!="Firm F", dmu(firm) gen(sub)}{p_end}
{phang2}{cmd:. assert r(n)==4}{p_end}
{phang2}{cmd:. assert missing(sub_indiv) if firm=="Firm A" | firm=="Firm F"}{p_end}

{pstd}Inspect full-precision diagnostics{p_end}
{phang2}{cmd:. dea_cooperative capital labor=output, dmu(firm) gen(audit) matrix}{p_end}
{phang2}{cmd:. assert audit_shapley==audit_indiv & audit_contrib==0}{p_end}
{phang2}{cmd:. matrix list r(lp_diagnostics)}{p_end}
{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker references}{...}
{title:References}

{marker SHAPLEY1953}{...}
{phang}
Shapley, L. S. 1953. A value for n-person games. In H. W. Kuhn and
A. W. Tucker, eds., {it:Contributions to the Theory of Games, Volume II},
307--317. Princeton, NJ: Princeton University Press.
{browse "https://doi.org/10.1515/9781400881970-018":doi:10.1515/9781400881970-018}

{marker NT2006}{...}
{phang}
Nakabayashi, K., and K. Tone. 2006. Egoist's dilemma: A DEA game.
{it:Omega} 34(2): 135--148.
{browse "https://doi.org/10.1016/j.omega.2004.08.003":doi:10.1016/j.omega.2004.08.003}

{marker LOZ2012}{...}
{phang}
Lozano, S. 2012. Information sharing in DEA: A cooperative game theory approach.
{it:European Journal of Operational Research} 222(3): 558--565.
{browse "https://doi.org/10.1016/j.ejor.2012.05.014":doi:10.1016/j.ejor.2012.05.014}

{marker LL2010}{...}
{phang}
Li, Y., and L. Liang. 2010. A Shapley value index on the importance of variables
in DEA models. {it:Expert Systems with Applications} 37(9): 6287--6292.
{browse "https://doi.org/10.1016/j.eswa.2010.02.093":doi:10.1016/j.eswa.2010.02.093}

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
