{smcl}
{* Help revision 2026-10-22.2: duplicate author email removed; self-contained examples retained.}{...}
{* version 1.1.1  02may2026}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_nash" "help dea_nash"}{...}
{vieweralsosee "dea_crossgame" "help dea_crossgame"}{...}
{title:Title}

{phang}
{bf:dea_cooperative} {hline 2} Cooperative Game DEA (Shapley Value)

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

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt inputs(varlist)}}input variables{p_end}
{synopt:{opt outputs(varlist)}}output variables{p_end}

{syntab:Optional}
{synopt:{opt dmu(varname)}}DMU identifier variable{p_end}
{synopt:{opt rts(string)}}{cmd:crs} (default) or {cmd:vrs}{p_end}
{synopt:{opt generate(string)}}prefix for generated variables; default {cmd:_coop}{p_end}
{synopt:{opt tolerance(#)}}LP solver tolerance; default {cmd:1e-7}{p_end}
{synopt:{opt detail}}display coalition-level details{p_end}
{synopt:{opt debug}}display debug information{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_cooperative} applies cooperative game theory to DEA by ranking each
DMU's marginal contribution to coalitions using a Shapley-style allocation
({help dea_cooperative##SHAPLEY1953:Shapley 1953}).  The implemented
characteristic function is a pragmatic DEA coalition-value convention; related
DEA game formulations may define coalition values differently
({help dea_cooperative##NT2006:Nakabayashi and Tone 2006}).

{pstd}
{ul:Cooperative game concepts:}

{phang2}
{bf:Characteristic function v(S)}: the maximum aggregate efficiency achievable
by a coalition S of DMUs when they pool their resources and share a common frontier.

{phang2}
{bf:Shapley value φ_i}: the fair share of coalition efficiency attributable to
DMU {it:i}, based on its average marginal contribution across all possible
coalition orderings:

{p 8 17 2}
φ_i = Σ_{S⊆N\{i}} [|S|!(n−|S|−1)!/n!] × [v(S∪{i}) − v(S)]

{phang2}
{bf:Core}: The set of stable allocations where no sub-coalition would prefer
to defect (i.e., achieve more on their own than their allocated share).

{pstd}
{ul:Interpretation:}

{pstd}
The Shapley value tells each DMU its expected marginal contribution to any
coalition it joins, averaged over all coalition orderings. A high Shapley value
implies a DMU is a valuable coalition partner; a low value suggests it
contributes little relative to what other DMUs can achieve alone.

{pstd}
Related DEA cooperative game applications include information sharing among
organizations ({help dea_cooperative##LOZ2012:Lozano 2012}) and Shapley-value
measures of variable importance in DEA models
({help dea_cooperative##LL2010:Li and Liang 2010}).

{pstd}
{ul:Computational note:}

{pstd}
The exact Shapley value requires evaluating all 2^n coalitions.
{cmd:dea_cooperative} warns when n > 15 and stops when n > 20 to prevent
computational overflow.

{marker options}{...}
{title:Options}

{phang}
{opt rts(string)} specifies returns to scale: {cmd:crs} (default) or {cmd:vrs}.

{phang}
{opt generate(string)} sets the prefix for generated variables (default {cmd:_coop}).
The following variables are created:

{p2colset 9 26 28 2}
{p2col:{it:prefix}{cmd:_indiv}}individual DEA efficiency (v({i})){p_end}
{p2col:{it:prefix}{cmd:_shapley}}Shapley value (φ_i){p_end}
{p2col:{it:prefix}{cmd:_contrib}}Shapley allocation minus individual efficiency{p_end}
{p2col:{it:prefix}{cmd:_rank}}rank by Shapley value (1 = highest){p_end}

{marker results}{...}
{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n)}}number of DMUs{p_end}
{synopt:{cmd:r(v_grand)}}characteristic function value of grand coalition{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(rts)}}returns to scale used{p_end}

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

{pstd}Shapley value allocation for six firms{p_end}
{phang2}{cmd:. dea_cooperative, inputs(capital labor) outputs(output) dmu(firm)}{p_end}

{pstd}VRS cooperative DEA with detailed coalition output{p_end}
{phang2}{cmd:. dea_cooperative capital labor = output, dmu(firm) rts(vrs) detail}{p_end}

{pstd}Rank DMUs by coalition value{p_end}
{phang2}{cmd:. dea_cooperative capital labor = output, dmu(firm)}{p_end}
{phang2}{cmd:. gsort _coop_rank}{p_end}
{phang2}{cmd:. list firm _coop_indiv _coop_shapley _coop_contrib _coop_rank}{p_end}

{pstd}
{ul:Note:} For datasets with more than 15 DMUs, computation time may be substantial
(2^n coalitions). Consider running on a subsample or using {helpb dea_nash}
as a computationally lighter alternative.


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
