{smcl}
{* Help revision 2026-10-22.2: duplicate author email removed; self-contained examples retained.}{...}
{* version 1.1.1  02may2026}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_cooperative" "help dea_cooperative"}{...}
{vieweralsosee "dea_crossgame" "help dea_crossgame"}{...}
{vieweralsosee "dea_stackelberg" "help dea_stackelberg"}{...}
{title:Title}

{phang}
{bf:dea_nash} {hline 2} Nash Bargaining DEA

{marker syntax}{...}
{title:Syntax}

{pstd}
Form A: direct input-output syntax

{p 8 17 2}
{cmd:dea_nash} {it:inputvars} {cmd:=} {it:outputvars}{cmd:,}
[{it:options}]

{pstd}
Form B: option syntax

{p 8 17 2}
{cmd:dea_nash} {ifin}{cmd:,}
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
{synopt:{opt generate(string)}}prefix for generated variables; default {cmd:_nash}{p_end}
{synopt:{opt tolerance(#)}}LP solver tolerance; default {cmd:1e-7}{p_end}
{synopt:{opt detail}}accepted for interface consistency{p_end}
{synopt:{opt debug}}display debug information{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_nash} implements a symmetric Nash-bargaining-style DEA allocation
inspired by bargaining DEA applications such as
{help dea_nash##DLLCZ2011:Du et al. (2011)} and the classical Nash bargaining
framework ({help dea_nash##NASH1950:Nash 1950}). The command computes DEA
disagreement points, cooperative-efficiency proxies, and a bargained allocation
of efficiency gains across DMUs.

{pstd}
{ul:Nash bargaining problem:}

{pstd}
Each DMU {it:i} has a disagreement point d_i - its efficiency under
non-cooperation (standard DEA). When DMUs cooperate (e.g., by pooling resources),
they can potentially achieve higher joint efficiency u_i. In the Nash
bargaining framework, the solution maximizes the product of net gains:

{p 8 17 2}
max Π_i (u_i − d_i)^α_i

{pstd}
where α_i is DMU {it:i}'s bargaining power. The current implementation assumes
{it:symmetric bargaining} and allocates cooperative gains using
inefficiency-based bargaining weights.

{pstd}
{ul:Key concepts:}

{p2colset 9 22 24 2}
{p2col:{bf:Disagreement point (d_i)}}the efficiency DMU {it:i} can guarantee
alone (standard DEA score) if negotiations fail.{p_end}
{p2col:{bf:Cooperative efficiency (u_i)}}the potential efficiency when DMU {it:i}
participates in the cooperative allocation.{p_end}
{p2col:{bf:Bargaining gain}}u_i - d_i: the net improvement from cooperation.{p_end}

{pstd}
{ul:Comparison with related approaches:}

{p2colset 9 26 28 2}
{p2col:{helpb dea_cooperative}}uses Shapley value - fair allocation based on
marginal contributions to coalitions; computationally intensive (2^n){p_end}
{p2col:{helpb dea_nash}}uses a Nash-bargaining-style gain allocation; scalable
to larger datasets{p_end}
{p2col:{helpb dea_crossgame}}game over weight choices (competitive/cooperative);
no explicit resource reallocation{p_end}

{marker options}{...}
{title:Options}

{phang}
{opt rts(string)} specifies returns to scale: {cmd:crs} (default) or {cmd:vrs}.

{phang}
{opt detail} is accepted for interface consistency. The current implementation
displays the bargaining table by default.

{phang}
{opt generate(string)} sets the prefix for generated variables (default {cmd:_nash}).
The following variables are created:

{p2colset 9 26 28 2}
{p2col:{it:prefix}{cmd:_disagree}}disagreement point (standard DEA efficiency){p_end}
{p2col:{it:prefix}{cmd:_coop}}cooperative-efficiency proxy{p_end}
{p2col:{it:prefix}{cmd:_bargain}}bargained allocation value{p_end}
{p2col:{it:prefix}{cmd:_gain}}bargaining gain (coop - disagree){p_end}

{marker results}{...}
{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n)}}number of DMUs{p_end}
{synopt:{cmd:r(nash_product)}}Nash product value at solution{p_end}
{synopt:{cmd:r(mean_gain)}}mean bargaining gain{p_end}

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


{pstd}Nash bargaining DEA for resource reallocation{p_end}

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
{phang2}{cmd:. dea_nash budget staff = output, dmu(agency) rts(vrs) detail}{p_end}

{pstd}Compare disagreement and cooperative efficiency{p_end}
{phang2}{cmd:. dea_nash budget staff = output, dmu(agency)}{p_end}
{phang2}{cmd:. list agency _nash_disagree _nash_coop _nash_bargain _nash_gain}{p_end}
{phang2}{cmd:. sum _nash_gain}{p_end}


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
