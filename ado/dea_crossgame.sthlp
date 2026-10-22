{smcl}
{* Help revision 2026-10-22.2: duplicate author email removed; self-contained examples retained.}{...}
{* version 2.2.0  03may2026}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_crosseff" "help dea_crosseff"}{...}
{vieweralsosee "dea_nash" "help dea_nash"}{...}
{vieweralsosee "dea_cooperative" "help dea_cooperative"}{...}
{title:Title}

{phang}
{bf:dea_crossgame} {hline 2} Cross-Efficiency Game DEA

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:dea_crossgame} {ifin}{cmd:,}
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
{synopt:{opt strategy(string)}}game strategy: {cmd:aggressive}, {cmd:benevolent}, {cmd:nash}, or {cmd:all} (default){p_end}
{synopt:{opt maxiter(#)}}maximum iterations for Nash convergence; default {cmd:100}{p_end}
{synopt:{opt generate(string)}}prefix for generated variables; default {cmd:_cg}{p_end}
{synopt:{opt tolerance(#)}}convergence tolerance; default {cmd:1e-6}{p_end}
{synopt:{opt matrix}}display full cross-efficiency matrix{p_end}
{synopt:{opt detail}}display detailed output{p_end}
{synopt:{opt debug}}display debug information{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_crossgame} extends the cross-efficiency framework of
{help dea_crossgame##DG1994:Doyle and Green (1994)}
by incorporating strategic behavior among DMUs. In the standard cross-efficiency
model, each DMU chooses optimal weights primarily for self-appraisal. The game
DEA approach ({help dea_crossgame##LWCZ2008:Liang et al. 2008}) models how DMUs would choose weights if they
anticipated — and responded to — other DMUs' choices.

{pstd}
{ul:Three competitive strategies:}

{phang2}
{bf:Aggressive}: Each DMU chooses weights that maximize its own self-appraisal
score {it:and}, as a secondary goal, minimizes the cross-efficiency scores of
competitors. This is the fully competitive ("zero-sum") strategy.

{phang2}
{bf:Benevolent}: Each DMU chooses weights that maximize its own score
{it:and}, as a secondary goal, maximizes others' cross-efficiency scores.
This is the cooperative ("win-win") strategy.

{phang2}
{bf:Nash}: The iterative game cross-efficiency algorithm of
{help dea_crossgame##LWCZ2008:Liang et al. (2008)} solves best-response
multiplier programs subject to preserving each peer DMU's current game
cross-efficiency score. Convergence is monitored up to {opt maxiter()}
iterations.

{phang2}
{bf:All}: Compute all three strategies and report results together (default).

{pstd}
{ul:Comparison:}

{pstd}
The three strategies typically yield different efficiency rankings.
Aggressive scores are lower (harder to achieve), benevolent scores are higher,
and Nash scores lie in between. Comparing them reveals the sensitivity of
rankings to competitive assumptions.

{marker options}{...}
{title:Options}

{phang}
{opt strategy(string)} specifies which game strategy to use.
{cmd:aggressive}, {cmd:benevolent}, {cmd:nash}, or {cmd:all} (default).

{phang}
{opt maxiter(#)} sets the maximum number of iterations for Nash equilibrium
computation. Default is 100. If convergence is not achieved within {opt maxiter()}
iterations, the last iterate is returned.

{phang}
{opt matrix} displays the full n×n cross-efficiency matrix for each strategy.

{phang}
{opt generate(string)} sets the prefix for generated variables (default {cmd:_cg}).
The following variables are created depending on {opt strategy()}:

{p2colset 9 22 24 2}
{p2col:{it:prefix}{cmd:_agg}}aggressive cross-efficiency score{p_end}
{p2col:{it:prefix}{cmd:_ben}}benevolent cross-efficiency score{p_end}
{p2col:{it:prefix}{cmd:_nash}}Nash equilibrium cross-efficiency score{p_end}

{marker results}{...}
{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n)}}number of DMUs{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(strategy)}}strategy used{p_end}

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


{pstd}Setup{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str6 dmu capital labor output}{p_end}
{phang2}{cmd:"Firm A" 12  8  95}{p_end}
{phang2}{cmd:"Firm B" 18 12 125}{p_end}
{phang2}{cmd:"Firm C" 15 10 118}{p_end}
{phang2}{cmd:"Firm D" 22 15 132}{p_end}
{phang2}{cmd:"Firm E" 10  7  90}{p_end}
{phang2}{cmd:"Firm F" 20 13 128}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Run all three strategies (default){p_end}
{phang2}{cmd:. dea_crossgame, inputs(capital labor) outputs(output) dmu(dmu)}{p_end}

{pstd}Nash equilibrium only{p_end}
{phang2}{cmd:. dea_crossgame, inputs(capital labor) outputs(output) dmu(dmu) strategy(nash) maxiter(200)}{p_end}

{pstd}Compare aggressive vs. benevolent{p_end}
{phang2}{cmd:. dea_crossgame, inputs(capital labor) outputs(output) dmu(dmu) strategy(all)}{p_end}
{phang2}{cmd:. list dmu _cg_agg _cg_ben _cg_nash}{p_end}

{pstd}Display cross-efficiency matrix for aggressive strategy{p_end}
{phang2}{cmd:. dea_crossgame, inputs(capital labor) outputs(output) dmu(dmu) strategy(aggressive) matrix}{p_end}


{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker references}{...}
{title:References}

{marker DG1994}{...}
{phang}
Doyle, J., and R. Green. 1994. Efficiency and cross-efficiency in DEA: Derivations,
meanings and uses. {it:Journal of the Operational Research Society} 45(5): 567–578.
{browse "https://doi.org/10.1057/jors.1994.84":doi:10.1057/jors.1994.84}

{marker LWCZ2008}{...}
{phang}
Liang, L., J. Wu, W.D. Cook, and J. Zhu. 2008. The DEA game cross-efficiency model
and its Nash equilibrium. {it:Operations Research} 56(5): 1278–1288.
{browse "https://doi.org/10.1287/opre.1070.0487":doi:10.1287/opre.1070.0487}

{marker WLY2009}{...}
{phang}
Wu, J., L. Liang, and F. Yang. 2009. Achievement and benchmarking of countries
at the Summer Olympics using cross efficiency evaluation method.
{it:European Journal of Operational Research} 197(2): 722–730.
{browse "https://doi.org/10.1016/j.ejor.2008.06.030":doi:10.1016/j.ejor.2008.06.030}

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
