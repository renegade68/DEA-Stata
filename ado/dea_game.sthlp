{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.1.0  02may2026}{...}
{viewerjumpto "Overview" "dea_game##overview"}{...}
{viewerjumpto "Commands" "dea_game##commands"}{...}
{viewerjumpto "Comparison" "dea_game##comparison"}{...}
{viewerjumpto "Examples" "dea_game##examples"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col:{bf:Game-theoretic DEA} {hline 2}}Overview of game-theoretic DEA extensions{p_end}
{p2colreset}{...}


{marker overview}{...}
{title:Overview}

{pstd}
{help dea_game} is an overview help file for DEA commands motivated by
game-theoretic ideas. The commands listed here are intended to open and
organize discussion about strategic interaction in DEA: competition over
weights, bargaining intuition, leader-follower ordering, and cooperative
coalitions.

{pstd}
These commands should not be read as a claim that every listed procedure is a
fully developed structural game-theoretic estimator. Several commands are
pragmatic implementations that translate an idea from the literature into a
usable DEA diagnostic. They are provided to share useful starting points and to
encourage deeper, command-specific development in future work.

{pstd}
For the revision validation record, the Nakabayashi and Tone (2006)
Egoist's dilemma Shapley examples are reproduced separately as paper examples.
Those examples use a benefit-only DEA game characteristic function and should
not be interpreted as exact validation of the current diagnostic commands
listed below.

{pstd}
Game-theoretic DEA ideas are useful when:

{p2colset 9 12 14 2}{...}
{p2col:-}DMUs compete or cooperate with each other{p_end}
{p2col:-}Resources need to be fairly allocated{p_end}
{p2col:-}Hierarchical (leader-follower) relationships exist{p_end}
{p2col:-}Coalition formation affects efficiency{p_end}
{p2colreset}{...}


{marker commands}{...}
{title:Available Commands}

{synoptset 20 tabbed}{...}
{synopthdr:Command}
{synoptline}
{synopt:{helpb dea_crossgame}}Cross-efficiency game and strategy comparison{p_end}
{synopt:{helpb dea_nash}}Nash-bargaining-style allocation diagnostic{p_end}
{synopt:{helpb dea_stackelberg}}Stackelberg-style two-stage ordering diagnostic{p_end}
{synopt:{helpb dea_cooperative}}Cooperative-game coalition and Shapley-style analysis{p_end}
{synoptline}


{marker comparison}{...}
{title:Model Comparison}

{col 5}{bf:Model}{col 25}{bf:Game Type}{col 45}{bf:Key Concept}{col 65}{bf:Application}
{col 5}{hline 75}
{col 5}Cross-Efficiency{col 25}Non-cooperative{col 45}Peer evaluation{col 65}Ranking, benchmarking
{col 5}Nash-style{col 25}Cooperative{col 45}Gain allocation{col 65}Resource-sharing ideas
{col 5}Stackelberg-style{col 25}Non-cooperative{col 45}Leader-follower order{col 65}Two-stage diagnostics
{col 5}Shapley-style{col 25}Cooperative{col 45}Marginal contribution{col 65}Coalition analysis
{col 5}{hline 75}

{pstd}
{bf:Choosing the Right Model:}

{p2colset 9 30 32 2}{...}
{p2col:{bf:Cross-Efficiency}}When you want to examine how peer evaluation and
secondary weight-selection goals affect rankings. Use aggressive, benevolent,
or Nash-style strategies as sensitivity diagnostics.{p_end}

{p2col:{bf:Nash-style allocation}}When DMUs can negotiate or share resources
and you want a simple way to compare disagreement points, cooperative-efficiency
proxies, and allocated gains. Interpret the output as a Nash-bargaining-style
diagnostic rather than a full axiomatic Nash solution.{p_end}

{p2col:{bf:Stackelberg-style ordering}}When a two-stage process has an upstream
or downstream stage that may be treated as leading the interpretation. The
current command reports stage efficiencies under a selected ordering; it is not
a full bilevel backward-induction solver.{p_end}

{p2col:{bf:Shapley-style coalition analysis}}When you want to examine how each
DMU contributes to coalitions and to summarize those contributions. Computation
can be expensive because the number of coalitions grows as 2^n.{p_end}
{p2colreset}{...}


{marker examples}{...}
{title:Quick Examples}

{pstd}Copy the command lines below (omit leading dots) from {cmd:preserve}
through {cmd:restore} into the Do-file Editor and run them as one block.
Examples use synthetic
data or Stata's bundled {cmd:auto} data and require no downloaded datasets.
The block restores the data in memory. Run dependent lines in order.{p_end}
{* DEAS_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}


{pstd}Setup data:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input byte dmu_id double(x1 x2 z1 y1 y2)}{p_end}
{phang2}{cmd:  1  2  5  3.8  1.0  1.5}{p_end}
{phang2}{cmd:  2  3  2  4.5  3.0  1.0}{p_end}
{phang2}{cmd:  3  1  3  3.6  2.0  2.0}{p_end}
{phang2}{cmd:  4  4  3  5.2  2.0  1.5}{p_end}
{phang2}{cmd:  5  2  4  4.1  1.0  2.0}{p_end}
{phang2}{cmd:  6  5  4  6.0  3.5  2.2}{p_end}
{phang2}{cmd:  end}{p_end}

{pstd}{bf:1. Cross-Efficiency Game} (compare strategies):{p_end}
{phang2}{cmd:. dea_crossgame, inputs(x1 x2) outputs(y1 y2) dmu(dmu_id) strategy(all)}{p_end}

{pstd}{bf:2. Nash-style allocation} (resource-sharing intuition):{p_end}
{phang2}{cmd:. dea_nash, inputs(x1 x2) outputs(y1 y2) dmu(dmu_id)}{p_end}

{pstd}{bf:3. Shapley-style coalition analysis} (marginal contribution):{p_end}
{phang2}{cmd:. dea_cooperative, inputs(x1 x2) outputs(y1 y2) dmu(dmu_id)}{p_end}

{pstd}{bf:4. Stackelberg-style ordering} (two-stage diagnostic):{p_end}
{phang2}{cmd:. dea_stackelberg, inputs(x1) intermediates(z1) outputs(y1) dmu(dmu_id) leader(stage1)}{p_end}


{title:Theoretical Background}

{pstd}
The concepts below provide vocabulary for interpreting the commands. They are
not a guarantee that every command solves the complete mathematical game named
in the heading.

{pstd}
{bf:Non-cooperative games:} Each player maximizes its own payoff given others'
strategies. In this package, this idea appears most directly in
cross-efficiency strategy comparisons.

{pstd}
{bf:Cooperative games:} Players can form agreements and share gains. In this
package, this idea appears in coalition analysis and in simpler bargaining-style
allocation diagnostics.

{pstd}
{bf:Stackelberg games:} Sequential games with a leader and follower. The
provided Stackelberg-style command uses this vocabulary to organize two-stage
DEA scores; it does not claim to solve a full nested Stackelberg program.



{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{title:References}

{phang}
Doyle, J., and R. Green. 1994. Efficiency and cross-efficiency in DEA:
Derivations, meanings and uses. {it:Journal of the Operational Research Society}
45(5): 567-578.
{browse "https://doi.org/10.1057/jors.1994.84":doi:10.1057/jors.1994.84}

{phang}
Liang, L., J. Wu, W. D. Cook, and J. Zhu. 2008. The DEA game cross-efficiency
model and its Nash equilibrium. {it:Operations Research} 56(5): 1278-1288.
{browse "https://doi.org/10.1287/opre.1070.0487":doi:10.1287/opre.1070.0487}

{phang}
Nash, J. F. 1950. The bargaining problem. {it:Econometrica} 18(2): 155-162.
{browse "https://doi.org/10.2307/1907266":doi:10.2307/1907266}

{phang}
Nakabayashi, K., and K. Tone. 2006. Egoist's dilemma: A DEA game. 
{it:Omega} 34(2): 135-148.
{browse "https://doi.org/10.1016/j.omega.2004.08.003":doi:10.1016/j.omega.2004.08.003}

{phang}
Shapley, L. S. 1953. A value for n-person games. In H. W. Kuhn and A. W.
Tucker, eds., {it:Contributions to the Theory of Games, Volume II}, 307-317.
Princeton, NJ: Princeton University Press.
{browse "https://doi.org/10.1515/9781400881970-018":doi:10.1515/9781400881970-018}

{phang}
von Stackelberg, H. 1934. {it:Marktform und Gleichgewicht}. Vienna: Julius Springer.


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
