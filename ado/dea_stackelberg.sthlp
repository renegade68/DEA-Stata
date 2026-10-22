{smcl}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* version 1.1.1  02may2026}{...}
{viewerjumpto "Related article" "dea_stackelberg##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_network" "help dea_network"}{...}
{vieweralsosee "dea_crossgame" "help dea_crossgame"}{...}
{title:Title}

{phang}
{bf:dea_stackelberg} {hline 2} Stackelberg-style two-stage DEA

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:dea_stackelberg} {ifin}{cmd:,}
{cmd:inputs(}{it:varlist}{cmd:)}
{cmd:intermediates(}{it:varlist}{cmd:)}
{cmd:outputs(}{it:varlist}{cmd:)}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt inputs(varlist)}}Stage 1 input variables (X){p_end}
{synopt:{opt intermediates(varlist)}}intermediate variables linking stages (Z){p_end}
{synopt:{opt outputs(varlist)}}Stage 2 output variables (Y){p_end}

{syntab:Optional}
{synopt:{opt dmu(varname)}}DMU identifier variable{p_end}
{synopt:{opt leader(string)}}leader stage: {cmd:stage1} (default) or {cmd:stage2}{p_end}
{synopt:{opt rts(string)}}{cmd:crs} (default) or {cmd:vrs}{p_end}
{synopt:{opt generate(string)}}prefix for generated variables; default {cmd:_stack}{p_end}
{synopt:{opt tolerance(#)}}LP solver tolerance; default {cmd:1e-7}{p_end}
{synopt:{opt detail}}accepted for interface consistency{p_end}
{synopt:{opt debug}}display debug information{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_stackelberg} implements a leader-follower-style two-stage DEA
decomposition inspired by supply-chain and network DEA models such as
Liang et al. (2006), Li et al. (2012), and Tavana and Khalili-Damghani
(2014). The command reports stage-specific DEA efficiencies under a selected
leadership ordering. It is intended as a pragmatic Stackelberg-style
diagnostic, not as a full fuzzy, alpha-cut, Monte-Carlo, or bilevel
backward-induction optimization model.

{pstd}
{ul:Production structure:}

{p 8 17 2}
Stage 1: X → Z  (inputs to intermediates){break}
Stage 2: Z → Y  (intermediates to outputs)

{pstd}
{ul:Computation:}

{phang2}{bf:1.} Stage 1 efficiency is computed from X to Z.{p_end}
{phang2}{bf:2.} Stage 2 efficiency is computed from Z to Y.{p_end}
{phang2}{bf:3.} {opt leader(stage1)} or {opt leader(stage2)} determines which
stage score is reported as the leader and which as the follower.{p_end}

{pstd}
{ul:Two leader configurations:}

{p2colset 9 20 22 2}
{p2col:{cmd:stage1}}Stage 1 is reported as the leader (e.g., upstream
production is treated as the leading stage); Stage 2 is reported as the
follower.{p_end}
{p2col:{cmd:stage2}}Stage 2 is reported as the leader (e.g., downstream
demand is treated as the leading stage); Stage 1 is reported as the
follower.{p_end}

{pstd}
Results are shown together with independent stage-efficiency benchmarks.
Because the current implementation does not solve a nested leader-follower
program, large differences should be interpreted as stage-ordering diagnostics
rather than as structural estimates of a Stackelberg equilibrium.

{marker options}{...}
{title:Options}

{phang}
{opt leader(string)} specifies which stage leads. {cmd:stage1} (default) or
{cmd:stage2}.

{phang}
{opt detail} is accepted for interface consistency. The current implementation
displays the stage comparison table by default.

{phang}
{opt generate(string)} sets the prefix for generated variables (default {cmd:_stack}).
The following variables are created:

{p2colset 9 26 28 2}
{p2col:{it:prefix}{cmd:_leader}}efficiency of the stage reported as leader{p_end}
{p2col:{it:prefix}{cmd:_follower}}efficiency of the stage reported as follower{p_end}
{p2col:{it:prefix}{cmd:_overall}}overall efficiency (leader × follower){p_end}
{p2col:{it:prefix}{cmd:_indep1}}Stage 1 efficiency under independent benchmark{p_end}
{p2col:{it:prefix}{cmd:_indep2}}Stage 2 efficiency under independent benchmark{p_end}

{marker results}{...}
{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n)}}number of DMUs{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(leader)}}leader stage used{p_end}
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


{pstd}Supply chain: manufacturer (Stage 1) leads, retailer (Stage 2) follows{p_end}

{pstd}Setup{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str6 unit raw labor parts sales}{p_end}
{phang2}{cmd:"UnitA" 40 12 68 102}{p_end}
{phang2}{cmd:"UnitB" 55 15 76 116}{p_end}
{phang2}{cmd:"UnitC" 48 13 72 108}{p_end}
{phang2}{cmd:"UnitD" 62 18 84 126}{p_end}
{phang2}{cmd:"UnitE" 36 11 61  94}{p_end}
{phang2}{cmd:"UnitF" 70 20 91 138}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Stage 1 as leader{p_end}
{phang2}{cmd:. dea_stackelberg, inputs(raw labor) intermediates(parts) ///}{p_end}
{phang2}{cmd:    outputs(sales) dmu(unit) leader(stage1)}{p_end}

{pstd}Reverse: downstream stage leads{p_end}
{phang2}{cmd:. dea_stackelberg, inputs(raw labor) intermediates(parts) ///}{p_end}
{phang2}{cmd:    outputs(sales) dmu(unit) leader(stage2) rts(vrs)}{p_end}

{pstd}Compare Stackelberg vs. independent efficiency{p_end}
{phang2}{cmd:. dea_stackelberg, inputs(raw labor) intermediates(parts) ///}{p_end}
{phang2}{cmd:    outputs(sales) dmu(unit)}{p_end}
{phang2}{cmd:. list unit _stack_leader _stack_follower _stack_overall _stack_indep1 _stack_indep2}{p_end}


{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker references}{...}
{title:References}

{phang}
Liang, L., F. Yang, W.D. Cook, and J. Zhu. 2006. DEA models for supply chain
efficiency evaluation. {it:Annals of Operations Research} 145(1): 35–49.
{browse "https://doi.org/10.1007/s10479-006-0026-7":doi:10.1007/s10479-006-0026-7}

{phang}
Li, Y., Y. Chen, L. Liang, and J. Xie. 2012. DEA models for extended two-stage
network structures. {it:Omega} 40(5): 611–618.
{browse "https://doi.org/10.1016/j.omega.2011.11.007":doi:10.1016/j.omega.2011.11.007}

{phang}
Tavana, M., and K. Khalili-Damghani. 2014. A new two-stage Stackelberg fuzzy
data envelopment analysis model. {it:Measurement} 53: 277–296.
{browse "https://doi.org/10.1016/j.measurement.2014.03.030":doi:10.1016/j.measurement.2014.03.030}

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
