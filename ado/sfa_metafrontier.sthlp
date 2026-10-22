{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.0.1 20jul2026}{...}
{title:Title}

{p2colset 5 25 27 2}{...}
{p2col:{cmd:sfa_metafrontier} {hline 2}}Two-step stochastic metafrontier comparison{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:sfa_metafrontier} {it:ln_output ln_inputs} {ifin}{cmd:,}
{cmd:group(}{it:varname}{cmd:)} [{cmd:method(huang|bro)} {it:options}]

{title:Description}

{pstd}The variables are the logged output and logged production-function regressors.
The command first estimates one stochastic frontier per group. {cmd:method(bro)}
fits the second-stage deterministic enveloping LP of Battese, Rao, and O'Donnell
(2004). {cmd:method(huang)}, the default, estimates the second stage as a
stochastic frontier following Huang, Huang, and Liu (2014).{p_end}

{pstd}The comparison isolates Huang et al.'s central criticism: the mixed LP
approach treats estimated group frontiers as error-free. The stochastic second
stage instead separates group-frontier estimation error from a nonnegative
technology-gap component.{p_end}

{pstd}This is a companion sensitivity command, not a replacement for the
nonparametric {cmd:dea_metafrontier} model. DEA and SFA TGRs arise from
different technologies and noise assumptions and should not be presented as
the same estimand. Huang et al.'s criticism concerns the BRO SFA-plus-
programming estimator, not the pooled DEA calculation.{p_end}

{title:Options}

{synoptset 28 tabbed}{...}
{synopt:{cmd:method(huang|bro)}}second-stage estimator; default {cmd:huang}{p_end}
{synopt:{cmd:distribution(hnormal|tnormal)}}inefficiency distribution; default {cmd:hnormal}{p_end}
{synopt:{cmd:cm(varlist)}}conditional-mean model for truncated-normal Huang stage{p_end}
{synopt:{cmd:vce(robust|oim|opg)}}second-stage VCE; default {cmd:robust}{p_end}
{synopt:{cmd:generate(prefix)}}generated-variable prefix; default {cmd:_smf}{p_end}
{synopt:{cmd:dmu(varname)}}identifier used with {cmd:detail}{p_end}
{synopt:{cmd:detail}}list observation-level results{p_end}
{synopt:{cmd:replace}}replace existing generated variables{p_end}

{pstd}{cmd:if} and {cmd:in} define one common estimation sample for every
group and the second stage. Numeric and string variables are accepted by
{cmd:group()}. With {cmd:replace}, existing result variables are changed only
after all group and meta-frontier estimations finish successfully.{p_end}

{title:Important Stata limitation}

{pstd}Stata's official {cmd:frontier} command rejects {cmd:vce(robust)} with
{cmd:distribution(tnormal)}. Huang et al. require a sandwich covariance matrix
because the second-stage dependent variable is estimated. Therefore the default
is the half-normal special case with robust QML covariance. A truncated-normal
run is allowed only with {cmd:vce(oim)} or {cmd:vce(opg)} and must not be
described as Huang-correct robust inference.{p_end}

{title:Generated variables}

{pstd}{it:prefix}{cmd:_groupfrontier}, {it:prefix}{cmd:_groupte},
{it:prefix}{cmd:_metafrontier}, {it:prefix}{cmd:_tgr}, and
{it:prefix}{cmd:_mte}. The identity is MTE = group TE x TGR.{p_end}

{title:Examples}

{pstd}Copy the command lines below (omit leading dots) from {cmd:preserve}
through {cmd:restore} into the Do-file Editor and run them as one block.
Examples use synthetic
data or Stata's bundled {cmd:auto} data and require no downloaded datasets.
The block restores the data in memory. Run dependent lines in order.{p_end}
{* DEAS_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set seed 210714}{p_end}
{phang2}{cmd:. set obs 120}{p_end}
{phang2}{cmd:. generate byte region = 1 + (_n > 60)}{p_end}
{phang2}{cmd:. generate double lnx = rnormal()}{p_end}
{phang2}{cmd:. generate double lny = 1 + .6*lnx - .25*(region==2) + .04*rnormal() - abs(.12*rnormal())}{p_end}
{phang2}{cmd:. sfa_metafrontier lny lnx, group(region) method(huang) generate(h)}{p_end}
{phang2}{cmd:. assert abs(h_mte-h_groupte*h_tgr) < 1e-10}{p_end}
{phang2}{cmd:. sfa_metafrontier lny lnx, group(region) method(bro) generate(b)}{p_end}
{phang2}{cmd:. assert abs(b_mte-b_groupte*b_tgr) < 1e-10}{p_end}
{phang2}{cmd:. summarize h_mte b_mte}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{title:References}

{pstd}Battese, G. E., Rao, D. S. P., and O'Donnell, C. J. (2004).
Journal of Productivity Analysis 21:91-103.{p_end}
{pstd}Huang, C. J., Huang, T.-H., and Liu, N.-H. (2014).
Journal of Productivity Analysis 42:241-254.{p_end}
