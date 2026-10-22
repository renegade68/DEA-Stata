{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.0.0  20jul2026}{...}
{title:dea_reference_sensitivity — Reference classification stability over tolerances}

{title:Syntax}
{p 8 12 2}
{cmd:dea_reference_sensitivity} {ifin}, {cmd:inputs(}{it:varlist}{cmd:)}
{cmd:outputs(}{it:varlist}{cmd:)}
[{cmd:optgrid(}{it:numlist}{cmd:)} {cmd:peergrid(}{it:numlist}{cmd:)}
{cmd:saving(}{it:filename}{cmd:)} {cmd:replace} and core model options]

{title:Description}
{pstd}
The command reruns {help dea_reference} over the Cartesian product of optimal-
set and positive-peer tolerances. For each evaluated-DMU--peer pair it reports
the fraction of scenarios in which the peer is possible, necessary, or
optional, and flags any change in the 0/1/2 classification. Defaults are
{cmd:optgrid(1e-6 1e-5 1e-4)} and {cmd:peergrid(1e-8 1e-7 1e-6)}.

{pstd}
This is a numerical-identification sensitivity analysis, not statistical
sampling uncertainty. A robust necessary peer is necessary in every scenario;
a robust impossible peer is impossible in every scenario. The grid must be
reported with substantive results.

{title:Stored results}
{p2colset 9 28 30 2}{...}
{p2col:{cmd:r(scenarios)}}number of tolerance combinations{p_end}
{p2col:{cmd:r(unstable_pairs)}}pairs whose classification changes{p_end}
{p2col:{cmd:r(robust_necessary)}}pairs necessary in every scenario{p_end}
{p2col:{cmd:r(robust_possible)}}always possible but not always necessary pairs{p_end}
{p2col:{cmd:r(robust_impossible)}}pairs impossible in every scenario{p_end}
{p2col:{cmd:r(stability)}}frequency and change matrix{p_end}
{p2colreset}{...}

{title:Example}

{pstd}Copy the command lines below (omit leading dots) from {cmd:preserve}
through {cmd:restore} into the Do-file Editor and run them as one block.
Examples use synthetic
data or Stata's bundled {cmd:auto} data and require no downloaded datasets.
The block restores the data in memory. Run dependent lines in order.{p_end}
{* DEAS_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input id x1 x2 y1 y2 region}{p_end}
{phang2}{cmd:. 1 1 2 1 1 1}{p_end}
{phang2}{cmd:. 2 2 4 1 1 1}{p_end}
{phang2}{cmd:. 3 3 6 1 1 2}{p_end}
{phang2}{cmd:. 4 4 8 1 1 2}{p_end}
{phang2}{cmd:. end}{p_end}
{phang2}{cmd:. generate priority = id}{p_end}
{phang2}{cmd:. generate distance = id}{p_end}
{phang2}{cmd:. generate cost = id}{p_end}
{phang2}{cmd:. generate byte review = id == 2}{p_end}
{phang2}{cmd:. matrix D = (1,2,3,4)}{p_end}
{phang2}{cmd:. tempfile peer_bounds peer_audit}{p_end}
{phang2}{cmd:. tempfile peer_stability}{p_end}

{phang2}{cmd:. dea_reference_sensitivity, inputs(x1 x2) outputs(y1 y2) dmu(id) optgrid(1e-7 1e-6 1e-5) peergrid(1e-9 1e-8 1e-7) saving("`peer_stability'") replace}{p_end}


{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{title:Also see}
{psee}{help dea_reference}, {help dea_reference_batch}
