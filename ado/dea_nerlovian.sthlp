{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.0.0  20jul2026}{...}
{vieweralsosee "dea" "help dea"}{...}
{title:Title}

{p2colset 5 26 28 2}{...}
{p2col:{cmd:dea_nerlovian} {hline 2}}Nerlovian directional profit inefficiency{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:dea_nerlovian} {it:inputvars} {cmd:=} {it:outputvars}
{ifin}{cmd:,} {cmd:w(}{it:inputprices}{cmd:)}
{cmd:p(}{it:outputprices}{cmd:)} [{it:options}]

{p 8 16 2}{cmd:dea_nerlovian} {it:outputvars} {ifin}{cmd:,}
{cmd:inputs(}{it:varlist}{cmd:)} {cmd:w(}{it:inputprices}{cmd:)}
{cmd:p(}{it:outputprices}{cmd:)} [{it:options}]

{synoptset 31 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:w(}{it:varlist}{cmd:)}}strictly positive input prices; required{p_end}
{synopt:{cmd:p(}{it:varlist}{cmd:)}}strictly positive output prices; required{p_end}
{synopt:{cmd:gx(}{it:varlist}{cmd:)}}nonnegative input direction; default observed inputs{p_end}
{synopt:{cmd:gy(}{it:varlist}{cmd:)}}nonnegative output direction; default observed outputs{p_end}
{synopt:{cmd:rts(vrs)}}finite unconditional profit technology; the only supported RTS{p_end}
{synopt:{cmd:engine(auto|simplex|ipm)}}engine for the directional technical component{p_end}
{synopt:{cmd:generate(}{it:prefix}{cmd:)}}result prefix; default {cmd:_npi}{p_end}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier used in detailed/saved results{p_end}
{synopt:{cmd:saving(}{it:filename}{cmd:)}}save estimation sample and results{p_end}
{synopt:{cmd:replace}}replace an existing results file{p_end}
{synopt:{cmd:detail}}list observation-level decomposition{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:dea_nerlovian} computes normalized profit loss without dividing by
observed profit. Consequently zero and negative observed profits are valid.
For direction ({it:gx},{it:gy}),

{p 8 12 2}
NPI = [maximum profit - observed profit] / [{it:p gy + w gx}]
    = directional technical inefficiency + allocative inefficiency.{p_end}

{pstd}
The default ({it:gx},{it:gy})=({it:x},{it:y}) normalizes by observed revenue
plus observed cost. Version 1.0.0 requires VRS: maximizing profit over its
convex technology is finite, while an unrestricted CRS technology can have an
unbounded profit function whenever it contains a positive-profit ray.

{pstd}
Generated suffixes are {cmd:_profit}, {cmd:_maxprofit},
{cmd:_normalization}, {cmd:_npi}, {cmd:_technical}, {cmd:_allocative}, and
{cmd:_identity_gap}. This additive loss measure is distinct from
{cmd:dea, model(profit)}, which returns an observed/maximum profit ratio and
requires strictly positive observed profit.

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
{phang2}{cmd:. input str2 firm x1 x2 y1 y2 w1 w2 p1 p2}{p_end}
{phang2}{cmd:. "A" 1 2 1 1 10 2 1 1}{p_end}
{phang2}{cmd:. "B" 2 3 2 1 10 2 1 1}{p_end}
{phang2}{cmd:. "C" 3 4 2 2 10 2 1 1}{p_end}
{phang2}{cmd:. "D" 4 5 3 2 10 2 1 1}{p_end}
{phang2}{cmd:. end}{p_end}


{phang2}{cmd:. dea_nerlovian x1 x2 = y1 y2, w(w1 w2) p(p1 p2) dmu(firm) generate(npi)}{p_end}
{phang2}{cmd:. list firm npi_profit npi_npi npi_technical npi_allocative}{p_end}


{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{title:Stored results}

{p2colset 9 31 33 2}{...}
{p2col:{cmd:r(N)}}evaluated DMUs{p_end}
{p2col:{cmd:r(N_nonpositive_profit)}}DMUs with zero or negative observed profit{p_end}
{p2col:{cmd:r(mean_npi)}}mean Nerlovian profit inefficiency{p_end}
{p2col:{cmd:r(mean_technical)}}mean directional technical component{p_end}
{p2col:{cmd:r(mean_allocative)}}mean allocative component{p_end}
{p2col:{cmd:r(max_identity_gap)}}largest decomposition residual{p_end}
{p2colreset}{...}

{title:Reference}

{phang}Chambers, R. G., Y. Chung, and R. Färe. 1998. Profit, directional
distance functions, and Nerlovian efficiency. {it:Journal of Optimization
Theory and Applications} 98: 351--364.{p_end}

{title:Author}

{pstd}Choonjoo Lee, Korea National Defense University.{p_end}
