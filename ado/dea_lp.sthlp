{smcl}
{* *! version 1.0.0  01nov2026}{...}
{title:dea_lp — compatibility entry point for IPM DEA}

{pstd}{cmd:dea_lp} accepts the current {help dea} syntax in both forms:
{p 8 16 2}{cmd:dea_lp} {it:inputs} {cmd:=} {it:outputs} {ifin} [{cmd:,} {it:dea_options}]
{p 8 16 2}{cmd:dea_lp} {it:outputs} {ifin}{cmd:, inputs(}{it:inputs}{cmd:)} [{it:dea_options}]

{pstd}Requires Stata 16 or later. This standalone entry point works in a fresh
session. It preserves the former {cmd:dea_lp} 6.1.0 defaults:
{cmd:rts(vrs)} and {cmd:engine(ipm)}. An explicit {cmd:rts()} overrides the
technology default. Explicit {cmd:engine()} or {cmd:method()} overrides the
solver default; {cmd:engine()} takes precedence when both are supplied.
For example, an explicit {cmd:method(simplex)} now requests simplex.
Ordinary {cmd:dea} retains its own CRS/simplex defaults.

{pstd}The remaining behavior follows current {cmd:dea}, including
{cmd:stage(2)}, output protection, {cmd:replacevars}, sample selection, and
failure reporting. This is a migration entry point, not an exact copy of
every historical option or numerical result. In particular, the former
{cmd:scale()} option is not accepted by current {cmd:dea}; use its documented
score convention. Check {cmd:r(rts)}, {cmd:r(engine)}, and {cmd:r(model)}.
The results of {cmd:dea} are returned with {cmd:r(cmd)} set to {cmd:dea_lp}.

{pstd}For a general linear program with objective and constraint coefficients,
use {help lp}. For a matrix-based mixed-integer problem, use {help milp_bb}.

{title:Examples}
{pstd}Run this complete block; it restores the original data.
{* DEAS_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input double id x y}{p_end}
{phang2}{cmd:1 1 1}{p_end}
{phang2}{cmd:2 2 3}{p_end}
{phang2}{cmd:3 4 4}{p_end}
{phang2}{cmd:4 5 3}{p_end}
{phang2}{cmd:. end}{p_end}
{phang2}{cmd:. dea_lp x=y, dmu(id) gen(vrs_score)}{p_end}
{phang2}{cmd:. assert "`r(rts)'"=="vrs"}{p_end}
{phang2}{cmd:. assert "`r(engine)'"=="ipm"}{p_end}
{phang2}{cmd:. assert abs(vrs_score-cond(id==4,.4,1))<1e-6}{p_end}
{phang2}{cmd:. dea_lp y, inputs(x) dmu(id) rts(crs) engine(simplex) gen(crs_score)}{p_end}
{phang2}{cmd:. assert abs(crs_score-y/(1.5*x))<1e-6}{p_end}
{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}
