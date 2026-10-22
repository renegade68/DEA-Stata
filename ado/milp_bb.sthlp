{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.2.0 20jul2026}{...}
{title:Title}

{phang}{bf:milp_bb} {hline 2} solve a small mixed-integer linear program by branch-and-bound

{title:Syntax}

{p 8 16 2}{cmd:milp_bb}, {cmd:c(}{it:C}{cmd:)}
{cmd:integer(}{it:numlist}{cmd:)}
[{cmd:sense(max|min)} {cmd:aeq(}{it:Aeq}{cmd:)} {cmd:beq(}{it:beq}{cmd:)}
{cmd:aub(}{it:Aub}{cmd:)} {cmd:bub(}{it:bub}{cmd:)}
{cmd:lb(}{it:lb}{cmd:)} {cmd:ub(}{it:ub}{cmd:)}
{cmd:maxnodes(}{it:#}{cmd:)} {cmd:inttol(}{it:#}{cmd:)}
{cmd:lptol(}{it:#}{cmd:)} {cmd:lpengine(ipm|simplex|auto)}
{cmd:lpscaling(none|row)} {cmd:trace}]

{title:Description}

{pstd}{cmd:milp_bb} optimizes c*x subject to Aeq*x=beq, Aub*x<=bub,
and lb<=x<=ub.  {cmd:integer()} contains one-based column indices.  Variables
not listed there remain continuous.  Defaults are maximization, lower bounds
zero, no finite upper bounds, 10,000 nodes, integrality tolerance 1e-7, and
LP tolerance 1e-10, {cmd:lpengine(ipm)}, and {cmd:lpscaling(none)}. Stata 16
or newer is required.

{pstd}The implementation uses deterministic depth-first branch-and-bound. Node
relaxations pass through the shared internal LP contract. {cmd:lpengine(ipm)}
uses Mata's {cmd:LinearProgram()} class and preserves the version 1.0 behavior;
{cmd:lpengine(simplex)} uses the package's revised-simplex core; and
{cmd:lpengine(auto)} tries simplex first and falls back to IPM only after a
numerical failure. It is intended for teaching, reproducibility, and small or
medium experiments; it is not a replacement for a production MILP solver.

{pstd}{cmd:lpscaling(row)} divides each equality and inequality row by the
largest absolute coefficient or right-hand-side magnitude in that row and
scales the objective by its largest absolute coefficient. The reported
solution, objective, bounds, and feasibility diagnostics remain in the
original units. Scaling is opt-in so that existing numerical behavior remains
reproducible.

{title:Stored results}

{synoptset 22 tabbed}{...}
{synopt:{cmd:r(x)}}incumbent solution row vector{p_end}
{synopt:{cmd:r(objective)}}incumbent objective{p_end}
{synopt:{cmd:r(status)}}0 optimal; 1 infeasible; 2 unbounded relaxation;
3 node limit with incumbent; 4 node limit without incumbent; 5 error{p_end}
{synopt:{cmd:r(status_text)}}status description{p_end}
{synopt:{cmd:r(nodes)}}processed nodes{p_end}
{synopt:{cmd:r(lp_solves)}}LP relaxations solved{p_end}
{synopt:{cmd:r(simplex_solves)}}simplex adapter calls{p_end}
{synopt:{cmd:r(ipm_solves)}}IPM adapter calls{p_end}
{synopt:{cmd:r(lp_fallbacks)}}simplex-to-IPM fallbacks{p_end}
{synopt:{cmd:r(objective_scale)}}objective divisor used by row scaling{p_end}
{synopt:{cmd:r(row_scale_min)}}smallest constraint-row multiplier{p_end}
{synopt:{cmd:r(row_scale_max)}}largest constraint-row multiplier{p_end}
{synopt:{cmd:r(rootbound)}}root LP objective{p_end}
{synopt:{cmd:r(hasinc)}}whether an incumbent exists{p_end}
{synopt:{cmd:r(integrality_gap)}}zero after an exhaustive optimal search,
missing otherwise{p_end}
{synopt:{cmd:r(cmd)}}{cmd:milp_bb}{p_end}
{synopt:{cmd:r(sense)}}requested objective direction{p_end}
{synopt:{cmd:r(engine)}}relaxation-engine family{p_end}
{synopt:{cmd:r(lpengine)}}requested LP adapter policy{p_end}
{synopt:{cmd:r(lpscaling)}}requested LP scaling policy{p_end}
{synopt:{cmd:r(integer)}}declared integer-variable columns{p_end}
{synopt:{cmd:r(maxnodes)}}configured node limit{p_end}
{synopt:{cmd:r(inttol)}}integrality tolerance{p_end}
{synopt:{cmd:r(lptol)}}LP convergence tolerance{p_end}

{title:Example}

{pstd}Copy the command lines below (omit leading dots) from {cmd:preserve}
through {cmd:restore} into the Do-file Editor and run them as one block.
Examples use synthetic
data or Stata's bundled {cmd:auto} data and require no downloaded datasets.
The block restores the data in memory. Run dependent lines in order.{p_end}
{* DEAS_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}


{phang2}{cmd:. matrix C=(3,2)}{p_end}
{phang2}{cmd:. matrix A=(2,1 \ 1,2)}{p_end}
{phang2}{cmd:. matrix b=(4 \ 4)}{p_end}
{phang2}{cmd:. milp_bb, c(C) aub(A) bub(b) integer(1 2) sense(max)}{p_end}
{phang2}{cmd:. milp_bb, c(C) aub(A) bub(b) integer(1 2) sense(max) lpengine(auto)}{p_end}
{phang2}{cmd:. matrix list r(x)}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

