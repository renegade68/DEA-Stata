{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 3.0.1  21jul2026}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_crosseff" "help dea_crosseff"}{...}
{vieweralsosee "dea_malmquist" "help dea_malmquist"}{...}
{vieweralsosee "dea_window" "help dea_window"}{...}
{viewerjumpto "Syntax" "dea_bootstrap##syntax"}{...}
{viewerjumpto "Description" "dea_bootstrap##description"}{...}
{viewerjumpto "Options" "dea_bootstrap##options"}{...}
{viewerjumpto "Remarks" "dea_bootstrap##remarks"}{...}
{viewerjumpto "Examples" "dea_bootstrap##examples"}{...}
{viewerjumpto "Stored results" "dea_bootstrap##stored_results"}{...}
{viewerjumpto "References" "dea_bootstrap##references"}{...}
{viewerjumpto "Authors" "dea_bootstrap##authors"}{...}
{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{cmd:dea_bootstrap} {hline 2}}Bootstrap bias-corrected DEA efficiency scores{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea_bootstrap}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}
[{cmd:,} {it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_bootstrap}
{ifin}{cmd:,}
{cmdab:in:puts(}{it:varlist}{cmd:)}
{cmdab:out:puts(}{it:varlist}{cmd:)}
[{it:options}]


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmdab:in:puts(}{it:varlist}{cmd:)}}input variable(s){p_end}
{synopt:{cmdab:out:puts(}{it:varlist}{cmd:)}}output variable(s){p_end}

{syntab:Model}
{synopt:{cmdab:r:ts(}{it:rtstype}{cmd:)}}returns to scale; default is {cmd:rts(crs)}{p_end}
{synopt:{cmdab:or:ientation(}{it:type}{cmd:)}}orientation; default is {cmd:orientation(input)}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}

{syntab:Bootstrap}
{synopt:{cmd:reps(}{it:#}{cmd:)}}number of bootstrap replications; default is {cmd:reps(200)}{p_end}
{synopt:{cmd:level(}{it:#}{cmd:)}}confidence level; default is {cmd:level(95)}{p_end}
{synopt:{cmd:seed(}{it:#}{cmd:)}}random-number seed for reproducibility{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:newvar}{cmd:)}}name for original efficiency variable; default is {cmd:_bs_eff}{p_end}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable (string or numeric){p_end}
{synopt:{cmd:detail}}display DMU-level bootstrap details{p_end}

{syntab:Technical}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmd:debug}}display debugging information{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:rtstype}: {cmd:crs} (constant returns to scale) | {cmd:vrs} (variable returns to scale){p_end}
{p 4 6 2}
{it:type}: {cmd:input} | {cmd:output}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_bootstrap} implements the
{help dea_bootstrap##SW1998:Simar and Wilson (1998)} bootstrap procedure
for DEA efficiency scores. It provides bias-corrected efficiency estimates and
confidence intervals that account for the finite-sample bias inherent in standard
DEA estimates.

{pstd}
The command estimates original DEA efficiency scores, then uses the smooth
bootstrap of {help dea_bootstrap##SW1998:Simar and Wilson (1998)} to construct a reference distribution.
Bias-corrected scores, bootstrap standard errors, and confidence intervals are
reported for each DMU.

{pstd}
Internally, {cmd:dea_bootstrap} solves the DEA linear programs with the
package-supplied Two-Phase Revised Simplex solver in {cmd:_dea_simplex.ado}.
The command does not expose an {cmd:engine()} option; this keeps the bootstrap
replication path fixed and reproducible.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{cmdab:r:ts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
{cmd:rts(crs)} assumes constant returns to scale
({help dea_bootstrap##CCR1978:Charnes, Cooper, and Rhodes 1978}).
{cmd:rts(vrs)} assumes variable returns to scale
({help dea_bootstrap##BCC1984:Banker, Charnes, and Cooper 1984}).
Default is {cmd:rts(crs)}.

{phang}
{cmdab:or:ientation(}{it:type}{cmd:)} specifies the direction of optimization.
{cmd:orientation(input)} minimizes inputs for given outputs.
{cmd:orientation(output)} maximizes outputs for given inputs.
Default is {cmd:orientation(input)}.

{dlgtab:Bootstrap}

{phang}
{cmd:reps(}{it:#}{cmd:)} specifies the number of bootstrap replications.
A minimum of 200 is recommended; 2000 or more is preferred for publication-quality
confidence intervals. Default is {cmd:reps(200)}.

{phang}
{cmd:level(}{it:#}{cmd:)} specifies the confidence level for the bootstrap
confidence intervals. Default is {cmd:level(95)} for 95% confidence intervals.

{phang}
{cmd:seed(}{it:#}{cmd:)} sets the random-number seed for reproducible results.
Default is {cmd:seed(0)} (no seed set).

{dlgtab:Output}

{phang}
{cmdab:gen:erate(}{it:newvar}{cmd:)} specifies the name for the generated
original efficiency variable. The command also creates variables
{it:newvar}{cmd:_bc} (bias-corrected efficiency),
{it:newvar}{cmd:_se} (bootstrap standard error),
{it:newvar}{cmd:_lo} (lower confidence bound), and
{it:newvar}{cmd:_hi} (upper confidence bound).
The default is {cmd:generate(_bs_eff)}, which creates {cmd:_bs_eff},
{cmd:_bs_eff_bc}, {cmd:_bs_eff_se}, {cmd:_bs_eff_lo}, and {cmd:_bs_eff_hi}.

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies a variable identifying each DMU.
Accepts both string and numeric variables.

{dlgtab:Technical}

{phang}
The LP solver is fixed to the package-supplied Two-Phase Revised Simplex
implementation in {cmd:_dea_simplex.ado}; {cmd:engine()} is not an option for
this wrapper.


{marker remarks}{...}
{title:Remarks}

{pstd}
The {help dea_bootstrap##SW1998:Simar and Wilson (1998)} algorithm uses a smooth bootstrap based on kernel
density estimation of the original efficiency scores. The method requires
{cmd:rts(crs)} or {cmd:rts(vrs)}; other returns-to-scale assumptions are not
supported.

{pstd}
Note that the confidence intervals produced by {cmd:dea_bootstrap} are
bias-corrected confidence intervals. The original efficiency score serves as the
lower bound when using input orientation (scores ≤ 1), since the bootstrap
distribution reflects sampling variation around the true frontier.


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
{phang2}{cmd:"Firm G" 16 11 110}{p_end}
{phang2}{cmd:"Firm H" 24 16 140}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Basic bootstrap DEA with 500 replications

{phang2}{cmd:. dea_bootstrap capital labor = output, reps(500) seed(12345)}{p_end}

{pstd}VRS with higher confidence level

{phang2}{cmd:. dea_bootstrap capital labor = output, rts(vrs) level(90) reps(1000)}{p_end}

{pstd}String DMU identifiers, save results

{phang2}{cmd:. dea_bootstrap capital labor = output, dmu(firm) reps(500) generate(bs)}{p_end}
{phang2}{cmd:. list firm bs bs_bc bs_se bs_lo bs_hi}{p_end}



{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea_bootstrap} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:r(N)}}number of DMUs{p_end}
{synopt:{cmd:r(reps)}}number of bootstrap replications{p_end}
{synopt:{cmd:r(level)}}confidence level{p_end}
{synopt:{cmd:r(mean_eff)}}mean original efficiency{p_end}
{synopt:{cmd:r(mean_bc)}}mean bias-corrected efficiency{p_end}
{synopt:{cmd:r(mean_se)}}mean bootstrap standard error{p_end}
{syntab:Macros}
{synopt:{cmd:r(rts)}}returns-to-scale assumption{p_end}
{synopt:{cmd:r(orientation)}}model orientation{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker BCC1984}{...}
{phang}
Banker, R. D., A. Charnes, and W. W. Cooper. 1984.
Some models for estimating technical and scale inefficiencies in
data envelopment analysis.
{it:Management Science} 30(9): 1078--1092.
{browse "https://doi.org/10.1287/mnsc.30.9.1078":doi:10.1287/mnsc.30.9.1078}

{marker CCR1978}{...}
{phang}
Charnes, A., W. W. Cooper, and E. Rhodes. 1978.
Measuring the efficiency of decision making units.
{it:European Journal of Operational Research} 2(6): 429--444.
{browse "https://doi.org/10.1016/0377-2217(78)90138-8":doi:10.1016/0377-2217(78)90138-8}

{marker SW1998}{...}
{phang}
Simar, L., and P. W. Wilson. 1998.
Sensitivity analysis of efficiency scores: How to bootstrap in nonparametric
frontier models.
{it:Management Science} 44(1): 49--61.
{browse "https://doi.org/10.1287/mnsc.44.1.49":doi:10.1287/mnsc.44.1.49}


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
