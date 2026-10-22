{smcl}
{* Help revision 2026-10-22.2: duplicate author email removed; self-contained examples retained.}{...}
{* *! version 2.0.0  26apr2026}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_sbmundesirable" "help dea_sbmundesirable"}{...}
{viewerjumpto "Syntax" "dea_congestion##syntax"}{...}
{viewerjumpto "Description" "dea_congestion##description"}{...}
{viewerjumpto "Options" "dea_congestion##options"}{...}
{viewerjumpto "Remarks" "dea_congestion##remarks"}{...}
{viewerjumpto "Examples" "dea_congestion##examples"}{...}
{viewerjumpto "Stored results" "dea_congestion##stored_results"}{...}
{viewerjumpto "References" "dea_congestion##references"}{...}
{viewerjumpto "Authors" "dea_congestion##authors"}{...}
{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{cmd:dea_congestion} {hline 2}}Congestion DEA (Cooper et al. 2004){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea_congestion}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}
[{cmd:,} {it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_congestion}
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

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for generated variables; default is {cmd:_cong}{p_end}

{syntab:Technical}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable{p_end}
{synopt:{cmd:engine(}{it:string}{cmd:)}}solver engine; default is {cmd:engine(ipm)}{p_end}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_congestion} identifies input congestion using the two-step procedure
of Cooper et al. (2004).  Congestion occurs when an input is used
in excess --- reducing the input would actually {it:increase} output --- a
perverse input-output relationship that standard DEA efficiency scores conceal.

{pstd}
The command implements the Cooper et al. (2004) procedure:

{phang2}
Step 1.  Standard VRS (BCC) DEA is solved under strong disposability,
yielding the VRS technical efficiency score {it:theta_VRS}.

{phang2}
Step 2.  A modified LP imposes weak input disposability through input
equality constraints, with outputs also held at equality, yielding the
congestion-model efficiency {it:theta_C}.

{pstd}
The congestion measure is {it:theta_VRS / theta_C}.  A value of 1 indicates
no congestion; values below 1 indicate congestion.

{pstd}
Five variables are created with the specified prefix:

{p2colset 9 32 34 2}
{p2col:{it:prefix}{cmd:_crs}}CRS technical efficiency{p_end}
{p2col:{it:prefix}{cmd:_vrs}}VRS technical efficiency{p_end}
{p2col:{it:prefix}{cmd:_cong_eff}}congestion-model efficiency ({it:theta_C}){p_end}
{p2col:{it:prefix}{cmd:_congestion}}congestion measure ({it:theta_VRS / theta_C}){p_end}
{p2col:{it:prefix}{cmd:_se}}scale efficiency ({it:theta_CRS / theta_VRS}){p_end}
{p2colreset}{...}

{pstd}
{cmd:dea_congestion} is a convenience wrapper around
{cmd:dea ..., model(congestion)}.


{marker options}{...}
{title:Options}

{dlgtab:Output}

{phang}
{cmdab:gen:erate(}{it:prefix}{cmd:)} specifies the prefix for all generated
variables.  Default is {cmd:_cong}.

{dlgtab:Technical}

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies the variable identifying each DMU.
Accepts both string and numeric variables.

{phang}
{cmd:engine(}{it:string}{cmd:)} specifies the solver engine.
Currently only {cmd:ipm} (Interior Point Method) is supported.

{phang}
{cmdab:tol:erance(}{it:#}{cmd:)} sets the numerical tolerance for the LP solver.
Default is {cmd:1e-7}.


{marker remarks}{...}
{title:Remarks}

{pstd}
Congestion arises when continued application of an input reduces rather
than increases output.  This can occur in situations of resource overuse:
too many workers on a fixed production line, excess fertilizer on a crop,
or overcrowding in a facility.

{pstd}
Standard radial DEA efficiency scores are insensitive to congestion
because they assume strong disposability of inputs --- an implicit
assumption that any input can be freely reduced without harm.
The Cooper et al. (2004) method detects congestion by comparing
efficiency under strong versus weak output disposability.

{pstd}
Mathematically, if the output equality constraint (Step 2) produces
a higher efficiency than the inequality constraint (Step 1), then
the DMU is producing {it:less} than it could at its current input
level, which is the signature of congestion.

{pstd}
The model requires Stata 16 or later (for the {cmd:LinearProgram()}
class used in the interior-point solver).


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
{phang2}{cmd:. input str10 hospital beds nurses patients}{p_end}
{phang2}{cmd:"Hospital A" 120 55 5100}{p_end}
{phang2}{cmd:"Hospital B" 150 70 6200}{p_end}
{phang2}{cmd:"Hospital C" 180 95 6400}{p_end}
{phang2}{cmd:"Hospital D" 110 50 5000}{p_end}
{phang2}{cmd:"Hospital E" 210 120 6600}{p_end}
{phang2}{cmd:"Hospital F" 160 75 6100}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Basic congestion analysis

{phang2}{cmd:. dea_congestion beds nurses = patients, dmu(hospital)}{p_end}

{pstd}Form B syntax

{phang2}{cmd:. dea_congestion, inputs(beds nurses) outputs(patients) dmu(hospital)}{p_end}

{pstd}Custom prefix for output variables

{phang2}{cmd:. dea_congestion beds nurses = patients, generate(hosp) dmu(hospital)}{p_end}
{phang2}{cmd:. list hospital hosp_vrs hosp_congestion if hosp_congestion < 0.99}{p_end}

{pstd}Equivalently using dea directly

{phang2}{cmd:. dea beds nurses = patients, model(congestion) generate(hosp) dmu(hospital)}{p_end}



{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea_congestion} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:r(N)}}number of DMUs{p_end}
{synopt:{cmd:r(n_congested)}}number of congested DMUs{p_end}
{synopt:{cmd:r(mean_vrs)}}mean VRS technical efficiency{p_end}
{synopt:{cmd:r(mean_efficiency)}}mean congestion measure{p_end}
{synopt:{cmd:r(n_efficient)}}number of non-congested DMUs{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{phang}
Cooper, W. W., H. Deng, L. M. Seiford, and J. Zhu. 2004.
Congestion: Its identification and management with DEA.
In {it:Handbook on Data Envelopment Analysis},
ed. W. W. Cooper, L. M. Seiford, and J. Zhu, 177--201.
Boston: Kluwer Academic Publishers.

{phang}
Cooper, W. W., L. M. Seiford, and K. Tone. 2007.
{it:Data Envelopment Analysis: A Comprehensive Text with Models, Applications,
References and DEA-Solver Software}. 2nd ed. New York: Springer.


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
