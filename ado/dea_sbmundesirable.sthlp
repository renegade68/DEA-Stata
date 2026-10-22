{smcl}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* *! version 1.0.1  02may2026}{...}
{viewerjumpto "Related article" "dea_sbmundesirable##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_congestion" "help dea_congestion"}{...}
{viewerjumpto "Syntax" "dea_sbmundesirable##syntax"}{...}
{viewerjumpto "Description" "dea_sbmundesirable##description"}{...}
{viewerjumpto "Options" "dea_sbmundesirable##options"}{...}
{viewerjumpto "Remarks" "dea_sbmundesirable##remarks"}{...}
{viewerjumpto "Examples" "dea_sbmundesirable##examples"}{...}
{viewerjumpto "Stored results" "dea_sbmundesirable##stored_results"}{...}
{viewerjumpto "References" "dea_sbmundesirable##references"}{...}
{viewerjumpto "Authors" "dea_sbmundesirable##authors"}{...}
{title:Title}

{p2colset 5 28 30 2}{...}
{p2col :{cmd:dea_sbmundesirable} {hline 2}}SBM DEA with undesirable outputs (Tone 2004){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea_sbmundesirable}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}{cmd:,}
{cmd:bad(}{it:varlist}{cmd:)}
[{it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_sbmundesirable}
{ifin}{cmd:,}
{cmdab:in:puts(}{it:varlist}{cmd:)}
{cmdab:out:puts(}{it:varlist}{cmd:)}
{cmd:bad(}{it:varlist}{cmd:)}
[{it:options}]


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmdab:in:puts(}{it:varlist}{cmd:)}}input variable(s){p_end}
{synopt:{cmdab:out:puts(}{it:varlist}{cmd:)}}desirable (good) output variable(s){p_end}
{synopt:{cmd:bad(}{it:varlist}{cmd:)}}undesirable (bad) output variable(s){p_end}

{syntab:Model}
{synopt:{cmdab:r:ts(}{it:rtstype}{cmd:)}}returns to scale; default is {cmd:rts(crs)}{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:name}{cmd:)}}name of generated efficiency variable; default is {cmd:sbmu_eff}{p_end}

{syntab:Technical}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable{p_end}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:rtstype}: {cmd:crs} | {cmd:vrs}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_sbmundesirable} computes efficiency scores for production processes
that simultaneously generate desirable (good) outputs and undesirable (bad)
outputs, using the Slacks-Based Measure (SBM) framework of
{help dea_sbmundesirable##TONE2004:Tone (2004)}.

{pstd}
Undesirable outputs --- such as pollutants, waste, or by-products --- cannot
be treated like regular outputs, because more is not necessarily better.
Standard DEA models that include bad outputs as ordinary outputs would
reward DMUs for producing more pollution.
{cmd:dea_sbmundesirable} handles this by imposing a {it:weak disposability}
(equality) constraint on bad outputs, reflecting the cost of reducing them.

{pstd}
The model solves the following fractional program
({help dea_sbmundesirable##TONE2004:Tone 2004}):

{p 8 12 2}
min  (1 - (1/m)*Sum(s_i^-/x_i0)) / (1 + (1/(s+b))*(Sum(s_r^+/y_r0) + Sum(s_k^b/b_k0)))

{pstd}subject to:

{p2colset 9 40 42 2}
{p2col:{it:X lambda + s^- = x_0}}input constraints{p_end}
{p2col:{it:Y lambda - s^+ = y_0}}good output constraints (free disposal){p_end}
{p2col:{it:B lambda + s^b = b_0}}bad output constraints (weak disposability){p_end}
{p2col:{it:Sum(lambda) = 1}}convexity (VRS only){p_end}
{p2col:{it:lambda, s^-, s^+, s^b >= 0}}non-negativity{p_end}
{p2colreset}{...}

{pstd}
where m is the number of inputs, s is the number of good outputs, and b is
the number of bad outputs.  The model is solved via the
Charnes-Cooper linearization and the Interior Point Method engine.

{pstd}
{cmd:dea_sbmundesirable} is a convenience wrapper around
{cmd:dea ..., model(sbmundesirable) bad(...)}.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{cmd:bad(}{it:varlist}{cmd:)} specifies the undesirable (bad) output variables.
These are treated under weak disposability: reducing bad outputs requires
proportional reduction in inputs or good outputs.

{dlgtab:Model}

{phang}
{cmdab:r:ts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
{cmd:crs} (constant, default) or {cmd:vrs} (variable).

{dlgtab:Output}

{phang}
{cmdab:gen:erate(}{it:name}{cmd:)} specifies the name for the generated
efficiency variable.  Default is {cmd:sbmu_eff}.

{dlgtab:Technical}

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies the DMU identifier variable.

{phang}
{cmdab:tol:erance(}{it:#}{cmd:)} sets the numerical tolerance.
Default is {cmd:1e-7}.


{marker remarks}{...}
{title:Remarks}

{pstd}
The SBM-Undesirable model extends the standard SBM efficiency measure
({help dea:dea model(sbm)}; {help dea_sbmundesirable##TONE2001:Tone 2001})
to accommodate pollution and other
undesirable by-products.  The key distinction from
{cmd:dea model(undesirable)}, which uses the Directional Distance Function
(DDF) approach of {help dea_sbmundesirable##CCF1998:Chambers, Chung, and Färe (1998)},
is that
{cmd:dea_sbmundesirable} uses a non-radial, slacks-based objective that
accounts for {it:all} input and output inefficiencies simultaneously.

{pstd}
Efficiency scores range from 0 to 1, with 1 indicating full efficiency.
A score below 1 reflects slack in inputs, underproduction of good outputs,
or overproduction of bad outputs.

{pstd}
The Charnes-Cooper linearization transforms the fractional program into
a linear program by substituting {it:t = 1/denominator} and
{it:z = t * lambda}.  The model requires Stata 16 or later.

{pstd}
For the DDF-based approach to undesirable outputs, see
{cmd:dea model(undesirable)} in {helpb dea}.


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


{pstd}Setup: power plants producing electricity (good) and SO2 (bad)

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str7 plant capital labor fuel electricity so2}{p_end}
{phang2}{cmd:"Plant A" 120 45 210  95 32}{p_end}
{phang2}{cmd:"Plant B" 160 52 260 125 45}{p_end}
{phang2}{cmd:"Plant C" 140 48 235 118 38}{p_end}
{phang2}{cmd:"Plant D" 190 60 300 138 58}{p_end}
{phang2}{cmd:"Plant E" 110 42 205  92 30}{p_end}
{phang2}{cmd:"Plant F" 175 56 275 132 50}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}SBM efficiency with undesirable output

{phang2}{cmd:. dea_sbmundesirable capital labor fuel = electricity, bad(so2) dmu(plant)}{p_end}

{pstd}VRS assumption

{phang2}{cmd:. dea_sbmundesirable capital labor fuel = electricity, bad(so2) rts(vrs) dmu(plant)}{p_end}

{pstd}Form B syntax with custom variable name

{phang2}{cmd:. dea_sbmundesirable, inputs(capital labor fuel) outputs(electricity) bad(so2) ///}{p_end}
{phang2}{cmd:    rts(vrs) generate(green_eff) dmu(plant)}{p_end}
{phang2}{cmd:. list plant green_eff}{p_end}

{pstd}Equivalently using dea directly

{phang2}{cmd:. dea capital labor fuel = electricity, model(sbmundesirable) bad(so2) dmu(plant)}{p_end}



{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea_sbmundesirable} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:r(N)}}number of DMUs{p_end}
{synopt:{cmd:r(N_eff)}}number of efficient DMUs, inherited from {cmd:dea}{p_end}
{synopt:{cmd:r(mean_eff)}}mean efficiency score, inherited from {cmd:dea}{p_end}
{synopt:{cmd:r(n_efficient)}}number of efficient DMUs (score = 1){p_end}
{synopt:{cmd:r(mean_efficiency)}}mean efficiency score{p_end}
{synopt:{cmd:r(min_eff)}}minimum efficiency score{p_end}
{synopt:{cmd:r(max_eff)}}maximum efficiency score{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker TONE2004}{...}
{phang}
Tone, K. 2004.
Dealing with undesirable outputs in DEA: A Slacks-Based Measure (SBM) approach.
GRIPS Discussion Paper I-2003-0005.
{browse "https://doi.org/10.24545/00000955":doi:10.24545/00000955}

{marker TONE2001}{...}
{phang}
Tone, K. 2001.
A slacks-based measure of efficiency in data envelopment analysis.
{it:European Journal of Operational Research} 130(3): 498--509.
{browse "https://doi.org/10.1016/S0377-2217(99)00407-5":doi:10.1016/S0377-2217(99)00407-5}

{marker CCF1998}{...}
{phang}
Chambers, R. G., Y. Chung, and R. Färe. 1998.
Profit, directional distance functions, and Nerlovian efficiency.
{it:Journal of Optimization Theory and Applications} 98(2): 351--364.
{browse "https://doi.org/10.1023/A:1022637501082":doi:10.1023/A:1022637501082}

{marker CST2007}{...}
{phang}
Cooper, W. W., L. M. Seiford, and K. Tone. 2007.
{it:Data Envelopment Analysis: A Comprehensive Text with Models, Applications,
References and DEA-Solver Software}. 2nd ed. New York: Springer.
{browse "https://doi.org/10.1007/978-0-387-45283-8":doi:10.1007/978-0-387-45283-8}


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
