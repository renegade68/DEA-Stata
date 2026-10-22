{smcl}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* *! version 1.1.2  03may2026}{...}
{viewerjumpto "Related article" "dea_network##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_bootstrap" "help dea_bootstrap"}{...}
{vieweralsosee "dea_crosseff" "help dea_crosseff"}{...}
{vieweralsosee "dea_malmquist" "help dea_malmquist"}{...}
{vieweralsosee "dea_window" "help dea_window"}{...}
{viewerjumpto "Syntax" "dea_network##syntax"}{...}
{viewerjumpto "Description" "dea_network##description"}{...}
{viewerjumpto "Options" "dea_network##options"}{...}
{viewerjumpto "Remarks" "dea_network##remarks"}{...}
{viewerjumpto "Examples" "dea_network##examples"}{...}
{viewerjumpto "Stored results" "dea_network##stored_results"}{...}
{viewerjumpto "References" "dea_network##references"}{...}
{viewerjumpto "Authors" "dea_network##authors"}{...}
{title:Title}

{p2colset 5 21 23 2}{...}
{p2col :{cmd:dea_network} {hline 2}}Network DEA (two-stage model){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea_network}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}{cmd:,}
{cmdab:int:ermediates(}{it:varlist}{cmd:)}
[{it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_network}
{ifin}{cmd:,}
{cmdab:in:puts(}{it:varlist}{cmd:)}
{cmdab:int:ermediates(}{it:varlist}{cmd:)}
{cmdab:out:puts(}{it:varlist}{cmd:)}
[{it:options}]


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmdab:in:puts(}{it:varlist}{cmd:)}}Stage 1 input variables (X){p_end}
{synopt:{cmdab:int:ermediates(}{it:varlist}{cmd:)}}intermediate variables linking Stage 1 and Stage 2 (Z){p_end}
{synopt:{cmdab:out:puts(}{it:varlist}{cmd:)}}Stage 2 output variables (Y){p_end}

{syntab:Model}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable (string or numeric){p_end}
{synopt:{cmdab:r:ts(}{it:rtstype}{cmd:)}}returns to scale; default is {cmd:rts(crs)}{p_end}
{synopt:{cmdab:or:ientation(}{it:type}{cmd:)}}orientation; default is {cmd:orientation(input)}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}
{synopt:{cmd:method(}{it:type}{cmd:)}}decomposition method; default is {cmd:method(multiplicative)}{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for generated efficiency variables; default is {cmd:_net_eff}{p_end}
{synopt:{cmd:detail}}accepted for compatibility; stage-by-stage results are displayed by default{p_end}

{syntab:Technical}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmd:debug}}display debugging information{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:rtstype}: {cmd:crs} | {cmd:vrs}{p_end}
{p 4 6 2}
{it:type} (orientation): {cmd:input} | {cmd:output}{p_end}
{p 4 6 2}
{it:type} (method): {cmd:multiplicative} | {cmd:additive} | {cmd:independent}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_network} computes Network DEA efficiency scores for two-stage production
processes where intermediate products connect the two stages.  The two-stage structure is:

{p 8 17 2}
Stage 1: X {hline 2}> Z  (external inputs to intermediate outputs){break}
Stage 2: Z {hline 2}> Y  (intermediate inputs to final outputs)

{pstd}
Overall efficiency is decomposed into stage-level efficiencies, which allows
identification of which stage is driving overall performance.

{pstd}
Three generated variables are created with the specified prefix:

{p2colset 9 32 34 2}
{p2col:{it:prefix}{cmd:_overall}}overall efficiency score{p_end}
{p2col:{it:prefix}{cmd:_stage1}}Stage 1 efficiency (X {hline 2}> Z){p_end}
{p2col:{it:prefix}{cmd:_stage2}}Stage 2 efficiency (Z {hline 2}> Y){p_end}
{p2colreset}{...}


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{cmdab:in:puts(}{it:varlist}{cmd:)} specifies the Stage 1 input variables (X).

{phang}
{cmdab:int:ermediates(}{it:varlist}{cmd:)} specifies the intermediate variables (Z)
that are outputs of Stage 1 and inputs to Stage 2.

{phang}
{cmdab:out:puts(}{it:varlist}{cmd:)} specifies the Stage 2 output variables (Y).

{dlgtab:Model}

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies the variable identifying each DMU.
Accepts both string and numeric variables.  If not specified, observation
numbers are used.

{phang}
{cmdab:r:ts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
Allowed values are {cmd:crs} (constant) and {cmd:vrs} (variable).
Default is {cmd:rts(crs)}.

{phang}
{cmdab:or:ientation(}{it:type}{cmd:)} specifies the direction of optimization.
Default is {cmd:orientation(input)}.

{phang}
{cmd:method(}{it:type}{cmd:)} specifies how overall efficiency is decomposed
from stage efficiencies:

{p2colset 9 26 28 2}
{p2col:{cmd:multiplicative}}Kao-Hwang shared intermediate-weight relational
multiplier model under {cmd:rts(crs)} and {cmd:orientation(input)}; default{p_end}
{p2col:{cmd:additive}}overall = 0.5 {it:x} Stage1 + 0.5 {it:x} Stage2 (equal weights){p_end}
{p2col:{cmd:independent}}stages computed independently with no linking constraint;
overall = Stage1 {it:x} Stage2{p_end}
{p2colreset}{...}

{dlgtab:Output}

{phang}
{cmdab:gen:erate(}{it:prefix}{cmd:)} specifies the prefix for generated variables.
Default is {cmd:_net_eff}.

{phang}
{cmd:detail} is accepted for compatibility; stage-level output for each DMU is
displayed by default.


{marker remarks}{...}
{title:Remarks}

{pstd}
Network DEA generalizes the standard DEA model by explicitly modeling
the internal structure of the production process
({help dea_network##FG2000:Färe and Grosskopf 2000}). Whereas standard DEA
treats each DMU as a black box that converts inputs into outputs
({help dea_network##CST2007:Cooper, Seiford, and Tone 2007}),
{cmd:dea_network} models the intermediate flow of products between stages.

{pstd}
The {cmd:multiplicative} method implements the shared intermediate-weight
relational multiplier model of {help dea_network##KH2008:Kao and Hwang (2008)}
for CRS input-oriented two-stage systems.  The same multipliers on intermediate
measures are used as Stage 1 output weights and Stage 2 input weights, and the
reported decomposition satisfies
E{sub:overall} = E{sub:1} {it:x} E{sub:2}.  For VRS or output-oriented calls,
the command retains the legacy stage-wise fallback.

{pstd}
The {cmd:additive} method assigns equal weight to each stage.
The {cmd:independent} method is equivalent to solving two separate
standard DEA problems and taking the product, with no linking constraint
on the intermediate products.

{pstd}
All models are estimated using the Interior Point Method engine via
Stata's {cmd:LinearProgram()} class, requiring Stata 16 or later.


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


{pstd}Setup: bank data where deposits are processed into loans (Stage 1),
which are then converted into interest income (Stage 2).

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input byte bankid deposits staff loans nonloans interest_income profit}{p_end}
{phang2}{cmd:1 120 18  95 28  8.5 12.0}{p_end}
{phang2}{cmd:2 160 24 130 34 11.4 15.5}{p_end}
{phang2}{cmd:3 140 21 118 31 10.2 14.0}{p_end}
{phang2}{cmd:4 190 30 145 40 12.0 16.2}{p_end}
{phang2}{cmd:5 110 16  88 25  7.9 11.1}{p_end}
{phang2}{cmd:6 175 27 138 37 11.8 15.9}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Basic two-stage network DEA (Form A syntax)

{phang2}{cmd:. dea_network deposits = interest_income, intermediates(loans) dmu(bankid)}{p_end}

{pstd}Form B syntax with explicit options

{phang2}{cmd:. dea_network, inputs(deposits) intermediates(loans) outputs(interest_income) dmu(bankid) rts(vrs) method(multiplicative)}{p_end}

{pstd}Multiple inputs and intermediates

{phang2}{cmd:. dea_network deposits staff = profit, intermediates(loans nonloans) dmu(bankid) method(additive)}{p_end}

{pstd}Save results with custom prefix

{phang2}{cmd:. dea_network deposits = interest_income, intermediates(loans) dmu(bankid) generate(bank)}{p_end}
{phang2}{cmd:. list bankid bank_stage1 bank_stage2 bank_overall}{p_end}



{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea_network} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:r(n)}}number of DMUs{p_end}
{synopt:{cmd:r(mean_overall)}}mean overall efficiency{p_end}
{synopt:{cmd:r(mean_stage1)}}mean Stage 1 efficiency{p_end}
{synopt:{cmd:r(mean_stage2)}}mean Stage 2 efficiency{p_end}

{syntab:Macros}
{synopt:{cmd:r(method)}}decomposition method used{p_end}
{synopt:{cmd:r(rts)}}returns to scale used{p_end}
{synopt:{cmd:r(orientation)}}orientation used{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker FG2000}{...}
{phang}
Färe, R., and S. Grosskopf. 2000.
Network DEA.
{it:Socio-Economic Planning Sciences} 34(1): 35--49.
{browse "https://doi.org/10.1016/S0038-0121(99)00012-9":doi:10.1016/S0038-0121(99)00012-9}

{marker KH2008}{...}
{phang}
Kao, C., and S.-N. Hwang. 2008.
Efficiency decomposition in two-stage data envelopment analysis: An
application to non-life insurance companies in Taiwan.
{it:European Journal of Operational Research} 185(1): 418--429.
{browse "https://doi.org/10.1016/j.ejor.2006.11.041":doi:10.1016/j.ejor.2006.11.041}

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
