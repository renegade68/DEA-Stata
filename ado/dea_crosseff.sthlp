{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.4.1  21jul2026}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_bootstrap" "help dea_bootstrap"}{...}
{vieweralsosee "dea_malmquist" "help dea_malmquist"}{...}
{vieweralsosee "dea_window" "help dea_window"}{...}
{viewerjumpto "Syntax" "dea_crosseff##syntax"}{...}
{viewerjumpto "Description" "dea_crosseff##description"}{...}
{viewerjumpto "Options" "dea_crosseff##options"}{...}
{viewerjumpto "Remarks" "dea_crosseff##remarks"}{...}
{viewerjumpto "Examples" "dea_crosseff##examples"}{...}
{viewerjumpto "Stored results" "dea_crosseff##stored_results"}{...}
{viewerjumpto "References" "dea_crosseff##references"}{...}
{viewerjumpto "Authors" "dea_crosseff##authors"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:dea_crosseff} {hline 2}}Cross-efficiency DEA{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea_crosseff}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}
[{cmd:,} {it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_crosseff}
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
{synopt:{cmd:method(}{it:methodname}{cmd:)}}secondary objective method; default is {cmd:method(arbitrary)}{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:newvar}{cmd:)}}name for cross-efficiency variable; default is {cmd:_ce_eff}{p_end}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable (string or numeric){p_end}
{synopt:{cmd:matrix}}display full cross-efficiency matrix{p_end}
{synopt:{cmd:detail}}display detailed results{p_end}

{syntab:Technical}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmd:debug}}display debugging information{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:rtstype}: {cmd:crs} | {cmd:vrs}{p_end}
{p 4 6 2}
{it:type}: {cmd:input} | {cmd:output}{p_end}
{p 4 6 2}
{it:methodname}: {cmd:arbitrary} | {cmd:benevolent} | {cmd:aggressive}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_crosseff} computes cross-efficiency scores for DEA.
In standard DEA, each DMU selects multiplier weights that maximize its own
efficiency, which may result in extreme weighting schemes and poor discrimination
among efficient DMUs. Cross-efficiency evaluation uses each DMU's optimal weights
to also evaluate all other DMUs, producing a cross-efficiency matrix.

{pstd}
The mean cross-efficiency score for each DMU — the column average of the
cross-efficiency matrix — provides a peer-evaluation ranking that is generally
better at discriminating among efficient units than self-appraisal DEA scores.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{cmdab:r:ts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
{cmd:crs} assumes constant returns to scale; {cmd:vrs} assumes variable returns
to scale. Default is {cmd:rts(crs)}. The original cross-efficiency formulation
is CRS-based; {cmd:rts(vrs)} implements a later VRS extension.

{phang}
{cmdab:or:ientation(}{it:type}{cmd:)} specifies whether to minimize inputs
({cmd:input}) or maximize outputs ({cmd:output}). Default is
{cmd:orientation(input)}.

{phang}
{cmd:method(}{it:methodname}{cmd:)} specifies the secondary objective when
multiple optimal weight sets exist.
{cmd:method(arbitrary)} uses the standard DEA solution without a secondary
objective ({help dea_crosseff##SSH1986:Sexton et al. 1986}).
{cmd:method(benevolent)} maximizes the cross-efficiency scores of other DMUs
({help dea_crosseff##DG1994:Doyle and Green 1994}).
{cmd:method(aggressive)} minimizes the cross-efficiency scores of other DMUs
({help dea_crosseff##DG1994:Doyle and Green 1994}).
Default is {cmd:method(arbitrary)}.

{dlgtab:Output}

{phang}
{cmdab:gen:erate(}{it:newvar}{cmd:)} specifies the name of the saved mean
cross-efficiency variable. The command also saves {it:newvar}{cmd:_self}
and {it:newvar}{cmd:_maverick}. Default is {cmd:_ce_eff}.

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies a variable identifying each DMU.
Accepts both string and numeric variables.

{phang}
{cmd:matrix} displays the full n×n cross-efficiency matrix.


{marker remarks}{...}
{title:Remarks}

{pstd}
The original cross-efficiency formulation is CRS-based.
{help dea_crosseff##SSH1986:Sexton et al. (1986)} introduced the
cross-efficiency matrix as a peer-evaluation ranking method, and the benevolent
and aggressive secondary objectives ({help dea_crosseff##DG1994:Doyle and Green 1994})
resolve non-uniqueness of the optimal multiplier solution. The command also
allows VRS cross-efficiency evaluation, following later work that extended
cross-efficiency analysis to variable returns to scale
({help dea_crosseff##LZ2015:Lim and Zhu 2015}).


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
{phang2}{cmd:end}{p_end}

{pstd}Basic cross-efficiency with arbitrary secondary objective

{phang2}{cmd:. dea_crosseff capital labor = output, dmu(firm)}{p_end}

{pstd}Benevolent cross-efficiency

{phang2}{cmd:. dea_crosseff capital labor = output, dmu(firm) method(benevolent) matrix}{p_end}

{pstd}VRS cross-efficiency

{phang2}{cmd:. dea_crosseff capital labor = output, dmu(firm) rts(vrs)}{p_end}

{pstd}String DMU identifiers

{phang2}{cmd:. dea_crosseff capital labor = output, dmu(firm) method(aggressive)}{p_end}



{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea_crosseff} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:r(n)}}number of DMUs{p_end}
{synopt:{cmd:r(mean_ce)}}mean cross-efficiency across all DMUs{p_end}
{synopt:{cmd:r(mean_self)}}mean self-appraisal score{p_end}
{syntab:Macros}
{synopt:{cmd:r(method)}}secondary objective used{p_end}
{synopt:{cmd:r(rts)}}returns-to-scale assumption{p_end}
{synopt:{cmd:r(orientation)}}model orientation{p_end}
{syntab:Matrices}
{synopt:{cmd:r(ce_matrix)}}cross-efficiency matrix; stored when {cmd:matrix} is specified{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker DG1994}{...}
{phang}
Doyle, J., and R. Green. 1994.
Efficiency and cross-efficiency in DEA: Derivations, meanings and uses.
{it:Journal of the Operational Research Society} 45(5): 567--578.
{browse "https://doi.org/10.1057/jors.1994.84":doi:10.1057/jors.1994.84}

{marker LZ2015}{...}
{phang}
Lim, S., and J. Zhu. 2015.
DEA cross-efficiency evaluation under variable returns to scale.
{it:Journal of the Operational Research Society} 66(3): 476--487.
{browse "https://doi.org/10.1057/jors.2014.13":doi:10.1057/jors.2014.13}

{marker SSH1986}{...}
{phang}
Sexton, T. R., R. H. Silkman, and A. J. Hogan. 1986.
Data envelopment analysis: Critique and extensions.
{it:New Directions for Program Evaluation} 1986(32): 73--105.
{browse "https://doi.org/10.1002/ev.1441":doi:10.1002/ev.1441}


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
