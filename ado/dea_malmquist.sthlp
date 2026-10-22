{smcl}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* *! version 1.3.1  28apr2026}{...}
{viewerjumpto "Related article" "dea_malmquist##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_bootstrap" "help dea_bootstrap"}{...}
{vieweralsosee "dea_window" "help dea_window"}{...}
{viewerjumpto "Syntax" "dea_malmquist##syntax"}{...}
{viewerjumpto "Description" "dea_malmquist##description"}{...}
{viewerjumpto "Options" "dea_malmquist##options"}{...}
{viewerjumpto "Remarks" "dea_malmquist##remarks"}{...}
{viewerjumpto "Examples" "dea_malmquist##examples"}{...}
{viewerjumpto "Stored results" "dea_malmquist##stored_results"}{...}
{viewerjumpto "References" "dea_malmquist##references"}{...}
{viewerjumpto "Authors" "dea_malmquist##authors"}{...}
{title:Title}

{p2colset 5 23 25 2}{...}
{p2col :{cmd:dea_malmquist} {hline 2}}Malmquist productivity index{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea_malmquist}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}{cmd:,}
{cmdab:dmu(}{it:varname}{cmd:)}
{cmdab:time(}{it:varname}{cmd:)}
[{it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_malmquist}
{ifin}{cmd:,}
{cmdab:dmu(}{it:varname}{cmd:)}
{cmdab:time(}{it:varname}{cmd:)}
{cmdab:in:puts(}{it:varlist}{cmd:)}
{cmdab:out:puts(}{it:varlist}{cmd:)}
[{it:options}]


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable (string or numeric){p_end}
{synopt:{cmd:time(}{it:varname}{cmd:)}}time period variable (numeric){p_end}
{synopt:{cmdab:in:puts(}{it:varlist}{cmd:)}}input variable(s){p_end}
{synopt:{cmdab:out:puts(}{it:varlist}{cmd:)}}output variable(s){p_end}

{syntab:Model}
{synopt:{cmdab:r:ts(}{it:rtstype}{cmd:)}}returns to scale; default is {cmd:rts(crs)}{p_end}
{synopt:{cmdab:or:ientation(}{it:type}{cmd:)}}orientation; default is {cmd:orientation(input)}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}
{synopt:{cmd:global}}use a pooled-frontier global Malmquist index{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for result variables{p_end}
{synopt:{cmd:detail}}accepted for compatibility; DMU-level results are displayed by default{p_end}

{syntab:Technical}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmd:debug}}display debugging information{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:rtstype}: {cmd:crs} | {cmd:vrs}{p_end}
{p 4 6 2}
{it:type}: {cmd:input} | {cmd:output}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_malmquist} computes the Malmquist productivity index (MPI) for panel
data. The MPI decomposes total factor productivity change between consecutive
time periods into two components:

{phang2}
{bf:Efficiency change (EC)}: The catch-up effect — how much closer a DMU moved
to the production frontier between periods t and t+1. EC > 1 indicates
efficiency improvement.

{phang2}
{bf:Technical change (TC)}: The frontier-shift effect — how much the best-practice
frontier itself moved between periods. TC > 1 indicates technological progress.

{pstd}
The overall Malmquist productivity index is MPI = EC × TC.
MPI > 1 indicates productivity growth; MPI < 1 indicates productivity decline.

{pstd}
The command computes four distance functions for each adjacent pair of periods
(t, t+1): D^t(x^t, y^t), D^{t+1}(x^{t+1}, y^{t+1}), D^t(x^{t+1}, y^{t+1}),
and D^{t+1}(x^t, y^t). The geometric mean of the two technology ratio components
constitutes the technical change index.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies the variable that identifies each
decision-making unit. Accepts both string and numeric variables.

{phang}
{cmd:time(}{it:varname}{cmd:)} specifies the numeric variable indicating the
time period. The command processes consecutive pairs of periods.

{dlgtab:Model}

{phang}
{cmdab:r:ts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
{cmd:crs} assumes constant returns to scale; {cmd:vrs} assumes variable returns
to scale. Default is {cmd:rts(crs)}.

{phang}
{cmdab:or:ientation(}{it:type}{cmd:)} specifies the direction of optimization.
Default is {cmd:orientation(input)}.


{marker remarks}{...}
{title:Remarks}

{pstd}
The data must be in long format with one observation per DMU per period.
Only DMUs that appear in both consecutive periods are included in the
Malmquist calculations for that period pair.

{pstd}
The geometric mean formulation of
{help dea_malmquist##CCD1982:Caves, Christensen, and Diewert (1982)} is used,
which avoids the need to choose a reference technology. The DEA decomposition
into efficiency change and technical change follows
{help dea_malmquist##FGNZ1994:Färe et al. (1994)}.

{pstd}
With {cmd:global}, the command estimates a single pooled frontier from all
periods and computes a global Malmquist index. This is useful when adjacent
period frontiers are infeasible or when a circular productivity index is desired.


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
{phang2}{cmd:. input str10 country year expenditure teachers graduation}{p_end}
{phang2}{cmd:"Aland"   2020 120 30  88}{p_end}
{phang2}{cmd:"Bordia"  2020 150 38 105}{p_end}
{phang2}{cmd:"Cyrenia" 2020 135 34  98}{p_end}
{phang2}{cmd:"Dameria" 2020 170 45 118}{p_end}
{phang2}{cmd:"Eldora"  2020 110 28  82}{p_end}
{phang2}{cmd:"Aland"   2021 125 31  94}{p_end}
{phang2}{cmd:"Bordia"  2021 154 39 112}{p_end}
{phang2}{cmd:"Cyrenia" 2021 138 35 103}{p_end}
{phang2}{cmd:"Dameria" 2021 176 46 126}{p_end}
{phang2}{cmd:"Eldora"  2021 113 29  87}{p_end}
{phang2}{cmd:"Aland"   2022 130 32 101}{p_end}
{phang2}{cmd:"Bordia"  2022 160 40 119}{p_end}
{phang2}{cmd:"Cyrenia" 2022 143 36 108}{p_end}
{phang2}{cmd:"Dameria" 2022 182 47 134}{p_end}
{phang2}{cmd:"Eldora"  2022 116 30  92}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Malmquist index with CRS

{phang2}{cmd:. dea_malmquist expenditure teachers = graduation, dmu(country) time(year)}{p_end}

{pstd}VRS with string DMU identifiers

{phang2}{cmd:. dea_malmquist expenditure teachers = graduation, dmu(country) time(year) rts(vrs)}{p_end}

{pstd}Global pooled-frontier Malmquist index

{phang2}{cmd:. dea_malmquist expenditure teachers = graduation, dmu(country) time(year) global}{p_end}

{pstd}Save decomposition results

{phang2}{cmd:. dea_malmquist expenditure teachers = graduation, dmu(country) time(year) generate(mpi)}{p_end}
{phang2}{cmd:. list country year mpi_ec mpi_tc mpi_mpi}{p_end}



{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea_malmquist} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:r(n_dmus)}}number of common DMUs in the last period pair processed{p_end}
{synopt:{cmd:r(n_periods)}}number of time periods{p_end}
{synopt:{cmd:r(mean_ec)}}mean efficiency change{p_end}
{synopt:{cmd:r(mean_tc)}}mean technical change{p_end}
{synopt:{cmd:r(mean_mpi)}}mean Malmquist productivity index{p_end}
{syntab:Macros}
{synopt:{cmd:r(malmquist_type)}}{cmd:adjacent} or {cmd:global}{p_end}
{synopt:{cmd:r(rts)}}returns to scale used{p_end}
{synopt:{cmd:r(orientation)}}orientation used{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker CCD1982}{...}
{phang}
Caves, D. W., L. R. Christensen, and W. E. Diewert. 1982.
The economic theory of index numbers and the measurement of input, output,
and productivity.
{it:Econometrica} 50(6): 1393--1414.
{browse "https://doi.org/10.2307/1913388":doi:10.2307/1913388}

{marker FGNZ1994}{...}
{phang}
Färe, R., S. Grosskopf, M. Norris, and Z. Zhang. 1994.
Productivity growth, technical progress, and efficiency change in
industrialized countries.
{it:American Economic Review} 84(1): 66--83.


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
