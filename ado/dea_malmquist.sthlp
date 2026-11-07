{smcl}
{* Help revision 2026-11-02: validated panel LPs, pair summaries and missing components.}{...}
{* Help revision 2026-11-01: validate selected DMU-period keys before output replacement.}{...}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* *! version 1.3.7  02nov2026}{...}
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
[{cmdab:dmu(}{it:varname}{cmd:)}]
{cmdab:time(}{it:varname}{cmd:)}
[{it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_malmquist}
{ifin}{cmd:,}
[{cmdab:dmu(}{it:varname}{cmd:)}]
{cmdab:time(}{it:varname}{cmd:)}
{cmdab:in:puts(}{it:varlist}{cmd:)}
{cmdab:out:puts(}{it:varlist}{cmd:)}
[{it:options}]


{pstd}
Both forms support {cmd:if} and {cmd:in}. Equality tests such as
{cmd:if group==1}, quoted strings, and functions such as
{cmd:if inlist(group,1,2)} are parsed as sample restrictions.
The restrictions apply to both evaluated observations and the reference data.
Input and output varlists may use standard variable ranges and wildcards.

{pstd}
For compatibility, output variables may also precede {cmd:if}/{cmd:in}, with
{cmd:inputs()} after the comma. If that positional varlist is omitted,
{cmd:outputs()} is required; unrelated variables are not inferred as outputs.
Input and output variables must be numeric; {cmd:dmu()} may be a string.

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmd:time(}{it:varname}{cmd:)}}time period variable (numeric){p_end}
{synopt:{cmdab:in:puts(}{it:varlist}{cmd:)}}input variable(s){p_end}
{synopt:{cmdab:out:puts(}{it:varlist}{cmd:)}}output variable(s){p_end}

{syntab:Model}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}optional stable DMU identifier (string or numeric); see remarks for positional IDs{p_end}
{synopt:{cmdab:r:ts(}{it:rtstype}{cmd:)}}returns to scale; default is {cmd:rts(crs)}{p_end}
{synopt:{cmdab:or:ientation(}{it:type}{cmd:)}}orientation; default is {cmd:orientation(input)}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}
{synopt:{cmd:global}}use a pooled-frontier global Malmquist index{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for result variables{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}
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
The command uses the requested RTS technology for every distance. With
{cmd:rts(vrs)}, the result is a VRS distance-based index; it is not the CRS TFP
index decomposed into pure efficiency, scale efficiency, and technical change.
Those additional decompositions are not computed. The selected technology
therefore matters for interpretation as well as feasibility.

{pstd}
The command computes four distance functions for each adjacent pair of periods
(t, t+1): D^t(x^t, y^t), D^{t+1}(x^{t+1}, y^{t+1}), D^t(x^{t+1}, y^{t+1}),
and D^{t+1}(x^t, y^t). The geometric mean of the two technology ratio components
constitutes the technical change index.

{pstd}
The displayed {cmd:E_tt}, {cmd:E_t1t1}, {cmd:E_t_t1}, and {cmd:E_t1_t} are
oriented efficiency measures: {it:theta} for input orientation and {it:1/phi}
for output orientation. The first period subscript identifies the reference
technology, and the second identifies the observation. Input {it:theta} is the
reciprocal of the input Shephard distance. Ratios are formed so that improvement
corresponds to values above 1 in either orientation. Cross-period efficiency
can exceed 1; the output expansion factor may be below 1.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{cmd:time(}{it:varname}{cmd:)} specifies the numeric variable indicating the
time period. The command processes consecutive pairs of periods.

{dlgtab:Model}

{phang}
{cmd:dmu(}{it:varname}{cmd:)} optionally specifies a stable identifier for each
DMU, using a string or numeric variable. If omitted, selected within-period
positions are used; see {help dea_malmquist##remarks:Remarks} for the identity assumption.

{phang}
{cmdab:r:ts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
{cmd:crs} assumes constant returns to scale; {cmd:vrs} assumes variable returns
to scale. Default is {cmd:rts(crs)}.

{phang}
{cmdab:or:ientation(}{it:type}{cmd:)} specifies the direction of optimization.
Default is {cmd:orientation(input)}.

{phang}
{cmd:tolerance(}{it:#}{cmd:)} must be positive and finite. From version 1.3.7,
it is passed to the native Stata {cmd:LinearProgram()} solver and to the
constraint/bound checks. Columns are rescaled internally, including the
evaluated observation when balancing cross-period problems. Native convergence
uses min(tolerance,1e-8); a converged candidate rejected by the primal checks
receives at most one stricter solve of the same LP. Unavailable or rejected
distances remain missing. Tolerance is a feasibility threshold, not a bound on
index error or a certificate of exact optimality. Earlier versions ignored
this option during optimization and accepted invalid values.


{marker remarks}{...}
{title:Remarks}

{pstd}
The effective analysis sample must contain at most one observation for each
DMU-period key. Repeated keys return {cmd:r(459)} before result variables are
created or replaced, even when the repeated observations have identical
inputs/outputs and even with {cmd:replacevars}. Resolve repeated keys using
a substantively appropriate rule before estimation; the command does not
choose one row or aggregate observations automatically.

{pstd}
The check applies after {cmd:if}/{cmd:in} restrictions and exclusion of missing
inputs, outputs, time values, or DMU identifiers. Repeated keys outside this
sample do not cause rejection. Nonmissing string identifiers are compared
exactly, including whitespace; numeric keys retain their stored precision.
Unbalanced panels and gaps between numeric time values remain supported.
Consecutive periods refer to the ordered observed periods, not necessarily
consecutive calendar dates. {cmd:time()} must be numeric.

{pstd}
Supply a stable {cmd:dmu()} identifier for panel analysis. If omitted, the
command numbers selected observations 1, 2, ... within each period, preserving
the incoming row order. Selection applies {cmd:if}/{cmd:in} and excludes missing
input, output, and time values before numbering. Excluded rows do not consume
an ID. The original row order is restored before computation.

{pstd}
Positional IDs assume that the same selected position represents the same
real DMU across periods. They do not infer identity: changing within-period
order or removing different DMUs across periods can change the matches. Use
an explicit stable {cmd:dmu()} for such panels. Artificial sequential IDs
cannot detect repeated observations of a real DMU. A note is displayed when
{cmd:dmu()} is omitted, and {cmd:r(dmu_method)} records the method used.
Recompute earlier results obtained without {cmd:dmu()}; earlier versions
counted excluded rows and did not fix within-period order.


{pstd}
The data must be in long format with one observation per DMU per period.
Only DMUs that appear in both consecutive periods are included in the
Malmquist calculations for that period pair.

{pstd}
Each period's reference frontier uses all selected complete observations in
that period, including units not matched to the other period. Nonnegative
quantities are required, with a positive input and a positive output in every
selected DMU-period. Zero components are permitted; no positive substitution
is made. Invalid quantities are rejected before replacing result columns.

{pstd}
VRS cross-period LPs can be genuinely infeasible even with strictly positive
data. Each computable component is retained: for example, EC can remain
available when cross-period distances prevent TC and MPI. Missing results
are counted separately from unmatched DMUs and first-period rows, which are
not attempted transitions. The counts alone do not distinguish mathematical
infeasibility from numerical failure. A period pair with no common DMUs has
zero counts and missing means, including when it is the final pair.

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

{pstd}
The global index is the ratio of later to earlier efficiency against the same
all-period pooled frontier. EC uses the two own-period frontiers, and the stored
TC is MPI/EC: the change in the gap between own-period and pooled best practice.
It differs from adjacent-period geometric technical change. Circularity holds
for the same fixed pooled reference data; changing the selected periods or
observations can change that frontier. Pooling contains each selected
observation and avoids cross-period infeasibility in these self-evaluations,
but does not guarantee numerical convergence.


{pstd}
{bf:Protecting existing data.} All output names are checked before results are
created. An existing result variable causes an error unless {cmd:replacevars}
is specified. The option applies only to variables produced by the current
command; it does not remove other variables sharing the prefix. Replacement
recreates whole result columns, including their storage types and labels;
observations not scored by the command do not retain old scores. Variables
used as model data or identifiers cannot be replaced, even with
{cmd:replacevars}. Invalid or duplicate output names are rejected before any
existing result column is changed. This name check does not promise rollback
of later computation or file-write errors.

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

{phang2}{cmd:. dea_malmquist expenditure teachers = graduation, dmu(country) time(year) rts(vrs) replacevars}{p_end}

{pstd}Global pooled-frontier Malmquist index

{phang2}{cmd:. dea_malmquist expenditure teachers = graduation, dmu(country) time(year) global replacevars}{p_end}

{pstd}Save decomposition results

{phang2}{cmd:. dea_malmquist expenditure teachers = graduation, dmu(country) time(year) generate(mpi)}{p_end}
{phang2}{cmd:. list country year mpi_ec mpi_tc mpi_mpi}{p_end}



{pstd}Equivalent sample restrictions in both syntax forms

{phang2}{cmd:. dea_malmquist expenditure teachers = graduation if country!="Aland" in 1/l, dmu(country) time(year) generate(qa)}{p_end}
{phang2}{cmd:. dea_malmquist if country!="Aland" in 1/l, inputs(expenditure teachers) outputs(graduation) dmu(country) time(year) generate(qb)}{p_end}
{phang2}{cmd:. assert qa_mpi==qb_mpi}{p_end}


{pstd}Documented minimum option abbreviations

{phang2}{cmd:. dea_malmquist, inputs(expenditure teachers) outputs(graduation) dmu(country) time(year) rts(crs) orientation(output) generate(fullopt)}{p_end}
{phang2}{cmd:. dea_malmquist, in(expenditure teachers) out(graduation) dmu(country) time(year) r(crs) or(output) gen(shortopt)}{p_end}
{phang2}{cmd:. assert fullopt_mpi==shortopt_mpi}{p_end}

{pstd}Detect repeated panel keys before estimation

{phang2}{cmd:. tempfile unique_panel}{p_end}
{phang2}{cmd:. save `unique_panel'}{p_end}
{phang2}{cmd:. expand 2 if _n==1}{p_end}
{phang2}{cmd:. capture noisily dea_malmquist expenditure teachers = graduation, dmu(country) time(year) gen(keycheck)}{p_end}
{phang2}{cmd:. assert _rc==459}{p_end}
{phang2}{cmd:. use `unique_panel', clear}{p_end}

{pstd}Selected positional IDs with consistently ordered DMUs

{phang2}{cmd:. tempfile positional_panel}{p_end}
{phang2}{cmd:. save `positional_panel'}{p_end}
{phang2}{cmd:. expand 2 if _n==1, generate(extra)}{p_end}
{phang2}{cmd:. sort year country extra}{p_end}
{phang2}{cmd:. dea_malmquist expenditure teachers = graduation if extra==0, time(year) dmu(country) generate(idref)}{p_end}
{phang2}{cmd:. dea_malmquist expenditure teachers = graduation if extra==0, time(year) generate(idpos)}{p_end}
{phang2}{cmd:. assert "`r(dmu_method)'"=="position"}{p_end}
{phang2}{cmd:. assert idref_mpi==idpos_mpi}{p_end}
{phang2}{cmd:. use `positional_panel', clear}{p_end}

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
{synopt:{cmd:r(mean_ec)}}arithmetic mean of available EC values in the last pair{p_end}
{synopt:{cmd:r(mean_tc)}}arithmetic mean of available TC values in the last pair{p_end}
{synopt:{cmd:r(mean_mpi)}}arithmetic mean of available MPI values in the last pair{p_end}
{synopt:{cmd:r(n_pairs)}}number of ordered adjacent period pairs{p_end}
{synopt:{cmd:r(n_failed)}}matched DMUs with any unavailable component in the last pair{p_end}
{synopt:{cmd:r(n_failed_ec)}}matched DMUs with missing EC in the last pair{p_end}
{synopt:{cmd:r(n_failed_tc)}}matched DMUs with missing TC in the last pair{p_end}
{synopt:{cmd:r(n_failed_mpi)}}matched DMUs with missing MPI in the last pair{p_end}
{synopt:{cmd:r(n_failed_total)}}sum of any-component failures over all pairs; counts transitions{p_end}
{synopt:{cmd:r(tolerance)}}requested feasibility tolerance{p_end}
{syntab:Matrices}
{synopt:{cmd:r(pair_summary)}}one row per pair: period0, period1, N, n_failed, n_failed_ec, n_failed_tc, n_failed_mpi, mean_ec, mean_tc, mean_mpi{p_end}
{syntab:Macros}
{synopt:{cmd:r(dmu_method)}}{cmd:variable} when {cmd:dmu()} is supplied; otherwise {cmd:position}{p_end}
{synopt:{cmd:r(malmquist_type)}}{cmd:adjacent} or {cmd:global}{p_end}
{synopt:{cmd:r(rts)}}returns to scale used{p_end}
{synopt:{cmd:r(orientation)}}orientation used{p_end}

{pstd}
Component means use their own available observations and are arithmetic means;
mean(MPI) need not equal mean(EC)*mean(TC). From version 1.3.7, final-pair returns
are reset even when the pair has no matched DMUs. Earlier versions could retain
the preceding pair's summary. Per-DMU generated variables are stored on the
later observation of each matched transition.
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
