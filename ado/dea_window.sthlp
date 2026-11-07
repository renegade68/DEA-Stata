{smcl}
{* Help revision 2026-11-02: validated window LPs and explicit failure statistics.}{...}
{* Help revision 2026-11-01: validate selected DMU-period keys before output replacement.}{...}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* *! version 1.3.8  02nov2026}{...}
{viewerjumpto "Related article" "dea_window##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_bootstrap" "help dea_bootstrap"}{...}
{vieweralsosee "dea_malmquist" "help dea_malmquist"}{...}
{viewerjumpto "Syntax" "dea_window##syntax"}{...}
{viewerjumpto "Description" "dea_window##description"}{...}
{viewerjumpto "Options" "dea_window##options"}{...}
{viewerjumpto "Remarks" "dea_window##remarks"}{...}
{viewerjumpto "Examples" "dea_window##examples"}{...}
{viewerjumpto "Stored results" "dea_window##stored_results"}{...}
{viewerjumpto "References" "dea_window##references"}{...}
{viewerjumpto "Authors" "dea_window##authors"}{...}
{title:Title}

{p2colset 5 20 22 2}{...}
{p2col :{cmd:dea_window} {hline 2}}Window analysis DEA{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea_window}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}{cmd:,}
[{cmdab:dmu(}{it:varname}{cmd:)}]
{cmdab:time(}{it:varname}{cmd:)}
[{it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_window}
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
{synopt:{cmd:width(}{it:#}{cmd:)}}window width in observed periods; default is min(3,T){p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for saved window-efficiency variables; default is {cmd:_win}{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}
{synopt:{cmd:detail}}accepted for compatibility; window-by-window results are displayed by default{p_end}

{syntab:Technical}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmd:debug}}display debugging information{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:rtstype}: {cmd:crs} | {cmd:vrs} | {cmd:drs} | {cmd:irs}{p_end}
{p 4 6 2}
{it:type}: {cmd:input} | {cmd:output}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_window} implements window analysis for DEA panel data. Window analysis
treats each DMU in each period as a distinct entity within a rolling window of
{it:w} consecutive periods, pooling the data across those periods to form the
reference technology ({help dea_window##CCGSS1985:Charnes et al. 1985}).

{pstd}
For a panel with T periods and window width w, the number of windows is T - w + 1.
Each DMU-period is evaluated against all selected observations in every window
containing that period.
Because each DMU appears multiple times across windows, efficiency may vary,
providing information about the stability of efficiency over time.

{pstd}
For each DMU, the command reports:

{phang2}
Efficiency scores for each window in which the DMU appears.

{phang2}
The mean and standard deviation of efficiency scores across windows, indicating
the stability of relative efficiency over time.


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{cmd:time(}{it:varname}{cmd:)} specifies the numeric variable for time periods.

{dlgtab:Model}

{phang}
{cmd:dmu(}{it:varname}{cmd:)} optionally specifies a stable identifier for each
DMU, using a string or numeric variable. If omitted, selected within-period
positions are used; see {help dea_window##remarks:Remarks} for the identity assumption.

{phang}
{cmdab:r:ts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
Allowed values are {cmd:crs} (constant), {cmd:vrs} (variable), {cmd:drs}
(decreasing), and {cmd:irs} (increasing). Default is {cmd:rts(crs)}.

{phang}
{cmdab:or:ientation(}{it:type}{cmd:)} specifies the direction of optimization.
Default is {cmd:orientation(input)}.

{phang}
{cmd:width(}{it:#}{cmd:)} specifies the window width — the number of consecutive
periods included in each window. Must be at least 2 and at most T (the total
number of periods). Omitting the option or specifying 0 selects min(3,T).

{phang}
{cmd:tolerance(}{it:#}{cmd:)} must be positive and finite. From version 1.3.8,
the native Stata {cmd:LinearProgram()} path applies the requested tolerance to
normalized constraint/bound checks and uses min(tolerance,1e-8) for convergence.
A converged candidate rejected by the primal checks receives at most one
stricter solve of the same LP. Columns are rescaled without changing the model.
Unavailable or rejected scores remain missing. Tolerance does not bound the
error in an efficiency score or certify exact optimality. Earlier versions
ignored this option during optimization and accepted invalid values.

{dlgtab:Output}

{phang}
{cmdab:gen:erate(}{it:prefix}{cmd:)} specifies the prefix for saved
window-efficiency variables. Scores are saved as
{it:prefix}{cmd:_w}{it:W}{cmd:_p}{it:P}, where {it:W} is the window number and
{it:P} is the observation's period position within that window. The command also
saves {it:prefix}{cmd:_mean} and {it:prefix}{cmd:_sd}. Default is {cmd:_win}.


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
Quantities must be nonnegative, with a positive input and a positive output
in every selected DMU-period. Individual zero components are allowed.
Validation precedes output replacement. Every selected observation is included
in each applicable window's reference set; unavailable scores can therefore
reflect numerical failure even though self-evaluation provides a feasible point.

{pstd}
An observation's generated mean and SD summarize its available scores across
windows. A single available score has row SD 0; no available score leaves both
summaries missing. The displayed DMU summary pools all available period-window
scores for that DMU. Overall returns pool all available scores, so observations
contained in more overlapping windows contribute more often. These are not
equal-weight averages of the per-observation means. SD uses the sample divisor;
the overall SD is missing with fewer than two available scores.

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
The data must be in long format with one observation per DMU per period. DMUs
need not be balanced across all periods; missing observations within a window
are handled automatically.

{pstd}
Window analysis was introduced by {help dea_window##CCGSS1985:Charnes et al. (1985)}
as a way to analyze panel data without requiring parametric assumptions about
the production technology over time. Standard DEA background is provided in
{help dea_window##CST2007:Cooper, Seiford, and Tone (2007)}.


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
{phang2}{cmd:"Aland"   2023 136 33 108}{p_end}
{phang2}{cmd:"Bordia"  2023 166 41 126}{p_end}
{phang2}{cmd:"Cyrenia" 2023 148 37 114}{p_end}
{phang2}{cmd:"Dameria" 2023 188 48 140}{p_end}
{phang2}{cmd:"Eldora"  2023 120 31  98}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Window analysis with default width 3

{phang2}{cmd:. dea_window expenditure teachers = graduation, dmu(country) time(year)}{p_end}

{pstd}Window width 2 with VRS

{phang2}{cmd:. dea_window expenditure teachers = graduation, dmu(country) time(year) width(2) rts(vrs) replacevars}{p_end}

{pstd}Output-oriented window analysis

{phang2}{cmd:. dea_window expenditure teachers = graduation, dmu(country) time(year) width(2) orientation(output) replacevars}{p_end}

{pstd}Save window scores with a custom prefix

{phang2}{cmd:. dea_window expenditure teachers = graduation, dmu(country) time(year) width(2) generate(win)}{p_end}
{phang2}{cmd:. list country year win_w1_p1 win_w1_p2 win_w2_p1 win_w2_p2 win_mean win_sd}{p_end}



{pstd}Equivalent sample restrictions in both syntax forms

{phang2}{cmd:. dea_window expenditure teachers = graduation if country!="Aland" in 1/l, dmu(country) time(year) width(2) generate(qa)}{p_end}
{phang2}{cmd:. dea_window if country!="Aland" in 1/l, inputs(expenditure teachers) outputs(graduation) dmu(country) time(year) width(2) generate(qb)}{p_end}
{phang2}{cmd:. assert qa_mean==qb_mean}{p_end}


{pstd}Documented minimum option abbreviations

{phang2}{cmd:. dea_window, inputs(expenditure teachers) outputs(graduation) dmu(country) time(year) width(2) rts(vrs) orientation(output) generate(fullopt)}{p_end}
{phang2}{cmd:. dea_window, in(expenditure teachers) out(graduation) dmu(country) time(year) width(2) r(vrs) or(output) gen(shortopt)}{p_end}
{phang2}{cmd:. assert fullopt_mean==shortopt_mean}{p_end}

{pstd}Detect repeated panel keys before estimation

{phang2}{cmd:. tempfile unique_panel}{p_end}
{phang2}{cmd:. save `unique_panel'}{p_end}
{phang2}{cmd:. expand 2 if _n==1}{p_end}
{phang2}{cmd:. capture noisily dea_window expenditure teachers = graduation, dmu(country) time(year) gen(keycheck)}{p_end}
{phang2}{cmd:. assert _rc==459}{p_end}
{phang2}{cmd:. use `unique_panel', clear}{p_end}

{pstd}Selected positional IDs with consistently ordered DMUs

{phang2}{cmd:. tempfile positional_panel}{p_end}
{phang2}{cmd:. save `positional_panel'}{p_end}
{phang2}{cmd:. expand 2 if _n==1, generate(extra)}{p_end}
{phang2}{cmd:. sort year country extra}{p_end}
{phang2}{cmd:. dea_window expenditure teachers = graduation if extra==0, time(year) dmu(country) generate(idref)}{p_end}
{phang2}{cmd:. dea_window expenditure teachers = graduation if extra==0, time(year) generate(idpos)}{p_end}
{phang2}{cmd:. assert "`r(dmu_method)'"=="position"}{p_end}
{phang2}{cmd:. assert idref_mean==idpos_mean}{p_end}
{phang2}{cmd:. use `positional_panel', clear}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea_window} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:r(n_dmus)}}number of DMUs{p_end}
{synopt:{cmd:r(n_periods)}}number of time periods{p_end}
{synopt:{cmd:r(width)}}window width used{p_end}
{synopt:{cmd:r(n_windows)}}number of windows{p_end}
{synopt:{cmd:r(mean_eff)}}mean window efficiency across reported scores{p_end}
{synopt:{cmd:r(sd_eff)}}standard deviation of window efficiency scores{p_end}
{synopt:{cmd:r(n_obs)}}number of window efficiency scores summarized{p_end}
{synopt:{cmd:r(n_attempted)}}number of observed DMU-period-window evaluations{p_end}
{synopt:{cmd:r(n_failed)}}attempted evaluations with unavailable scores{p_end}
{synopt:{cmd:r(n_failed_obs)}}selected observations with no available window score{p_end}
{synopt:{cmd:r(tolerance)}}requested feasibility tolerance{p_end}
{syntab:Macros}
{synopt:{cmd:r(dmu_method)}}{cmd:variable} when {cmd:dmu()} is supplied; otherwise {cmd:position}{p_end}
{synopt:{cmd:r(rts)}}returns to scale used{p_end}
{synopt:{cmd:r(orientation)}}orientation used{p_end}

{pstd}
{cmd:r(n_attempted)} = {cmd:r(n_obs)} + {cmd:r(n_failed)}. Structural blank
cells outside an observation's windows or for absent DMU-period combinations
are not failures. All-missing runs return zero successful scores, missing
overall mean/SD, and explicit failure counts.
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker CCGSS1985}{...}
{phang}
Charnes, A., W. W. Cooper, B. Golany, L. Seiford, and J. Stutz. 1985.
Foundations of data envelopment analysis for Pareto-Koopmans efficient empirical
production functions.
{it:Journal of Econometrics} 30(1--2): 91--107.
{browse "https://doi.org/10.1016/0304-4076(85)90133-2":doi:10.1016/0304-4076(85)90133-2}

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
