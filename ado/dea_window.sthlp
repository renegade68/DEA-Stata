{smcl}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* *! version 1.3.0  02may2026}{...}
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
{cmdab:dmu(}{it:varname}{cmd:)}
{cmdab:time(}{it:varname}{cmd:)}
[{it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_window}
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
{synopt:{cmd:width(}{it:#}{cmd:)}}window width in periods; default is {cmd:width(3)}{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for saved window-efficiency variables; default is {cmd:_win}{p_end}
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
DMU j in period t is evaluated against all observations in the window [t, t+w-1].
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
{cmd:dmu(}{it:varname}{cmd:)} specifies the variable identifying each DMU.
Accepts both string and numeric variables.

{phang}
{cmd:time(}{it:varname}{cmd:)} specifies the numeric variable for time periods.

{dlgtab:Model}

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
number of periods). Default is {cmd:width(3)}.

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
The data must be in long format with one observation per DMU per period. DMUs
need not be balanced across all periods; missing observations within a window
are handled automatically.

{pstd}
Window analysis was introduced by {help dea_window##CCGSS1985:Charnes et al. (1985)}
as a way to analyze panel data without requiring parametric assumptions about
the production technology over time. Standard DEA background is provided in
{help dea_window##CST2007:Cooper, Seiford, and Tone (2007)}.


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

{phang2}{cmd:. dea_window expenditure teachers = graduation, dmu(country) time(year) width(2) rts(vrs)}{p_end}

{pstd}Output-oriented window analysis

{phang2}{cmd:. dea_window expenditure teachers = graduation, dmu(country) time(year) width(2) orientation(output)}{p_end}

{pstd}Save window scores with a custom prefix

{phang2}{cmd:. dea_window expenditure teachers = graduation, dmu(country) time(year) width(2) generate(win)}{p_end}
{phang2}{cmd:. list country year win_w1_p1 win_w1_p2 win_w2_p1 win_w2_p2 win_mean win_sd}{p_end}



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
{syntab:Macros}
{synopt:{cmd:r(rts)}}returns to scale used{p_end}
{synopt:{cmd:r(orientation)}}orientation used{p_end}
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
