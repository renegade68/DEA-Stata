{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.0.0 20jul2026}{...}
{title:Title}

{phang}{bf:dea_dynamic} {hline 2} Dynamic slacks-based DEA with carry-over links

{title:Syntax}

{p 8 16 2}
{cmd:dea_dynamic} {ifin}, {cmd:inputs(}{it:varlist}{cmd:)}
{cmd:outputs(}{it:varlist}{cmd:)} {cmd:dmu(}{it:varname}{cmd:)}
{cmd:time(}{it:varname}{cmd:)} {cmd:carryover(}{it:varlist}{cmd:)}
[{it:options}]

{title:Description}

{pstd}
{cmd:dea_dynamic} implements the input-oriented, output-oriented, and
non-oriented dynamic SBM of Tone and Tsutsui (2010). One intensity vector is
estimated for each period, and equality constraints connect projected
carry-over flows in adjacent periods. The data must form a strictly balanced
DMU-by-time panel and all model quantities must be strictly positive.

{title:Options}

{phang}{cmd:carrytype(}{it:typelist}{cmd:)} supplies one type for each variable
in {cmd:carryover()}. Types are {cmd:desirable}, {cmd:undesirable},
{cmd:discretionary}, and {cmd:nondiscretionary}. Synonyms {cmd:good},
{cmd:bad}, {cmd:free}, and {cmd:fixed} are accepted. The default is desirable
for every carry-over.

{phang}{cmd:orientation(input|output|nonoriented)} selects the dynamic SBM
orientation. {cmd:ort()} is an alias. The default is {cmd:nonoriented}.

{phang}{cmd:rts(crs|vrs)} selects returns to scale. The default is {cmd:vrs}.

{phang}{cmd:generate(}{it:prefix}{cmd:)} creates
{it:prefix}{cmd:_overall} and {it:prefix}{cmd:_period}. The default prefix is
{cmd:_dyn}.

{phang}{cmd:tolerance(}{it:#}{cmd:)} sets the LP tolerance; default
{cmd:1e-7}.

{phang}{cmd:detail} lists DMU-period results. {cmd:debug} reports failed LP
status codes.

{title:Examples}

{pstd}Copy the command lines below (omit leading dots) from {cmd:preserve}
through {cmd:restore} into the Do-file Editor and run them as one block.
Examples use synthetic
data or Stata's bundled {cmd:auto} data and require no downloaded datasets.
The block restores the data in memory. Run dependent lines in order.{p_end}
{* DEAS_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input id year x y z}{p_end}
{phang2}{cmd:. 1 1 1 1 1}{p_end}
{phang2}{cmd:. 1 2 1 1 1}{p_end}
{phang2}{cmd:. 2 1 2 1 1}{p_end}
{phang2}{cmd:. 2 2 2 1 1}{p_end}
{phang2}{cmd:. end}{p_end}
{phang2}{cmd:. dea_dynamic, inputs(x) outputs(y) dmu(id) time(year) carryover(z) orientation(input) generate(dyn)}{p_end}
{phang2}{cmd:. assert r(n_failed) == 0}{p_end}
{phang2}{cmd:. list id year dyn_overall dyn_period}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{title:Stored results}

{synoptset 20 tabbed}{...}
{synopt:{cmd:r(n_dmus)}}number of DMUs{p_end}
{synopt:{cmd:r(T)}}number of periods{p_end}
{synopt:{cmd:r(mean_eff)}}mean overall efficiency{p_end}
{synopt:{cmd:r(n_failed)}}number of failed DMU problems{p_end}
{synopt:{cmd:r(rts)}}returns-to-scale specification{p_end}
{synopt:{cmd:r(orientation)}}orientation{p_end}
{synopt:{cmd:r(carrytype)}}normalized carry-over type list{p_end}

{title:Reference}

{pstd}
Tone, K., and M. Tsutsui. 2010. Dynamic DEA: A slacks-based measure
approach. {it:Omega} 38:145-156. doi:10.1016/j.omega.2009.07.003.
