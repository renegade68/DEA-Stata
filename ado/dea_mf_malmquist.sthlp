{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.0.0  20jul2026}{...}
{vieweralsosee "dea_metafrontier" "help dea_metafrontier"}{...}
{vieweralsosee "dea_malmquist" "help dea_malmquist"}{...}
{title:Title}

{p2colset 5 28 30 2}{...}
{p2col:{cmd:dea_mf_malmquist} {hline 2}}Global metafrontier Malmquist index{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:dea_mf_malmquist} {it:inputvars} {cmd:=} {it:outputvars}
{ifin}{cmd:,} {cmd:group(}{it:varname}{cmd:)}
{cmd:dmu(}{it:varname}{cmd:)} {cmd:time(}{it:varname}{cmd:)} [{it:options}]

{p 8 16 2}{cmd:dea_mf_malmquist} {it:outputvars} {ifin}{cmd:,}
{cmd:inputs(}{it:varlist}{cmd:)} {cmd:group(}{it:varname}{cmd:)}
{cmd:dmu(}{it:varname}{cmd:)} {cmd:time(}{it:varname}{cmd:)} [{it:options}]

{synoptset 31 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:group(}{it:varname}{cmd:)}}technology-group identifier; required{p_end}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier; required{p_end}
{synopt:{cmd:time(}{it:varname}{cmd:)}}numeric period identifier; required{p_end}
{synopt:{cmd:rts(crs|vrs)}}returns to scale; default {cmd:crs}{p_end}
{synopt:{cmd:orientation(input|output)}}orientation; default {cmd:input}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}
{synopt:{cmd:engine(auto|simplex|ipm)}}LP engine passed to {cmd:dea}{p_end}
{synopt:{cmd:generate(}{it:prefix}{cmd:)}}result prefix; default {cmd:_mfm}{p_end}
{synopt:{cmd:tolerance(}{it:#}{cmd:)}}positive numerical tolerance; default 1e-7{p_end}
{synopt:{cmd:saving(}{it:filename}{cmd:)}}save the estimation sample and results{p_end}
{synopt:{cmd:replace}}replace an existing {cmd:saving()} file{p_end}
{synopt:{cmd:detail}}list transition-level results{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:dea_mf_malmquist} constructs a global technology for each group by
pooling all selected periods within that group. It also constructs one global
metafrontier by pooling all groups and periods. For each DMU and adjacent
period pair it generates group productivity change (GMPI), technology-gap
change (TGCH), and metafrontier productivity change (MMPI), satisfying

{p 12 16 2}{it:MMPI = GMPI x TGCH}.{p_end}

{pstd}
The common global benchmarks give a circular productivity index. Values above
one indicate improvement and values below one indicate decline. The command
requires a balanced panel, unique {cmd:dmu()}--{cmd:time()} observations, and
constant group membership within each DMU. The first observation of every DMU
has missing change indexes because it has no preceding period.

{pstd}
Generated variables are {it:prefix}{cmd:_group_eff},
{it:prefix}{cmd:_meta_eff}, {it:prefix}{cmd:_tgr},
{it:prefix}{cmd:_gmpi}, {it:prefix}{cmd:_tgch},
{it:prefix}{cmd:_mmpi}, and {it:prefix}{cmd:_identity_gap}.

{title:Example}

{pstd}Copy the command lines below (omit leading dots) from {cmd:preserve}
through {cmd:restore} into the Do-file Editor and run them as one block.
Examples use synthetic
data or Stata's bundled {cmd:auto} data and require no downloaded datasets.
The block restores the data in memory. Run dependent lines in order.{p_end}
{* DEAS_HELP_EXAMPLES_BEGIN}{...}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input firm region year labor land output}{p_end}
{phang2}{cmd:. 1 1 1 1 2 1}{p_end}
{phang2}{cmd:. 1 1 2 1 2 1.1}{p_end}
{phang2}{cmd:. 2 1 1 2 4 1.5}{p_end}
{phang2}{cmd:. 2 1 2 2 4 1.7}{p_end}
{phang2}{cmd:. 3 2 1 3 6 1}{p_end}
{phang2}{cmd:. 3 2 2 3 6 1.1}{p_end}
{phang2}{cmd:. 4 2 1 4 8 1.2}{p_end}
{phang2}{cmd:. 4 2 2 4 8 1.4}{p_end}
{phang2}{cmd:. end}{p_end}


{phang2}{cmd:. dea_mf_malmquist labor land = output, group(region) dmu(firm) time(year) rts(crs) generate(mfm)}{p_end}
{phang2}{cmd:. list firm year mfm_gmpi mfm_tgch mfm_mmpi}{p_end}


{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{title:Stored results}

{p2colset 9 30 32 2}{...}
{p2col:{cmd:r(N)}}selected observations{p_end}
{p2col:{cmd:r(N_groups)}}technology groups{p_end}
{p2col:{cmd:r(N_periods)}}periods{p_end}
{p2col:{cmd:r(N_transitions)}}DMU-period transitions{p_end}
{p2col:{cmd:r(mean_group_mpi)}}mean group productivity change{p_end}
{p2col:{cmd:r(mean_tgch)}}mean technology-gap change{p_end}
{p2col:{cmd:r(mean_meta_mpi)}}mean metafrontier productivity change{p_end}
{p2col:{cmd:r(max_identity_gap)}}largest absolute decomposition residual{p_end}
{p2col:{cmd:r(index)}}{cmd:global_metafrontier}{p_end}
{p2col:{cmd:r(rts)}}returns-to-scale specification{p_end}
{p2col:{cmd:r(orientation)}}orientation{p_end}
{p2colreset}{...}

{title:References}

{phang}O'Donnell, C. J., D. S. P. Rao, and G. E. Battese. 2008.
Metafrontier frameworks for the study of firm-level efficiencies and
technology ratios. {it:Empirical Economics} 34: 231--255.{p_end}

{phang}Oh, D.-h. 2010. A global Malmquist--Luenberger productivity index.
{it:Journal of Productivity Analysis} 34: 183--197.{p_end}

{title:Author}

{pstd}Choonjoo Lee, Korea National Defense University.{p_end}
