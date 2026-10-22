{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.4.0  20jul2026}{...}
{vieweralsosee "dea" "help dea"}{...}
{title:Title}

{p2colset 5 26 28 2}{...}
{p2col:{cmd:dea_metafrontier} {hline 2}}Nonparametric DEA-based metafrontier{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 16 2}{cmd:dea_metafrontier} {it:inputvars} {cmd:=} {it:outputvars}
{ifin}{cmd:,} {cmd:group(}{it:varname}{cmd:)} [{it:options}]

{p 8 16 2}{cmd:dea_metafrontier} {it:outputvars} {ifin}{cmd:,}
{cmd:inputs(}{it:varlist}{cmd:)} {cmd:group(}{it:varname}{cmd:)} [{it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{cmd:group(}{it:varname}{cmd:)}}group identifier; required{p_end}
{synopt:{cmd:method(dea)}}nonparametric DEA method; currently the only method{p_end}
{synopt:{cmd:rts(crs|vrs|drs|irs)}}returns to scale; default {cmd:vrs}{p_end}
{synopt:{cmd:orientation(input|output)}}orientation; default {cmd:input}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}
{synopt:{cmd:generate(}{it:prefix}{cmd:)}}output prefix; default {cmd:_mf}{p_end}
{synopt:{cmd:engine(auto|simplex|ipm)}}LP engine passed to {cmd:dea}{p_end}
{synopt:{cmd:tolerance(}{it:#}{cmd:)}}positive numerical tolerance; default 1e-7{p_end}
{synopt:{cmd:minobs(}{it:#}{cmd:)}}require at least # observations in every group{p_end}
{synopt:{cmd:pooltest(}{it:#}{cmd:)}}run a group-label permutation test with # replications{p_end}
{synopt:{cmd:breps(}{it:#}{cmd:)}}run at least 50 joint smooth-bootstrap replications; default 0{p_end}
{synopt:{cmd:level(}{it:#}{cmd:)}}confidence level for bootstrap TGR intervals; default 95{p_end}
{synopt:{cmd:seed(}{it:#}{cmd:)}}random-number seed for resampling options{p_end}
{synopt:{cmd:saving(}{it:filename}{cmd:)}}save the estimation sample and results{p_end}
{synopt:{cmd:replace}}replace an existing {cmd:saving()} file{p_end}
{synopt:{cmd:detail}}display group means{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:dea_metafrontier} estimates the input- or output-oriented DEA frontier
separately within each group and then for the pooled sample.  It generates
{it:prefix}{cmd:_group}, {it:prefix}{cmd:_meta}, {it:prefix}{cmd:_tgr}, and
{it:prefix}{cmd:_gap}.  The technology gap ratio is
theta_meta/theta_group and the gap is 1-TGR.  No functional form or random
noise component is specified.

{pstd}
The command reports the smallest group size and the descriptive
3*(number of inputs + number of outputs) rule.  The rule is a diagnostic, not
a validity theorem.  {cmd:minobs()} converts a researcher-selected minimum
group size into an enforceable data contract.

{pstd}
{cmd:pooltest()} tests the exchangeability implication of a common technology.
It uses mean(1-TGR) as the test statistic and randomly reallocates the observed
group labels while preserving every group size.  The reported Monte Carlo
p-value is (1 + number of permuted statistics at least as large as observed) /
(1 + replications).  This is a conditional randomization diagnostic; dependence,
clustered sampling, or nonexchangeable covariates require a research-specific
resampling design.  The command restores the caller's random-number state.

{pstd}
{cmd:breps()} applies the Simar--Wilson smooth-bootstrap construction within
each group. In every replication, group pseudo frontiers are generated jointly
and their union forms the corresponding pooled pseudo metafrontier. The same
replication therefore supplies both efficiency components of each TGR. The
command generates {it:prefix}{cmd:_tgr_bc}, {it:prefix}{cmd:_tgr_se},
{it:prefix}{cmd:_tgr_lo}, and {it:prefix}{cmd:_tgr_hi}. The interval is the
pivot interval truncated to the theoretical [0,1] range. At least 50
replications are required for testing; 1,000 or more are recommended for
reported empirical work. Dependence across observations requires a
research-specific block or cluster bootstrap and is not handled automatically.

{pstd}
{cmd:saving()} writes a compact Stata dataset containing the requested DMU
identifier (when supplied), group, input and output variables, and the four
generated result variables. The caller's dataset, observation order, and
variable set are restored after saving. {cmd:replace} is allowed only with
{cmd:saving()}.

{pstd}
This command implements the DEA-based metafrontier, not the deterministic
coefficient-envelopment SFA model of Battese, Rao, and O'Donnell (2004) and
not the stochastic metafrontier of Huang, Huang, and Liu (2014).

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
{phang2}{cmd:. input region labor land production}{p_end}
{phang2}{cmd:. 1 1 2 1}{p_end}
{phang2}{cmd:. 1 2 3 1.5}{p_end}
{phang2}{cmd:. 1 3 5 2}{p_end}
{phang2}{cmd:. 1 4 6 2.3}{p_end}
{phang2}{cmd:. 2 2 3 1}{p_end}
{phang2}{cmd:. 2 3 4 1.3}{p_end}
{phang2}{cmd:. 2 4 6 1.6}{p_end}
{phang2}{cmd:. 2 5 8 2}{p_end}
{phang2}{cmd:. end}{p_end}


{phang2}{cmd:. dea_metafrontier labor land = production, group(region) rts(vrs) ort(input) generate(mf) detail}
{phang2}{cmd:. dea_metafrontier production, inputs(labor land) group(region) rts(crs) orientation(output) generate(mfout)}
{phang2}{cmd:. dea_metafrontier labor land = production, group(region) breps(100) level(95) seed(12345) generate(mfbs)}


{pstd}The small sample and 100 bootstrap repetitions illustrate syntax only; use adequate data and repetitions for inference.{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{title:Stored results}

{p2colset 9 28 30 2}{...}
{p2col:{cmd:r(N)}}number of evaluated DMUs{p_end}
{p2col:{cmd:r(N_groups)}}number of groups{p_end}
{p2col:{cmd:r(mean_group)}}mean group-frontier efficiency{p_end}
{p2col:{cmd:r(mean_meta)}}mean pooled-frontier efficiency{p_end}
{p2col:{cmd:r(mean_tgr)}}mean technology gap ratio{p_end}
{p2col:{cmd:r(min_group_N)}}smallest group sample size{p_end}
{p2col:{cmd:r(recommended_N)}}descriptive 3*(inputs+outputs) threshold{p_end}
{p2col:{cmd:r(pool_stat)}}observed mean technology-gap statistic{p_end}
{p2col:{cmd:r(p_pool)}}permutation p-value; missing when not requested{p_end}
{p2col:{cmd:r(pool_reps)}}number of requested permutations{p_end}
{p2col:{cmd:r(bootstrap_reps)}}number of smooth-bootstrap replications{p_end}
{p2col:{cmd:r(level)}}bootstrap confidence level{p_end}
{p2col:{cmd:r(mean_tgr_bc)}}mean bias-corrected TGR; missing when not requested{p_end}
{p2col:{cmd:r(mean_tgr_se)}}mean observation-level bootstrap standard error{p_end}
{p2col:{cmd:r(saving)}}requested results filename; empty when not requested{p_end}
{p2col:{cmd:r(group_summary)}}group means and group sample sizes{p_end}
{p2col:{cmd:r(method)}}{cmd:dea}{p_end}
{p2col:{cmd:r(rts)}}returns-to-scale specification{p_end}
{p2col:{cmd:r(orientation)}}orientation{p_end}
{p2colreset}{...}

{title:References}

{phang}O'Donnell, C. J., D. S. P. Rao, and G. E. Battese. 2008.
Metafrontier frameworks for the study of firm-level efficiencies and
technology ratios. {it:Empirical Economics} 34: 231-255.{p_end}

{title:Author}

{pstd}Choonjoo Lee, Korea National Defense University.{p_end}
