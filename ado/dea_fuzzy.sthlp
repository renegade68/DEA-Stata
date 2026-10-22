{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.0.1  03may2026}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_bootstrap" "help dea_bootstrap"}{...}
{vieweralsosee "dea_window" "help dea_window"}{...}
{viewerjumpto "Syntax" "dea_fuzzy##syntax"}{...}
{viewerjumpto "Description" "dea_fuzzy##description"}{...}
{viewerjumpto "Options" "dea_fuzzy##options"}{...}
{viewerjumpto "Remarks" "dea_fuzzy##remarks"}{...}
{viewerjumpto "Examples" "dea_fuzzy##examples"}{...}
{viewerjumpto "Stored results" "dea_fuzzy##stored_results"}{...}
{viewerjumpto "References" "dea_fuzzy##references"}{...}
{viewerjumpto "Authors" "dea_fuzzy##authors"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{cmd:dea_fuzzy} {hline 2}}Fuzzy data envelopment analysis (α-cut approach){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:dea_fuzzy}
{ifin}{cmd:,}
{cmdab:in:puts(}{it:varlist}{cmd:)}
{cmdab:out:puts(}{it:varlist}{cmd:)}
[{it:options}]


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{cmdab:in:puts(}{it:varlist}{cmd:)}}center (modal) values of input variables{p_end}
{synopt:{cmdab:out:puts(}{it:varlist}{cmd:)}}center (modal) values of output variables{p_end}

{syntab:Fuzzy bounds}
{synopt:{cmd:xlo(}{it:varlist}{cmd:)}}lower-bound variables for inputs at α = 0; default uses {cmd:inputs()}{p_end}
{synopt:{cmd:xhi(}{it:varlist}{cmd:)}}upper-bound variables for inputs at α = 0; default uses {cmd:inputs()}{p_end}
{synopt:{cmd:ylo(}{it:varlist}{cmd:)}}lower-bound variables for outputs at α = 0; default uses {cmd:outputs()}{p_end}
{synopt:{cmd:yhi(}{it:varlist}{cmd:)}}upper-bound variables for outputs at α = 0; default uses {cmd:outputs()}{p_end}

{syntab:Model}
{synopt:{cmdab:alp:ha(}{it:numlist}{cmd:)}}α-cut levels ∈ [0,1]; default is {cmd:alpha(0 0.25 0.5 0.75 1)}{p_end}
{synopt:{cmdab:r:ts(}{it:rtstype}{cmd:)}}returns to scale; default is {cmd:rts(crs)}{p_end}
{synopt:{cmdab:or:ientation(}{it:type}{cmd:)}}orientation; default is {cmd:orientation(input)}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for efficiency variables; default is {cmd:fdea_eff}{p_end}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable{p_end}
{synopt:{cmd:detail}}display DMU-level efficiency table{p_end}

{syntab:Technical}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmd:debug}}display debugging information{p_end}
{synoptline}

{p 4 6 2}
{it:rtstype}: {cmd:crs} {cmd:ccr} {cmd:vrs} {cmd:bcc} {cmd:drs} {cmd:irs}{p_end}
{p 4 6 2}
{it:type}: {cmd:input} {cmd:i} {cmd:output} {cmd:o}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_fuzzy} computes Fuzzy Data Envelopment Analysis (FDEA) efficiency
scores using the α-cut (possibility) approach
({help dea_fuzzy##KL2000:Kao and Liu 2000};
{help dea_fuzzy##LFJ2003:Lertworasirikul et al. 2003}).
When inputs and/or outputs are not precisely known but can be
characterized as triangular fuzzy numbers, this command evaluates the
efficiency of each DMU at each specified α-level, producing a pair of
bounds — pessimistic and optimistic — that together trace the
{it:fuzzy efficiency membership function}.

{pstd}
A triangular fuzzy number for variable {it:x} is defined by three values
({it:x_lo}, {it:x_center}, {it:x_hi}) where {it:x_center} is the modal
(most likely) value, {it:x_lo} is the lower support, and {it:x_hi} is the
upper support.  At α-cut level α ∈ [0,1], the fuzzy number is
represented by the crisp interval [{it:x_L}(α), {it:x_U}(α)] where

{p 8 8 2}
{it:x_L}(α) = {it:x_lo} + α × ({it:x_center} − {it:x_lo})

{p 8 8 2}
{it:x_U}(α) = {it:x_hi} − α × ({it:x_hi} − {it:x_center})

{pstd}
At α = 1 both bounds collapse to {it:x_center}, so the result is identical
to standard (crisp) DEA.  At α = 0 the full support [{it:x_lo}, {it:x_hi}]
is used, representing maximum uncertainty.

{pstd}
For each α-level, two LPs are solved per DMU:

{phang2}
{it:Pessimistic LP.}  DMU {it:k} is evaluated at its worst-case values
(high inputs, low outputs).  Peer DMUs enter the reference set at favorable
bounds (low inputs, high outputs), while DMU {it:k} enters the reference set
at the same worst-case realization.  This gives the {it:lower efficiency bound}
θ_L(α).

{phang2}
{it:Optimistic LP.}  DMU {it:k} is evaluated at its best-case values
(low inputs, high outputs).  Peer DMUs enter the reference set at unfavorable
bounds (high inputs, low outputs), while DMU {it:k} enters the reference set
at the same best-case realization.  This gives the {it:upper efficiency bound}
θ_U(α).

{pstd}
The pair [θ_L(α), θ_U(α)] is the α-cut of the fuzzy efficiency score.
Plotting θ_L and θ_U against α traces the membership function of the
fuzzy efficiency score.

{pstd}
For output-oriented models, the raw expansion factor φ (≥ 1) is
converted to 1/φ ∈ (0,1] for consistency with the input-oriented scale,
following the {cmd:dea} package convention for imprecise DEA (IDEA).


{marker options}{...}
{title:Options}

{dlgtab:Fuzzy bounds}

{phang}
{cmd:xlo(}{it:varlist}{cmd:)} specifies the lower-bound variables for inputs
at α = 0 (the full lower support of the fuzzy number).  The number of
variables must equal the number of variables in {cmd:inputs()}, and each
corresponds to the matching input variable in order.  If not specified,
the center values in {cmd:inputs()} are used, making those inputs crisp.

{phang}
{cmd:xhi(}{it:varlist}{cmd:)} specifies the upper-bound variables for inputs
at α = 0.  Same correspondence rules as {cmd:xlo()}.

{phang}
{cmd:ylo(}{it:varlist}{cmd:)} specifies the lower-bound variables for outputs
at α = 0.  If not specified, the center values in {cmd:outputs()} are used.

{phang}
{cmd:yhi(}{it:varlist}{cmd:)} specifies the upper-bound variables for outputs
at α = 0.  If not specified, the center values in {cmd:outputs()} are used.

{pstd}
Fuzzy bounds can be specified for any subset of inputs and outputs.
Variables without bounds are treated as crisp (the lower and upper
bounds both equal the center value).  If none of {cmd:xlo()},
{cmd:xhi()}, {cmd:ylo()}, {cmd:yhi()} are specified, all variables are
treated as crisp and the results are identical to standard DEA at every
α-level.

{dlgtab:Model}

{phang}
{cmd:alpha(}{it:numlist}{cmd:)} specifies the α-cut levels at which to
evaluate fuzzy efficiency.  Each value must be in [0,1].  Multiple values
are supplied as a space-separated list; the list is automatically sorted.
The default is {cmd:alpha(0 0.25 0.5 0.75 1)}.  Use {cmd:alpha(1)} to
obtain only the standard DEA result.

{phang}
{cmd:rts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
{cmd:crs} (constant, the default; synonym {cmd:ccr}) imposes no constraint
on the intensity variables λ.  {cmd:vrs} (variable; synonym {cmd:bcc})
adds sum(λ) = 1.  {cmd:drs} (decreasing) adds sum(λ) ≤ 1.
{cmd:irs} (increasing) adds sum(λ) ≥ 1.  The same RTS assumption is
applied to both the pessimistic and optimistic LPs at every α-level.

{phang}
{cmd:orientation(}{it:type}{cmd:)} or {cmd:ort(}{it:type}{cmd:)} specifies
the model orientation.  {cmd:input} (the default) minimizes the radial
contraction factor θ, giving input efficiency scores in (0,1].
{cmd:output} maximizes the radial expansion factor φ ≥ 1, stored as
1/φ ∈ (0,1].  Abbreviations {cmd:i} and {cmd:o} are accepted.

{dlgtab:Output}

{phang}
{cmd:generate(}{it:prefix}{cmd:)} specifies the variable name prefix for
the efficiency results.  The default prefix is {cmd:fdea_eff}.  For each
α-level, two variables are created:

{p2colset 9 42 44 2}{...}
{p2col :{it:prefix}{cmd:_pess_a}{it:N}}pessimistic efficiency at α = N/100{p_end}
{p2col :{it:prefix}{cmd:_opt_a}{it:N}}optimistic efficiency at α = N/100{p_end}
{p2colreset}{...}

{pstd}
where {it:N} = round(α × 100).  For example, {cmd:generate(eff)} with
{cmd:alpha(0 0.5 1)} creates {cmd:eff_pess_a0}, {cmd:eff_opt_a0},
{cmd:eff_pess_a50}, {cmd:eff_opt_a50}, {cmd:eff_pess_a100}, and
{cmd:eff_opt_a100}.

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies a string or numeric variable
containing DMU identifiers used in output displays.  If not specified,
DMUs are labelled {cmd:dmu1}, {cmd:dmu2}, etc.

{phang}
{cmd:detail} displays a DMU-level table of pessimistic and optimistic
efficiency scores at every α-level.

{dlgtab:Technical}

{phang}
{cmd:tolerance(}{it:#}{cmd:)} specifies the numerical tolerance for LP
convergence.  The default is {cmd:1e-7}.  This value is passed to the
Stata 16+ {cmd:LinearProgram()} interior point solver.

{phang}
{cmd:debug} displays additional status messages during execution.


{marker remarks}{...}
{title:Remarks}

{pstd}
{it:Relationship to IDEA.}
{cmd:dea_fuzzy} generalizes the Imprecise DEA ({cmd:model(imprecise)}) in
{cmd:dea} to a parametric family of interval problems indexed by α.
At a single α-level the LP structure is identical to IDEA: pessimistic
uses (X_U[k], Y_L[k]) for DMU {it:k} and (X_L, Y_U) for the frontier,
while optimistic reverses the roles.  The novelty of {cmd:dea_fuzzy} is
the α-cut parameterization, which reduces the full fuzzy number to a
linear interpolation between the support bounds and the modal value.

{pstd}
{it:Interpreting efficiency intervals.}
For input-oriented models, the efficiency interval [θ_L(α), θ_U(α)] has
the following interpretation:

{p2colset 9 36 38 2}{...}
{p2col :{it:Condition}}{it:Interpretation}{p_end}
{p2col :{hline 60}}{p_end}
{p2col :θ_L(α) = 1}DMU is {it:certainly efficient} at level α{p_end}
{p2col :θ_L(α) < 1, θ_U(α) ≥ 1}DMU is {it:possibly efficient} at level α{p_end}
{p2col :θ_U(α) < 1}DMU is {it:certainly inefficient} at level α{p_end}
{p2colreset}{...}

{pstd}
The same classification applies for output-oriented models using 1/φ in
place of θ.

{pstd}
{it:Membership function.}
Plotting θ_L(α) and θ_U(α) against α traces the lower and upper branches
of the fuzzy efficiency membership function.  At α = 1 both branches
meet at the crisp DEA score.  At α = 0 the interval is widest, reflecting
maximum data uncertainty.

{pstd}
{it:LP formulation — input-oriented.}
At each α-level, {cmd:dea_fuzzy} solves the following two LPs for each
DMU {it:k}:

{p 8 8 2}
{it:Pessimistic:} {space 2} min θ{break}
subject to  X_W(α)′ λ − θ × x_{k,U}(α) ≤ 0  {space 2} (m input constraints){break}
{space 12}  −Y_W(α)′ λ {space 12} ≤ −y_{k,L}(α)  (s output constraints){break}
{space 12}  [convexity constraint per {cmd:rts()}]{break}
{space 12}  λ ≥ 0,  θ ≥ 0

{p 8 8 2}
{it:Optimistic:} {space 3} min θ{break}
subject to  X_B(α)′ λ − θ × x_{k,L}(α) ≤ 0{break}
{space 12}  −Y_B(α)′ λ {space 12} ≤ −y_{k,U}(α){break}
{space 12}  [convexity constraint per {cmd:rts()}]{break}
{space 12}  λ ≥ 0,  θ ≥ 0

{pstd}
where X_W/Y_W use favorable bounds for peers and DMU {it:k}'s unfavorable
realization for DMU {it:k}; X_B/Y_B use unfavorable bounds for peers and
DMU {it:k}'s favorable realization for DMU {it:k}. This is the
Kao-Liu/deaR convention for best/worst alpha-cut evaluation.

{pstd}
{it:LP formulation — output-oriented.}
The output-oriented model maximizes φ (stored as 1/φ):

{p 8 8 2}
{it:Pessimistic:} {space 2} max φ  (equivalently, min −φ){break}
subject to  X_W(α)′ λ {space 14} ≤ x_{k,U}(α){break}
{space 12}  −Y_W(α)′ λ + φ × y_{k,L}(α) ≤ 0{break}
{space 12}  [convexity constraint per {cmd:rts()}]{break}
{space 12}  λ ≥ 0,  φ ≥ 1

{p 8 8 2}
{it:Optimistic:} {space 3} max φ{break}
subject to  X_B(α)′ λ {space 14} ≤ x_{k,L}(α){break}
{space 12}  −Y_B(α)′ λ + φ × y_{k,U}(α) ≤ 0{break}
{space 12}  [convexity constraint per {cmd:rts()}]{break}
{space 12}  λ ≥ 0,  φ ≥ 1

{pstd}
{it:Solver.}
All LPs are solved using Stata 16's {cmd:LinearProgram()} interior point
method.  The {cmd:dea_fuzzy} command therefore requires Stata 16.0 or later.

{pstd}
{it:Number of LPs solved.}
For {it:n} DMUs, {it:K} α-levels, and two orientations per level
(pessimistic and optimistic), the total number of LPs solved is
2 × {it:n} × {it:K}.  For example, with 50 DMUs and 5 α-levels,
500 LPs are solved.

{pstd}
{it:Verifying results at α = 1.}
Because α = 1 collapses all bounds to the center value, results at
α = 1 are identical for pessimistic and optimistic, and should
match standard DEA results from {cmd:dea} with the same {cmd:rts()} and
{cmd:orientation()}.  This provides a useful cross-check:

{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) ort(i) gen(crs_eff)}{p_end}
{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) alpha(1) gen(f)}{p_end}
{phang2}{cmd:. assert abs(crs_eff - f_pess_a100) < 1e-6}{p_end}


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


    {title:Setup}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str7 dmu i_employee i_area o_sales o_profits}{p_end}
{phang2}{cmd:. "store A" 10 20 70 6}{p_end}
{phang2}{cmd:. "store B" 15 15 100 3}{p_end}
{phang2}{cmd:. "store C" 20 30 80 5}{p_end}
{phang2}{cmd:. "store D" 25 15 100 2}{p_end}
{phang2}{cmd:. "store E" 12 9 90 8}{p_end}
{phang2}{cmd:. "store F" 16 18 75 5}{p_end}
{phang2}{cmd:. end}{p_end}

{phang2}{cmd:. // Variables: dmu, i_employee, i_area, o_sales, o_profits}{p_end}
{phang2}{cmd:. // Create ±10% triangular fuzzy spreads around observed values}{p_end}
{phang2}{cmd:. foreach v of varlist i_employee i_area o_sales o_profits {c -(}}{p_end}
{phang2}{cmd:. {space 4}gen `v'_lo = `v' * 0.90}{p_end}
{phang2}{cmd:. {space 4}gen `v'_hi = `v' * 1.10}{p_end}
{phang2}{cmd:. {c )-}}{p_end}

    {title:Basic fuzzy DEA: CRS, input-oriented, 5 α-levels}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) rts(crs) ort(input) gen(crs_f)}{p_end}

    {title:VRS fuzzy DEA with DMU labels and detailed table}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) rts(vrs) ort(input) dmu(dmu) gen(vrs_f) detail}{p_end}

    {title:Custom α-levels}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) alpha(0 0.2 0.4 0.6 0.8 1) gen(f6)}{p_end}

    {title:Output-oriented (stored as 1/φ)}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) rts(crs) ort(output) gen(out_f)}{p_end}

    {title:Fuzzy inputs only (crisp outputs)}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) gen(fi)}{p_end}

    {title:Verification: α=1 should equal standard DEA}

{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) ort(i) gen(theta)}{p_end}
{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) alpha(1) gen(chk)}{p_end}
{phang2}{cmd:. list dmu theta chk_pess_a100 chk_opt_a100, clean}{p_end}

    {title:Compute a finer membership grid}

{phang2}{cmd:. dea_fuzzy, inputs(i_employee i_area) outputs(o_sales o_profits) xlo(i_employee_lo i_area_lo) xhi(i_employee_hi i_area_hi) ylo(o_sales_lo o_profits_lo) yhi(o_sales_hi o_profits_hi) alpha(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1) gen(mf) dmu(dmu)}{p_end}
{phang2}{cmd:. list dmu mf_pess_a0 mf_opt_a0 mf_pess_a100 mf_opt_a100, clean}{p_end}



{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea_fuzzy} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of DMUs{p_end}
{synopt:{cmd:r(m_in)}}number of inputs{p_end}
{synopt:{cmd:r(m_out)}}number of outputs{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}}{cmd:dea_fuzzy}{p_end}
{synopt:{cmd:r(orientation)}}orientation ({cmd:input} or {cmd:output}){p_end}
{synopt:{cmd:r(rts)}}returns to scale ({cmd:crs}, {cmd:vrs}, {cmd:drs}, or {cmd:irs}){p_end}
{synopt:{cmd:r(alpha)}}α-levels used{p_end}
{synopt:{cmd:r(generate)}}generate prefix used{p_end}

{pstd}
Generated variables (with prefix {it:P} = {cmd:r(generate)}, α-level N = round(α×100)):

{synoptset 28 tabbed}{...}
{synopt:{it:P}{cmd:_pess_a}{it:N}}pessimistic efficiency at α = N/100{p_end}
{synopt:{it:P}{cmd:_opt_a}{it:N}}optimistic efficiency at α = N/100{p_end}
{synoptline}


{marker references}{...}
{title:References}

{pstd}
{it:Fuzzy DEA}

{marker KL2000}{...}
{phang}
Kao, C., and S.-T. Liu. 2000.
Fuzzy efficiency measures in data envelopment analysis.
{it:Fuzzy Sets and Systems} 113: 427-437.

{marker LFJ2003}{...}
{phang}
Lertworasirikul, S., S.-C. Fang, J. A. Joines, and H. L. W. Nuttle. 2003.
Fuzzy data envelopment analysis (DEA): a possibility approach.
{it:Fuzzy Sets and Systems} 139: 379-394.

{marker GT2001}{...}
{phang}
Guo, P., and H. Tanaka. 2001.
Fuzzy DEA: a perceptual evaluation method.
{it:Fuzzy Sets and Systems} 119: 149-160.

{pstd}
{it:Imprecise DEA (related)}

{marker CPY1999}{...}
{phang}
Cooper, W. W., K. S. Park, and G. Yu. 1999.
IDEA and AR-IDEA: Models for dealing with imprecise data in DEA.
{it:Management Science} 45(4): 597-607.
{browse "https://doi.org/10.1287/mnsc.45.4.597":doi:10.1287/mnsc.45.4.597}

{pstd}
{it:Foundational DEA}

{marker CCR1978}{...}
{phang}
Charnes, A., W. W. Cooper, and E. Rhodes. 1978.
Measuring the efficiency of decision making units.
{it:European Journal of Operational Research} 2: 429-444.

{marker BCC1984}{...}
{phang}
Banker, R. D., A. Charnes, and W. W. Cooper. 1984.
Some models for estimating technical and scale inefficiencies in data
envelopment analysis.
{it:Management Science} 30: 1078-1092.

{pstd}
{it:Software}

{marker JL2010}{...}
{phang}
Ji, Y.-B., and C. Lee. 2010.
Data envelopment analysis.
{it:Stata Journal} 10(2): 267-280.


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
{synopt:{helpb dea_fuzzy}}Fuzzy DEA{p_end}
{synopt:{helpb dea_game}}Game-theoretic DEA overview{p_end}
{synopt:{helpb dea_malmquist}}Malmquist productivity index{p_end}
{synopt:{helpb dea_nash}}Nash bargaining DEA{p_end}
{synopt:{helpb dea_network}}Network DEA{p_end}
{synopt:{helpb dea_sbmundesirable}}SBM with undesirable outputs{p_end}
{synopt:{helpb dea_stackelberg}}Stackelberg game DEA{p_end}
{synopt:{helpb dea_sw2007}}Simar-Wilson two-stage DEA{p_end}
{synopt:{helpb dea_window}}Window analysis DEA{p_end}
{synoptline}
