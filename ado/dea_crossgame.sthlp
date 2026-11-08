{smcl}
{* Help revision 2026-11-02: validated CRS LPs, complete Nash states, matrix/detail diagnostics.}{...}
{* Help revision 2026-11-02: singleton sample rules by strategy.}{...}
{* Help revision 2026-11-01: if/in defines both evaluation and reference samples.}{...}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* version 2.4.0  08nov2026}{...}
{viewerjumpto "Related article" "dea_crossgame##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_crosseff" "help dea_crosseff"}{...}
{vieweralsosee "dea_nash" "help dea_nash"}{...}
{vieweralsosee "dea_cooperative" "help dea_cooperative"}{...}
{title:Title}

{phang}
{bf:dea_crossgame} {hline 2} Cross-Efficiency Game DEA

{marker syntax}{...}
{title:Syntax}

{pstd}Form A{p_end}

{p 8 17 2}
{cmd:dea_crossgame} {it:inputvars} {cmd:=} {it:outputvars}
{ifin} [{cmd:,} {it:options}]

{pstd}Form B{p_end}

{p 8 17 2}
{cmd:dea_crossgame} {ifin}{cmd:,}
{cmd:inputs(}{it:varlist}{cmd:)}
{cmd:outputs(}{it:varlist}{cmd:)}
[{it:options}]

{pstd}
{cmd:if} and {cmd:in} restrict both the DMUs evaluated and the reference
set. Observations with missing input/output values or a missing explicit DMU identifier
are excluded. Quantities must be numeric and nonnegative, with at least one
positive input and one positive output per selected observation. Unrelated
variables do not affect sample selection. Generated results are missing outside
the selected sample. An empty estimation sample returns {cmd:r(2000)}.
{cmd:strategy(aggressive)}, {cmd:strategy(benevolent)} and {cmd:strategy(all)}
(including the default) require at least two selected observations; a single
remaining DMU returns {cmd:r(2001)} before creating or replacing result columns.
When {cmd:replacevars} is used, it replaces entire result columns, including
setting observations outside the selected sample to missing.{p_end}

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt inputs(varlist)}}input variables{p_end}
{synopt:{opt outputs(varlist)}}output variables{p_end}

{syntab:Optional}
{synopt:{opt dmu(varname)}}DMU identifier variable{p_end}
{synopt:{opt strategy(string)}}game strategy: {cmd:aggressive}, {cmd:benevolent}, {cmd:nash}, or {cmd:all} (default){p_end}
{synopt:{opt maxiter(#)}}maximum iterations for Nash convergence; default {cmd:100}{p_end}
{synopt:{opt generate(string)}}prefix for generated variables; default {cmd:_cg}{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}
{synopt:{opt tolerance(#)}}positive finite LP feasibility and Nash stopping tolerance; default {cmd:1e-6}{p_end}
{synopt:{opt matrix}}display full cross-efficiency matrix{p_end}
{synopt:{opt detail}}display detailed output{p_end}
{synopt:{opt debug}}display debug information{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_crossgame} extends the cross-efficiency framework of
{help dea_crossgame##DG1994:Doyle and Green (1994)}
by incorporating strategic behavior among DMUs. In the standard cross-efficiency
model, each DMU chooses optimal weights primarily for self-appraisal. The game
DEA approach ({help dea_crossgame##LWCZ2008:Liang et al. 2008}) models how DMUs would choose weights if they
anticipated — and responded to — other DMUs' choices.

{pstd}
{ul:Three competitive strategies:}

{phang2}
{bf:Aggressive}: After primary input-oriented CRS efficiency is computed,
minimize sum over j!=i of u_i*y_j, subject to v_i*x_i=1,
u_i*y_i=self_i, frontier constraints u_i*y_j<=v_i*x_j, and nonnegative
weights. This is the command's existing weighted-output surrogate.

{phang2}
{bf:Benevolent}: Maximize the same weighted-output sum on the same optimal
self-efficiency face. Neither secondary objective optimizes the sum of peer
ratios u_i*y_j/(v_i*x_j), since other DMUs' denominators are not normalized.
These strategies are diagnostics; they are not the exact fractional secondary
models or a zero-sum/cooperative structural game.

{phang2}
{bf:Nash}: Implements the CRS input model (5) and synchronous iteration of
{help dea_crossgame##LWCZ2008:Liang et al. (2008)}. For each protected peer d
and target j, maximize u*y_j subject to v*x_j=1, u*Y<=v*X,
theta_d*v*x_d-u*y_d<=0, and u,v>=0. Weights can differ for every (d,j).
Initial theta is the column mean of the ordinary primary cross-efficiency
matrix, including its diagonal. Each new theta averages all n game LP
objectives for the target, including d=j.

{phang2}
{bf:All}: Compute all three strategies (default). A single strategy uses the
same primary solution and procedure as its corresponding part of {cmd:all}.

{pstd}
There is no implemented guarantee that aggressive scores are below benevolent
scores or that Nash lies between them. The secondary objectives and averaging
sets differ. Alternative optimal weights may change peer matrices and ranks.
Nash's numerical stopping rule does not certify an exact equilibrium.

{marker options}{...}
{title:Options}

{phang}
{opt strategy(string)} specifies which game strategy to use.
{cmd:aggressive}, {cmd:benevolent}, {cmd:nash}, or {cmd:all} (default).
Aggressive and benevolent scores average evaluations by other DMUs, excluding
the evaluated DMU itself, so they require at least two observations after
{cmd:if}/{cmd:in} and missing-value exclusions. The count is based on selected
observations, not distinct ID labels or distinct input/output vectors.
The existing {cmd:strategy(nash)} calculation includes self-evaluations and
remains available for one DMU; interpret that case as self-appraisal.

{phang}
{opt maxiter(#)} sets the maximum number of iterations for Nash equilibrium
computation. Default is 100. If convergence is not achieved within {opt maxiter()}
iterations, the last complete state is returned with {cmd:r(nash_status)=1}.
The value must be a positive integer, even when only a secondary strategy is
requested. Invalid values return error 198 before result replacement.

{phang}
{opt matrix} displays and returns the full n×n matrix for each requested
strategy. Rows and columns follow selected observation order. For aggressive
and benevolent matrices, row i is an evaluator and column j is a target.
In the Nash matrix, row d is the protected peer and column j is the target;
each cell has its own optimal weights. The initial matrix, selected ordinary
and secondary weights, statuses, peer counts and Nash history are also returned.
Weights use original data units. Stata matrix-size limits apply.

{phang}
{opt detail} displays secondary stage statuses and available-peer counts, and
the Nash trace when requested. {opt debug} additionally prints each completed
Nash iteration's maximum change.

{phang}
{opt tolerance(#)} must be positive and finite. LPs are column/evaluated-row
reparameterized and constraint rows positively rescaled. The native solver uses
min(tolerance,1e-12) for convergence and the requested tolerance for scaled
primal/bound checks. The same requested value sets the Nash stopping criterion:
max absolute change of all scores between two complete states < tolerance.
This is not a bound on score error or a proof of exact optimality. No solver
engine option is provided by this command.

{phang}
{opt generate(string)} sets the prefix for generated variables (default {cmd:_cg}).
The following variables are created depending on {opt strategy()}:

{p2colset 9 22 24 2}
{p2col:{it:prefix}{cmd:_agg}}aggressive cross-efficiency score{p_end}
{p2col:{it:prefix}{cmd:_ben}}benevolent cross-efficiency score{p_end}
{p2col:{it:prefix}{cmd:_nash}}Nash equilibrium cross-efficiency score{p_end}

{marker recovery}{...}
{title:Certified recovery after a native failure}

{pstd}
Version 2.4.0 preserves accepted native results and attempts a bounded recovery
only when a secondary or game LP is unavailable. The common simplex proposes a
basis; exact arithmetic verifies the primal constraints, dual signs and objective
before a recovery is accepted. This certification uses the original binary64
input values. Exported double weights are checked against the requested tolerance.
An ordinary native success is not relabeled as an exact certificate.{p_end}

{pstd}
For a failed aggressive/benevolent solve, the original-data primary LP is
certified and its optimal face is used for the secondary problem. A rounded
primary efficiency is not inserted again as an equality RHS. The original
reported self score is retained. The diagnostic fields {cmd:self_certified} and
{cmd:self_delta} show the rounded certified target and its displayed difference
from {cmd:self_initial}; a zero displayed difference need not imply exact equality.
The weighted-output secondary objective is unchanged. For a failed game LP,
z=v*x_d and theta_d*z-u*y_d<=0 avoid rounding theta_d*x_d before solving.
The same theta and original constraints are retained.{p_end}

{pstd}
Recovery is limited to 32 selected observations, 16 variables (including the
additional game variable), and 256 recovery attempts per command. Existing exact
solver work, limb and pivot limits apply. Unresolved results retain the primary
fallback or last complete Nash state; no constraint relaxation or tolerance increase
is used. A native warning about infeasibility describes that solver attempt.
The final recovery diagnostics and Nash status determine the returned result;
an unresolved recovery is not a mathematical infeasibility certificate.{p_end}

{pstd}
Recovery counts are returned even without {cmd:matrix} or {cmd:detail}. When
records exist and {cmd:c(matsize)>=32}, {cmd:r(recovery_diagnostics)} stores the
first min(256,c(matsize)) eligible failed-call records, with selected-observation
indices and original observation numbers. It never stores an arbitrary numeric
encoding of a string DMU label. Counts include records omitted by the attempt or
matrix-size limit. With smaller matrix sizes only aggregate diagnostics are
returned. No diagnostic files or command-specific permanent Mata counters are created.{p_end}

{pstd}
Additional scalars: {cmd:r(recovery_version)=11401}, {cmd:r(n_recovery_attempted)},
{cmd:r(n_recovered)}, {cmd:r(n_recovery_unresolved)}, {cmd:r(n_recovery_limit)},
{cmd:r(n_native_failed_agg)}, {cmd:r(n_native_failed_ben)}, {cmd:r(n_native_failed_game)},
{cmd:r(n_recovered_agg)}, {cmd:r(n_recovered_ben)}, {cmd:r(n_recovered_game)},
{cmd:r(n_primary_target_changed)}, {cmd:r(max_primary_target_change)},
{cmd:r(max_recovery_residual)}, {cmd:r(n_recovery_rows)}, and
{cmd:r(n_recovery_rows_omitted)}. Unresolved includes attempts skipped by the
256-attempt limit. Target-change counts cover accepted secondary recoveries with
a nonzero displayed delta. The maximum residual covers accepted recovered points.
Game recovery counts include solves in an iteration that may subsequently be
discarded; {cmd:iteration_accepted} discloses that distinction.{p_end}

{pstd}
Diagnostic columns, in order:{break}
{cmd:iteration evaluator observation mode peer peer_observation}{break}
{cmd:native_converged native_error native_value native_rc}{break}
{cmd:self_initial self_certified self_delta theta recovered reason}{break}
{cmd:primary_cert primary_reason secondary_cert secondary_reason}{break}
{cmd:primary_ops secondary_ops max_limbs primary_pivots secondary_pivots}{break}
{cmd:export_eq export_ineq export_gap export_bounds export_error}{break}
{cmd:objective_value iteration_accepted}.{p_end}

{pstd}
{cmd:mode}: 1 aggressive, 2 benevolent, 3 game. Secondary iteration is 0.
The {cmd:evaluator} index denotes the evaluator in secondary solves and target j
in a game solve; {cmd:peer} is the protected peer d in the game matrix.
{cmd:native_value} is the optimizer's objective before command sign conversion;
{cmd:native_rc}=1 means unavailable or rejected by the native acceptance checks.
{cmd:recovered}: 1 accepted certified recovery, 0 unresolved.
{cmd:reason}: 0 accepted; 1 domain/dimension limit; 2 exact mapping unavailable;
3 primary/game certificate unavailable; 4 exact secondary objective build unavailable;
5 secondary certificate unavailable; 6 point export/arithmetic/residual unavailable;
7 objective export/range unavailable. Certificate fields equal 1 only when verified;
0/-1 mean unverified/unavailable, not infeasible. Additional reason/operation/pivot
fields describe the bounded certificate procedure; unused fields are missing.
For game records {cmd:iteration_accepted} is 1 only if every pair in that iteration
was accepted; 0 denotes a discarded update. It is missing for secondary records.
{cmd:export_eq}, {cmd:export_ineq}, and {cmd:export_gap} are scaled exact residuals
of the actual exported double point, and {cmd:export_bounds} reports nonnegativity.
{cmd:objective_value} uses the secondary minimization sign (negative weighted-output
sum for benevolent) or the positive game score.{p_end}

{pstd}
With {cmd:matrix}, {cmd:r(source_agg)} and {cmd:r(source_ben)} distinguish
0 primary unavailable, 1 native secondary accepted, 2 certified secondary recovery,
and 3 primary fallback. Existing {cmd:status_agg}/{cmd:status_ben} retain their
0/2/3 meanings. Optional matrices exist only for requested strategies. New
recovery records are command-local, so subsequent calls do not inherit them.{p_end}

{pstd}
All companion ado files must be installed together. Dependencies are checked
before result replacement. Solver state is restored on load or compute errors;
existing result columns are replaced only after computation completes. An error
or user Break during computation leaves existing result columns intact.{p_end}

{marker results}{...}
{title:Stored results and failure handling}

{pstd}
From version 2.3.0 the fixed 1e-6 positive weight floor is removed; all weights
have lower bound zero. This prevents a unit-dependent change in the model.
Historical weights and peer scores may change. A finite nonnegative weighted
output, a finite positive weighted input, and a finite ratio <=1+tolerance
are required for ordinary and secondary peer cells. Small positive denominators
are retained. No ratio clipping or epsilon imputation is used. Individual zeros
and redundant all-zero columns are allowed; 0/0 remains missing.

{pstd}
Aggressive/benevolent means use available off-diagonal cells. If a secondary
LP remains unavailable after the bounded recovery, the already accepted primary weights are retained with
status 3. This is a fallback, not a successful secondary optimization. An
unavailable primary leaves that evaluator's weights missing, although the DMU
can still receive evaluations from other evaluators. Peer counts disclose
changing averaging sets.

{pstd}
Nash requires a complete initial matrix and complete n-by-n iterations.
An incomplete initial matrix returns missing Nash scores (status 2). If any
pair LP fails in an iteration, the entire attempted update is discarded and
the last complete state is returned (status 3). That state can be the initial
ordinary cross-efficiency vector when zero game iterations succeeded. Neither
status 1 nor 3 should be interpreted as an equilibrium estimate. Primary and
secondary results remain available. An unavailable solve is not a certificate
of mathematical infeasibility.

{synoptset 27 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:r(n)}}selected DMUs{p_end}
{synopt:{cmd:r(tolerance)}}requested tolerance{p_end}
{synopt:{cmd:r(maxiter)}}iteration limit{p_end}
{synopt:{cmd:r(n_failed_primary)}}unavailable primary solves{p_end}
{synopt:{cmd:r(n_failed_agg)}, {cmd:r(n_failed_ben)}}secondary fallbacks; requested strategies only{p_end}
{synopt:{cmd:r(n_missing_agg)}, {cmd:r(n_missing_ben)}}missing off-diagonal cells; requested strategies only{p_end}
{synopt:{cmd:r(nash_status)}}0 converged; 1 iteration limit; 2 initial matrix incomplete; 3 iteration failed{p_end}
{synopt:{cmd:r(nash_converged)}}1 only for status 0{p_end}
{synopt:{cmd:r(nash_iterations)}}number of completed game iterations{p_end}
{synopt:{cmd:r(nash_attempts)}}number of attempted game iterations{p_end}
{synopt:{cmd:r(nash_failed_pairs)}}unavailable pairs in the failed attempted iteration; otherwise 0{p_end}
{synopt:{cmd:r(nash_initial_missing)}}missing cells in initial n-by-n matrix{p_end}
{synopt:{cmd:r(nash_last_change)}}last completed iteration's maximum change; missing if none{p_end}
{syntab:Macros}
{synopt:{cmd:r(strategy)}}requested strategy{p_end}
{syntab:Matrices: require matrix option}
{synopt:{cmd:r(self_eff)}}n-by-1 ordinary self scores{p_end}
{synopt:{cmd:r(weights_input)}, {cmd:r(weights_output)}}n-by-m/n-by-s primary weights{p_end}
{synopt:{cmd:r(ce_initial)}}ordinary primary cross-efficiency matrix{p_end}
{synopt:{cmd:r(ce_agg)}, {cmd:r(ce_ben)}}requested secondary peer matrices{p_end}
{synopt:{cmd:r(weights_agg_input)}, {cmd:r(weights_agg_output)}}selected aggressive weights{p_end}
{synopt:{cmd:r(weights_ben_input)}, {cmd:r(weights_ben_output)}}selected benevolent weights{p_end}
{synopt:{cmd:r(status_agg)}, {cmd:r(status_ben)}}n-by-1: 0 primary unavailable; 2 secondary accepted; 3 primary fallback{p_end}
{synopt:{cmd:r(count_agg)}, {cmd:r(count_ben)}}n-by-1 available off-diagonal counts{p_end}
{synopt:{cmd:r(ce_nash)}}last complete state matrix; all missing for status 2{p_end}
{synopt:{cmd:r(nash_history)}}score vectors by row, initial state then completed iterations; n columns{p_end}
{synopt:{cmd:r(nash_trace)}}initial state then all attempts: iteration, max_change, missing_pairs, complete{p_end}
{p2colreset}{...}

{pstd}
Nash scalars/matrices exist only when requested. Secondary-specific results
exist only for their requested strategy. At trace iteration 0, missing_pairs
counts initial missing ratios and complete indicates whether initialization
succeeded. The trace retains a failed attempt, while the score history does not.

{pstd}
{bf:Protecting existing data.} All output names are checked before results are
created. An existing result variable causes an error unless {cmd:replacevars}
is specified. The option applies only to variables produced by the current
command; it does not remove other variables sharing the prefix. Replacement
recreates whole result columns, including their storage types and labels;
observations not scored by the command do not retain old scores. Variables
used as model data or identifiers cannot be replaced, even with
{cmd:replacevars}. Invalid or duplicate output names are rejected before any
existing result column is changed. Name checks run before computation. Result columns are committed only after
computation succeeds; existing outputs survive errors during that computation.

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


{pstd}Setup{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str6 dmu capital labor output}{p_end}
{phang2}{cmd:"Firm A" 12  8  95}{p_end}
{phang2}{cmd:"Firm B" 18 12 125}{p_end}
{phang2}{cmd:"Firm C" 15 10 118}{p_end}
{phang2}{cmd:"Firm D" 22 15 132}{p_end}
{phang2}{cmd:"Firm E" 10  7  90}{p_end}
{phang2}{cmd:"Firm F" 20 13 128}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Run all three strategies (default){p_end}
{phang2}{cmd:. dea_crossgame, inputs(capital labor) outputs(output) dmu(dmu)}{p_end}

{pstd}Nash equilibrium only{p_end}
{phang2}{cmd:. dea_crossgame, inputs(capital labor) outputs(output) dmu(dmu) strategy(nash) maxiter(200) replacevars}{p_end}

{pstd}Compare aggressive vs. benevolent{p_end}
{phang2}{cmd:. dea_crossgame, inputs(capital labor) outputs(output) dmu(dmu) strategy(all) replacevars}{p_end}
{phang2}{cmd:. list dmu _cg_agg _cg_ben _cg_nash}{p_end}
{phang2}{cmd:. return list}{p_end}
{phang2}{cmd:. assert r(recovery_version)==11401}{p_end}
{phang2}{cmd:. if r(n_recovery_rows)>0 matrix list r(recovery_diagnostics)}{p_end}

{pstd}Display cross-efficiency matrix for aggressive strategy{p_end}
{phang2}{cmd:. dea_crossgame, inputs(capital labor) outputs(output) dmu(dmu) strategy(aggressive) matrix detail replacevars}{p_end}


{pstd}Inspect a bounded Nash run and its stopping status{p_end}
{phang2}{cmd:. dea_crossgame capital labor = output, dmu(dmu) strategy(nash) gen(short) maxiter(1) matrix detail}{p_end}
{phang2}{cmd:. assert inlist(r(nash_status),0,1,2,3)}{p_end}
{phang2}{cmd:. matrix list r(nash_trace)}{p_end}

{pstd}Restrict both evaluation and reference DMUs{p_end}
{phang2}{cmd:. dea_crossgame capital labor = output if inlist(dmu,"Firm B","Firm C","Firm D","Firm E"), dmu(dmu) gen(sub)}{p_end}
{phang2}{cmd:. assert r(n)==4}{p_end}
{phang2}{cmd:. assert missing(sub_agg) if dmu=="Firm A" | dmu=="Firm F"}{p_end}

{pstd}Single-DMU behavior by strategy:{p_end}
{phang2}{cmd:. capture noisily dea_crossgame capital labor = output if dmu=="Firm A", strategy(all) gen(one)}{p_end}
{phang2}{cmd:. assert _rc==2001}{p_end}
{phang2}{cmd:. capture confirm variable one_agg, exact}{p_end}
{phang2}{cmd:. assert _rc==111}{p_end}
{phang2}{cmd:. dea_crossgame capital labor = output if dmu=="Firm A", strategy(nash) gen(one)}{p_end}
{phang2}{cmd:. assert r(n)==1}{p_end}
{phang2}{cmd:. assert abs(one_nash-1)<1e-6 if dmu=="Firm A"}{p_end}
{phang2}{cmd:. assert missing(one_nash) if dmu!="Firm A"}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker references}{...}
{title:References}

{marker DG1994}{...}
{phang}
Doyle, J., and R. Green. 1994. Efficiency and cross-efficiency in DEA: Derivations,
meanings and uses. {it:Journal of the Operational Research Society} 45(5): 567–578.
{browse "https://doi.org/10.1057/jors.1994.84":doi:10.1057/jors.1994.84}

{marker LWCZ2008}{...}
{phang}
Liang, L., J. Wu, W.D. Cook, and J. Zhu. 2008. The DEA game cross-efficiency model
and its Nash equilibrium. {it:Operations Research} 56(5): 1278–1288.
{browse "https://doi.org/10.1287/opre.1070.0487":doi:10.1287/opre.1070.0487}

{marker WLY2009}{...}
{phang}
Wu, J., L. Liang, and F. Yang. 2009. Achievement and benchmarking of countries
at the Summer Olympics using cross efficiency evaluation method.
{it:European Journal of Operational Research} 197(2): 722–730.
{browse "https://doi.org/10.1016/j.ejor.2008.06.030":doi:10.1016/j.ejor.2008.06.030}

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
