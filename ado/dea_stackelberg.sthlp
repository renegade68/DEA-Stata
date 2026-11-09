{smcl}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* version 1.3.0  09nov2026}{...}
{viewerjumpto "Related article" "dea_stackelberg##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_network" "help dea_network"}{...}
{vieweralsosee "dea_crossgame" "help dea_crossgame"}{...}
{title:Title}

{phang}
{bf:dea_stackelberg} {hline 2} Two-stage DEA with optional CRS leader-follower optimization

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:dea_stackelberg} {ifin}{cmd:,}
{cmd:inputs(}{it:varlist}{cmd:)}
{cmd:intermediates(}{it:varlist}{cmd:)}
{cmd:outputs(}{it:varlist}{cmd:)}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt inputs(varlist)}}Stage 1 input variables (X){p_end}
{synopt:{opt intermediates(varlist)}}intermediate variables linking stages (Z){p_end}
{synopt:{opt outputs(varlist)}}Stage 2 output variables (Y){p_end}

{syntab:Optional}
{synopt:{opt dmu(varname)}}DMU identifier variable{p_end}
{synopt:{opt leader(string)}}leader stage: {cmd:stage1} (default) or {cmd:stage2}{p_end}
{synopt:{cmd:method(}{it:string}{cmd:)}}{cmd:independent} (default) or {cmd:noncooperative}{p_end}
{synopt:{cmd:matrix}}return component, status, weight and reference matrices{p_end}
{synopt:{opt rts(string)}}{cmd:crs} (default) or {cmd:vrs}{p_end}
{synopt:{opt generate(string)}}prefix for generated variables; default {cmd:_stack}{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}
{synopt:{opt tolerance(#)}}positive finite LP residual tolerance; default {cmd:1e-7}{p_end}
{synopt:{opt detail}}show benchmark and stage availability diagnostics{p_end}
{synopt:{opt debug}}display debug information{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_stackelberg} evaluates the two-stage structure X -> Z -> Y.
All first-stage outputs are the only second-stage inputs. Two explicitly
selected methods are available.

{pstd}
{cmd:method(independent)} preserves the historical default: compute input-oriented
DEA separately for X -> Z and Z -> Y, label the selected stage as leader, and
report the product. Changing {cmd:leader()} changes labels only. These scores
are independent-stage diagnostics and do not estimate a leader-follower game.
CRS and VRS are supported for this method.

{pstd}
{cmd:method(noncooperative)} implements the CRS leader-follower DEA approach of
Liang, Cook, and Zhu (2008). First obtain the selected leader's independent
optimal efficiency. Then hold that efficiency fixed and maximize the follower's
efficiency using common intermediate weights. The follower may choose among
weights that preserve the leader optimum; an arbitrary leader weight vector
is not frozen. The two leader choices can yield different decompositions.
This option currently requires {cmd:rts(crs)}; VRS requests return error 198.

{pstd}
If only one intermediate column is positive anywhere in the selected sample,
the CRS model reduces to independent stage optima. The command uses this exact
model reduction before attempting a full joint LP, obtains shared weights from
two multiplier LPs, and checks them against the leader-fixed model. This is a
structural calculation route, not a replacement after an optimization failure.
Zero intermediate columns have zero weights in this route. With multiple active
intermediates, the full joint LP is used. The returned route identifies which
calculation was attempted.

{pstd}
With nonnegative weights (v,w,u), all reference DMUs satisfy wZ <= vX and
uY <= wZ. If Stage 1 leads, normalize wZ0=1, impose wZ0=e1*vX0, and maximize
uY0. If Stage 2 leads, normalize vX0=1, impose uY0=e2*wZ0, and maximize wZ0.
Here e1 or e2 is the corresponding independent leader optimum. These positive
normalizations are equivalent to the CRS ratio models. No intercept-based VRS
extension or fuzzy/bilevel-general-purpose optimization is implied.

{pstd}
Independent benchmarks are always returned. For the noncooperative method, a
failed follower LP leaves the leader and available benchmarks intact, while
follower and overall scores are missing. No independent follower score is
substituted. If the leader score is unavailable or nonpositive, no follower LP
is attempted. An unavailable independent follower benchmark does not block a
successful joint LP when the leader benchmark is available.

{pstd}
All quantities must be nonnegative, with at least one positive variable in
each of X, Z, and Y for every selected observation. Missing quantities and
explicitly missing IDs are excluded. The selected sample is also the reference
sample. An empty sample returns 2000. Invalid data, tolerance, or unsupported
method/RTS combinations are checked before replacing existing results.

{pstd}
Positive column and constraint scaling preserve the LP. Native solutions must
converge and pass constraint/bound checks at {cmd:tolerance()}; noncooperative benchmarks
and single-intermediate multiplier solves use an internal cap of 1e-12.
Independent-mode and full follower LPs use a cap of 1e-8, always no larger than
the requested tolerance. This native check alone is not an exact optimality certificate.
Scores are not clipped. The exact verification path and retained native
provenance are reported explicitly below.

{marker options}{...}
{title:Options}

{phang}
{opt leader(string)} specifies which stage leads. {cmd:stage1} (default) or
{cmd:stage2}.

{phang}
{opt detail} adds independent benchmarks, component availability, and follower
solve status to the stage table. {cmd:debug} prints per-observation scores.

{phang}
{opt generate(string)} sets the prefix for generated variables (default {cmd:_stack}).
The following variables are created:

{p2colset 9 26 28 2}
{p2col:{it:prefix}{cmd:_leader}}efficiency of the stage reported as leader{p_end}
{p2col:{it:prefix}{cmd:_follower}}efficiency of the stage reported as follower{p_end}
{p2col:{it:prefix}{cmd:_overall}}overall efficiency (leader × follower){p_end}
{p2col:{it:prefix}{cmd:_indep1}}Stage 1 efficiency under independent benchmark{p_end}
{p2col:{it:prefix}{cmd:_indep2}}Stage 2 efficiency under independent benchmark{p_end}

{marker results}{...}

{pstd}
{bf:Numerical quality and retained results.} Native {cmd:LinearProgram()}
results are checked first. A bounded exact-arithmetic path then verifies the
model formed from the original stored double data, its optimality certificate,
and the precision of exported scores and weights or reference coefficients.
Certification is attempted for at most 32 selected DMUs and 16 multiplier
variables (including two variables for a VRS intercept). Larger models or
unsuccessful certification retain the native result and label it unverified;
this does not prove infeasibility. No tolerance relaxation or score clipping
is applied. Exact certification concerns the stored binary inputs, not unknown
unrounded source quantities.

{pstd}
For independent formulations, certified stage objectives are combined as exact
fractions before a single conversion to double. If a positive product rounds
to zero, or the final product/mean fails the requested relative precision,
the overall score is missing and a warning is printed. Available stage scores
and reference coefficients are retained. If either stage is unverified, the
aggregate of the available stage values is retained as unverified. Linked
relational and noncooperative scores keep their own joint-model verification.
A genuine certified zero is distinct from a positive value lost at export.

{pstd}
With {cmd:matrix}, {cmd:r(score_quality)} and {cmd:r(score_reason)} have columns
overall, stage1, stage2. Quality codes are 0 for unavailable, 1 for finite but
unverified, 2 for certified with export checks passed, and 3 for an overall
value withheld after aggregation/export checks. {cmd:r(status)} continues to
report only numeric availability. Reason 0 means certification succeeded;
201 and 202 identify the DMU and multiplier-size limits; other stage reasons
identify model construction, certificate, arithmetic, or export failures.
Aggregate reason 1001 means at least one stage is unverified, 1002 invalid
internal composition arguments, 1101 positive underflow, 1102 invalid/range
export, 1103 insufficient relative precision, and 1200 plus an arithmetic
error identifies bounded-arithmetic failure. These reasons are not proofs
that the DEA model has no feasible solution.

{pstd}
{cmd:r(native_efficiencies)} preserves the checked native scores before exact
verification. For every nonrelational formulation,
{cmd:r(independent_efficiencies)}, {cmd:r(independent_quality)}, and
{cmd:r(independent_reason)} describe the final independent Stage 1 and Stage 2
benchmarks. {cmd:r(native_independent)}, {cmd:r(native_lambda_stage1)}, and
{cmd:r(native_lambda_stage2)} preserve their native provenance. The final
{cmd:r(lambda_stage1)} and {cmd:r(lambda_stage2)} refer to the final independent
benchmarks; in a noncooperative model they are not the linked follower weights.
Native and certified reference sets can differ when alternative optima exist.

{pstd}
For independent formulations, {cmd:r(aggregate_diagnostics)} has columns
quality, reason, value, rel_error, operations, maxlimbs, arithmetic_error,
and legacy_value. The last column is the ordinary double aggregate of the
final stage scores before the exact aggregate check. Withheld values have a
missing value column. Linked models do not return this matrix.

{pstd}
The scalars {cmd:r(n_certified_overall)}, {cmd:r(n_certified_stage1)},
{cmd:r(n_certified_stage2)} and the corresponding {cmd:r(n_unverified_*)}
count quality 2 and quality 1 results; unavailable results are not counted as
finite unverified values. {cmd:r(n_withheld_overall)} counts quality 3.
Nonrelational models also return {cmd:r(n_certified_indep1)},
{cmd:r(n_certified_indep2)}, {cmd:r(n_unverified_indep1)}, and
{cmd:r(n_unverified_indep2)}. These counts are available without {cmd:matrix}.
{cmd:r(quality_policy)} is {cmd:certify_or_retain_unverified};
{cmd:r(aggregate_policy)} is {cmd:exact_or_withhold}, conditional on both
independent stages being certified. Internal diagnostics named t120, t121,
and t122 are private and not a stable interface.

{title:Stored results}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(n)}}number of DMUs{p_end}
{synopt:{cmd:r(tolerance)}}requested residual tolerance{p_end}
{synopt:{cmd:r(n_failed_overall)}}missing overall scores{p_end}
{synopt:{cmd:r(n_failed_stage1)}}missing reported Stage 1 scores{p_end}
{synopt:{cmd:r(n_failed_stage2)}}missing reported Stage 2 scores{p_end}
{synopt:{cmd:r(n_failed_indep1)}}missing independent Stage 1 benchmarks{p_end}
{synopt:{cmd:r(n_failed_indep2)}}missing independent Stage 2 benchmarks{p_end}
{synopt:{cmd:r(n_failed_follower_lp)}}attempted but unavailable follower solves, including the reduced route{p_end}
{synopt:{cmd:r(n_single_link)}}DMUs attempting the single-intermediate reduction{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(method)}}independent or noncooperative{p_end}
{synopt:{cmd:r(formulation)}}independent_ordering or noncooperative_crs{p_end}
{synopt:{cmd:r(orientation)}}input{p_end}
{synopt:{cmd:r(leader)}}leader stage used{p_end}
{synopt:{cmd:r(rts)}}returns to scale used{p_end}

{pstd}
With {cmd:matrix}, {cmd:r(efficiencies)} reports overall, stage1, stage2 in
selected observation order. {cmd:r(status)} uses the same columns: 1 available,
0 missing. {cmd:r(lambda_stage1)} and {cmd:r(lambda_stage2)} hold the independent
benchmark reference intensities (evaluator rows, reference-DMU columns).
For {cmd:method(noncooperative)}, {cmd:r(weights_input)},
{cmd:r(weights_intermediate)}, and {cmd:r(weights_output)} contain the joint
(v,w,u) in original units. They are missing if the joint solve is unavailable.
{cmd:r(follower_status)} is 0 if no positive finite leader score is available,
1 if the follower LP was attempted but unavailable, and 2 if accepted.
{cmd:r(follower_route)} is 0 when not attempted, 1 for the full joint LP, or
2 for the single-intermediate reduction, including unsuccessful attempts.
Weights need not be unique even if the optimal stage scores are unique.

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


{pstd}Supply chain: manufacturer (Stage 1) leads, retailer (Stage 2) follows{p_end}

{pstd}Setup{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str6 unit raw labor parts sales}{p_end}
{phang2}{cmd:"UnitA" 40 12 68 102}{p_end}
{phang2}{cmd:"UnitB" 55 15 76 116}{p_end}
{phang2}{cmd:"UnitC" 48 13 72 108}{p_end}
{phang2}{cmd:"UnitD" 62 18 84 126}{p_end}
{phang2}{cmd:"UnitE" 36 11 61  94}{p_end}
{phang2}{cmd:"UnitF" 70 20 91 138}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Stage 1 as leader{p_end}
{phang2}{cmd:. dea_stackelberg, inputs(raw labor) intermediates(parts) ///}{p_end}
{phang2}{cmd:    outputs(sales) dmu(unit) leader(stage1)}{p_end}

{pstd}Reverse: downstream stage leads{p_end}
{phang2}{cmd:. dea_stackelberg, inputs(raw labor) intermediates(parts) ///}{p_end}
{phang2}{cmd:    outputs(sales) dmu(unit) leader(stage2) rts(vrs) replacevars}{p_end}

{pstd}Compare Stackelberg vs. independent efficiency{p_end}
{phang2}{cmd:. dea_stackelberg, inputs(raw labor) intermediates(parts) ///}{p_end}
{phang2}{cmd:    outputs(sales) dmu(unit) replacevars}{p_end}
{phang2}{cmd:. list unit _stack_leader _stack_follower _stack_overall _stack_indep1 _stack_indep2}{p_end}


{pstd}CRS shared-weight leader-follower optimization{p_end}
{phang2}{cmd:. dea_stackelberg, inputs(raw labor) intermediates(parts) outputs(sales) dmu(unit) method(noncooperative) leader(stage1) gen(game1) matrix detail}{p_end}
{phang2}{cmd:. matrix list r(follower_status)}{p_end}
{phang2}{cmd:. dea_stackelberg, inputs(raw labor) intermediates(parts) outputs(sales) dmu(unit) method(noncooperative) leader(stage2) gen(game2) matrix detail}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker references}{...}
{title:References}

{phang}
Liang, L., W. D. Cook, and J. Zhu. 2008. DEA models for two-stage processes:
Game approach and efficiency decomposition. {it:Naval Research Logistics}
55(7): 643--653.
{browse "https://doi.org/10.1002/nav.20308":doi:10.1002/nav.20308}

{pstd}The following references provide broader supply-chain/network or fuzzy-model
background; they do not define every option implemented here.{p_end}

{phang}
Liang, L., F. Yang, W.D. Cook, and J. Zhu. 2006. DEA models for supply chain
efficiency evaluation. {it:Annals of Operations Research} 145(1): 35–49.
{browse "https://doi.org/10.1007/s10479-006-0026-7":doi:10.1007/s10479-006-0026-7}

{phang}
Li, Y., Y. Chen, L. Liang, and J. Xie. 2012. DEA models for extended two-stage
network structures. {it:Omega} 40(5): 611–618.
{browse "https://doi.org/10.1016/j.omega.2011.11.007":doi:10.1016/j.omega.2011.11.007}

{phang}
Tavana, M., and K. Khalili-Damghani. 2014. A new two-stage Stackelberg fuzzy
data envelopment analysis model. {it:Measurement} 53: 277–296.
{browse "https://doi.org/10.1016/j.measurement.2014.03.030":doi:10.1016/j.measurement.2014.03.030}

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
