{smcl}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-10-22.3: dea journal attribution and software-version notice added.}{...}
{* *! version 1.3.0  09nov2026}{...}
{viewerjumpto "Related article" "dea_network##publication"}{...}
{vieweralsosee "dea" "help dea"}{...}
{vieweralsosee "dea_bootstrap" "help dea_bootstrap"}{...}
{vieweralsosee "dea_crosseff" "help dea_crosseff"}{...}
{vieweralsosee "dea_malmquist" "help dea_malmquist"}{...}
{vieweralsosee "dea_window" "help dea_window"}{...}
{viewerjumpto "Syntax" "dea_network##syntax"}{...}
{viewerjumpto "Description" "dea_network##description"}{...}
{viewerjumpto "Options" "dea_network##options"}{...}
{viewerjumpto "Remarks" "dea_network##remarks"}{...}
{viewerjumpto "Examples" "dea_network##examples"}{...}
{viewerjumpto "Stored results" "dea_network##stored_results"}{...}
{viewerjumpto "References" "dea_network##references"}{...}
{viewerjumpto "Authors" "dea_network##authors"}{...}
{title:Title}

{p2colset 5 21 23 2}{...}
{p2col :{cmd:dea_network} {hline 2}}Network DEA (two-stage model){p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea_network}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}{cmd:,}
{cmdab:int:ermediates(}{it:varlist}{cmd:)}
[{it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea_network}
{ifin}{cmd:,}
{cmdab:in:puts(}{it:varlist}{cmd:)}
{cmdab:int:ermediates(}{it:varlist}{cmd:)}
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
{synopt:{cmdab:in:puts(}{it:varlist}{cmd:)}}Stage 1 input variables (X){p_end}
{synopt:{cmdab:int:ermediates(}{it:varlist}{cmd:)}}intermediate variables linking Stage 1 and Stage 2 (Z){p_end}
{synopt:{cmdab:out:puts(}{it:varlist}{cmd:)}}Stage 2 output variables (Y){p_end}

{syntab:Model}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable (string or numeric){p_end}
{synopt:{cmdab:r:ts(}{it:rtstype}{cmd:)}}returns to scale; default is {cmd:rts(crs)}{p_end}
{synopt:{cmdab:or:ientation(}{it:type}{cmd:)}}orientation; default is {cmd:orientation(input)}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}
{synopt:{cmd:method(}{it:type}{cmd:)}}decomposition method; default is {cmd:method(multiplicative)}{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:prefix}{cmd:)}}prefix for generated efficiency variables; default is {cmd:_net_eff}{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}
{synopt:{cmd:detail}}show component availability and decomposition status{p_end}
{synopt:{cmdab:mat:rix}}return component, status, weight or reference matrices{p_end}

{syntab:Technical}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}positive finite LP residual tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmd:debug}}display debugging information{p_end}
{synoptline}
{p2colreset}{...}

{p 4 6 2}
{it:rtstype}: {cmd:crs} | {cmd:vrs}{p_end}
{p 4 6 2}
{it:type} (orientation): {cmd:input} | {cmd:output}{p_end}
{p 4 6 2}
{it:type} (method): {cmd:multiplicative} | {cmd:additive} | {cmd:independent}{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea_network} computes Network DEA efficiency scores for two-stage production
processes where intermediate products connect the two stages.  The two-stage structure is:

{p 8 17 2}
Stage 1: X {hline 2}> Z  (external inputs to intermediate outputs){break}
Stage 2: Z {hline 2}> Y  (intermediate inputs to final outputs)

{pstd}
Overall efficiency is decomposed into stage-level efficiencies, which allows
identification of which stage is driving overall performance.

{pstd}
Three generated variables are created with the specified prefix:

{p2colset 9 32 34 2}
{p2col:{it:prefix}{cmd:_overall}}overall efficiency score{p_end}
{p2col:{it:prefix}{cmd:_stage1}}Stage 1 efficiency (X {hline 2}> Z){p_end}
{p2col:{it:prefix}{cmd:_stage2}}Stage 2 efficiency (Z {hline 2}> Y){p_end}
{p2colreset}{...}


{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{cmdab:in:puts(}{it:varlist}{cmd:)} specifies the Stage 1 input variables (X).

{phang}
{cmdab:int:ermediates(}{it:varlist}{cmd:)} specifies the intermediate variables (Z)
that are outputs of Stage 1 and inputs to Stage 2.

{phang}
{cmdab:out:puts(}{it:varlist}{cmd:)} specifies the Stage 2 output variables (Y).

{dlgtab:Model}

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies the variable identifying each DMU.
Accepts both string and numeric variables.  If not specified, observation
numbers are used.

{phang}
{cmdab:r:ts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
Allowed values are {cmd:crs} (constant) and {cmd:vrs} (variable).
Default is {cmd:rts(crs)}.

{phang}
{cmdab:or:ientation(}{it:type}{cmd:)} specifies the direction of optimization.
Default is {cmd:orientation(input)}.

{phang}
{cmd:method(}{it:type}{cmd:)} specifies how overall efficiency is decomposed
from stage efficiencies:

{p2colset 9 26 28 2}
{p2col:{cmd:multiplicative}}Kao-Hwang shared intermediate-weight relational
multiplier model under {cmd:rts(crs)} and {cmd:orientation(input)}; default{p_end}
{p2col:{cmd:additive}}overall = 0.5 {it:x} Stage1 + 0.5 {it:x} Stage2 (equal weights){p_end}
{p2col:{cmd:independent}}stages computed independently with no linking constraint;
overall = Stage1 {it:x} Stage2{p_end}
{p2colreset}{...}

{dlgtab:Output}

{phang}
{cmdab:gen:erate(}{it:prefix}{cmd:)} specifies the prefix for generated variables.
Default is {cmd:_net_eff}.

{phang}
{cmd:detail} adds component availability and decomposition status to the score table.

{phang}
{cmd:matrix} returns matrices described below. Rows follow the selected observation
order. Multiplier weights are expressed in the original variable units.


{marker remarks}{...}
{title:Remarks}

{pstd}
Network DEA generalizes the standard DEA model by explicitly modeling
the internal structure of the production process
({help dea_network##FG2000:Färe and Grosskopf 2000}). Whereas standard DEA
treats each DMU as a black box that converts inputs into outputs
({help dea_network##CST2007:Cooper, Seiford, and Tone 2007}),
{cmd:dea_network} models the intermediate flow of products between stages.

{pstd}
The {cmd:multiplicative} method implements the shared intermediate-weight
relational multiplier model of {help dea_network##KH2008:Kao and Hwang (2008)}
for CRS input-oriented two-stage systems.  The same multipliers on intermediate
measures are used as Stage 1 output weights and Stage 2 input weights, and the
reported decomposition satisfies
E{sub:overall} = E{sub:1} {it:x} E{sub:2}.  For VRS or output-oriented calls,
the command retains the legacy independent-stage product and prints a notice.
{cmd:r(formulation)} identifies the actual model; {cmd:r(legacy_fallback)} is 1
for this option combination. It is not a shared-weight VRS or output network model.

{pstd}
The {cmd:additive} method averages two independently optimized stage scores.
It is not the linked additive network DEA model in the literature.
The {cmd:independent} method is equivalent to solving two separate
standard DEA problems and taking the product, with no linking constraint
on the intermediate products.

{pstd}
For the relational model, nonnegative (v,w,u) satisfy wZ <= vX and uY <= wZ
for every selected reference DMU. The command first maximizes uY at vX=1,
then holds that system optimum fixed and maximizes wZ. Thus the secondary
choice favors Stage 1 among alternative optimal decompositions. If that
secondary LP is unavailable, accepted primary weights are retained and the
fallback is reported. Stage ratios are computed from those weights; an
undefined zero denominator remains missing. The system score can remain
available when a component is undefined.

{pstd}
All quantities must be nonnegative. Every selected observation must have at
least one positive variable in each of X, Z, and Y. Missing quantities or an
explicitly missing DMU identifier are excluded before validation. The selected
sample is also the reference sample. An empty sample returns error 2000.
Invalid tolerance or data are rejected before existing results are replaced.

{pstd}
Column and constraint scaling preserve the specified LP. Native solutions must
converge and satisfy the submitted constraints and bounds at the requested
tolerance; relational LPs use an internal convergence cap of 1e-12 and independent
stages use 1e-8, always no larger than the requested tolerance. This is a
numerical acceptance test, not an exact optimality or infeasibility certificate.
Scores are not clipped to [0,1]. Independent output-oriented stages report
1/phi. A failed stage remains missing while an available stage is retained.

{pstd}
The native calculation uses Stata's {cmd:LinearProgram()} class. The exact
verification path described below requires the bundled helpers. Stata 16 or
later is required.


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


{pstd}Setup: bank data where deposits are processed into loans (Stage 1),
which are then converted into interest income (Stage 2).

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input byte bankid deposits staff loans nonloans interest_income profit}{p_end}
{phang2}{cmd:1 120 18  95 28  8.5 12.0}{p_end}
{phang2}{cmd:2 160 24 130 34 11.4 15.5}{p_end}
{phang2}{cmd:3 140 21 118 31 10.2 14.0}{p_end}
{phang2}{cmd:4 190 30 145 40 12.0 16.2}{p_end}
{phang2}{cmd:5 110 16  88 25  7.9 11.1}{p_end}
{phang2}{cmd:6 175 27 138 37 11.8 15.9}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}Basic two-stage network DEA (Form A syntax)

{phang2}{cmd:. dea_network deposits = interest_income, intermediates(loans) dmu(bankid)}{p_end}

{pstd}Form B syntax with explicit options

{phang2}{cmd:. dea_network, inputs(deposits) intermediates(loans) outputs(interest_income) dmu(bankid) rts(vrs) method(multiplicative) replacevars}{p_end}

{pstd}Multiple inputs and intermediates

{phang2}{cmd:. dea_network deposits staff = profit, intermediates(loans nonloans) dmu(bankid) method(additive) replacevars}{p_end}

{pstd}Save results with custom prefix

{phang2}{cmd:. dea_network deposits = interest_income, intermediates(loans) dmu(bankid) generate(bank)}{p_end}
{phang2}{cmd:. list bankid bank_stage1 bank_stage2 bank_overall}{p_end}



{pstd}Equivalent sample restrictions in both syntax forms

{phang2}{cmd:. dea_network deposits = interest_income if inrange(_n,2,5), intermediates(loans) dmu(bankid) generate(qa)}{p_end}
{phang2}{cmd:. dea_network if inrange(_n,2,5), inputs(deposits) outputs(interest_income) intermediates(loans) dmu(bankid) generate(qb)}{p_end}
{phang2}{cmd:. assert qa_overall==qb_overall}{p_end}


{pstd}Documented minimum option abbreviations

{phang2}{cmd:. dea_network, inputs(deposits) outputs(interest_income) intermediates(loans) dmu(bankid) rts(vrs) orientation(output) generate(fullopt)}{p_end}
{phang2}{cmd:. dea_network, in(deposits) out(interest_income) int(loans) dmu(bankid) r(vrs) or(output) gen(shortopt)}{p_end}
{phang2}{cmd:. assert fullopt_overall==shortopt_overall}{p_end}

{pstd}Inspect shared weights and decomposition status (CRS/input){p_end}
{phang2}{cmd:. dea_network deposits staff = profit, intermediates(loans nonloans) dmu(bankid) generate(linked) matrix detail}{p_end}
{phang2}{cmd:. matrix list r(weights_intermediate)}{p_end}
{phang2}{cmd:. matrix list r(decomposition_status)}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}


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

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea_network} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:r(n)}}number of DMUs{p_end}
{synopt:{cmd:r(tolerance)}}requested residual tolerance{p_end}
{synopt:{cmd:r(n_failed_overall)}}missing overall scores{p_end}
{synopt:{cmd:r(n_failed_stage1)}}missing Stage 1 scores{p_end}
{synopt:{cmd:r(n_failed_stage2)}}missing Stage 2 scores{p_end}
{synopt:{cmd:r(n_decomp_fallback)}}relational secondary failures using primary weights{p_end}
{synopt:{cmd:r(legacy_fallback)}}1 if multiplicative uses independent stages, otherwise 0{p_end}
{synopt:{cmd:r(mean_overall)}}mean overall efficiency{p_end}
{synopt:{cmd:r(mean_stage1)}}mean Stage 1 efficiency{p_end}
{synopt:{cmd:r(mean_stage2)}}mean Stage 2 efficiency{p_end}

{syntab:Macros}
{synopt:{cmd:r(method)}}requested method{p_end}
{synopt:{cmd:r(formulation)}}relational_crs, independent_product, or independent_mean{p_end}
{synopt:{cmd:r(rts)}}returns to scale used{p_end}
{synopt:{cmd:r(orientation)}}orientation used{p_end}
{p2colreset}{...}


{pstd}
With {cmd:matrix}, {cmd:r(efficiencies)} has columns overall, stage1, stage2;
{cmd:r(status)} has the same columns with 1 for available and 0 for missing.
For {cmd:relational_crs}, {cmd:r(weights_input)},
{cmd:r(weights_intermediate)}, and {cmd:r(weights_output)} contain (v,w,u).
{cmd:r(decomposition_status)} is 0 for primary failure, 1 for primary-weight
fallback, or 2 for accepted secondary optimization. Undefined component ratios
are also reflected in {cmd:r(status)}. For independent formulations,
{cmd:r(lambda_stage1)} and {cmd:r(lambda_stage2)} contain evaluator rows and
reference-DMU columns. Missing solve rows remain missing.

{marker references}{...}
{title:References}

{marker FG2000}{...}
{phang}
Färe, R., and S. Grosskopf. 2000.
Network DEA.
{it:Socio-Economic Planning Sciences} 34(1): 35--49.
{browse "https://doi.org/10.1016/S0038-0121(99)00012-9":doi:10.1016/S0038-0121(99)00012-9}

{marker KH2008}{...}
{phang}
Kao, C., and S.-N. Hwang. 2008.
Efficiency decomposition in two-stage data envelopment analysis: An
application to non-life insurance companies in Taiwan.
{it:European Journal of Operational Research} 185(1): 418--429.
{browse "https://doi.org/10.1016/j.ejor.2006.11.041":doi:10.1016/j.ejor.2006.11.041}

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
