{smcl}
{* Help revision 2026-11-04: bounded exact duplicate recovery after numerical failure.}{...}
{* Help revision 2026-11-04: original-column reconstruction after failed Phase-1 pivots.}{...}
{* Help revision 2026-11-04: certified price recovery after repeated Phase-2 states.}{...}
{* Help revision 2026-11-03: checked simplex inverse identities and conservative fast acceptance.}{...}
{* Help revision 2026-11-03: explicit convex joint weak-disposal DDF and diagnostics.}{...}
{* Help revision 2026-11-02: SBM zero means, unit scaling and oriented super-SBM.}{...}
{* Help revision 2026-11-02: early allocative price validation and economic-model conventions.}{...}
{* Help revision 2026-11-02: consistent interval realizations and legacy tolerance delivery.}{...}
{* Help revision 2026-11-02: TE numerical diagnostics and small positive score preservation.}{...}
{* Help revision 2026-11-02: congestion/SBM-undesirable missing-result statistics and saving.}{...}
{* Help revision 2026-11-02: TE slack columns identify the expanded input/output varlists.}{...}
{* Help revision 2026-11-02: literal DMU labels in score and reference displays.}{...}
{* Help revision 2026-11-02: explicit IPM Stage-2 radial recovery and saved projection targets.}{...}
{* Help revision 2026-11-01: bounded TE/VRS native retries and per-DMU evidence.}{...}
{* Help revision 2026-11-01: simplex stage routing and stage-specific reference results.}{...}
{* Help revision 2026-11-01: numeric bound lists connected to the legacy mixed-data transformation.}{...}
{* Help revision 2026-11-01: FDH/IDEA result contracts, effective model and saved columns.}{...}
{* Help revision 2026-11-01: observation-aligned VRS matrices and r(vrsobs) source-row mapping.}{...}
{* Help revision 2026-11-01: corrected undesirable DDF, observed directions, sample/domain checks and numerical scaling.}{...}
{* Help revision 2026-11-01: model-specific score counts/ranks and strict super-efficiency count; undesirable LP limitation.}{...}
{* Help revision 2026-11-01: shared single-stage radial bootstrap, actual SE/CI, explicit unsupported options.}{...}
{* Help revision 2026-11-01: validate output names and require replacevars for existing results.}{...}
{* Help revision 2026-11-01: persistent reference weights, final-engine saving, automatic ID protection.}{...}
{* Public option abbreviations and sequential examples validated; legacy using alias documented.}{...}
{* *! version 2.0.73 06nov2026}{...}
{viewerjumpto "Related article" "dea##publication"}{...}
{vieweralsosee "[R] frontier" "help frontier"}{...}
{vieweralsosee "[R] teradial" "help teradial"}{...}
{viewerjumpto "Syntax" "dea##syntax"}{...}
{viewerjumpto "Description" "dea##description"}{...}
{viewerjumpto "Options" "dea##options"}{...}
{viewerjumpto "Remarks" "dea##remarks"}{...}
{viewerjumpto "Examples" "dea##examples"}{...}
{viewerjumpto "Ordinary TE data domain" "dea##te_domain"}{...}
{viewerjumpto "Simplex recovery log" "dea##recoverylog"}{...}
{viewerjumpto "Stored results" "dea##stored_results"}{...}
{viewerjumpto "References" "dea##references"}{...}
{viewerjumpto "Authors" "dea##authors"}{...}
{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{cmd:dea} {hline 2}}Data envelopment analysis{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Form A (legacy syntax)

{p 8 16 2}
{cmd:dea}
{it:inputvars} {cmd:=} {it:outputvars}
{ifin}
[{cmd:using} {it:filename}]
[{cmd:,} {it:options}]


{pstd}
Form B (extended syntax)

{p 8 16 2}
{cmd:dea}
{it:outputvars}
{ifin}{cmd:,}
{cmdab:i:nputs(}{it:varlist}{cmd:)}
[{it:options}]


{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{synopt:{cmdab:mod:el(}{it:modelname}{cmd:)}}DEA model; default is {cmd:model(te)}{p_end}
{synopt:{cmdab:r:ts(}{it:rtstype}{cmd:)}}returns to scale; default is {cmd:rts(crs)}{p_end}
{synopt:{cmdab:or:ientation(}{it:type}{cmd:)}}orientation; default is {cmd:orientation(input)}{p_end}
{synopt:{cmd:ort(}{it:type}{cmd:)}}synonym for {cmd:orientation()}{p_end}
{synopt:{cmd:stage(1|2)}}optimization stages; default is {cmd:stage(2)}{p_end}
{synopt:{cmd:fdh}}shorthand for {cmd:model(fdh)}{p_end}

{syntab:Prices, directions, and undesirable outputs}
{synopt:{cmd:w(}{it:varlist}{cmd:)}}input prices for {cmd:model(cost)}, {cmd:model(profit)}, or {cmd:model(allocative)}{p_end}
{synopt:{cmd:p(}{it:varlist}{cmd:)}}output prices for {cmd:model(revenue)}, {cmd:model(profit)}, or {cmd:model(allocative)}{p_end}
{synopt:{cmd:gx(}{it:varlist}{cmd:)}}input direction vector for {cmd:model(ddf)}{p_end}
{synopt:{cmd:gy(}{it:varlist}{cmd:)}}output direction vector for {cmd:model(ddf)}{p_end}
{synopt:{cmd:ddfbeta(}{it:newvar}{cmd:)}}variable name for DDF beta score{p_end}
{synopt:{cmd:bad(}{it:varlist}{cmd:)}}undesirable (bad) output variables for {cmd:model(undesirable)} and {cmd:model(sbmundesirable)}{p_end}
{synopt:{cmdab:disp:osability(}{it:equality|weak}{cmd:)}}DDF technology; default {cmd:equality}{p_end}
{synopt:{cmdab:weakmat:rix}}return active/idle intensities, directional targets, slacks and diagnostics with {cmd:disposability(weak)}{p_end}

{syntab:Imprecise DEA (IDEA) options}
{synopt:{cmd:xlo(}{it:varlist}{cmd:)}}lower bounds for interval inputs{p_end}
{synopt:{cmd:xhi(}{it:varlist}{cmd:)}}upper bounds for interval inputs{p_end}
{synopt:{cmd:ylo(}{it:varlist}{cmd:)}}lower bounds for interval outputs{p_end}
{synopt:{cmd:yhi(}{it:varlist}{cmd:)}}upper bounds for interval outputs{p_end}
{synopt:{cmd:xord(}{it:varlist}{cmd:)}}ordinal input variables{p_end}
{synopt:{cmd:yord(}{it:varlist}{cmd:)}}ordinal output variables{p_end}
{synopt:{cmd:bounded(}{it:varlist}{cmd:)}}bounded data variables{p_end}
{synopt:{cmd:ordinal(}{it:varlist}{cmd:)}}ordinal data variables{p_end}

{syntab:Output}
{synopt:{cmdab:gen:erate(}{it:newvar}{cmd:)}}name for efficiency variable; default is {cmd:dea_efficiency}{p_end}
{synopt:{cmd:replacevars}}explicitly replace existing result variables; model data and identifiers remain protected{p_end}
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable{p_end}
{synopt:{cmdab:ref:erence}}display reference sets{p_end}
{synopt:{cmdab:saver:ef}}save reference weights as variables{p_end}
{synopt:{cmdab:savew:eights}}save multiplier weights{p_end}
{synopt:{cmdab:sav:ing(}{it:filename}{cmd:)}}save detailed results to dataset{p_end}
{synopt:{cmd:replace}}overwrite existing dataset when using {cmd:saving()}{p_end}

{syntab:Bootstrap (Simar-Wilson)}
{synopt:{cmdab:boo:tstrap}}smoothed bootstrap for radial technical efficiency{p_end}
{synopt:{cmd:reps1(}{it:#}{cmd:)}}bootstrap replications, at least 50; default is {cmd:100}{p_end}
{synopt:{cmd:level(}{it:#}{cmd:)}}confidence level (%); default is {cmd:95}{p_end}

{syntab:Technical}
{synopt:{cmdab:eng:ine(}{it:enginename}{cmd:)}}solver engine; {cmd:simplex} (default), {cmd:ipm}, or {cmd:auto}{p_end}
{synopt:{cmdab:met:hod(}{it:enginename}{cmd:)}}synonym for {cmd:engine()}{p_end}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}positive finite numerical tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmdab:radialt:olerance(}{it:#}{cmd:)}}explicit IPM TE Stage-2 recovery allowance in [0,1); default 0 disables adjustment{p_end}
{synopt:{cmdab:recoveryl:og}}record per-LP simplex recovery context; supported models and codes below{p_end}
{synopt:{cmd:trace}}display iteration details{p_end}
{synopt:{cmd:debug}}display debugging information{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea} computes Data Envelopment Analysis (DEA) efficiency scores for a set
of decision-making units (DMUs) that use multiple inputs to produce multiple
outputs.  Version 2.0 provides sixteen models within a unified interface
and two LP solver engines with model-specific domains and backend limits.

{pstd}
Two syntax forms are supported.  Form A uses the intuitive
{it:inputs} {cmd:=} {it:outputs} formulation from the original {cmd:dea}
command ({help dea##JL2010:Ji and Lee 2010}; {help dea##LJ2009:Lee and Ji 2009}).  Form B specifies outputs as the
main varlist and inputs via the {cmd:inputs()} option.

{pstd}
In addition to the sixteen models in {cmd:dea}, fifteen extension commands
are currently available (see {help dea##alsosee:Also see}).


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{cmd:model(}{it:modelname}{cmd:)} specifies the DEA model to compute.
The available models are:

{p2colset 9 25 27 2}{...}
{p2col :{cmd:te}}Technical efficiency (radial envelopment); this is the
default.  Minimizes (input) or maximizes (output) proportional
adjustment.  Two-stage slack maximization is performed when
{cmd:stage(2)} is specified (the default).{p_end}

{p2col :{cmd:multiplier}}Multiplier (dual) form.  Returns virtual input/output
weights and dual efficiency scores.  Use {cmd:saveweights} to store the
optimal weight vectors and intercept.  Under the command's signed-intercept
convention, the VRS multiplier {it:u0} is positive for DRS, negative for IRS,
and zero for CRS at an unambiguous optimum.{p_end}

{p2col :{cmd:cost}}Cost efficiency.  Requires input prices in {cmd:w()}.
Returns the ratio of minimum cost to observed cost.{p_end}

{p2col :{cmd:revenue}}Revenue efficiency.  Requires output prices in {cmd:p()}.
Returns the ratio of observed revenue to maximum revenue.{p_end}

{p2col :{cmd:profit}}Profit efficiency
({help dea##CCS2007:Cooper, Seiford, and Tone 2007}, Ch. 8).
Requires both input prices in {cmd:w()} and output prices in {cmd:p()}.
Maximizes profit (revenue minus cost) over the technology set subject
to producing at least observed outputs using at most observed inputs.
Returns the ratio of observed profit to maximum profit.  Bounded under
all RTS assumptions because the input constraint X{it:lambda} <= x_o
limits the scale of activity.{p_end}
All price variables and every evaluated DMU's observed profit must be strictly
positive. Zero or negative profit requires a separately specified
directional/Nerlovian measure and is not silently converted to a ratio.{p_end}

{p2col :{cmd:ddf}}Directional distance function
({help dea##CCF1996:Chambers, Chung, and Fare 1996};
{help dea##CRS1998:1998}).
Maximizes simultaneous input contraction and output expansion along
direction vectors specified by {cmd:gx()} and {cmd:gy()}.  Default
direction is the observed values.  The stored beta measures inefficiency:
zero is the directional frontier, smaller values rank better.{p_end}

{p2col :{cmd:sbm}}Slacks-based measure ({help dea##T2001:Tone 2001}).
Non-oriented efficiency incorporates relative slacks on both sides, solved via
{help dea##CC1962:Charnes-Cooper (1962)} transformation. Explicit input or output
orientation evaluates relative slacks on the chosen side. A score of 1 in an
oriented model need not imply zero slacks on the other side.{p_end}

{pstd}
{cmd:model(sbm)} requires nonnegative quantities and at least one positive input
and one positive output in each selected DMU. A zero component is omitted from
its relative-slack average, and the divisor is the number of positive components
on that side for that DMU; all input/output balance constraints remain in force.
For zero outputs this adopts Tone (2001), section 6.2, case 1: the DMU does not
have the function to produce that output. An observed zero caused by unused
production potential requires a separate modeling decision. The command does
not replace zeros with a small positive number or handle negative data by
translation. Entirely zero columns are allowed if each DMU still has a positive
component on each side.

{p2col :{cmd:additive}}Additive model
({help dea##CCD1985:Charnes et al. 1985}).
Maximizes the sum of input and output slacks.  A DMU is efficient if
and only if all optimal slacks are zero.{p_end}

{p2col :{cmd:supereff}}Radial super-efficiency
({help dea##AP1993:Andersen and Petersen 1993}).  Excludes the
evaluated DMU from its own reference set, allowing scores above
unity.  Use {cmd:reference} or {cmd:saveref} to display or store the
reference sets (lambda weights).  In output orientation the stored score is
{it:1/phi}, so larger values rank better in both orientations.  Scores within
tolerance of 1 remain efficient; {cmd:r(N_super)} separately counts scores
strictly greater than {it:1+tolerance}.{p_end}

{p2col :{cmd:supersbm}}Super-efficiency SBM ({help dea##T2002:Tone 2002}).
Super-efficiency extension of the slacks-based measure.  Efficient DMUs
(ordinary SBM score >= 1-tolerance) receive a super-efficiency score via Tone's
projected-quantity formulation with self-exclusion. DMUs below that threshold
retain their ordinary SBM score in the requested orientation. All selected
inputs and outputs must be strictly positive. No automatic zero substitution
is performed.{p_end}

{pstd}
The default super-SBM is non-oriented. Input orientation fixes the projected
outputs at the observed outputs and minimizes the average projected-input
ratio. Output orientation fixes the projected inputs at the observed inputs,
maximizes the average projected-output ratio, and stores its reciprocal.
Both orientations therefore use larger scores for better performance.
Oriented super-SBM can be infeasible under non-CRS technologies, even with
strictly positive data; an empty reference set also has no score. Such results
remain missing and are included in {cmd:r(n_failed)}. That count records
unavailable results, without by itself distinguishing mathematical
infeasibility from numerical solver failure.

{p2col :{cmd:allocative}}Allocative efficiency.  Decomposes economic efficiency
into technical and allocative components following
{help dea##F1957:Farrell (1957)} and
{help dea##FL1978:Fare and Lovell (1978)}.  Requires {cmd:w()} for input-oriented
or {cmd:p()} for output-oriented analysis.{p_end}

{pstd}
For {cmd:model(allocative)}, input orientation requires one {cmd:w()} variable
per input, and output orientation requires one {cmd:p()} variable per output.
Every supplied price list must match the corresponding quantity list in length
and order. From version 2.0.57, a missing required list or incorrect length is
rejected before any result variable is created or replaced. Earlier versions
could fail during computation or accept a short list and return missing scores.{p_end}

{pstd}
Allocative efficiency is {it:CE/TE} in input orientation and {it:RE/TE} in
output orientation, using this command's stored technical score in (0,1].
When comparing with software that reports output expansion {it:phi >= 1},
the output decomposition is {it:AE = RE*phi}. The profit model above restricts
improvements to {it:X*lambda <= x_o} and {it:Y*lambda >= y_o}; it is not the
unrestricted profit maximization returned by R {cmd:Benchmarking::profit.opt}.
Different feasible sets can produce different profit scores.{p_end}

{p2col :{cmd:fdh}}Free disposal hull
({help dea##DPT2006:Deprins, Simar, and Tulkens 2006}).
Replaces convexity with integer constraints on lambda, producing a
staircase frontier.  Computed by enumeration without LP.{p_end}

{p2col :{cmd:imprecise}}Imprecise DEA (IDEA)
({help dea##CPY1999:Cooper, Park, and Yu 1999}).  Handles interval,
ordinal, and bounded data through separate interval and legacy
transformation paths. See {it:Imprecise DEA options} below for their
output definitions and current limitations.{p_end}

{p2col :{cmd:undesirable}}Directional distance function with undesirable
(bad) outputs.  Extends {cmd:model(ddf)} by incorporating bad outputs
that are contracted while good outputs are expanded.  Requires
{cmd:bad(}{it:varlist}{cmd:)}.  Direction controlled by {cmd:direction()}:
{cmd:yb} (default, expand good outputs and contract bad outputs),
{cmd:xyb} (contract inputs and bad outputs while expanding good outputs),
or {cmd:x} (contract inputs only).{p_end}

{p2col :{cmd:congestion}}Legacy radial equality-ratio diagnostic (default), or
explicit {cmd:congmethod(cooper)} slack decomposition.
For background on congestion methods, see {help dea##CDSZ2004:Cooper et al. (2004)}.
Retains a legacy equality-restricted radial diagnostic. The ratio alone
does not identify input congestion; see {helpb dea_congestion}. Reports the ratio
({it:theta_VRS / theta_C}; below 1 indicates an equality-restriction gap), plus CRS and
VRS efficiency and scale efficiency.  Additional variables are created
with suffixes {cmd:_crs}, {cmd:_vrs}, {cmd:_cong_eff}, {cmd:_congestion},
and {cmd:_se}. With {cmd:congmethod(cooper)}, the primary result instead is mean
relative input congestion (0 means none); suffixes are {cmd:_phi}, {cmd:_eff},
{cmd:_total}, {cmd:_technical} and {cmd:_status}. Options {cmd:congweight(raw|normalized)},
{cmd:congmatrix} and {cmd:congranges} control slack objectives and optional
conditional optimal-face diagnostics. Simplex uses bounded exact certification
for strictly positive data with {cmd:2*N+inputs+outputs+3<=64}; ipm retains the
native numerical path. Cooper auto is rejected. Simplex-only options
{cmd:conganchor(A)}, {cmd:congpolicy(stage1|stage2)}, {cmd:congreference(M)} and
{cmd:congselectlimit(16|128)} control explicit reference selection; defaults keep
the baseline point and 16 exact basis exchanges. References follow the selected
sample order, with no automatic ID matching. Both engines reject unrepresentable
objective scaling. The wrapper's default engine is ipm, while dea defaults to
simplex. See {cmd:r(actual_engine)} and {cmd:r(actual_backend)}. See {help dea_congestion##cooper:Cooper method}
for equations, strict-positive domain, phase failures and stored results.
Always uses VRS; the {cmd:rts()} option is ignored.{p_end}

{p2col :{cmd:sbmundesirable}}SBM with undesirable outputs
({help dea##T2004:Tone 2004}).  Non-radial efficiency incorporating
slacks for inputs, good outputs, and bad outputs simultaneously.
Bad-output reduction uses independent nonnegative slacks; B lambda+s^b=b0
is equivalent to B lambda<=b0, not a joint-abatement restriction.
Requires {cmd:bad(}{it:varlist}{cmd:)}.
See also {helpb dea_sbmundesirable}.{p_end}
{p2colreset}{...}

{pstd}
From version 2.0.53, legacy congestion and SBM-undesirable failures remain missing.
The following ratio-specific counts apply to {cmd:congmethod(legacy)}.
Their {cmd:r(N_eff)} counts require a finite main score within tolerance of 1;
means and extrema use nonmissing main scores. For SBM-undesirable,
{cmd:r(n_failed)} counts missing scores. For congestion it counts the union of
rows with unavailable CRS, VRS, congestion-model or derived ratio results;
{cmd:r(n_failed_primary)} counts only missing main measures. Separate component
counts are described in {helpb dea_congestion}. A failed component does not erase
other available components. A ratio with a denominator no greater than
{cmd:tolerance()} remains missing and is never assigned 1.

{pstd}
For these models, {cmd:saving()} retains failed rows in the selected sample.
Congestion files also include all five auxiliary columns; SBM-undesirable
files include the {cmd:bad()} variables. These additions apply to the core
command; the two convenience wrappers do not accept {cmd:saving()}.

{pstd}
{bf:Correction in version 2.0.34:} {cmd:undesirable} now uses the correct
good-output expansion constraint in both engine routes.  Recompute results
obtained with {cmd:direction(yb)} or {cmd:direction(xyb)} in earlier versions;
even identical positive CRS DMUs could previously receive beta=1 instead of 0.
This correction is separate from {cmd:model(sbmundesirable)}.{p_end}


{marker te_domain}{...}
{phang}
{bf:Quantity domain for ordinary TE}

{pstd}
From {cmd:dea} 2.0.63, ordinary {cmd:model(te)} requires nonnegative inputs and
outputs and at least one strictly positive input and one strictly positive
output in every selected DMU. Individual zero components and entire zero
columns are permitted if this per-DMU condition holds. The check uses original
values, without a tolerance band, translation, clipping or a small positive
replacement. It applies to both syntaxes, all RTS choices, both orientations,
both stages, and {cmd:engine(simplex|ipm|auto)}. It also applies to the ordinary
TE sample before bootstrap and with {cmd:radialtolerance()}.

{pstd}
Selection by {cmd:if}/{cmd:in} and exclusion of observations with missing
required quantities precede this validation. Invalid values in an excluded
observation do not stop the command. An empty complete sample retains error
{cmd:r(2000)}. A selected negative quantity or an all-zero input or output
vector causes {cmd:r(459)} before existing result variables, reference variables,
automatic IDs or saved files are replaced. With {cmd:recoverylog}, the rejected
command's recording is discarded. The existing DMU-label selection convention
is unchanged; a missing label alone is not a missing input/output quantity.

{pstd}
Earlier ordinary TE calls could accept these invalid quantity vectors and
return scores. Such scores require review under the intended DEA data domain.
This change does not alter the LP or its handling of valid data, small positive
scores, numerical failures or alternative optima. FDH, IDEA and other models
retain their own domain rules; in particular, natural VRS super-efficiency
infeasibility is not reclassified by this ordinary-TE check.

{phang}
{cmd:rts(}{it:rtstype}{cmd:)} specifies the returns-to-scale assumption.
{cmd:crs} (constant, the default; synonym {cmd:ccr}) imposes no constraint on
lambda.  {cmd:vrs} (variable; synonym {cmd:bcc}) adds the convexity constraint
sum(lambda) = 1.  {cmd:nirs} (non-increasing) adds sum(lambda) <= 1;
the legacy synonym {cmd:drs} is retained.  {cmd:ndrs} (non-decreasing) adds
sum(lambda) >= 1; the legacy synonym {cmd:irs} is retained.  NIRS and NDRS
name technology restrictions.  They should not be confused with the DRS and
IRS classifications reported by the VRS frontier analysis.
When {cmd:rts(vrs)} is specified with {cmd:model(te)}, scale efficiency and
returns-to-scale classification are automatically computed following
{help dea##FGL1993:Fare, Grosskopf, and Lovell (1993)}.  This VRS
Frontier Analysis table is available under both the simplex and IPM
engines.

{phang}
{cmd:orientation(}{it:type}{cmd:)} or {cmd:ort(}{it:type}{cmd:)} specifies
{cmd:input} (the default) for input-oriented models that minimize input usage,
or {cmd:output} for output-oriented models that maximize output production.
Abbreviations {cmd:i} and {cmd:o} are accepted.
Not applicable to {cmd:model(additive)}, {cmd:model(cost)},
{cmd:model(revenue)}, or {cmd:model(profit)}.
For {cmd:model(sbm)} and {cmd:model(supersbm)}, omitting the option requests the non-oriented model;
{cmd:orientation(input)} and {cmd:orientation(output)} request the corresponding
oriented models.

{pstd}
From {cmd:dea} 2.0.58 (simplex 2.0.14, IPM 2.0.16), super-SBM actually applies
the explicit orientation; earlier versions returned the non-oriented score
even when an orientation was specified. Recompute those results when orientation
matters. Ordinary SBM now reduces the averaging divisor when omitting zero
terms. Both models rescale columns internally to limit dependence on measurement
units, preserving physical constraints and the original units of any slacks.
Earlier extreme-unit results could be missing or incorrect. The result display
and {cmd:r(orientation)} now report {cmd:non-oriented} when the SBM orientation
option is omitted.

{phang}
{cmd:stage(1|2)} selects radial-only ({cmd:1}) or two-stage ({cmd:2},
the default) optimization for technical efficiency. Stage 1 optimizes the
radial score and, when requested, returns that LP's reference weights.
It does not maximize slacks or return public slack matrices. Stage 2 holds
the radial score fixed and maximizes slacks. A successful Stage 2 selects a
slack-maximizing projection; the projection and reference weights need not
be unique. Both engines maximize the sum of slacks in their normalized
input/output columns and reports slacks in the original units.

{pstd}
From {cmd:dea} 2.0.41 with simplex 2.0.8, {cmd:stage(1)} skips the second LP.
Earlier simplex versions also ran Stage 2 for {cmd:stage(1)}. Radial scores
are unchanged on successful solves, but Stage-1 reference weights can now
differ. Recompute earlier {cmd:stage(1) reference} or {cmd:saveref} results
when the distinction matters. Auxiliary CRS/VRS/NIRS fits for the VRS
decomposition use only radial scores and skip slack maximization.

{pstd}
From {cmd:dea} 2.0.52, both TE engines retain an available radial score
but leave reference weights and slacks missing when the requested Stage 2
is unavailable. Earlier simplex versions could return Stage-1 references
and slacks after Stage-2 failure; those were not certified slack maxima.
If the simplex Stage-1 primal fails but the dual succeeds, {cmd:stage(1)}
retains the score with missing references; {cmd:stage(2)} attempts to recover
references and slacks in its second LP. Other models retain their own
stage conventions described below.

{pstd}
From {cmd:dea} 2.0.46 with IPM 2.0.8, IPM technical-efficiency Stage 2 uses
the same RTS technology as Stage 1: no lambda-sum restriction for CRS,
sum(lambda)=1 for VRS, <=1 for DRS/NIRS, and >=1 for IRS/NDRS.
On success, {cmd:r(lambda)} and the two slack matrices describe the Stage-2
projection. Earlier IPM versions returned Stage-1 weights and did not return
the slack matrices; recompute those reference results when Stage 2 matters.

{pstd}
For {cmd:model(te)} with either engine, {cmd:r(stage_status)} reports two
columns for each selected observation. {cmd:radial} is 0 for an available
score (including dual recovery) and 1 for a missing score. {cmd:stage2} is
0 for success, positive for failure, or missing if not attempted.
Simplex codes are 1 (infeasible), 2 (unbounded), and 3 (numerical failure or
nonconvergence). IPM uses native error codes, including 3 for nonconvergence
or exhausted retries, 90 for residual rejection after unsuccessful recovery,
and 91 for contradictory native reports. These codes describe solver outcomes;
a numerical failure does not prove mathematical infeasibility.
Stage 2 is skipped for {cmd:stage(1)}, when both radial formulations fail,
or in IPM when its primal remains unavailable after retries.

{pstd}
{cmd:r(solution_source)} has two columns: {cmd:radial} is 0 for unavailable,
1 for the primal score, and 2 for dual recovery; {cmd:reference} is 0 for
unavailable, 1 for Stage-1 references, and 2 for Stage-2 references.
Provenance is reported even without {cmd:reference} or {cmd:saveref}.
A finite score alone does not imply an available projection.

{pstd}
For both TE engines, {cmd:r(n_failed)} counts a row once if its radial score
is missing or its requested Stage 2 is unavailable. {cmd:r(n_radial_failed)}
counts missing scores; {cmd:r(n_stage2_failed)} counts unavailable requested
Stage-2 projections, including ones not attempted after a primal failure.
It is zero for {cmd:stage(1)}. Score summaries and ranks still use available
radial scores. The status and projection returns survive VRS decomposition,
bootstrap, and an automatic retry to IPM.
The existing {cmd:engine(auto)} choice and retry policy is unchanged: it may
select IPM for heterogeneous scales or retry missing simplex radial scores.
Stage-2 failure alone does not trigger an engine switch.

{pstd}
From {cmd:dea} 2.0.54, {cmd:r(numerical_diagnostics)} reports a 23-column
matrix for ordinary TE with either engine, even without {cmd:reference}.
Rows follow {cmd:r(efficiency)}. It audits the returned Stage-2 vector;
it does not replace the solver status or trigger a new optimization.
The columns are:

{p2colset 5 34 36 2}{...}
{p2col:{cmd:radial_initial}}factor used to fix the original Stage-2 LP{p_end}
{p2col:{cmd:radial_stage2}}factor of the returned projection; missing if unavailable{p_end}
{p2col:{cmd:radial_adjustment}}returned factor minus initial factor{p_end}
{p2col:{cmd:stage2_status}}same final code as {cmd:r(stage_status)}{p_end}
{p2col:{cmd:checked}}1 if the returned projection could be checked, otherwise 0{p_end}
{p2col:{cmd:within_tolerance}}1 if the maximum scaled residual is <= {cmd:tolerance()},
0 if above it, missing if not checked{p_end}
{p2col:{cmd:input_abs output_abs}}maximum absolute equality residual in original input/output units{p_end}
{p2col:{cmd:rts_abs}}lambda-sum restriction violation; zero for CRS{p_end}
{p2col:{cmd:lambda_sum lambda_gap}}sum(lambda), and sum(lambda)-1 (signed, also shown for CRS){p_end}
{p2col:{cmd:lambda_bound}}maximum negative-lambda violation{p_end}
{p2col:{cmd:slack_in_bound slack_out_bound}}maximum negative-slack violations in original units{p_end}
{p2col:{cmd:input_scaled output_scaled}}maximum scaled input/output residual using column-normalized data and slacks{p_end}
{p2col:{cmd:rts_scaled}}scaled violation of the original RTS equality/inequality{p_end}
{p2col:{cmd:lambda_scaled}}scaled lambda-bound violation{p_end}
{p2col:{cmd:slack_in_scaled slack_out_scaled}}scaled bounds on normalized slacks{p_end}
{p2col:{cmd:max_scaled}}maximum of the preceding six scaled measures{p_end}
{p2col:{cmd:input_normalized output_normalized}}maximum absolute input/output residual in column-normalized units{p_end}
{p2colreset}{...}

{pstd}
Factors are theta for input orientation and phi for output orientation;
the output efficiency score is 1/phi. Original and normalized absolute
residuals are different quantities. For each normalized constraint row,
the scaled denominator is max(abs(b), sum(abs(a_j)*max(1,abs(z_j))));
zero/zero contributes zero. Nonnegativity bounds use max(1,abs(z_j)).
RTS diagnostics use the original lambda-sum restriction, without an
internal simplex convexity-slack variable. Absolute maxima can mix variable
units within each input/output group and should be interpreted accordingly.

{pstd}
{cmd:stage(1)} returns the initial factor and missing Stage-2 diagnostics;
the diagnostic audit itself performs no extra optimization. The simplex
stage-link guard below can retry Stage 1 without running Stage 2. Missing or arithmetically uncheckable
projections have {cmd:checked}=0 and missing residuals. Counts are
{cmd:r(n_projection_checked)}, {cmd:r(n_projection_outside_tol)}, and
{cmd:r(n_projection_unchecked)} (the last counts available Stage-2 results
that could not be checked). A warning identifies the latter two cases.
These are diagnostic counts; they do not redefine {cmd:r(n_failed)}.
The matrix is returned in memory; {cmd:saving()} does not add its columns.

{pstd}
An original fixed-radial LP can be exactly infeasible while its returned
vector passes the requested numerical criterion. Neither a zero displayed
residual nor {cmd:within_tolerance}=1 proves exact feasibility or optimality.
Already accepted projections are not reoptimized; {cmd:radialtolerance()}
keeps its default zero and its existing explicit-IPM, failed-projection scope.
The diagnostics show the changed target when that recovery succeeds.
Native IPM statuses and residual checks remain separately visible in
{cmd:r(ipm_trace)} and, for changed targets, {cmd:r(projection_trace)}.

{pstd}
Simplex no longer rounds a positive radial efficiency below 1e-8 to zero.
The same score is retained after a requested Stage-2 failure. Reference/slack
cleanup is retained only if the cleaned vector still passes the original LP's
scaled residual test. Otherwise the original solver-accepted components are
returned, without a new optimization or radial adjustment. This avoids erasing
small weights needed to satisfy the constraints. Tiny negative components may
remain within numerical tolerance; bound diagnostics report their size.
The returned-vector audit includes this extraction step. Solver pivot/convergence
tolerances are unchanged; very small scores are not guaranteed exact by this correction.

{pstd}
For TE, {cmd:saving()} also writes {cmd:_dea_radial_status},
{cmd:_dea_stage2_status}, {cmd:_dea_score_source}, and {cmd:_dea_ref_stage}
in selected-observation order, matching these two matrices. These four columns
are added only to the saved dataset. A name conflict with a retained input,
output, ID, or score is rejected before computation; unrelated source columns
are still omitted by {cmd:saving()}.

{pstd}
From {cmd:dea} 2.0.47 / IPM 2.0.9, {cmd:tolerance()} is passed to native
TE radial, dual, Stage-2 and auxiliary VRS solves. Successful candidates must
also pass the shared scaled primal-residual check described in {help lp}.
TE uses column-normalized input/output data; Stage-2 slack variables are in
those normalized units during checking and are returned in source units.
IPM checked slacks are no longer rounded to zero. The simplex core applies the
same criterion to the real-variable standard form after Phase I and Phase II.

{pstd}
A rejected IPM radial primal candidate may recover a score from the checked
dual solution, without primal reference weights. Output dual recovery returns
{it:1/phi}, correcting earlier versions that returned {it:phi} as efficiency.
The auxiliary VRS analysis uses the same data normalization and numerical
settings as the main TE path. Stage 2 keeps the computed radial score fixed;
no implicit tolerance band or score adjustment is added.

{pstd}
Residual acceptance certifies only numerical primal feasibility within the
requested criterion. It is not an independent optimality or infeasibility
certificate. Native boundary failures remain possible; convergence settings
can change both last digits and availability of projections. Other native IPM
model kernels retain their individual numerical policies. General optimality
and infeasibility certification remains separate work.

{pstd}
From {cmd:dea} 2.0.48 / IPM 2.0.10, the native TE primal, Stage-2, dual-recovery
and auxiliary CRS/VRS/DRS LPs use bounded native retries after a failed or
rejected initial solve. The three variants are the shared policy in {help lp}:
replace equalities by both inequality directions; scale rows/objective in that
representation; then tighten convergence to min(requested tolerance,1e-12) in
the original representation. Inapplicable variants are skipped. There are at
most three extra native solves per LP. Every candidate must pass the original
LP's residual check at the requested tolerance. A successful initial solution
is returned unchanged. These formulation retries do not move the radial target.

{pstd}
An infeasible report may recover through a checked native solution. An unbounded
report contradicted by an infeasible report or an accepted finite candidate
produces final status 91. Any unbounded report in the primal history prevents
using dual fallback to conceal that unresolved report. Explicit {cmd:engine(ipm)}
uses native IPM throughout. A checked dual recovery supplies only the radial
score; Stage-2 projections still require a recovered primal solution. The
separate top-level {cmd:engine(auto)} policy remains as documented below.

{pstd}
{cmd:r(ipm_trace)} has one row per selected DMU and 75 columns: three 25-column
blocks with prefixes {cmd:p} (radial primal), {cmd:s} (Stage 2), and {cmd:d}
(radial dual). Each block starts with {it:prefix}{cmd:_status}, its final LP
status, followed by four six-column attempt groups. Names such as
{cmd:p0_native} and {cmd:s1_resid} use variant 0 for the initial solve and
1/2/3 for the retries above. Each group contains {cmd:native}, {cmd:checked},
{cmd:tol}, {cmd:iter}, {cmd:resid}, and {cmd:obj}. The raw native code is retained;
checked code 0 means the candidate passed the primal check, and 90 means it was
rejected. A checked code 0 may still lead to final status 91 when reports conflict.
Unused attempts/phases are missing. Stage(1) leaves the whole {cmd:s} block missing.

{pstd}
The trace records actual convergence settings, iteration counts, maximum scaled
residuals, and LP minimization objectives. Objectives are in the normalized
LP's coordinates and sign convention; they are not always efficiency scores.
For example, the output primal objective is -phi and Stage 2 minimizes the
negative sum of normalized slacks. {cmd:r(vrs_ipm_trace)} has 150 columns for
six analogous blocks: {cmd:crs_p}, {cmd:crs_d}, {cmd:vrs_p}, {cmd:vrs_d},
{cmd:drs_p}, {cmd:drs_d}. Rows follow {cmd:r(vrsobs)}. It is returned only for
an IPM auxiliary VRS decomposition. These matrices survive saving, bootstrap,
and top-level auto fallback; successful simplex calls do not return them.

{pstd}
{cmd:r(n_ipm_retries)}, {cmd:r(n_ipm_recovered)}, and {cmd:r(n_ipm_conflicts)}
count extra native solves, LP phases recovered after initial failure, and
conflicted LP phases in the main TE calculation. These are LP counts, not DMU
counts. The corresponding auxiliary counts use {cmd:n_vrs_ipm_} as prefix.
They do not replace {cmd:r(n_failed)} or the two DMU-level failure counts.
The common TE function also supplies radial scores for the multiplier and
allocative wrappers, but their other optimization kernels and public returns
retain their existing policies.

{phang}
{cmd:radialtolerance(}{it:#}{cmd:)} permits a separate recovery attempt when an
IPM TE Stage-2 projection remains unavailable. The value must be finite and
in [0,1); the default 0 makes no adjustment and preserves the previous results.
A positive value requires explicit {cmd:engine(ipm)} (or {cmd:method(ipm)}),
{cmd:model(te)} and {cmd:stage(2)}, without FDH or imprecise-data options.
It is rejected with {cmd:engine(simplex)}, {cmd:engine(auto)}, Stage 1, and other
models, including VRS super-efficiency. Selected complete observations must
have nonnegative quantities and at least one positive input and one positive
output per DMU. This option does not change bootstrap's original radial scores.

{pstd}
Let {it:r} be the original radial factor: theta for input orientation, phi=1/score
for output orientation. After ordinary fixed-radial Stage 2 and its native
retries fail, the new target is {it:r+d*(1+abs(r))} for input orientation or
{it:max(1,r-d*(1+abs(r)))} for output orientation, where {it:d} is the declared
allowance. This is one conservative endpoint attempt, not a search for the
smallest adjustment. A nonfinite, unchanged or incorrectly directed endpoint
is not attempted. An unavailable primal or any unresolved unbounded/conflicting
Stage-2 report also prevents adjustment. The changed LP uses the same RTS,
normalization, native retry policy and requested residual tolerance.

{pstd}
A successful recovery supplies reference weights and slacks at the reported
new radial factor. The generated efficiency variable, {cmd:r(efficiency)},
score summaries/ranks and auxiliary VRS scores retain their original values.
Already accepted fixed-radial projections remain unchanged. If recovery fails,
references and slacks remain missing. The allowance does not guarantee a
solution or establish exact feasibility/optimality. In particular, an already
accepted approximate projection is not reclassified by an exact mathematical
test. The structural infeasibility of VRS super-efficiency is a separate issue.

{pstd}
With a positive allowance, {cmd:r(projection)} has eight columns:
{cmd:radial_initial}, {cmd:radial_stage2}, {cmd:radial_adjustment},
{cmd:radial_target}, {cmd:status_initial}, {cmd:status_retry},
{cmd:status_final}, and {cmd:adjusted}. Rows follow {cmd:r(efficiency)}.
The second and third columns are missing when no projection is available;
otherwise the third is the second minus the first. The attempted target and
retry status are missing if no adjusted LP was attempted. {cmd:adjusted}=1
only for a successful changed-target solve. Input adjustments are positive;
output adjustments in phi units are negative. {cmd:r(stage_status)} describes
the final available projection, including this recovery.

{pstd}
The original fixed-radial history remains in {cmd:r(ipm_trace)}. Its {cmd:s_status}
may differ from final Stage-2 status when adjustment is enabled.
{cmd:r(projection_trace)} is a separate {it:N} x 25 matrix for the changed LP,
with prefix {cmd:a} and the same final-status/attempt fields as each IPM trace
block. The existing {cmd:n_ipm_} counters exclude this separate recovery.
{cmd:r(n_projection_attempted)} counts changed-target LPs,
{cmd:r(n_projection_adjusted)} counts successful changed projections, and
{cmd:r(n_projection_solves)} counts all native calls for changed LPs, including
formulation retries. Original and changed histories survive saving and bootstrap.

{pstd}
Three double variables are also generated: {cmd:_dea_radial},
{cmd:_dea_radial2}, and {cmd:_dea_radial_adj}, corresponding to the first three
projection columns. They are missing outside the selected sample and included
in {cmd:saving()} alongside the original efficiency score and any saved
references. They follow the same {cmd:replacevars} and protected-input rules as
other result variables. No additional variables or projection matrices are
created when the allowance is zero. Output-oriented factors are phi, while
the efficiency score remains 1/phi; use the recorded projection factor when
checking the reference/slack identities.

{phang}
{cmd:fdh} is shorthand for {cmd:model(fdh)}.  When specified, the
returned RTS label is VRS for both spellings. Default/CRS and VRS labels
are accepted; NIRS/NDRS (including DRS/IRS aliases) are rejected. The
nonconvex dominance comparison does not use an LP solver or a convex VRS
constraint. {cmd:engine()} and {cmd:stage()} do not change FDH scores.

{dlgtab:Prices and directions}

{phang}
{cmd:w(}{it:varlist}{cmd:)} specifies input price variables for
{cmd:model(cost)}, {cmd:model(profit)}, and input-oriented {cmd:model(allocative)}.
The number of price variables must equal the number of inputs.

{phang}
{cmd:p(}{it:varlist}{cmd:)} specifies output price variables for
{cmd:model(revenue)}, {cmd:model(profit)}, and output-oriented {cmd:model(allocative)}.
The number of price variables must equal the number of outputs.

{pstd}
For {cmd:model(ddf)}, observations missing any explicitly supplied {cmd:gx()}
or {cmd:gy()} value are excluded from both evaluation and reference samples.
For {cmd:model(undesirable)} and {cmd:model(sbmundesirable)}, the same rule
applies to {cmd:bad()}. Ordinary input/output missing values are also excluded.
These rules apply after {cmd:if}/{cmd:in}; {cmd:r(N)} and saved results reflect
the selected complete cases. Result variables are missing outside that sample.
An empty sample stops before result replacement or file writing.

{phang}
{cmd:gx(}{it:varlist}{cmd:)} specifies the input direction vector for
{cmd:model(ddf)}.  The default uses observed input values (proportional
direction).

{phang}
{cmd:gy(}{it:varlist}{cmd:)} specifies the output direction vector for
{cmd:model(ddf)}.  The default uses observed output values (proportional
direction).

{phang}
{cmd:ddfbeta(}{it:newvar}{cmd:)} specifies the variable name for the DDF
inefficiency score beta.  It is allowed only with {cmd:model(ddf)} and is an
alias for {cmd:generate()}; do not specify both options.

{phang}
{cmd:bad(}{it:varlist}{cmd:)} specifies the undesirable (bad) output variables.
Required for {cmd:model(undesirable)} and {cmd:model(sbmundesirable)}.
The technologies differ: DDF uses a bad-output equality without a slack;
SBM uses independent bad-output reduction slacks. Neither alone imposes
null-jointness or estimates an abatement cost.

{phang}
{cmd:direction(}{it:type}{cmd:)} specifies the contraction/expansion direction
for {cmd:model(undesirable)}.  Options are:

{p2colset 9 16 18 2}
{p2col:{cmd:yb}}expand good outputs and contract bad outputs; hold inputs fixed (default){p_end}
{p2col:{cmd:xyb}}contract inputs and bad outputs simultaneously while expanding good outputs{p_end}
{p2col:{cmd:x}}contract inputs only{p_end}
{p2colreset}{...}

{pstd}
For {cmd:model(undesirable)} with default {cmd:disposability(equality)},
the LP maximizes beta subject to
{it:X'lambda + beta*gx <= x0},
{it:Y'lambda >= y0 + beta*gy}, and
{it:B'lambda + beta*gb = b0}, with nonnegative lambda and beta.
The observed directions (gx,gy,gb) are (0,y0,b0) for {cmd:yb},
(x0,y0,b0) for {cmd:xyb}, and (x0,0,0) for {cmd:x}.
{cmd:direction()} determines these movements; {cmd:orientation()} and
{cmd:stage()} do not change this single LP.
The chosen RTS restriction on sum(lambda) is retained: none for CRS,
1 for VRS, at most 1 for NIRS, and at least 1 for NDRS.

{pstd}
The bad-equality DDF technology with VRS or NDRS is not generally closed
under proportional joint output reduction: sum(lambda)=1 or >=1 can prevent
scaling all activities down. Thus this is the documented observed-direction
equality model, not a general convex weak-disposal VRS model with
activity-specific abatement factors. Under CRS/NIRS, nonnegative inputs
permit proportional activity reduction. No additional abatement variables
are introduced. With all reference bad-output sums positive, B'lambda=0
forces lambda=0; zero-bad reference producers can invalidate that
null-jointness implication. The program does not impose it as an extra rule.

{pstd}
For the default paths of these specialized models, {cmd:r(formulation)} identifies the actual
specification and {cmd:r(lp_backend)} is {cmd:stata_linearprogram}.
{cmd:r(engine)} retains the shared engine selection; it does not imply that
this model used the package simplex solver. Effective orientations are
{cmd:directional} for undesirable DDF, {cmd:non-oriented} for SBM-undesirable,
and {cmd:input} for the legacy congestion ratio. For the default path, the parsed request remains
in {cmd:r(orientation_requested)}. For DDF, {cmd:r(direction)} is the actual
observed-direction choice. SBM-undesirable requires strictly positive
quantities and supports all four RTS restrictions; see its specialized help.

{pstd}
For the default {cmd:disposability(equality)} path,
observations missing any input, good output, or {cmd:bad()} variable are
excluded after applying {cmd:if}/{cmd:in}.  Selected values must be nonnegative,
and every selected DMU must have at least one strictly positive component in
its chosen direction.  Individual zero values and zero columns are allowed;
an all-zero direction is rejected before any existing result is replaced.
These checks do not impose null-jointness or guarantee a bounded LP for every
nonnegative dataset.  An unsolved or unbounded DMU has a missing score and is
included in {cmd:r(n_failed)}.

{pstd}
Both engine routes use Stata's {cmd:LinearProgram()} for this model.
Variables and their observed directions are scaled together by column maxima,
and native convergence uses {it:min(tolerance,1e-12)}. Returned solutions
must satisfy the common scaled primal and bound checks at the requested
{cmd:tolerance()}; rejected solutions remain missing. This option also controls
score classification and ranking. Independent LP comparisons and unit-change
tests are provided in {cmd:tests/test_undesirable_ddf.py} and
{cmd:tests/test_nonte_tolerance.py}.

{phang}
{cmd:disposability(weak)} selects a convex technology with separate nonnegative
active intensities z and idle/abatement intensities w. It maximizes beta >= 0:

{p 8 8 2}
X'(z+w) + beta*gx <= x0; Y'z >= y0 + beta*gy;
B'z + beta*gb = b0.

{pstd}
The RTS restriction applies to the {it:total} intensity z+w: none for CRS,
sum(z+w)=1 for VRS, <=1 for NIRS, and >=1 for NDRS. The VRS technology follows
the activity-specific abatement construction of Kuosmanen (2005);
the NIRS/NDRS variants here explicitly extend the restriction on total intensity.
This is an input-commitment technology without a materials-balance constraint.
It is not a single common abatement factor or a by-production model.
The observed directions and stored beta use the conventions above.
Zero beta is directionally efficient; a larger beta permits greater improvement.
A zero distance does not rule out non-directional slacks.

{pstd}
Joint weak disposability follows directly: for 0 <= a <= 1, replace z by a*z
and w by w+(1-a)*z. Total intensity and input use stay fixed, while both outputs
contract together. For nonnegative inputs the CRS/NIRS technology is equivalent
to the default equality technology; for VRS/NDRS it contains that technology,
so the mathematical optimum beta cannot decrease. Numerical solver failures
can differ. Null-jointness is conditional on the reference data: every producer
of positive good output must have a positive total bad output. The command
allows zero-bad producers and therefore does not guarantee null-jointness.

{pstd}
The weak option requires nonnegative data with at least one strictly positive
input and one strictly positive good output in every selected DMU. Individual
zeros and all-zero columns are allowed. Missing input/good/bad/ID values are
excluded after {cmd:if}/{cmd:in}; sample order determines matrix rows and columns.
These conditions ensure a bounded self-referencing LP in exact arithmetic.
Floating-point availability is checked separately; unavailable results remain
missing rather than being called infeasible or replaced with an efficient score.

{pstd}
Both {cmd:stage(1)} and {cmd:stage(2)} solve the same single DDF LP.
Use {cmd:direction(yb|xyb|x)} and omit {cmd:ort()}/{cmd:orientation()}.
FDH, bootstrap, imprecise data, prices, gx()/gy(), reference/saveref,
saveweights and trace cannot be combined with this option.
Unlike the default equality path, the new worker honors the actual solver:
{cmd:engine(simplex)} uses the package simplex; {cmd:engine(ipm)} uses Stata
LinearProgram; {cmd:engine(auto)} starts with simplex and tries native only
for unavailable DMUs. The default is simplex. {cmd:r(engine)} records the
selected label. {cmd:r(lp_backend)} is {cmd:package_simplex},
{cmd:stata_linearprogram}, or {cmd:simplex_then_native}, respectively.
The auto label describes the policy; per-DMU diagnostics identify the last
attempted backend. Existing solver implementations are reused unchanged.

{pstd}
{cmd:generate(name)} stores beta in {it:name} and status in {it:name}_status:
0 is an accepted LP result; 1 is unavailable. Both fields are missing outside
the sample. Existing result fields require {cmd:replacevars}; inputs, outputs,
bad outputs and IDs cannot be overwritten. {cmd:saving()} (or Form A {cmd:using})
saves selected ID, input/good/bad and result columns. {cmd:replace} authorizes
file replacement. A failed save preserves the caller's data.

{pstd}
Returns include {cmd:r(N)}, {cmd:r(N_complete)}, {cmd:r(n_failed)},
{cmd:r(N_eff)} (finite abs(beta)<=tolerance), {cmd:r(mean_beta)},
{cmd:r(mean_eff)} (alias), {cmd:r(n_native_retry)}, {cmd:r(n_native_recovered)},
{cmd:r(n_simplex_retry)}, {cmd:r(n_simplex_recovered)}, {cmd:r(n_native_attempted)},
{cmd:r(n_simplex)} and {cmd:r(n_native)} (accepted scores by backend),
{cmd:r(tolerance)}, and {cmd:r(acceptance_tolerance)}.
{cmd:r(formulation)} is {cmd:convex_joint_weak_disposal} and
{cmd:r(disposability)} is {cmd:weak}.

{pstd}
{cmd:weakmatrix} additionally returns N-by-N {cmd:r(active)}, {cmd:r(idle)},
and {cmd:r(intensity)} (= active+idle). Each row evaluates one DMU; each column
refers to a sample observation, mapped by the N-by-1 {cmd:r(obs_index)}.
{cmd:r(target_inputs)}, {cmd:r(target_outputs)}, and {cmd:r(target_bad)} are
{it:x0-beta*gx, y0+beta*gy, b0-beta*gb}, in original units. These are directional
targets, not slack-adjusted strongly efficient projections.
{cmd:r(input_slack)} and {cmd:r(output_slack)} are the remaining nonnegative
inequality slacks at the selected optimum; optima/intensities need not be unique.
No secondary slack maximization or component range is performed.

{pstd}
The N-by-10 {cmd:r(diagnostics)} columns are status, backend, simplex_retry, native_retry,
eq_abs, ineq_abs, bound_abs, eq_scaled, ineq_scaled, bound_scaled.
Backend is 1 (package simplex) or 2 (native), including failed attempts.
Each retry field is 0 (no retry), 1 (unavailable) or 2 (recovered).
The six residuals refer to the original LP after column-max normalization;
the first three are absolute residuals in those normalized coordinates.
The scaled checks use the common row/bound denominators described above.
Acceptance and simplex tolerance use min(tolerance,1e-9); native convergence
uses min(tolerance,1e-12). On simplex failure, one equivalent total/active
formulation (lambda=z+w, mu=z, 0<=mu<=lambda) is tried, then transformed back
and checked in the original active/idle LP. Existing bounded native retries
follow initial native failure. First accepted solutions are retained;
there is no RHS relaxation or score clipping. Small floating-point negatives
may remain within the acceptance tolerance. Residual checks are not an exact
optimality certificate. {cmd:weakmatrix} requires matsize >= max(N,m,s,q,10);
without it N-by-N intensity matrices are not stored.

{phang}
{bf:Weak-disposal gap recovery.}
{pstd}
For joint weak-disposal VRS with {cmd:direction(x)}, after unsuccessful primal
representations the simplex route may validate a dual upper bound and the
exactly feasible self point. It returns beta=0 only if that conservative bound
is within {it:min(tolerance,1e-9)} and the self point passes the original scaled
constraints. This narrow gap recovery does not relax the technology or treat
an arbitrary missing score as zero. {cmd:r(n_dual_attempted)} and
{cmd:r(n_dual_recovered)} count these DMUs. With {cmd:weakmatrix},
{cmd:r(dual_recovery)} has columns {cmd:state upper_bound dual_solver_status
 dual_margin_lower normalization_lower corrections dual_attempts}.
State 0 = not attempted, 1 = no validated bound, 2 = accepted gap recovery,
3 = validated bound too wide. The dual may retry with reversed peer order.
The bound uses conservative floating-point arithmetic; it is distinct from an
independent exact-rational certification of a particular dataset.

{pstd}
For example, the first observation below has beta=1/3 under VRS joint weak
disposal, versus beta=0 under the default VRS equality technology:

{* DEAS_HELP_EXAMPLES_BEGIN *}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input id x y bad}{p_end}
{phang2}{cmd:1 1 1 1}{p_end}
{phang2}{cmd:2 1 4 2}{p_end}
{phang2}{cmd:3 2 1 2}{p_end}
{phang2}{cmd:end}{p_end}
{phang2}{cmd:. dea x = y, model(undesirable) bad(bad) dmu(id) rts(vrs) disposability(weak) weakmatrix generate(wd)}{p_end}
{phang2}{cmd:. dea x = y, model(undesirable) bad(bad) dmu(id) rts(vrs) generate(eq)}{p_end}
{phang2}{cmd:. assert abs(wd-1/3)<1e-8 in 1}{p_end}
{phang2}{cmd:. assert abs(eq)<1e-8 in 1}{p_end}
{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END *}

{pstd}
Reference: Kuosmanen, T. 2005. Weak disposability in nonparametric production
analysis with undesirable outputs. {it:American Journal of Agricultural Economics}
87(4): 1077-1082. doi:10.1111/j.1467-8276.2005.00788.x.
Independent multiplier-dual and total/active-intensity LP comparisons,
closure witnesses, zeros, unit/order changes and syntax/output checks are in
{cmd:tests/test_weak_disposability.py}.

{dlgtab:Imprecise DEA (IDEA) options}

{pstd}
Imprecise DEA options use separate interval and legacy transformation paths
({help dea##CPY1999:Cooper, Park, and Yu 1999};
{help dea##Z2003:Zhu 2003}). {cmd:xlo()}, {cmd:xhi()}, {cmd:ylo()},
{cmd:yhi()}, {cmd:bounded()}, {cmd:ordinal()}, {cmd:xord()}, or {cmd:yord()}
activate IDEA when the model is omitted. Explicit {cmd:model(te)} also permits
this automatic selection; {cmd:model(imprecise)} requires at least one of these
options. Combining IDEA options with {cmd:fdh} or another explicit model is
rejected before results are changed. {cmd:r(model)} then reports {cmd:imprecise}, and
{cmd:r(model_requested)} retains the normalized model option (default {cmd:te}).
The ordinary radial VRS decomposition is not run after IDEA.

{phang}
{cmd:xlo(}{it:varlist}{cmd:)} and {cmd:xhi(}{it:varlist}{cmd:)} specify
lower and upper interval bounds for inputs, in the same order as the input
variables. {cmd:ylo(}{it:varlist}{cmd:)} and {cmd:yhi(}{it:varlist}{cmd:)}
specify the corresponding output bounds. When one bound list is omitted,
the main input or output variables supply that bound.

{pstd}
Both interval-only and mixed IDEA calls require finite endpoints with lower
<= upper for each selected observation. Missing endpoints (including extended
missing values) and reversed bounds return {cmd:r(459)} before creating IDs,
replacing results, or writing files. An omitted endpoint is checked using its
nominal model variable. When both endpoints are supplied, an unused nominal
value may be missing; it is still required if used in an ordinary model role.
Restrictions and missing ordinary/ordinal model values determine the sample
before this check. Invalid endpoints outside that sample are ignored.

{phang}
{cmd:bounded(}{it:varlist}{cmd:)} and {cmd:ordinal(}{it:varlist}{cmd:)}
select the legacy transformation, including when numeric bound lists are present.
{cmd:bounded()} accepts legacy string intervals such as {cmd:"[1,2]"}; a
numeric bounded variable supplies the same value for both endpoints.
{cmd:ordinal()} identifies numeric model variables transformed by their ordering.
Every bounded/ordinal variable must belong to the input or output list. A
transformed variable must have exactly one input/output role. Distinct bounded
and ordinal variables can be used together; one variable cannot be both.

{pstd}
String intervals are accepted for both inputs and outputs in both public
syntaxes. Each string model variable must be explicitly named in {cmd:bounded()};
all other model variables must be numeric. The string must contain exactly
{cmd:[lower,upper]} with two finite, ordered endpoints ({it:lower <= upper}).
Signed decimal and scientific notation, and whitespace around the interval
or its endpoints, are accepted. A numeric bounded value must be nonmissing.
Malformed strings, missing endpoints and reversed bounds stop with {cmd:r(459)}
before automatic IDs, result replacement or file writes. Validation applies
after {cmd:if}/{cmd:in} selection and exclusion of missing ordinary model
variables; invalid bounded values outside that sample are ignored. These
checks validate interval syntax and ordering, not feasibility of every DEA LP.

{phang}
{cmd:xord(}{it:varlist}{cmd:)} and {cmd:yord(}{it:varlist}{cmd:)} identify
ordinal inputs and outputs in the legacy path. Each variable in {cmd:xord()}
must be an input and each variable in {cmd:yord()} must be an output. Both
options activate IDEA on their own. Their variables are combined with
{cmd:ordinal()} into one unique list, so a variable named in multiple ordinal
options is transformed once.

{pstd}
{bf:Recompute earlier xord()/yord()-only results.} Before version 2.0.37,
these options alone could silently run ordinary TE. Explicit
{cmd:model(imprecise)} or another legacy activation option avoided that
particular routing defect. Repeated ordinal names could also apply the same
transformation twice; they now have the same meaning as a single occurrence.

{pstd}
{bf:Recompute earlier ordinal-output results.} Before version 2.0.39, an
output named in {cmd:ordinal()} or {cmd:yord()} could produce different scores
after observations were reordered. This also affected calls combining such
outputs with {cmd:bounded()}. The evaluated DMU's original comparison value is
now held in a double-precision scalar before transforming the data. The
existing strict comparison for inputs and inclusive comparison for outputs,
including treatment of ties, are retained. This corrects the in-place
transformation error; it does not establish equivalence with the separate
numeric interval formulation.

{pstd}
Numeric lists {cmd:xlo()/xhi()/ylo()/yhi()} can be combined with
{cmd:bounded()}, {cmd:ordinal()}, {cmd:xord()}, or {cmd:yord()} through the
legacy transformation. Lists retain the complete model-variable order. At
positions belonging to an ordinal or {cmd:bounded()} variable, each supplied
list must repeat that variable's own name as a placeholder. This also permits
string {cmd:bounded()} placeholders on the same side as numeric bounds.
Specifying a different endpoint variable at such a position is ambiguous and
stops with {cmd:r(198)}. Each numeric bounded dimension must occur once in the
input/output lists. Other endpoint variables must be numeric.

{pstd}
For example, {cmd:dea x q = y, xlo(x_low q) xhi(x_high q) ordinal(q)}
uses numeric bounds for {cmd:x} and the original ordering of {cmd:q}.
The evaluated DMU uses its lower input and upper output bounds, while reference
DMUs use upper input and lower output bounds, exactly as with string
{cmd:bounded()} data. All endpoint columns are copied before any dimension is
transformed, so a bound may safely refer to another numeric model variable.
The ordinal comparisons use the fixed values described above.

{pstd}
The mixed path requires finite, ordered numeric endpoints in each selected
observation. An omitted endpoint uses the nominal model variable; when both
endpoints are supplied, an unused nominal value may be missing. Ordinary and
ordinal model-variable missing values are excluded before this check. Invalid
selected endpoints stop with {cmd:r(459)} before automatic IDs, result
replacement or file writes; invalid endpoints outside that sample are ignored.
String bounds follow the same selection and rejection policy described above.
{cmd:r(idea_numeric_vars)} lists numeric bounded dimensions in model order.
Calls with numeric bound lists alone keep the separate interval calculation.

{pstd}
Earlier integrated versions could silently ignore ordinal options in numeric
mixed calls; versions 2.0.37 through 2.0.39 rejected the combination. Recompute
earlier mixed results using the legacy definition now specified for this path.

{pstd}
For interval IDEA, {cmd:generate(}{it:newvar}{cmd:)} creates {it:newvar}
(the pessimistic score) and {it:newvar}{cmd:_opt} (the optimistic score).
It creates neither {cmd:_pess} nor {cmd:_mid}.
{cmd:r(idea_method)} is {cmd:interval}. This path uses Mata's
{cmd:LinearProgram} regardless of {cmd:engine()}; selecting {cmd:auto}
does not trigger a second identical interval calculation. Interval LPs use
convergence tolerance min(tolerance(),1e-12) to reduce endpoint ordering errors;
the original-constraint acceptance threshold remains tolerance().

{pstd}
From version 2.0.56, numeric interval IDEA computes worst/best relative
efficiencies over independently varying input/output intervals. For the
pessimistic endpoint, the evaluated DMU uses upper inputs and lower outputs;
all other DMUs use lower inputs and upper outputs. For the optimistic endpoint,
these choices are reversed. In each LP the evaluated DMU also appears in the
reference set at its own chosen values. This is the consistent-realization
construction discussed by {help dea##SMD2006:Smirlis, Maragos, and Despotis (2006)},
also used for alpha-cut bounds by {help dea_fuzzy:dea_fuzzy}.
Input scores are theta; output scores are 1/phi, so a larger reported score
always means greater efficiency. This is not a fixed-frontier interval method.

{pstd}
For strictly positive, finite, ordered bounds, both endpoint LPs admit the
self-reference solution, and the mathematical scores satisfy
0 < pessimistic <= optimistic <= 1. A singleton has both endpoints equal to 1;
point intervals reduce to ordinary radial DEA. Widening the independent
intervals cannot narrow the efficiency range. These are model properties,
subject to numerical solution accuracy. Bounds describe data uncertainty,
not statistical confidence intervals. Dependence restrictions between uncertain
values are not imposed by these options.

{pstd}
For the legacy bounded/ordinal path, including numeric mixed bounds, {it:newvar} is the transformed DEA
score. Existing {it:newvar}{cmd:_opt} and {it:newvar}{cmd:_mid} variables
are retained for compatibility, but both are copies of that same score.
They are not independently computed optimistic bounds or interval midpoints.
{cmd:r(idea_method)} is {cmd:legacy}; the transformed DEA calls use the
selected engine. {cmd:r(idea_vars)} lists all generated IDEA score columns,
and {cmd:saving()} includes them.

{pstd}
For both paths, {cmd:r(mean_eff)}, {cmd:r(min_eff)}, {cmd:r(max_eff)}, and
{cmd:r(N_eff)} summarize the main generated score. {cmd:r(mean_opt)},
{cmd:r(min_opt)}, {cmd:r(max_opt)}, and {cmd:r(N_eff_opt)} summarize
{it:newvar}{cmd:_opt}. Each efficiency count uses finite scores within
{cmd:tolerance()} of 1; it is not a classification of every optimistic
score above 1. Statistics exclude missing values in the corresponding column.
{cmd:r(n_failed)} counts selected rows with either main or {cmd:_opt} score
missing, counting each row once. Separate {cmd:r(n_failed_primary)} and
{cmd:r(n_failed_opt)} counts identify the affected columns. If all values in
a column are missing, its mean/minimum/maximum are missing and its efficiency
count is zero. These return-value conventions do not establish equivalence
between the interval and legacy methodologies.
Recompute numeric interval IDEA results from versions before 2.0.56: both
endpoints used an inconsistent self-reference row. Earlier optimistic scores
could exceed 1 or be missing even with positive, ordered data, and pessimistic
scores could also be incorrect. No score is clipped or substituted with 1;
an unsuccessful numerical solve still returns missing and contributes to the
failure counts. The legacy bounded/ordinal transformation is unchanged, but
from 2.0.56 its internal DEA calls also receive the requested {cmd:tolerance()};
earlier versions silently used the default within those calls.

{dlgtab:Output}

{phang}
{cmd:generate(}{it:newvar}{cmd:)} specifies the variable name for storing
efficiency scores. The default name is {cmd:dea_efficiency}. Interval IDEA
also creates {it:newvar}{cmd:_opt}; the legacy bounded/ordinal path creates
{it:newvar}{cmd:_opt} and {it:newvar}{cmd:_mid} as copies of the main score.
See {it:Imprecise DEA options} for the distinction.

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies a string or numeric variable
containing DMU identifiers for display purposes. If omitted, {cmd:_dmu} is
created from observation order. An existing {cmd:_dmu} is never overwritten:
specify {cmd:dmu(_dmu)} to use it or supply another identifier variable.

{pstd}
For core {cmd:dea}, each selected observation is evaluated separately even if
its DMU label is repeated or missing.  In version 2.0.35, the VRS decomposition
matrix follows the same selected-observation order as {cmd:r(efficiency)};
duplicate labels are not averaged and missing labels do not remove rows.
Its row names are {cmd:r1}, {cmd:r2}, ..., not DMU labels.  Use
{cmd:r(vrsobs)} to map each matrix row to the observation number in the dataset
at the time of the call.  This number is not a permanent identifier after
sorting or changing the dataset.  Code that previously matched VRS matrix
rows by DMU row names should use this mapping instead.

{pstd}
From version 2.0.50, score and reference tables display string DMU labels
literally, including quotation marks and text resembling Stata macros or
formatting directives. This also applies to the FDH and numeric interval
IDEA score tables. The original identifier values are preserved. Numeric
identifiers retain their general-format display.
Saved reference names replace backticks as needed for valid Stata names;
the original peer ID remains in the {cmd:dea_peer_id} variable characteristic.

{phang}
{cmd:reference} displays the reference set (peer DMUs and their weights) for
each DMU in the results table. Display values are local to the command;
existing user scalars and matrices are preserved, and no public scratch
objects are left behind. Data variables named {cmd:_dmu_k}, {cmd:_dmu_j},
or {cmd:_lam_val} do not alter the displayed peer labels or weights.
Earlier versions could overwrite objects with these names or display incorrect
reference lists when same-named data variables were present. The correction
in 2.0.43 changes display handling, not the computed scores or weights.

{phang}
{cmd:saveref} requires {cmd:reference} and is supported with radial
{cmd:model(te)} or {cmd:model(supereff)}, without FDH or imprecise-data options.
It saves the final
engine's reference weights as persistent double variables. Names start with
{cmd:_ref_} followed by a Stata-compatible form of the peer DMU identifier;
duplicate names after conversion/truncation receive numeric suffixes.
Existing reference variables are never overwritten; explicitly remove only
reference variables that you intend to replace. If {cmd:generate()} itself
conflicts with a planned reference name, choose a different score variable.
The {cmd:replace} option controls the output file, not variables in memory.
Variables are missing outside the evaluation sample. Their labels and
{cmd:dea_peer_id}, {cmd:dea_peer_index}, and {cmd:dea_peer_dmu} characteristics
record the mapping; {cmd:r(refvars)} gives variable names in lambda-column
order. With {cmd:saving()}, the saved dataset contains DMU, input, output,
score, and these persistent reference variables, without dispatcher tempvars.
This adds one variable per peer DMU and can require substantial memory for
large evaluation samples. Optimal peer weights need not be unique.

{phang}
{cmd:saveweights} saves optimal multiplier weights as new variables after
{cmd:stage(2)}. Input and output weights are named {cmd:_dea_v#} and
{cmd:_dea_u#}; the signed intercept is {cmd:_dea_u0}.  For input orientation,
the Cooper-signed supporting inequalities are {it:u'y_j-v'x_j-u0} <= 0 and
the objective identity is {it:efficiency=u'y_o-u0}.  For output orientation,
they are
{it:v'x_j-u'y_j+u0} >= 0 and
{it:1/efficiency=v'x_o+u0}. Available with {cmd:model(multiplier)}. With
{cmd:stage(1)}, only the efficiency score is reported.

{phang}
{cmd:saving(}{it:filename}{cmd:)} saves the evaluation sample's DMU identifier,
inputs, outputs, and final efficiency score to a Stata dataset. The file is
written once after any automatic engine retry. Persistent reference variables
are included with {cmd:saveref}; multiplier weight variables are included with
{cmd:model(multiplier) stage(2) saveweights}. All IDEA score columns listed
in {cmd:r(idea_vars)} are also included. Other user variables and internal
temporary variables are excluded. Slack matrices, when available, are returned
in {cmd:r(slack_in)} and {cmd:r(slack_out)} and are not automatically saved as
dataset variables. A file-write failure does not roll back generated variables
in memory.

{phang}
In the {it:inputvars} {cmd:=} {it:outputvars} syntax, legacy
{cmd:using} {it:filename} is an alias for {cmd:saving()}. Specify one form only.

{phang}
{cmd:replace} permits {cmd:saving()} to overwrite an existing dataset.

{dlgtab:Bootstrap (Simar-Wilson)}

{phang}
{cmd:bootstrap} computes smoothed-bootstrap inference for radial
{cmd:model(te)} using the same calculation routines as {helpb dea_bootstrap},
following the pseudo-frontier approach of
{help dea##SW1998:Simar and Wilson (1998)}. Each replication constructs
pseudo inputs or outputs and solves DEA again against that reference set.
Both orientations and CRS/VRS/NIRS/NDRS are supported. FDH, imprecise data,
and other model() choices are not supported with this option.

{phang}
{cmd:reps1(}{it:#}{cmd:)} is the number of bootstrap replications used for
both bias estimation and confidence intervals. The default is {cmd:100};
at least 50 are required. Increase it for empirical work. The legacy name
{cmd:reps1()} does not imply a second bootstrap stage.

{phang}
{cmd:level(}{it:#}{cmd:)} is the confidence level as a percentage, strictly
between 0 and 100; decimals are allowed. The default is {cmd:95}.
Intervals use empirical inverse-CDF quantiles of the bootstrap deviations.
This option affects the interval endpoints, not the bias or random draws.

{pstd}
With {cmd:generate(eff)}, the additional variables are {cmd:eff_bc}
(bias-corrected score), {cmd:eff_bias} (original minus corrected score),
{cmd:eff_se} (bootstrap standard error), and {cmd:eff_lo}/{cmd:eff_hi}
(pivot confidence limits). They are included in {cmd:saving()}.
All results use the stored efficiency scale: input theta, or output 1/phi.
Corrections and interval endpoints are not clipped to [0,1]. They may fall
outside that range; this is not a new feasible DEA point estimate.

{pstd}
{cmd:if}/{cmd:in} restrict evaluation and reference DMUs together.
Bootstrap result variables are missing outside this sample. Existing output
columns require {cmd:replacevars}, including all five new suffixes.
Use {cmd:set seed} before the call for reproducible draws.

{pstd}
{cmd:engine()} selects the initial DEA point-score engine. Bootstrap
pseudo-frontier LPs always use the package simplex engine, as does
{cmd:dea_bootstrap}. Original point scores and reference weights are retained.
The bootstrap requires positive, nonmissing initial scores for every selected DMU.

{pstd}
With {cmd:bootstrap envars(Z)}, the command instead runs sequential
Simar-Wilson (2007) Algorithm 2 through {helpb dea_sw2007}.
{cmd:reps1()} controls the score-correction bootstrap (minimum 50), and
{cmd:reps2()} controls the subsequent regression bootstrap (default 2000;
minimum 2). Both estimation and reference samples exclude missing Z values.
The original engine scores are retained, and the only additional score column
is {cmd:generate()}_bc, obtained by correcting inverse efficiency and then
inverting back. It is included in {cmd:saving()}; {cmd:r(bootstrap)} has one
column named bias_corrected. Coefficients, bootstrap covariance and intervals
are returned in {cmd:e()}, with {cmd:e(cmd)} equal to {cmd:dea_sw2007}.
The normal {cmd:r()} DEA results remain available. Without {cmd:envars()},
{cmd:reps2()} is an error. Use {cmd:dea_sw2007} directly to select the legacy
{cmd:procedure(joint)} or the public {cmd:backend(simarwilson)}.

{pstd}
These calculations use a smoothed radial bootstrap. Numerical agreement of
LP solutions does not establish nominal confidence-interval coverage in a
particular data-generating process; modeling assumptions and sample size
remain relevant to empirical inference.

{pstd}
The common radial smoother does not account for heterogeneous efficiency
distributions or dependence within clusters. The SW2007 path uses independent,
constant-variance truncated-normal errors and a separable common technology;
its period-specific panel extension is not a cluster bootstrap. Bias correction
can increase RMSE for some DMUs, and larger bootstrap counts do not guarantee
nominal coverage. See {help dea_bootstrap##remarks:dea_bootstrap remarks} and
{help dea_sw2007##remarks:dea_sw2007 remarks} for the relevant assumptions.

{pstd}
The smoother samples the efficiency scores and their mirror images around 1,
adds Gaussian noise, and corrects the mirrored distribution's variance around
its mean of 1 before reflecting values back below 1. This ordering preserves
support at the efficiency boundary. Earlier versions rescaled after reflection,
which artificially separated bootstrap draws from 1 and could severely distort
bias corrections and interval coverage. Inference results from those versions
should be recomputed.

{dlgtab:Technical}

{phang}
{cmd:engine(}{it:enginename}{cmd:)} specifies the LP solver engine.
{cmd:simplex} (the default) uses the package-supplied Two-Phase Revised
Simplex method implemented in Mata.  {cmd:ipm} uses Stata's {cmd:LinearProgram()}
interior point method, available since Stata 16.  {cmd:auto} normally tries
simplex first and re-solves all DMUs with IPM if any score is missing.  For
radial technical efficiency, it instead selects IPM before solving when the
ratio between the largest and smallest nonzero variable scales exceeds
{cmd:1e8}.  Explicit {cmd:engine(simplex)} is never overridden.  Model-specific domains and backend limits still apply. Cooper auto is rejected
and its certified simplex path has a bounded domain. Unusually heterogeneous
units may require an explicit engine choice or rescaling.

{phang}
{cmd:method(}{it:enginename}{cmd:)} is a synonym for {cmd:engine()}.
If both are specified, {cmd:engine()} takes precedence.

{phang}
{cmd:tolerance(}{it:#}{cmd:)} specifies the numerical tolerance for
convergence and efficiency classification. The default is {cmd:1e-7}.
The effective value must be positive and finite. Zero, negative values, and
Stata missing values ({cmd:.}, {cmd:.a}, ..., {cmd:.z}) return {cmd:r(198)}
before sample evaluation, automatic IDs, result replacement, or file writing.
Both public syntax forms and all core model paths use this validation.
For the TE residual/convergence policy, see {cmd:stage()} above.
From version 2.0.55, core non-TE LP paths also receive this value: economic
and allocative efficiency, DDF, additive, SBM and super-SBM, super-efficiency,
multiplier weights, numeric IDEA, congestion, and SBM with undesirable outputs.
Native non-TE convergence uses {it:min(tolerance,1e-8)}; undesirable DDF retains
its stricter {it:min(tolerance,1e-12)} cap. These are convergence settings,
not separate user-adjustable pivot thresholds. Simplex receives the requested
value and retains its internal pivot stability rules.

{pstd}
A finite solver result must also pass equality, inequality and bound checks
on the submitted LP, before simplex standardization or native optimization.
For each row, the residual is divided by
{it:max(|b|,sum_j |a_j| max(1,|x_j|))}; zero rows with zero RHS have zero
residual. Bounds use {it:max(1,|bound|,|x_j|)}. The maximum scaled violation
must not exceed {cmd:tolerance()}. A native solve that reports convergence
but fails these checks gets at most one restart of the identical LP, with
convergence tightened to {it:min(tolerance,1e-12)}, only if that is stricter.
There is no restart here after native nonconvergence; the existing engine
fallback remains separate. This native residual-restart rule does not
reoptimize accepted initial solutions. The simplex successful-basis and
stage-link checks described separately have their own bounded acceptance rules.
Rejected scores stay missing and enter the
existing failure counts; automatic fallback retains its documented scope.
This primal test does not certify optimality or guarantee an objective-value
error smaller than {cmd:tolerance()}. A failure does not prove mathematical
infeasibility. In particular, VRS super-efficiency can be truly infeasible.
Positive data are not rounded to zero using this option; SBM-undesirable
normalizes columns and corresponding slack units before forming its ratios.
Small positive multiplier weights are also retained: setting them to zero can
break normalization when data are large. Unavailable native multiplier weights
remain missing, even when the separately computed radial score is available;
{cmd:r(n_failed)} counts missing main scores, not missing weight vectors.

{pstd}
For the common LP model display, a finite score is classified as efficient
when its distance from the frontier is less than or equal to the tolerance:
0 for {cmd:ddf}, {cmd:undesirable}, and {cmd:additive}; 1 for the other ordinary
models.  For {cmd:supereff} and {cmd:supersbm}, {cmd:r(N_eff)} instead counts
finite scores greater than or equal to {it:1-tolerance}, including ties at 1;
{cmd:r(N_super)} counts scores strictly greater than {it:1+tolerance}.
These are score-based classifications; they do not replace checks for radial
slacks or establish strong Pareto efficiency.

{pstd}
The common display ranks the three zero-frontier models from smallest to
largest score, and the other models from largest to smallest.  Scores within
tolerance of the frontier share a ranking key at that frontier; elsewhere
the key is rounded to the tolerance grid.  Competition ranks are used
(for example, 1, 1, 3).  Missing scores have no rank and are excluded from
the efficiency counts.  These operations never change the stored scores.
The separate FDH, imprecise, congestion, and SBM-undesirable displays are
outside this common-display rule.

{pstd}
Legacy Form A accepts {cmd:tol1()} only as a deprecated alias for
{cmd:tolerance()}; its supplied value, including {cmd:1e-14}, is used literally.
Specifying both options returns {cmd:r(198)}, even when their values match.
{cmd:tol2()} also returns {cmd:r(198)} with migration guidance: the integrated
engine no longer implements the old second threshold. Earlier releases silently
ignored {cmd:tol2()} and could override an explicit {cmd:tolerance(1e-7)} based
on {cmd:tol1()}. Use {cmd:tolerance()} alone in new scripts and inspect
{cmd:r(tolerance)} for the effective value. Form B supports the modern option.

{phang}
{cmd:trace} displays iteration-level details from the LP solver.

{phang}
{cmd:debug} displays additional debugging information including problem
formulation details.


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:dea} version 2.0 features a dual-engine architecture.  The engines cover the common formulation family, with model-specific domains
and backend limits. The explicit Cooper path has the separate certification,
selection and auto restrictions described above.

{pstd}
{it:Two-Phase Revised Simplex engine.}  The default engine is a
complete Mata implementation of the revised simplex method, requiring no
external dependencies beyond base Mata.  It employs Dantzig pricing
with automatic Bland's anti-cycling rule activation, Harris two-pass
ratio test for numerical stability, eta-vector product-form basis
updates with periodic LU refactorization, and a general-purpose LP
converter that transforms inequality/equality formulations to
standard form.  Non-radial models (cost, revenue, profit, DDF, SBM,
additive, super-efficiency, super-SBM, allocative) are solved through the
converter infrastructure.

{pstd}
Simplex 2.0.10 checks basis inverses and propagates failed pivots or
refactorizations as numerical failures. Remaining dependent Phase-I rows
are removed without arbitrary replacement of basic variables. Bland uses
minimum-ratio ties; Harris selection limits the basic-feasibility loss.
One bounded retry after numerical failure rebuilds the initial basis,
equilibrates rows, columns and the objective, and starts with Bland pricing.
Recovered variables and the objective are returned in original units and
the original constraints are checked. This retry keeps the fixed radial
target. It requires no external solver. {cmd:trace} includes pivot choices
and the restart notice. Floating-point acceptance is not an exact
feasibility or optimality certificate.

{pstd}
{it:Interior Point Method (IPM) engine.}  The alternative engine leverages
Stata 16's {cmd:LinearProgram()} class, which implements a primal-dual
interior point method.  IPM offers polynomial-time worst-case complexity,
efficient handling of degeneracy common in DEA, and automatic
primal-dual fallback when the primal formulation encounters numerical
difficulties.

{pstd}
Cross-engine checks compare efficiency scores, projection feasibility and
slack objectives. Solver status and projection availability can differ on
numerical boundary cases. Optimal reference weights need not be identical
when the LP has multiple optima. A returned radial score alone does not
establish that the requested stage-2 projection was completed; inspect
the documented fallback behavior and available status returns.

{pstd}
{it:Automatic engine selection.} Cooper rejects auto; this paragraph describes
the ordinary model routes. When {cmd:engine(auto)} is specified,
the command normally first runs all DMUs with the simplex engine.  If any DMU
produces a missing efficiency score (e.g., due to numerical difficulty),
the entire problem is re-solved from scratch with the IPM engine.  This
provides the speed advantage of the simplex engine (typically
22--57{it:x} faster) while retaining the robustness of IPM as a
safety net.  The {cmd:r(engine)} macro reports which engine produced
the final results.  For radial {cmd:model(te)}, version 2.0.22 also compares
the sample maxima of the input and output variables.  A dimensional scale
ratio above {cmd:1e8} selects IPM before the first solve to protect unit
invariance.  The decision is reported in {cmd:r(auto_preemptive_ipm)}.

{pstd}
Beginning with version 2.0.28, radial {cmd:model(te)} envelopment calculations
also scale each input and output column internally before either solver is
called.  This positive change of units leaves radial efficiency scores and
reference intensities unchanged, protects explicit {cmd:engine(simplex)} and
{cmd:engine(ipm)} calls when variable magnitudes differ greatly, and converts
stage-2 slacks back to the source variables' units.  The input and output
variables in memory are not modified.  This statement applies to the radial
technical-efficiency path; other model-specific solvers retain their documented
implementations.

{pstd}
{it:Backward compatibility.} The standalone {help dea_lp} entry point works
in a fresh session and calls {cmd:dea} with the former VRS/IPM defaults.
Explicit RTS and engine/method options override those defaults; other behavior
uses the current command, including {cmd:replacevars}. Ordinary {cmd:dea}
retains CRS/simplex defaults. This alias is for DEA; independent general LP
and MILP problems use {help lp} and {help milp_bb}.

{pstd}
{it:Model selection guide:}

{p2colset 9 32 34 2}{...}
{p2col :{it:Goal}}{it:Recommended model}{p_end}
{p2col :{hline 52}}{p_end}
{p2col :Radial efficiency}model(te){p_end}
{p2col :Multiplier weights}model(multiplier){p_end}
{p2col :Cost minimization}model(cost) with w(){p_end}
{p2col :Revenue maximization}model(revenue) with p(){p_end}
{p2col :Profit maximization}model(profit) with w() and p(){p_end}
{p2col :Directional distance}model(ddf) with gx() gy(){p_end}
{p2col :Non-radial efficiency}model(sbm){p_end}
{p2col :Slack analysis}model(additive){p_end}
{p2col :Radial super-efficiency}model(supereff){p_end}
{p2col :Non-radial super-efficiency}model(supersbm){p_end}
{p2col :Allocative decomposition}model(allocative) with w() or p(){p_end}
{p2col :No convexity}model(fdh){p_end}
{p2col :Interval/ordinal data}model(imprecise) with xlo() xhi(){p_end}
{p2col :Undesirable outputs (DDF)}model(undesirable) with bad(){p_end}
{p2col :Input congestion}model(congestion){p_end}
{p2col :Undesirable outputs (SBM)}model(sbmundesirable) with bad(){p_end}
{p2colreset}{...}


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

{pstd}
The check includes VRS columns {cmd:CRS_TE}, {cmd:VRS_TE}, {cmd:SCALE},
{cmd:RTS_NUM}, {cmd:RTS}, requested reference/multiplier weights, and model-specific
suffix columns. {cmd:replace} continues to control only {cmd:saving()}; it does
not replace variables. Use {cmd:replacevars} to explicitly refresh existing
results, and {cmd:replace} as well if an existing saved dataset must be replaced.
An existing {cmd:_dmu} is never replaced automatically: specify {cmd:dmu(_dmu)}
to reuse it, or provide another identifier.

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


    {title:Quick start with a familiar Stata dataset}

{pstd}{cmd:auto} is useful for command practice, not a substantive automobile
production study.  The selected quantities are strictly positive.{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. summarize price weight mpg trunk}{p_end}
{phang2}{cmd:. dea price weight = mpg trunk}{p_end}
{phang2}{cmd:. return list}{p_end}

{pstd}With no options, the Results window reports the effective defaults:
{cmd:model(te)}, {cmd:rts(crs)}, {cmd:orientation(input)}, {cmd:stage(2)},
{cmd:engine(simplex)}, and {cmd:tolerance(1e-7)}.  It also reports that a DMU
identifier named {cmd:_dmu} was generated from observation order.  In saved
research code, specify the model explicitly and use a verified identifier.{p_end}

{pstd}The following repeated analyses specify {cmd:dmu()} explicitly to reuse
their identifiers without replacing an automatically generated {cmd:_dmu}.{p_end}

    {title:Setup}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str7 dmu i_employee i_area o_sales o_profits}{p_end}
{phang2}{cmd:"store A" 10 20  70 6}{p_end}
{phang2}{cmd:"store B" 15 15 100 3}{p_end}
{phang2}{cmd:"store C" 20 30  80 5}{p_end}
{phang2}{cmd:"store D" 25 15 100 2}{p_end}
{phang2}{cmd:"store E" 12  9  90 8}{p_end}
{phang2}{cmd:"store F" 16 18  75 5}{p_end}
{phang2}{cmd:end}{p_end}

    {title:Basic radial efficiency}

{pstd}CRS input-oriented (CCR model){p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) ort(i) gen(theta) dmu(dmu)}{p_end}

{pstd}Invalid tolerances preserve existing results; valid values are returned{p_end}
{phang2}{cmd:. clonevar theta_before = theta}{p_end}
{phang2}{cmd:. capture noisily dea i_employee i_area = o_sales o_profits, tolerance(0) gen(theta) replacevars dmu(dmu)}{p_end}
{phang2}{cmd:. assert _rc==198}{p_end}
{phang2}{cmd:. assert theta==theta_before}{p_end}
{phang2}{cmd:. capture noisily dea o_sales o_profits, inputs(i_employee i_area) tolerance(.z) gen(theta) replacevars dmu(dmu)}{p_end}
{phang2}{cmd:. assert _rc==198}{p_end}
{phang2}{cmd:. assert theta==theta_before}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, tolerance(1e-8) gen(theta_tol) dmu(dmu)}{p_end}
{phang2}{cmd:. assert r(tolerance)==1e-8}{p_end}
{phang2}{cmd:. assert !missing(theta_tol)}{p_end}
{phang2}{cmd:. drop theta_before}{p_end}

{pstd}VRS input-oriented (BCC model) with reference sets{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(vrs) ort(i) gen(te_vrs) dmu(dmu) reference}{p_end}
{phang2}{cmd:. matrix list r(vrsfrontier)}{p_end}
{phang2}{cmd:. matrix list r(vrsobs)}{p_end}

{pstd}Output-oriented{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(vrs) ort(o) gen(te_out) dmu(dmu) replacevars}{p_end}

    {title:Non-radial models}

{pstd}Slacks-based measure{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(sbm) rts(crs) gen(sbm_eff) dmu(dmu)}{p_end}

{pstd}Additive model{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(additive) rts(vrs) gen(add) dmu(dmu)}{p_end}

    {title:Economic efficiency}

{pstd}Cost efficiency with input prices{p_end}
{phang2}{cmd:. gen w1 = 10}{p_end}
{phang2}{cmd:. gen w2 = 2}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(cost) w(w1 w2) rts(vrs) gen(ce) dmu(dmu)}{p_end}

{pstd}Revenue efficiency with output prices{p_end}
{phang2}{cmd:. gen p1 = 5}{p_end}
{phang2}{cmd:. gen p2 = 8}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(revenue) p(p1 p2) gen(re) dmu(dmu)}{p_end}

{pstd}Profit efficiency with input and output prices{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(profit) w(w1 w2) p(p1 p2) rts(vrs) gen(pe) dmu(dmu)}{p_end}

{pstd}Allocative efficiency decomposition{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(allocative) w(w1 w2) rts(vrs) gen(ae) dmu(dmu)}{p_end}

    {title:Directional distance function}

{pstd}Default direction (proportional){p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(ddf) rts(vrs) gen(beta) dmu(dmu)}{p_end}

{pstd}Custom direction vector{p_end}
{phang2}{cmd:. gen gx1 = 1}{p_end}
{phang2}{cmd:. gen gx2 = 1}{p_end}
{phang2}{cmd:. gen gy1 = 1}{p_end}
{phang2}{cmd:. gen gy2 = 1}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(ddf) gx(gx1 gx2) gy(gy1 gy2) gen(b_unit) dmu(dmu)}{p_end}

{pstd}Missing custom directions exclude the same evaluation and reference rows{p_end}
{phang2}{cmd:. clonevar gx1_before = gx1}{p_end}
{phang2}{cmd:. replace gx1 = .a in 1}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(ddf) gx(gx1 gx2) gy(gy1 gy2) gen(dir_complete) dmu(dmu)}{p_end}
{phang2}{cmd:. assert r(N)==_N-1}{p_end}
{phang2}{cmd:. assert missing(dir_complete) in 1}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits if !missing(gx1), model(ddf) gx(gx1 gx2) gy(gy1 gy2) gen(dir_manual) dmu(dmu)}{p_end}
{phang2}{cmd:. assert dir_complete==dir_manual}{p_end}
{phang2}{cmd:. replace gx1 = gx1_before}{p_end}
{phang2}{cmd:. drop gx1_before}{p_end}

    {title:Super-efficiency}

{pstd}Radial super-efficiency{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(supereff) rts(crs) ort(i) gen(se) dmu(dmu)}{p_end}
{phang2}{cmd:. display "Efficient / strictly super-efficient: " r(N_eff) " / " r(N_super)}{p_end}

{pstd}Super-SBM{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(supersbm) rts(vrs) gen(ssbm) dmu(dmu)}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(supersbm) ort(input) rts(crs) gen(ssbm_in) dmu(dmu)}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(supersbm) ort(output) rts(crs) gen(ssbm_out) dmu(dmu)}{p_end}

    {title:Multiplier (dual) model}

{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(multiplier) rts(crs) ort(i) gen(mult) saveweights dmu(dmu)}{p_end}

    {title:Free disposal hull}

{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(fdh) gen(fdh_eff) dmu(dmu)}{p_end}

{pstd}Shorthand using {cmd:fdh} flag{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, fdh gen(fdh_eff2) dmu(dmu)}{p_end}

    {title:Imprecise DEA (IDEA)}

{pstd}Interval data with known lower/upper bounds on inputs{p_end}
{phang2}{cmd:. gen emp_lo = i_employee * 0.9}{p_end}
{phang2}{cmd:. gen emp_hi = i_employee * 1.1}{p_end}
{phang2}{cmd:. gen area_lo = i_area * 0.9}{p_end}
{phang2}{cmd:. gen area_hi = i_area * 1.1}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(imprecise) xlo(emp_lo area_lo) xhi(emp_hi area_hi) gen(idea) dmu(dmu)}{p_end}

{pstd}Supplied endpoints permit an unused nominal missing value{p_end}
{phang2}{cmd:. clonevar employee_before = i_employee}{p_end}
{phang2}{cmd:. replace i_employee = .z in 1}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, xlo(emp_lo area_lo) xhi(emp_hi area_hi) gen(idea_nominal) dmu(dmu)}{p_end}
{phang2}{cmd:. assert r(N)==_N}{p_end}
{phang2}{cmd:. assert idea_nominal==idea & idea_nominal_opt==idea_opt}{p_end}
{phang2}{cmd:. replace i_employee = employee_before}{p_end}
{phang2}{cmd:. drop employee_before}{p_end}

{pstd}Invalid selected endpoints leave existing results intact{p_end}
{phang2}{cmd:. clonevar emp_lo_before = emp_lo}{p_end}
{phang2}{cmd:. replace emp_lo = .a in 1}{p_end}
{phang2}{cmd:. capture noisily dea i_employee i_area = o_sales o_profits, xlo(emp_lo area_lo) xhi(emp_hi area_hi) gen(idea) replacevars dmu(dmu)}{p_end}
{phang2}{cmd:. assert _rc==459}{p_end}
{phang2}{cmd:. assert idea==idea_nominal & idea_opt==idea_nominal_opt}{p_end}
{phang2}{cmd:. replace emp_lo = emp_lo_before}{p_end}
{phang2}{cmd:. drop emp_lo_before}{p_end}

{pstd}Interval data on both inputs and outputs{p_end}
{phang2}{cmd:. gen sales_lo = o_sales * 0.95}{p_end}
{phang2}{cmd:. gen sales_hi = o_sales * 1.05}{p_end}
{phang2}{cmd:. gen prof_lo = o_profits * 0.9}{p_end}
{phang2}{cmd:. gen prof_hi = o_profits * 1.1}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, xlo(emp_lo area_lo) xhi(emp_hi area_hi) ylo(sales_lo prof_lo) yhi(sales_hi prof_hi) gen(idea2) dmu(dmu)}{p_end}

{pstd}Legacy string intervals on inputs and outputs, with a distinct ordinal
input. Both syntax forms produce the same score columns.{p_end}
{phang2}{cmd:. gen strL emp_bound = "[" + strtrim(string(emp_lo,"%21.17g")) + "," + strtrim(string(emp_hi,"%21.17g")) + "]"}{p_end}
{phang2}{cmd:. gen strL sales_bound = "[" + strtrim(string(sales_lo,"%21.17g")) + "," + strtrim(string(sales_hi,"%21.17g")) + "]"}{p_end}
{phang2}{cmd:. dea emp_bound i_area = sales_bound o_profits, bounded(emp_bound sales_bound) ordinal(i_area) gen(bound_a) dmu(dmu)}{p_end}
{phang2}{cmd:. dea sales_bound o_profits, inputs(emp_bound i_area) bounded(emp_bound sales_bound) ordinal(i_area) gen(bound_b) dmu(dmu)}{p_end}
{phang2}{cmd:. assert bound_a == bound_b}{p_end}

{pstd}Equivalent numeric bounds with the ordinal input repeated as a positional
placeholder. The main score and its two compatibility copies agree.{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, xlo(emp_lo i_area) xhi(emp_hi i_area) ylo(sales_lo o_profits) yhi(sales_hi o_profits) ordinal(i_area) gen(mix_numeric) dmu(dmu)}{p_end}
{phang2}{cmd:. assert mix_numeric == bound_a}{p_end}
{phang2}{cmd:. assert mix_numeric == mix_numeric_opt & mix_numeric == mix_numeric_mid}{p_end}
{phang2}{cmd:. return list}{p_end}

{pstd}Ordinal inputs activate IDEA without an explicit model option{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, xord(i_employee) gen(ord_in) dmu(dmu)}{p_end}

{pstd}Ordinal inputs and outputs together, using standard syntax. Reversing
observation order retains scores for the same DMUs within numerical tolerance.{p_end}
{phang2}{cmd:. dea o_sales o_profits, inputs(i_employee i_area) xord(i_employee) yord(o_profits) rts(vrs) gen(ord_both) dmu(dmu)}{p_end}
{phang2}{cmd:. gen long ord_visit = _n}{p_end}
{phang2}{cmd:. gsort -ord_visit}{p_end}
{phang2}{cmd:. dea o_sales o_profits, inputs(i_employee i_area) xord(i_employee) yord(o_profits) rts(vrs) gen(ord_reverse) dmu(dmu)}{p_end}
{phang2}{cmd:. assert missing(ord_both) == missing(ord_reverse)}{p_end}
{phang2}{cmd:. assert abs(ord_both - ord_reverse) < 1e-7 if !missing(ord_both)}{p_end}
{phang2}{cmd:. sort ord_visit}{p_end}
{phang2}{cmd:. drop ord_visit}{p_end}

    {title:Undesirable outputs (DDF approach)}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input capital labor electricity so2}{p_end}
{phang2}{cmd:. 120 34 400 42}{p_end}
{phang2}{cmd:. 150 40 470 55}{p_end}
{phang2}{cmd:. 135 37 430 48}{p_end}
{phang2}{cmd:. 180 48 520 66}{p_end}
{phang2}{cmd:. 110 32 385 40}{p_end}
{phang2}{cmd:. 165 43 500 58}{p_end}
{phang2}{cmd:. end}{p_end}
{phang2}{cmd:. gen long plant = _n}{p_end}

{pstd}Contract bad outputs (SO2) while expanding good outputs{p_end}
{phang2}{cmd:. dea capital labor = electricity, model(undesirable) bad(so2) rts(vrs) gen(ddf_u) dmu(plant)}{p_end}

{pstd}Contract inputs and bad outputs simultaneously{p_end}
{phang2}{cmd:. dea capital labor = electricity, model(undesirable) bad(so2) direction(xyb) gen(ddf_xyb) dmu(plant)}{p_end}

{pstd}Contract inputs while holding good and bad outputs fixed:{p_end}
{phang2}{cmd:. dea capital labor = electricity, model(undesirable) bad(so2) direction(x) gen(ddf_x) dmu(plant)}{p_end}

    {title:SBM with undesirable outputs (Tone 2004)}

{pstd}Non-radial efficiency with bad outputs{p_end}
{phang2}{cmd:. dea capital labor = electricity, model(sbmundesirable) bad(so2) rts(vrs) gen(sbmu) dmu(plant)}{p_end}

{pstd}Missing bad outputs are excluded before solving{p_end}
{phang2}{cmd:. clonevar so2_before = so2}{p_end}
{phang2}{cmd:. replace so2 = .z in 1}{p_end}
{phang2}{cmd:. dea capital labor = electricity, model(sbmundesirable) bad(so2) rts(vrs) gen(sbmu_complete) dmu(plant)}{p_end}
{phang2}{cmd:. assert r(N)==_N-1}{p_end}
{phang2}{cmd:. assert missing(sbmu_complete) in 1}{p_end}
{phang2}{cmd:. replace so2 = so2_before}{p_end}
{phang2}{cmd:. drop so2_before}{p_end}

    {title:Congestion analysis (Cooper et al. 2004)}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str10 hospital beds nurses patients}{p_end}
{phang2}{cmd:. "Hospital A" 120 55 5100}{p_end}
{phang2}{cmd:. "Hospital B" 150 70 6200}{p_end}
{phang2}{cmd:. "Hospital C" 180 95 6400}{p_end}
{phang2}{cmd:. "Hospital D" 110 50 5000}{p_end}
{phang2}{cmd:. "Hospital E" 210 120 6600}{p_end}
{phang2}{cmd:. "Hospital F" 160 75 6100}{p_end}
{phang2}{cmd:. end}{p_end}

{pstd}Compute the historical VRS equality-restriction ratio{p_end}
{phang2}{cmd:. dea beds nurses = patients, model(congestion) gen(cong) dmu(hospital)}{p_end}
{phang2}{cmd:. list hospital cong_vrs cong_congestion if cong_congestion < 0.99}{p_end}

    {title:Super-efficiency with reference sets}

{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input str7 dmu i_employee i_area o_sales o_profits}{p_end}
{phang2}{cmd:. "store A" 10 20 70 6}{p_end}
{phang2}{cmd:. "store B" 15 15 100 3}{p_end}
{phang2}{cmd:. "store C" 20 30 80 5}{p_end}
{phang2}{cmd:. "store D" 25 15 100 2}{p_end}
{phang2}{cmd:. "store E" 12 9 90 8}{p_end}
{phang2}{cmd:. "store F" 16 18 75 5}{p_end}
{phang2}{cmd:. end}{p_end}

{pstd}Radial super-efficiency with reference set display{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(supereff) rts(vrs) ort(i) gen(se_vrs) reference dmu(dmu)}{p_end}
{pstd}In this VRS leave-one-out example, stores B and E cannot be dominated
by a convex combination of the remaining stores. Their scores are missing
and the command reports two failed DMUs; this illustrates infeasibility,
not a complete ranking. Inspect missing scores before interpreting results.{p_end}
{phang2}{cmd:. assert missing(se_vrs) if inlist(dmu, "store B", "store E")}{p_end}
{phang2}{cmd:. assert !missing(se_vrs) if !inlist(dmu, "store B", "store E")}{p_end}


{pstd}Save super-efficiency reference weights as variables{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(supereff) rts(crs) ort(i) gen(se_crs) reference saveref dmu(dmu)}{p_end}

    {title:Bootstrap bias correction}

{pstd}Use a fixed seed for this small demonstration. Increase replications
for substantive inference. See {help dea_sw2007} for a complete two-stage
example with environmental variables.{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. set seed 613}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(vrs) ort(i) gen(te_bc) bootstrap reps1(100) dmu(dmu)}{p_end}
{phang2}{cmd:. list dmu te_bc te_bc_bc te_bc_se te_bc_lo te_bc_hi}{p_end}
{phang2}{cmd:. assert te_bc_lo<=te_bc_hi}{p_end}

    {title:Engine and syntax options}

{pstd}Using the Simplex engine{p_end}
{pstd}Force the package-supplied Two-Phase Revised Simplex engine instead of
automatic engine selection.{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) ort(i) engine(simplex) gen(theta_s) dmu(dmu)}{p_end}

{pstd}Using {cmd:method()} synonym{p_end}
{pstd}{cmd:method(simplex)} is accepted as a synonym for {cmd:engine(simplex)}
for backward compatibility.{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) ort(i) method(simplex) gen(theta_m) dmu(dmu)}{p_end}

{pstd}Extended syntax form (Form B){p_end}
{phang2}{cmd:. dea o_sales o_profits, inputs(i_employee i_area) rts(vrs) gen(te_vrs2) dmu(dmu) replacevars}{p_end}

{pstd}Single-stage (no slack maximization){p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) stage(1) engine(simplex) gen(theta_s1) dmu(dmu) reference}{p_end}
{phang2}{cmd:. matrix list r(stage_linkage)}{p_end}
{phang2}{cmd:. display r(n_stage_link_adjusted)}{p_end}
{phang2}{cmd:. matrix stage1_refs=r(lambda)}{p_end}
{phang2}{cmd:. capture matrix list r(slack_in)}{p_end}
{phang2}{cmd:. assert _rc==111}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) stage(2) engine(simplex) gen(theta_s2) dmu(dmu) reference}{p_end}
{phang2}{cmd:. assert theta_s1==theta_s2}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) stage(2) engine(ipm) gen(theta_ipm) dmu(dmu) reference}{p_end}
{phang2}{cmd:. matrix list r(stage_status)}{p_end}
{phang2}{cmd:. display r(n_radial_failed) " radial failures; " r(n_stage2_failed) " unavailable stage-2 projections"}{p_end}
{phang2}{cmd:. matrix list r(slack_in)}{p_end}
{phang2}{cmd:. assert "`:colnames r(slack_in)'"=="i_employee i_area"}{p_end}
{phang2}{cmd:. assert "`:colnames r(slack_out)'"=="o_sales o_profits"}{p_end}

{pstd}Explicit recovery of an unavailable IPM Stage-2 projection{p_end}
{pstd}The original efficiency score is preserved. A successful original
projection needs no adjustment; inspect the reported targets and statuses.{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) stage(2) engine(ipm) radialtolerance(1e-7) gen(theta_recover) dmu(dmu) reference}{p_end}
{phang2}{cmd:. assert theta_recover==theta_ipm}{p_end}
{phang2}{cmd:. matrix list r(projection)}{p_end}
{phang2}{cmd:. list dmu theta_recover _dea_radial _dea_radial2 _dea_radial_adj}{p_end}

{pstd}Save results to a temporary dataset{p_end}
{phang2}{cmd:. tempfile dea_results}{p_end}

{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(vrs) gen(te) saving("`dea_results'") replace dmu(dmu) replacevars}{p_end}



{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{title:Simplex numerical checks}

{pstd}
The package simplex solver checks both inverse identities when accepting a
new basis inverse. A conservative rounding-error bound permits a fast check
when there is sufficient margin; other cases use compensated residuals.
The inverse check retains its internal 1e-8 threshold and does not change
{cmd:tolerance()} or the original-LP acceptance checks.

{pstd}
On difficult data, a result may become available or be rejected when an older
inverse check was inaccurate. Different reference weights may represent
alternative optimal solutions. These checks do not guarantee a solution for
every numerical problem or change Stata's native IPM implementation.
VRS super-efficiency can remain mathematically infeasible.

{pstd}
After the established simplex paths fail numerically, the solver can retry an
artificial-free basis and refine its solution using the original LP coefficients.
This candidate must pass the original primal checks and additional quality checks
that account for arithmetic rounding error. Previously accepted results proceed to the successful-basis check below.
The extra checks have finite dimension and exponent limits; an uncertain candidate
remains unavailable. A failed transition can propose at most two real-column bases while preserving
all original rows. Column selection uses bounded QR workspace for small problems
and a thin Householder method for wider ones; original feasibility and quality
checks still decide acceptance. The precision screen requires {it:n+2m<=4096}
for the standard-form problem and retains its finite arithmetic range.
A detected repeat of a Phase-2 Bland pivot state ends that numerical path early.
Existing restarts retain priority. If every preceding path is still rejected,
the solver may retry with price intervals that include basis-solve error.
An uncertain price is not treated as proof of an improving variable or an optimum.
This bounded repeat check does not detect every possible cycle.
If all earlier paths still fail after a blocked Phase-1 pivot, a final numerical
retry can construct a feasible basis from the original columns. It preserves
the original constraint rows and requires the existing original-unit primal/dual
quality checks. Earlier accepted scores and reference weights keep priority.
All repeated work is included in the iteration total; a failed retry can add work.
This recovery does not make {cmd:tolerance()} a guarantee of
relative objective error or certify exact optimality for all LPs. No new option is required.

{pstd}
Successful simplex results receive a basis-quality check. Already-certified
values are retained. A verified defect may be corrected only when the refined
point passes the original primal and conservative primal-dual checks. If
certification is unavailable or refinement fails, the original success is
retained; this does not certify every returned optimum. The public LP adapter
checks original coefficients and bounds, including lower-bound objective
constants. Exact zero rows are handled separately. A conservative fast check
falls back to compensated arithmetic when necessary. No new public option is
required, and native IPM calculations are unchanged. Earlier-result priority
in the failure-only descriptions is subject to this additional check.

{pstd}
The simplex engine also supports bounded exact selection of independent
original rows, including the right-hand side, and analytic checks for a trivial
empty basis. After every earlier recovery path fails, it may solve a reduced
row system once; acceptance still requires checks against all original rows,
columns and bounds. Near dependence alone never authorizes row removal.
The selector is limited to256 rows/columns,rank32,32 expansion components and
one million work units, with additional arithmetic-range guards. Unsupported
cases retain the existing outcome. No new public option or default strict
policy is introduced. These finite checks do not guarantee exact optimality
for every dataset. The packaged ROW_RECOVERY.txt note describes the regression scope.

{marker stagelink}{...}
{title:Simplex stage-link diagnostics}

{pstd}
Ordinary simplex TE now checks the radial factor passed from Stage 1 to
Stage 2. A bounded retry is permitted only after arithmetic error intervals
prove a negative component in the original-row basis solution. A corrected
factor must pass the change budget and original primal/dual quality checks.
The retry uses min(tolerance(),1e-10) for that DMU; the public default tolerance
is unchanged. Stage 2 rebuilds its LP from the same original-data-derived,
column-normalized inputs and outputs, with a fresh basis and the accepted
factor. A failed corrected Stage 2 retains the score and withholds references.
Stage(1) uses the same radial decision and skips slack optimization.

{pstd}
r(stage_linkage) records the before/after factors, reason, retry tolerance and
Stage-2 quality for the main ordinary simplex TE calculation. It is available
without recoverylog or reference; VRS auxiliary calls do not replace its rows.
r(numerical_diagnostics)'s radial_initial is the accepted factor used to build
Stage 2, so use stage_linkage's radial_before for the value before this retry.
Factors are theta (input) and phi (output); output efficiency is 1/phi, subject
to the existing near-one presentation cleanup. The matrix is not added to
saving() files. The IPM projection-recovery contract remains separate.

{pstd}
The basis check is bounded by m<=32,n+2m<=4096,m*n<=16384 and arithmetic-range
guards. Here m and n are the original Stage-1 LP row and real-variable counts.
It requires a complete original-row basis; unsupported/uncertain cases retain
the prior result. The dual fallback is not covered by this correction. Shared radial TE
subproblems also use the guard, but only ordinary model(te) exposes this
public diagnostic matrix.
One fresh retry is allowed; existing solver-internal recovery caps still apply.
This does not guarantee a feasible exact fixed target for every dataset, nor
remove all tolerance sensitivity or nonuniqueness of reference weights.

{pstd}
{cmd:r(stage_linkage)} has N rows in the same order as {cmd:r(efficiency)}
and the following 30 columns. It is absent for IPM and other public models.
Missing certificate fields mean not evaluated or unavailable, never a pass.

{p2colset 5 37 39 2}{...}
{p2col:{cmd:obs output stage policy}}original observation number; output=1/input=0; requested stage; internal policy (normally 2){p_end}
{p2col:{cmd:gate old_basis negative}}residual trigger passed; old basis certified; number of proven negative basic components{p_end}
{p2col:{cmd:old_upper old_error}}minimum old component upper bound; absolute basis-solve error bound{p_end}
{p2col:{cmd:retry_tol retry_status}}requested retry tolerance and retry solver status; missing if not attempted{p_end}
{p2col:{cmd:new_basis new_lower new_error}}new basis certified; minimum lower bound; absolute solve-error bound{p_end}
{p2col:{cmd:applied}}1 if the bounded radial repair was adopted, otherwise 0{p_end}
{p2col:{cmd:radial_before radial_after}}factors before/after the guard; equal when unchanged{p_end}
{p2col:{cmd:relative_change}}abs(after-before)/(1+abs(before)); bounded by tolerance() for an adopted repair{p_end}
{p2col:{cmd:residual_lower residual_upper}}old Stage-1 equality residual interval, maximum over rows{p_end}
{p2col:{cmd:residual_after safe_radial}}trial residual upper bound and trial safe factor; a trial may be rejected{p_end}
{p2col:{cmd:stage2_status reason}}final Stage-2 status and stage-link decision code{p_end}
{p2col:{cmd:old_iters retry_iters}}iterations in the first Stage-1 solve and additional retry, including their internal recovery attempts{p_end}
{p2col:{cmd:effective_tol rts}}accepted solve/extraction tolerance; RTS code CRS=1,VRS=2,DRS=3,IRS=4{p_end}
{p2col:{cmd:stage2_quality}}extra corrected-target check: 1 preserved,2 refined,0 failed,-1 unavailable; missing if not performed{p_end}
{p2col:{cmd:requested_tol}}the command's tolerance(){p_end}
{p2colreset}{...}

{pstd}
The stage-link residual denominator is 1+abs(b_i)+sum(abs(A_ij*x_j));
these interval bounds use stored binary64 normalized coefficients and differ
from the existing returned-projection diagnostic scale. The trigger is
64*(m+n+4)*machine epsilon. A small-residual screen retains an existing
result; it is not a feasibility or optimality certificate.

{pstd}
Reason codes: 0 off/unchanged (private testing); 1 unsupported range;
2 residual trigger not met; 3 old basis uncertified; 4 no proven negative
component; 5 audit-only policy (private testing); 6 retry failed; 7 new basis
uncertified or nonnegativity unproved; 8 radial/change budget failed;
9 primal/residual/quality check failed; 10 repair adopted;
11 fast small-residual screen; 12 no primal Stage-1 result to check
(dual fallback or both solves failed). Reason 10 describes Stage 1; inspect
stage2_status and stage2_quality separately for the returned reference.
{cmd:r(n_stage_link_adjusted)} and {cmd:r(n_stage_link_retries)} count only
main-analysis DMUs. {cmd:r(recovery_log)} additionally includes auxiliary
LP calls when enabled. Internal policies are not new public syntax options.

{phang2}{cmd:. matrix list r(stage_linkage)}{p_end}
{phang2}{cmd:. display r(n_stage_link_adjusted)}{p_end}

{marker recoverylog}{...}
{title:Simplex recovery log}

{pstd}
After all earlier simplex recovery paths fail, one bounded retry may remove
exact positive power-of-two duplicate rows or columns. Row comparisons include
the right-hand side; column comparisons include the objective and every original
constraint and bound row. The solution is restored to the original variables
and must pass the original primal and primal/dual quality checks. Earlier
accepted common-simplex solutions and reference weights retain priority. All attempted
iterations are included in the reported total, with one aggregate recovery
record per common simplex call.

{pstd}
This retry uses exact stored-coefficient comparisons, not tolerance-based row
deletion. Supported ratios are powers of two with exponents from -400 to 400;
finite-arithmetic checks can further restrict this range. The original
standard-form size limit is n+2m <= 4096, and duplicate comparisons have a
1,000,000-element budget. Unverified or unsupported candidates remain numerical
failures; the retry does not publish a new infeasible or unbounded certificate.
It does not guarantee recovery of every degenerate LP or improve every already
accepted result. User tolerances and native IPM calculations are unchanged.

{pstd}
With {cmd:engine(auto)}, a new simplex recovery can replace a previously
successful native fallback. Scores, reference weights, iteration totals and the
reported backend can therefore change across versions. Earlier-result priority
within the common simplex solver does not promise identical automatic-engine
results across versions.



{phang}
{cmd:recoverylog} records the package simplex LP calls made by this command.
It is optional and off by default. It supports ordinary {cmd:model(te)},
{cmd:model(sbm)}, {cmd:model(supersbm)}, and
{cmd:model(undesirable) disposability(weak)} with either command syntax.
It rejects bootstrap, IDEA/bounded/ordinal, FDH, {cmd:contractcheck}, and other
models. It does not change tolerances, solving policy, scores or references.

{pstd}
A row represents one common simplex call, including its optional factorized
retry, rather than one DMU. The initial eta solver retries only numerical or
iteration failure (status 3). The retry uses compensated, uncached basis
factorization and must pass the submitted standard-form LP residual check.
An infeasible or unbounded initial classification is not retried by this rule.
Later model checks can still reject a core success. This is numerical evidence,
not an exact feasibility or optimality certificate.

{pstd}
{cmd:r(recovery_stats)} is 1 x 8 with columns {cmd:lp_calls attempted recovered
eta_iters factor_iters total_iters omitted stored}. Counts cover every recorded
simplex call, including VRS auxiliary calls. The final two columns count omitted
and stored detail rows, so {it:lp_calls=omitted+stored}. Iteration totals add the
initial and retry iterations. At most {it:min(c(matsize),10000)} first rows are
stored; counts continue after that limit. {cmd:r(recovery_log)} is returned only
when at least one row is stored. A native-only run returns zero counts and no
log matrix; native attempts remain in their existing diagnostics.
{cmd:r(recovery_scope)} is {cmd:package_simplex_calls}.

{pstd}
The 25 columns of {cmd:r(recovery_log)} are:

{p2colset 8 31 33 2}{...}
{p2col:{cmd:obsno}}Original observation number in the current dataset, after applying
selection; not a numeric conversion of {cmd:dmu()}. Copy the log before sorting
or replacing the data if you need to join it to string or duplicate labels.{p_end}
{p2col:{cmd:scope}}1 = requested primary calculation; 2 = auxiliary VRS frontier
CRS/VRS/NIRS calculations.{p_end}
{p2col:{cmd:model}}1 = TE; 2 = ordinary SBM (including super-SBM classification);
3 = super-SBM; 4 = joint weak-disposal DDF.{p_end}
{p2col:{cmd:rts}}1 = CRS; 2 = VRS; 3 = DRS/NIRS; 4 = IRS/NDRS.{p_end}
{p2col:{cmd:stage}}TE: 1 = radial, 2 = slack. SBM and weak DDF: 1 = single stage.{p_end}
{p2col:{cmd:form}}TE: 1 = primal, 2 = dual. SBM: 1. Weak DDF: 1 = active/idle,
2 = total/active reformulation, 3 = dual bound, 4 = reversed-row dual bound.{p_end}
{p2col:{cmd:lp_seq}}Call sequence within this command.{p_end}
{p2col:{cmd:attempted recovered}}Factorized retry attempted and accepted by the common
solver (0/1). Acceptance is not a promise of a final score or reference.{p_end}
{p2col:{cmd:eta_status factor_status}}Initial and retry solver statuses: 0 = optimal,
1 = infeasible, 2 = unbounded, 3 = numerical/iteration failure. An unattempted
retry has missing status.{p_end}
{p2col:{cmd:eta_iters factor_iters total_iters}}Initial, retry and combined iteration
counts; an unattempted retry has zero iterations.{p_end}
{p2col:{cmd:core_status}}Status returned by the common solver after recovery policy.{p_end}
{p2col:{cmd:eta_reason factor_reason}}Numerical reason codes described below;
unattempted retry is missing.{p_end}
{p2col:{cmd:n_rows n_cols}}Rows and original-variable columns submitted to the common
standard-form solver, before redundant-row removal. These need not equal the
DEA input/output dimensions.{p_end}
{p2col:{cmd:caller_status}}Status observed by the calling adapter. Non-TE code 4
means its original submitted constraints or finite objective check rejected the
candidate. Further weak-DDF checks can still reject this result.{p_end}
{p2col:{cmd:model_score model_ref}}Availability in the matching observation/scope/model/RTS
calculation. Score is 0/1. Reference is 0/1 for TE and missing for other models.
These group flags repeat across that calculation's LP rows. For TE/SBM they
describe the simplex model call before the main command's native fallback; for
weak DDF they include its internal native fallback and gap recovery.{p_end}
{p2col:{cmd:final_available final_score final_source}}Main command's final score for
this observation: available (0/1), its value, and source (0 = unavailable,
1 = package simplex, 2 = native, 3 = weak-DDF feasible-self/dual-gap recovery).{p_end}
{p2colreset}{...}

{pstd}
Reason codes are 0 = no mapped numerical reason; 1 = invalid initial basis;
2 = invalid basic solution; 3 = price computation unavailable; 4 = pivot
unavailable; 5 = Phase-1 transition inverse; 6 = artificial-basis pivot;
7 = artificial-basis transition; 8 = Phase-2 basis inverse;
9 = residual check; 10 = iteration limit; 11 = lost feasibility;
12 = refactorization. These describe detected failure points, not a proven
underlying cause.

{pstd}
Final score fields repeat the main command result even on auxiliary scope-2
rows or ordinary-SBM classification rows. They do not say that a particular LP
provided the chosen score. For example, a failed simplex row can coexist with
a final native score after {cmd:engine(auto)}. A successful TE radial score can
coexist with unavailable stage-2 references. Use {cmd:r(stage_status)},
{cmd:r(solution_source)} and the existing 23-column
{cmd:r(numerical_diagnostics)} for those TE contracts. Their columns are unchanged.
A subsequent {cmd:dea} command resets the recording; a nonzero Stata command
return code aborts and discards it. A completed command can retain missing
scores and their failed LP records.
The option does not add diagnostic variables to {cmd:saving()} datasets.


{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of DMUs{p_end}
{synopt:{cmd:r(N_eff)}}number of efficient DMUs, based on the main score{p_end}
{synopt:{cmd:r(n_failed)}}selected rows with missing results; TE also counts unavailable requested Stage 2; congestion counts unavailable requested components or ratios once per row; IDEA counts either main or {cmd:_opt} missing, once per row{p_end}
{synopt:{cmd:r(n_radial_failed)}}missing radial scores, for {cmd:model(te)}{p_end}
{synopt:{cmd:r(n_stage2_failed)}}unavailable requested Stage-2 projections, for {cmd:model(te)}; zero with {cmd:stage(1)}{p_end}
{synopt:{cmd:r(N_eff_opt)}}IDEA {cmd:_opt} scores within tolerance of 1{p_end}
{synopt:{cmd:r(mean_opt)}}IDEA mean of nonmissing {cmd:_opt} scores{p_end}
{synopt:{cmd:r(min_opt)}}IDEA minimum {cmd:_opt} score{p_end}
{synopt:{cmd:r(max_opt)}}IDEA maximum {cmd:_opt} score{p_end}
{synopt:{cmd:r(n_failed_primary)}}IDEA rows with missing main score{p_end}
{synopt:{cmd:r(n_failed_opt)}}IDEA rows with missing {cmd:_opt} score{p_end}
{synopt:{cmd:r(N_super)}}number with score > 1+tolerance; returned for {cmd:supereff} and {cmd:supersbm}{p_end}
{synopt:{cmd:r(mean_eff)}}mean stored score (inefficiency for DDF/additive){p_end}
{synopt:{cmd:r(min_eff)}}minimum stored score{p_end}
{synopt:{cmd:r(max_eff)}}maximum stored score{p_end}
{synopt:{cmd:r(stage)}}number of optimization stages{p_end}
{synopt:{cmd:r(tolerance)}}numerical tolerance used{p_end}
{synopt:{cmd:r(defaults_applied)}}1 when the call supplied no options and all public defaults were applied{p_end}
{synopt:{cmd:r(generated_dmu)}}1 when {cmd:dea} generated {cmd:_dmu} from observation order{p_end}
{synopt:{cmd:r(auto_scale_ratio)}}dimensional scale ratio inspected by {cmd:engine(auto)}{p_end}
{synopt:{cmd:r(auto_scale_threshold)}}scale-ratio threshold ({cmd:1e8}){p_end}
{synopt:{cmd:r(auto_preemptive_ipm)}}1 if scale protection selected IPM before solving{p_end}
{synopt:{cmd:r(n_crs)}}number classified as CRS; returned after {cmd:model(te) rts(vrs)}{p_end}
{synopt:{cmd:r(n_drs)}}number classified as DRS; returned after {cmd:model(te) rts(vrs)}{p_end}
{synopt:{cmd:r(n_irs)}}number classified as IRS; returned after {cmd:model(te) rts(vrs)}{p_end}
{synopt:{cmd:r(bootstrap_reps)}}replications after {cmd:bootstrap}{p_end}
{synopt:{cmd:r(bootstrap_level)}}confidence level after {cmd:bootstrap}{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}}{cmd:dea}{p_end}
{synopt:{cmd:r(cmdline)}}full command as typed{p_end}
{synopt:{cmd:r(model)}}effective model actually dispatched{p_end}
{synopt:{cmd:r(model_requested)}}normalized model option; {cmd:te} when omitted{p_end}
{synopt:{cmd:r(idea_method)}}{cmd:interval} or {cmd:legacy}, after IDEA{p_end}
{synopt:{cmd:r(idea_numeric_vars)}}numeric bounded dimensions in the legacy mixed path; empty otherwise{p_end}
{synopt:{cmd:r(idea_vars)}}all generated IDEA score-variable names{p_end}
{synopt:{cmd:r(rts)}}returns to scale{p_end}
{synopt:{cmd:r(orientation)}}orientation; {cmd:non-oriented} for SBM/super-SBM when omitted{p_end}
{synopt:{cmd:r(engine)}}engine selection ({cmd:ipm} or {cmd:simplex}), after retries;
FDH and interval IDEA use their dedicated backends independently of this selector{p_end}
{synopt:{cmd:r(inputs)}}input variable names{p_end}
{synopt:{cmd:r(outputs)}}output variable names{p_end}
{synopt:{cmd:r(generate)}}name of generated efficiency variable{p_end}
{synopt:{cmd:r(refvars)}}persistent reference variables in lambda-column order when {cmd:saveref} is used{p_end}

{synopt:{cmd:r(bootstrap_engine)}}{cmd:simplex}, after {cmd:bootstrap}{p_end}
{synopt:{cmd:r(bootstrap_vars)}}five additional result-variable names, after {cmd:bootstrap}{p_end}

{synopt:{cmd:r(radial_tolerance)}}declared recovery allowance, default 0{p_end}
{synopt:{cmd:r(n_projection_attempted)}}changed-target Stage-2 LPs attempted, positive allowance only{p_end}
{synopt:{cmd:r(n_projection_adjusted)}}successful changed-target projections{p_end}
{synopt:{cmd:r(n_projection_solves)}}native calls for changed-target LPs, including their retries{p_end}
{synopt:{cmd:r(projection_vars)}}three persistent radial metadata variables, positive allowance only{p_end}
{synopt:{cmd:r(n_ipm_retries)}}extra native solves in main TE phases{p_end}
{synopt:{cmd:r(n_ipm_recovered)}}main TE LP phases recovered by native retries{p_end}
{synopt:{cmd:r(n_ipm_conflicts)}}main TE LP phases with final status 91{p_end}
{synopt:{cmd:r(n_vrs_ipm_retries)}}extra native solves in auxiliary VRS analysis{p_end}
{synopt:{cmd:r(n_vrs_ipm_recovered)}}auxiliary LP phases recovered by native retries{p_end}
{synopt:{cmd:r(n_vrs_ipm_conflicts)}}auxiliary LP phases with final status 91{p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:r(efficiency)}}efficiency scores ({it:N} x 1){p_end}
{synopt:{cmd:r(bootstrap)}}{it:N} x 5: bias_corrected, bias, se, ci_lo, ci_hi;
rows follow the same selected-observation order as {cmd:r(efficiency)}{p_end}
{synopt:{cmd:r(lambda)}}reference weights ({it:N} x {it:N}); returned when
{cmd:reference} requests them and the selected model supports them{p_end}
{synopt:{cmd:r(slack_in)}}input slacks ({it:N} x {it:m}); returned for
technical-efficiency calls with {cmd:stage(2)}, with either engine;
columns have the expanded input variable names, in the specified order{p_end}
{synopt:{cmd:r(slack_out)}}output slacks ({it:N} x {it:s}); returned for
technical-efficiency calls with {cmd:stage(2)}, with either engine;
columns have the expanded output variable names, in the specified order{p_end}
{pstd}From version 2.0.51, these slack column names identify the model variables
instead of internal temporary variables. Rows retain the selected-observation
order of {cmd:r(efficiency)}. The names are preserved after automatic engine
fallback, VRS decomposition, bootstrap and {cmd:saving()}. Slack values remain
in the original units; {cmd:saving()} does not add slack variables to the saved
dataset. {cmd:stage(1)} does not return these matrices.{p_end}
{synopt:{cmd:r(stage_status)}}TE, either engine: {it:N} x 2, columns {cmd:radial} and
{cmd:stage2}; same row order as {cmd:r(efficiency)}. Codes described under
{cmd:stage()}; missing Stage-2 status means not attempted{p_end}
{synopt:{cmd:r(solution_source)}}TE, either engine: {it:N} x 2, columns {cmd:radial}
and {cmd:reference}; score and reference provenance, codes under {cmd:stage()}{p_end}
{synopt:{cmd:r(recovery_stats)}}with {cmd:recoverylog}: 1 x 8 simplex call/retry/iteration and row-limit counts{p_end}
{synopt:{cmd:r(recovery_log)}}with {cmd:recoverylog} and stored calls: LP x 25 context, solver and final-result fields{p_end}
{synopt:{cmd:r(recovery_scope)}}with {cmd:recoverylog}: {cmd:package_simplex_calls}{p_end}
{synopt:{cmd:r(stage_linkage)}}ordinary TE using simplex: N x 30 before/after radial and retry diagnostics; see {help dea##stagelink:stage-link diagnostics}{p_end}
{synopt:{cmd:r(n_stage_link_adjusted)}}main DMUs with an adopted radial repair{p_end}
{synopt:{cmd:r(n_stage_link_retries)}}main DMUs with a stage-link retry{p_end}
{synopt:{cmd:r(numerical_diagnostics)}}TE, either engine: {it:N} x 23, radial factors and returned-projection residuals; see {cmd:stage()}{p_end}
{synopt:{cmd:r(n_projection_checked)}}TE: number of Stage-2 projections checked{p_end}
{synopt:{cmd:r(n_projection_outside_tol)}}TE: checked projections above the scaled tolerance{p_end}
{synopt:{cmd:r(n_projection_unchecked)}}TE: available Stage-2 results not checkable{p_end}
{synopt:{cmd:r(projection)}}positive allowance only: {it:N} x 8 targets, adjustments and initial/retry/final statuses{p_end}
{synopt:{cmd:r(projection_trace)}}positive allowance only: {it:N} x 25 changed-target native attempt history{p_end}
{synopt:{cmd:r(ipm_trace)}}IPM TE only: {it:N} x 75, per-phase native attempt evidence;
see {cmd:stage()}{p_end}
{synopt:{cmd:r(vrs_ipm_trace)}}IPM VRS decomposition only: {it:N} x 150,
auxiliary primal/dual attempt evidence{p_end}
{synopt:{cmd:r(vrsfrontier)}}CRS TE, VRS TE, NIRS TE (legacy column name
{cmd:DRS_TE}), scale efficiency, and RTS code
after {cmd:model(te) rts(vrs)} ({it:N} x 5){p_end}
{synopt:{cmd:r(vrsobs)}}source observation numbers ({it:N} x 1, column
{cmd:obsno}) corresponding to {cmd:r(vrsfrontier)} rows; returned with the
VRS decomposition. Both matrices use selected-observation order and row names
{cmd:r1}, ..., {cmd:rN}{p_end}
{synopt:{cmd:r(multiplier_u0)}}signed multiplier intercept ({it:N} x 1);
returned after {cmd:model(multiplier) stage(2) saveweights}{p_end}


{marker references}{...}
{title:References}

{pstd}
{it:DEA Foundational Models}

{marker F1957}{...}
{phang}
Farrell, M. J. 1957.
The measurement of productive efficiency.
{it:Journal of the Royal Statistical Society, Series A} 120: 253-290.

{marker CCR1978}{...}
{phang}
Charnes, A., W. W. Cooper, and E. Rhodes. 1978.
Measuring the efficiency of decision making units.
{it:European Journal of Operational Research} 2: 429-444.

{marker BCC1984}{...}
{phang}
Banker, R. D., A. Charnes, and W. W. Cooper. 1984.
Some models for estimating technical and scale inefficiencies in data
envelopment analysis.  {it:Management Science} 30: 1078-1092.

{marker CCS2007}{...}
{phang}
Cooper, W. W., L. M. Seiford, and K. Tone. 2007.
{it:Data Envelopment Analysis: A Comprehensive Text with Models, Applications,}
{it:References and DEA-Solver Software}. 2nd ed. New York: Springer.

{pstd}
{it:Extended Models}

{marker CCD1985}{...}
{phang}
Charnes, A., W. W. Cooper, B. Golany, L. Seiford, and J. Stutz. 1985.
Foundations of data envelopment analysis for Pareto-Koopmans efficient
empirical production functions.
{it:Journal of Econometrics} 30: 91-107.

{marker AP1993}{...}
{phang}
Andersen, P., and N. C. Petersen. 1993.
A procedure for ranking efficient units in data envelopment analysis.
{it:Management Science} 39: 1261-1264.

{marker CCF1996}{...}
{phang}
Chambers, R. G., Y. Chung, and R. Fare. 1996.
Benefit and distance functions.
{it:Journal of Economic Theory} 70: 407-419.

{marker CRS1998}{...}
{phang}
Chambers, R. G., Y. Chung, and R. Fare. 1998.
Profit, directional distance functions, and Nerlovian efficiency.
{it:Journal of Optimization Theory and Applications} 98: 351-364.

{marker T2001}{...}
{phang}
Tone, K. 2001.
A slacks-based measure of efficiency in data envelopment analysis.
{it:European Journal of Operational Research} 130: 498-509.

{marker T2002}{...}
{phang}
Tone, K. 2002.
A slacks-based measure of super-efficiency in data envelopment analysis.
{it:European Journal of Operational Research} 143: 32-41.

{marker T2004}{...}
{phang}
Tone, K. 2004.
Dealing with undesirable outputs in DEA: A Slacks-Based Measure (SBM) approach.
Presentation at NAPW III, Toronto, pp. 44--45.

{marker CDSZ2004}{...}
{phang}
Cooper, W. W., H. Deng, L. M. Seiford, and J. Zhu. 2004.
Congestion: Its identification and management with DEA.
In {it:Handbook on Data Envelopment Analysis},
ed. W. W. Cooper, L. M. Seiford, and J. Zhu, 177--201.
Boston: Kluwer Academic Publishers.

{marker DPT2006}{...}
{phang}
Deprins, D., L. Simar, and H. Tulkens. 2006.
Measuring labor-efficiency in post offices.
In {it:Public Goods, Environmental Externalities and Fiscal Competition},
ed. P. Chander, J. Drèze, C. K. Lovell, and J. Mintz.
Boston, MA: Springer.
{browse "https://doi.org/10.1007/978-0-387-25534-7_16":doi:10.1007/978-0-387-25534-7_16}

{marker CPY1999}{...}
{phang}
Cooper, W. W., K. S. Park, and G. Yu. 1999.
IDEA and AR-IDEA: Models for dealing with imprecise data in DEA.
{it:Management Science} 45(4): 597--607.
{browse "https://doi.org/10.1287/mnsc.45.4.597":doi:10.1287/mnsc.45.4.597}

{marker Z2003}{...}
{phang}
Zhu, J. 2003.
Imprecise data envelopment analysis (IDEA): A review and improvement
with an application.
{it:European Journal of Operational Research} 144: 513-529.

{marker SMD2006}{...}
{phang}
Smirlis, Y. G., E. K. Maragos, and D. K. Despotis. 2006.
Data envelopment analysis with missing values: An interval DEA approach.
{it:Applied Mathematics and Computation} 177: 1--10.
{browse "https://doi.org/10.1016/j.amc.2005.10.028":doi:10.1016/j.amc.2005.10.028}

{pstd}
{it:Decomposition and Scale}

{marker FL1978}{...}
{phang}
Fare, R., and C. A. K. Lovell. 1978.
Measuring the technical efficiency of production.
{it:Journal of Economic Theory} 19: 150-162.

{marker FGL1993}{...}
{phang}
Fare, R., S. Grosskopf, and C. A. K. Lovell. 1993.
{it:Production Frontiers}. Cambridge: Cambridge University Press.

{marker FGNZ1994}{...}
{phang}
Fare, R., S. Grosskopf, M. Norris, and Z. Zhang. 1994.
Productivity growth, technical progress, and efficiency change in
industrialized countries.
{it:American Economic Review} 84: 66-83.

{marker CC1962}{...}
{phang}
Charnes, A., and W. W. Cooper. 1962.
Programming with linear fractional functionals.
{it:Naval Research Logistics Quarterly} 9: 181-186.

{pstd}
{it:Statistical Inference}

{marker SW1998}{...}
{phang}
Simar, L., and P. W. Wilson. 1998.
Sensitivity analysis of efficiency scores: How to bootstrap in
nonparametric frontier models.
{it:Management Science} 44: 49-61.

{marker SW2007}{...}
{phang}
Simar, L., and P. W. Wilson. 2007.
Estimation and inference in two-stage, semi-parametric models of
production processes.
{it:Journal of Econometrics} 136: 31-64.

{pstd}
{it:Original Software}

{marker JL2010}{...}
{phang}
Ji, Y.-B., and C. Lee. 2010.
Data envelopment analysis.
{it:Stata Journal} 10(2): 267-280.
{browse "https://doi.org/10.1177/1536867X1001000207"}

{marker LJ2009}{...}
{phang}
Lee, C., and Y.-B. Ji. 2009.
Data envelopment analysis in Stata.
{it:Stata Conference DC 2009 Proceedings}.
{browse "https://ideas.repec.org/p/boc/dcon09/4.html"}


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
