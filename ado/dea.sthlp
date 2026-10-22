{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 2.0.28  05aug2026}{...}
{vieweralsosee "[R] frontier" "help frontier"}{...}
{vieweralsosee "[R] teradial" "help teradial"}{...}
{viewerjumpto "Syntax" "dea##syntax"}{...}
{viewerjumpto "Description" "dea##description"}{...}
{viewerjumpto "Options" "dea##options"}{...}
{viewerjumpto "Remarks" "dea##remarks"}{...}
{viewerjumpto "Examples" "dea##examples"}{...}
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
{synopt:{cmd:dmu(}{it:varname}{cmd:)}}DMU identifier variable{p_end}
{synopt:{cmdab:ref:erence}}display reference sets{p_end}
{synopt:{cmdab:saver:ef}}save reference weights as variables{p_end}
{synopt:{cmdab:savew:eights}}save multiplier weights{p_end}
{synopt:{cmdab:sav:ing(}{it:filename}{cmd:)}}save detailed results to dataset{p_end}
{synopt:{cmd:replace}}overwrite existing dataset when using {cmd:saving()}{p_end}

{syntab:Bootstrap (Simar-Wilson)}
{synopt:{cmdab:boo:tstrap}}perform Simar-Wilson bootstrap{p_end}
{synopt:{cmdab:env:ars(}{it:varlist}{cmd:)}}environmental variables for second-stage regression{p_end}
{synopt:{cmd:reps1(}{it:#}{cmd:)}}first-stage bootstrap replications; default is {cmd:100}{p_end}
{synopt:{cmd:reps2(}{it:#}{cmd:)}}second-stage bootstrap replications; default is {cmd:2000}{p_end}
{synopt:{cmd:level(}{it:#}{cmd:)}}confidence level (%); default is {cmd:95}{p_end}

{syntab:Technical}
{synopt:{cmdab:eng:ine(}{it:enginename}{cmd:)}}solver engine; {cmd:simplex} (default), {cmd:ipm}, or {cmd:auto}{p_end}
{synopt:{cmdab:met:hod(}{it:enginename}{cmd:)}}synonym for {cmd:engine()}{p_end}
{synopt:{cmdab:tol:erance(}{it:#}{cmd:)}}numerical tolerance; default is {cmd:1e-7}{p_end}
{synopt:{cmd:trace}}display iteration details{p_end}
{synopt:{cmd:debug}}display debugging information{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:dea} computes Data Envelopment Analysis (DEA) efficiency scores for a set
of decision-making units (DMUs) that use multiple inputs to produce multiple
outputs.  Version 2.0 provides sixteen models within a unified interface
and two LP solver engines with full model parity.

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
direction is the observed values.{p_end}

{p2col :{cmd:sbm}}Slacks-based measure ({help dea##T2001:Tone 2001}).
Non-radial efficiency incorporating all slacks, solved via
{help dea##CC1962:Charnes-Cooper (1962)} transformation.{p_end}

{p2col :{cmd:additive}}Additive model
({help dea##CCD1985:Charnes et al. 1985}).
Maximizes the sum of input and output slacks.  A DMU is efficient if
and only if all optimal slacks are zero.{p_end}

{p2col :{cmd:supereff}}Radial super-efficiency
({help dea##AP1993:Andersen and Petersen 1993}).  Excludes the
evaluated DMU from its own reference set, allowing scores above
unity.  Use {cmd:reference} or {cmd:saveref} to display or store the
reference sets (lambda weights).{p_end}

{p2col :{cmd:supersbm}}Super-efficiency SBM ({help dea##T2002:Tone 2002}).
Super-efficiency extension of the slacks-based measure.  Efficient DMUs
(SBM score = 1) receive a super-efficiency score via Tone's [LP]
formulation with self-exclusion from the reference set.  Inefficient
DMUs retain their standard SBM score.{p_end}

{p2col :{cmd:allocative}}Allocative efficiency.  Decomposes economic efficiency
into technical and allocative components following
{help dea##F1957:Farrell (1957)} and
{help dea##FL1978:Fare and Lovell (1978)}.  Requires {cmd:w()} for input-oriented
or {cmd:p()} for output-oriented analysis.{p_end}

{p2col :{cmd:fdh}}Free disposal hull
({help dea##DPT2006:Deprins, Simar, and Tulkens 2006}).
Replaces convexity with integer constraints on lambda, producing a
staircase frontier.  Computed by enumeration without LP.{p_end}

{p2col :{cmd:imprecise}}Imprecise DEA (IDEA)
({help dea##CPY1999:Cooper, Park, and Yu 1999}).  Handles interval,
ordinal, and bounded data.  Computes pessimistic and optimistic
efficiency bounds.  See {it:Imprecise DEA options} below.{p_end}

{p2col :{cmd:undesirable}}Directional distance function with undesirable
(bad) outputs.  Extends {cmd:model(ddf)} by incorporating bad outputs
that are contracted while good outputs are expanded.  Requires
{cmd:bad(}{it:varlist}{cmd:)}.  Direction controlled by {cmd:direction()}:
{cmd:yb} (default, contract bad outputs only),
{cmd:xyb} (contract inputs and bad outputs while expanding good outputs),
or {cmd:x} (contract inputs only).{p_end}

{p2col :{cmd:congestion}}Congestion DEA
({help dea##CDSZ2004:Cooper et al. 2004}).
Detects input congestion by comparing VRS efficiency under strong versus
weak output disposability.  Reports the congestion measure
({it:theta_VRS / theta_C}; below 1 indicates congestion), plus CRS and
VRS efficiency and scale efficiency.  Additional variables are created
with suffixes {cmd:_crs}, {cmd:_vrs}, {cmd:_cong_eff}, {cmd:_congestion},
and {cmd:_se}.  Always uses VRS; the {cmd:rts()} option is ignored.{p_end}

{p2col :{cmd:sbmundesirable}}SBM with undesirable outputs
({help dea##T2004:Tone 2004}).  Non-radial efficiency incorporating
slacks for inputs, good outputs, and bad outputs simultaneously.
Bad outputs are treated under weak disposability (equality constraint).
Requires {cmd:bad(}{it:varlist}{cmd:)}.
See also {helpb dea_sbmundesirable}.{p_end}
{p2colreset}{...}

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
For {cmd:model(sbm)}, omitting the option requests the non-oriented SBM;
{cmd:orientation(input)} and {cmd:orientation(output)} request the corresponding
oriented SBM models.

{phang}
{cmd:stage(1|2)} specifies single-stage ({cmd:1}) or two-stage ({cmd:2},
the default) optimization.  In Stage 2, slacks are maximized while holding
the efficiency score fixed, ensuring a unique Pareto-Koopmans efficient
projection.

{phang}
{cmd:fdh} is shorthand for {cmd:model(fdh)}.  When specified, the
returns-to-scale is automatically set to VRS.

{dlgtab:Prices and directions}

{phang}
{cmd:w(}{it:varlist}{cmd:)} specifies input price variables for
{cmd:model(cost)}, {cmd:model(profit)}, and input-oriented {cmd:model(allocative)}.
The number of price variables must equal the number of inputs.

{phang}
{cmd:p(}{it:varlist}{cmd:)} specifies output price variables for
{cmd:model(revenue)}, {cmd:model(profit)}, and output-oriented {cmd:model(allocative)}.
The number of price variables must equal the number of outputs.

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
Bad outputs are treated under weak disposability --- they cannot be freely
reduced, reflecting the real cost of abatement.

{phang}
{cmd:direction(}{it:type}{cmd:)} specifies the contraction/expansion direction
for {cmd:model(undesirable)}.  Options are:

{p2colset 9 16 18 2}
{p2col:{cmd:yb}}contract bad outputs only (default){p_end}
{p2col:{cmd:xyb}}contract inputs and bad outputs simultaneously while expanding good outputs{p_end}
{p2col:{cmd:x}}contract inputs only{p_end}
{p2colreset}{...}

{dlgtab:Imprecise DEA (IDEA) options}

{pstd}
The following options activate imprecise DEA
({help dea##CPY1999:Cooper, Park, and Yu 1999};
{help dea##Z2003:Zhu 2003}).  Specifying any of these options
automatically invokes {cmd:model(imprecise)}.  IDEA computes both
pessimistic (worst-case) and optimistic (best-case) efficiency bounds.

{phang}
{cmd:xlo(}{it:varlist}{cmd:)} and {cmd:xhi(}{it:varlist}{cmd:)} specify
the lower and upper bounds for interval-valued inputs.  Each variable in
{cmd:xlo()} corresponds to the matching variable in {cmd:inputs()}, and
similarly for {cmd:xhi()}.  The main input variables serve as the center
(modal) values.

{phang}
{cmd:ylo(}{it:varlist}{cmd:)} and {cmd:yhi(}{it:varlist}{cmd:)} specify
the lower and upper bounds for interval-valued outputs.

{phang}
{cmd:xord(}{it:varlist}{cmd:)} and {cmd:yord(}{it:varlist}{cmd:)} specify
ordinal input and output variables.  Ordinal data preserve rank ordering
but not cardinal magnitude.

{phang}
{cmd:bounded(}{it:varlist}{cmd:)} specifies variables with known upper and
lower bounds.

{phang}
{cmd:ordinal(}{it:varlist}{cmd:)} specifies variables with ordinal
(rank-order) data.

{pstd}
IDEA generates three efficiency variables with the suffixes {cmd:_pess}
(pessimistic bound), {cmd:_opt} (optimistic bound), and {cmd:_mid}
(midpoint).  A DMU is {it:certainly efficient} if the pessimistic bound
equals 1, {it:possibly efficient} if the pessimistic bound is below 1
but the optimistic bound is at least 1, and {it:certainly inefficient}
if the optimistic bound is below 1.

{dlgtab:Output}

{phang}
{cmd:generate(}{it:newvar}{cmd:)} specifies the variable name for storing
efficiency scores.  The default name is {cmd:dea_efficiency}.  For {cmd:model(imprecise)},
three variables are created: {it:newvar}{cmd:_pess}, {it:newvar}{cmd:_opt},
and {it:newvar}{cmd:_mid}.

{phang}
{cmd:dmu(}{it:varname}{cmd:)} specifies a string or numeric variable
containing DMU identifiers for display purposes.

{phang}
{cmd:reference} displays the reference set (peer DMUs and their weights) for
each DMU in the results table.

{phang}
{cmd:saveref} saves reference weights (lambda values) as new variables named
{cmd:_ref_}{it:dmuname}.

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
{cmd:saving(}{it:filename}{cmd:)} saves detailed results including efficiency
scores, lambda weights, and slacks to a Stata dataset.

{phang}
{cmd:replace} permits {cmd:saving()} to overwrite an existing dataset.

{dlgtab:Bootstrap (Simar-Wilson)}

{pstd}
The following options implement the Simar-Wilson
({help dea##SW1998:1998}; {help dea##SW2007:2007}) bootstrap procedure for bias-corrected
efficiency scores and confidence intervals.

{phang}
{cmd:bootstrap} activates the Simar-Wilson two-stage bootstrap
procedure.  Bias-corrected efficiency scores and bootstrap confidence
intervals are computed.

{phang}
{cmd:envars(}{it:varlist}{cmd:)} specifies environmental (contextual)
variables for the second-stage truncated regression in the Simar-Wilson
procedure.

{phang}
{cmd:reps1(}{it:#}{cmd:)} specifies the number of first-stage bootstrap
replications used to estimate bias.  The default is {cmd:100}.

{phang}
{cmd:reps2(}{it:#}{cmd:)} specifies the number of second-stage bootstrap
replications for confidence intervals.  The default is {cmd:2000}.

{phang}
{cmd:level(}{it:#}{cmd:)} specifies the confidence level as a percentage.
The default is {cmd:95}.

{dlgtab:Technical}

{phang}
{cmd:engine(}{it:enginename}{cmd:)} specifies the LP solver engine.
{cmd:simplex} (the default) uses the package-supplied Two-Phase Revised
Simplex method implemented in Mata.  {cmd:ipm} uses Stata's {cmd:LinearProgram()}
interior point method, available since Stata 16.  {cmd:auto} normally tries
simplex first and re-solves all DMUs with IPM if any score is missing.  For
radial technical efficiency, it instead selects IPM before solving when the
ratio between the largest and smallest nonzero variable scales exceeds
{cmd:1e8}.  Explicit {cmd:engine(simplex)} is never overridden.  Both engines
support all sixteen models; unusually heterogeneous units may require IPM or
rescaling.

{phang}
{cmd:method(}{it:enginename}{cmd:)} is a synonym for {cmd:engine()}.
If both are specified, {cmd:engine()} takes precedence.

{phang}
{cmd:tolerance(}{it:#}{cmd:)} specifies the numerical tolerance for
convergence and efficiency classification.  The default is {cmd:1e-7}.
A DMU is classified as efficient when its score is within this tolerance
of the efficient frontier (1 for radial models, 0 for DDF).

{phang}
{cmd:trace} displays iteration-level details from the LP solver.

{phang}
{cmd:debug} displays additional debugging information including problem
formulation details.


{marker remarks}{...}
{title:Remarks}

{pstd}
{cmd:dea} version 2.0 features a dual-engine architecture.  Both engines
achieve full model parity across all sixteen DEA formulations.

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
{it:Interior Point Method (IPM) engine.}  The alternative engine leverages
Stata 16's {cmd:LinearProgram()} class, which implements a primal-dual
interior point method.  IPM offers polynomial-time worst-case complexity,
efficient handling of degeneracy common in DEA, and automatic
primal-dual fallback when the primal formulation encounters numerical
difficulties.

{pstd}
Both engines are validated to produce identical results within numerical
tolerance (typically < 1e-6) across all model-RTS-orientation
combinations, including the CRS/VRS/DRS decomposition used for VRS
Frontier Analysis.

{pstd}
{it:Automatic engine selection.}  When {cmd:engine(auto)} is specified,
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
{it:Backward compatibility.}  The former {cmd:dea_lp} command is now
integrated into {cmd:dea}.  Typing {cmd:dea_lp} automatically redirects
to {cmd:dea} with a notification message.

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
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) ort(i) gen(theta)}{p_end}

{pstd}VRS input-oriented (BCC model) with reference sets{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(vrs) ort(i) gen(te_vrs) dmu(dmu) reference}{p_end}

{pstd}Output-oriented{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(vrs) ort(o) gen(te_out)}{p_end}

    {title:Non-radial models}

{pstd}Slacks-based measure{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(sbm) rts(crs) gen(sbm_eff)}{p_end}

{pstd}Additive model{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(additive) rts(vrs) gen(add)}{p_end}

    {title:Economic efficiency}

{pstd}Cost efficiency with input prices{p_end}
{phang2}{cmd:. gen w1 = 10}{p_end}
{phang2}{cmd:. gen w2 = 2}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(cost) w(w1 w2) rts(vrs) gen(ce)}{p_end}

{pstd}Revenue efficiency with output prices{p_end}
{phang2}{cmd:. gen p1 = 5}{p_end}
{phang2}{cmd:. gen p2 = 8}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(revenue) p(p1 p2) gen(re)}{p_end}

{pstd}Profit efficiency with input and output prices{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(profit) w(w1 w2) p(p1 p2) rts(vrs) gen(pe)}{p_end}

{pstd}Allocative efficiency decomposition{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(allocative) w(w1 w2) rts(vrs) gen(ae)}{p_end}

    {title:Directional distance function}

{pstd}Default direction (proportional){p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(ddf) rts(vrs) gen(beta)}{p_end}

{pstd}Custom direction vector{p_end}
{phang2}{cmd:. gen gx1 = 1}{p_end}
{phang2}{cmd:. gen gx2 = 1}{p_end}
{phang2}{cmd:. gen gy1 = 1}{p_end}
{phang2}{cmd:. gen gy2 = 1}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(ddf) gx(gx1 gx2) gy(gy1 gy2) gen(b_unit)}{p_end}

    {title:Super-efficiency}

{pstd}Radial super-efficiency{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(supereff) rts(crs) ort(i) gen(se)}{p_end}

{pstd}Super-SBM{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(supersbm) rts(vrs) gen(ssbm)}{p_end}

    {title:Multiplier (dual) model}

{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(multiplier) rts(crs) ort(i) gen(mult) saveweights}{p_end}

    {title:Free disposal hull}

{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(fdh) gen(fdh_eff)}{p_end}

{pstd}Shorthand using {cmd:fdh} flag{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, fdh gen(fdh_eff2)}{p_end}

    {title:Imprecise DEA (IDEA)}

{pstd}Interval data with known lower/upper bounds on inputs{p_end}
{phang2}{cmd:. gen emp_lo = i_employee * 0.9}{p_end}
{phang2}{cmd:. gen emp_hi = i_employee * 1.1}{p_end}
{phang2}{cmd:. gen area_lo = i_area * 0.9}{p_end}
{phang2}{cmd:. gen area_hi = i_area * 1.1}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(imprecise) xlo(emp_lo area_lo) xhi(emp_hi area_hi) gen(idea)}{p_end}

{pstd}Interval data on both inputs and outputs{p_end}
{phang2}{cmd:. gen sales_lo = o_sales * 0.95}{p_end}
{phang2}{cmd:. gen sales_hi = o_sales * 1.05}{p_end}
{phang2}{cmd:. gen prof_lo = o_profits * 0.9}{p_end}
{phang2}{cmd:. gen prof_hi = o_profits * 1.1}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, xlo(emp_lo area_lo) xhi(emp_hi area_hi) ylo(sales_lo prof_lo) yhi(sales_hi prof_hi) gen(idea2)}{p_end}

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

{pstd}Contract bad outputs (SO2) while expanding good outputs{p_end}
{phang2}{cmd:. dea capital labor = electricity, model(undesirable) bad(so2) rts(vrs) gen(ddf_u)}{p_end}

{pstd}Contract inputs and bad outputs simultaneously{p_end}
{phang2}{cmd:. dea capital labor = electricity, model(undesirable) bad(so2) direction(xyb) gen(ddf_xyb)}{p_end}

    {title:SBM with undesirable outputs (Tone 2004)}

{pstd}Non-radial efficiency with bad outputs{p_end}
{phang2}{cmd:. dea capital labor = electricity, model(sbmundesirable) bad(so2) rts(vrs) gen(sbmu)}{p_end}

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

{pstd}Detect input congestion with VRS decomposition{p_end}
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
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(supereff) rts(vrs) ort(i) gen(se_vrs) reference}{p_end}
{pstd}In this VRS leave-one-out example, stores B and E cannot be dominated
by a convex combination of the remaining stores. Their scores are missing
and the command reports two failed DMUs; this illustrates infeasibility,
not a complete ranking. Inspect missing scores before interpreting results.{p_end}
{phang2}{cmd:. assert missing(se_vrs) if inlist(dmu, "store B", "store E")}{p_end}
{phang2}{cmd:. assert !missing(se_vrs) if !inlist(dmu, "store B", "store E")}{p_end}


{pstd}Save super-efficiency reference weights as variables{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, model(supereff) rts(crs) ort(i) gen(se_crs) reference saveref}{p_end}

    {title:Bootstrap bias correction}

{pstd}Use a fixed seed for this small demonstration. Increase replications
for substantive inference. See {help dea_sw2007} for a complete two-stage
example with environmental variables.{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(vrs) ort(i) gen(te_bc) bootstrap reps1(100)}{p_end}

    {title:Engine and syntax options}

{pstd}Using the Simplex engine{p_end}
{pstd}Force the package-supplied Two-Phase Revised Simplex engine instead of
automatic engine selection.{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) ort(i) engine(simplex) gen(theta_s)}{p_end}

{pstd}Using {cmd:method()} synonym{p_end}
{pstd}{cmd:method(simplex)} is accepted as a synonym for {cmd:engine(simplex)}
for backward compatibility.{p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) ort(i) method(simplex) gen(theta_m)}{p_end}

{pstd}Extended syntax form (Form B){p_end}
{phang2}{cmd:. dea o_sales o_profits, inputs(i_employee i_area) rts(vrs) gen(te_vrs2)}{p_end}

{pstd}Single-stage (no slack maximization){p_end}
{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(crs) stage(1) gen(theta_s1)}{p_end}

{pstd}Save results to a temporary dataset{p_end}
{phang2}{cmd:. tempfile dea_results}{p_end}

{phang2}{cmd:. dea i_employee i_area = o_sales o_profits, rts(vrs) gen(te) saving("`dea_results'") replace}{p_end}



{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{marker stored_results}{...}
{title:Stored results}

{pstd}
{cmd:dea} stores the following in {cmd:r()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of DMUs{p_end}
{synopt:{cmd:r(N_eff)}}number of efficient DMUs{p_end}
{synopt:{cmd:r(mean_eff)}}mean efficiency score{p_end}
{synopt:{cmd:r(min_eff)}}minimum efficiency score{p_end}
{synopt:{cmd:r(max_eff)}}maximum efficiency score{p_end}
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

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}}{cmd:dea}{p_end}
{synopt:{cmd:r(cmdline)}}full command as typed{p_end}
{synopt:{cmd:r(model)}}model specification{p_end}
{synopt:{cmd:r(rts)}}returns to scale{p_end}
{synopt:{cmd:r(orientation)}}orientation{p_end}
{synopt:{cmd:r(engine)}}solver engine used ({cmd:ipm} or {cmd:simplex}){p_end}
{synopt:{cmd:r(inputs)}}input variable names{p_end}
{synopt:{cmd:r(outputs)}}output variable names{p_end}
{synopt:{cmd:r(generate)}}name of generated efficiency variable{p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:r(efficiency)}}efficiency scores ({it:N} x 1){p_end}
{synopt:{cmd:r(lambda)}}reference weights ({it:N} x {it:N}); returned when
{cmd:reference} requests them and the selected model supports them{p_end}
{synopt:{cmd:r(slack_in)}}input slacks ({it:N} x {it:m}); returned for
simplex technical-efficiency calls with {cmd:stage(2)}{p_end}
{synopt:{cmd:r(slack_out)}}output slacks ({it:N} x {it:s}); returned for
simplex technical-efficiency calls with {cmd:stage(2)}{p_end}
{synopt:{cmd:r(vrsfrontier)}}CRS TE, VRS TE, NIRS TE (legacy column name
{cmd:DRS_TE}), scale efficiency, and RTS code
after {cmd:model(te) rts(vrs)} ({it:N} x 5){p_end}
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
{synopt:{helpb dea_game}}Game-theoretic DEA overview{p_end}
{synopt:{helpb dea_malmquist}}Malmquist productivity index{p_end}
{synopt:{helpb dea_nash}}Nash bargaining DEA{p_end}
{synopt:{helpb dea_network}}Network DEA{p_end}
{synopt:{helpb dea_sbmundesirable}}SBM with undesirable outputs{p_end}
{synopt:{helpb dea_stackelberg}}Stackelberg game DEA{p_end}
{synopt:{helpb dea_sw2007}}Simar-Wilson two-stage DEA{p_end}
{synopt:{helpb dea_window}}Window analysis DEA{p_end}
{synopt:{helpb milp_bb}}Mixed-integer linear programming by branch-and-bound{p_end}
{synoptline}
