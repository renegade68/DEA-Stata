{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.2.0  20jul2026}{...}
{title:dea_reference_batch — Restartable batched reference diagnostics}

{title:Syntax}

{p 8 12 2}
{cmd:dea_reference_batch} {ifin}, {cmd:inputs(}{it:varlist}{cmd:)}
{cmd:outputs(}{it:varlist}{cmd:)} {cmd:checkpoint(}{it:filename}{cmd:)}
[{cmd:journal(}{it:filename}{cmd:)} {cmd:batchsize(}{it:#}{cmd:)} {cmd:maxbatches(}{it:#}{cmd:)}
{cmd:resume} {cmd:replace} and {help dea_reference:dea_reference options}]

{title:Description}

{pstd}
{cmd:dea_reference_batch} preserves the complete marked-sample frontier but
splits evaluated DMUs into batches. After every successful batch it atomically
replaces a long-format checkpoint containing all completed DMU--peer results.
{cmd:resume} skips completed evaluated-DMU indices. The data, marked sample,
and observation order must remain unchanged between runs.

{pstd}
Version 1.1 stores both a marked-data fingerprint and a canonical model-option
fingerprint. Resume fails before solving if an input, output, selector,
identifier, ordering/policy variable, policy matrix, or result-defining option
has changed. Recovery controls such as {cmd:batchsize()}, {cmd:maxbatches()},
and {cmd:maxauxlp()} may change between invocations.

{pstd}
{cmd:batchsize()} defaults to 25. {cmd:maxbatches()} deliberately stops after
the specified number of new batches; default -1 completes all remaining work.
Use it to fit a run into a time window and then invoke the same command with
{cmd:resume}. {cmd:replace} starts a new checkpoint. The two options are
mutually exclusive.

{pstd}
All model and selection options are passed to {cmd:dea_reference}. The
checkpoint adds {cmd:checkpoint_N}, {cmd:checkpoint_version}, and
{cmd:checkpoint_batch}. A checkpoint is updated only after a complete batch.

{pstd}
If {cmd:journal()} is supplied, every attempted batch adds one row containing
its DMU-index range, start/stop time, elapsed seconds, return code, status, and
continuous/sparse LP effort. Failed attempts leave the checkpoint unchanged.
With {cmd:replace}, an existing journal is also restarted.

{title:Stored results}

{p2colset 9 28 30 2}{...}
{p2col:{cmd:r(N_requested)}}eligible evaluated DMUs{p_end}
{p2col:{cmd:r(N_completed)}}completed evaluated DMUs{p_end}
{p2col:{cmd:r(N_remaining)}}remaining evaluated DMUs{p_end}
{p2col:{cmd:r(batches_run)}}new batches run in this invocation{p_end}
{p2col:{cmd:r(complete)}}1 if no eligible DMU remains{p_end}
{p2col:{cmd:r(aux_lp_solves)}}continuous LP calls in this invocation{p_end}
{p2col:{cmd:r(sparse_lp_solves)}}sparse B&B LP calls in this invocation{p_end}
{p2col:{cmd:r(total_lp_solves)}}total LP calls in this invocation{p_end}
{p2col:{cmd:r(checkpoint)}}checkpoint filename{p_end}
{p2col:{cmd:r(journal)}}journal filename{p_end}
{p2col:{cmd:r(data_signature)}}current marked-data fingerprint{p_end}
{p2col:{cmd:r(model_signature)}}canonical model fingerprint{p_end}
{p2colreset}{...}

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
{phang2}{cmd:. input id x1 x2 y1 y2 region}{p_end}
{phang2}{cmd:. 1 1 2 1 1 1}{p_end}
{phang2}{cmd:. 2 2 4 1 1 1}{p_end}
{phang2}{cmd:. 3 3 6 1 1 2}{p_end}
{phang2}{cmd:. 4 4 8 1 1 2}{p_end}
{phang2}{cmd:. end}{p_end}
{phang2}{cmd:. generate priority = id}{p_end}
{phang2}{cmd:. generate distance = id}{p_end}
{phang2}{cmd:. generate cost = id}{p_end}
{phang2}{cmd:. generate byte review = id == 2}{p_end}
{phang2}{cmd:. matrix D = (1,2,3,4)}{p_end}
{phang2}{cmd:. tempfile peer_bounds peer_audit}{p_end}
{phang2}{cmd:. tempfile peer_work peer_journal}{p_end}
{phang2}{cmd:. dea_reference_batch, inputs(x1 x2) outputs(y1 y2) dmu(id) checkpoint("`peer_work'") journal("`peer_journal'") batchsize(2) maxbatches(1) replace}{p_end}
{phang2}{cmd:. dea_reference_batch_status using "`peer_work'", journal("`peer_journal'")}{p_end}
{phang2}{cmd:. assert r(N_remaining) == 2}{p_end}
{phang2}{cmd:. dea_reference_batch, inputs(x1 x2) outputs(y1 y2) dmu(id) checkpoint("`peer_work'") journal("`peer_journal'") batchsize(2) resume}{p_end}
{phang2}{cmd:. dea_reference_batch_status using "`peer_work'", journal("`peer_journal'")}{p_end}
{phang2}{cmd:. assert r(N_remaining) == 0}{p_end}

{phang2}{cmd:. restore}{p_end}
{* DEAS_HELP_EXAMPLES_END}{...}

{title:Also see}

{psee}{help dea_reference_batch_status}, {help dea_reference}, {help milp_bb}
