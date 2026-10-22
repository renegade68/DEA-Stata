{smcl}
{* Help revision 2026-10-22.1: self-contained examples; command version retained.}{...}
{* *! version 1.1.0  20jul2026}{...}
{title:dea_reference_batch_status — Inspect a reference checkpoint}

{title:Syntax}
{p 8 12 2}{cmd:dea_reference_batch_status using} {it:checkpoint}
[{cmd:, journal(}{it:filename}{cmd:)}]

{title:Description}
{pstd}
This read-only command reports completed, requested, and remaining evaluated
DMUs without solving an LP. With {cmd:journal()}, it also reports attempts,
successes, failures, the last return code, and cumulative LP and B&B effort.
The dataset in memory is preserved.

{pstd}
Version 1.1 estimates remaining LP calls, B&B nodes, seconds, and completion
time from successful journal work per evaluated DMU. These are throughput
projections, not solver guarantees. ETA is missing until positive successful
elapsed time is available; completed work returns zero remaining effort.

{title:Stored results}
{p2colset 9 31 33 2}{...}
{p2col:{cmd:r(N_completed)}}completed evaluated DMUs{p_end}
{p2col:{cmd:r(N_remaining)}}remaining evaluated DMUs{p_end}
{p2col:{cmd:r(batches_completed)}}highest committed batch number{p_end}
{p2col:{cmd:r(attempts)}}journal attempts{p_end}
{p2col:{cmd:r(successes)}}successful attempts{p_end}
{p2col:{cmd:r(failures)}}failed attempts{p_end}
{p2col:{cmd:r(last_rc)}}last journal return code{p_end}
{p2col:{cmd:r(aux_lp_solves)}}cumulative successful continuous LP calls{p_end}
{p2col:{cmd:r(sparse_lp_solves)}}cumulative successful sparse LP calls{p_end}
{p2col:{cmd:r(completed_fraction)}}completed share from 0 to 1{p_end}
{p2col:{cmd:r(successful_evaluated)}}DMUs represented by successful attempts{p_end}
{p2col:{cmd:r(avg_seconds_per_dmu)}}observed successful seconds per DMU{p_end}
{p2col:{cmd:r(estimated_remaining_total_lp)}}projected remaining LP calls{p_end}
{p2col:{cmd:r(estimated_remaining_sparse_nodes)}}projected remaining B&B nodes{p_end}
{p2col:{cmd:r(estimated_remaining_seconds)}}projected remaining seconds{p_end}
{p2col:{cmd:r(eta_datetime)}}projected completion datetime in Stata milliseconds{p_end}
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

