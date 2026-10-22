*! version 1.1.0  20jul2026
*! Read-only status inspector for dea_reference_batch checkpoints
program define dea_reference_batch_status, rclass
    version 16.0
    syntax using/ [, JOURNAL(string)]
    preserve
    quietly use `"`using'"', clear
    foreach v in dmu_index checkpoint_N checkpoint_requested checkpoint_batch {
        capture confirm variable `v'
        if _rc {
            restore
            di as error "dea_reference_batch_status: incompatible checkpoint; `v' is missing"
            exit 459
        }
    }
    quietly summarize checkpoint_N, meanonly
    local N=r(max)
    quietly summarize checkpoint_requested, meanonly
    local requested=r(max)
    quietly summarize checkpoint_batch, meanonly
    local batches=r(max)
    quietly keep dmu_index
    quietly duplicates drop
    local completed=_N
    local remaining=`requested'-`completed'
    local rows=.
    restore

    local attempts=0
    local successes=0
    local failures=0
    local last_rc=.
    local aux=0
    local sparse=0
    local nodes=0
    local successful_dmus=0
    local successful_seconds=0
    if "`journal'"!="" {
        preserve
        quietly use `"`journal'"', clear
        local attempts=_N
        quietly count if return_code==0
        local successes=r(N)
        local failures=`attempts'-`successes'
        quietly summarize attempt, meanonly
        local lastattempt=r(max)
        quietly summarize return_code if attempt==`lastattempt', meanonly
        local last_rc=r(max)
        quietly summarize aux_lp_solves, meanonly
        local aux=r(sum)
        quietly summarize sparse_lp_solves, meanonly
        local sparse=r(sum)
        quietly summarize sparse_nodes, meanonly
        local nodes=r(sum)
        quietly summarize N_evaluated if return_code==0, meanonly
        local successful_dmus=r(sum)
        quietly summarize seconds if return_code==0, meanonly
        local successful_seconds=r(sum)
        restore
    }
    local fraction=`completed'/`requested'
    local avg_aux=.
    local avg_sparse=.
    local avg_nodes=.
    local avg_seconds=.
    local est_aux=.
    local est_sparse=.
    local est_nodes=.
    local est_seconds=.
    local eta=.
    if `successful_dmus'>0 {
        local avg_aux=`aux'/`successful_dmus'
        local avg_sparse=`sparse'/`successful_dmus'
        local avg_nodes=`nodes'/`successful_dmus'
        local est_aux=`remaining'*`avg_aux'
        local est_sparse=`remaining'*`avg_sparse'
        local est_nodes=`remaining'*`avg_nodes'
        if `successful_seconds'>0 {
            local avg_seconds=`successful_seconds'/`successful_dmus'
            local est_seconds=`remaining'*`avg_seconds'
            local eta=clock("`c(current_date)' `c(current_time)'","DMYhms")+1000*`est_seconds'
        }
    }
    if `remaining'==0 {
        local est_aux=0
        local est_sparse=0
        local est_nodes=0
        local est_seconds=0
        local eta=clock("`c(current_date)' `c(current_time)'","DMYhms")
    }
    return scalar N=`N'
    return scalar N_requested=`requested'
    return scalar N_completed=`completed'
    return scalar N_remaining=`remaining'
    return scalar complete=(`remaining'==0)
    return scalar batches_completed=`batches'
    return scalar attempts=`attempts'
    return scalar successes=`successes'
    return scalar failures=`failures'
    return scalar last_rc=`last_rc'
    return scalar aux_lp_solves=`aux'
    return scalar sparse_lp_solves=`sparse'
    return scalar sparse_nodes=`nodes'
    return scalar completed_fraction=`fraction'
    return scalar successful_evaluated=`successful_dmus'
    return scalar avg_aux_lp_per_dmu=`avg_aux'
    return scalar avg_sparse_lp_per_dmu=`avg_sparse'
    return scalar avg_sparse_nodes_per_dmu=`avg_nodes'
    return scalar avg_seconds_per_dmu=`avg_seconds'
    return scalar estimated_remaining_aux_lp=`est_aux'
    return scalar estimated_remaining_sparse_lp=`est_sparse'
    return scalar estimated_remaining_sparse_nodes=`est_nodes'
    return scalar estimated_remaining_total_lp=`est_aux'+`est_sparse'
    return scalar estimated_remaining_seconds=`est_seconds'
    return scalar eta_datetime=`eta'
    return local checkpoint `"`using'"'
    return local journal `"`journal'"'
    return local cmd "dea_reference_batch_status"
    di as text "Reference batch checkpoint: " as result `"`using'"'
    di as text "  completed/requested/remaining: " as result ///
        `completed' "/" `requested' "/" `remaining'
    if "`journal'"!="" di as text "  attempts/successes/failures: " ///
        as result `attempts' "/" `successes' "/" `failures'
    di as text "  progress: " as result %6.2f (100*`fraction') as text "%"
    if "`journal'"!="" & `remaining'>0 & `successful_dmus'>0 {
        di as text "  estimated remaining LP calls: " as result ///
            %12.1f (`est_aux'+`est_sparse')
        if `est_seconds'<. di as text "  estimated remaining seconds: " ///
            as result %12.1f `est_seconds'
        else di as text "  ETA unavailable: successful journal time is zero"
    }
end
