*! version 1.2.0  20jul2026
*! Restartable batch runner for dea_reference
program define dea_reference_batch, rclass
    version 16.0
    syntax [if] [in], INPUTS(varlist numeric min=1) OUTPUTS(varlist numeric min=1) ///
        CHECKPoint(string) [JOURNAL(string) BATCHSize(integer 25) MAXBatches(integer -1) ///
        RESUME REPLACE RTS(string) ORIentation(string) ENGINE(string) SCALING(string) ///
        DMU(varname) EVALuate(varname numeric) CANDidates(varname numeric) ///
        LEXicographic(string) LEXOrder(varname numeric) ///
        POLICYWeight(varname numeric) POLICYMatrix(name) POLICYSense(string) ///
        SPARSE BROAD RELativeinterior BASIC AUDIT ///
        MAXNodes(integer 10000) MAXAUXLP(integer -1) INTTol(real 1e-7) ///
        TOLerance(real 1e-8) PEERTOLerance(real 1e-7) OPTTOLerance(real -1) TRACE]

    if `batchsize'<1 | (`maxbatches'==0 | `maxbatches'<-1) {
        di as error "dea_reference_batch: batchsize() must be positive and maxbatches() positive or -1"
        exit 198
    }
    if "`resume'"!="" & "`replace'"!="" {
        di as error "dea_reference_batch: resume and replace may not be combined"
        exit 198
    }
    capture confirm file `"`checkpoint'"'
    local exists = !_rc
    if `exists' & "`resume'`replace'"=="" {
        di as error "dea_reference_batch: checkpoint exists; specify resume or replace"
        exit 602
    }
    if !`exists' & "`resume'"!="" {
        di as error "dea_reference_batch: resume checkpoint not found"
        exit 601
    }
    if "`replace'"!="" & "`journal'"!="" capture erase `"`journal'"'

    marksample touse, novarlist
    markout `touse' `inputs' `outputs'
    if "`dmu'"!="" {
        capture confirm string variable `dmu'
        if !_rc markout `touse' `dmu', strok
        else markout `touse' `dmu'
    }
    tempvar techindex eligible completed batchflag originalorder
    quietly generate long `originalorder'=_n
    quietly generate long `techindex'=sum(`touse') if `touse'
    quietly generate byte `eligible'=`touse'
    if "`evaluate'"!="" quietly replace `eligible'=`touse' & !missing(`evaluate') & `evaluate'!=0
    quietly count if `touse'
    local N=r(N)
    quietly count if `eligible'
    local requested=r(N)
    if `N'==0 | `requested'==0 {
        di as error "dea_reference_batch: no complete technology or evaluation observations"
        exit 2000
    }
    quietly generate byte `completed'=0

    * Fingerprint every data column that can affect the technology, selection,
    * ordering, policy objective, or saved identifiers.
    local sigvars `inputs' `outputs' `dmu' `evaluate' `candidates' `lexorder' `policyweight'
    local sigvars : list uniq sigvars
    tempvar rowgroup orderfinger
    quietly egen long `rowgroup'=group(`sigvars') if `touse', missing
    quietly generate double `orderfinger'=`rowgroup'*`techindex' if `touse'
    quietly _datasignature `sigvars' `orderfinger' if `touse', nonames
    local datasignature `"`r(datasignature)'"'
    local matrixsignature "none"
    if "`policymatrix'"!="" {
        preserve
        quietly clear
        quietly svmat double `policymatrix', names(__pm)
        quietly _datasignature __pm*, nonames
        local matrixsignature `"`r(datasignature)'"'
        restore
    }
    if "`rts'"=="" local rts "vrs"
    if "`orientation'"=="" local orientation "input"
    if "`engine'"=="" local engine "auto"
    if "`scaling'"=="" local scaling "row"
    if "`policysense'"=="" local policysense "min"
    local modelsignature ///
        "rts=`rts'|orientation=`orientation'|engine=`engine'|scaling=`scaling'|dmu=`dmu'|evaluate=`evaluate'|candidates=`candidates'|lexicographic=`lexicographic'|lexorder=`lexorder'|policyweight=`policyweight'|policymatrix=`policymatrix'|matrixsig=`matrixsignature'|policysense=`policysense'|sparse=`sparse'|broad=`broad'|relativeinterior=`relativeinterior'|basic=`basic'|audit=`audit'|maxnodes=`maxnodes'|inttol=`inttol'|tolerance=`tolerance'|peertolerance=`peertolerance'|opttolerance=`opttolerance'"

    tempfile done batchout combined
    local batchbase=0
    if "`resume'"!="" {
        preserve
        quietly use `"`checkpoint'"', clear
        foreach v in dmu_index peer_index checkpoint_N checkpoint_version checkpoint_data_signature checkpoint_model_signature {
            capture confirm variable `v'
            if _rc {
                restore
                di as error "dea_reference_batch: incompatible checkpoint; `v' is missing"
                exit 459
            }
        }
        quietly summarize checkpoint_N, meanonly
        if r(min)!=`N' | r(max)!=`N' {
            restore
            di as error "dea_reference_batch: checkpoint technology size differs from current sample"
            exit 459
        }
        quietly count if checkpoint_data_signature!=`"`datasignature'"'
        if r(N)>0 {
            restore
            di as error "dea_reference_batch: checkpoint data fingerprint differs from current marked sample"
            exit 459
        }
        quietly count if checkpoint_model_signature!=`"`modelsignature'"'
        if r(N)>0 {
            restore
            di as error "dea_reference_batch: checkpoint model fingerprint differs from current options"
            exit 459
        }
        quietly summarize checkpoint_batch, meanonly
        local batchbase=r(max)
        quietly keep dmu_index
        quietly duplicates drop
        quietly rename dmu_index `techindex'
        quietly generate byte `completed'=1
        quietly save "`done'", replace
        restore
        tempvar donemerge
        quietly merge m:1 `techindex' using "`done'", keep(master match) generate(`donemerge')
        quietly replace `completed'=1 if `donemerge'==3
        drop `donemerge'
        quietly sort `originalorder'
    }

    local forward ""
    foreach o in rts orientation engine scaling {
        if "``o''"!="" local forward `"`forward' `o'(``o'')"'
    }
    local forward `"`forward' maxnodes(`maxnodes') maxauxlp(`maxauxlp') inttol(`inttol')"'
    local forward `"`forward' tolerance(`tolerance') peertolerance(`peertolerance') opttolerance(`opttolerance')"'
    foreach o in dmu candidates lexorder policyweight {
        if "``o''"!="" local forward `"`forward' `o'(``o'')"'
    }
    foreach o in lexicographic policymatrix {
        if "``o''"!="" local forward `"`forward' `o'(``o'')"'
    }
    if "`policyweight'`policymatrix'"!="" ///
        local forward `"`forward' policysense(`policysense')"'
    foreach o in sparse broad relativeinterior basic audit trace {
        if "``o''"!="" local forward `"`forward' `o'"'
    }

    local batches=0
    local aux_total=0
    local sparse_lp_total=0
    local sparse_nodes_total=0
    while 1 {
        quietly count if `eligible' & !`completed'
        if r(N)==0 continue, break
        if `maxbatches'>0 & `batches'>=`maxbatches' continue, break
        quietly summarize `techindex' if `eligible' & !`completed', meanonly
        local first=r(min)
        quietly generate byte `batchflag'=`eligible' & !`completed' & ///
            `techindex'>=`first' & `techindex'<`first'+`batchsize'
        quietly count if `batchflag'
        if r(N)==0 {
            drop `batchflag'
            continue, break
        }
        quietly summarize `techindex' if `batchflag', meanonly
        local last=r(max)
        quietly count if `batchflag'
        local batch_n=r(N)
        local started=clock("`c(current_date)' `c(current_time)'","DMYhms")
        capture quietly dea_reference if `touse', inputs(`inputs') outputs(`outputs') ///
            evaluate(`batchflag') `forward' saving("`batchout'") replace
        local batchrc=_rc
        local batchaux=0
        local batchsparse=0
        local batchnodes=0
        if `batchrc'==0 {
            local batchaux=r(aux_lp_solves)
            local batchsparse=r(sparse_lp_solves)
            local batchnodes=r(sparse_nodes)
        }
        local stopped=clock("`c(current_date)' `c(current_time)'","DMYhms")
        if "`journal'"!="" quietly _dea_reference_batch_log, journal(`"`journal'"') ///
            batch(`=`batchbase'+`batches'+1') first(`first') last(`last') ///
            neval(`batch_n') rc(`batchrc') started(`started') stopped(`stopped') ///
            aux(`batchaux') sparselp(`batchsparse') nodes(`batchnodes')
        if `batchrc'!=0 {
            drop `batchflag'
            di as error "dea_reference_batch: batch failed; checkpoint unchanged"
            exit `batchrc'
        }
        local aux_total=`aux_total'+`batchaux'
        local sparse_lp_total=`sparse_lp_total'+`batchsparse'
        local sparse_nodes_total=`sparse_nodes_total'+`batchnodes'
        preserve
        quietly use "`batchout'", clear
        quietly generate long checkpoint_N=`N'
        quietly generate byte checkpoint_version=2
        quietly generate long checkpoint_requested=`requested'
        quietly generate long checkpoint_batch=`batchbase'+`batches'+1
        quietly generate str80 checkpoint_data_signature=`"`datasignature'"'
        quietly generate str2045 checkpoint_model_signature=`"`modelsignature'"'
        if `exists' | `batches'>0 {
            quietly append using `"`checkpoint'"'
        }
        quietly sort dmu_index peer_index
        quietly save "`combined'", replace
        quietly copy "`combined'" `"`checkpoint'"', replace
        restore
        local exists=1
        quietly replace `completed'=1 if `batchflag'
        drop `batchflag'
        local ++batches
        di as text "dea_reference_batch: completed batch `batches'; checkpoint updated"
    }

    quietly count if `eligible' & `completed'
    local ncompleted=r(N)
    quietly count if `eligible' & !`completed'
    local nremaining=r(N)
    return scalar N=`N'
    return scalar N_requested=`requested'
    return scalar N_completed=`ncompleted'
    return scalar N_remaining=`nremaining'
    return scalar batches_run=`batches'
    return scalar complete=(`nremaining'==0)
    return scalar aux_lp_solves=`aux_total'
    return scalar sparse_lp_solves=`sparse_lp_total'
    return scalar sparse_nodes=`sparse_nodes_total'
    return scalar total_lp_solves=`aux_total'+`sparse_lp_total'
    return local checkpoint `"`checkpoint'"'
    return local journal `"`journal'"'
    return local data_signature `"`datasignature'"'
    return local model_signature `"`modelsignature'"'
    return local cmd "dea_reference_batch"
end

program define _dea_reference_batch_log
    version 16.0
    syntax, JOURNAL(string) BATCH(integer) FIRST(integer) LAST(integer) ///
        NEVAL(integer) RC(integer) STARTED(real) STOPPED(real) ///
        AUX(integer) SPARSELP(integer) NODES(integer)
    tempfile entry merged
    preserve
    quietly clear
    quietly set obs 1
    quietly generate long attempt=1
    capture confirm file `"`journal'"'
    if !_rc {
        quietly use `"`journal'"', clear
        quietly summarize attempt, meanonly
        local next=r(max)+1
        quietly clear
        quietly set obs 1
        quietly generate long attempt=`next'
    }
    quietly generate long checkpoint_batch=`batch'
    quietly generate long first_dmu_index=`first'
    quietly generate long last_dmu_index=`last'
    quietly generate long N_evaluated=`neval'
    quietly generate int return_code=`rc'
    quietly generate str8 status=cond(`rc'==0,"success","failed")
    quietly generate double started=`started'
    quietly generate double stopped=`stopped'
    format started stopped %tcCCYY-NN-DD_HH:MM:SS
    quietly generate double seconds=(`stopped'-`started')/1000
    quietly generate long aux_lp_solves=`aux'
    quietly generate long sparse_lp_solves=`sparselp'
    quietly generate long sparse_nodes=`nodes'
    quietly save "`entry'", replace
    capture confirm file `"`journal'"'
    if !_rc quietly append using `"`journal'"'
    quietly sort attempt
    quietly save "`merged'", replace
    quietly copy "`merged'" `"`journal'"', replace
    restore
end
