*! version 1.0.0  20jul2026
*! Nerlovian directional profit inefficiency

capture program drop dea_nerlovian
program define dea_nerlovian, rclass
    version 16.0

    local raw0 `"`0'"'
    if strpos(`"`raw0'"',"=")>0 {
        local _inputs ""
        gettoken _w 0 : 0, parse("=,")
        while "`_w'"!="=" & "`_w'"!="" {
            if "`_w'"=="," {
                di as error "syntax error: comma before '='"
                exit 198
            }
            local _inputs `_inputs' `_w'
            gettoken _w 0 : 0, parse("=,")
        }
        local _cpos=strpos(`"`0'"',",")
        if `_cpos'>0 {
            local _pre=substr(`"`0'"',1,`_cpos')
            local _post=substr(`"`0'"',`_cpos'+1,.)
            local 0 `"`_pre' inputs(`_inputs') `_post'"'
        }
        else local 0 `"`0', inputs(`_inputs')"'
    }

    syntax [varlist(numeric)] [if] [in], W(varlist numeric) P(varlist numeric) ///
        [Inputs(varlist numeric) Outputs(varlist numeric) GX(varlist numeric) ///
        GY(varlist numeric) DMU(varname) RTS(string) ENGINE(string) ///
        GENerate(name) TOLerance(real 1e-7) SAVing(string asis) REPLACE DETail]

    if "`outputs'"=="" local outputs `varlist'
    if "`inputs'"=="" | "`outputs'"=="" {
        di as error "specify inputs and outputs, or use: dea_nerlovian inputs = outputs"
        exit 198
    }
    if "`rts'"=="" local rts "vrs"
    if "`engine'"=="" local engine "auto"
    if "`generate'"=="" local generate "_npi"
    local rts=lower("`rts'")
    if inlist("`rts'","v","bcc") local rts "vrs"
    if "`rts'"!="vrs" {
        di as error "dea_nerlovian 1.0.0 requires rts(vrs) to ensure finite unconditional maximum profit"
        exit 198
    }
    if !inlist("`engine'","auto","simplex","ipm") {
        di as error "engine() must be auto, simplex, or ipm"
        exit 198
    }
    if `tolerance'<=0 {
        di as error "tolerance() must be positive"
        exit 198
    }
    if "`replace'"!="" & `"`saving'"'=="" {
        di as error "replace requires saving()"
        exit 198
    }

    local m : word count `inputs'
    local s : word count `outputs'
    local nw : word count `w'
    local np : word count `p'
    if `nw'!=`m' {
        di as error "w() must contain one price for each input"
        exit 198
    }
    if `np'!=`s' {
        di as error "p() must contain one price for each output"
        exit 198
    }
    if "`gx'"!="" {
        local ngx : word count `gx'
        if `ngx'!=`m' {
            di as error "gx() must contain one direction variable for each input"
            exit 198
        }
    }
    if "`gy'"!="" {
        local ngy : word count `gy'
        if `ngy'!=`s' {
            di as error "gy() must contain one direction variable for each output"
            exit 198
        }
    }

    local v_obs "`generate'_profit"
    local v_max "`generate'_maxprofit"
    local v_den "`generate'_normalization"
    local v_npi "`generate'_npi"
    local v_te  "`generate'_technical"
    local v_ae  "`generate'_allocative"
    local v_gap "`generate'_identity_gap"
    foreach v in `v_obs' `v_max' `v_den' `v_npi' `v_te' `v_ae' `v_gap' {
        confirm name `v'
        capture confirm new variable `v'
        if _rc {
            di as error "variable `v' already exists"
            exit 110
        }
    }

    marksample touse, novarlist
    markout `touse' `inputs' `outputs' `w' `p'
    if "`gx'"!="" markout `touse' `gx'
    if "`gy'"!="" markout `touse' `gy'
    if "`dmu'"!="" markout `touse' `dmu', strok
    quietly count if `touse'
    if r(N)==0 error 2000
    local N=r(N)
    foreach v of varlist `inputs' `outputs' `w' `p' {
        quietly count if `touse' & `v'<=0
        if r(N) {
            di as error "quantities and prices must be strictly positive; check `v'"
            exit 459
        }
    }
    local dirvars `gx' `gy'
    if "`dirvars'"!="" {
        foreach v of varlist `dirvars' {
            quietly count if `touse' & `v'<0
            if r(N) {
                di as error "direction variables must be nonnegative; check `v'"
                exit 459
            }
        }
    }

    tempvar beta
    local ddfdir ""
    if "`gx'"!="" local ddfdir "`ddfdir' gx(`gx')"
    if "`gy'"!="" local ddfdir "`ddfdir' gy(`gy')"
    quietly dea `outputs' if `touse', inputs(`inputs') model(ddf) rts(vrs) ///
        `ddfdir' engine(`engine') tolerance(`tolerance') generate(`beta')
    quietly generate double `v_obs'=.
    quietly generate double `v_max'=.
    quietly generate double `v_den'=.
    quietly generate double `v_npi'=.
    quietly generate double `v_te'=`beta' if `touse'
    quietly generate double `v_ae'=.
    quietly generate double `v_gap'=.

    mata: _dea_nerlovian_compute("`inputs'","`outputs'","`w'","`p'", ///
        "`gx'","`gy'","`touse'","`v_obs'","`v_max'","`v_den'","`v_npi'")
    quietly replace `v_ae'=`v_npi'-`v_te' if `touse'
    quietly replace `v_ae'=0 if `touse' & abs(`v_ae')<=max(`tolerance',1e-8)
    quietly replace `v_gap'=abs(`v_npi'-`v_te'-`v_ae') if `touse'

    quietly count if `touse' & (`v_den'<=0 | missing(`v_npi') | missing(`v_te'))
    if r(N) {
        di as error "one or more Nerlovian problems were not solved"
        exit 430
    }
    quietly count if `touse' & `v_ae'<-max(`tolerance',1e-6)
    if r(N) {
        di as error "Nerlovian allocative component is negative beyond tolerance"
        exit 459
    }
    quietly summarize `v_gap' if `touse', meanonly
    local max_gap=r(max)
    quietly summarize `v_npi' if `touse', meanonly
    local mean_npi=r(mean)
    quietly summarize `v_te' if `touse', meanonly
    local mean_te=r(mean)
    quietly summarize `v_ae' if `touse', meanonly
    local mean_ae=r(mean)
    quietly count if `touse' & `v_obs'<=0
    local N_nonpositive=r(N)

    label variable `v_obs' "Observed profit"
    label variable `v_max' "Maximum VRS profit"
    label variable `v_den' "Nerlovian direction value"
    label variable `v_npi' "Nerlovian profit inefficiency"
    label variable `v_te' "Directional technical inefficiency"
    label variable `v_ae' "Nerlovian allocative inefficiency"
    label variable `v_gap' "Nerlovian decomposition residual"

    di as text "{hline 72}"
    di as text "Nerlovian directional profit inefficiency"
    di as text "DMUs: " as result `N' as text "   Nonpositive observed profits: " as result `N_nonpositive'
    di as text "Mean overall NPI:  " as result %9.6f `mean_npi'
    di as text "Mean technical:    " as result %9.6f `mean_te'
    di as text "Mean allocative:   " as result %9.6f `mean_ae'
    di as text "Identity: NPI = technical + allocative; max residual = " as result %9.2e `max_gap'
    di as text "{hline 72}"
    if "`detail'"!="" {
        local listvars `v_obs' `v_max' `v_npi' `v_te' `v_ae'
        if "`dmu'"!="" local listvars `dmu' `listvars'
        list `listvars' if `touse', noobs separator(0)
    }

    if `"`saving'"'!="" {
        preserve
        quietly keep if `touse'
        local savevars `inputs' `outputs' `w' `p' `gx' `gy' ///
            `v_obs' `v_max' `v_den' `v_npi' `v_te' `v_ae' `v_gap'
        if "`dmu'"!="" local savevars `dmu' `savevars'
        keep `savevars'
        if "`replace'"!="" save `saving', replace
        else save `saving'
        restore
    }

    return clear
    return scalar N=`N'
    return scalar N_nonpositive_profit=`N_nonpositive'
    return scalar mean_npi=`mean_npi'
    return scalar mean_technical=`mean_te'
    return scalar mean_allocative=`mean_ae'
    return scalar max_identity_gap=`max_gap'
    return local rts "vrs"
    return local engine "`engine'"
    return local inputs "`inputs'"
    return local outputs "`outputs'"
    return local generate "`generate'"
    return local saving `"`saving'"'
end

mata:
void _dea_nerlovian_compute(string scalar inputs,string scalar outputs,
    string scalar wvars,string scalar pvars,string scalar gxvars,
    string scalar gyvars,string scalar touse,string scalar obsvar,
    string scalar maxvar,string scalar denvar,string scalar npivar)
{
    real matrix X,Y,W,P,GX,GY
    real colvector obs
    real scalar n,i,j,op,mp,den
    X=st_data(.,tokens(inputs),touse); Y=st_data(.,tokens(outputs),touse)
    W=st_data(.,tokens(wvars),touse); P=st_data(.,tokens(pvars),touse)
    GX=(gxvars=="" ? X : st_data(.,tokens(gxvars),touse))
    GY=(gyvars=="" ? Y : st_data(.,tokens(gyvars),touse))
    obs=selectindex(st_data(.,touse)); n=rows(X)
    for (i=1;i<=n;i++) {
        op=sum(P[i,.]:*Y[i,.])-sum(W[i,.]:*X[i,.])
        mp=.
        for (j=1;j<=n;j++) {
            if (mp==. | sum(P[i,.]:*Y[j,.])-sum(W[i,.]:*X[j,.])>mp)
                mp=sum(P[i,.]:*Y[j,.])-sum(W[i,.]:*X[j,.])
        }
        den=sum(P[i,.]:*GY[i,.])+sum(W[i,.]:*GX[i,.])
        st_store(obs[i],st_varindex(obsvar),op)
        st_store(obs[i],st_varindex(maxvar),mp)
        st_store(obs[i],st_varindex(denvar),den)
        if (den>0) st_store(obs[i],st_varindex(npivar),(mp-op)/den)
    }
}
end
