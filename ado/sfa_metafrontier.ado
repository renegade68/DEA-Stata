*! version 1.0.1  20jul2026
*! Stochastic metafrontier: Battese-Rao-O'Donnell mixed LP and Huang et al. QML

capture program drop sfa_metafrontier
program define sfa_metafrontier, rclass
    version 16.0
    syntax varlist(numeric min=2) [if] [in], GROUP(varname) ///
        [METHOD(string) DISTribution(string) CM(varlist numeric) ///
         GENerate(name) VCE(string) DMU(varname) DETail REPLACE]

    gettoken depvar indepvars : varlist
    if "`method'" == "" local method huang
    if "`distribution'" == "" local distribution hnormal
    if "`generate'" == "" local generate _smf
    if "`vce'" == "" local vce robust
    local method = lower("`method'")
    local distribution = lower("`distribution'")
    local vce = lower("`vce'")

    if !inlist("`method'","huang","bro") {
        di as error "method() must be huang or bro"
        exit 198
    }
    if !inlist("`distribution'","hnormal","tnormal") {
        di as error "distribution() must be hnormal or tnormal"
        exit 198
    }
    if "`method'" == "bro" & "`cm'" != "" {
        di as error "cm() applies only to method(huang) distribution(tnormal)"
        exit 198
    }
    if "`distribution'" == "hnormal" & "`cm'" != "" {
        di as error "cm() requires distribution(tnormal)"
        exit 198
    }
    if "`method'" == "huang" & "`distribution'" == "tnormal" & ///
       inlist("`vce'","robust","cluster") {
        di as error "Stata frontier does not allow sandwich/cluster VCE with distribution(tnormal)"
        di as error "Use distribution(hnormal) vce(robust), or tnormal with vce(oim|opg) and label inference non-QML."
        exit 198
    }
    if !inlist("`vce'","robust","oim","opg") {
        di as error "vce() must be robust, oim, or opg"
        exit 198
    }

    marksample touse
    markout `touse' `group', strok
    if "`dmu'" != "" markout `touse' `dmu', strok
    quietly count if `touse'
    if r(N)==0 error 2000
    local N=r(N)

    local out_gf `generate'_groupfrontier
    local out_te `generate'_groupte
    local out_mf `generate'_metafrontier
    local out_tgr `generate'_tgr
    local out_mte `generate'_mte
    foreach x in `out_gf' `out_te' `out_mf' `out_tgr' `out_mte' {
        confirm name `x'
        capture confirm new variable `x'
        if _rc & "`replace'"=="" {
            di as error "variable `x' already exists; specify replace"
            exit 110
        }
    }

    tempvar gid xb_g te_g gf te mf tgr mte
    quietly egen long `gid'=group(`group') if `touse', label
    quietly summarize `gid' if `touse', meanonly
    local G=r(max)
    quietly generate double `gf'=.
    quietly generate double `te'=.
    local k : word count `indepvars'
    forvalues g=1/`G' {
        quietly count if `touse' & `gid'==`g'
        if r(N)<=`k'+2 {
            di as error "group `g' has too few observations for its stochastic frontier"
            exit 2001
        }
        quietly frontier `depvar' `indepvars' if `touse' & `gid'==`g', ///
            distribution(`distribution')
        quietly predict double `xb_g' if e(sample), xb
        quietly predict double `te_g' if e(sample), te
        quietly replace `gf'=`xb_g' if `gid'==`g' & `touse'
        quietly replace `te'=`te_g' if `gid'==`g' & `touse'
        drop `xb_g' `te_g'
    }

    quietly generate double `mf'=.
    quietly generate double `tgr'=.
    if "`method'"=="huang" {
        tempvar xb_m te_m
        local cmopt
        if "`cm'"!="" local cmopt cm(`cm')
        quietly frontier `gf' `indepvars' if `touse', ///
            distribution(`distribution') `cmopt' vce(`vce')
        quietly predict double `xb_m' if e(sample), xb
        quietly predict double `te_m' if e(sample), te
        quietly replace `mf'=`xb_m' if e(sample)
        quietly replace `tgr'=`te_m' if e(sample)
    }
    else {
        tempname bbro
        mata: _sfa_mf_bro("`gf'", "`indepvars'", "`touse'", "`bbro'")
        matrix colnames `bbro' = _cons `indepvars'
        quietly replace `mf'=`bbro'[1,1] if `touse'
        local j=1
        foreach x of varlist `indepvars' {
            local ++j
            quietly replace `mf'=`mf'+`bbro'[1,`j']*`x' if `touse'
        }
        quietly replace `tgr'=exp(`gf'-`mf') if `touse'
        quietly replace `tgr'=1 if `touse' & `tgr'>1 & `tgr'<1+1e-7
    }
    quietly generate double `mte'=`te'*`tgr' if `touse'

    label variable `gf' "Estimated group stochastic frontier (log)"
    label variable `te' "Group technical efficiency"
    label variable `mf' "Estimated stochastic metafrontier (log)"
    label variable `tgr' "Technology gap ratio"
    label variable `mte' "Metafrontier technical efficiency"

    quietly summarize `tgr' if `touse', meanonly
    local mean_tgr=r(mean)
    quietly summarize `te' if `touse', meanonly
    local mean_te=r(mean)
    quietly summarize `mte' if `touse', meanonly
    local mean_mte=r(mean)

    * Commit results only after every estimation step succeeds. This keeps
    * replace atomic if a group or second-stage frontier fails to converge.
    foreach x in `out_gf' `out_te' `out_mf' `out_tgr' `out_mte' {
        capture drop `x'
    }
    quietly generate double `out_gf'=`gf'
    quietly generate double `out_te'=`te'
    quietly generate double `out_mf'=`mf'
    quietly generate double `out_tgr'=`tgr'
    quietly generate double `out_mte'=`mte'
    label variable `out_gf' "Estimated group stochastic frontier (log)"
    label variable `out_te' "Group technical efficiency"
    label variable `out_mf' "Estimated stochastic metafrontier (log)"
    label variable `out_tgr' "Technology gap ratio"
    label variable `out_mte' "Metafrontier technical efficiency"

    di as text "{hline 72}"
    di as text "Stochastic metafrontier: " as result upper("`method'")
    di as text "Groups: " as result `G' as text "   Observations: " as result `N'
    di as text "Stage 2: " as result cond("`method'"=="huang", ///
        "stochastic frontier (Huang et al. 2014)", "deterministic LP (Battese et al. 2004)")
    di as text "Mean group TE: " as result %9.6f `mean_te'
    di as text "Mean TGR:      " as result %9.6f `mean_tgr'
    di as text "Mean meta TE:  " as result %9.6f `mean_mte'
    if "`method'"=="huang" & "`distribution'"=="tnormal" {
        di as error "Inference note: Stata prohibits robust VCE for truncated-normal frontier; reported VCE is `vce'."
    }
    if "`detail'"!="" {
        local listvars `group' `out_te' `out_tgr' `out_mte'
        if "`dmu'"!="" local listvars `dmu' `listvars'
        list `listvars' if `touse', noobs separator(0)
    }

    return clear
    return scalar N=`N'
    return scalar N_groups=`G'
    return scalar mean_group_te=`mean_te'
    return scalar mean_tgr=`mean_tgr'
    return scalar mean_mte=`mean_mte'
    return local cmd "sfa_metafrontier"
    return local method "`method'"
    return local distribution "`distribution'"
    return local vce "`vce'"
    return local depvar "`depvar'"
    return local indepvars "`indepvars'"
    return local generate "`generate'"
end

capture mata: mata drop _sfa_mf_bro()
mata:
void _sfa_mf_bro(string scalar gf, string scalar xs, string scalar tousename,
                 string scalar outname)
{
    real colvector y, use, bub
    real matrix X, A
    real rowvector c, lb, ub, b
    real scalar rc
    class LinearProgram scalar S
    use=st_data(.,tousename)
    y=select(st_data(.,gf),use)
    X=select(st_data(.,tokens(xs)),use)
    X=J(rows(X),1,1),X
    c=colsum(X)
    A=-X; bub=-y
    lb=J(1,cols(X),.); ub=J(1,cols(X),.)
    S=LinearProgram(); S.setMaxOrMin("min"); S.setCoefficients(c)
    S.setInequality(A,bub); S.setBounds(lb,ub); S.setTol(1e-10)
    rc=S.optimize()
    if (S.errorcode()!=0 | !S.converged()) _error(430)
    b=S.parameters()
    if (rows(b)>1) b=b'
    st_matrix(outname,b)
}
end
