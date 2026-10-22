*! version 1.2.0  20jul2026
program define milp_bb, rclass
    version 16.0
    syntax, C(name) INTEGER(numlist integer min=1) ///
        [SENSE(string) Aeq(name) beq(name) Aub(name) bub(name) ///
         LB(name) UB(name) MAXNodes(integer 10000) ///
         INTTol(real 1e-7) LPTol(real 1e-10) LPENGINE(string) ///
         LPSCALING(string) FORCELPFALLBACK TRACE]

    if "`sense'" == "" local sense "max"
    local sense = lower("`sense'")
    if !inlist("`sense'", "min", "max") {
        di as error "milp_bb: sense() must be min or max"
        exit 198
    }
    if `maxnodes' < 1 {
        di as error "milp_bb: maxnodes() must be positive"
        exit 198
    }
    if `inttol' <= 0 | `lptol' <= 0 {
        di as error "milp_bb: tolerances must be positive"
        exit 198
    }
    if "`lpengine'" == "" local lpengine "ipm"
    local lpengine = lower("`lpengine'")
    if !inlist("`lpengine'", "simplex", "ipm", "auto") {
        di as error "milp_bb: lpengine() must be simplex, ipm, or auto"
        exit 198
    }
    if "`lpscaling'" == "" local lpscaling "none"
    local lpscaling = lower("`lpscaling'")
    if !inlist("`lpscaling'", "none", "row") {
        di as error "milp_bb: lpscaling() must be none or row"
        exit 198
    }
    if "`forcelpfallback'" != "" & "`lpengine'" != "auto" {
        di as error "milp_bb: forcelpfallback is a test hook for lpengine(auto)"
        exit 198
    }
    local engine_label "shared LP contract"
    if "`lpengine'" == "ipm" local engine_label "LinearProgram"
    foreach M in c aeq beq aub bub lb ub {
        if "``M''" != "" {
            capture confirm matrix ``M''
            if _rc {
                di as error "milp_bb: matrix ``M'' not found"
                exit 111
            }
        }
    }
    if ("`aeq'" == "") != ("`beq'" == "") {
        di as error "milp_bb: aeq() and beq() must be specified together"
        exit 198
    }
    if ("`aub'" == "") != ("`bub'" == "") {
        di as error "milp_bb: aub() and bub() must be specified together"
        exit 198
    }

    tempname ints
    local intcsv : subinstr local integer " " ",", all
    matrix `ints' = (`intcsv')

    mata: milp_bb_solve("`c'", "`sense'", "`aeq'", "`beq'", ///
        "`aub'", "`bub'", "`lb'", "`ub'", "`ints'", ///
        `maxnodes', `inttol', `lptol', "`lpengine'", "`lpscaling'", ///
        "`forcelpfallback'", "`trace'")

    local status = r(status)
    local nodes = r(nodes)
    local lps = r(lp_solves)
    local simplex_lps = r(simplex_solves)
    local ipm_lps = r(ipm_solves)
    local fallbacks = r(lp_fallbacks)
    local objscale = r(objective_scale)
    local rowscalemin = r(row_scale_min)
    local rowscalemax = r(row_scale_max)
    local hasinc = r(hasinc)
    local objective = r(objective)
    local rootbound = r(rootbound)
    local gap = r(integrality_gap)
    tempname answer
    matrix `answer' = r(x)
    if `status' == 0 local status_text "optimal"
    else if `status' == 1 local status_text "infeasible"
    else if `status' == 2 local status_text "unbounded LP relaxation"
    else if `status' == 3 local status_text "node limit; incumbent available"
    else if `status' == 4 local status_text "node limit; no incumbent"
    else local status_text "input or LP solver error"
    return scalar status = `status'
    return scalar nodes = `nodes'
    return scalar lp_solves = `lps'
    return scalar simplex_solves = `simplex_lps'
    return scalar ipm_solves = `ipm_lps'
    return scalar lp_fallbacks = `fallbacks'
    return scalar objective_scale = `objscale'
    return scalar row_scale_min = `rowscalemin'
    return scalar row_scale_max = `rowscalemax'
    return scalar hasinc = `hasinc'
    return scalar objective = `objective'
    return scalar rootbound = `rootbound'
    return scalar integrality_gap = `gap'
    return matrix x = `answer'
    return local status_text "`status_text'"
    return local cmd "milp_bb"
    return local sense "`sense'"
    return local engine "`engine_label'"
    return local lpengine "`lpengine'"
    return local lpscaling "`lpscaling'"
    return local integer "`integer'"
    return scalar maxnodes = `maxnodes'
    return scalar inttol = `inttol'
    return scalar lptol = `lptol'
end

capture mata: mata drop milp_bb_solve()
mata:
void milp_bb_solve(string scalar cn, string scalar sense,
    string scalar aeqn, string scalar beqn,
    string scalar aubn, string scalar bubn,
    string scalar lbn, string scalar ubn, string scalar intn,
    real scalar maxnodes, real scalar inttol, real scalar lptol,
    string scalar lpengine, string scalar lpscaling,
    string scalar forcefallback_s, string scalar trace)
{
    real rowvector c, lb, ub, ilist, nlb, nub, x, bestx, leftlb, leftub, rightlb, rightub
    real colvector beq, bub
    real matrix Aeq, Aub, q
    real scalar n, nodes, lps, simplex_lps, ipm_lps, fallbacks
    real scalar nodefallback, objscale, nodeobjscale
    real scalar rowscalemin, rowscalemax, noderowscalemin, noderowscalemax
    real scalar hasinc, best, rootbound, status
    real scalar i, j, v, obj, integral, prune, maxfrac, f
    struct _DEA_LPResult scalar lpres

    c = st_matrix(cn)
    if (rows(c)>1) c = c'
    n = cols(c)
    ilist = st_matrix(intn)
    if (rows(ilist)>1) ilist = ilist'
    lb = (lbn=="" ? J(1,n,0) : st_matrix(lbn))
    ub = (ubn=="" ? J(1,n,.) : st_matrix(ubn))
    if (rows(lb)>1) lb=lb'
    if (rows(ub)>1) ub=ub'
    Aeq = (aeqn=="" ? J(0,n,.) : st_matrix(aeqn))
    beq = (beqn=="" ? J(0,1,.) : st_matrix(beqn))
    Aub = (aubn=="" ? J(0,n,.) : st_matrix(aubn))
    bub = (bubn=="" ? J(0,1,.) : st_matrix(bubn))

    status=0; nodes=0; lps=0; simplex_lps=0; ipm_lps=0; fallbacks=0
    objscale=1; rowscalemin=.; rowscalemax=0
    hasinc=0; rootbound=.; best=.
    bestx=J(1,n,.); q=J(0,2*n,.)
    if (n<1 | cols(lb)!=n | cols(ub)!=n |
        (rows(Aeq)>0 & cols(Aeq)!=n) | rows(Aeq)!=rows(beq) |
        (rows(Aub)>0 & cols(Aub)!=n) | rows(Aub)!=rows(bub) |
        any(ilist:<1) | any(ilist:>n) | any(lb:>ub)) status=5
    else q = (lb,ub)

    while (rows(q)>0 & status==0 & nodes<maxnodes) {
        nlb=q[rows(q),1..n]; nub=q[rows(q),(n+1)..(2*n)]
        q = (rows(q)==1 ? J(0,2*n,.) : q[|1,1 \ rows(q)-1,2*n|])
        nodes++; lps++

        if (lpengine=="ipm") ipm_lps++
        else simplex_lps++
        lpres=_dea_lpcontract_dispatch(c,sense,Aeq,beq,Aub,bub,
            nlb,nub,lpengine,lptol,lpscaling,forcefallback_s!="",trace,
            nodefallback,nodeobjscale,noderowscalemin,noderowscalemax)
        if (nodefallback) {
            fallbacks++; ipm_lps++
        }
        objscale=nodeobjscale
        rowscalemin=(rowscalemin==. ? noderowscalemin : min((rowscalemin,noderowscalemin)))
        rowscalemax=max((rowscalemax,noderowscalemax))
        if (lpres.status==1) {
            if (trace!="") printf("{txt}node %g: infeasible\n",nodes)
            continue
        }
        if (lpres.status==2) {
            status=2
            continue
        }
        if (lpres.status!=0) {
            status=5
            continue
        }
        obj=lpres.objective; x=lpres.x
        if (nodes==1) rootbound=obj

        prune=0
        if (hasinc) {
            if (sense=="max") prune=(obj<=best+inttol)
            else prune=(obj>=best-inttol)
        }
        if (prune) {
            if (trace!="") printf("{txt}node %g: bound pruned (%g)\n",nodes,obj)
            continue
        }

        integral=1; j=.; maxfrac=-1
        for (i=1;i<=cols(ilist);i++) {
            v=x[ilist[i]]; f=abs(v-round(v))
            if (f>inttol) {
                integral=0
                /* Stable most-fractional branching; declared order breaks ties. */
                if (f>maxfrac+inttol) {
                    maxfrac=f
                    j=ilist[i]
                }
            }
        }
        if (integral) {
            if (!hasinc | (sense=="max" & obj>best+inttol) |
                (sense=="min" & obj<best-inttol)) {
                hasinc=1; best=obj; bestx=x
                for (i=1;i<=cols(ilist);i++) bestx[ilist[i]]=round(bestx[ilist[i]])
                if (trace!="") printf("{res}node %g: incumbent %g\n",nodes,best)
            }
            continue
        }

        v=x[j]
        /* Push left then right so deterministic LIFO visits the right first. */
        rightlb=nlb; rightub=nub; rightlb[j]=ceil(v)
        leftlb=nlb; leftub=nub; leftub[j]=floor(v)
        if (leftlb[j]<=leftub[j]) q=q\(leftlb,leftub)
        if (rightlb[j]<=rightub[j]) q=q\(rightlb,rightub)
    }

    if (status==0 & rows(q)>0) status=(hasinc ? 3 : 4)
    else if (status==0 & !hasinc) status=1
    st_numscalar("r(status)",status)
    st_numscalar("r(nodes)",nodes); st_numscalar("r(lp_solves)",lps)
    st_numscalar("r(simplex_solves)",simplex_lps)
    st_numscalar("r(ipm_solves)",ipm_lps)
    st_numscalar("r(lp_fallbacks)",fallbacks)
    st_numscalar("r(objective_scale)",objscale)
    st_numscalar("r(row_scale_min)",(rowscalemin==. ? 1 : rowscalemin))
    st_numscalar("r(row_scale_max)",(rowscalemax==0 ? 1 : rowscalemax))
    st_numscalar("r(hasinc)",hasinc); st_numscalar("r(objective)",best)
    st_numscalar("r(rootbound)",rootbound)
    st_numscalar("r(integrality_gap)",(status==0 & hasinc ? 0 : .))
    st_matrix("r(x)",bestx)
}
end
