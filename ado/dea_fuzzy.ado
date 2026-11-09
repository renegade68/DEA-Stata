*! version 1.2.0  09nov2026
*! F11: disambiguate rounded alpha suffixes; return sorted alpha/output mapping.
*! 2026-11-01: accept the minimum option abbreviations advertised in help.
*! 2026-11-01: apply solver tolerance; check feasibility, bounds, and interval ordering.
*! 2026-11-01: correct NIRS -> DRS and NDRS -> IRS alias mapping.
*! 2026-11-01: honor if/in for evaluation and reference samples; reject empty samples.
*! 2026-11-01: preflight output names; explicit replacevars; protect model data.
*! dea_fuzzy.ado — Triangular alpha-cut DEA, Kao-Liu best/worst convention
*! 2026-11-03: complete-cohort summary; gap-aware nesting; exact center endpoints; safe ID display.
*! Requires: Stata 16.0 or later (LinearProgram class)
*!
*! AUTHORS:
*! --------
*!   Choonjoo Lee (corresponding author)
*!     Korea National Defense University
*!     Seoul, Korea
*!     bloom.rampike@gmail.com
*!
*! METHOD:
*! -------
*! Implements triangular alpha-cut DEA with Cartesian-product joint cuts.
*! For each α ∈ [0,1], triangular fuzzy inputs/outputs are sliced into
*! crisp intervals [x_L(α), x_U(α)], and two parametric LPs are solved
*! per DMU:
*!
*!   Pessimistic LP: DMU k uses unfavorable bounds (high inputs, low outputs);
*!                   peer frontier uses favorable bounds, with DMU k itself
*!                   entering the reference set at its same unfavorable bounds.
*!                   → lower bound on efficiency  θ_L(α)
*!
*!   Optimistic LP:  DMU k uses favorable bounds (low inputs, high outputs);
*!                   peer frontier uses unfavorable bounds, with DMU k itself
*!                   entering the reference set at its same favorable bounds.
*!                   → upper bound on efficiency  θ_U(α)
*!
*! α-cut bounds for triangular fuzzy number (x_lo, x_center, x_hi):
*!   x_L(α) = x_lo    + α * (x_center − x_lo)      [lo at α=0, center at α=1]
*!   x_U(α) = x_hi    − α * (x_hi     − x_center)  [hi at α=0, center at α=1]
*!
*! At α = 1: x_L = x_U = x_center  →  collapses to standard (crisp) DEA.
*! At α = 0: full fuzzy spread      →  maximum parametric uncertainty.
*!
*! For output-oriented models, φ (≥ 1) is converted to 1/φ ∈ (0,1] for
*! display, consistent with dea.ado conventions for IDEA output-oriented.
*!
*! REFERENCES:
*! -----------
*! Kao, C., & Liu, S.-T. (2000). Fuzzy efficiency measures in data
*!   envelopment analysis. Fuzzy Sets and Systems, 113(3), 427-437.
*!   → α-cut LP formulation (pessimistic/optimistic)
*!
*! Lertworasirikul, S., Fang, S.-C., Joines, J. A., & Nuttle, H. L. W.
*!   (2003). Fuzzy data envelopment analysis (DEA): a possibility approach.
*!   Fuzzy Sets and Systems, 139(2), 379-394.
*!   → Background: a different possibilistic DEA formulation
*!
*! Guo, P., & Tanaka, H. (2001). Fuzzy DEA: a perceptual evaluation method.
*!   Fuzzy Sets and Systems, 119(1), 149-160.
*!   → Background: a different fuzzy DEA formulation
*!
*! Cooper, W. W., Park, K. S., & Yu, G. (1999). IDEA and AR-IDEA: Models
*!   for dealing with imprecise data in DEA. Management Science, 45(4),
*!   597-607.
*!   → IDEA LP structure (pessimistic/optimistic); basis for this extension
*!
*! Charnes, A., Cooper, W. W., & Rhodes, E. (1978). Measuring the efficiency
*!   of decision making units. European Journal of Operational Research,
*!   2, 429-444.
*!
*! Banker, R. D., Charnes, A., & Cooper, W. W. (1984). Some models for
*!   estimating technical and scale inefficiencies in data envelopment
*!   analysis. Management Science, 30, 1078-1092.
*!
*! Ji, Y.-B., & Lee, C. (2010). Data envelopment analysis. Stata Journal,
*!   10(2), 267-280.
*!
*! Changelog:
*!   1.0.1  03may2026  Align best/worst reference-set convention with
*!                      deaR::modelfuzzy_kaoliu (evaluated DMU enters the
*!                      reference set at the same alpha-cut realization).
*!   1.0.0  07apr2026  Initial packaged version; based on earlier teaching/seminar versions


// ============================================================================
// Compile Mata function when this file is (re)loaded
// ============================================================================
capture mata: mata drop _dea_fuzzy_solve()
capture mata: mata drop _dea_fuzzy_lp()
capture mata: mata drop _dea_fuzzy_gridcheck()
capture mata: mata drop _dea_fuzzy_summary()
capture mata: mata drop _dea_fuzzy_display()

mata:
mata set matastrict on

// Solve and validate a single LP. A convergence flag alone is insufficient.
// Return: score, retry count, infeasible, unbounded, numerical failure, residual.
real rowvector _dea_fuzzy_lp(real rowvector c, real matrix A, real colvector b,
    real matrix E, real colvector d, real rowvector lb, real rowvector ub,
    real scalar tol, string scalar orientation)
{
    class LinearProgram scalar S
    real rowvector x, ix
    real scalar attempt, use_tol, status, first_status, value, score, residual
    real scalar check_tol
    check_tol=1e-8
    first_status=0
    for (attempt=1; attempt<=2; attempt++) {
        use_tol=tol
        if (attempt==2) {
            use_tol=min((tol/100,1e-12))
            if (use_tol<=0) use_tol=tol
        }
        S=LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A,b)
        if (rows(E)>0) S.setEquality(E,d)
        S.setBounds(lb,ub)
        S.setTrace("off")
        S.setTol(use_tol)
        value=S.optimize()
        status=S.errorcode()
        if (S.converged() & status==0 & value<.) {
            x=S.parameters()
            status=3
            if (cols(x)==cols(c) & !hasmissing(x)) {
                residual=max((0,max((A*x'-b):/(1:+abs(b):+abs(A)*abs(x)'))))
                if (rows(E)>0) residual=max((residual,max(abs(E*x'-d):/(1:+abs(d):+abs(E)*abs(x)'))))
                ix=selectindex(lb:<.)
                if (cols(ix)>0) residual=max((residual,max((lb[ix]-x[ix]):/(1:+abs(lb[ix])))))
                ix=selectindex(ub:<.)
                if (cols(ix)>0) residual=max((residual,max((x[ix]-ub[ix]):/(1:+abs(ub[ix])))))
                residual=max((residual,abs(c*x'-value)/(1+abs(value))))
                score=value
                if (orientation=="output") score=(value<0 ? -1/value : .)
                if (residual<=check_tol & score>=0 & score<=1+check_tol) {
                    return((score,attempt-1,0,0,0,residual))
                }
            }
        }
        else if (status==0) status=3
        if (attempt==1) first_status=status
    }
    // Preserve missing-score reporting for repeated native status 1/2 reports.
    // Native status is not a mathematical infeasibility/unboundedness certificate.
    // Other failures must abort before permanent result variables are changed.
    if (status==first_status & (status==1 | status==2)) {
        return((.,1,status==1,status==2,0,.))
    }
    return((.,1,0,0,1,.))
}

// ----------------------------------------------------------------------------
// _dea_fuzzy_solve()
//
// Solves pessimistic and optimistic Fuzzy DEA LPs for all DMUs at one
// α-cut level.
//
// Arguments:
//   inputs      - space-separated center input variable names
//   outputs     - space-separated center output variable names
//   xlo_vars    - space-separated lower-bound input variable names
//   xhi_vars    - space-separated upper-bound input variable names
//   ylo_vars    - space-separated lower-bound output variable names
//   yhi_vars    - space-separated upper-bound output variable names
//   alpha       - α-cut level ∈ [0,1]
//   orientation - "input" or "output"
//   rts         - "crs", "vrs", "drs", or "irs"
//   tol         - requested solver convergence tolerance; retried more tightly if needed
//   pess_var    - variable name to store pessimistic efficiency scores
//   opt_var     - variable name to store optimistic efficiency scores
//   touse       - sample marker variable name
//
// Pessimistic input-oriented LP:
//   min θ  s.t.  X_L(α)' λ  ≤  θ · x_{k,U}(α)    [frontier best, DMU worst]
//               −Y_U(α)' λ  ≤ −y_{k,L}(α)
//                         λ  ≥  0
//
// Optimistic input-oriented LP:
//   min θ  s.t.  X_U(α)' λ  ≤  θ · x_{k,L}(α)    [frontier worst, DMU best]
//               −Y_L(α)' λ  ≤ −y_{k,U}(α)
//                         λ  ≥  0
//
// Output-oriented: analogously; result stored as 1/φ ∈ (0,1].
// ----------------------------------------------------------------------------

void _dea_fuzzy_solve(
    string scalar inputs,
    string scalar outputs,
    string scalar xlo_vars,
    string scalar xhi_vars,
    string scalar ylo_vars,
    string scalar yhi_vars,
    real scalar   alpha,
    string scalar orientation,
    string scalar rts,
    real scalar   tol,
    string scalar pess_var,
    string scalar opt_var,
    string scalar touse,
    string scalar stats_name, | string scalar native_name)
{
    // ----------------------------------------------------------------
    // Load data
    // ----------------------------------------------------------------
    real matrix X_c, X_lo, X_hi, Y_c, Y_lo, Y_hi
    real matrix X_L, X_U, Y_L, Y_U, X_ref, Y_ref
    real colvector obs_indices
    real scalar n, m, s, k, obs_k

    X_c  = st_data(., tokens(inputs),   touse)
    X_lo = st_data(., tokens(xlo_vars), touse)
    X_hi = st_data(., tokens(xhi_vars), touse)
    Y_c  = st_data(., tokens(outputs),  touse)
    Y_lo = st_data(., tokens(ylo_vars), touse)
    Y_hi = st_data(., tokens(yhi_vars), touse)

    obs_indices = selectindex(st_data(., touse))
    n = rows(X_c)
    m = cols(X_c)
    s = cols(Y_c)

    real rowvector xs, ys
    xs=colmax(abs(X_c)\abs(X_lo)\abs(X_hi))
    ys=colmax(abs(Y_c)\abs(Y_lo)\abs(Y_hi))
    for (k=1;k<=cols(xs);k++) if (xs[k]==0) xs[k]=1
    for (k=1;k<=cols(ys);k++) if (ys[k]==0) ys[k]=1
    X_c=X_c:/xs; X_lo=X_lo:/xs; X_hi=X_hi:/xs
    Y_c=Y_c:/ys; Y_lo=Y_lo:/ys; Y_hi=Y_hi:/ys

    // ----------------------------------------------------------------
    // Compute α-cut bounds
    //   x_L(α) = x_lo + α*(x_center − x_lo)   [lo → center as α: 0→1]
    //   x_U(α) = x_hi − α*(x_hi − x_center)   [hi → center as α: 0→1]
    // ----------------------------------------------------------------
    // Convex combinations avoid cancellation; alpha=1 uses the identical
    // center matrix in all four paths, not separately rounded differences.
    X_L = (1-alpha):*X_lo + alpha:*X_c
    X_U = (1-alpha):*X_hi + alpha:*X_c
    Y_L = (1-alpha):*Y_lo + alpha:*Y_c
    Y_U = (1-alpha):*Y_hi + alpha:*Y_c
    if (alpha==0) {
        X_L=X_lo; X_U=X_hi; Y_L=Y_lo; Y_U=Y_hi
    }
    if (alpha==1) {
        X_L=X_c; X_U=X_c; Y_L=Y_c; Y_U=Y_c
    }

    // ----------------------------------------------------------------
    // Result vectors
    // ----------------------------------------------------------------
    real colvector eff_pess, eff_opt
    eff_pess = J(n, 1, .)
    eff_opt  = J(n, 1, .)

    // ----------------------------------------------------------------
    // LP building blocks
    // ----------------------------------------------------------------
    real rowvector c, lb, ub
    real matrix    A_ineq, A_eq
    real colvector b_ineq, b_eq, x0, y0
    real rowvector fit, lpstats
    real matrix nativefits
    lpstats=J(1,6,0)
    nativefits=J(n,12,.)

    for (k = 1; k <= n; k++) {

        obs_k = obs_indices[k]

        if (orientation == "input") {
            // ============================================================
            // INPUT-ORIENTED
            //
            // Variables: [λ₁, …, λₙ, θ]    (n+1 variables)
            //
            // Pessimistic:
            //   DMU k worst → x0 = X_U[k], y0 = Y_L[k]
            //   Peer frontier best; DMU k enters reference at worst case
            //   min θ  s.t.  X_L'λ − θ·X_U[k] ≤ 0   (m rows)
            //                −Y_U'λ             ≤ −Y_L[k]   (s rows)
            //                λ ≥ 0, θ ≥ 0
            //
            // Optimistic:
            //   DMU k best  → x0 = X_L[k], y0 = Y_U[k]
            //   Peer frontier worst; DMU k enters reference at best case
            //   min θ  s.t.  X_U'λ − θ·X_L[k] ≤ 0   (m rows)
            //                −Y_L'λ             ≤ −Y_U[k]   (s rows)
            //                λ ≥ 0, θ ≥ 0
            // ============================================================

            c  = J(1, n, 0), 1
            lb = J(1, n+1, 0)
            ub = J(1, n+1, .)

            // --- PESSIMISTIC ---
            x0 = X_U[k, ]'
            y0 = Y_L[k, ]'
            X_ref = X_L
            Y_ref = Y_U
            X_ref[k, ] = X_U[k, ]
            Y_ref[k, ] = Y_L[k, ]

            A_ineq = (X_ref', -x0) \ (-Y_ref', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)

            A_eq = J(0, n+1, .)
            b_eq = J(0, 1, .)
            if (rts == "vrs") {
                A_eq = (J(1, n, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n, -1), 0)
                b_ineq = b_ineq \ (-1)
            }

            fit=_dea_fuzzy_lp(c,A_ineq,b_ineq,A_eq,b_eq,lb,ub,tol,orientation)
            eff_pess[k]=fit[1]
            nativefits[k,1..6]=fit
            lpstats[1..5]=lpstats[1..5]+(1,fit[2..5])
            if (fit[6]<.) lpstats[6]=max((lpstats[6],fit[6]))

            // --- OPTIMISTIC ---
            x0 = X_L[k, ]'
            y0 = Y_U[k, ]'
            X_ref = X_U
            Y_ref = Y_L
            X_ref[k, ] = X_L[k, ]
            Y_ref[k, ] = Y_U[k, ]

            A_ineq = (X_ref', -x0) \ (-Y_ref', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)

            A_eq = J(0, n+1, .)
            b_eq = J(0, 1, .)
            if (rts == "vrs") {
                A_eq = (J(1, n, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n, -1), 0)
                b_ineq = b_ineq \ (-1)
            }

            fit=_dea_fuzzy_lp(c,A_ineq,b_ineq,A_eq,b_eq,lb,ub,tol,orientation)
            eff_opt[k]=fit[1]
            nativefits[k,7..12]=fit
            lpstats[1..5]=lpstats[1..5]+(1,fit[2..5])
            if (fit[6]<.) lpstats[6]=max((lpstats[6],fit[6]))
        }
        else {
            // ============================================================
            // OUTPUT-ORIENTED
            //
            // Variables: [λ₁, …, λₙ, φ]    (n+1 variables)
            //
            // Standard output-oriented LP (min −φ form):
            //   c       = [0,…,0, −1]
            //   A_ineq  = (X_ref', J(m,1,0)) row-block
            //             \ (−Y_ref', y0)     row-block
            //   b_ineq  = x0 \ J(s,1,0)
            //   lb      = J(1,n,0), 1     (φ ≥ 1)
            //   result  stored as 1/φ ∈ (0,1]  (consistent with IDEA)
            //
            // Pessimistic:
            //   DMU k worst → x0 = X_U[k], y0 = Y_L[k]
            //   Peer frontier best; DMU k enters reference at worst case
            //   → large φ (less efficient); stored as 1/φ_pess (small)
            //
            // Optimistic:
            //   DMU k best  → x0 = X_L[k], y0 = Y_U[k]
            //   Peer frontier worst; DMU k enters reference at best case
            //   → small φ (more efficient); stored as 1/φ_opt (large)
            // ============================================================

            c  = J(1, n, 0), -1
            lb = J(1, n, 0), 1
            ub = J(1, n+1, .)

            // --- PESSIMISTIC ---
            x0 = X_U[k, ]'
            y0 = Y_L[k, ]'
            X_ref = X_L
            Y_ref = Y_U
            X_ref[k, ] = X_U[k, ]
            Y_ref[k, ] = Y_L[k, ]

            A_ineq = (X_ref', J(m, 1, 0)) \ (-Y_ref', y0)
            b_ineq = x0 \ J(s, 1, 0)

            A_eq = J(0, n+1, .)
            b_eq = J(0, 1, .)
            if (rts == "vrs") {
                A_eq = (J(1, n, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n, -1), 0)
                b_ineq = b_ineq \ (-1)
            }

            fit=_dea_fuzzy_lp(c,A_ineq,b_ineq,A_eq,b_eq,lb,ub,tol,orientation)
            eff_pess[k]=fit[1]
            nativefits[k,1..6]=fit
            lpstats[1..5]=lpstats[1..5]+(1,fit[2..5])
            if (fit[6]<.) lpstats[6]=max((lpstats[6],fit[6]))

            // --- OPTIMISTIC ---
            x0 = X_L[k, ]'
            y0 = Y_U[k, ]'
            X_ref = X_U
            Y_ref = Y_L
            X_ref[k, ] = X_L[k, ]
            Y_ref[k, ] = Y_U[k, ]

            A_ineq = (X_ref', J(m, 1, 0)) \ (-Y_ref', y0)
            b_ineq = x0 \ J(s, 1, 0)

            A_eq = J(0, n+1, .)
            b_eq = J(0, 1, .)
            if (rts == "vrs") {
                A_eq = (J(1, n, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n, -1), 0)
                b_ineq = b_ineq \ (-1)
            }

            fit=_dea_fuzzy_lp(c,A_ineq,b_ineq,A_eq,b_eq,lb,ub,tol,orientation)
            eff_opt[k]=fit[1]
            nativefits[k,7..12]=fit
            lpstats[1..5]=lpstats[1..5]+(1,fit[2..5])
            if (fit[6]<.) lpstats[6]=max((lpstats[6],fit[6]))
        }
    }

    // ----------------------------------------------------------------
    // Store results back to Stata
    // ----------------------------------------------------------------
    if(args()>=15 & native_name!="") st_matrix(native_name,nativefits)
    st_matrix(stats_name,lpstats)
    st_store(obs_indices, st_varindex(pess_var), eff_pess)
    st_store(obs_indices, st_varindex(opt_var),  eff_opt)
}


// Check all observed levels, including comparisons across a missing middle cut.
// The final alpha=1 interval must collapse whenever both endpoints are known.
real scalar _dea_fuzzy_gridcheck(string scalar pvars,string scalar ovars,
    string scalar touse,string scalar alphas,real scalar tol)
{
    real matrix P,O
    real rowvector a
    real scalar i,j,lastp,lasto,bad
    P=st_data(.,tokens(pvars),touse);O=st_data(.,tokens(ovars),touse)
    a=strtoreal(tokens(alphas));bad=0
    for(i=1;i<=rows(P);i++) {
        lastp=.;lasto=.
        for(j=1;j<=cols(P);j++) {
            if(P[i,j]<. & O[i,j]<.) {
                bad=bad+(P[i,j]>O[i,j]+tol)
                if(a[j]==1) bad=bad+(abs(P[i,j]-O[i,j])>tol)
            }
            if(P[i,j]<.) {
                if(lastp<.) bad=bad+(lastp>P[i,j]+tol)
                lastp=(lastp<. ? max((lastp,P[i,j])) : P[i,j])
            }
            if(O[i,j]<.) {
                if(lasto<.) bad=bad+(O[i,j]>lasto+tol)
                lasto=(lasto<. ? min((lasto,O[i,j])) : O[i,j])
            }
        }
    }
    return(bad)
}

// Summary uses one common cohort: rows with both endpoints at every requested
// alpha. Per-cut availability counts remain visible; this is not defuzzification.
void _dea_fuzzy_summary(string scalar pvars,string scalar ovars,
    string scalar touse,string scalar alphas,string scalar dest,real scalar matrices)
{
    real matrix P,O,T
    string rowvector destinations
    real colvector complete,ii,pp,oo
    real rowvector a
    real scalar j,nc
    P=st_data(.,tokens(pvars),touse);O=st_data(.,tokens(ovars),touse)
    a=strtoreal(tokens(alphas));complete=(rowmissing((P,O)):==0)
    ii=selectindex(complete);nc=rows(ii);T=J(cols(P),11,.)
    for(j=1;j<=cols(P);j++) {
        T[j,1..5]=(a[j],sum(P[,j]:<.),sum(O[,j]:<.),sum((P[,j]:<.):&(O[,j]:<.)),nc)
        if(nc) {
            pp=P[ii,j];oo=O[ii,j]
            T[j,6..11]=(mean(pp),mean(oo),min(pp),max(pp),min(oo),max(oo))
        }
    }
    destinations=tokens(dest)
    st_matrix(destinations[1],T)
    if(matrices) {
        st_matrix(destinations[2],P);st_matrix(destinations[3],O)
        st_matrix(destinations[4],complete)
        st_matrix(destinations[5],selectindex(st_data(.,touse)))
    }
}

void _dea_fuzzy_display(string scalar pvars,string scalar ovars,
    string scalar touse,string scalar dmu,string scalar alphas)
{
    real matrix P,O
    string colvector labels
    real rowvector a
    real scalar i,j
    P=st_data(.,tokens(pvars),touse);O=st_data(.,tokens(ovars),touse)
    a=strtoreal(tokens(alphas))
    if(substr(st_vartype(dmu),1,3)=="str") labels=st_sdata(.,dmu,touse)
    else labels=strtrim(strofreal(st_data(.,dmu,touse),"%21.0g"))
    printf("\nDMU-level fuzzy efficiency\n%21s","DMU")
    for(j=1;j<=cols(P);j++) printf("  alpha=%g (pess/opt)",a[j])
    printf("\n")
    for(i=1;i<=rows(P);i++) {
        printf("%21s",labels[i])
        for(j=1;j<=cols(P);j++) printf("  %9.6g/%9.6g",P[i,j],O[i,j])
        printf("\n")
    }
}

end   // end mata


// ============================================================================
// Main Stata program
// ============================================================================
capture program drop dea_fuzzy
program define dea_fuzzy, rclass
    version 16.0

    // =========================================================================
    // Syntax
    // =========================================================================
    syntax [if] [in] , ///
        INputs(varlist numeric)  ///
        OUTputs(varlist numeric) ///
        [XLO(varlist numeric)    ///
         XHI(varlist numeric)    ///
         YLO(varlist numeric)    ///
         YHI(varlist numeric)    ///
         ALPha(numlist >=0 <=1 sort) ///
         Rts(string)             ///
         ORT(string)             ///
         ORientation(string)     ///
         DMU(varname)            ///
         GENerate(string)        ///
         REPLACEVARS             ///
         TOLerance(real 1e-10)    ///
         MATrix                  ///
         DETail                  ///
         DEBUG]

    // =========================================================================
    // Defaults
    // =========================================================================
    if "`alpha'"       == "" local alpha "0 0.25 0.5 0.75 1"
    if "`rts'"         == "" local rts "crs"
    if "`generate'"    == "" local generate "fdea_eff"

    // Resolve orientation (ort() is an alias)
    if "`ort'" != "" & "`orientation'" != "" {
        di as error "Specify only one of: ort() or orientation()"
        exit 198
    }
    if "`ort'" != "" local orientation "`ort'"
    if "`orientation'" == "" local orientation "input"

    // Normalize RTS
    local rts = lower("`rts'")
    if "`rts'" == "ccr"  local rts "crs"
    if "`rts'" == "bcc"  local rts "vrs"
    // NIRS: sum(lambda) <= 1; NDRS: sum(lambda) >= 1.
    if "`rts'" == "nirs" local rts "drs"
    if "`rts'" == "ndrs" local rts "irs"
    if !inlist("`rts'", "crs", "vrs", "drs", "irs") {
        di as error "rts() must be one of: crs (ccr), vrs (bcc), nirs (drs), ndrs (irs)"
        exit 198
    }

    // Normalize orientation
    local orientation = lower("`orientation'")
    if strpos("input", "`orientation'") == 1 local orientation "input"
    if strpos("output", "`orientation'") == 1 local orientation "output"
    if !inlist("`orientation'", "input", "output") {
        di as error "orientation() must be 'input' or 'output'"
        exit 198
    }

    if missing(`tolerance') | `tolerance' <= 0 {
        di as error "tolerance() must be positive and finite"
        exit 198
    }

    // Count inputs / outputs
    local m_in  : word count `inputs'
    local m_out : word count `outputs'

    // =========================================================================
    // Validate bounds: if provided, must match count of inputs/outputs
    // =========================================================================
    if "`xlo'" != "" {
        local n_xlo : word count `xlo'
        if `n_xlo' != `m_in' {
            di as error "xlo() must contain `m_in' variable(s) (same as inputs())"
            exit 198
        }
    }
    if "`xhi'" != "" {
        local n_xhi : word count `xhi'
        if `n_xhi' != `m_in' {
            di as error "xhi() must contain `m_in' variable(s) (same as inputs())"
            exit 198
        }
    }
    if "`ylo'" != "" {
        local n_ylo : word count `ylo'
        if `n_ylo' != `m_out' {
            di as error "ylo() must contain `m_out' variable(s) (same as outputs())"
            exit 198
        }
    }
    if "`yhi'" != "" {
        local n_yhi : word count `yhi'
        if `n_yhi' != `m_out' {
            di as error "yhi() must contain `m_out' variable(s) (same as outputs())"
            exit 198
        }
    }

    // Keep historical suffixes when unique. Disambiguate every member of a
    // colliding rounded bucket in sorted alpha order; reject repeated levels.
    // Use the syntax-parsed numlist, which may itself round very close values.
    local _previous_alpha=-1
    foreach a of numlist `alpha' {
        if `a' == `_previous_alpha' {
            di as error "alpha() contains repeated levels after numlist parsing: `a'"
            exit 198
        }
        local _previous_alpha `a'
        local asuf = string(round(`a' * 100))
        if "`_bucket_`asuf''" == "" local _bucket_`asuf'=0
        local ++_bucket_`asuf'
    }
    local _planned ""
    local _pessvars ""
    local _optvars ""
    local _alpha_collisions=0
    foreach a of numlist `alpha' {
        local asuf = string(round(`a' * 100))
        if `_bucket_`asuf'' > 1 {
            local _alpha_collisions=1
            if "`_index_`asuf''" == "" local _index_`asuf'=0
            local ++_index_`asuf'
            local asuf "`asuf'_`_index_`asuf''"
        }
        local _planned `_planned' `generate'_pess_a`asuf' `generate'_opt_a`asuf'
        local _pessvars `_pessvars' `generate'_pess_a`asuf'
        local _optvars `_optvars' `generate'_opt_a`asuf'
    }
    // Validate every output before any existing variable can be replaced.
    _dea_output_check, names(`_planned') ///
        protect(`inputs' `outputs' `xlo' `xhi' `ylo' `yhi' `dmu') `replacevars'

    // Warn if no fuzzy bounds supplied (degenerates to standard DEA)
    if "`xlo'`xhi'`ylo'`yhi'" == "" {
        di as text "Note: No fuzzy bounds specified; all variables treated as crisp."
        di as text "      Results will be identical to standard DEA at all α-levels."
    }

    // =========================================================================
    // Build sample
    // =========================================================================
    marksample touse, novarlist
    markout `touse' `inputs' `outputs' `xlo' `xhi' `ylo' `yhi'

    local _generated_dmu 0
    if "`dmu'" == "" {
        capture confirm new variable _dmu
        if _rc {
            di as error "_dmu already exists; specify dmu(_dmu) to use it or supply another dmu() variable"
            exit 110
        }
        tempvar _automatic_dmu
        qui gen str32 `_automatic_dmu' = "dmu" + string(_n) if `touse'
        local dmu "`_automatic_dmu'"
        local _generated_dmu 1
    }

    qui count if `touse'
    if r(N) == 0 {
        di as error "No observations in estimation sample"
        exit 2000
    }
    local n_dmus = r(N)

    // =========================================================================
    // Resolve bound variable lists
    // (use center values where bounds not specified → crisp for that variable)
    // =========================================================================
    local xlo_vars "`xlo'"
    local xhi_vars "`xhi'"
    local ylo_vars "`ylo'"
    local yhi_vars "`yhi'"

    if "`xlo_vars'" == "" local xlo_vars "`inputs'"
    if "`xhi_vars'" == "" local xhi_vars "`inputs'"
    if "`ylo_vars'" == "" local ylo_vars "`outputs'"
    if "`yhi_vars'" == "" local yhi_vars "`outputs'"

    foreach side in x y {
        if "`side'" == "x" {
            local _center "`inputs'"
            local _lo "`xlo_vars'"
            local _hi "`xhi_vars'"
        }
        else {
            local _center "`outputs'"
            local _lo "`ylo_vars'"
            local _hi "`yhi_vars'"
        }
        local _j=0
        foreach v of local _center {
            local ++_j
            local l : word `_j' of `_lo'
            local h : word `_j' of `_hi'
            quietly count if `touse' & (`l'<0 | `l'>`v' | `v'>`h')
            if r(N) {
                di as error "Fuzzy data require 0 <= lower <= center <= upper: `v'"
                exit 459
            }
        }
    }

    // At the smallest requested alpha, each lower-cut row must have a positive
    // input and output. Nesting then guarantees the same for every selected cut
    // and every coordinatewise realization in it. Individual zero entries remain valid.
    local _first_alpha : word 1 of `alpha'
    foreach side in x y {
        local centers "`inputs'"
        local lowers "`xlo_vars'"
        if "`side'"=="y" {
            local centers "`outputs'"
            local lowers "`ylo_vars'"
        }
        tempvar _positive
        quietly gen byte `_positive'=0 if `touse'
        local j=0
        foreach v of local centers {
            local ++j
            local l : word `j' of `lowers'
            quietly replace `_positive'=1 if `touse' & (1-`_first_alpha')*`l'+`_first_alpha'*`v'>0
        }
        quietly count if `touse' & `_positive'==0
        if r(N) {
            di as error "Each selected lower alpha-cut row requires a positive input and a positive output"
            di as error "All-zero input/output rows are outside the finite radial fuzzy model"
            exit 459
        }
    }

    // =========================================================================
    // Parse alpha list
    // =========================================================================
    local n_alpha : word count `alpha'
    local _has_center=0
    foreach a of numlist `alpha' {
        if `a'==1 local _has_center=1
    }
    if `_has_center' {
        quietly findfile _dea_fuzzy_load.ado
        local _fuzzy_loader "`r(fn)'"
        capture program drop _dea_fuzzy_load
        quietly run "`_fuzzy_loader'"
        quietly _dea_fuzzy_load
    }


    // =========================================================================
    // Orientation label for display
    // =========================================================================
    if "`orientation'" == "input" {
        local ort_label "Input-oriented"
        local eff_label "θ (input efficiency)"
    }
    else {
        local ort_label "Output-oriented"
        local eff_label "1/φ (output efficiency)"
    }

    // RTS label
    if "`rts'" == "crs" local rts_label "CRS"
    if "`rts'" == "vrs" local rts_label "VRS"
    if "`rts'" == "drs" local rts_label "DRS"
    if "`rts'" == "irs" local rts_label "IRS"

    // =========================================================================
    // Header
    // =========================================================================
    di as text "{hline `c(linesize)'}"
    di as text "dea_fuzzy v1.2 — Triangular alpha-cut DEA (Kao-Liu best/worst convention)"
    di as text "  `ort_label', `rts_label'"
    di as text "  DMUs: `n_dmus'  |  Inputs: `m_in'  |  Outputs: `m_out'"
    di as text "  α-levels: `alpha'"
    di as text "{hline `c(linesize)'}"

    // =========================================================================
    // Loop over α-levels: create temp vars, call Mata solver
    // =========================================================================
    tempname _lpstats _native _center_diag _center_counts _fctx
    local _center_attempts=0
    matrix `_center_counts'=J(1,4,0)
    local _lp_count=0
    local _lp_retries=0
    local _lp_infeasible=0
    local _lp_unbounded=0
    local _max_residual=0
    local _check_tol=1e-8
    local _use_tol=`tolerance'
    local _grid_retries=0
    forvalues _pass=1/2 {
        local ia=0
        local _pess_temp ""
        local _opt_temp ""
        foreach a of numlist `alpha' {
            local ++ia
            tempvar _pess_`ia' _opt_`ia'
            qui gen double `_pess_`ia''=.
            qui gen double `_opt_`ia''=.
            local _pess_temp `_pess_temp' `_pess_`ia''
            local _opt_temp `_opt_temp' `_opt_`ia''
            quietly {
                mata: _dea_fuzzy_solve(                           ///
                    "`inputs'", "`outputs'",                       ///
                    "`xlo_vars'", "`xhi_vars'",                   ///
                    "`ylo_vars'", "`yhi_vars'",                   ///
                    `a', "`orientation'", "`rts'", `_use_tol',     ///
                    "`_pess_`ia''", "`_opt_`ia''",                 ///
                    "`touse'", "`_lpstats'", "`_native'")
            }
            if `a'==1 {
                mata: `_fctx'=_c57_context_save()
                capture noisily {
                    mata: _dea_rec_abort();_rsm_audit_reset(0);_rsm_success_reset(1,0);_dea_link_reset(2,0)
                    mata: _f135_apply("`inputs'","`outputs'","`touse'","`_pess_`ia''","`_opt_`ia''", ///
                        "`rts'","`orientation'",`tolerance',"`_native'","`_lpstats'","`_center_diag'","`_center_counts'")
                }
                local _center_rc=_rc
                mata: _c57_context_restore(`_fctx')
                mata: mata drop `_fctx'
                if `_center_rc' exit `_center_rc'
                local _center_attempts=`_center_attempts'+`n_dmus'
            }
            local _lp_count=`_lp_count'+`_lpstats'[1,1]
            local _lp_retries=`_lp_retries'+`_lpstats'[1,2]
            local _lp_infeasible=`_lp_infeasible'+`_lpstats'[1,3]
            local _lp_unbounded=`_lp_unbounded'+`_lpstats'[1,4]
            local _max_residual=max(`_max_residual',`_lpstats'[1,6])
            if `_lpstats'[1,5]>0 {
                di as error "Fuzzy DEA numerical validation failed at alpha=`a' after retry"
                di as error "No permanent result variables were changed."
                exit 498
            }
            if "`debug'"!="" di as text "  [alpha=`a'] solver tolerance `_use_tol'; retries " `_lpstats'[1,2]
        }

        // Compare every available level, bridging missing interior cuts.
        tempname _grid_bad
        mata: st_numscalar("`_grid_bad'",_dea_fuzzy_gridcheck("`_pess_temp'","`_opt_temp'","`touse'","`alpha'",`_check_tol'))
        local _bad=scalar(`_grid_bad')
        if `_bad'==0 continue, break
        if `_pass'==2 {
            di as error "Fuzzy efficiency intervals fail ordering/nesting checks after retry"
            di as error "No permanent result variables were changed."
            exit 498
        }
        local _use_tol=min(`tolerance'/100,1e-12)
        if `_use_tol'<=0 local _use_tol=`tolerance'
        local _grid_retries=1
        di as text "Recomputing alpha grid at tighter tolerance after interval checks..."
    }
    local _lp_missing=0
    forvalues ia=1/`n_alpha' {
        foreach _kind in pess opt {
            qui count if `touse' & missing(`_`_kind'_`ia'')
            local _lp_missing=`_lp_missing'+r(N)
        }
    }
    if `_lp_missing'>0 {
        di as text "Note: `_lp_missing' scores unavailable; repeated native status 1/2 is not a mathematical certificate."
    }

    if `_has_center' {
        if `_center_counts'[1,4]>0 di as text "Note: alpha=1 scores recovered for " `_center_counts'[1,4] " DMUs; native failure reports remain in center_diagnostics."
        if `_center_counts'[1,2]>0 di as text "Warning: " `_center_counts'[1,2] " alpha=1 DMUs retain numerical scores without exact certification; inspect center_diagnostics."
        if `_center_counts'[1,1]>0 di as text "Warning: " `_center_counts'[1,1] " alpha=1 DMUs remain unavailable; certification limits do not prove infeasibility."
    }

    tempname _summary _summary_pess _summary_opt _summary_complete _summary_obs
    mata: _dea_fuzzy_summary("`_pess_temp'","`_opt_temp'","`touse'","`alpha'","`_summary' `_summary_pess' `_summary_opt' `_summary_complete' `_summary_obs'","`matrix'"!="")
    matrix colnames `_summary'=alpha N_pess N_opt N_pair N_complete mean_pess mean_opt min_pess max_pess min_opt max_opt
    local _n_complete=`_summary'[1,5]

    // =========================================================================
    // Save results to permanent named variables
    // =========================================================================
    local ia = 0
    foreach a of numlist `alpha' {
        local ++ia
        local vp : word `ia' of `_pessvars'
        local vo : word `ia' of `_optvars'

        // Publish the validated outputs; replacing existing columns requires replacevars.
        foreach _drop_output in `vp' {
            capture confirm variable `_drop_output', exact
            if !_rc quietly drop `_drop_output'
        }
        foreach _drop_output in `vo' {
            capture confirm variable `_drop_output', exact
            if !_rc quietly drop `_drop_output'
        }

        qui gen double `vp' = `_pess_`ia'' if `touse'
        qui gen double `vo' = `_opt_`ia''  if `touse'

        label var `vp' "Fuzzy DEA pessimistic efficiency (α=`a', `rts_label', `ort_label')"
        label var `vo' "Fuzzy DEA optimistic efficiency  (α=`a', `rts_label', `ort_label')"
    }

    di as text ""
    di as text "Fuzzy efficiency summary: `_n_complete'/`n_dmus' DMUs complete at every requested alpha"
    di as text "Means and ranges use this same cohort; counts show each cut's available endpoints."
    di as text " alpha       N_pess N_opt N_pair    Mean_pess    Mean_opt"
    forvalues ia=1/`n_alpha' {
        di as result %9.6g `_summary'[`ia',1] ///
            "  " %6.0f `_summary'[`ia',2] " " %6.0f `_summary'[`ia',3] ///
            " " %6.0f `_summary'[`ia',4] "  " %11.7g `_summary'[`ia',6] ///
            " " %11.7g `_summary'[`ia',7]
    }
    di as text "Ranges and full-precision summaries: r(summary). No alpha-weighted aggregate is computed."
    di as text "Variables created: `generate'_pess_a* (pessimistic)"
    di as text "                   `generate'_opt_a*  (optimistic)"
    if `_alpha_collisions' | "`detail'" != "" {
        di as text "Alpha / result variable mapping (sorted alpha order):"
        local ia=0
        foreach a of numlist `alpha' {
            local ++ia
            local vp : word `ia' of `_pessvars'
            local vo : word `ia' of `_optvars'
            di as text "  alpha=`a': `vp'  `vo'"
        }
    }

    if "`detail'"!="" {
        mata: _dea_fuzzy_display("`_pessvars'","`_optvars'","`touse'","`dmu'","`alpha'")
    }

    // =========================================================================
    // Return values
    // =========================================================================
    if `_generated_dmu' rename `_automatic_dmu' _dmu

    return local cmd         "dea_fuzzy"
    return local orientation "`orientation'"
    return local rts         "`rts'"
    return local alpha       "`alpha'"
    return local generate    "`generate'"
    return local pessvars    "`_pessvars'"
    return local optvars     "`_optvars'"
    return local formulation "triangular_kaoliu_basic"
    return local summary_sample "complete_alpha_grid"
    return scalar N_complete = `_n_complete'
    return scalar N_partial = `n_dmus'-`_n_complete'
    return matrix summary = `_summary'
    if "`matrix'"!="" {
        local _columns ""
        forvalues j=1/`n_alpha' {
            local _columns `_columns' alpha`j'
        }
        matrix colnames `_summary_pess'=`_columns'
        matrix colnames `_summary_opt'=`_columns'
        matrix colnames `_summary_complete'=complete
        matrix colnames `_summary_obs'=observation
        return matrix eff_pess=`_summary_pess'
        return matrix eff_opt=`_summary_opt'
        return matrix complete=`_summary_complete'
        return matrix obs_index=`_summary_obs'
    }
    return scalar n_center_unavailable=`_center_counts'[1,1]
    return scalar n_center_unverified=`_center_counts'[1,2]
    return scalar n_center_certified=`_center_counts'[1,3]
    return scalar n_center_recovered=`_center_counts'[1,4]
    return scalar center_attempts=`_center_attempts'
    return local center_scope "alpha1_original_data"
    if `_has_center' {
        matrix colnames `_center_diag'=observation quality reason recovered native_pess_score native_pess_retries native_pess_infeasible native_pess_unbounded native_pess_numerical native_pess_residual native_opt_score native_opt_retries native_opt_infeasible native_opt_unbounded native_opt_numerical native_opt_residual exact_ok phase cert_status cert_reason reserved5 reserved6 cert_ops reserved8 cert_limbs weight_obj_error weight_eq_error weight_ineq_error weight_obj_guard lambda_eq_error lambda_ineq_error lambda_error score_error arithmetic_error reserved19 guard_ops exact_reason map_reason zero_columns repair_stop repair_steps rts_code input_orientation
        return matrix center_diagnostics=`_center_diag'
    }
    return scalar n_alpha    = `n_alpha'
    return scalar N          = `n_dmus'
    return scalar m_in       = `m_in'
    return scalar m_out      = `m_out'
    return scalar tolerance  = `tolerance'
    return scalar check_tolerance = `_check_tol'
    return scalar lp_count = `_lp_count'
    return scalar lp_retries = `_lp_retries'
    return scalar grid_retries = `_grid_retries'
    return scalar lp_missing = `_lp_missing'
    return scalar lp_infeasible = `_lp_infeasible'
    return scalar lp_unbounded = `_lp_unbounded'
    return scalar max_primal_violation = `_max_residual'

end
