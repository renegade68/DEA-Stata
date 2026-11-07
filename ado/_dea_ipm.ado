*! version 2.0.16  02nov2026
*! SBM: zero-term means, unit rescaling, and oriented super-SBM.
*! Interval IDEA: evaluate each DMU against the same realized dataset, including itself.
*! Forward non-TE tolerance and validate original LP constraints before returning results.
*! Return projection residual diagnostics without changing IPM acceptance.
*! Report radial and reference provenance without changing IPM solves.
*! Explicit failed-projection radial recovery; preserve original efficiency scores.
*! Bounded native retries and per-DMU evidence for TE and VRS auxiliaries.
*! Correct TE stage-2 RTS, projection weights, and incomplete-result reporting.
*! Correct the good-output beta sign in the undesirable-output DDF.
*! _dea_ipm.ado — IPM Engine for dea v2.0
*! Interior Point Method engine using LinearProgram() class (Stata 16+)
*!
*! AUTHORS:
*! --------
*!   Choonjoo Lee (corresponding author)
*!     Korea National Defense University
*!     Seoul, Korea
*!     bloom.rampike@gmail.com
*!
*!
*! NOTE: This IPM engine module is entirely new in v2.0 (no code from
*!   prior versions). It uses Stata 16+ LinearProgram() class with
*!   primal-dual fallback. Prior versions used ldea.do (Revised Simplex).
*!
*! MODEL REFERENCES (see dea.ado header for complete list):
*! --------------------------------------------------------
*!   [2]  Charnes, Cooper, & Rhodes (1978) — CCR (CRS)
*!   [3]  Banker, Charnes, & Cooper (1984) — BCC (VRS)
*!   [4]  Cooper, Seiford, & Tone (2007) — envelopment LP, cost/revenue/profit
*!   [6]  Andersen & Petersen (1993) — super-efficiency
*!   [7]  Chambers, Chung, & Färe (1996) — DDF
*!   [9]  Tone (2001) — SBM
*!   [10] Tone (2002) — Super-SBM
*!   [11] Deprins, Simar, & Tulkens (2006) — FDH
*!   [12] Cooper, Park, & Yu (1999) — IDEA
*!   [13] Färe, Grosskopf, & Lovell (1993) — RTS determination
*!   [15] Charnes & Cooper (1962) — Charnes-Cooper transformation
*! 
*! Changelog:
*!   2.0.6  05aug2026  Scale TE input/output columns before LinearProgram
*!                      stage 1/2 solves and restore slacks to source units.
*!   2.0.5  20jul2026  Align input multiplier u0 sign and DRS/IRS bounds with Cooper
*!   2.0.4  20jul2026  Return a common signed multiplier intercept
*!   2.0.3  20jul2026  Fix DRS/IRS stage-2 scalar min/max Mata syntax
*!   2.0.1  08feb2026  Fix: supereff output-oriented φ lower bound 1 → 1e-8
*!                      (φ≥1 prevented super-efficient DMU scores < 1)
*!   2.0.0  05feb2026  Initial packaged version; based on earlier teaching/seminar versions
*!
*! Features:
*!   - Unified envelopment: CRS, VRS, DRS, IRS
*!   - Stage 2 slack maximization
*!   - Primal → Dual (Multiplier) fallback
*!   - Lambda (reference weights) extraction
*!   - VRS Frontier analysis (SCALE, RTS)

/**
 * MATA FUNCTION INDEX:
 * ---------------------------------------------------------------------------
 * Core Solver:
 *   _dea_ipm_envelopment()     Unified envelopment: all RTS, all orientations
 *                               Returns: theta/phi, lambda, convergence info
 *   _dea_ipm_slacks()          Stage 2: slack maximization (theta fixed)
 *                               Returns: input slacks, output slacks
 *
 * Dual Fallback:
 *   _dea_ipm_dual_input()      Multiplier form, input-oriented
 *   _dea_ipm_dual_output()     Multiplier form, output-oriented
 *
 * VRS Frontier:
 *   _dea_ipm_vrsfrontier()     CRS→VRS→DRS + SCALE + RTS
 *
 * Utilities:
 *   _dea_ipm_count_results()   Count converged/failed DMUs
 * ---------------------------------------------------------------------------
 */


// ============================================================================
// Clear existing functions to allow clean reload
// ============================================================================
// Load the shared acceptance criterion before compiling solver functions.
capture mata: assert(_dea_lpcheck_ok(J(1,6,0),1e-8))
if _rc {
    quietly findfile _dea_lpcheck.ado
    quietly run "`r(fn)'"
}

capture mata: assert(_dea_nonte_lp_ready()==1)
if _rc {
    quietly findfile _dea_nonte_lp.ado
    quietly run "`r(fn)'"
}

capture mata: assert(_dea_teaudit_ready()==1)
if _rc {
    quietly findfile _dea_teaudit.ado
    quietly run "`r(fn)'"
}

capture mata: assert(_dea_native_retry_ready())
if _rc {
    quietly findfile _dea_native_retry.ado
    quietly run "`r(fn)'"
}
capture mata: mata drop _dea_ipm_audit_ready()
capture mata: mata drop _dea_ipm_finish()
capture mata: mata drop _dea_ipm_trace_store()
capture mata: mata drop _dea_ipm_trace_counts()
capture mata: mata drop _dea_ipm_envelopment()
capture mata: mata drop _dea_ipm_slacks()
capture mata: mata drop _dea_ipm_dual_input()
capture mata: mata drop _dea_ipm_dual_output()
capture mata: mata drop _dea_ipm_vrsfrontier()
capture mata: mata drop _dea_ipm_count_results()
capture mata: mata drop _dea_ipm_supereff()
capture mata: mata drop _dea_ipm_ddf()
capture mata: mata drop _dea_ipm_ddf_undesirable()
capture mata: mata drop _dea_ipm_additive()
capture mata: mata drop _dea_ipm_cost()
capture mata: mata drop _dea_ipm_revenue()
capture mata: mata drop _dea_ipm_profit()
capture mata: mata drop _dea_ipm_sbm()
capture mata: mata drop _dea_ipm_sbm_oriented()
capture mata: mata drop _dea_ipm_multiplier()
capture mata: mata drop _dea_ipm_supersbm()
capture mata: mata drop _dea_fdh()
capture mata: mata drop _dea_idea()

// Helper function for bootstrap CI computation
capture mata: mata drop _dea_bootstrap_ci()
mata:
real scalar _dea_ipm_audit_ready() return(4)
void _dea_bootstrap_ci(string scalar matname, 
                       real scalar lo_pct, 
                       real scalar hi_pct,
                       string scalar lo_var,
                       string scalar hi_var)
{
    real matrix B
    real scalar n, reps, i, lo_idx, hi_idx
    real rowvector sorted_row
    real colvector ci_lo, ci_hi
    
    B = st_matrix(matname)
    n = rows(B)
    reps = cols(B)
    
    lo_idx = max((1, ceil(reps * lo_pct / 100)))
    hi_idx = min((reps, floor(reps * hi_pct / 100)))
    
    ci_lo = J(n, 1, .)
    ci_hi = J(n, 1, .)
    
    for (i = 1; i <= n; i++) {
        sorted_row = sort(B[i, .]', 1)'
        ci_lo[i] = sorted_row[lo_idx]
        ci_hi[i] = sorted_row[hi_idx]
    }
    
    st_store(., st_varindex(lo_var), ci_lo)
    st_store(., st_varindex(hi_var), ci_hi)
}
end

mata:
mata set matastrict on

// ============================================================================
// 1. UNIFIED ENVELOPMENT SOLVER
// ============================================================================
/**
 * _dea_ipm_envelopment()
 *
 * Solves the DEA envelopment (primal) form for all DMUs with dual fallback.
 *
 * Supports all four RTS assumptions:
 *   CRS: no convexity constraint
 *   VRS: sum(lambda) = 1
 *   DRS: sum(lambda) <= 1  (Non-Increasing RTS)
 *   IRS: sum(lambda) >= 1  (Non-Decreasing RTS)
 *
 * Input-oriented:  min theta  s.t. X*lambda <= theta*x0, Y*lambda >= y0
 * Output-oriented: max phi    s.t. X*lambda <= x0,       Y*lambda >= phi*y0
 *
 * Arguments:
 *   inputs       - string: input variable names
 *   outputs      - string: output variable names
 *   eff_var      - string: efficiency score variable name
 *   lambda_var   - string: lambda variable names (space-separated), "" if none
 *   slack_in_var - string: input slack variable names, "" if none
 *   slack_out_var- string: output slack variable names, "" if none
 *   touse        - string: estimation sample marker
 *   orientation  - string: "input" or "output"
 *   rts          - string: "crs", "vrs", "drs", "irs"
 *   do_stage2    - real:   1 = run Stage 2 slack maximization, 0 = skip
 *   tol          - real:   numerical tolerance
 */
// Finish one native LP without changing its first accepted solution. The
// history retains all four variant slots; skipped attempts remain missing.
// Final 90 is residual rejection; 91 is contradictory native evidence.
real scalar _dea_ipm_finish(real scalar status, real scalar native,
    real scalar tol, real scalar iterations, real rowvector c,
    real matrix AE, real colvector BE, real matrix AU, real colvector BU,
    real rowvector lb, real rowvector ub, real rowvector sol,
    real rowvector audit)
{
    struct _DEA_NativeRetry scalar retry
    real matrix H
    real colvector codes
    H=J(4,6,.)
    H[1,.]=_dea_native_record(native,status,tol,iterations,sol,c,AE,BE,AU,BU,lb,ub)
    if (status!=0) {
        retry=_dea_native_retry(c,"min",AE,BE,AU,BU,lb,ub,tol,"")
        H[2..4,.]=retry.history
        codes=select(H[.,2],H[.,2]:<.)
        if (any(H[.,1]:==2) & (any(H[.,1]:==1) | retry.status==0)) status=91
        else if (retry.status==0) { sol=retry.x;status=0; }
        else if (any(codes:==90)) status=90
        else if (any(codes:==3)) status=3
        // Identical failure reports retain their original native status.
    }
    audit=(status,H[1,.],H[2,.],H[3,.],H[4,.])
    return(status)
}

// One row per selected DMU; a phase uses 25 columns (final + 4 x 6).
void _dea_ipm_trace_store(string scalar name, real matrix T,
    string rowvector prefixes)
{
    string rowvector names, fields
    real scalar p,v,j
    if (name=="") return
    names=J(1,0,"");fields=("native","checked","tol","iter","resid","obj")
    for (p=1;p<=cols(prefixes);p++) {
        names=(names,prefixes[p]+"_status")
        for (v=0;v<=3;v++) {
            for (j=1;j<=6;j++) names=(names,prefixes[p]+strofreal(v)+"_"+fields[j])
        }
    }
    st_matrix(name,T)
    st_matrixcolstripe(name,(J(cols(names),1,""),names'))
}

real rowvector _dea_ipm_trace_counts(real matrix T)
{
    real scalar p,v,b,retries,recovered,conflicts
    retries=0;recovered=0;conflicts=0
    for (p=1;p<=cols(T)/25;p++) {
        b=(p-1)*25
        for (v=1;v<=3;v++) retries=retries+sum(T[.,b+2+6*v]:<.)
        recovered=recovered+sum((T[.,b+1]:==0):&(T[.,b+3]:!=0))
        conflicts=conflicts+sum(T[.,b+1]:==91)
    }
    return((retries,recovered,conflicts))
}

void _dea_ipm_envelopment(string scalar inputs,
                           string scalar outputs,
                           string scalar eff_var,
                           string scalar lambda_var,
                           string scalar slack_in_var,
                           string scalar slack_out_var,
                           string scalar touse,
                           string scalar orientation,
                           string scalar rts,
                           real scalar do_stage2,
                           real scalar tol,
                           | string scalar stage1_var,
                           string scalar stage2_var,
                           string scalar trace_mat,
                           real scalar radialtol,
                           string scalar projection_mat,
                           string scalar projection_trace_mat,
                           string scalar score_source_var,
                           string scalar ref_stage_var,
                           string scalar diagnostics_mat)
{
    // ----------------------------------------------------------------
    // Data loading
    // ----------------------------------------------------------------
    real matrix X, Y, XR, YR, D
    real rowvector x_scale, y_scale
    real colvector obs_indices
    real scalar n, m, s, k, obs_k, j
    if (args()<12) stage1_var = ""
    if (args()<13) stage2_var = ""
    if (args()<14) trace_mat = ""
    if (args()<15) radialtol = 0
    if (args()<16) projection_mat = ""
    if (args()<17) projection_trace_mat = ""
    if (args()<18) score_source_var = ""
    if (args()<19) ref_stage_var = ""
    if (args()<20) diagnostics_mat = ""
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    if (diagnostics_mat!="") { XR=X; YR=Y; }
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)

    // Positive column scaling leaves radial DEA scores and lambdas unchanged,
    // while protecting LinearProgram stage 1 and stage 2 from heterogeneous
    // financial-statement magnitudes.  Stored slacks are restored below.
    x_scale = J(1, m, 1)
    y_scale = J(1, s, 1)
    for (j = 1; j <= m; j++) {
        x_scale[j] = max(abs(X[., j]))
        if (x_scale[j] <= 0 | x_scale[j] >= .) x_scale[j] = 1
        X[., j] = X[., j] / x_scale[j]
    }
    for (j = 1; j <= s; j++) {
        y_scale[j] = max(abs(Y[., j]))
        if (y_scale[j] <= 0 | y_scale[j] >= .) y_scale[j] = 1
        Y[., j] = Y[., j] / y_scale[j]
    }
    
    // Lambda variable list (may be empty)
    string rowvector lam_varlist, si_varlist, so_varlist
    real scalar has_lambda, has_slacks
    
    has_lambda = (lambda_var != "")
    has_slacks = (slack_in_var != "" & slack_out_var != "")
    
    if (has_lambda) lam_varlist = tokens(lambda_var)
    if (has_slacks) {
        si_varlist = tokens(slack_in_var)
        so_varlist = tokens(slack_out_var)
    }
    
    // ----------------------------------------------------------------
    // Solve each DMU
    // ----------------------------------------------------------------
    real colvector x0, y0
    real rowvector c, lb, ub, sol, lam_k
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar rc, eff_score, primal_ok
    real scalar phi_raw, stage2_status, primal_status
    real rowvector p_trace, s_trace, d_trace, adjusted_trace
    real scalar radial0, radialused, radialtry, step, original_stage2, adjusted_status
    real matrix T, P, PT
    class LinearProgram scalar S
    
    // Stage 2 variables
    real rowvector slacks_in, slacks_out
    
    D=J((diagnostics_mat!="" ? n : 0),23,.)
    T=J((trace_mat!="" ? n : 0),75,.)
    P=J((projection_mat!="" ? n : 0),8,.)
    PT=J((projection_trace_mat!="" ? n : 0),25,.)
    for (k = 1; k <= n; k++) {
        
        obs_k = obs_indices[k]
        x0 = X[k, ]'   // column vector (m × 1)
        y0 = Y[k, ]'   // column vector (s × 1)
        primal_ok = 0
        eff_score = .
        lam_k = J(1, n, .)
        slacks_in = J(1, m, .)
        slacks_out = J(1, s, .)
        stage2_status = .
        sol=J(1,n+1,.)
        p_trace=s_trace=d_trace=J(1,25,.)
        
        // ============================================================
        // PRIMAL (ENVELOPMENT)
        // ============================================================
        if (orientation == "input") {
            // ---------------------------------------------------------
            // Input-oriented: min theta
            //   s.t. X'*lambda - theta*x0 <= 0   (m rows)
            //        -Y'*lambda           <= -y0  (s rows)
            //        [convexity constraint]
            //        lambda >= 0, theta >= 0
            // Variables: [lambda(1..n), theta]
            // ---------------------------------------------------------
            
            c = J(1, n, 0), 1    // min theta
            
            A_ineq = (X', -x0) \ (-Y', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)
            
            lb = J(1, n+1, 0)
            ub = J(1, n+1, .)
            
            // Convexity constraint based on RTS
            A_eq = J(0, n+1, .)
            b_eq = J(0, 1, .)
            
            if (rts == "vrs") {
                // sum(lambda) = 1
                A_eq = (J(1, n, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                // sum(lambda) <= 1 → add as inequality
                A_ineq = A_ineq \ (J(1, n, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                // sum(lambda) >= 1 → -sum(lambda) <= -1
                A_ineq = A_ineq \ (J(1, n, -1), 0)
                b_ineq = b_ineq \ (-1)
            }
            // CRS: no convexity constraint
            
            // Solve
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) {
                S.setEquality(A_eq, b_eq)
            }
            S.setBounds(lb, ub)
            S.setTrace("off")
            S.setTol(tol)
            
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                sol = S.parameters()
                if (_dea_lpcheck_ok(_dea_lpcheck(sol,A_eq,b_eq,A_ineq,b_ineq,lb,ub),tol)) {
                    eff_score = c*sol'
                    lam_k = sol[1..n]
                    primal_ok = (eff_score<.)
                }
            }
        }
        else {
            // ---------------------------------------------------------
            // Output-oriented: max phi  (= min -phi)
            //   s.t. X'*lambda           <= x0   (m rows)
            //        -Y'*lambda + phi*y0 <= 0    (s rows)
            //        [convexity constraint]
            //        lambda >= 0, phi >= 1
            // Variables: [lambda(1..n), phi]
            // ---------------------------------------------------------
            
            c = J(1, n, 0), -1   // min(-phi) = max(phi)
            
            A_ineq = (X', J(m, 1, 0)) \ (-Y', y0)
            b_ineq = x0 \ J(s, 1, 0)
            
            lb = J(1, n, 0), 1   // phi >= 1
            ub = J(1, n+1, .)
            
            // Convexity constraint
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
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) {
                S.setEquality(A_eq, b_eq)
            }
            S.setBounds(lb, ub)
            S.setTrace("off")
            S.setTol(tol)
            
            rc = S.optimize()
            
            if (S.converged() && S.errorcode() == 0) {
                sol = S.parameters()
                if (_dea_lpcheck_ok(_dea_lpcheck(sol,A_eq,b_eq,A_ineq,b_ineq,lb,ub),tol)) {
                    phi_raw = -c*sol'
                    if (phi_raw>0 & phi_raw<.) {
                        eff_score = 1 / phi_raw
                        lam_k = sol[1..n]
                        primal_ok = 1
                    }
                }
            }
        }
        
        primal_status=(primal_ok ? 0 : (S.errorcode()!=0 ? S.errorcode() :
            (S.converged() ? 90 : 3)))
        primal_status=_dea_ipm_finish(primal_status,S.errorcode(),tol,S.iterations(),
            c,A_eq,b_eq,A_ineq,b_ineq,lb,ub,sol,p_trace)
        if (!primal_ok & primal_status==0) {
            eff_score=c*sol'
            if (orientation=="output") eff_score=(eff_score<0 ? -1/eff_score : .)
            primal_ok=(eff_score<.)
            if (primal_ok) lam_k=sol[1..n]
            else { primal_status=90;p_trace[1]=90; }
        }

        // ============================================================
        // DUAL FALLBACK (if primal failed)
        // ============================================================
        if (!primal_ok & !any(p_trace[(2,8,14,20)]:==2)) {
            if (orientation == "input") {
                eff_score = _dea_ipm_dual_input(X, Y, k, rts, tol, d_trace)
            }
            else {
                eff_score = _dea_ipm_dual_output(X, Y, k, rts, tol, d_trace)
            }
            lam_k = J(1, n, .)  // lambda unavailable from dual
        }
        
        // ============================================================
        // STAGE 2: SLACK MAXIMIZATION (if requested and Stage 1 OK)
        // ============================================================
        if (do_stage2 && primal_ok && eff_score != .) {
            stage2_status = _dea_ipm_slacks(X, Y, k, orientation, rts,
                            eff_score, lam_k,
                            slacks_in, slacks_out, tol, s_trace)
        }
        // Preserve the original fixed-radial history independently of recovery.
        if (trace_mat!="") T[k,.]=(p_trace,s_trace,d_trace)
        original_stage2=stage2_status
        radial0=(orientation=="input" ? eff_score : 1/eff_score)
        radialused=(stage2_status==0 ? radial0 : .)
        radialtry=.;adjusted_status=.;adjusted_trace=J(1,25,.)
        // A user-declared radial band permits one changed RHS after ordinary
        // native retries fail. Never conceal contradictory/unbounded reports.
        if (radialtol>0 & do_stage2 & primal_ok & stage2_status!=0 &
            stage2_status<. & !any(s_trace[(2,8,14,20)]:==2) &
            stage2_status!=91 & radial0>0 & radial0<.) {
            step=radialtol*(1+abs(radial0))
            if (step>0 & step<.) {
                radialtry=(orientation=="input" ? radial0+step : max((1,radial0-step)))
                if (radialtry>0 & radialtry<. &
                    ((orientation=="input" & radialtry>radial0) |
                     (orientation=="output" & radialtry<radial0))) {
                    adjusted_status=_dea_ipm_slacks(X,Y,k,orientation,rts,
                        eff_score,lam_k,slacks_in,slacks_out,tol,adjusted_trace,radialtry)
                    stage2_status=adjusted_status
                    if (adjusted_status==0) radialused=radialtry
                }
                else radialtry=.
            }
        }
        if (diagnostics_mat!="") D[k,.]=_dea_teaudit(X,Y,XR,YR,k,orientation,rts,
            x_scale,y_scale,lam_k,slacks_in,slacks_out,radial0,radialused,stage2_status,tol)
        if (projection_mat!="") P[k,.]=(radial0,radialused,
            radialused-radial0,radialtry,original_stage2,adjusted_status,
            stage2_status,(adjusted_status==0))
        if (projection_trace_mat!="") PT[k,.]=adjusted_trace
        // A requested stage-2 projection must not expose stage-1 references.
        // Retain the radial score independently of projection availability.
        if (do_stage2 && stage2_status != 0) lam_k = J(1, n, .)
        if (stage1_var != "") st_store(obs_k, st_varindex(stage1_var), eff_score>=.)
        if (stage2_var != "") st_store(obs_k, st_varindex(stage2_var), stage2_status)
        if (score_source_var!="") st_store(obs_k,st_varindex(score_source_var),
            (eff_score>=. ? 0 : (primal_ok ? 1 : 2)))
        if (ref_stage_var!="") st_store(obs_k,st_varindex(ref_stage_var),
            (do_stage2 ? (stage2_status==0 ? 2 : 0) : (primal_ok ? 1 : 0)))
        
        // ============================================================
        // STORE RESULTS
        // ============================================================
        st_store(obs_k, st_varindex(eff_var), eff_score)
        
        if (has_lambda) {
            for (j = 1; j <= n; j++) {
                st_store(obs_k, st_varindex(lam_varlist[j]), lam_k[j])
            }
        }
        
        if (has_slacks && do_stage2) {
            for (j = 1; j <= m; j++) {
                st_store(obs_k, st_varindex(si_varlist[j]), slacks_in[j] * x_scale[j])
            }
            for (j = 1; j <= s; j++) {
                st_store(obs_k, st_varindex(so_varlist[j]), slacks_out[j] * y_scale[j])
            }
        }
    }
    
    _dea_teaudit_store(diagnostics_mat,D)
    _dea_ipm_trace_store(trace_mat,T,("p","s","d"))
    if (projection_mat!="") {
        st_matrix(projection_mat,P)
        st_matrixcolstripe(projection_mat,(J(8,1,""),
            ("radial_initial","radial_stage2","radial_adjustment","radial_target",
             "status_initial","status_retry","status_final","adjusted")'))
    }
    _dea_ipm_trace_store(projection_trace_mat,PT,"a")

    // Store convergence counts
    _dea_ipm_count_results(eff_var, touse, n)
}


// ============================================================================
// 2. STAGE 2: SLACK MAXIMIZATION
// ============================================================================
/**
 * _dea_ipm_slacks()
 *
 * Stage 2: Given theta* from Stage 1, maximize total slacks.
 *
 * Input-oriented:
 *   max sum(s_in) + sum(s_out)
 *   s.t. X*lambda + s_in  = theta* * x_k    (inputs)
 *        Y*lambda - s_out = y_k              (outputs)
 *        [convexity constraint]
 *        lambda, s_in, s_out >= 0
 *
 * Output-oriented:
 *   max sum(s_in) + sum(s_out)
 *   s.t. X*lambda + s_in  = x_k             (inputs)
 *        Y*lambda - s_out = phi* * y_k       (outputs)
 *        [convexity constraint]
 *        lambda, s_in, s_out >= 0
 *
 * Variables: [lambda(n), s_in(m), s_out(s)]
 *
 * Modifies lam_k and slacks in place only on success; returns the native
 * error code (3 for nonconvergence without a nonzero native error code).
 */
real scalar _dea_ipm_slacks(real matrix X,
                      real matrix Y,
                      real scalar k,
                      string scalar orientation,
                      string scalar rts,
                      real scalar theta_star,
                      real rowvector lam_k,
                      real rowvector slacks_in,
                      real rowvector slacks_out,
                      real scalar tol, | real rowvector audit,
                      real scalar radial_override)
{
    real scalar n, m, s, total_vars, j
    real rowvector c, lb, ub, sol
    real matrix A_eq, A_ineq
    real colvector b_eq, b_ineq, x0, y0
    real scalar rc, status
    class LinearProgram scalar S2
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    total_vars = n + m + s  // lambda(n) + s_in(m) + s_out(s)
    
    x0 = X[k, ]'
    y0 = Y[k, ]'
    
    // Objective: max sum(s_in) + sum(s_out)
    //         = min -sum(s_in) - sum(s_out)
    c = (J(1, n, 0), -J(1, m, 1), -J(1, s, 1))
    
    // Equality constraints
    if (orientation == "input") {
        // X'*lambda + I*s_in = theta* * x_k
        // Y'*lambda - I*s_out = y_k
        A_eq = (X', I(m), J(m, s, 0)) \
               (Y', J(s, m, 0), -I(s))
        b_eq = ((args()>=12 ? radial_override : theta_star) * x0) \ y0
    }
    else {
        // X'*lambda + I*s_in = x_k
        // Y'*lambda - I*s_out = (1/theta_star) * y_k  [phi* = 1/theta_star]
        real scalar phi_star
        phi_star = (args()>=12 ? radial_override : 1/theta_star)
        
        A_eq = (X', I(m), J(m, s, 0)) \
               (Y', J(s, m, 0), -I(s))
        b_eq = x0 \ (phi_star * y0)
    }
    
    // Preserve the stage-1 technology: CRS has no lambda-sum restriction.
    A_ineq = J(0, total_vars, 0)
    b_ineq = J(0, 1, 0)
    if (rts == "vrs") {
        A_eq = A_eq \ (J(1, n, 1), J(1, m + s, 0))
        b_eq = b_eq \ 1
    }
    else if (rts == "drs") {
        A_ineq = (J(1, n, 1), J(1, m + s, 0))
        b_ineq = 1
    }
    else if (rts == "irs") {
        A_ineq = (J(1, n, -1), J(1, m + s, 0))
        b_ineq = -1
    }
    
    // Bounds: all >= 0
    lb = J(1, total_vars, 0)
    ub = J(1, total_vars, .)
    
    // Solve
    S2 = LinearProgram()
    S2.setMaxOrMin("min")
    S2.setCoefficients(c)
    S2.setEquality(A_eq, b_eq)
    if (rows(A_ineq)>0) S2.setInequality(A_ineq, b_ineq)
    S2.setBounds(lb, ub)
    S2.setTrace("off")
    S2.setTol(tol)
    
    rc = S2.optimize()
    sol=J(1,total_vars,.)
    status=(S2.errorcode()!=0 ? S2.errorcode() : 3)
    if (S2.converged() && S2.errorcode() == 0) {
        sol = S2.parameters()
        status=(_dea_lpcheck_ok(_dea_lpcheck(sol,A_eq,b_eq,A_ineq,b_ineq,lb,ub),tol)
            & c*sol'<. ? 0 : 90)
    }
    status=_dea_ipm_finish(status,S2.errorcode(),tol,S2.iterations(),
        c,A_eq,b_eq,A_ineq,b_ineq,lb,ub,sol,audit)
    if (status==0) {
        lam_k=sol[1..n]
        slacks_in=sol[(n+1)..(n+m)]
        slacks_out=sol[(n+m+1)..(n+m+s)]
    }
    else {
        lam_k=J(1,n,.)
        slacks_in=J(1,m,.)
        slacks_out=J(1,s,.)
    }
    return(status)
}


// ============================================================================
// 3. DUAL FALLBACK SOLVERS
// ============================================================================

/**
 * _dea_ipm_dual_input()
 *
 * Multiplier form, input-oriented:
 *   max u'y0 [+ u0]  s.t. v'x0 = 1, u'Y - v'X [+ u0] <= 0
 *   theta = u'y0 [+ u0]
 *
 * Supports CRS/VRS/DRS/IRS:
 *   CRS: no u0 term
 *   VRS: u0 free (unrestricted)
 *   DRS: u0 <= 0  (non-positive)
 *   IRS: u0 >= 0  (non-negative)
 */
real scalar _dea_ipm_dual_input(real matrix X,
                                real matrix Y,
                                real scalar k,
                                string scalar rts, | real scalar tol, real rowvector audit)
{
    real scalar n, m_in, s_out, total_vars, eff, status
    real rowvector x0, y0, c, lb, ub, sol
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    class LinearProgram scalar lp
    
    if (args()<5) tol=1e-8
    n = rows(X)
    m_in = cols(X)
    s_out = cols(Y)
    x0 = X[k, .]
    y0 = Y[k, .]
    
    lp = LinearProgram()
    lp.setTrace("off")
    lp.setTol(tol)
    
    if (rts == "crs") {
        // Variables: [u(s), v(m)]
        total_vars = s_out + m_in
        c = (-y0, J(1, m_in, 0))        // min -u'y0 = max u'y0
        A_eq = (J(1, s_out, 0), x0)     // v'x0 = 1
        A_ineq = (Y, -X)                // u'Y - v'X <= 0
        lb = J(1, total_vars, 0)
    }
    else {
        // Variables: [u(s), v(m), u0]
        total_vars = s_out + m_in + 1
        c = (-y0, J(1, m_in, 0), -1)    // min -(u'y0 + u0)
        A_eq = (J(1, s_out, 0), x0, 0)  // v'x0 = 1
        A_ineq = (Y, -X, J(n, 1, 1))   // u'Y - v'X + u0 <= 0
        
        // Bounds on u0 depend on RTS
        if (rts == "vrs") {
            lb = (J(1, s_out + m_in, 0), .)  // u0 free
        }
        else if (rts == "drs") {
            lb = (J(1, s_out + m_in, 0), .)  // u0 free (but u0 <= 0)
            // DRS dual: u0 <= 0, implemented via upper bound
        }
        else if (rts == "irs") {
            lb = (J(1, s_out + m_in, 0), 0)  // u0 >= 0
        }
    }
    b_eq = 1
    b_ineq = J(n, 1, 0)
    ub = J(1, total_vars, .)
    
    // DRS: u0 <= 0
    if (rts == "drs") {
        ub[total_vars] = 0
    }
    
    lp.setCoefficients(c)
    lp.setMaxOrMin("min")
    lp.setEquality(A_eq, b_eq)
    lp.setInequality(A_ineq, b_ineq)
    lp.setBounds(lb, ub)
    
    sol=J(1,total_vars,.)
    (void) lp.optimize()
    
    if (lp.converged() == 1 & lp.errorcode() == 0) {
        sol=lp.parameters()
        eff=(_dea_lpcheck_ok(_dea_lpcheck(sol,A_eq,b_eq,A_ineq,b_ineq,lb,ub),tol) ? -c*sol' : .)
    }
    else {
        eff = .
    }
    
    status=(eff<. ? 0 : (lp.errorcode()!=0 ? lp.errorcode() :
        (lp.converged() ? 90 : 3)))
    status=_dea_ipm_finish(status,lp.errorcode(),tol,lp.iterations(),
        c,A_eq,b_eq,A_ineq,b_ineq,lb,ub,sol,audit)
    if (status==0 & eff>=.) {
        eff=-c*sol'
        if (eff>=.) audit[1]=90
    }
    if (status!=0) eff=.
    return(eff)
}


/**
 * _dea_ipm_dual_output()
 *
 * Multiplier form, output-oriented:
 *   min v'x0 [- u0]  s.t. u'y0 = 1, u'Y - v'X [+ u0] <= 0
 *   phi = v'x0 [- u0]   →  efficiency = 1/phi
 *
 * Returns efficiency in [0,1] scale.
 */
real scalar _dea_ipm_dual_output(real matrix X,
                                  real matrix Y,
                                  real scalar k,
                                  string scalar rts, | real scalar tol, real rowvector audit)
{
    real scalar n, m_in, s_out, total_vars, eff, status
    real rowvector x0, y0, c, lb, ub, sol
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    class LinearProgram scalar lp
    
    if (args()<5) tol=1e-8
    n = rows(X)
    m_in = cols(X)
    s_out = cols(Y)
    x0 = X[k, .]
    y0 = Y[k, .]
    
    lp = LinearProgram()
    lp.setTrace("off")
    lp.setTol(tol)
    
    if (rts == "crs") {
        total_vars = s_out + m_in
        c = (J(1, s_out, 0), x0)        // min v'x0
        A_eq = (y0, J(1, m_in, 0))      // u'y0 = 1
        A_ineq = (Y, -X)                // u'Y - v'X <= 0
        lb = J(1, total_vars, 0)
    }
    else {
        total_vars = s_out + m_in + 1
        c = (J(1, s_out, 0), x0, -1)    // min v'x0 - u0
        A_eq = (y0, J(1, m_in, 0), 0)   // u'y0 = 1
        A_ineq = (Y, -X, J(n, 1, 1))   // u'Y - v'X + u0 <= 0
        
        if (rts == "vrs") {
            lb = (J(1, s_out + m_in, 0), .)  // u0 free
        }
        else if (rts == "drs") {
            lb = (J(1, s_out + m_in, 0), .)  // u0 free, capped at 0
        }
        else if (rts == "irs") {
            lb = (J(1, s_out + m_in, 0), 0)  // u0 >= 0
        }
    }
    b_eq = 1
    b_ineq = J(n, 1, 0)
    ub = J(1, total_vars, .)
    
    if (rts == "drs") {
        ub[total_vars] = 0  // u0 <= 0
    }
    
    lp.setCoefficients(c)
    lp.setMaxOrMin("min")
    lp.setEquality(A_eq, b_eq)
    lp.setInequality(A_ineq, b_ineq)
    lp.setBounds(lb, ub)
    
    sol=J(1,total_vars,.)
    (void) lp.optimize()
    
    if (lp.converged() == 1 & lp.errorcode() == 0) {
        sol=lp.parameters()
        eff=(_dea_lpcheck_ok(_dea_lpcheck(sol,A_eq,b_eq,A_ineq,b_ineq,lb,ub),tol) ? c*sol' : .)
        // The dual optimum equals phi, not 1/phi.
        eff=(eff>0 & eff<. ? 1/eff : .)
    }
    else {
        eff = .
    }
    
    status=(eff<. ? 0 : (lp.errorcode()!=0 ? lp.errorcode() :
        (lp.converged() ? 90 : 3)))
    status=_dea_ipm_finish(status,lp.errorcode(),tol,lp.iterations(),
        c,A_eq,b_eq,A_ineq,b_ineq,lb,ub,sol,audit)
    if (status==0 & eff>=.) {
        eff=c*sol'
        eff=(eff>0 & eff<. ? 1/eff : .)
        if (eff>=.) audit[1]=90
    }
    if (status!=0) eff=.
    return(eff)
}


// ============================================================================
// 4. VRS FRONTIER ANALYSIS
// ============================================================================
/**
 * _dea_ipm_vrsfrontier()
 *
 * Runs CRS → VRS → DRS sequentially for VRS Frontier analysis.
 * Computes:
 *   CRS_TE  = CRS efficiency
 *   VRS_TE  = VRS efficiency
 *   SCALE   = CRS_TE / VRS_TE
 *   RTS     = based on DRS_TE comparison:
 *             if CRS_TE == VRS_TE → CRS (0)
 *             else if DRS_TE == VRS_TE → DRS (-1)
 *             else → IRS (1)
 *
 * Stores results in 5 Stata variables: crs_te, vrs_te, drs_te, scale, rts_type
 *
 * Alternative RTS method (sum of lambda):
 *   sum(lambda_CRS) < 1 → IRS
 *   sum(lambda_CRS) > 1 → DRS
 *   sum(lambda_CRS) = 1 → CRS
 */
void _dea_ipm_vrsfrontier(string scalar inputs,
                           string scalar outputs,
                           string scalar crs_var,
                           string scalar vrs_var,
                           string scalar drs_var,
                           string scalar scale_var,
                           string scalar rts_var,
                           string scalar lambda_sum_var,
                           string scalar touse,
                           string scalar orientation,
                           real scalar tol, | string scalar trace_mat)
{
    real matrix X, Y, T
    real colvector obs, x0, y0, BE, BU
    real rowvector c, lb, ub, sol, scores, lam, p_trace, d_trace
    real matrix AE, AU
    real scalar n,m,q,k,j,t,rc,eff,scale_eff,rts_type,sum_lambda,den,status
    string scalar rts
    class LinearProgram scalar S
    if (args()<12) trace_mat=""
    X=st_data(.,tokens(inputs),touse);Y=st_data(.,tokens(outputs),touse)
    obs=selectindex(st_data(.,touse));n=rows(X);m=cols(X);q=cols(Y)
    // Match the main radial path's units, convergence and acceptance settings.
    for (j=1;j<=m;j++) {
        den=max(abs(X[.,j]));if (den>0 & den<.) X[.,j]=X[.,j]/den
    }
    for (j=1;j<=q;j++) {
        den=max(abs(Y[.,j]));if (den>0 & den<.) Y[.,j]=Y[.,j]/den
    }
    T=J((trace_mat!="" ? n : 0),150,.)
    for (k=1;k<=n;k++) {
        x0=X[k,.]';y0=Y[k,.]';scores=J(1,3,.);sum_lambda=.
        for (t=1;t<=3;t++) {
            rts=(t==1 ? "crs" : (t==2 ? "vrs" : "drs"))
            if (orientation=="input") {
                c=(J(1,n,0),1);lb=J(1,n+1,0)
                AU=(X',-x0) \ (-Y',J(q,1,0));BU=J(m,1,0) \ (-y0)
            }
            else {
                c=(J(1,n,0),-1);lb=(J(1,n,0),1)
                AU=(X',J(m,1,0)) \ (-Y',y0);BU=x0 \ J(q,1,0)
            }
            ub=J(1,n+1,.);AE=J(0,n+1,0);BE=J(0,1,0)
            if (rts=="vrs") { AE=(J(1,n,1),0);BE=1; }
            if (rts=="drs") { AU=AU\(J(1,n,1),0);BU=BU\1; }
            S=LinearProgram();S.setMaxOrMin("min");S.setCoefficients(c)
            S.setInequality(AU,BU)
            if (rows(AE)>0) S.setEquality(AE,BE)
            S.setBounds(lb,ub);S.setTrace("off");S.setTol(tol)
            sol=J(1,n+1,.);p_trace=d_trace=J(1,25,.)
            rc=S.optimize();eff=.
            if (S.converged() & S.errorcode()==0) {
                sol=S.parameters()
                if (_dea_lpcheck_ok(_dea_lpcheck(sol,AE,BE,AU,BU,lb,ub),tol)) {
                    eff=c*sol'
                    if (orientation=="output") eff=(eff<0 ? -1/eff : .)
                    if (eff<. & t==1) sum_lambda=sum(sol[1..n])
                }
            }
            status=(eff<. ? 0 : (S.errorcode()!=0 ? S.errorcode() :
                (S.converged() ? 90 : 3)))
            status=_dea_ipm_finish(status,S.errorcode(),tol,S.iterations(),
                c,AE,BE,AU,BU,lb,ub,sol,p_trace)
            if (eff>=. & status==0) {
                eff=c*sol'
                if (orientation=="output") eff=(eff<0 ? -1/eff : .)
                if (eff<. & t==1) sum_lambda=sum(sol[1..n])
                if (eff>=.) { status=90;p_trace[1]=90; }
            }
            if (eff>=. & !any(p_trace[(2,8,14,20)]:==2)) {
                eff=(orientation=="input" ? _dea_ipm_dual_input(X,Y,k,rts,tol,d_trace) :
                    _dea_ipm_dual_output(X,Y,k,rts,tol,d_trace))
            }
            if (trace_mat!="") T[k,((t-1)*50+1)..(t*50)]=(p_trace,d_trace)
            scores[t]=eff
        }
        scale_eff=(scores[1]<. & scores[2]>0 & scores[2]<. ? scores[1]/scores[2] : .)
        rts_type=.
        if (scores[1]<. & scores[2]<.) {
            if (abs(scores[1]-scores[2])<tol) rts_type=0
            else if (scores[3]<. & abs(scores[3]-scores[2])<tol) rts_type=-1
            else if (scores[3]<.) rts_type=1
        }
        st_store(obs[k],st_varindex(crs_var),scores[1])
        st_store(obs[k],st_varindex(vrs_var),scores[2])
        st_store(obs[k],st_varindex(drs_var),scores[3])
        st_store(obs[k],st_varindex(scale_var),scale_eff)
        st_store(obs[k],st_varindex(rts_var),rts_type)
        if (lambda_sum_var!="") st_store(obs[k],st_varindex(lambda_sum_var),sum_lambda)
    }
    _dea_ipm_trace_store(trace_mat,T,("crs_p","crs_d","vrs_p","vrs_d","drs_p","drs_d"))
}


// ============================================================================
// 5. UTILITY: COUNT RESULTS
// ============================================================================
void _dea_ipm_count_results(string scalar eff_var,
                             string scalar touse,
                             real scalar n)
{
    real colvector eff_data
    real scalar n_conv
    
    eff_data = st_data(., eff_var, touse)
    n_conv = rows(eff_data) - colmissing(eff_data)
    st_numscalar("r(n_converged)", n_conv)
    st_numscalar("r(n_failed)", n - n_conv)
}


// ============================================================================
// 6. SUPER-EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_ipm_supereff()
 *
 * Computes super-efficiency by excluding the evaluated DMU from its reference set.
 * Input-oriented:  min theta  s.t. X_{-o}*lambda <= theta*x0, Y_{-o}*lambda >= y0
 *   theta > 1 → super-efficient (stored as theta)
 * Output-oriented: max phi    s.t. X_{-o}*lambda <= x0,       Y_{-o}*lambda >= phi*y0
 *   phi < 1 → super-efficient (stored as 1/phi > 1)
 *
 * Arguments:
 *   lambda_var   - string: lambda variable names (n, space-separated), "" if none
 *                  LP solves over n-1 lambdas (excluding DMU k);
 *                  results are mapped back to n positions with lambda_k = 0.
 */
void _dea_ipm_supereff(string scalar inputs,
                        string scalar outputs,
                        string scalar eff_var,
                        string scalar lambda_var,
                        string scalar touse,
                        string scalar orientation,
                        string scalar rts,
    | real scalar tol)
{
    real matrix X, Y, X_ref, Y_ref
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, sol, lam_ref, lam_k
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc, eff_score, n_ref, j
    real scalar phi_raw, has_lambda
    class LinearProgram scalar S
    
    // Lambda variable handling
    string rowvector lam_varlist
    has_lambda = (lambda_var != "")
    if (has_lambda) lam_varlist = tokens(lambda_var)
    
    if (args()<8) tol=1e-7

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        eff_score = .
        lam_k = J(1, n, 0)  // Initialize: all zeros (k-th is always 0)
        
        // Create reference set excluding DMU k
        n_ref = n - 1
        if (k == 1) {
            X_ref = X[2..n, ]
            Y_ref = Y[2..n, ]
        }
        else if (k == n) {
            X_ref = X[1..n-1, ]
            Y_ref = Y[1..n-1, ]
        }
        else {
            X_ref = X[1..k-1, ] \ X[k+1..n, ]
            Y_ref = Y[1..k-1, ] \ Y[k+1..n, ]
        }
        
        if (orientation == "input") {
            // min theta s.t. X_ref'*lambda <= theta*x0, Y_ref'*lambda >= y0
            // Variables: [lambda(1..n_ref), theta]
            c = J(1, n_ref, 0), 1
            A_ineq = (X_ref', -x0) \ (-Y_ref', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)
            lb = J(1, n_ref+1, 0)
            ub = J(1, n_ref+1, .)
            
            A_eq = J(0, n_ref+1, .)
            b_eq = J(0, 1, .)
            
            if (rts == "vrs") {
                A_eq = (J(1, n_ref, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n_ref, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n_ref, -1), 0)
                b_ineq = b_ineq \ (-1)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) {
                S.setEquality(A_eq, b_eq)
            }
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = _dea_nonte_optimize(S,tol)
            
            if (rc == 0) {
                sol = S.parameters()
                eff_score = S.value()
                lam_ref = sol[1..n_ref]
                
                // Map n-1 lambdas back to n positions (k-th = 0)
                if (k == 1) {
                    lam_k[2..n] = lam_ref
                }
                else if (k == n) {
                    lam_k[1..n-1] = lam_ref
                }
                else {
                    lam_k[1..k-1] = lam_ref[1..k-1]
                    lam_k[k+1..n] = lam_ref[k..n_ref]
                }
            }
        }
        else {
            // Output-oriented: max phi
            // Super-efficiency allows phi < 1 (super-efficient DMUs)
            // Variables: [lambda(1..n_ref), phi]
            c = J(1, n_ref, 0), -1
            A_ineq = (X_ref', J(m, 1, 0)) \ (-Y_ref', y0)
            b_ineq = x0 \ J(s, 1, 0)
            lb = J(1, n_ref, 0), 1e-8
            ub = J(1, n_ref+1, .)
            
            A_eq = J(0, n_ref+1, .)
            b_eq = J(0, 1, .)
            
            if (rts == "vrs") {
                A_eq = (J(1, n_ref, 1), 0)
                b_eq = 1
            }
            else if (rts == "drs") {
                A_ineq = A_ineq \ (J(1, n_ref, 1), 0)
                b_ineq = b_ineq \ 1
            }
            else if (rts == "irs") {
                A_ineq = A_ineq \ (J(1, n_ref, -1), 0)
                b_ineq = b_ineq \ (-1)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) {
                S.setEquality(A_eq, b_eq)
            }
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = _dea_nonte_optimize(S,tol)
            
            if (rc == 0) {
                sol = S.parameters()
                phi_raw = -S.value()
                eff_score = 1 / phi_raw
                lam_ref = sol[1..n_ref]
                
                // Map n-1 lambdas back to n positions (k-th = 0)
                if (k == 1) {
                    lam_k[2..n] = lam_ref
                }
                else if (k == n) {
                    lam_k[1..n-1] = lam_ref
                }
                else {
                    lam_k[1..k-1] = lam_ref[1..k-1]
                    lam_k[k+1..n] = lam_ref[k..n_ref]
                }
            }
        }
        
        // Store efficiency score
        st_store(obs_k, st_varindex(eff_var), eff_score)
        
        // Store lambda weights (mapped to original n DMU positions)
        if (has_lambda) {
            for (j = 1; j <= n; j++) {
                st_store(obs_k, st_varindex(lam_varlist[j]), lam_k[j])
            }
        }
    }
}


// ============================================================================
// 7. DIRECTIONAL DISTANCE FUNCTION (DDF)
// ============================================================================
/**
 * _dea_ipm_ddf()
 *
 * Solves the Directional Distance Function model.
 *
 * max beta s.t. X*lambda <= x0 - beta*gx
 *               Y*lambda >= y0 + beta*gy
 *               [convexity constraint]
 *               lambda >= 0, beta >= 0
 */
void _dea_ipm_ddf(string scalar inputs,
                  string scalar outputs,
                  string scalar eff_var,
                  string scalar gx_var,
                  string scalar gy_var,
                  string scalar touse,
                  string scalar rts,
    | real scalar tol)
{
    real matrix X, Y, GX, GY
    real colvector obs_indices, x0, y0, gx0, gy0
    real rowvector c, lb, ub, sol
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc, beta
    class LinearProgram scalar S
    
    if (args()<8) tol=1e-7

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    // Load direction vectors (or use observed values as default)
    if (gx_var != "") {
        GX = st_data(., tokens(gx_var), touse)
    }
    else {
        GX = X  // Default: direction = observed input values
    }
    
    if (gy_var != "") {
        GY = st_data(., tokens(gy_var), touse)
    }
    else {
        GY = Y  // Default: direction = observed output values
    }
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        gx0 = GX[k, ]'
        gy0 = GY[k, ]'
        beta = .
        
        // max beta => min -beta
        // Variables: [lambda(1..n), beta]
        // Constraints:
        //   X'*lambda + beta*gx0 <= x0  =>  X'*lambda + beta*gx0 - x0 <= 0
        //   Y'*lambda - beta*gy0 >= y0  => -Y'*lambda + beta*gy0 <= -y0
        
        c = J(1, n, 0), -1   // min(-beta) = max(beta)
        
        A_ineq = (X', gx0) \ (-Y', gy0)
        b_ineq = x0 \ (-y0)
        
        lb = J(1, n+1, 0)
        ub = J(1, n+1, .)
        
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
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) {
            S.setEquality(A_eq, b_eq)
        }
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = _dea_nonte_optimize(S,tol)
        
        if (rc == 0) {
            sol = S.parameters()
            beta = sol[n+1]
        }
        
        st_store(obs_k, st_varindex(eff_var), beta)
    }
}


// ============================================================================
// 7b. DDF WITH UNDESIRABLE (BAD) OUTPUTS — Wang et al. (2022)
// ============================================================================
/**
 * _dea_ipm_ddf_undesirable()
 *
 * Directional Distance Function with undesirable outputs (bad-output equality).
 * VRS/NDRS equality technology need not be closed under joint output contraction.
 * Based on Chung et al. (1997), Färe et al. (2007), Wang et al. (2022).
 *
 * max β
 * s.t.  X'λ + β·gx ≤ x₀     (inputs, if gx ≠ 0)
 *       Y'λ - β·gy ≥ y₀     (good outputs: (-Y',+gy)·[λ;β] ≤ -y₀)
 *       B'λ + β·gb = b₀     (bad outputs: equality, no reduction slack)
 *       [Σλ = 1  VRS]
 *       λ ≥ 0, β ≥ 0
 *
 * Direction vectors (direction()):
 *   "yb"  (default): gx=0,  gy=y₀, gb=b₀  — adjust good/bad outputs only
 *   "xyb":           gx=x₀, gy=y₀, gb=b₀  — adjust all three
 *   "x":             gx=x₀, gy=0,  gb=0   — input-oriented (bad outputs fixed)
 *
 * Good outputs expand in direction gy:
 *   Y'λ - β·gy ≥ y₀, equivalently Y'λ ≥ y₀ + β·gy.
 */
void _dea_ipm_ddf_undesirable(string scalar inputs,
                               string scalar outputs,
                               string scalar bad_outputs,
                               string scalar eff_var,
                               string scalar direction,
                               string scalar touse,
                               string scalar rts,
    | real scalar tol)
{
    real matrix X, Y, B
    real colvector obs_indices, x0, y0, b0, gx0, gy0, gb0
    real rowvector c, lb, ub, sol, xscale, yscale, bscale
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, q, k, obs_k, rc, beta
    class LinearProgram scalar S

    if (args()<8) tol=1e-7

    X = st_data(., tokens(inputs),      touse)
    Y = st_data(., tokens(outputs),     touse)
    B = st_data(., tokens(bad_outputs), touse)
    obs_indices = selectindex(st_data(., touse))

    n = rows(X)
    m = cols(X)    // number of inputs
    s = cols(Y)    // number of good outputs
    q = cols(B)    // number of bad outputs

    // Observed directions scale with their variables. Positive column scaling
    // preserves beta/lambda and avoids a dependence on measurement units.
    xscale = colmax(abs(X)); xscale = xscale :+ (xscale:==0)
    yscale = colmax(abs(Y)); yscale = yscale :+ (yscale:==0)
    bscale = colmax(abs(B)); bscale = bscale :+ (bscale:==0)
    X = X :/ xscale
    Y = Y :/ yscale
    B = B :/ bscale

    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0  = X[k, .]'
        y0  = Y[k, .]'
        b0  = B[k, .]'
        beta = .

        // --- Direction vectors ---
        if (direction == "xyb") {
            gx0 = x0
            gy0 = y0
            gb0 = b0
        }
        else if (direction == "x") {
            gx0 = x0
            gy0 = J(s, 1, 0)
            gb0 = J(q, 1, 0)
        }
        else {  // "yb" (default)
            gx0 = J(m, 1, 0)
            gy0 = y0
            gb0 = b0
        }

        // Variables: [λ(n), β];  objective: min -β (= max β)
        c  = (J(1, n, 0), -1)
        lb = J(1, n+1, 0)
        ub = J(1, n+1, .)

        // Inequality constraints
        // X'λ + β·gx ≤ x₀ :  (X', gx0) · [λ;β] ≤ x₀
        // Y'λ - β·gy ≥ y₀ :  (-Y', +gy0) · [λ;β] ≤ -y₀
        A_ineq = (X', gx0) \ (-Y', gy0)
        b_ineq = x0 \ (-y0)

        // Equality constraints: bad outputs + (optional) convexity
        // B'λ + β·gb = b₀ :  (B', gb0) · [λ;β] = b₀
        A_eq = (B', gb0)
        b_eq = b0

        if (rts == "vrs") {
            A_eq = A_eq \ (J(1, n, 1), 0)
            b_eq = b_eq \ 1
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ (J(1, n, 1), 0)
            b_ineq = b_ineq \ 1
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ (J(1, n, -1), 0)
            b_ineq = b_ineq \ (-1)
        }

        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A_ineq, b_ineq)
        S.setEquality(A_eq, b_eq)
        S.setBounds(lb, ub)
        S.setTrace("off")

        rc = _dea_nonte_optimize(S,tol,1e-12)

        if (rc == 0) {
            sol  = S.parameters()
            beta = sol[n+1]
        }

        st_store(obs_k, st_varindex(eff_var), beta)
    }
}


// ============================================================================
// 8. ADDITIVE MODEL
// ============================================================================
/**
 * _dea_ipm_additive()
 *
 * Maximizes the sum of input and output slacks.
 *
 * max sum(s-) + sum(s+)
 * s.t. X*lambda + s- = x0
 *      Y*lambda - s+ = y0
 *      [convexity constraint]
 *      lambda, s-, s+ >= 0
 */
void _dea_ipm_additive(string scalar inputs,
                        string scalar outputs,
                        string scalar eff_var,
                        string scalar slack_in_var,
                        string scalar slack_out_var,
                        string scalar touse,
                        string scalar rts,
    | real scalar tol)
{
    real matrix X, Y
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, sol
    real matrix A_eq
    real colvector b_eq
    real matrix A_ineq
    real colvector b_ineq
    real scalar n, m, s, k, obs_k, rc, add_score
    real rowvector slack_in, slack_out
    class LinearProgram scalar S
    
    if (args()<8) tol=1e-7

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    string rowvector si_varlist, so_varlist
    real scalar has_slacks, j
    
    has_slacks = (slack_in_var != "" & slack_out_var != "")
    if (has_slacks) {
        si_varlist = tokens(slack_in_var)
        so_varlist = tokens(slack_out_var)
    }
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        add_score = .
        slack_in = J(1, m, 0)
        slack_out = J(1, s, 0)
        
        // Variables: [lambda(n), s-(m), s+(s)]
        // Total: n + m + s variables
        // Objective: max sum(s-) + sum(s+) => min -sum(s-) - sum(s+)
        
        c = J(1, n, 0), J(1, m, -1), J(1, s, -1)
        
        // Equality constraints:
        // X'*lambda + I*s- = x0  => X'*lambda + I*s- - x0 = 0
        // Y'*lambda - I*s+ = y0  => Y'*lambda - I*s+ - y0 = 0
        
        A_eq = (X', I(m), J(m, s, 0)) \ (Y', J(s, m, 0), -I(s))
        b_eq = x0 \ y0
        
        // RTS constraints
        A_ineq = J(0, n+m+s, .)
        b_ineq = J(0, 1, .)
        
        if (rts == "vrs") {
            A_eq = A_eq \ (J(1, n, 1), J(1, m+s, 0))
            b_eq = b_eq \ 1
        }
        else if (rts == "drs") {
            A_ineq = (J(1, n, 1), J(1, m+s, 0))
            b_ineq = 1
        }
        else if (rts == "irs") {
            A_ineq = (J(1, n, -1), J(1, m+s, 0))
            b_ineq = -1
        }
        
        lb = J(1, n+m+s, 0)
        ub = J(1, n+m+s, .)
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setEquality(A_eq, b_eq)
        if (rows(A_ineq) > 0) {
            S.setInequality(A_ineq, b_ineq)
        }
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = _dea_nonte_optimize(S,tol)
        
        if (rc == 0) {
            sol = S.parameters()
            slack_in = sol[n+1..n+m]
            slack_out = sol[n+m+1..n+m+s]
            add_score = sum(slack_in) + sum(slack_out)
        }
        
        // Store efficiency (0 if efficient, >0 if inefficient)
        st_store(obs_k, st_varindex(eff_var), add_score)
        
        // Store slacks if requested
        if (has_slacks) {
            for (j = 1; j <= m; j++) {
                st_store(obs_k, st_varindex(si_varlist[j]), slack_in[j])
            }
            for (j = 1; j <= s; j++) {
                st_store(obs_k, st_varindex(so_varlist[j]), slack_out[j])
            }
        }
    }
}


// ============================================================================
// 9. COST EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_ipm_cost()
 *
 * Computes cost efficiency given input prices.
 *
 * min w'*x
 * s.t. X*lambda <= x
 *      Y*lambda >= y0
 *      [convexity constraint]
 *      lambda, x >= 0
 *
 * Cost Efficiency = min_cost / observed_cost = (w'*x*) / (w'*x0)
 */
void _dea_ipm_cost(string scalar inputs,
                   string scalar outputs,
                   string scalar eff_var,
                   string scalar w_var,
                   string scalar touse,
                   string scalar rts,
    | real scalar tol)
{
    real matrix X, Y, W
    real colvector obs_indices, x0, y0, w0
    real rowvector c, lb, ub, sol
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc
    real scalar min_cost, obs_cost, cost_eff
    class LinearProgram scalar S
    
    if (args()<7) tol=1e-7

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    W = st_data(., tokens(w_var), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        w0 = W[k, ]'
        cost_eff = .
        
        // Observed cost
        obs_cost = sum(w0 :* x0)
        
        // Variables: [lambda(n), x(m)]
        // Objective: min w'*x
        // Constraints:
        //   X'*lambda - x <= 0
        //   -Y'*lambda <= -y0
        //   [convexity]
        
        c = J(1, n, 0), w0'
        
        A_ineq = (X', -I(m)) \ (-Y', J(s, m, 0))
        b_ineq = J(m, 1, 0) \ (-y0)
        
        lb = J(1, n+m, 0)
        ub = J(1, n+m, .)
        
        A_eq = J(0, n+m, .)
        b_eq = J(0, 1, .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), J(1, m, 0))
            b_eq = 1
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ (J(1, n, 1), J(1, m, 0))
            b_ineq = b_ineq \ 1
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ (J(1, n, -1), J(1, m, 0))
            b_ineq = b_ineq \ (-1)
        }
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) {
            S.setEquality(A_eq, b_eq)
        }
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = _dea_nonte_optimize(S,tol)
        
        if (rc == 0) {
            min_cost = S.value()
            if (obs_cost > 0) {
                cost_eff = min_cost / obs_cost
            }
        }
        
        st_store(obs_k, st_varindex(eff_var), cost_eff)
    }
}


// ============================================================================
// 10. REVENUE EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_ipm_revenue()
 *
 * Computes revenue efficiency given output prices.
 *
 * max p'*y
 * s.t. X*lambda <= x0
 *      Y*lambda >= y
 *      [convexity constraint]
 *      lambda, y >= 0
 *
 * Revenue Efficiency = observed_revenue / max_revenue = (p'*y0) / (p'*y*)
 */
void _dea_ipm_revenue(string scalar inputs,
                      string scalar outputs,
                      string scalar eff_var,
                      string scalar p_var,
                      string scalar touse,
                      string scalar rts,
    | real scalar tol)
{
    real matrix X, Y, P
    real colvector obs_indices, x0, y0, p0
    real rowvector c, lb, ub, sol
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc
    real scalar max_rev, obs_rev, rev_eff
    class LinearProgram scalar S
    
    if (args()<7) tol=1e-7

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    P = st_data(., tokens(p_var), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        p0 = P[k, ]'
        rev_eff = .
        
        // Observed revenue
        obs_rev = sum(p0 :* y0)
        
        // Variables: [lambda(n), y(s)]
        // Objective: max p'*y => min -p'*y
        // Constraints:
        //   X'*lambda <= x0
        //   Y'*lambda - y >= 0 => -Y'*lambda + y <= 0
        //   [convexity]
        
        c = J(1, n, 0), -p0'
        
        A_ineq = (X', J(m, s, 0)) \ (-Y', I(s))
        b_ineq = x0 \ J(s, 1, 0)
        
        lb = J(1, n+s, 0)
        ub = J(1, n+s, .)
        
        A_eq = J(0, n+s, .)
        b_eq = J(0, 1, .)
        
        if (rts == "vrs") {
            A_eq = (J(1, n, 1), J(1, s, 0))
            b_eq = 1
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ (J(1, n, 1), J(1, s, 0))
            b_ineq = b_ineq \ 1
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ (J(1, n, -1), J(1, s, 0))
            b_ineq = b_ineq \ (-1)
        }
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) {
            S.setEquality(A_eq, b_eq)
        }
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = _dea_nonte_optimize(S,tol)
        
        if (rc == 0) {
            max_rev = -S.value()
            if (max_rev > 0) {
                rev_eff = obs_rev / max_rev
            }
        }
        
        st_store(obs_k, st_varindex(eff_var), rev_eff)
    }
}


// ============================================================================
// 10a. PROFIT EFFICIENCY MODEL
// ============================================================================
/**
 * _dea_ipm_profit()
 *
 * Computes profit efficiency (Cooper, Seiford & Tone, Ch.8).
 *
 * [Profit] eq. 8.8:
 *   py* - cx* = max  p'(Yλ) - c'(Xλ)
 *               λ
 *   s.t. x = Xλ ≤ x_o       (inputs bounded by observed)
 *        y = Yλ ≥ y_o        (outputs meet observed)
 *        [convexity]
 *        λ ≥ 0
 *
 * Variables: λ(n) only — x and y are determined by λ.
 * The LP is bounded under CRS because Xλ ≤ x_o constrains λ.
 *
 * Profit Efficiency = observed_profit / max_profit  (eq. 8.9)
 *   where observed_profit = p'y_o - c'x_o
 *         max_profit = p'(Yλ*) - c'(Xλ*)
 */
void _dea_ipm_profit(string scalar inputs,
                     string scalar outputs,
                     string scalar eff_var,
                     string scalar w_var,
                     string scalar p_var,
                     string scalar touse,
                     string scalar rts,
    | real scalar tol)
{
    real matrix X, Y, W, P
    real colvector obs_indices, x0, y0, w0, p0
    real rowvector c, lb, ub
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc
    real scalar max_profit, obs_profit, profit_eff
    class LinearProgram scalar S
    
    if (args()<8) tol=1e-7

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    W = st_data(., tokens(w_var), touse)
    P = st_data(., tokens(p_var), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        w0 = W[k, ]'
        p0 = P[k, ]'
        profit_eff = .
        
        // Observed profit
        obs_profit = sum(p0 :* y0) - sum(w0 :* x0)
        
        // Variables: λ(n) only
        // Objective: max p0'(Y'λ) - w0'(X'λ) = (Yp0 - Xw0)'λ
        //         => min (Xw0 - Yp0)'λ
        c = (X * w0 - Y * p0)'
        
        // Constraints:
        //   X'λ ≤ x0
        //   -Y'λ ≤ -y0
        A_ineq = X' \ (-Y')
        b_ineq = x0 \ (-y0)
        
        lb = J(1, n, 0)
        ub = J(1, n, .)
        
        A_eq = J(0, n, .)
        b_eq = J(0, 1, .)
        
        if (rts == "vrs") {
            A_eq = J(1, n, 1)
            b_eq = 1
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ J(1, n, 1)
            b_ineq = b_ineq \ 1
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ J(1, n, -1)
            b_ineq = b_ineq \ (-1)
        }
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setInequality(A_ineq, b_ineq)
        if (rows(A_eq) > 0) {
            S.setEquality(A_eq, b_eq)
        }
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = _dea_nonte_optimize(S,tol)
        
        if (rc == 0) {
            max_profit = -S.value()
            if (max_profit > 0) {
                profit_eff = obs_profit / max_profit
            }
        }
        
        st_store(obs_k, st_varindex(eff_var), profit_eff)
    }
}


// ============================================================================
// 11. SBM (SLACKS-BASED MEASURE) MODEL
// ============================================================================
/**
 * _dea_ipm_sbm()
 *
 * Computes Slacks-Based Measure efficiency (Tone, 2001).
 * Uses Charnes-Cooper transformation to convert fractional program to LP.
 *
 * Original fractional program:
 *   min rho = [1 - (1/m)*sum(s-_i/x_io)] / [1 + (1/s)*sum(s+_r/y_ro)]
 *   s.t. x_o = X*lambda + s-
 *        y_o = Y*lambda - s+
 *        [convexity]
 *        lambda, s-, s+ >= 0
 *
 * After Charnes-Cooper transformation (t = 1/denominator):
 *   min t - (1/m)*sum(S-_i/x_io)
 *   s.t. t + (1/s)*sum(S+_r/y_ro) = 1
 *        t*x_o = X*Lambda + S-
 *        t*y_o = Y*Lambda - S+
 *        [convexity: sum(Lambda) = t]
 *        Lambda, S-, S+, t >= 0
 */
void _dea_ipm_sbm(string scalar inputs,
                  string scalar outputs,
                  string scalar eff_var,
                  string scalar slack_in_var,
                  string scalar slack_out_var,
                  string scalar touse,
                  string scalar rts,
    | real scalar tol)
{
    real matrix X, Y
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, sol
    real matrix A_eq
    real colvector b_eq
    real matrix A_ineq
    real colvector b_ineq
    real scalar n, m, s, k, obs_k, rc
    real scalar sbm_eff, t_val
    real rowvector slack_in, slack_out, coef_sin, coef_sout
    class LinearProgram scalar S
    real scalar j
    
    if (args()<8) tol=1e-7

    real rowvector x_scale, y_scale
    real scalar mp, sp

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    // Column rescaling preserves the SBM ratios and balances the LP.
    x_scale = colmax(X); y_scale = colmax(Y)
    for (j=1; j<=cols(X); j++) if (x_scale[j]==0) x_scale[j]=1
    for (j=1; j<=cols(Y); j++) if (y_scale[j]==0) y_scale[j]=1
    X = X :/ x_scale; Y = Y :/ y_scale
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    string rowvector si_varlist, so_varlist
    real scalar has_slacks
    
    has_slacks = (slack_in_var != "" & slack_out_var != "")
    if (has_slacks) {
        si_varlist = tokens(slack_in_var)
        so_varlist = tokens(slack_out_var)
    }
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        sbm_eff = .
        slack_in = J(1, m, 0)
        slack_out = J(1, s, 0)
        
        // Variables: [Lambda(n), S-(m), S+(s), t]
        // Total: n + m + s + 1 variables
        
        // Coefficients for input slacks: -(1/m)/x_io
        // Tone (2001), section 6: omitted zero terms also reduce the mean.
        // Zero outputs use the structural-zero interpretation (no perturbation).
        mp = sum(x0:>0); sp = sum(y0:>0)
        if (mp==0 | sp==0) {
            st_store(obs_k, st_varindex(eff_var), .)
            continue
        }
        coef_sin = J(1, m, 0)
        for (j = 1; j <= m; j++) {
            if (x0[j] > 0) {
                coef_sin[j] = -1 / (mp * x0[j])
            }
        }
        
        // Coefficients for output slacks in normalizing constraint: (1/s)/y_ro
        coef_sout = J(1, s, 0)
        for (j = 1; j <= s; j++) {
            if (y0[j] > 0) {
                coef_sout[j] = 1 / (sp * y0[j])
            }
        }
        
        // Objective: min t + sum(coef_sin * S-)
        c = J(1, n, 0), coef_sin, J(1, s, 0), 1
        
        // Equality constraints:
        // (1) Normalizing: t + sum(coef_sout * S+) = 1
        // (2) Input:  X'*Lambda + I*S- - t*x0 = 0
        // (3) Output: Y'*Lambda - I*S+ - t*y0 = 0
        
        // (1) Normalizing constraint
        A_eq = (J(1, n, 0), J(1, m, 0), coef_sout, 1)
        b_eq = 1
        
        // (2) Input constraints: X'*Lambda + I*S- = t*x0
        //     => X'*Lambda + I*S- - x0*t = 0
        A_eq = A_eq \ (X', I(m), J(m, s, 0), -x0)
        b_eq = b_eq \ J(m, 1, 0)
        
        // (3) Output constraints: Y'*Lambda - I*S+ = t*y0
        //     => Y'*Lambda - I*S+ - y0*t = 0
        A_eq = A_eq \ (Y', J(s, m, 0), -I(s), -y0)
        b_eq = b_eq \ J(s, 1, 0)
        
        // RTS constraint
        A_ineq = J(0, n+m+s+1, .)
        b_ineq = J(0, 1, .)
        
        if (rts == "vrs") {
            // sum(Lambda) = t => sum(Lambda) - t = 0
            A_eq = A_eq \ (J(1, n, 1), J(1, m+s, 0), -1)
            b_eq = b_eq \ 0
        }
        else if (rts == "drs") {
            // sum(Lambda) <= t => sum(Lambda) - t <= 0
            A_ineq = (J(1, n, 1), J(1, m+s, 0), -1)
            b_ineq = 0
        }
        else if (rts == "irs") {
            // sum(Lambda) >= t => -sum(Lambda) + t <= 0
            A_ineq = (J(1, n, -1), J(1, m+s, 0), 1)
            b_ineq = 0
        }
        // CRS: sum(Lambda) is free (no constraint added)
        
        lb = J(1, n+m+s+1, 0)
        ub = J(1, n+m+s+1, .)
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setEquality(A_eq, b_eq)
        if (rows(A_ineq) > 0) {
            S.setInequality(A_ineq, b_ineq)
        }
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = _dea_nonte_optimize(S,tol)
        
        if (rc == 0) {
            sol = S.parameters()
            t_val = sol[n+m+s+1]
            if (t_val > 0) {
                // Convert back: s- = S-/t, s+ = S+/t
                slack_in = sol[n+1..n+m] / t_val
                slack_out = sol[n+m+1..n+m+s] / t_val
                sbm_eff = S.value()
            }
        }
        
        // Public/internal slack requests remain in the original data units.
        slack_in = slack_in :* x_scale
        slack_out = slack_out :* y_scale
        st_store(obs_k, st_varindex(eff_var), sbm_eff)
        
        if (has_slacks) {
            for (j = 1; j <= m; j++) {
                st_store(obs_k, st_varindex(si_varlist[j]), slack_in[j])
            }
            for (j = 1; j <= s; j++) {
                st_store(obs_k, st_varindex(so_varlist[j]), slack_out[j])
            }
        }
    }
}


// ============================================================================
// 12. MULTIPLIER (DUAL) MODEL
// ============================================================================
/**
 * _dea_ipm_multiplier()
 *
 * Solves the multiplier (dual) form of DEA.
 *
 * Input-oriented (max virtual output ratio):
 *   max u'*y0 - u0  (for VRS)
 *   s.t. v'*x0 = 1
 *        u'*Y - v'*X - u0 <= 0  (for each DMU)
 *        u, v >= epsilon
 *
 * Output-oriented (min virtual input ratio):
 *   min v'*x0 + v0  (for VRS)
 *   s.t. u'*y0 = 1
 *        u'*Y - v'*X - v0 >= 0
 *        u, v >= epsilon
 */
void _dea_ipm_multiplier(string scalar inputs,
                          string scalar outputs,
                          string scalar eff_var,
                          string scalar v_var,
                          string scalar u_var,
                          string scalar u0_var,
                          string scalar touse,
                          string scalar orientation,
                          string scalar rts,
    | real scalar tol)
{
    real matrix X, Y
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, sol
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc, eff_score, u0_value
    real scalar epsilon, has_u0
    real rowvector v_weights, u_weights
    class LinearProgram scalar S
    
    if (args()<10) tol=1e-7

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    epsilon = 0     // Standard multiplier form permits zero weights
    has_u0 = (rts == "vrs" | rts == "drs" | rts == "irs")
    
    string rowvector v_varlist, u_varlist
    real scalar has_weights, j
    
    has_weights = (v_var != "" & u_var != "")
    if (has_weights) {
        v_varlist = tokens(v_var)
        u_varlist = tokens(u_var)
    }
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        eff_score = .
        v_weights = J(1, m, .)
        u_weights = J(1, s, .)
        u0_value = .
        
        if (orientation == "input") {
            // Variables: [v(m), u(s), u0] if VRS, else [v(m), u(s)]
            // Objective: max u'*y0 - u0 => min -u'*y0 + u0
            
            if (has_u0) {
                c = J(1, m, 0), -y0', 1
                
                // Normalizing: v'*x0 = 1
                A_eq = (x0', J(1, s, 0), 0)
                b_eq = 1
                
                // Constraints: u'*y_j - v'*x_j - u0 <= 0 for all j
                A_ineq = (-X, Y, J(n, 1, -1))
                b_ineq = J(n, 1, 0)
                
                lb = J(1, m, epsilon), J(1, s, epsilon), .
                ub = J(1, m+s+1, .)
                
                // Cooper signed intercept: DRS u0 >= 0; IRS u0 <= 0
                if (rts == "drs") {
                    lb[m+s+1] = 0
                }
                else if (rts == "irs") {
                    ub[m+s+1] = 0
                }
            }
            else {
                // CRS
                c = J(1, m, 0), -y0'
                
                A_eq = (x0', J(1, s, 0))
                b_eq = 1
                
                A_ineq = (-X, Y)
                b_ineq = J(n, 1, 0)
                
                lb = J(1, m, epsilon), J(1, s, epsilon)
                ub = J(1, m+s, .)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setEquality(A_eq, b_eq)
            S.setInequality(A_ineq, b_ineq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = _dea_nonte_optimize(S,tol)
            
            if (rc == 0) {
                sol = S.parameters()
                eff_score = -S.value()
                v_weights = sol[1..m]
                u_weights = sol[m+1..m+s]
                u0_value = (has_u0 ? sol[m+s+1] : 0)
            }
        }
        else {
            // Output-oriented
            // Variables: [v(m), u(s), v0] if VRS
            // Objective: min v'*x0 + v0
            
            if (has_u0) {
                c = x0', J(1, s, 0), 1
                
                // Normalizing: u'*y0 = 1
                A_eq = (J(1, m, 0), y0', 0)
                b_eq = 1
                
                // Constraints: v'*x_j - u'*y_j + v0 >= 0 => -v'*x_j + u'*y_j - v0 <= 0
                A_ineq = (-X, Y, J(n, 1, -1))
                b_ineq = J(n, 1, 0)
                
                lb = J(1, m, epsilon), J(1, s, epsilon), .
                ub = J(1, m+s+1, .)
                
                if (rts == "drs") {
                    lb[m+s+1] = 0
                }
                else if (rts == "irs") {
                    ub[m+s+1] = 0
                }
            }
            else {
                c = x0', J(1, s, 0)
                
                A_eq = (J(1, m, 0), y0')
                b_eq = 1
                
                A_ineq = (-X, Y)
                b_ineq = J(n, 1, 0)
                
                lb = J(1, m, epsilon), J(1, s, epsilon)
                ub = J(1, m+s, .)
            }
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setEquality(A_eq, b_eq)
            S.setInequality(A_ineq, b_ineq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = _dea_nonte_optimize(S,tol)
            
            if (rc == 0) {
                sol = S.parameters()
                eff_score = 1 / S.value()  // Convert to standard efficiency
                v_weights = sol[1..m]
                u_weights = sol[m+1..m+s]
                u0_value = (has_u0 ? sol[m+s+1] : 0)
            }
        }
        
        st_store(obs_k, st_varindex(eff_var), eff_score)
        
        if (has_weights) {
            for (j = 1; j <= m; j++) {
                st_store(obs_k, st_varindex(v_varlist[j]), v_weights[j])
            }
            for (j = 1; j <= s; j++) {
                st_store(obs_k, st_varindex(u_varlist[j]), u_weights[j])
            }
            if (u0_var != "") {
                st_store(obs_k, st_varindex(u0_var), u0_value)
            }
        }
    }
}


// ============================================================================
// 13. SUPER-SBM MODEL (Tone, 2002)
// ============================================================================
/**
 * _dea_ipm_supersbm()
 *
 * Tone (2002): min mean(x_bar/x_o) / mean(y_bar/y_o), with
 * x_bar >= X_{-o}'lambda, x_bar >= x_o,
 * y_bar <= Y_{-o}'lambda, 0 <= y_bar <= y_o, lambda >= 0.
 * Ordinary SBM scores below 1-tol are retained in the chosen orientation.
 * Input orientation fixes y_bar=y_o; output fixes x_bar=x_o and returns
 * the reciprocal maximum mean(y_bar/y_o). Scores equal to 1 are possible.
 */
void _dea_ipm_supersbm(string scalar inputs,
                        string scalar outputs,
                        string scalar eff_var,
                        string scalar touse,
                        string scalar rts,
    | real scalar tol, string scalar orientation)
{
    real matrix X, Y, X_ref, Y_ref
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, coef_xbar, coef_ybar
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, obs_k, rc, n_ref, j
    real scalar sbm_eff, sbm_score
    class LinearProgram scalar S
    
    real colvector ordinary_scores
    real rowvector x_scale, y_scale
    string scalar ort

    if (args()<6) tol=1e-7

    if (args()<7) orientation=""
    ort = strlower(orientation)
    if (ort=="in") ort="input"
    if (ort=="out") ort="output"
    if (ort=="inout") ort=""
    // Classify and retain ordinary SBM scores using the requested orientation.
    if (ort=="") _dea_ipm_sbm(inputs, outputs, eff_var, "", "", touse, rts, tol)
    else _dea_ipm_sbm_oriented(inputs, outputs, eff_var, touse, rts, ort, tol)
    ordinary_scores = st_data(., eff_var, touse)

    // No peer exists even if rounding puts the ordinary score just below 1.
    if (rows(ordinary_scores)<2) {
        st_store(., eff_var, touse, J(rows(ordinary_scores),1,.))
        return
    }

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    x_scale = colmax(X); y_scale = colmax(Y)
    X = X :/ x_scale; Y = Y :/ y_scale
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k, ]'
        y0 = Y[k, ]'
        sbm_eff = .
        
        sbm_score = ordinary_scores[k]

        // If SBM-inefficient, return SBM score directly
        if (sbm_score == . || sbm_score < 1 - tol) {
            st_store(obs_k, st_varindex(eff_var), sbm_score)
            continue
        }
        
        // ==============================================================
        // Step 2: Super-SBM for efficient DMUs — Tone (2002) eq. 14
        //   Charnes-Cooper linearization:
        //   Variables: [Lambda(n_ref), x_tilde(m), y_tilde(s), t]
        //
        //   min  (1/m) * sum(x_tilde_i / x_i0)
        //   s.t. (1/s) * sum(y_tilde_r / y_r0) = 1
        //        x_tilde >= X_ref' * Lambda
        //        y_tilde <= Y_ref' * Lambda
        //        x_tilde >= t * x_0
        //        y_tilde <= t * y_0
        //        VRS: sum(Lambda) = t
        //        Lambda >= 0, y_tilde >= 0, t > 0
        // ==============================================================
        
        n_ref = n - 1
        if (n_ref < 1) {
            st_store(obs_k, st_varindex(eff_var), .)
            continue
        }
        
        if (k == 1) {
            X_ref = X[2..n, ]
            Y_ref = Y[2..n, ]
        }
        else if (k == n) {
            X_ref = X[1..n-1, ]
            Y_ref = Y[1..n-1, ]
        }
        else {
            X_ref = X[1..k-1, ] \ X[k+1..n, ]
            Y_ref = Y[1..k-1, ] \ Y[k+1..n, ]
        }
        
        coef_xbar = J(1, m, 0)
        for (j = 1; j <= m; j++) {
            if (x0[j] > 0) coef_xbar[j] = 1 / (m * x0[j])
        }
        coef_ybar = J(1, s, 0)
        for (j = 1; j <= s; j++) {
            if (y0[j] > 0) coef_ybar[j] = 1 / (s * y0[j])
        }
        
        // Non-oriented: fractional model. Oriented: t=1 and fix the
        // opposite-side projection. Output score is the reciprocal maximum.
        c = J(1, n_ref, 0), coef_xbar, J(1, s, 0), 0
        if (ort=="") {
            A_eq = (J(1, n_ref+m, 0), coef_ybar, 0)
            b_eq = 1
        }
        else {
            A_eq = (J(1, n_ref+m+s, 0), 1)
            b_eq = 1
            if (ort=="input") {
                A_eq = A_eq \ (J(s,n_ref+m,0), I(s), -y0)
                b_eq = b_eq \ J(s,1,0)
            }
            else {
                c = J(1,n_ref+m,0), -coef_ybar, 0
                A_eq = A_eq \ (J(m,n_ref,0), I(m), J(m,s,0), -x0)
                b_eq = b_eq \ J(m,1,0)
            }
        }

        // Inequality 1: x_tilde >= X_ref' * Lambda
        A_ineq = (X_ref', -I(m), J(m, s, 0), J(m, 1, 0))
        b_ineq = J(m, 1, 0)
        
        // Inequality 2: y_tilde <= Y_ref' * Lambda
        A_ineq = A_ineq \ (-Y_ref', J(s, m, 0), I(s), J(s, 1, 0))
        b_ineq = b_ineq \ J(s, 1, 0)
        
        // Inequality 3: x_tilde >= t * x_0
        A_ineq = A_ineq \ (J(m, n_ref, 0), -I(m), J(m, s, 0), x0)
        b_ineq = b_ineq \ J(m, 1, 0)
        
        // Inequality 4: y_tilde <= t * y_0
        A_ineq = A_ineq \ (J(s, n_ref, 0), J(s, m, 0), I(s), -y0)
        b_ineq = b_ineq \ J(s, 1, 0)
        
        // RTS constraints
        if (rts == "vrs") {
            A_eq = A_eq \ (J(1, n_ref, 1), J(1, m+s, 0), -1)
            b_eq = b_eq \ 0
        }
        else if (rts == "drs") {
            A_ineq = A_ineq \ (J(1, n_ref, 1), J(1, m+s, 0), -1)
            b_ineq = b_ineq \ 0
        }
        else if (rts == "irs") {
            A_ineq = A_ineq \ (J(1, n_ref, -1), J(1, m+s, 0), 1)
            b_ineq = b_ineq \ 0
        }
        
        lb = J(1, n_ref + m + s + 1, 0)
        ub = J(1, n_ref + m + s + 1, .)
        
        S = LinearProgram()
        S.setMaxOrMin("min")
        S.setCoefficients(c)
        S.setEquality(A_eq, b_eq)
        S.setInequality(A_ineq, b_ineq)
        S.setBounds(lb, ub)
        S.setTrace("off")
        
        rc = _dea_nonte_optimize(S,tol)
        
        if (rc == 0) {
            if (ort!="output") sbm_eff = S.value()
            else if (S.value()<0) sbm_eff = -1/S.value()
        }
        
        st_store(obs_k, st_varindex(eff_var), sbm_eff)
    }
}


// ============================================================================
// 14. FDH (FREE DISPOSAL HULL) - No LP Required
// ============================================================================
/**
 * _dea_fdh()
 *
 * Free Disposal Hull efficiency using dominance comparison.
 * No convexity assumption - only actual observed DMUs as references.
 *
 * Input-oriented:
 *   theta_k = min_j max_i (x_ji / x_ki), over references with y_j >= y_k.
 *   This is the radial Farrell input score relative to the non-convex FDH.
 *
 * Output-oriented:
 *   phi_k = max_j min_r (y_jr / y_kr), over references with x_j <= x_k.
 *   Efficiency is stored as 1/phi, matching dea.ado output-score convention.
 */
void _dea_fdh(string scalar inputs,
              string scalar outputs,
              string scalar eff_var,
              string scalar touse,
              string scalar orientation)
{
    real matrix X, Y
    real colvector obs_indices
    real scalar n, m, s, k, j, i, r, obs_k
    real scalar theta, phi, ratio, min_ratio, max_ratio
    real scalar feasible_ref
    real colvector x_k, y_k, x_j, y_j
    
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        x_k = X[k, ]'
        y_k = Y[k, ]'
        
        if (orientation == "input") {
            // Input-oriented radial FDH:
            // among output-feasible references, minimize the maximum input ratio.
            theta = 1  // default: efficient
            
            for (j = 1; j <= n; j++) {
                x_j = X[j, ]'
                y_j = Y[j, ]'
                
                feasible_ref = 1
                for (r = 1; r <= s; r++) {
                    if (y_j[r] < y_k[r]) {
                        feasible_ref = 0
                        break
                    }
                }
                
                if (feasible_ref) {
                    max_ratio = 0
                    for (i = 1; i <= m; i++) {
                        if (x_k[i] > 0) {
                            ratio = x_j[i] / x_k[i]
                            if (ratio > max_ratio) {
                                max_ratio = ratio
                            }
                        }
                    }
                    if (max_ratio < theta) {
                        theta = max_ratio
                    }
                }
            }
            
            st_store(obs_k, st_varindex(eff_var), theta)
        }
        else {
            // Output-oriented radial FDH:
            // among input-feasible references, maximize the minimum output ratio.
            phi = 1  // default: efficient (phi=1 means eff=1)
            
            for (j = 1; j <= n; j++) {
                x_j = X[j, ]'
                y_j = Y[j, ]'
                
                feasible_ref = 1
                for (i = 1; i <= m; i++) {
                    if (x_j[i] > x_k[i]) {
                        feasible_ref = 0
                        break
                    }
                }
                
                if (feasible_ref) {
                    min_ratio = .
                    for (r = 1; r <= s; r++) {
                        if (y_k[r] > 0) {
                            ratio = y_j[r] / y_k[r]
                            if (min_ratio == . || ratio < min_ratio) {
                                min_ratio = ratio
                            }
                        }
                    }
                    if (min_ratio != . && min_ratio > phi) {
                        phi = min_ratio
                    }
                }
            }
            
            // Output efficiency = 1/phi (or just store phi for expansion factor)
            st_store(obs_k, st_varindex(eff_var), 1/phi)
        }
    }
}


end  // end of mata block


// ============================================================================
// 15. IDEA (IMPRECISE DEA) - Interval/Bounded Data
// ============================================================================
/**
 * IDEA computes efficiency under data uncertainty.
 * 
 * For interval data [xlo, xhi] and [ylo, yhi]:
 *   - Pessimistic: DMU k at worst (high inputs, low outputs)
 *                  Other DMUs at best (low inputs, high outputs)
 *   - Optimistic:  DMU k at best (low inputs, high outputs)
 *                  Other DMUs at worst (high inputs, low outputs)
 * In each LP, the reference row for k equals its evaluated realization.
 * These are worst/best relative efficiencies over independent data intervals,
 * not scores against one fixed frontier. Output scores are reciprocals of phi.
 *
 * Uses IPM solver for LP.
 */

mata:

void _dea_idea(string scalar inputs,
               string scalar outputs,
               string scalar xlo_vars,
               string scalar xhi_vars,
               string scalar ylo_vars,
               string scalar yhi_vars,
               string scalar eff_pess_var,
               string scalar eff_opt_var,
               string scalar touse,
               string scalar orientation,
               string scalar rts,
               real scalar tol)
{
    real matrix X, Y, Xlo, Xhi, Ylo, Yhi
    real matrix X_pess, Y_pess, X_opt, Y_opt
    real colvector obs_indices, x0, y0
    real rowvector c, lb, ub, sol, lam_k
    real matrix A_ineq, A_eq
    real colvector b_ineq, b_eq
    real scalar n, m, s, k, j, obs_k, rc
    real scalar eff_pess, eff_opt
    class LinearProgram scalar S
    
    // Load base data
    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    obs_indices = selectindex(st_data(., touse))
    
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    
    // Load interval bounds (if provided)
    if (xlo_vars != "") {
        Xlo = st_data(., tokens(xlo_vars), touse)
    } else {
        Xlo = X
    }
    if (xhi_vars != "") {
        Xhi = st_data(., tokens(xhi_vars), touse)
    } else {
        Xhi = X
    }
    if (ylo_vars != "") {
        Ylo = st_data(., tokens(ylo_vars), touse)
    } else {
        Ylo = Y
    }
    if (yhi_vars != "") {
        Yhi = st_data(., tokens(yhi_vars), touse)
    } else {
        Yhi = Y
    }
    
    // Tight convergence keeps near-unit endpoint ordering stable; the requested
    // tolerance still governs original-LP residual acceptance. No clipping.
    for (k = 1; k <= n; k++) {
        obs_k = obs_indices[k]
        eff_pess = .
        eff_opt = .
        
        // ================================================================
        // PESSIMISTIC EFFICIENCY
        // DMU k: worst case (Xhi, Ylo)
        // Other DMUs: best case (Xlo, Yhi); own row stays at (Xhi[k], Ylo[k]).
        // ================================================================
        X_pess = Xlo
        Y_pess = Yhi
        X_pess[k, ] = Xhi[k, ]
        Y_pess[k, ] = Ylo[k, ]
        if (orientation == "input") {
            // min theta over the complete worst-case realized dataset.
            x0 = Xhi[k, ]'
            y0 = Ylo[k, ]'
            
            c = J(1, n, 0), 1
            A_ineq = (X_pess', -x0) \ (-Y_pess', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)
            
            lb = J(1, n+1, 0)
            ub = J(1, n+1, .)
            
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
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) S.setEquality(A_eq, b_eq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = _dea_nonte_optimize(S,tol,1e-12)
            if (rc == 0) {
                eff_pess = S.value()
            }
        }
        else {
            // Output-oriented pessimistic
            x0 = Xhi[k, ]'
            y0 = Ylo[k, ]'
            
            c = J(1, n, 0), -1
            A_ineq = (X_pess', J(m, 1, 0)) \ (-Y_pess', y0)
            b_ineq = x0 \ J(s, 1, 0)
            
            lb = J(1, n, 0), 1
            ub = J(1, n+1, .)
            
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
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) S.setEquality(A_eq, b_eq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = _dea_nonte_optimize(S,tol,1e-12)
            if (rc == 0) {
                eff_pess = 1 / (-S.value())
            }
        }
        
        // ================================================================
        // OPTIMISTIC EFFICIENCY
        // DMU k: best case (Xlo, Yhi)
        // Other DMUs: worst case (Xhi, Ylo); own row stays at (Xlo[k], Yhi[k]).
        // ================================================================
        X_opt = Xhi
        Y_opt = Ylo
        X_opt[k, ] = Xlo[k, ]
        Y_opt[k, ] = Yhi[k, ]
        if (orientation == "input") {
            x0 = Xlo[k, ]'
            y0 = Yhi[k, ]'
            
            c = J(1, n, 0), 1
            A_ineq = (X_opt', -x0) \ (-Y_opt', J(s, 1, 0))
            b_ineq = J(m, 1, 0) \ (-y0)
            
            lb = J(1, n+1, 0)
            ub = J(1, n+1, .)
            
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
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) S.setEquality(A_eq, b_eq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = _dea_nonte_optimize(S,tol,1e-12)
            if (rc == 0) {
                eff_opt = S.value()
            }
        }
        else {
            // Output-oriented optimistic
            x0 = Xlo[k, ]'
            y0 = Yhi[k, ]'
            
            c = J(1, n, 0), -1
            A_ineq = (X_opt', J(m, 1, 0)) \ (-Y_opt', y0)
            b_ineq = x0 \ J(s, 1, 0)
            
            lb = J(1, n, 0), 1
            ub = J(1, n+1, .)
            
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
            
            S = LinearProgram()
            S.setMaxOrMin("min")
            S.setCoefficients(c)
            S.setInequality(A_ineq, b_ineq)
            if (rows(A_eq) > 0) S.setEquality(A_eq, b_eq)
            S.setBounds(lb, ub)
            S.setTrace("off")
            
            rc = _dea_nonte_optimize(S,tol,1e-12)
            if (rc == 0) {
                eff_opt = 1 / (-S.value())
            }
        }
        
        // Store results
        st_store(obs_k, st_varindex(eff_pess_var), eff_pess)
        st_store(obs_k, st_varindex(eff_opt_var), eff_opt)
    }
}

// Oriented SBM (Cooper, Seiford, and Tone 2007, eqs. 4.78--4.79).
// Input orientation maximizes normalized input slack and reports 1-z.
// Output orientation maximizes normalized output slack and reports 1/(1+z).
void _dea_ipm_sbm_oriented(string scalar inputs,
                           string scalar outputs,
                           string scalar eff_var,
                           string scalar touse,
                           string scalar rts,
                           string scalar orientation,
    | real scalar tol)
{
    real matrix X, Y, A_eq, A_ineq
    real colvector obs_indices, x0, y0, b_eq, b_ineq
    real rowvector c, lb, ub, sol, coef_sin, coef_sout
    real scalar n, m, s, k, j, obs_k, nv, obj, eff, rc
    class LinearProgram scalar L

    if (args()<7) tol=1e-7

    real rowvector x_scale, y_scale
    real scalar mp, sp

    X = st_data(., tokens(inputs), touse)
    Y = st_data(., tokens(outputs), touse)
    // Column rescaling preserves the SBM ratios and balances the LP.
    x_scale = colmax(X); y_scale = colmax(Y)
    for (j=1; j<=cols(X); j++) if (x_scale[j]==0) x_scale[j]=1
    for (j=1; j<=cols(Y); j++) if (y_scale[j]==0) y_scale[j]=1
    X = X :/ x_scale; Y = Y :/ y_scale
    obs_indices = selectindex(st_data(., touse))
    n = rows(X); m = cols(X); s = cols(Y); nv = n + m + s

    for (k=1; k<=n; k++) {
        obs_k = obs_indices[k]
        x0 = X[k,.]'; y0 = Y[k,.]'
        // Tone (2001), section 6: omitted zero terms also reduce the mean.
        // Zero outputs use the structural-zero interpretation (no perturbation).
        mp = sum(x0:>0); sp = sum(y0:>0)
        if (mp==0 | sp==0) {
            st_store(obs_k, st_varindex(eff_var), .)
            continue
        }
        coef_sin = J(1,m,0); coef_sout = J(1,s,0)
        for (j=1; j<=m; j++) if (x0[j]>0) coef_sin[j] = 1/(mp*x0[j])
        for (j=1; j<=s; j++) if (y0[j]>0) coef_sout[j] = 1/(sp*y0[j])

        if (orientation=="input") c = J(1,n,0), -coef_sin, J(1,s,0)
        else c = J(1,n,0), J(1,m,0), -coef_sout

        A_eq = (X', I(m), J(m,s,0))
        b_eq = x0
        A_eq = A_eq \ (Y', J(s,m,0), -I(s))
        b_eq = b_eq \ y0
        A_ineq = J(0,nv,0); b_ineq = J(0,1,0)

        if (rts=="vrs") {
            A_eq = A_eq \ (J(1,n,1), J(1,m+s,0))
            b_eq = b_eq \ 1
        }
        else if (rts=="drs") {
            A_ineq = (J(1,n,1), J(1,m+s,0)); b_ineq = 1
        }
        else if (rts=="irs") {
            A_ineq = (J(1,n,-1), J(1,m+s,0)); b_ineq = -1
        }

        lb = J(1,nv,0); ub = J(1,nv,.)
        L = LinearProgram()
        L.setMaxOrMin("min")
        L.setCoefficients(c)
        L.setEquality(A_eq,b_eq)
        if (rows(A_ineq)>0) L.setInequality(A_ineq,b_ineq)
        L.setBounds(lb,ub); L.setTrace("off")
        rc = _dea_nonte_optimize(L,tol)
        eff = .
        if (rc == 0) {
            obj = L.value()
            if (orientation=="input") eff = 1 + obj
            else if (1-obj>0) eff = 1/(1-obj)
        }
        st_store(obs_k,st_varindex(eff_var),eff)
    }
}

end
