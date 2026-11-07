*! version 1.3.1 08nov2026
*! Cross-efficiency simplex, original-data optimal faces and recovery diagnostics.
version 16.0
capture mata: mata drop _c57_simplex_ready()
capture mata: mata drop _c57_solve()
capture mata: mata drop _c57_multiplier()
capture mata: mata drop _c57_compute()
capture mata: mata drop _c58_multiplier()
capture mata: mata drop _c57_dmu_read()
capture mata: mata drop _c57_lzratio()
capture mata: mata drop _c57_rank()
mata:
mata set matastrict on
real scalar _c57_simplex_ready() return(8402)
// Split each free variable into two nonnegative columns, then check the
// reconstructed solution against the submitted, unsplit constraints.
// Diagnostics: accepted-status, solver-status, objective, iterations,
// maximum scaled primal residual, maximum absolute primal residual.
real scalar _c57_solve(real rowvector c,real matrix E,real colvector b,
    real matrix U,real colvector h,real rowvector lo,real rowvector hi,
    real scalar tol,real rowvector x,real scalar obj,real rowvector diag)
{
 real matrix T
 real rowvector z,quality
 real scalar n,j
 struct RSM_Result scalar R
 n=cols(c);T=I(n);x=J(1,n,.);obj=.;diag=(91,.,.,.,.,.)
 if(hasmissing(c) | tol<=0 | tol>=. | cols(lo)!=n | cols(hi)!=n) return(91)
 for(j=1;j<=n;j++) {
  if(hi[j]<. | (lo[j]<. & lo[j]!=0)) return(91)
  if(lo[j]>=.) T=(T,-I(n)[.,j])
 }
 R=_rsm_solve_lp(U*T,h,E*T,b,c*T,tol)
 diag[2]=R.status;diag[4]=R.iters
 if(R.status!=0 | rows(R.x)!=cols(T) | hasmissing(R.x)) return(91)
 x=(T*R.x)';obj=c*x'
 quality=_dea_lpcheck(x,E,b,U,h,lo,hi)
 diag[3]=obj;diag[5]=max(quality[4..6]);diag[6]=max(quality[1..3])
 if(obj>=. | !_dea_lpcheck_ok(quality,tol)) return(91)
 diag[1]=0
 return(0)
}
void _c57_dmu_read(string scalar  dmu_var,
                   string scalar  touse,
                   real colvector dmu_ids,
                   string colvector dmu_labels)
{
    real scalar n, i, j
    string colvector dmu_str, uniq

    if (substr(st_vartype(dmu_var), 1, 3) == "str") {
        dmu_str    = st_sdata(., dmu_var, touse)
        n          = rows(dmu_str)
        uniq       = sort(uniqrows(dmu_str), 1)
        dmu_ids    = J(n, 1, .)
        dmu_labels = dmu_str
        for (i = 1; i <= n; i++) {
            for (j = 1; j <= rows(uniq); j++) {
                if (dmu_str[i] == uniq[j]) {
                    dmu_ids[i] = j
                    break
                }
            }
        }
    }
    else {
        dmu_ids    = st_data(., dmu_var, touse)
        dmu_labels = strtrim(strofreal(dmu_ids,"%21.17g"))
    }
}
real scalar _c57_lzratio(real scalar vx, real scalar uy,
                                 real scalar w, real scalar tol)
{
    real scalar den, scale, score
    if (vx>=. | uy>=. | w>=. | vx<0 | uy<0) return(.)
    den=vx-w
    scale=max((abs(vx),abs(w),uy))
    if (den>=. | den<=0 | den<=tol*scale) return(.)
    score=uy/den
    if (score>=. | score<0 | score>1+tol) return(.)
    return(score)
}
real colvector _c57_rank(real colvector x)
{
    real scalar n, i, j, rank
    real colvector ranks

    n = rows(x)
    ranks = J(n, 1, .)

    for (i = 1; i <= n; i++) {
        if (x[i] == .) {
            ranks[i] = .
            continue
        }

        rank = 1
        for (j = 1; j <= n; j++) {
            if (j != i && x[j] != . && x[j] > x[i]) {
                rank++
            }
        }
        ranks[i] = rank
    }

    return(ranks)
}
void _c57_multiplier(real rowvector x_original,
                               real rowvector y_original,
                               real matrix X_original,
                               real matrix Y_original,
                               string scalar rts,
                               string scalar orientation,
                               string scalar method,
                               real scalar tol,
                               real rowvector v_opt,
                               real rowvector u_opt,
                               real scalar w_opt,
                               real scalar eff, real scalar stage_status, real rowvector diag)
{
    real scalar n, m, s, j
    real matrix X, Y
    real rowvector x0, y0, sx, sy
    real colvector row_scale
    // Reparameterize weights around the evaluated row, then balance each
    // peer constraint. This preserves the LP and keeps normalization usable
    // when DMUs have very different magnitudes within the same column.
    sx=x_original; sy=y_original
    for (j=1; j<=cols(sx); j++) if (sx[j]<=0) sx[j]=1
    for (j=1; j<=cols(sy); j++) if (sy[j]<=0) sy[j]=1
    X=X_original:/sx; Y=Y_original:/sy
    x0=x_original:/sx; y0=y_original:/sy
    n = rows(X)
    m = cols(X)
    s = cols(Y)

    real scalar is_vrs, is_output, nvars, rc
    real rowvector c, lb, ub
    real matrix    A_ineq, A_eq
    real colvector b_ineq, b_eq
    real rowvector sol, d1, d2
    real scalar obj
    diag=J(1,12,.)

    stage_status=0
    is_vrs    = (rts == "vrs")
    is_output = (orientation == "output")
    nvars     = m + s + is_vrs        // add 1 free variable w if VRS

    // ── Objective (Stage 1: self-efficiency) ─────────────────────────
    if (!is_output) {
        // Input: min −u·y₀ [− w]
        c = (J(1, m, 0), -y0)
        if (is_vrs) c = (c, -1)
    }
    else {
        // Output: min v·x₀ [+ w]
        c = (x0, J(1, s, 0))
        if (is_vrs) c = (c, 1)
    }

    // ── Equality constraint ──────────────────────────────────────────
    if (!is_output) {
        // Input: v·x₀ = 1
        A_eq = !is_vrs ? (x0, J(1, s, 0)) : (x0, J(1, s, 0), 0)
    }
    else {
        // Output: u·y₀ = 1
        A_eq = !is_vrs ? (J(1, m, 0), y0) : (J(1, m, 0), y0, 0)
    }
    b_eq = 1

    // ── Inequality constraints: [−X, Y {, ±1}]·z ≤ 0 ────────────────
    if (!is_vrs) {
        A_ineq = (-X, Y)                          // CRS: −v·xⱼ + u·yⱼ ≤ 0
    }
    else if (!is_output) {
        A_ineq = (-X, Y, J(n, 1,  1))            // Input VRS:  +w term (u·y+w ≤ v·x)
    }
    else {
        A_ineq = (-X, Y, J(n, 1, -1))            // Output VRS: −w term (u·y−w ≤ v·x)
    }
    row_scale=rowmax(abs(A_ineq))
    for (j=1; j<=n; j++) if (row_scale[j]==0) row_scale[j]=1
    A_ineq=A_ineq:/row_scale
    b_ineq = J(n, 1, 0)

    // ── Bounds: nonnegative v, u; w free ────────────────────────────────────
    lb = !is_vrs ? J(1, nvars, 0) : (J(1, m+s, 0), .)
    ub = J(1, nvars, .)

    // ── Stage 1: Solve self-efficiency ───────────────────────────────
    rc=_c57_solve(c,A_eq,b_eq,A_ineq,b_ineq,lb,ub,tol,sol,obj,d1)
    diag[1..6]=d1

    if (rc==0) {
        v_opt = sol[1..m]:/sx
        u_opt = sol[(m+1)..(m+s)]:/sy
        w_opt = (is_vrs ? sol[nvars] : 0)
        if (!is_output) {
            eff = -obj           // Input:  −(−u·y₀ [−w]) = u·y₀ [+w]
        }
        else {
            eff = (obj>0 & obj<. ? 1/obj : .)   // Output: 1/(v·x₀ [+w])
        }
    }
    else {
        v_opt = J(1, m, .)
        u_opt = J(1, s, .)
        w_opt = .
        eff   = .
        return
    }

    if (eff<=0 | eff>=.) {
        v_opt=J(1,m,.); u_opt=J(1,s,.); w_opt=.; eff=.
        return
    }
    stage_status=1

    // ── Stage 2: Benevolent / Aggressive secondary LP ────────────────
    // Preserve the existing surrogate objectives (not a sum of peer ratios).
    // Input uses aggregate peer slack without the intercept; output uses peer inputs.
    if (method != "arbitrary" && eff < .) {
        stage_status=3
        real scalar eff_prime
        real rowvector c2, A_eq2_row
        real matrix    A_eq2
        real colvector b_eq2
        real rowvector sol2

        eff_prime = eff

        // Second equality: lock self-efficiency at eff_prime
        if (!is_output) {
            // u·y₀ [+w] = eff_prime
            A_eq2_row = !is_vrs ? (J(1,m,0), y0) : (J(1,m,0), y0, 1)
            b_eq2 = (1 \ eff_prime)
        }
        else {
            // v·x₀ [+w] = 1/eff_prime
            A_eq2_row = !is_vrs ? (x0, J(1,s,0)) : (x0, J(1,s,0), 1)
            b_eq2 = (1 \ (1/eff_prime))
        }
        A_eq2 = (A_eq \ A_eq2_row)

        // Secondary objective
        real rowvector sum_y, sum_x
        sum_y = colsum(Y) - y0    // Σ_{j≠k} y_j  (1×s)
        sum_x = colsum(X) - x0    // Σ_{j≠k} x_j  (1×m)

        if (!is_output) {
            // The input secondary goal uses aggregate peer constraint slack
            // d_j = v*x_j - u*y_j.  Benevolent minimizes sum d_j; aggressive
            // maximizes it, while preserving the Stage-1 self-efficiency.
            c2 = !is_vrs ? (sum_x, -sum_y) : (sum_x, -sum_y, 0)
            if (method == "aggressive") c2 = -c2
        }
        else {
            // Benevolent: min v·sum_x  (minimise peer input → higher CE for peers)
            // Aggressive: max v·sum_x  → min −v·sum_x
            c2 = !is_vrs ? (sum_x, J(1,s,0)) : (sum_x, J(1,s,0), 0)
            if (method == "aggressive") c2 = -c2
        }

        rc=_c57_solve(c2,A_eq2,b_eq2,A_ineq,b_ineq,lb,ub,tol,sol2,obj,d2)
        diag[7..12]=d2

        if (rc==0) {
            v_opt = sol2[1..m]:/sx
            u_opt = sol2[(m+1)..(m+s)]:/sy
            w_opt = (is_vrs ? sol2[nvars] : 0)
            stage_status=2
            // eff unchanged — keeps Stage 1 self-efficiency
        }
        // If Stage 2 fails: keep Stage 1 solution (already in v_opt, u_opt, w_opt)
    }
}

// Bounded exact original-data face. The accepted T57 primary remains the
// reported score; no rounded objective is used as a new equality RHS.
void _c58_multiplier(real rowvector x0,real rowvector y0,real matrix X,real matrix Y,
 string scalar rts,string scalar orientation,string scalar method,real scalar tol,
 real rowvector v,real rowvector u,real scalar w,real scalar eff,real scalar stage,
 real rowvector diag,real matrix XR,real matrix YR,real scalar index,
 real rowvector sx,real rowvector sy,real rowvector link,
 real rowvector p_diag,real rowvector s_diag,real scalar map_reason,real rowvector p_trace,real rowvector s_trace,real rowvector objective_diag)
{
 real scalar n,m,s,d,free,out,obj,exacteff
 real rowvector c,c2,x,lo,hi,q
 real matrix E,U,T,G,F,G2
 real colvector b,h,fb,h2,ub
 struct _C58_Cert scalar primary,secondary
 transmorphic scalar exactC
 real scalar exact_exponent
 objective_diag=J(1,6,0)
 p_diag=(0,.,.,.,.,.,.,.,.);s_diag=p_diag;map_reason=.
 p_trace=s_trace=_c62_empty_trace()
 link=(0,.,0,0,.,.,.,.,.,.)
 _c57_multiplier(x0,y0,X,Y,rts,orientation,"arbitrary",tol,v,u,w,eff,stage,diag)
 if(stage==0 | method=="arbitrary") return
 stage=3;link[4]=2
 n=rows(XR);m=cols(XR);s=cols(YR);free=rts=="vrs";out=orientation=="output";d=m+s+free
 if(d+free>16 | n>32) {link[1]=-2;link[4]=6;p_diag[2]=-2;objective_diag=(1,-2,0,3,0,0);return;}
 c=(out ? (XR[index,.],J(1,s,0)) : (J(1,m,0),-YR[index,.]))
 E=(out ? (J(1,m,0),YR[index,.]) : (XR[index,.],J(1,s,0)))
 U=(-XR,YR)
 if(free) {c=(c,(out ? 1 : -1));E=(E,0);U=(U,J(n,1,(out ? -1 : 1)));}
 b=1;ub=J(n,1,0);lo=J(1,d,0);if(free) lo[d]=.;hi=J(1,d,.)
 T=_c58_map(E,U,c,free,map_reason)
 if(rows(T)==0) {link[1]=-1;p_diag[2]=-1;return;}
 _c58_fit(c*T,E*T,b,U*T,ub,tol,primary,G,h)
 p_trace=primary.trace62
 p_diag=_c60_certdiag(primary,rows(E))
 link[1]=primary.status;link[5]=primary.ops;link[7]=primary.solver_status;link[9]=primary.reason
 if(primary.status!=1) return
 obj=_c58_value(c*T,primary)
 exacteff=(out ? (obj>0 & obj<. ? 1/obj : .) : -obj)
 if(exacteff>=. | abs(exacteff-eff)>tol*max((1,abs(exacteff),abs(eff)))) {link[4]=3;return;}
 link[4]=4;F=E*T;fb=b
 if(rows(primary.face)) {F=F\G[primary.face,.];fb=fb\h[primary.face];}
 link[10]=rows(primary.face)
 _c70_build(XR,YR,index,out,method=="aggressive",T,c2,exactC,exact_exponent,objective_diag)
 if(objective_diag[2]!=1) return
 _c70_fit(c2,F,fb,U*T,ub,tol,exactC,exact_exponent,objective_diag,secondary,G2,h2)
 s_trace=secondary.trace62
 s_diag=_c60_certdiag(secondary,rows(F))
 link[2]=secondary.status;link[6]=secondary.ops;link[8]=secondary.solver_status
 diag[7..12]=(91,secondary.solver_status,.,.,.,.)
 if(secondary.status!=1) return
 objective_diag[4]=1
 x=_c58_point(T,secondary);q=_dea_lpcheck(x,E,b,U,ub,lo,hi)
 if(hasmissing(x) | !_dea_lpcheck_ok(q,tol)) {link[4]=5;return;}
 obj=c*x'
 exacteff=(out ? (obj>0 & obj<. ? 1/obj : .) : -obj)
 if(exacteff>=. | abs(exacteff-eff)>tol*max((1,abs(exacteff),abs(eff)))) {link[4]=5;return;}
 v=x[1..m]:*sx;u=x[m+1..m+s]:*sy;w=(free ? x[d] : 0)
 stage=2;link[3]=1;link[4]=1
 diag[7..12]=(0,secondary.solver_status,_c70_value(exactC,exact_exponent,secondary),.,max(q[4..6]),max(q[1..3]))
}

void _c57_compute(string scalar inputs,
                            string scalar outputs,
                            string scalar dmu_var,
                            string scalar touse,
                            string scalar eff_ce_var,
                            string scalar eff_self_var,
                            string scalar eff_maverick_var,
                            string scalar method,
                            string scalar rts,
                            string scalar orientation,
                            string scalar vrscross,
                            real scalar tol,
                            real scalar show_matrix,
                            real scalar show_detail,
                            real scalar show_debug,
                            string scalar ce_matrix_name,
                            string scalar v_name, string scalar u_name,
                            string scalar w_name, string scalar status_name,
                            string scalar peers_name, string scalar diag_name, string scalar link_name,
                            string scalar p_diag_name,string scalar s_diag_name,string scalar map_name,string scalar p_trace_name,string scalar s_trace_name,string scalar objective_name,string scalar raw_name,string scalar scaling_name,string scalar peer_name,|real scalar capture_optima)
{
    // Load data
    real matrix X, Y, XR, YR
    real rowvector sx, sy
    real scalar jj
    real colvector dmu_ids, obs_indices
    string colvector dmu_labels
    string rowvector input_names, output_names

    input_names  = tokens(inputs)
    output_names = tokens(outputs)

    X           = st_data(., input_names,  touse)
    Y           = st_data(., output_names, touse)
    obs_indices = selectindex(st_data(., touse))

    // Read DMU — handles string and numeric
    _c57_dmu_read(dmu_var, touse, dmu_ids, dmu_labels)

    real scalar n, m, s
    n = rows(X)
    m = cols(X)
    s = cols(Y)
    // A change of measurement units leaves the normalized LP unchanged.
    sx=colmax(X); sy=colmax(Y)
    for (jj=1; jj<=m; jj++) if (sx[jj]==0) sx[jj]=1
    for (jj=1; jj<=s; jj++) if (sy[jj]==0) sy[jj]=1
    XR=X;YR=Y
    X=X:/sx; Y=Y:/sy

    if (show_debug) {
        printf("  n=%g, m=%g (inputs), s=%g (outputs)\n", n, m, s)
    }

    // Step 1: Compute optimal weights for each DMU using multiplier model
    // For input-oriented CRS: max sum(u*y0) s.t. sum(v*x0)=1, sum(u*yj)-sum(v*xj)<=0

    real matrix V_raw,U_raw,RD,SC,PS
    real rowvector raw_x,rd_i,scaling_i,peer_diag
    real scalar raw_stage,raw_score,raw_supported
    V_raw=J(n,m,.);U_raw=J(n,s,.);RD=J(n,17,0)
    raw_supported=(n<=32 & m+s+2*(rts=="vrs")<=16)
    SC=J(n,(raw_supported ? m+s+(rts=="vrs")+n+2 : 0),.)
    PS=J((raw_supported | show_matrix ? n : 1),(raw_supported | show_matrix ? n : 1),0)
    real matrix V_opt, U_opt   // Optimal weights: V (inputs), U (outputs)
    real colvector W_opt       // VRS free intercept (0 for CRS)
    V_opt = J(n, m, .)
    U_opt = J(n, s, .)
    W_opt = J(n, 1, 0)

    real colvector self_eff, stage_status, peer_count
    self_eff = J(n, 1, .)
    stage_status=J(n,1,0); peer_count=J(n,1,0)

    real scalar i, j
    real rowvector x0, y0, v_i, u_i
    real scalar eff_i, w_i, stage_i
    real matrix recovery_p,recovery_s
    real rowvector pt_i,st_i
    recovery_p=recovery_s=J(n,36,0)
    real matrix cert_p,cert_s
    real rowvector pc_i,sc_i
    real colvector maps
    real scalar map_i
    cert_p=cert_s=J(n,9,.);maps=J(n,1,.)
    real matrix link_diag
    real rowvector link_i
    link_diag=J(n,10,.)
    real matrix solver_diag
    real rowvector d_i,objective_i
    real matrix objective_diag
    objective_diag=J(n,6,0)
    solver_diag=J(n,12,.)

    for (i = 1; i <= n; i++) {
        x0 = X[i, .]
        y0 = Y[i, .]

        // Solve multiplier model (Stage 1 + optional Stage 2 for benevolent/aggressive)
        _c58_multiplier(x0,y0,X,Y,rts,orientation,method,tol,v_i,u_i,w_i,eff_i,stage_i,d_i,XR,YR,i,sx,sy,link_i,pc_i,sc_i,map_i,pt_i,st_i,objective_i)
        V_raw[i,.]=v_i:/sx;U_raw[i,.]=u_i:/sy
        if(stage_i==0) {
            _c84_multiplier(XR,YR,i,rts=="vrs",orientation=="output",method,tol,raw_x,raw_score,raw_stage,rd_i,scaling_i,pc_i,sc_i,pt_i,st_i,link_i,objective_i)
            RD[i,.]=rd_i;if(raw_supported) SC[i,.]=scaling_i
            if(raw_stage>0) {
                V_raw[i,.]=raw_x[1..m];U_raw[i,.]=raw_x[m+1..m+s]
                w_i=(rts=="vrs" ? raw_x[m+s+1] : 0)
                eff_i=raw_score;stage_i=raw_stage
            }
        }
        cert_p[i,.]=pc_i;cert_s[i,.]=sc_i;maps[i]=map_i
        recovery_p[i,.]=pt_i;recovery_s[i,.]=st_i
        link_diag[i,.]=link_i
        solver_diag[i,.]=d_i
        objective_diag[i,.]=objective_i

        V_opt[i, .] = v_i
        U_opt[i, .] = u_i
        W_opt[i]    = w_i
        self_eff[i] = eff_i
        stage_status[i]=stage_i

        if (show_debug) {
            printf("  DMU %g: self-eff = %9.4f\n", dmu_ids[i], eff_i)
        }
    }

    // Step 2: Compute cross-efficiency matrix
    // CE[i,j] = efficiency of DMU j using weights of DMU i
    real matrix CE
    CE = J(n, n, .)

    real rowvector v_i2, u_i2
    real scalar num, denom

    // Legacy: u*y/(v*x). Lim-Zhu (2015), equation (4), input VRS:
    // u*y/(v*x-w), where our input multiplier w equals -xi in the paper.
    // Apply the chosen formula to the diagonal too; do not replace it by self.
    for (i = 1; i <= n; i++) {
        v_i2 = V_opt[i, .]
        u_i2 = U_opt[i, .]

        for (j = 1; j <= n; j++) {
            num   = u_i2 * Y[j, .]'
            denom = v_i2 * X[j, .]'
            // CRS permits genuinely small positive peer denominators. Under
            // the legacy VRS ratio convention retain its tolerance cutoff:
            // the intercept can carry the entire output-oriented self value.
            if (vrscross=="limzhu") CE[i,j]=_c57_lzratio(denom,num,W_opt[i],tol)
            else CE[i, j] = (denom>(rts=="vrs" ? tol : 0) & denom<. & num>=0 & num<. ? num/denom : .)
            if(RD[i,6]==1) {
                CE[i,j]=_c84_peer(V_raw[i,.],XR[j,.],U_raw[i,.],YR[j,.],W_opt[i],rts=="vrs",vrscross=="limzhu",tol,peer_diag)
                PS[i,j]=peer_diag[1]
            }
        }
    }

    // Step 3: Compute cross-efficiency scores (column means excluding self)
    real colvector ce_scores
    ce_scores = J(n, 1, .)

    real scalar sum_ce, count_ce

    for (j = 1; j <= n; j++) {
        sum_ce = 0
        count_ce = 0

        for (i = 1; i <= n; i++) {
            if (i != j && CE[i, j] < .) {
                sum_ce = sum_ce + CE[i, j]
                count_ce++
            }
        }

        peer_count[j]=count_ce
        if (count_ce > 0) {
            ce_scores[j] = sum_ce / count_ce
        }
    }

    // Step 4: Compute maverick index (self-appraisal vs cross-efficiency)
    real colvector maverick
    maverick = J(n,1,.)
    for (i=1; i<=n; i++) {
        if (self_eff[i]<. & ce_scores[i]>0 & ce_scores[i]<.) maverick[i]=self_eff[i]/ce_scores[i]
    }

    // Compute max label width for aligned display
    real scalar lw
    string scalar fmt_dmu
    lw = 10
    for (i = 1; i <= n; i++) {
        if (strlen(dmu_labels[i]) > lw) lw = strlen(dmu_labels[i])
    }
    fmt_dmu = "%" + strofreal(lw) + "s"

    // Display results — use string labels
    printf("\n")
    printf("Cross-Efficiency Results\n")
    printf("{hline 70}\n")
    printf(fmt_dmu + " | %10s %10s %10s %10s\n",
           "DMU", "Self", "Cross-Eff", "Maverick", "Rank")
    printf(fmt_dmu + " | %10s %10s %10s %10s\n",
           "-":*lw, "----------", "----------", "----------", "----------")

    // Rank by cross-efficiency
    real colvector ranks
    ranks = _c57_rank(ce_scores)

    for (i = 1; i <= n; i++) {
        printf(fmt_dmu + " | %10.4f %10.4f %10.4f %10.0f\n",
               dmu_labels[i], self_eff[i], ce_scores[i], maverick[i], ranks[i])
    }
    printf("{hline 70}\n")

    // Show cross-efficiency matrix if requested
    if (show_matrix) {
        printf("\nCross-Efficiency Matrix (row i evaluates column j)\n")
        printf("{hline 70}\n")

        // Header
        printf(fmt_dmu + " |", "Evaluator")
        for (j = 1; j <= n; j++) {
            printf(" %7s", dmu_labels[j])
        }
        printf("\n")

        printf(fmt_dmu + " |", "-":*lw)
        for (j = 1; j <= n; j++) {
            printf(" -------")
        }
        printf("\n")

        // Data
        for (i = 1; i <= n; i++) {
            printf(fmt_dmu + " |", dmu_labels[i])
            for (j = 1; j <= n; j++) {
                printf(" %7.4f", CE[i, j])
            }
            printf("\n")
        }
        printf("{hline 70}\n")
        st_matrix(objective_name,objective_diag)
        st_matrix(ce_matrix_name, CE)
        st_matrix(v_name,V_raw)
        st_matrix(u_name,U_raw)
        st_matrix(w_name,W_opt)
        st_matrix(status_name,stage_status)
        st_matrix(peers_name,peer_count)
        st_matrix(diag_name,solver_diag)
        st_matrix(link_name,link_diag)
        st_matrix(p_diag_name,cert_p);st_matrix(s_diag_name,cert_s);st_matrix(map_name,maps)
        st_matrix(p_trace_name,recovery_p[.,(1..13,16..24,26..34)]);st_matrix(s_trace_name,recovery_s[.,(1..13,16..24,26..34)])
    }

    if (!show_matrix & capture_optima==1) {
        st_matrix(v_name,V_raw);st_matrix(u_name,U_raw)
        st_matrix(w_name,W_opt);st_matrix(status_name,stage_status)
    }

    if (show_detail) {
        printf("\nWeight stages: 0=primary unavailable, 1=primary only, 2=secondary accepted, 3=primary fallback\n")
        for (i=1; i<=n; i++) printf(fmt_dmu+" | stage %g | available peers %g\n",dmu_labels[i],stage_status[i],peer_count[i])
    }
    if(show_matrix) {
        st_matrix(raw_name,(obs_indices,RD));st_matrix(scaling_name,(obs_indices,SC));st_matrix(peer_name,PS)
    }
    if(show_detail) {
        for(i=1;i<=n;i++) if(RD[i,1]) printf("  Raw recovery: observation %g, DMU %s, stop %g, stage %g; primary certificate %g, secondary certificate %g\n",obs_indices[i],dmu_labels[i],RD[i,11],stage_status[i],RD[i,4],RD[i,7])
    }
    st_numscalar("r(n_raw_attempted)",sum(RD[.,1]))
    st_numscalar("r(n_raw_recovered)",sum(RD[.,6]))
    st_numscalar("r(n_raw_unavailable)",sum(RD[.,1])-sum(RD[.,6]))
    st_numscalar("r(n_raw_primary_unknown)",sum(RD[.,11]:==7))
    st_numscalar("r(n_raw_export_obstructed)",sum(RD[.,17]))
    st_numscalar("r(n_raw_secondary_unknown)",sum((RD[.,11]:==3):|(RD[.,11]:==11):|(RD[.,11]:==12)))
    st_numscalar("r(n_raw_secondary_unbounded)",sum(RD[.,11]:==4))
    st_numscalar("r(n_raw_peers_evaluated)",sum(PS:>0))
    st_numscalar("r(n_raw_peers_unavailable)",sum(PS:>1))
    st_numscalar("r(n_raw_peers_range)",sum(PS:==3))
    st_numscalar("r(n_raw_peers_arithmetic)",sum(PS:==4))
    _c70_counts(objective_diag)
    if(show_detail) _c70_detail(objective_diag,obs_indices,dmu_labels)
    if(show_detail) _c60_cert_detail(cert_p,cert_s,obs_indices)
    if(show_detail) _c62_recovery_detail(recovery_p,recovery_s,obs_indices)
    st_numscalar("r(n_recovery_requested)",sum(recovery_p[.,1])+sum(recovery_s[.,1]))
    st_numscalar("r(n_recovery_finite)",sum(recovery_p[.,2]:==1)+sum(recovery_s[.,2]:==1))
    st_numscalar("r(n_recovery_unknown)",sum(recovery_p[.,2]:==3)+sum(recovery_s[.,2]:==3))
    st_numscalar("r(n_recovery_skipped)",sum(recovery_p[.,2]:==-1)+sum(recovery_s[.,2]:==-1))
    st_numscalar("r(n_recovery_equality)",sum(recovery_p[.,8])+sum(recovery_s[.,8]))
    st_numscalar("r(n_recovery_auxiliary)",sum(recovery_p[.,10])+sum(recovery_s[.,10]))
    st_numscalar("r(n_exact_unbounded)",sum(recovery_p[.,16])+sum(recovery_s[.,16]))
    st_numscalar("r(recovery_integer_ops)",sum(recovery_p[.,22]-recovery_p[.,17])+sum(recovery_s[.,22]-recovery_s[.,17]))
    st_numscalar("r(recovery_steps)",sum(recovery_p[.,13]-recovery_p[.,7])+sum(recovery_s[.,13]-recovery_s[.,7]))
    st_numscalar("r(certificate_integer_ops)",sum(recovery_p[.,22])+sum(recovery_s[.,22]))
    st_numscalar("r(n_cert_primary_attempted)",sum(cert_p[.,1]))
    st_numscalar("r(n_cert_primary_verified)",sum(cert_p[.,2]:==1))
    st_numscalar("r(n_cert_secondary_attempted)",sum(cert_s[.,1]))
    st_numscalar("r(n_cert_secondary_verified)",sum(cert_s[.,2]:==1))
    st_numscalar("r(n_cert_repair_eq_stops)",sum(cert_p[.,4]:==3)+sum(cert_s[.,4]:==3))
    st_numscalar("r(n_cert_repair_both_stops)",sum(cert_p[.,4]:==2)+sum(cert_s[.,4]:==2))
    st_numscalar("r(n_cert_mapping_failed)",sum((maps:>0):&(maps:<.)))
    st_numscalar("r(n_failed_primary)",sum(stage_status:==0))
    st_numscalar("r(n_failed_secondary)",sum(stage_status:==3))
    st_numscalar("r(n_failed_ce)",sum(ce_scores:>=.))
    st_numscalar("r(n_missing_peer_pairs)",n*(n-1)-sum(peer_count))
    st_numscalar("r(n_missing_diagonal)",sum(diagonal(CE):>=.))

    // Store results in Stata
    for (i = 1; i <= n; i++) {
        st_store(obs_indices[i], st_varindex(eff_ce_var), ce_scores[i])
        st_store(obs_indices[i], st_varindex(eff_self_var), self_eff[i])
        st_store(obs_indices[i], st_varindex(eff_maverick_var), maverick[i])
    }
}

end
