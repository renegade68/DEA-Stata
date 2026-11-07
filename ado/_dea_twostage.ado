*! version 1.0.0  02nov2026
*! Shared two-stage LP construction; CRS shared weights and independent stages.
version 16.0
capture mata: assert(_dea_nonte_lp_ready()==1)
if _rc {
    quietly findfile _dea_nonte_lp.ado
    quietly run "`r(fn)'"
}
capture mata: mata drop _dea_ts_ready()
capture mata: mata drop _dea_ts_stage()
capture mata: mata drop _dea_ts_joint()
capture mata: mata drop _dea_ts_link()
capture mata: mata drop _dea_ts_compute()
mata:
mata set matastrict on
real scalar _dea_ts_ready() return(1)

// Independent radial envelopment. No automatic switch to another engine.
real scalar _dea_ts_stage(real matrix X, real matrix Y, real scalar k,
    string scalar rts, string scalar ort, real scalar tol, real rowvector lambda,
    real scalar cap)
{
    real scalar n,m,s,eff
    real matrix A,E
    real colvector b,h,scale
    real rowvector c,lb,ub
    class LinearProgram scalar S
    n=rows(X);m=cols(X);s=cols(Y);lambda=J(1,n,.)
    E=J(0,n+1,0);h=J(0,1,0)
    if(rts=="vrs") {E=(J(1,n,1),0);h=1;}
    if(ort=="input") {
        c=(J(1,n,0),1);A=(X',-X[k,.]')\(-Y',J(s,1,0))
        b=J(m,1,0)\(-Y[k,.]');lb=J(1,n+1,0)
    }
    else {
        c=(J(1,n,0),-1);A=(X',J(m,1,0))\(-Y',Y[k,.]')
        b=X[k,.]'\J(s,1,0);lb=(J(1,n,0),1)
    }
    scale=rowmax(abs((A,b)));scale=scale:+(scale:==0)
    A=A:/scale;b=b:/scale;ub=J(1,n+1,.)
    S=LinearProgram();S.setMaxOrMin("min");S.setCoefficients(c)
    S.setInequality(A,b);if(rows(E)) S.setEquality(E,h)
    S.setBounds(lb,ub);S.setTrace("off")
    if(_dea_nonte_optimize(S,tol,cap)) return(.)
    eff=(ort=="input" ? S.value() : -S.value())
    if(eff>=. | eff<0) return(.)
    if(ort=="output") {if(eff<=0) return(.);eff=1/eff;}
    if(eff>1+tol) return(.)
    lambda=S.parameters();lambda=lambda[1..n]
    return(eff)
}

// With a single active intermediate, common weighting imposes no tradeoff
// between stage optima (Liang et al. 2008). Normalize w*Zk=1 and solve the
// two independent multiplier sides. This is an equivalent model reduction,
// not substitution after a failed joint LP. The caller still checks the full
// leader-locked constraints at the requested tolerance.
real scalar _dea_ts_link(real matrix X,real matrix Z,real matrix Y,
    real scalar k,real scalar q,real scalar tol,real rowvector w)
{
    real scalar n,m,d,s,j,den
    real rowvector v,u,zweight,c,active
    real colvector rhs,scale
    real matrix A
    class LinearProgram scalar S
    n=rows(X);m=cols(X);d=cols(Z);s=cols(Y)
    active=selectindex(colmax(Z):>0);j=active[1]
    zweight=J(1,d,0);zweight[j]=1/Z[k,j];rhs=Z*zweight'
    A=-X;scale=rowmax(abs((A,rhs)));scale=scale:+(scale:==0)
    S=LinearProgram();S.setMaxOrMin("min");S.setCoefficients(X[k,.])
    S.setInequality(A:/scale,-rhs:/scale);S.setBounds(J(1,m,0),J(1,m,.));S.setTrace("off")
    if(_dea_nonte_optimize(S,tol,1e-12)) return(1)
    v=S.parameters()
    A=Y;scale=rowmax(abs((A,rhs)));scale=scale:+(scale:==0)
    S=LinearProgram();S.setMaxOrMin("min");S.setCoefficients(-Y[k,.])
    S.setInequality(A:/scale,rhs:/scale);S.setBounds(J(1,s,0),J(1,s,.));S.setTrace("off")
    if(_dea_nonte_optimize(S,tol,1e-12)) return(1)
    u=S.parameters();w=(v,zweight,u)
    if(q==3) {den=v*X[k,.]';if(den<=0 | den>=.) return(1);w=w/den;}
    return(0)
}

// q=0: centralized; q=1: maximize stage 1 on the centralized optimal face;
// q=2/3: follower with stage 1/2 leader's independently optimal score fixed.
// Columns are (v,w,u). Only the two stage technologies are needed: the
// whole-system inequality follows by adding them. All weights are nonnegative.
real scalar _dea_ts_joint(real matrix X0,real matrix Z0,real matrix Y0,
    real scalar k,real scalar q,real scalar lock,real scalar tol,
    real rowvector weights,real scalar route)
{
    real scalar n,m,d,s,p
    real matrix X,Z,Y,A,E
    real colvector b,h,scale
    real rowvector sx,sz,sy,c,ex,ez,ey,lb,ub,candidate
    class LinearProgram scalar S
    n=rows(X0);m=cols(X0);d=cols(Z0);s=cols(Y0);p=m+d+s
    route=1
    weights=J(1,p,.)
    sx=X0[k,.];sz=Z0[k,.];sy=Y0[k,.]
    sx=sx:+(sx:==0);sz=sz:+(sz:==0);sy=sy:+(sy:==0)
    X=X0:/sx;Z=Z0:/sz;Y=Y0:/sy
    ex=(X[k,.],J(1,d+s,0));ez=(J(1,m,0),Z[k,.],J(1,s,0))
    ey=(J(1,m+d,0),Y[k,.])
    A=(-X,Z,J(n,s,0))\(J(n,m,0),-Z,Y);b=J(2*n,1,0)
    if(q==0) {c=-ey;E=ex;h=1;}
    else if(q==1) {c=-ez;E=ex\ey;h=1\lock;}
    else if(q==2) {c=-ey;E=ez\(ez-lock*ex);h=1\0;}
    else {c=-ez;E=ex\(ey-lock*ez);h=1\0;}
    scale=rowmax(abs((A,b)));scale=scale:+(scale:==0);A=A:/scale;b=b:/scale
    scale=rowmax(abs((E,h)));scale=scale:+(scale:==0);E=E:/scale;h=h:/scale
    lb=J(1,p,0);ub=J(1,p,.)
    if(q>=2 & sum(colmax(Z):>0)==1) {
        route=2
        if(_dea_ts_link(X,Z,Y,k,q,tol,candidate)) return(1)
        if(!_dea_lpcheck_ok(_dea_lpcheck(candidate,E,h,A,b,lb,ub),tol)) return(1)
        if(-c*candidate'<0 | -c*candidate'>1+tol | hasmissing(candidate)) return(1)
        weights=candidate:/(sx,sz,sy)
        return(0)
    }
    S=LinearProgram();S.setMaxOrMin("min");S.setCoefficients(c)
    S.setInequality(A,b);S.setEquality(E,h);S.setBounds(lb,ub);S.setTrace("off")
    // Keep the centralized optimum precise before locking its decomposition.
    // Follower solves use the requested residual test with a 1e-8 convergence
    // cap; a blanket 1e-12 cap can stall on the fixed leader face.
    if(_dea_nonte_optimize(S,tol,(q<2 ? 1e-12 : 1e-8))) return(1)
    if(-S.value()<0 | -S.value()>1+tol | S.value()>=.) return(1)
    weights=S.parameters():/(sx,sz,sy)
    return(0)
}

// mode: relational, product, mean, independent ordering, or noncooperative.
void _dea_ts_compute(string scalar xs,string scalar zs,string scalar ys,
    string scalar id,string scalar touse,string scalar gen,string scalar kind,
    string scalar mode,string scalar leader,string scalar rts,string scalar ort,
    real scalar tol,real scalar matrices,real scalar detail,real scalar debug,
    string scalar prefix)
{
    real matrix X,Z,Y,R,W,L1,L2,B,Status
    real rowvector sx,sz,sy,w,alt,lam,meta
    real colvector obs,decomp,joint,routes
    string scalar fmt
    string colvector labels
    real scalar n,m,d,s,i,rc,den,mid,num,e1,e2,overall,lead,follow,lw,route
    X=st_data(.,tokens(xs),touse);Z=st_data(.,tokens(zs),touse);Y=st_data(.,tokens(ys),touse)
    obs=selectindex(st_data(.,touse));n=rows(X);m=cols(X);d=cols(Z);s=cols(Y)
    if(substr(st_vartype(id),1,3)=="str") labels=st_sdata(.,id,touse)
    else labels=strtrim(strofreal(st_data(.,id,touse),"%21.0g"))
    sx=colmax(X);sz=colmax(Z);sy=colmax(Y)
    sx=sx:+(sx:==0);sz=sz:+(sz:==0);sy=sy:+(sy:==0)
    X=X:/sx;Z=Z:/sz;Y=Y:/sy
    R=J(n,3,.);B=J(n,2,.);W=J(n,m+d+s,.);L1=L2=J(n,n,.)
    decomp=joint=J(n,1,.);routes=J(n,1,0)
    for(i=1;i<=n;i++) {
        if(mode=="relational") {
            decomp[i]=0
            rc=_dea_ts_joint(X,Z,Y,i,0,.,tol,w,route)
            if(rc==0) {
                overall=w[(m+d+1)..(m+d+s)]*Y[i,.]'
                R[i,1]=overall;decomp[i]=1
                if(_dea_ts_joint(X,Z,Y,i,1,overall,tol,alt,route)==0) {w=alt;decomp[i]=2;}
                den=w[1..m]*X[i,.]';mid=w[(m+1)..(m+d)]*Z[i,.]'
                num=w[(m+d+1)..(m+d+s)]*Y[i,.]'
                e1=(den>0 & den<. ? mid/den : .)
                e2=(mid>0 & mid<. ? num/mid : .)
                if(e1>=0 & e1<=1+tol) R[i,2]=e1
                if(e2>=0 & e2<=1+tol) R[i,3]=e2
                W[i,.]=w
            }
        }
        else {
            B[i,1]=_dea_ts_stage(X,Z,i,rts,ort,tol,lam,(mode=="noncooperative" ? 1e-12 : 1e-8));L1[i,.]=lam
            B[i,2]=_dea_ts_stage(Z,Y,i,rts,ort,tol,lam,(mode=="noncooperative" ? 1e-12 : 1e-8));L2[i,.]=lam
            R[i,2..3]=B[i,.]
            if(mode=="noncooperative") {
                lead=(leader=="stage1" ? B[i,1] : B[i,2]);follow=.
                joint[i]=0
                if(lead>0 & lead<.) {
                    joint[i]=1
                    rc=_dea_ts_joint(X,Z,Y,i,(leader=="stage1" ? 2 : 3),lead,tol,w,route)
                    routes[i]=route
                    if(rc==0) {
                        den=w[1..m]*X[i,.]';mid=w[(m+1)..(m+d)]*Z[i,.]'
                        num=w[(m+d+1)..(m+d+s)]*Y[i,.]'
                        follow=(leader=="stage1" ? (mid>0 & mid<. ? num/mid : .) : (den>0 & den<. ? mid/den : .))
                        if(follow>=0 & follow<=1+tol) {joint[i]=2;W[i,.]=w;}
                        else follow=.
                    }
                }
                if(leader=="stage1") R[i,3]=follow
                else R[i,2]=follow
            }
            if(!hasmissing(R[i,2..3])) {
                R[i,1]=(mode=="mean" ? (R[i,2]+R[i,3])/2 : R[i,2]*R[i,3])
            }
        }
        if(debug) printf("  observation %g: overall=%g, stage1=%g, stage2=%g\n",obs[i],R[i,1],R[i,2],R[i,3])
    }
    Status=R:<.;lw=10
    for(i=1;i<=n;i++) lw=max((lw,strlen(labels[i])))
    fmt="%"+strofreal(lw)+"s"
    printf("\n"+fmt+" | %10s %10s %10s\n","DMU","Stage 1","Stage 2","Overall")
    for(i=1;i<=n;i++) {
        printf(fmt+" | %10.6f %10.6f %10.6f\n",labels[i],R[i,2],R[i,3],R[i,1])
        if(detail) {
            printf("  available (overall/stage1/stage2): %g/%g/%g",Status[i,1],Status[i,2],Status[i,3])
            if(mode=="relational") printf("; decomposition=%g",decomp[i])
            if(mode=="noncooperative") printf("; follower solve=%g, route=%g",joint[i],routes[i])
            if(kind=="stack") printf("; independent stage1=%g, stage2=%g",B[i,1],B[i,2])
            printf("\n")
        }
    }
    if(sum(R[,1]:>=.)) printf("Note: %g overall score(s) unavailable; available components retained.\n",sum(R[,1]:>=.))
    if(sum(decomp:==1)) printf("Note: %g decomposition LP(s) unavailable; primary weights retained.\n",sum(decomp:==1))
    if(sum(joint:==1)) printf("Note: %g follower LP(s) unavailable; no independent-score substitution.\n",sum(joint:==1))
    if(kind=="network") {
        st_store(obs,st_varindex(gen+"_overall"),R[,1]);st_store(obs,st_varindex(gen+"_stage1"),R[,2]);st_store(obs,st_varindex(gen+"_stage2"),R[,3])
    }
    else {
        st_store(obs,st_varindex(gen+"_overall"),R[,1])
        st_store(obs,st_varindex(gen+"_leader"),R[,(leader=="stage1" ? 2 : 3)])
        st_store(obs,st_varindex(gen+"_follower"),R[,(leader=="stage1" ? 3 : 2)])
        st_store(obs,st_varindex(gen+"_indep1"),B[,1]);st_store(obs,st_varindex(gen+"_indep2"),B[,2])
    }
    meta=(colsum(R:>=.),sum(decomp:==1),sum(joint:==1),colsum(B:>=.),sum(routes:==2))
    st_matrix(prefix+"_meta",meta)
    if(matrices) {
        st_matrix(prefix+"_scores",R);st_matrix(prefix+"_status",Status)
        if(mode=="relational" | mode=="noncooperative") {
            W=W:/(sx,sz,sy);st_matrix(prefix+"_v",W[,1..m]);st_matrix(prefix+"_w",W[,(m+1)..(m+d)]);st_matrix(prefix+"_u",W[,(m+d+1)..(m+d+s)])
            st_matrix(prefix+"_joint",(mode=="relational" ? decomp : joint))
        }
        if(mode=="noncooperative") st_matrix(prefix+"_route",routes)
        if(mode!="relational") {st_matrix(prefix+"_l1",L1);st_matrix(prefix+"_l2",L2);}
    }
}
end
