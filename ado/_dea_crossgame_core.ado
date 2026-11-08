*! version 1.0.0 08nov2026
*! Explicit core load before the recovery dispatcher.
version 16.0
capture mata: mata drop _g114_core_ready()
capture mata: mata drop _g114_ready()
// Reloads must replace old helper definitions in the same Stata session.
capture mata: mata drop _dea_cg_dmu_read()
capture mata: mata drop _dea_crossgame_compute()
capture mata: mata drop _dea_cg_multiplier()
capture mata: mata drop _dea_cg_aggressive()
capture mata: mata drop _dea_cg_benevolent()
capture mata: mata drop _dea_cg_nash()
capture mata: mata drop _dea_cg_lp()
capture mata: mata drop _dea_cg_peers()
capture mata: mata drop _dea_cg_average()
capture mata: mata drop _dea_cg_show()

mata:
mata set matastrict on
real scalar _g114_core_ready() return(11401)
// mode: 0 primary, 1 aggressive numerator surrogate, 2 benevolent surrogate,
    // 3 Liang et al. (2008) model (5), with protected peer d and score theta.
real scalar _dea_cg_lp(real matrix X0, real matrix Y0, real scalar j,
    real scalar mode, real scalar own, real scalar d, real scalar theta,
    real scalar tol, real rowvector v, real rowvector u, real scalar value, | real rowvector nd)
{
    real matrix X,Y,E,A
    real rowvector sx,sy,c,sol
    real colvector b,rs
    real scalar k,m,s,rc
    class LinearProgram scalar S
    nd=(.,.,.,1)
    m=cols(X0);s=cols(Y0);sx=X0[j,.];sy=Y0[j,.]
    v=J(1,m,.);u=J(1,s,.);value=.
    for(k=1;k<=m;k++) if(sx[k]<=0) sx[k]=1
    for(k=1;k<=s;k++) if(sy[k]<=0) sy[k]=1
    X=X0:/sx;Y=Y0:/sy
    if(hasmissing(X) | hasmissing(Y)) return(1)
    E=(X[j,.],J(1,s,0));b=1
    A=(-X,Y)
    if(mode==1 | mode==2) {
        if(own>=.) return(1)
        E=E\(J(1,m,0),Y[j,.]);b=b\own
        c=(J(1,m,0),colsum(Y)-Y[j,.])
        if(mode==2) c=-c
    }
    else c=(J(1,m,0),-Y[j,.])
    if(mode==3) {
        if(theta>=. | theta<0) return(1)
        A=A\(theta*X[d,.],-Y[d,.])
    }
    rs=rowmax(abs(A))
    for(k=1;k<=rows(rs);k++) if(rs[k]==0) rs[k]=1
    A=A:/rs
    S=LinearProgram();S.setMaxOrMin("min");S.setCoefficients(c)
    S.setEquality(E,b);S.setInequality(A,J(rows(A),1,0))
    S.setBounds(J(1,m+s,0),J(1,m+s,.));S.setTrace("off")
    rc=_dea_nonte_optimize(S,tol,1e-12)
    nd=(S.converged(),S.errorcode(),S.value(),1)
    if(rc) return(1)
    sol=S.parameters();value=S.value()
    if(mode==0 | mode==3) {
        value=-value
        if(value<0 | value>1+tol | value>=.) {
            value=.;return(1)
        }
    }
    nd[4]=0
    v=sol[1..m]:/sx;u=sol[(m+1)..(m+s)]:/sy
    return(0)
}

real matrix _dea_cg_peers(real matrix X,real matrix Y,real matrix V,
    real matrix U,real scalar tol)
{
    real matrix CE
    real scalar i,j,num,den,q,n
    n=rows(X);CE=J(n,n,.)
    for(i=1;i<=n;i++) for(j=1;j<=n;j++) {
        num=U[i,.]*Y[j,.]';den=V[i,.]*X[j,.]'
        if(num>=0 & num<. & den>0 & den<.) {
            q=num/den
            if(q<. & q<=1+tol) CE[i,j]=q
        }
    }
    return(CE)
}

real colvector _dea_cg_average(real matrix CE,real scalar exclude,
    real colvector counts)
{
    real colvector scores
    real scalar n,i,j,total
    n=rows(CE);counts=J(n,1,0);scores=J(n,1,.)
    for(j=1;j<=n;j++) {
        total=0
        for(i=1;i<=n;i++) if((!exclude | i!=j) & CE[i,j]<.) {
            total=total+CE[i,j];counts[j]=counts[j]+1
        }
        if(counts[j]) scores[j]=total/counts[j]
    }
    return(scores)
}

void _dea_cg_show(real matrix CE,string colvector labels,string scalar title)
{
    real scalar i,j
    printf("\n%s (row d, column j; selected observation order)\n",title)
    printf("%12s", "d / j")
    for(j=1;j<=rows(CE);j++) printf(" %9.0f",j)
    printf("\n")
    for(i=1;i<=rows(CE);i++) {
        printf("%12s",labels[i])
        for(j=1;j<=cols(CE);j++) printf(" %9.5f",CE[i,j])
        printf("\n")
    }
}

void _dea_crossgame_compute(string scalar inputs,string scalar outputs,
    string scalar dmu_var,string scalar touse,string scalar gen_a,string scalar gen_b,string scalar gen_n,
    string scalar strategy,real scalar maxiter,real scalar tol,
    real scalar show_matrix,real scalar show_detail,real scalar show_debug,
    string scalar prefix)
{
    real matrix X,Y,V,U,Va,Ua,Vb,Ub,C0,Ca,Cb,Cn,Cnew,hist,trace,XR,YR,D
    real rowvector sx,sy,v,u,meta,counts
    real colvector obs_indices,self,sa,sb,pa,pb,dummy,a,b,theta,next,sourcea,sourceb
    string colvector labels
    real scalar n,m,s,i,j,d,k,rc,value,ns,iterations,attempts,failed,delta,initmissing
    real scalar wanta,wantb,wantn,source,row0
    D=J(0,32,.);counts=J(1,13,0)
    X = st_data(., tokens(inputs), touse)
    Y=st_data(.,tokens(outputs),touse);obs_indices=selectindex(st_data(.,touse))
    if(substr(st_vartype(dmu_var),1,3)=="str") labels=st_sdata(.,dmu_var,touse)
    else labels=strtrim(strofreal(st_data(.,dmu_var,touse),"%21.17g"))
    n=rows(X);m=cols(X);s=cols(Y)
    sx=colmax(X);sy=colmax(Y)
    for(k=1;k<=m;k++) if(sx[k]==0) sx[k]=1
    for(k=1;k<=s;k++) if(sy[k]==0) sy[k]=1
    XR=X;YR=Y
    X=X:/sx;Y=Y:/sy
    wanta=(strategy=="all" | strategy=="aggressive")
    wantb=(strategy=="all" | strategy=="benevolent")
    wantn=(strategy=="all" | strategy=="nash")
    self=J(n,1,.);V=J(n,m,.);U=J(n,s,.)
    for(i=1;i<=n;i++) {
        rc=_g114_dispatch(X,Y,XR,YR,sx,sy,i,0,.,.,.,tol,v,u,value,0,obs_indices,D,counts,source)
        if(!rc) {
            self[i]=value;V[i,.]=v;U[i,.]=u
        }
    }
    C0=_dea_cg_peers(X,Y,V,U,tol)
    Va=V;Ua=U;Vb=V;Ub=U;sa=J(n,1,0);sb=J(n,1,0)
    a=J(n,1,.);b=a;pa=J(n,1,0);pb=pa;sourcea=sourceb=J(n,1,0)
    for(i=1;i<=n;i++) if(self[i]<.) {
        if(wanta) {
            rc=_g114_dispatch(X,Y,XR,YR,sx,sy,i,1,self[i],.,.,tol,v,u,value,0,obs_indices,D,counts,source)
            sa[i]=(rc ? 3 : 2);sourcea[i]=(rc ? 3 : source)
            if(!rc) {
                Va[i,.]=v;Ua[i,.]=u
            }
        }
        if(wantb) {
            rc=_g114_dispatch(X,Y,XR,YR,sx,sy,i,2,self[i],.,.,tol,v,u,value,0,obs_indices,D,counts,source)
            sb[i]=(rc ? 3 : 2);sourceb[i]=(rc ? 3 : source)
            if(!rc) {
                Vb[i,.]=v;Ub[i,.]=u
            }
        }
    }
    if(wanta) {
        Ca=_dea_cg_peers(X,Y,Va,Ua,tol);a=_dea_cg_average(Ca,1,pa)
    }
    if(wantb) {
        Cb=_dea_cg_peers(X,Y,Vb,Ub,tol);b=_dea_cg_average(Cb,1,pb)
    }
    ns=.;iterations=0;attempts=0;failed=0;delta=.;initmissing=sum(C0:>=.)
    theta=J(n,1,.);Cn=J(n,n,.)
    if(wantn) {
        // Never run a different game by silently averaging an incomplete matrix.
        ns=2
        if(!hasmissing(C0)) {
            Cn=C0;theta=colsum(C0)'/n;ns=1
        }
        hist=theta';trace=(0,.,initmissing,ns==1)
        if(ns==1) for(k=1;k<=maxiter;k++) {
            attempts=k;Cnew=J(n,n,.);row0=rows(D)
            for(d=1;d<=n;d++) for(j=1;j<=n;j++) {
                rc=_g114_dispatch(X,Y,XR,YR,sx,sy,j,3,.,d,theta[d],tol,v,u,value,k,obs_indices,D,counts,source)
                if(!rc) Cnew[d,j]=value
            }
            failed=sum(Cnew:>=.)
            if(failed) {
                ns=3;trace=trace\(k,.,failed,0);break
            }
            if(rows(D)>row0) D[(row0+1)..rows(D),32]=J(rows(D)-row0,1,1)
            next=colsum(Cnew)'/n;delta=max(abs(next-theta))
            iterations=k;Cn=Cnew;theta=next;hist=hist\theta'
            trace=trace\(k,delta,0,1)
            if(show_debug) printf("  Nash iteration %g: max change = %12.8g\n",k,delta)
            if(delta<tol) {
                ns=0;break
            }
        }
    }
    printf("\nCross-Efficiency Game Results (CRS, input orientation)\n")
    printf("%20s | %10s", "DMU", "Self")
    if(wanta) printf(" %10s","Aggressive")
    if(wantb) printf(" %10s","Benevolent")
    if(wantn) printf(" %10s","Nash state")
    printf("\n")
    for(i=1;i<=n;i++) {
        printf("%20s | %10.5f",labels[i],self[i])
        if(wanta) {
            printf(" %10.5f",a[i]);st_store(obs_indices[i],gen_a,a[i])
        }
        if(wantb) {
            printf(" %10.5f",b[i]);st_store(obs_indices[i],gen_b,b[i])
        }
        if(wantn) {
            printf(" %10.5f",theta[i]);st_store(obs_indices[i],gen_n,theta[i])
        }
        printf("\n")
    }
    if(sum(self:>=.)) printf("Note: %g primary weight solve(s) unavailable.\n",sum(self:>=.))
    if(wanta & sum(sa:==3)) printf("Note: %g aggressive secondary solve(s) unavailable; primary weights retained.\n",sum(sa:==3))
    if(wantb & sum(sb:==3)) printf("Note: %g benevolent secondary solve(s) unavailable; primary weights retained.\n",sum(sb:==3))
    if(wanta & sum(pa)<n*(n-1)) printf("Note: aggressive means use available peers; %g pairs unavailable.\n",n*(n-1)-sum(pa))
    if(wantb & sum(pb)<n*(n-1)) printf("Note: benevolent means use available peers; %g pairs unavailable.\n",n*(n-1)-sum(pb))
    if(wantn) {
        printf("Nash status %g: %g complete iteration(s), %g attempted.\n",ns,iterations,attempts)
        if(ns==0) printf("Successive complete states satisfy tolerance; this is a numerical stopping test.\n")
        if(ns==1) printf("Iteration limit reached; returning the last complete state, not confirmed convergence.\n")
        if(ns==2) printf("Initialization incomplete (%g missing cells); Nash scores unavailable.\n",initmissing)
        if(ns==3) printf("Iteration failed (%g pair solves); returning the last complete state, not an equilibrium.\n",failed)
    }
    if(show_detail) {
        printf("\nWeight stages: 0=primary unavailable, 2=secondary accepted, 3=primary fallback\n")
        for(i=1;i<=n;i++) {
            printf("%20s",labels[i])
            if(wanta) printf(" | aggressive stage %g, peers %g",sa[i],pa[i])
            if(wantb) printf(" | benevolent stage %g, peers %g",sb[i],pb[i])
            printf("\n")
        }
        if(wantn) {
            printf("Nash trace: iteration, max change, missing pairs, complete state\n")
            for(i=1;i<=rows(trace);i++) printf("%6.0f %14.8g %8.0f %3.0f\n",trace[i,1],trace[i,2],trace[i,3],trace[i,4])
        }
    }
    meta=(sum(self:>=.),(wanta ? sum(sa:==3) : .),(wantb ? sum(sb:==3) : .),
    ns,iterations,attempts,failed,delta,initmissing,
    (wanta ? n*(n-1)-sum(pa) : .),(wantb ? n*(n-1)-sum(pb) : .))
    st_matrix(prefix+"_meta",meta)
    _g114_report(D,counts,prefix,show_detail)
    if(show_matrix) {
        st_matrix(prefix+"_self",self);st_matrix(prefix+"_v",V:/sx);st_matrix(prefix+"_u",U:/sy)
        st_matrix(prefix+"_initial",C0)
        if(wanta) {
            _dea_cg_show(Ca,labels,"Aggressive cross-efficiency matrix")
            st_matrix(prefix+"_agg",Ca);st_matrix(prefix+"_va",Va:/sx);st_matrix(prefix+"_ua",Ua:/sy)
            st_matrix(prefix+"_sourcea",sourcea);st_matrix(prefix+"_sa",sa);st_matrix(prefix+"_pa",pa)
        }
        if(wantb) {
            _dea_cg_show(Cb,labels,"Benevolent cross-efficiency matrix")
            st_matrix(prefix+"_ben",Cb);st_matrix(prefix+"_vb",Vb:/sx);st_matrix(prefix+"_ub",Ub:/sy)
            st_matrix(prefix+"_sourceb",sourceb);st_matrix(prefix+"_sb",sb);st_matrix(prefix+"_pb",pb)
        }
        if(wantn) {
            _dea_cg_show(Cn,labels,"Nash last complete matrix (protected peer d, target j)")
            st_matrix(prefix+"_nash",Cn);st_matrix(prefix+"_history",hist);st_matrix(prefix+"_trace",trace)
        }
    }
}

end
