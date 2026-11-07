*! version 1.0.0  02nov2026
*! Diagnose returned TE stage-2 projections without changing their acceptance.
version 16.0
capture mata: mata drop _dea_teaudit_ready()
capture mata: mata drop _dea_teaudit()
capture mata: mata drop _dea_teaudit_store()
mata:
mata set matastrict on
real scalar _dea_teaudit_ready() return(1)

// X/Y are the normalized solver data; XR/YR are the original source data.
// si/so are the returned slacks in normalized units. These diagnostics inspect
// the returned vector (including any existing extraction cleanup), not an
// exact feasible/optimal certificate and not a new solver acceptance gate.
real rowvector _dea_teaudit(real matrix X, real matrix Y,
    real matrix XR, real matrix YR, real scalar k, string scalar ort,
    string scalar rts, real rowvector xs, real rowvector ys,
    real rowvector lam, real rowvector si, real rowvector so,
    real scalar radial0, real scalar radial2, real scalar status, real scalar tol)
{
    real rowvector out, dn, dr, z, zr, a, ar, bs, v
    real matrix empty
    real colvector rhs0
    real scalar n, m, s, j, side, rhs, rhsr, factor, total, err, den
    out=J(1,23,.)
    out[1..5]=(radial0,radial2,radial2-radial0,status,0)
    n=rows(X);m=cols(X);s=cols(Y)
    if (status!=0 | radial2>=. | hasmissing(lam) | hasmissing(si) |
        hasmissing(so) | cols(lam)!=n | cols(si)!=m | cols(so)!=s) return(out)
    out[7..23]=J(1,17,0)
    empty=J(0,n+1,0);rhs0=J(0,1,0)
    for (side=1;side<=2;side++) {
        factor=(side==1 ? (ort=="input" ? radial2 : 1) :
                         (ort=="output" ? radial2 : 1))
        for (j=1;j<=(side==1 ? m : s);j++) {
            if (side==1) {
                a=(X[.,j]',1);ar=(XR[.,j]',1)
                rhs=factor*X[k,j];rhsr=factor*XR[k,j]
                z=(lam,si[j]);zr=(lam,si[j]*xs[j])
            }
            else {
                a=(Y[.,j]',-1);ar=(YR[.,j]',-1)
                rhs=factor*Y[k,j];rhsr=factor*YR[k,j]
                z=(lam,so[j]);zr=(lam,so[j]*ys[j])
            }
            dn=_dea_lpcheck(z,a,rhs,empty,rhs0,J(1,n+1,.),J(1,n+1,.))
            dr=_dea_lpcheck(zr,ar,rhsr,empty,rhs0,J(1,n+1,.),J(1,n+1,.))
            // Do not report an arithmetic overflow as a zero residual.
            if (hasmissing(dn) | hasmissing(dr)) return((out[1..5],J(1,18,.)))
            out[6+side]=max((out[6+side],dr[1]))
            out[14+side]=max((out[14+side],dn[4]))
            out[21+side]=max((out[21+side],dn[1]))
        }
    }
    total=sum(lam)
    if (total>=.) return((out[1..5],J(1,18,.)))
    out[10]=total;out[11]=total-1
    if (rts!="crs") {
        err=(rts=="vrs" ? abs(total-1) :
            (rts=="drs" ? max((0,total-1)) : max((0,1-total))))
        den=max((1,sum(colmax((J(1,n,1)\abs(lam))))))
        out[9]=err;out[17]=err/den
    }
    // Bound diagnostics are separate because raw slack units differ by column.
    out[12]=max((0,-lam))
    out[13]=max((0,-si:*xs));out[14]=max((0,-so:*ys))
    for (side=1;side<=3;side++) {
        z=(side==1 ? lam : (side==2 ? si : so))
        v=colmax((J(1,cols(z),0)\(-z)))
        bs=v:/colmax((J(1,cols(z),1)\abs(z)))
        out[17+side]=max(bs)
    }
    out[21]=max(out[15..20])
    if (hasmissing(out[7..23])) return((out[1..5],J(1,18,.)))
    out[5]=1;out[6]=(out[21]<=tol)
    return(out)
}

void _dea_teaudit_store(string scalar name, real matrix D)
{
    if (name=="") return
    st_matrix(name,D)
    st_matrixcolstripe(name,(J(23,1,""),
        ("radial_initial","radial_stage2","radial_adjustment","stage2_status",
         "checked","within_tolerance","input_abs","output_abs","rts_abs",
         "lambda_sum","lambda_gap","lambda_bound","slack_in_bound","slack_out_bound",
         "input_scaled","output_scaled","rts_scaled","lambda_scaled",
         "slack_in_scaled","slack_out_scaled","max_scaled",
         "input_normalized","output_normalized")'))
}
end
