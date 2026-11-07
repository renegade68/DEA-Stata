*! version 1.0.0 06nov2026
*! Shared Cooper objective-weight representability; no LP engine dependency.
version 16.0
capture mata: mata drop _dea_cooper_weights_safe()
capture mata: mata drop _dea_cooper_weights_ready()
mata:
mata set matastrict on
real scalar _dea_cooper_weights_safe(real matrix X,real matrix Y,string scalar weight)
{
 real rowvector sx,sy,w
 real matrix xs,ys
 real scalar i
 sx=colmax(X);sy=colmax(Y);xs=X:/sx;ys=Y:/sy
 if(weight=="raw") {
  w=(sx,sy);w=w:/max(w)
  return(!hasmissing(w) & min(w)>0)
 }
 for(i=1;i<=rows(X);i++) {
  w=1:/(xs[i,.],ys[i,.])
  if(hasmissing(w) | min(w)<=0) return(0)
  w=w:/max(w);if(hasmissing(w) | min(w)<=0) return(0)
 }
 return(1)
}

real scalar _dea_cooper_weights_ready() return(5401)
end
