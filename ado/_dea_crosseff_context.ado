*! version 1.0.0 06nov2026
*! Shared solver context preservation for cross-efficiency.
version 16.0
capture mata: mata drop _c57_context_ready()
capture mata: mata drop _c57_context_names()
capture mata: mata drop _c57_context_save()
capture mata: mata drop _c57_context_defaults()
capture mata: mata drop _c57_context_restore()
mata:
mata set matastrict on
real scalar _c57_context_ready() return(57)
string rowvector _c57_context_names()
{
 return(tokens("_dea_rec_on _dea_rec_limit _dea_rec_scope_id _dea_rec_rows _dea_rec_counts _dea_rec_ctx _dea_link_mode _dea_link_record _dea_link_fast _dea_link_records _rsm_recording _rsm_records _rsm_counts _rsm_success_policy _rsm_success_recording _rsm_success_records _rsm_success_counts _rsm_success_fast _rsm_success_fastcounts"))
}
transmorphic scalar _c57_context_save()
{
 transmorphic scalar T
 pointer scalar p
 string rowvector names
 real scalar i
 T=asarray_create("string",1);names=_c57_context_names()
 for(i=1;i<=cols(names);i++) {
  p=findexternal(names[i]);asarray(T,names[i]+":existed",p!=NULL)
  if(p!=NULL) asarray(T,names[i],*p)
 }
 return(T)
}
// After loading and quiet reset, retain initialized defaults for globals that
// did not exist on entry. Removing them would leave loaded solvers uninitialized.
void _c57_context_defaults(transmorphic scalar T)
{
 pointer scalar p
 string rowvector names
 real scalar i
 names=_c57_context_names()
 for(i=1;i<=cols(names);i++) if(!asarray(T,names[i]+":existed")) {
  p=findexternal(names[i]);if(p!=NULL) asarray(T,names[i],*p)
 }
}
void _c57_context_restore(transmorphic scalar T)
{
 pointer scalar p
 string rowvector names
 real scalar i
 names=_c57_context_names()
 for(i=1;i<=cols(names);i++) if(asarray_contains(T,names[i])) {
  p=findexternal(names[i]);if(p==NULL) p=crexternal(names[i])
  *p=asarray(T,names[i])
 }
}

end
