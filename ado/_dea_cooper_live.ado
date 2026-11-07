*! version 1.0.0 06nov2026
*! Internal file-input Cooper pipeline interface; public commands use data arrays.
capture program drop _dea_cooper_live
program define _dea_cooper_live, rclass
 version 16.0
 syntax using/, N(integer) Inputs(integer) Outputs(integer) Index(integer) [POLicy(string) REFerence1(string) REFerence2(string) XSCALE(string) YSCALE(string) RANGes SELECTLimit(integer 16) ENGine(string) INputnames(string) OUTputnames(string) DMUnames(string) AUDit(string)]
 if "`policy'"=="" local policy "stage1"
 local policy=lower(strtrim("`policy'"))
 if !inlist("`policy'","stage1","stage2") exit 198
 if "`engine'"=="" local engine "simplex"
 local engine=lower(strtrim("`engine'"))
 if "`engine'"!="simplex" {
  di as error "this experimental pipeline supports engine(simplex) only"
  exit 198
 }
 if !inlist(`selectlimit',16,128) | `n'<1 | `inputs'<1 | `outputs'<1 | 2*`n'+`inputs'+`outputs'+3>64 | `index'<1 | `index'>`n' exit 198
 foreach pair in xscale yscale {
  local dim=cond("`pair'"=="xscale",`inputs',`outputs')
  if "``pair''"=="" {
   forvalues j=1/`dim' {
    local `pair' "``pair'' 1"
   }
  }
  local k: word count ``pair''
  if `k'!=`dim' exit 198
  foreach x of local `pair' {
   capture confirm number `x'
   if _rc exit 198
   if `x'<=0 | missing(`x') exit 198
  }
 }
 foreach kind in input output dmu {
  local dim=cond("`kind'"=="input",`inputs',cond("`kind'"=="output",`outputs',`n'))
  if "``kind'names'"=="" {
   forvalues j=1/`dim' {
    local `kind'names "``kind'names' `kind'`j'"
   }
  }
  local k: word count ``kind'names'
  local uniq: list uniq `kind'names
  local ku: word count `uniq'
  if `k'!=`dim' | `ku'!=`dim' exit 198
  foreach name of local `kind'names {
   capture confirm name `name'
   if _rc exit 198
  }
 }
 local bothnames "`inputnames' `outputnames'"
 local uniq: list uniq bothnames
 local ku: word count `uniq'
 if `ku'!=`inputs'+`outputs' exit 198
 foreach f in c E b U ub lo hi {
  confirm file "`using'/radial/`f'.bin"
 }
 confirm file "`using'/slack/c.bin"
 confirm file "`using'/decomposition/c.bin"
 foreach ref in reference1 reference2 {
  if "``ref''"!="" confirm file "``ref''"
 }
 capture mata: assert(_c52_ready()==52)
 if _rc {
  quietly findfile _dea_cooper_live_core.ado
  quietly run "`r(fn)'"
 }
 capture mata: assert(_c52_ready()==52)
 if _rc exit 499
 tempname l1 l2 sc it tt co pmin pmax cmin cmax info state
 local _c52_names "`l1' `l2' `sc' `it' `tt' `co' `pmin' `pmax' `cmin' `cmax' `info' `state'"
 mata: _c52_bridge(st_local("using"),strtoreal(st_local("n")),strtoreal(st_local("inputs")),strtoreal(st_local("outputs")),strtoreal(st_local("index")),st_local("policy"),st_local("reference1"),st_local("reference2"),strtoreal(tokens(st_local("xscale"))),strtoreal(tokens(st_local("yscale"))),st_local("ranges")!="",strtoreal(st_local("selectlimit")),tokens(st_local("_c52_names")))
 matrix colnames `l1'=`dmunames'
 matrix colnames `l2'=`dmunames'
 matrix colnames `sc'=congestion phi te total technical
 foreach mat in it tt co cmin cmax {
  matrix colnames ``mat''=`inputnames'
 }
 matrix colnames `pmin'=`bothnames'
 matrix colnames `pmax'=`bothnames'
 matrix colnames `state'=radial slack anchor decomposition projection congestion
 return matrix lambda_stage1=`l1'
 return matrix lambda_stage2=`l2'
 return matrix scores=`sc'
 return matrix input_total=`it'
 return matrix input_technical=`tt'
 return matrix input_congestion=`co'
 return matrix projection_min=`pmin'
 return matrix projection_max=`pmax'
 return matrix congestion_min=`cmin'
 return matrix congestion_max=`cmax'
 return matrix pipeline_state=`state'
 local fields "usable selection_status selection_applied conversion_status rounded_values base_certificate selection_certificate selection_maxpivots reference_distance ranges_complete range_conversion_status range_rounded_values"
 local j=0
 foreach field of local fields {
  local ++j
  return scalar `field'=`info'[1,`j']
 }
 return scalar ranges_requested=("`ranges'"!="")
 return local requested_policy "`policy'"
 return local result_source "`_c52_source'"
 return local selection_reason "`_c52_reason'"
 return local range_reason "`_c52_range_reason'"
 return local actual_engine "simplex"
 return local actual_backend "rsm_proposal_and_bounded_integer_basis_repair"
 return local certificate_engine "bounded_integer_primal_dual"
 return local component_units "lp_units_times_xscale"
 return local range_units "lp_units_times_xscale_yscale"
 return local certificate_scope "internal_exact_point"
 return local range_scope "projection_on_slack_face_congestion_on_original_technical_face"
 return local integration_status "isolated_prototype"
 return local cmd "_dea_cooper_live"
end
