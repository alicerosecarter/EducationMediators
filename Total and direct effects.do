
global PCs PC01 PC02 PC03 PC04 PC05 PC06 PC07 PC08 PC09 PC10
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth

drop if CVD_inc==.

********************************************************************************

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {

	putexcel set RESULTS_FILE, sheet(`out') modify
	putexcel A1="Method" B1="Scale" C1="Mediator" D1="Total Effect" E1="Total - LCI" F1="Total - UCI" G1="Direct Effect" H1="Direct - LCI" I1="Direct - UCI" J1="" ///
		K1="Exposure - Mediator" L1="Exp-Med - LCI" M1="Exp-Med - UCI" N1="Mediator - Outcome" O1="Med-Out - LCI" P1="Med-Out - UCI" 
}

gen combined = 1
lab var combined "Combined risk factors"

******************** Total Effects - Observational on RD scale *****************

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=1	
	
	foreach med in zbmi zsbp zcsi combined {
		
		local x=`x'+1
	
	putexcel set RESULTS_FILE, sheet(`out') modify

		regress `out' eduyears_scaled $max_covar
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="Observational" B`x'="Risk Difference" C`x'="`med_label'" D`x'=_b[eduyears_scaled] E`x'=`ll_2sls1' F`x'=`ul_2sls1'
	
	}	
}

******************* Total Effects - Observational on logOR scale ***************

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=5
	
	foreach med in zbmi zsbp zcsi combined {

		local x=`x'+1
	
	putexcel set RESULTS_FILE, sheet(`out') modify

		logit `out' eduyears_scaled $max_covar
	
	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="Observational" B`x'="Log OR" C`x'="`med_label'" D`x'=_b[eduyears_scaled] E`x'=`ll_2sls1' F`x'=`ul_2sls1'
	
	}	
}


foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=5
	
	foreach med in zbmi zsbp zcsi combined {

		local x=`x'+1
	
	putexcel set RESULTS_FILE, sheet(`out') modify

		logit `out' `med' $max_covar
	
	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[`med']-1.96*_se[`med']
	local ul_2sls1 = _b[`med']+1.96*_se[`med']
	
	putexcel AA`x'=_b[`med'] AB`x'=`ll_2sls1' AC`x'=`ul_2sls1'
	
	}	
}


******************** Total Effects - One-sample MR on RD scale *****************

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=11
		
	foreach med in zbmi zsbp zcsi combined {

		local x=`x'+1
		
	putexcel set RESULTS_FILE, sheet(`out') modify

		ivreg2 `out' (eduyears_scaled = ea_weighted) $PCs $max_covar
	
	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="MR - One-sample" B`x'="Risk Difference" C`x'="`med_label'" D`x'=_b[eduyears_scaled] E`x'=`ll_2sls1' F`x'=`ul_2sls1'
	
	}	
}


******************* Total Effects - One-sample MR on logOR scale ***************
capture drop education_MR
regress eduyears_scaled ea_weighted $max_covar $PCs
predict education_MR, xb

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=15
	
	foreach med in zbmi zsbp zcsi combined {

		local x=`x'+1
	
	putexcel set RESULTS_FILE, sheet(`out') modify

		logit `out' education_MR $PCs $max_covar
	
	
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[education_MR]-1.96*_se[education_MR]
	local ul_2sls1 = _b[education_MR]+1.96*_se[education_MR]
	
	putexcel A`x'="MR - One-sample" B`x'="Log OR" C`x'="`med_label'" D`x'=_b[education_MR] E`x'=`ll_2sls1' F`x'=`ul_2sls1'
	
	}	
}

************* Direct & Indirect Effects - Observational on RD scale ************

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=1	
	
	foreach med in zbmi zsbp zcsi {
		
		local x=`x'+1
	
	putexcel set RESULTS_FILE, sheet(`out') modify

	regress `out' eduyears_scaled `med' $max_covar
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	local ll_med = _b[`med']-1.96*_se[`med']
	local ul_med = _b[`med']+1.96*_se[`med']
	
	putexcel G`x'=_b[eduyears_scaled] H`x'=`ll_2sls1' I`x'= `ul_2sls1'
	putexcel N`x'=_b[`med'] O`x'=`ll_med' P`x'= `ul_med'

	}	
}

******* Direct Effects for Multiple Mediators - Observational on RD scale ******


foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=5	
	
	putexcel set RESULTS_FILE, sheet(`out') modify

	regress `out' eduyears_scaled zbmi zsbp zcsi $max_covar
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel G`x'=_b[eduyears_scaled] H`x'=`ll_2sls1' I`x'= `ul_2sls1'

	}	

**************** Direct & Indirect Effects - Observational on logOR scale ******

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=5
	
	foreach med in zbmi zsbp zcsi {

		local x=`x'+1
	
	putexcel set RESULTS_FILE, sheet(`out') modify

	logit `out' eduyears_scaled `med' $max_covar
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	local ll_med = _b[`med']-1.96*_se[`med']
	local ul_med = _b[`med']+1.96*_se[`med']
	
	putexcel G`x'=_b[eduyears_scaled] H`x'=`ll_2sls1' I`x'= `ul_2sls1'
	putexcel N`x'=_b[`med'] O`x'=`ll_med' P`x'= `ul_med'
	
	}	
}

******* Direct Effects for Multiple Mediators - Observational on RD scale ******


foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=9	
	
	putexcel set RESULTS_FILE, sheet(`out') modify

	logit `out' eduyears_scaled zbmi zsbp zcsi $max_covar
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel G`x'=_b[eduyears_scaled] H`x'=`ll_2sls1' I`x'= `ul_2sls1'

	}	

*************** Direct & Indirect Effects - One-sample MR on RD scale **********
							
										*BMI*
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=11

		local x=`x'+1
		
	putexcel set RESULTS_FILE, sheet(`out') modify

	ivreg2 `out' (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
	
	
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	local ll_med = _b[zbmi]-1.96*_se[zbmi]
	local ul_med = _b[zbmi]+1.96*_se[zbmi]
	
	putexcel G`x'=_b[eduyears_scaled] H`x'=`ll_2sls1' I`x'= `ul_2sls1'
	putexcel N`x'=_b[zbmi] O`x'=`ll_med' P`x'= `ul_med'

	}
								*Composite smoking*
			
								
local x=1
local y=5
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {


local x=`x'+1
local y=`y'+1
	putexcel set  zcsi_results, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or_education" D1="lci_education" E1="uci_education" F1="n" G1="beta_zcsi" H1="LCI_zcsi" I1="UCI_zcsi"

	ivreg2 `out' (eduyears_scaled zcsi = ea_weighted csi_score_1_weighted) $PCs $max_covar if sample==2
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	local ll_med = _b[zcsi]-1.96*_se[zcsi]
	local ul_med = _b[zcsi]+1.96*_se[zcsi]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[eduyears_scaled] D`x'=`ll_2sls1' E`x'=`ul_2sls1' F`x'=`e(N)' ///
		G`x'=_b[zcsi] H`x'=`ll_med' I`x'=`ul_med' 
	
	ivreg2 `out' (eduyears_scaled zcsi = ea_weighted csi_score_2_weighted) $PCs $max_covar if sample==1
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	local ll_med = _b[zcsi]-1.96*_se[zcsi]
	local ul_med = _b[zcsi]+1.96*_se[zcsi]
	
	putexcel A`y'="2" B`y'="`out_label'" C`y'=_b[eduyears_scaled] D`y'=`ll_2sls1' E`y'=`ul_2sls1' F`y'=`e(N)' ///	
		G`y'=_b[zcsi] H`y'=`ll_med' I`y'=`ul_med' 
}

									*SBP*
local x=1
local y=5
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {


local x=`x'+1
local y=`y'+1
	putexcel set  SBP_results, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or_education" D1="lci_education" E1="uci_education" F1="n" G1="beta_bp" H1="LCI_bp" I1="UCI_bp"

	ivreg2 `out' (eduyears_scaled zsbp = ea_weighted sbp_score_1_weighted) $PCs $max_covar if sample==2
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	local ll_med = _b[zsbp]-1.96*_se[zsbp]
	local ul_med = _b[zsbp]+1.96*_se[zsbp]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[eduyears_scaled] D`x'=`ll_2sls1' E`x'=`ul_2sls1' F`x'=`e(N)' ///
		G`x'=_b[zsbp] H`x'=`ll_med' I`x'=`ul_med' 
	
	ivreg2 `out' (eduyears_scaled zsbp = ea_weighted sbp_score_2_weighted) $PCs $max_covar if sample==1
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	local ll_med = _b[zsbp]-1.96*_se[zsbp]
	local ul_med = _b[zsbp]+1.96*_se[zsbp]
	
	putexcel A`y'="2" B`y'="`out_label'" C`y'=_b[eduyears_scaled] D`y'=`ll_2sls1' E`y'=`ul_2sls1' F`y'=`e(N)' ///	
		G`y'=_b[zsbp] H`y'=`ll_med' I`y'=`ul_med' 
}

****************************MULTIPLE MEDIATORS**********************************
							
local x=1
local y=5
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
local x=`x'+1
local y=`y'+1
		
	putexcel set SPLIT_SAMPLE_MULTIPLEMED, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or_education_controlled_for_med" D1="lci_education_controlled_for_med" E1="uci_education_controlled_for_med" F1="n" 

	ivreg2 `out' (eduyears_scaled zbmi zcsi zsbp = ea_weighted zbmi_weighted csi_score_1_weighted sbp_score_1_weighted) $PCs $max_covar if sample==2
	
	
	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[eduyears_scaled] D`x'=`ll_2sls1' E`x'=`ul_2sls1' F`x'=`e(N)' 
		
ivreg2 `out' (eduyears_scaled zbmi zcsi zsbp = ea_weighted zbmi_weighted csi_score_2_weighted sbp_score_2_weighted) $PCs $max_covar if sample==1

	local out_label : var label `out'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel C`y'=_b[eduyears_scaled] D`y'=`ll_2sls1' E`y'=`ul_2sls1' F`y'=`e(N)' 

	}
	
*********** Direct & Indirect Effects - One-sample MR on log OR scale **********

							*BMI*

capture drop dir_eabmi
capture drop dir_bmiea
regress eduyears_scaled ea_weighted zbmi_weighted $max_covar $PCs
predict dir_eabmi, xb
regress zbmi zbmi_weighted ea_weighted $max_covar $PCs
predict dir_bmiea, xb


foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {

	putexcel set RESULTS_FILE, sheet(`out') modify
	
	logit `out' dir_eabmi dir_bmiea $max_covar $PCs, vce(robust) // Natural Direct effect of Education
	
	local out_label : var label `out'
	local ll_2sls1 = _b[dir_eabmi]-1.96*_se[dir_eabmi]
	local ul_2sls1 = _b[dir_eabmi]+1.96*_se[dir_eabmi]
	local ll_med = _b[dir_bmiea]-1.96*_se[dir_bmiea]
	local ul_med = _b[dir_bmiea]+1.96*_se[dir_bmiea]
	
		
	putexcel G16=_b[dir_eabmi] H16=`ll_2sls1' I16= `ul_2sls1'
	putexcel N16=_b[dir_bmiea] O16=`ll_med' P16= `ul_med'
	
	
}	

							*zcsi*
							
capture drop dir_eazcsi_1
capture drop dir_zcsiea_1
regress eduyears_scaled ea_weighted csi_score_2_weighted $max_covar $PCs if sample==1
predict dir_eazcsi_1, xb
regress zcsi  csi_score_2_weighted ea_weighted $max_covar $PCs if sample==1
predict dir_zcsiea_1, xb

capture drop dir_eazcsi_2
capture drop dir_zcsiea_2
regress eduyears_scaled ea_weighted csi_score_1_weighted $max_covar $PCs if sample==2
predict dir_eazcsi_2, xb
regress zcsi  csi_score_1_weighted ea_weighted $max_covar $PCs if sample==2
predict dir_zcsiea_2, xb


local x=1
local y =5
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
	putexcel set zcsi_results_logOR, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or_education" D1="lci_education" E1="uci_education" F1="n" G1="beta_zcsi" H1="LCI_zcsi" I1="UCI_zcsi"
	
	local x=`x'+1
	local y=`y'+1
	
	logit `out' dir_eazcsi_1 dir_zcsiea_1 $max_covar $PCs if sample==1, vce(robust) // Natural Direct effect of Education

	
	local out_label : var label `out'
	local ll_2sls1 = _b[dir_eazcsi_1]-1.96*_se[dir_eazcsi_1]
	local ul_2sls1 = _b[dir_eazcsi_1]+1.96*_se[dir_eazcsi_1]
	local ll_med = _b[dir_zcsiea_1]-1.96*_se[dir_zcsiea_1]
	local ul_med = _b[dir_zcsiea_1]+1.96*_se[dir_zcsiea_1]
	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[dir_eazcsi_1] D`x'=`ll_2sls1' E`x'=`ul_2sls1'  ///
		G`x'=_b[dir_zcsiea_1] H`x'=`ll_med' I`x'=`ul_med' 
	
	logit `out' dir_eazcsi_2 dir_zcsiea_2 $max_covar $PCs if sample==2, vce(robust) // Effect of zcsi on outcomes
	
	local out_label : var label `out'
	local ll_2sls1 = _b[dir_eazcsi_2]-1.96*_se[dir_eazcsi_2]
	local ul_2sls1 = _b[dir_eazcsi_2]+1.96*_se[dir_eazcsi_2]
	local ll_med = _b[dir_zcsiea_2]-1.96*_se[dir_zcsiea_2]
	local ul_med = _b[dir_zcsiea_2]+1.96*_se[dir_zcsiea_2]
	
	putexcel A`y'="2" B`y'="`out_label'" C`y'=_b[dir_eazcsi_2] D`y'=`ll_2sls1' E`y'=`ul_2sls1'  ///	
		G`y'=_b[dir_zcsiea_2] H`y'=`ll_med' I`y'=`ul_med' 
		
	
}	

							*SBP*

capture drop dir_easbp_1
capture drop dir_sbpea_1
regress eduyears_scaled ea_weighted zsbp_score_2_weighted $max_covar $PCs if sample==1
predict dir_easbp_1, xb
regress zsbp  zsbp_score_2_weighted ea_weighted $max_covar $PCs if sample==1
predict dir_sbpea_1, xb

capture drop dir_easbp_2
capture drop dir_sbpea_2
regress eduyears_scaled ea_weighted zsbp_score_1_weighted $max_covar $PCs if sample==2
predict dir_easbp_2, xb
regress zsbp  zsbp_score_1_weighted ea_weighted $max_covar $PCs if sample==2
predict dir_sbpea_2, xb

local x=1
local y =5
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
	putexcel set SBP_results_logOR, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or_education" D1="lci_education" E1="uci_education" F1="n" G1="beta_sbp" H1="LCI_sbp" I1="UCI_sbp"
	
	local x=`x'+1
		
	logit `out' dir_easbp_1 dir_sbpea_1 $max_covar $PCs if sample==1, vce(robust) 
	local lci = _b[dir_easbp_1]-1.96*_se[dir_easbp_1]
	local uci = _b[dir_easbp_1]+1.96*_se[dir_easbp_1]
	local ll_med = _b[dir_sbpea_1]-1.96*_se[dir_sbpea_1]
	local ul_med = _b[dir_sbpea_1]+1.96*_se[dir_sbpea_1]
	local out_label : var label `out'
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[dir_easbp_1] D`x'=`lci' E`x'=`uci' F`x'=`e(N)' ///
		G`x'=_b[dir_sbpea_1] H`x'=`ll_med' I`x'=`ul_med' 

	local y = `y'+1
	
	logit `out' dir_easbp_2 dir_sbpea_2 $max_covar $PCs if sample==2, vce(robust) 
	local lci = _b[dir_easbp_2]-1.96*_se[dir_easbp_2]
	local uci = _b[dir_easbp_2]+1.96*_se[dir_easbp_2]
	local ll_med = _b[dir_sbpea_2]-1.96*_se[dir_sbpea_2]
	local ul_med = _b[dir_sbpea_2]+1.96*_se[dir_sbpea_2]
	local out_label : var label `out'
	putexcel A`y'="2" B`y'="`out_label'" C`y'=_b[dir_easbp_2] D`y'=`lci' E`y'=`uci' F`y'=`e(N)' ///
		G`y'=_b[dir_sbpea_2] H`y'=`ll_med' I`y'=`ul_med' 
}	
									*Multiple Mediators*
							
capture drop ea_1 bmi_1 csi_1 sbp_1
regress eduyears_scaled ea_weighted zbmi_weighted zsbp_score_2_weighted csi_score_2_weighted $max_covar $PCs if sample==1
predict ea_1, xb
regress zbmi ea_weighted zbmi_weighted zsbp_score_2_weighted csi_score_2_weighted $max_covar $PCs if sample==1
predict bmi_1, xb
regress zcsi ea_weighted zbmi_weighted zsbp_score_2_weighted csi_score_2_weighted $max_covar $PCs if sample==1
predict csi_1, xb
regress zsbp ea_weighted zbmi_weighted zsbp_score_2_weighted csi_score_2_weighted $max_covar $PCs if sample==1
predict sbp_1, xb

capture drop ea_2 bmi_2 csi_2 sbp_2
regress eduyears_scaled ea_weighted zbmi_weighted zsbp_score_1_weighted csi_score_1_weighted $max_covar $PCs if sample==2
predict ea_2, xb
regress zbmi ea_weighted zbmi_weighted zsbp_score_1_weighted csi_score_1_weighted $max_covar $PCs if sample==2
predict bmi_2, xb
regress zcsi ea_weighted zbmi_weighted zsbp_score_1_weighted csi_score_1_weighted $max_covar $PCs if sample==2
predict csi_2, xb
regress zsbp ea_weighted zbmi_weighted zsbp_score_1_weighted csi_score_1_weighted $max_covar $PCs if sample==2
predict sbp_2, xb


local x=1
local y =5
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
	putexcel set MultipledMed_logOR, sheet(all) modify
	putexcel A1="sample" B1="out" C1="or_education" D1="lci_education" E1="uci_education" F1="n" G1="beta_zcsi" H1="LCI_zcsi" I1="UCI_zcsi"
	
	local x=`x'+1
	local y=`y'+1
	
	logit `out' ea_1 bmi_1 csi_1 sbp_1 $max_covar $PCs if sample==1, vce(robust)

	
	local out_label : var label `out'
	local ll_2sls1 = _b[ea_1]-1.96*_se[ea_1]
	local ul_2sls1 = _b[ea_1]+1.96*_se[ea_1]

	
	putexcel A`x'="1" B`x'="`out_label'" C`x'=_b[ea_1] D`x'=`ll_2sls1' E`x'=`ul_2sls1'  
	
	logit `out' ea_2 bmi_2 csi_2 sbp_2 $max_covar $PCs if sample==2, vce(robust)
	
	local out_label : var label `out'
	local ll_2sls1 = _b[ea_2]-1.96*_se[ea_2]
	local ul_2sls1 = _b[ea_2]+1.96*_se[ea_2]

	
	putexcel A`y'="2" B`y'="`out_label'" C`y'=_b[ea_2] D`y'=`ll_2sls1' E`y'=`ul_2sls1' 
		
	
}
	
********* Total Effects Exposure-Mediator - Observational on RD scale **********

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=1	
	
	foreach med in zbmi zsbp zcsi {
		
		local x=`x'+1
	
	putexcel set RESULTS_FILE, sheet(`out') modify

		regress `med' eduyears_scaled $max_covar
	
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel K`x'=_b[eduyears_scaled] L`x'=`ll_2sls1' M`x'=`ul_2sls1'
	
	}	
}

******** Total Effects Exposure-Mediator - Observational on log OR scale *******

/* All mediators are continuous variables and as such the exposure-mediator
	pathway cannot be estimated on the log OR scale */

********** Total Effects Exposure-Mediator - One-sample MR on RD scale *********

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
	
		local x=11
		
	foreach med in zbmi zsbp zcsi {

		local x=`x'+1
		
	putexcel set RESULTS_FILE, sheet(`out') modify

		ivreg2 `med' (eduyears_scaled = ea_weighted) $PCs $max_covar
	
		
	local out_label : var label `out'
	local med_label : var label `med'
	local ll_2sls1 = _b[eduyears_scaled]-1.96*_se[eduyears_scaled]
	local ul_2sls1 = _b[eduyears_scaled]+1.96*_se[eduyears_scaled]
	
	putexcel K`x'=_b[eduyears_scaled] L`x'=`ll_2sls1' M`x'=`ul_2sls1'
	
	}	
}

******** Total Effects Exposure-Mediator - One-sample MR on log OR scale *******

/* All mediators are continuous variables and as such the exposure-mediator
	pathway cannot be estimated on the log OR scale */

*******************************************************************************
*Meta-analysing split samples

********************************************************************************

*zcsi on the RD scale

import excel "zcsi_results.xlsx", sheet("all") firstrow clear

destring sample, replace
lab var sample "Sample"
lab def sample 1 "Sample 1" 2 "Sample 2"
lab val sample sample

lab var n "N"

generate outcome = .
replace outcome = 1 if out=="Incident CVD Case"
replace outcome = 2 if out=="Incident Stroke Case"
replace outcome = 3 if out=="Incident AMI Case"
replace outcome = 4 if out=="Incident IHD Case"
destring outcome, replace


local i = 0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		putexcel set RESULTS_FILE, sheet(`out') modify
	
local i=`i'+1
local x=14

metan or_education lci_education uci_education if outcome==`i' , lcols(sample n ) effect(Risk difference) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) 
local estimate_education = r(ES)
	local lci_education = r(ci_low)
	local uci_education = r(ci_upp)	
	putexcel G`x'=`estimate_education' H`x'=`lci_education' I`x'=`uci_education'

metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==`i' , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
local beta_zcsi = r(ES)
	local LCI_zcsi = r(ci_low)
	local UCI_zcsi = r(ci_upp)
	
	putexcel N`x'=`beta_zcsi' O`x'=`LCI_zcsi' P`x'=`UCI_zcsi'
}

********************************************************************************

*SBP on the RD scale

import excel "SBP_results.xlsx", sheet("all") firstrow clear
destring sample, replace
lab var sample "Sample"
lab def sample 1 "Sample 1" 2 "Sample 2"
lab val sample sample

lab var n "N"

generate outcome = .
replace outcome = 1 if out=="Incident CVD Case"
replace outcome = 2 if out=="Incident Stroke Case"
replace outcome = 3 if out=="Incident AMI Case"
replace outcome = 4 if out=="Incident IHD Case"
destring outcome, replace


local i = 0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		putexcel set RESULTS_FILE, sheet(`out') modify

local i=`i'+1
local x = 13
metan or_education lci_education uci_education if outcome==`i', lcols(sample n ) effect(log OR) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6)
local estimate_education = r(ES)
	local lci_education = r(ci_low)
	local uci_education = r(ci_upp)	
	putexcel G`x'=`estimate_education' H`x'=`lci_education' I`x'=`uci_education'

metan beta_bp LCI_bp UCI_bp if outcome==`i', lcols(sample n ) effect(Risk difference) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6)
local estimate_sbp = r(ES)
	local lci_sbp = r(ci_low)
	local uci_sbp = r(ci_upp)
	
	putexcel N`x'=`estimate_sbp' O`x'=`lci_sbp' P`x'=`uci_sbp'
}

********************************************************************************

*zcsi on the logOR scale


import excel "zcsi_results_logOR.xlsx", sheet("all") firstrow clear

destring sample, replace
lab var sample "Sample"
lab def sample 1 "Sample 1" 2 "Sample 2"
lab val sample sample


generate outcome = .
replace outcome = 1 if out=="Incident CVD Case"
replace outcome = 2 if out=="Incident Stroke Case"
replace outcome = 3 if out=="Incident AMI Case"
replace outcome = 4 if out=="Incident IHD Case"
destring outcome, replace


lab var n "N"

local i = 0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		putexcel set RESULTS_FILE, sheet(`out') modify
	
local i=`i'+1
local x = 18
metan or_education lci_education uci_education if outcome==`i', lcols(sample n ) effect(log OR) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6)
local estimate_education = r(ES)
	local lci_education = r(ci_low)
	local uci_education = r(ci_upp)	
	putexcel G`x'=`estimate_education' H`x'=`lci_education' I`x'=`uci_education'

metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==`i', lcols(sample n ) effect(log OR) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6)
local beta_zcsi = r(ES)
	local LCI_zcsi = r(ci_low)
	local UCI_zcsi = r(ci_upp)
	
	putexcel N`x'=`beta_zcsi' O`x'=`LCI_zcsi' P`x'=`UCI_zcsi'
}

********************************************************************************
	
*SBP on the logOR scale


import excel "SBP_results_logOR.xlsx", sheet("all") firstrow clear

destring sample, replace
lab var sample "Sample"
lab def sample 1 "Sample 1" 2 "Sample 2"
lab val sample sample


generate outcome = .
replace outcome = 1 if out=="Incident CVD Case"
replace outcome = 2 if out=="Incident Stroke Case"
replace outcome = 3 if out=="Incident AMI Case"
replace outcome = 4 if out=="Incident IHD Case"
destring outcome, replace


lab var n "N"

local i = 0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		putexcel set RESULTS_FILE, sheet(`out') modify

local i=`i'+1
local x = 17
metan or_education lci_education uci_education if outcome==`i', lcols(sample n ) effect(log OR) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6)
local estimate_education = r(ES)
	local lci_education = r(ci_low)
	local uci_education = r(ci_upp)	
	putexcel G`x'=`estimate_education' H`x'=`lci_education' I`x'=`uci_education'

metan beta_sbp LCI_sbp UCI_sbp if outcome==`i', lcols(sample n ) effect(log OR) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6)
local estimate_sbp = r(ES)
	local lci_sbp = r(ci_low)
	local uci_sbp = r(ci_upp)
	
	putexcel N`x'=`estimate_sbp' O`x'=`lci_sbp' P`x'=`uci_sbp'
}
******* Direct Effects for Multiple Mediators - One-sample MR on RD scale ******

import excel "SPLIT_SAMPLE_MULTIPLEMED.xlsx", sheet("all") firstrow clear

destring sample, replace
lab var sample "Sample"
lab def sample 1 "Sample 1" 2 "Sample 2"
lab val sample sample


generate outcome = .
replace outcome = 1 if out=="Incident CVD Case"
replace outcome = 2 if out=="Incident Stroke Case"
replace outcome = 3 if out=="Incident AMI Case"
replace outcome = 4 if out=="Incident IHD Case"
destring outcome, replace


lab var n "N"

local i = 0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		putexcel set RESULTS_FILE, sheet(`out') modify

local i=`i'+1
local x = 15
metan or_education_controlled_for_med lci_education_controlled_for_med uci_education_controlled_for_med if outcome==`i', lcols(sample n ) effect(RD) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6)
local estimate_education = r(ES)
	local lci_education = r(ci_low)
	local uci_education = r(ci_upp)	
	putexcel G`x'=`estimate_education' H`x'=`lci_education' I`x'=`uci_education'

}




******* Direct Effects for Multiple Mediators - One-sample MR on log OR scale ******

import excel "MultipledMed_logOR.xlsx", sheet("all") firstrow clear

destring sample, replace
lab var sample "Sample"
lab def sample 1 "Sample 1" 2 "Sample 2"
lab val sample sample


generate outcome = .
replace outcome = 1 if out=="Incident CVD Case"
replace outcome = 2 if out=="Incident Stroke Case"
replace outcome = 3 if out=="Incident AMI Case"
replace outcome = 4 if out=="Incident IHD Case"
destring outcome, replace


lab var n "N"

local i = 0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		putexcel set RESULTS_FILE, sheet(`out') modify

local i=`i'+1
local x = 19
metan or lci uci if outcome==`i', lcols(sample n ) effect(RD) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6)
local estimate_education = r(ES)
	local lci_education = r(ci_low)
	local uci_education = r(ci_upp)	
	putexcel G`x'=`estimate_education' H`x'=`lci_education' I`x'=`uci_education'

}

