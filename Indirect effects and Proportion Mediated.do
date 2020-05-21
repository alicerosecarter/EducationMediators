
global PCs PC01 PC02 PC03 PC04 PC05 PC06 PC07 PC08 PC09 PC10
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth

set more off

***********************************************
* Observational mediators on the log OR scale *
***********************************************

capture program drop boot
program boot, rclass
foreach med in zbmi zsbp zcsi { 
        
     if (`med' == zbmi) {
        regress zbmi eduyears_scaled $max_covar
        scalar exp_zbmi = _b[eduyears_scaled]
        }
        if (`med' == zsbp) {
        regress zsbp eduyears_scaled $max_covar
        scalar exp_zsbp = _b[eduyears_scaled]
        }
        if (`med' == zcsi) {
        regress zcsi eduyears_scaled $max_covar
        scalar exp_zcsi = _b[eduyears_scaled]
        }

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {     

     
     if (`out' == CVD_inc) {
	logit CVD_inc eduyears_scaled $max_covar
	scalar total_CVD = _b[eduyears_scaled]
        logit CVD_inc `med' eduyears_scaled $max_covar 
        scalar med_CVD = _b[`med']
		scalar direct_CVD = _b[eduyears_scaled]
		scalar indirect_diff_`med'_CVD_inc = total_CVD-direct_CVD
        scalar indirect_`med'_CVD_inc = exp_`med'*med_CVD
	scalar prop_`med'_CVD_inc = (indirect_`med'_CVD_inc/total_CVD)*100
	scalar diff_prop_`med'_CVD_inc = (indirect_diff_`med'_CVD_inc/total_CVD)*100
        }
        if (`out' == Stroke_inc) {
	logit Stroke_inc eduyears_scaled $max_covar
	scalar total_Stroke = _b[eduyears_scaled]
        logit Stroke_inc `med' eduyears_scaled $max_covar 
        scalar med_Stroke = _b[`med']
		scalar direct_Stroke = _b[eduyears_scaled]
		scalar indirect_diff_`med'_Stroke_inc = total_Stroke-direct_Stroke
        scalar indirect_`med'_Stroke_inc = exp_`med'*med_Stroke
	scalar prop_`med'_Stroke_inc = (indirect_`med'_Stroke_inc/total_Stroke)*100
	scalar diff_prop_`med'_Stroke_inc = (indirect_diff_`med'_Stroke_inc/total_Stroke)*100
        }
        if (`out' == AMI_inc) {
	logit AMI_inc eduyears_scaled $max_covar
	scalar total_AMI = _b[eduyears_scaled]
        logit AMI_inc `med' eduyears_scaled $max_covar 
        scalar med_AMI = _b[`med']
		scalar direct_AMI = _b[eduyears_scaled]
		scalar indirect_diff_`med'_AMI_inc = total_AMI-direct_AMI
        scalar indirect_`med'_AMI_inc = exp_`med'*med_AMI
	scalar prop_`med'_AMI_inc = (indirect_`med'_AMI_inc/total_AMI)*100
	scalar diff_prop_`med'_AMI_inc = (indirect_diff_`med'_AMI_inc/total_AMI)*100
        }
        if (`out' == IHD_inc) {
	logit IHD_inc eduyears_scaled $max_covar
	scalar total_IHD = _b[eduyears_scaled]
        logit IHD_inc `med' eduyears_scaled $max_covar 
        scalar med_IHD = _b[`med']
		scalar direct_IHD = _b[eduyears_scaled]
		scalar indirect_diff_`med'_IHD_inc = total_IHD-direct_IHD
        scalar indirect_`med'_IHD_inc = exp_`med'*med_IHD
	scalar prop_`med'_IHD_inc = (indirect_`med'_IHD_inc/total_IHD)*100
	scalar diff_prop_`med'_IHD_inc = (indirect_diff_`med'_IHD_inc/total_IHD)*100
        }
}
}
 end     

 foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc{
 local x=5
 local y=5
 local z=5
 foreach med in zbmi zsbp zcsi { 
 
     putexcel set RESULTS_FILE, sheet(`out') modify
		reg `med'  eduyears_scaled $max_covar
         logit `out' `med' eduyears_scaled $max_covar 
 bootstrap indirect = indirect_`med'_`out', reps(5) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]
 local x=`x'+1

		putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI' 
  
		 logit `out' eduyears_scaled $max_covar
         reg `med'  eduyears_scaled $max_covar
         logit `out' `med' eduyears_scaled $max_covar  	
bootstrap prop_med = prop_`med'_`out', reps(5) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 local y = `y'+1 
 
      		putexcel X`y'=`proportion' Y`y'=`proportion_LCI' Z`y'= `proportion_UCI' 

	logit `out'  eduyears_scaled $max_covar
	logit `out'  eduyears_scaled `med' $max_covar
bootstrap diff_prop_`med'_`out', reps(5) nodrop : boot 
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 local z = `z'+1
	 
   		putexcel U`z'=`proportion' V`z'=`proportion_LCI' W`z'= `proportion_UCI'


 }
}

*******************************************************
* Observational multiple mediators on the LogOR scale *
*******************************************************

global med zbmi zsbp zcsi 
capture program drop boot
program boot, rclass

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {     

     
     if (`out' == CVD_inc) {
	logit CVD_inc eduyears_scaled $max_covar
	scalar total_CVD = _b[eduyears_scaled]
        logit CVD_inc $med eduyears_scaled $max_covar 
		scalar direct_CVD = _b[eduyears_scaled]
		scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar diff_prop_CVD_inc = (indirect_diff_CVD_inc/total_CVD)*100
        }
     if (`out' == Stroke_inc) {
	logit Stroke_inc eduyears_scaled $max_covar
	scalar total_Stroke = _b[eduyears_scaled]
        logit Stroke_inc $med eduyears_scaled $max_covar 
		scalar direct_Stroke = _b[eduyears_scaled]
		scalar indirect_diff_Stroke_inc = total_Stroke-direct_Stroke
	scalar diff_prop_Stroke_inc = (indirect_diff_Stroke_inc/total_Stroke)*100
        }
     if (`out' == AMI_inc) {
	logit AMI_inc eduyears_scaled $max_covar
	scalar total_AMI = _b[eduyears_scaled]
        logit AMI_inc $med eduyears_scaled $max_covar 
		scalar direct_AMI = _b[eduyears_scaled]
		scalar indirect_diff_AMI_inc = total_AMI-direct_AMI
	scalar diff_prop_AMI_inc = (indirect_diff_AMI_inc/total_AMI)*100
        }
     if (`out' == IHD_inc) {
	logit IHD_inc eduyears_scaled $max_covar
	scalar total_IHD = _b[eduyears_scaled]
        logit IHD_inc $med eduyears_scaled $max_covar 
		scalar direct_IHD = _b[eduyears_scaled]
		scalar indirect_diff_IHD_inc = total_IHD-direct_IHD
	scalar diff_prop_IHD_inc = (indirect_diff_IHD_inc/total_IHD)*100
        }
}

 end     

 foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc{
 local x=9
 local y=9
 local z=9
     	 putexcel set RESULTS_FILE, sheet(`out') modify
	 logit `out'  eduyears_scaled $max_covar
	logit `out'  eduyears_scaled $med $max_covar
 
 bootstrap indirect = indirect_diff_`out', reps(5) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]


     	
		 putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'
  	 
		
	logit `out'  eduyears_scaled $max_covar
	logit `out'  eduyears_scaled $med $max_covar
bootstrap diff_prop_`out', reps(5) nodrop : boot 
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]

	 
   		putexcel U`y'=`proportion' V`y'=`proportion_LCI' W`y'= `proportion_UCI'
 }


**************************************************
* Observational mediators on the Risk diff scale *
**************************************************

capture program drop boot
program boot, rclass
foreach med in zbmi zsbp zcsi { 
        
     if (`med' == zbmi) {
        regress zbmi eduyears_scaled $max_covar
        scalar exp_zbmi = _b[eduyears_scaled]
        }
        if (`med' == zsbp) {
        regress zsbp eduyears_scaled $max_covar
        scalar exp_zsbp = _b[eduyears_scaled]
        }
        if (`med' == zcsi) {
        regress zcsi eduyears_scaled $max_covar
        scalar exp_zcsi = _b[eduyears_scaled]
        }

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {     

     
     if (`out' == CVD_inc) {
	reg CVD_inc eduyears_scaled $max_covar
	scalar total_CVD = _b[eduyears_scaled]
        regress CVD_inc `med' eduyears_scaled $max_covar 
        scalar med_CVD = _b[`med']
		scalar direct_CVD = _b[eduyears_scaled]
		scalar indirect_diff_`med'_CVD_inc = total_CVD-direct_CVD
        scalar indirect_`med'_CVD_inc = exp_`med'*med_CVD
	scalar prop_`med'_CVD_inc = (indirect_`med'_CVD_inc/total_CVD)*100
	scalar diff_prop_`med'_CVD_inc = (indirect_diff_`med'_CVD_inc/total_CVD)*100
        }
        if (`out' == Stroke_inc) {
	reg Stroke_inc eduyears_scaled $max_covar
	scalar total_Stroke = _b[eduyears_scaled]
        regress Stroke_inc `med' eduyears_scaled $max_covar 
        scalar med_Stroke = _b[`med']
		scalar direct_Stroke = _b[eduyears_scaled]
		scalar indirect_diff_`med'_Stroke_inc = total_Stroke-direct_Stroke
        scalar indirect_`med'_Stroke_inc = exp_`med'*med_Stroke
	scalar prop_`med'_Stroke_inc = (indirect_`med'_Stroke_inc/total_Stroke)*100
	scalar diff_prop_`med'_Stroke_inc = (indirect_diff_`med'_Stroke_inc/total_Stroke)*100
        }
        if (`out' == AMI_inc) {
	reg AMI_inc eduyears_scaled $max_covar
	scalar total_AMI = _b[eduyears_scaled]
        regress AMI_inc `med' eduyears_scaled $max_covar 
        scalar med_AMI = _b[`med']
		scalar direct_AMI = _b[eduyears_scaled]
		scalar indirect_diff_`med'_AMI_inc = total_AMI-direct_AMI
        scalar indirect_`med'_AMI_inc = exp_`med'*med_AMI
	scalar prop_`med'_AMI_inc = (indirect_`med'_AMI_inc/total_AMI)*100
	scalar diff_prop_`med'_AMI_inc = (indirect_diff_`med'_AMI_inc/total_AMI)*100
        }
        if (`out' == IHD_inc) {
	reg IHD_inc eduyears_scaled $max_covar
	scalar total_IHD = _b[eduyears_scaled]
        regress IHD_inc `med' eduyears_scaled $max_covar 
        scalar med_IHD = _b[`med']
		scalar direct_IHD = _b[eduyears_scaled]
		scalar indirect_diff_`med'_IHD_inc = total_IHD-direct_IHD
	    scalar indirect_`med'_IHD_inc = exp_`med'*med_IHD
		scalar prop_`med'_IHD_inc = (indirect_`med'_IHD_inc/total_IHD)*100
		scalar diff_prop_`med'_IHD_inc = (indirect_diff_`med'_IHD_inc/total_IHD)*100

        }
}
}
 end     

 foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc{
 local x=1
 local y=1
 local z=1
 foreach med in zbmi zsbp zcsi { 
     putexcel set RESULTS_FILE, sheet(`out') modify
       reg `med'  eduyears_scaled $max_covar
         reg `out' `med' eduyears_scaled $max_covar 
 bootstrap indirect = indirect_`med'_`out', reps(5) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]
         local x=`x'+1

     	
		 putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'
  	 
	 reg `out' eduyears_scaled $max_covar
         reg `med'  eduyears_scaled $max_covar
         reg `out' `med' eduyears_scaled $max_covar 

bootstrap prop_med = prop_`med'_`out', reps(5) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
         local y=`y'+1

      		putexcel X`y'=`proportion' Y`y'=`proportion_LCI' Z`y'= `proportion_UCI' 
		
	reg `out'  eduyears_scaled $max_covar
	reg `out'  eduyears_scaled `med' $max_covar
bootstrap prop_med_diff = diff_prop_`med'_`out', reps(5) nodrop : boot 
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 local z=`z'+1
	 
   		putexcel U`z'=`proportion' V`z'=`proportion_LCI' W`z'=`proportion_UCI'
 }
}


***********************************************************
* Observational multiple mediators on the Risk diff scale *
***********************************************************

global med zbmi zsbp zcsi 
capture program drop boot
program boot, rclass

foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {     

     
     if (`out' == CVD_inc) {
	reg CVD_inc eduyears_scaled $max_covar
	scalar total_CVD = _b[eduyears_scaled]
        regress CVD_inc $med eduyears_scaled $max_covar 
		scalar direct_CVD = _b[eduyears_scaled]
		scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar diff_prop_CVD_inc = (indirect_diff_CVD_inc/total_CVD)*100
        }
     if (`out' == Stroke_inc) {
	reg Stroke_inc eduyears_scaled $max_covar
	scalar total_Stroke = _b[eduyears_scaled]
        regress Stroke_inc $med eduyears_scaled $max_covar 
		scalar direct_Stroke = _b[eduyears_scaled]
		scalar indirect_diff_Stroke_inc = total_Stroke-direct_Stroke
	scalar diff_prop_Stroke_inc = (indirect_diff_Stroke_inc/total_Stroke)*100
        }
     if (`out' == AMI_inc) {
	reg AMI_inc eduyears_scaled $max_covar
	scalar total_AMI = _b[eduyears_scaled]
        regress AMI_inc $med eduyears_scaled $max_covar 
		scalar direct_AMI = _b[eduyears_scaled]
		scalar indirect_diff_AMI_inc = total_AMI-direct_AMI
	scalar diff_prop_AMI_inc = (indirect_diff_AMI_inc/total_AMI)*100
        }
     if (`out' == IHD_inc) {
	reg IHD_inc eduyears_scaled $max_covar
	scalar total_IHD = _b[eduyears_scaled]
        regress IHD_inc $med eduyears_scaled $max_covar 
		scalar direct_IHD = _b[eduyears_scaled]
		scalar indirect_diff_IHD_inc = total_IHD-direct_IHD
	scalar diff_prop_IHD_inc = (indirect_diff_IHD_inc/total_IHD)*100
        }
}

 end     

 foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc{
 local x=5
 local y=5
 local z=5
     	 putexcel set RESULTS_FILE, sheet(`out') modify
         reg `out' eduyears_scaled $max_covar
         reg `out' $med eduyears_scaled $max_covar 
 
 bootstrap indirect = indirect_diff_`out', reps(5) nodrop : boot  // indirect effect
         matrix results = r(table)
         local indirect = results[1,1]
         local indirect_LCI = results[5,1]
         local indirect_UCI = results[6,1]

     	
		 putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'
  	 
		
	reg `out'  eduyears_scaled $max_covar
	reg `out'  eduyears_scaled $med $max_covar
bootstrap diff_prop_`out', reps(5) nodrop : boot 
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
	 
   		putexcel U`y'=`proportion' V`y'=`proportion_LCI' W`y'= `proportion_UCI'
 }

****************************************
* BMI MR mediation on the log OR scale *
****************************************

capture drop dir_eabmi
capture drop dir_bmiea
capture drop education
regress eduyears_scaled ea_weighted  $max_covar $PCs
predict education, xb
regress eduyears_scaled ea_weighted zbmi_weighted $max_covar $PCs
predict dir_eabmi, xb
regress zbmi zbmi_weighted ea_weighted $max_covar $PCs
predict dir_bmiea, xb


capture program drop boot
program boot, rclass
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		
	if (`out' == CVD_inc)  {
	logit CVD_inc education $max_covar $PCs 
	scalar total_CVD = _b[education]
	logit CVD_inc dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
	scalar med_CVD = _b[dir_bmiea]
	scalar direct_CVD = _b[dir_eabmi]
	scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar indirect_CVD_inc = exp_med*med_CVD
	scalar prop_bmi_CVD_inc = (indirect_CVD_inc/total_CVD)*100 
	scalar diff_prop_bmi_CVD_inc = (indirect_diff_CVD_inc/total_CVD)*100
	}
	if (`out' == Stroke_inc) {
	logit Stroke_inc education $max_covar $PCs 
	scalar total_Stroke = _b[education]
	logit Stroke_inc dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
	scalar med_Stroke = _b[dir_bmiea]
	scalar direct_stroke = _b[dir_eabmi]
	scalar indirect_diff_Stroke_inc = total_Stroke-direct_Stroke
	scalar indirect_Stroke_inc = exp_med*med_Stroke
	scalar prop_bmi_Stroke_inc = (indirect_Stroke_inc/total_Stroke)*100 
	scalar diff_prop_bmi_Stroke_inc = (indirect_diff_Stroke_inc/total_Stroke)*100

	}
	if (`out' == AMI_inc) {
	logit AMI_inc education $max_covar $PCs 
	scalar total_AMI = _b[education]
	logit AMI_inc dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
	scalar med_AMI = _b[dir_bmiea]
	scalar direct_AMI = _b[dir_eabmi]
	scalar indirect_diff_AMI_inc = total_AMI-direct_AMI
	scalar indirect_AMI_inc = exp_med*med_AMI
	scalar prop_bmi_AMI_inc = (indirect_AMI_inc/total_AMI)*100 
	scalar diff_prop_bmi_AMI_inc = (indirect_diff_AMI_inc/total_AMI)*100

	}
	if (`out' == IHD_inc) {
	logit IHD_inc education $max_covar $PCs 
	scalar total_IHD = _b[education]
	logit IHD_inc dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
	scalar med_IHD = _b[dir_bmiea]
	scalar direct_IHD = _b[dir_eabmi]
	scalar indirect_diff_IHD_inc = total_IHD-direct_IHD
	scalar indirect_IHD_inc = exp_med*med_IHD
	scalar prop_bmi_IHD_inc = (indirect_IHD_inc/total_IHD)*100
	scalar diff_prop_bmi_IHD_inc = (indirect_diff_IHD_inc/total_IHD)*100

	}
	}
	
end	
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
local x=16
local y=16
local z=16
	putexcel set RESULTS_FILE, sheet(`out') modify
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	logit `out' dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 
bootstrap indirect = indirect_`out', reps(5) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]

	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	logit `out' education $max_covar $PCs 
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	logit `out' dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 

bootstrap prop_med = prop_bmi_`out', reps(5) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	putexcel X`y'=`proportion' Y`y'=`proportion_LCI' Z`y'= `proportion_UCI'

	
	logit `out' education $max_covar $PCs 
	logit `out' dir_bmiea dir_eabmi  $max_covar $PCs, vce(robust) 

bootstrap diff_prop_med = diff_prop_bmi_`out', reps(5) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	putexcel U`z'=`proportion' V`z'=`proportion_LCI' W`z'= `proportion_UCI'
}

**************************************
*SBP MR mediation on the log OR scale*
**************************************

capture drop education
regress eduyears_scaled ea_weighted $PCs $max_covar
predict education, xb

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

append using "sbp_results.dta"

capture program drop boot
program boot, rclass
	ivreg2 zsbp (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		
	if (`out' == CVD_inc)  {
	logit CVD_inc education $max_covar $PCs 
	scalar total_CVD = _b[education]
	metan beta_sbp LCI_sbp UCI_sbp if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_CVD = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_CVD = `estimate_CVD'
	metan or_education lci_education uci_education if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_CVD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_CVD = `direct_CVD'
	scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar indirect_CVD_inc = exp_med*med_CVD
	scalar prop_sbp_CVD_inc = (indirect_CVD_inc/total_CVD)*100
	scalar diff_prop_sbp_CVD_inc = (indirect_diff_CVD_inc/total_CVD)*100
	}
	if (`out' == Stroke_inc) {
	logit Stroke_inc education $max_covar $PCs 
	scalar total_Stroke = _b[education]
	metan beta_sbp LCI_sbp UCI_sbp if outcome==2 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_Stroke = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_Stroke = `estimate_Stroke'
	metan or_education lci_education uci_education if outcome==2 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_Stroke = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_Stroke = `direct_Stroke'
	scalar indirect_diff_Stroke_inc = total_Stroke-direct_Stroke
	scalar indirect_Stroke_inc = exp_med*med_Stroke
	scalar prop_sbp_Stroke_inc = (indirect_Stroke_inc/total_Stroke)*100
	scalar diff_prop_sbp_Stroke_inc = (indirect_diff_Stroke_inc/total_Stroke)*100
	}
	if (`out' == AMI_inc) {
	logit AMI_inc education $max_covar $PCs 
	scalar total_AMI = _b[education]
	metan beta_sbp LCI_sbp UCI_sbp if outcome==3 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_AMI = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_AMI = `estimate_AMI'
	metan or_education lci_education uci_education if outcome==3 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_AMI = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_AMI = `direct_AMI'
	scalar indirect_diff_AMI_inc = total_AMI-direct_AMI
	scalar indirect_AMI_inc = exp_med*med_AMI
	scalar prop_sbp_AMI_inc = (indirect_AMI_inc/total_AMI)*100
	scalar diff_prop_sbp_AMI_inc = (indirect_diff_AMI_inc/total_AMI)*100
	}
	if (`out' == IHD_inc) {
	logit IHD_inc education $max_covar $PCs 
	scalar total_IHD = _b[education]
	metan beta_sbp LCI_sbp UCI_sbp if outcome==4 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_IHD = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_IHD = `estimate_IHD'
	metan or_education lci_education uci_education if outcome==4 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_IHD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_IHD = `direct_IHD'
	scalar indirect_diff_IHD_inc = total_IHD-direct_IHD
	scalar indirect_IHD_inc = exp_med*med_IHD
	scalar prop_sbp_IHD_inc = (indirect_IHD_inc/total_IHD)*100
	scalar diff_prop_sbp_IHD_inc = (indirect_diff_IHD_inc/total_IHD)*100
	}
	}
	
end	
local i=0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc{
local i=`i'+1
local x=17
local y=17
local z=17
	putexcel set RESULTS_FILE, sheet(`out') modify
	ivreg2 zsbp (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	metan beta_sbp LCI_sbp UCI_sbp if outcome==`i' , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
bootstrap indirect = indirect_`out', reps(5) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]

	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	logit `out' education $max_covar $PCs 
	ivreg2 zsbp (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	metan beta_sbp LCI_sbp UCI_sbp if outcome==`i' , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph

bootstrap prop_med = prop_sbp_`out', reps(5) nodrop : boot  // indirect effect
        matrix results = r(table)
        local proportion = results[1,1]
        local proportion_LCI = results[5,1]
        local proportion_UCI = results[6,1]
        
	putexcel X`y'=`proportion' Y`y'=`proportion_LCI' Z`y'= `proportion_UCI'

	
	logit `out' education $max_covar $PCs	
	metan or_education lci_education uci_education if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
bootstrap diff_prop_med = diff_prop_sbp_`out', reps(5) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	putexcel U`z'=`proportion' V`z'=`proportion_LCI' W`z'= `proportion_UCI'
}

**************************************
*CSI MR Mediation on the log OR scale*
**************************************
use "indirect_bootstrap_CIs_data.dta", clear
global PCs PC01 PC02 PC03 PC04 PC05 PC06 PC07 PC08 PC09 PC10
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth

set more off

capture drop education
regress eduyears_scaled ea_weighted $PCs $max_covar
predict education, xb

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

append using "CSI_results"

capture program drop boot
program boot, rclass
	ivreg2 zcsi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		
	if (`out' == CVD_inc)  {
	logit CVD_inc education $max_covar $PCs 
	scalar total_CVD = _b[education]
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_CVD = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_CVD = `estimate_CVD'
	metan or_education lci_education uci_education if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_CVD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_CVD = `direct_CVD'
	scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar indirect_CVD_inc = exp_med*med_CVD
	scalar prop_zcsi_CVD_inc = (indirect_CVD_inc/total_CVD)*100
	scalar diff_prop_zcsi_CVD_inc = (indirect_diff_CVD_inc/total_CVD)*100
	}
	if (`out' == Stroke_inc) {
	logit Stroke_inc education $max_covar $PCs 
	scalar total_Stroke = _b[education]
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==2 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_Stroke = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_Stroke = `estimate_Stroke'
	metan or_education lci_education uci_education if outcome==2 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_Stroke = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_Stroke = `direct_Stroke'
	scalar indirect_diff_Stroke_inc = total_Stroke-direct_Stroke
	scalar indirect_Stroke_inc = exp_med*med_Stroke
	scalar prop_zcsi_Stroke_inc = (indirect_Stroke_inc/total_Stroke)*100
	scalar diff_prop_zcsi_Stroke_inc = (indirect_diff_Stroke_inc/total_Stroke)*100
	}
	if (`out' == AMI_inc) {
	logit AMI_inc education $max_covar $PCs 
	scalar total_AMI = _b[education]
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==3 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_AMI = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_AMI = `estimate_AMI'
	metan or_education lci_education uci_education if outcome==3 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_AMI = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_AMI = `direct_AMI'
	scalar indirect_diff_AMI_inc = total_AMI-direct_AMI
	scalar indirect_AMI_inc = exp_med*med_AMI
	scalar prop_zcsi_AMI_inc = (indirect_AMI_inc/total_AMI)*100
	scalar diff_prop_zcsi_AMI_inc = (indirect_diff_AMI_inc/total_AMI)*100
	}
	if (`out' == IHD_inc) {
	logit IHD_inc education $max_covar $PCs 
	scalar total_IHD = _b[education]
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==4 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_IHD = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_IHD = `estimate_IHD'
	metan or_education lci_education uci_education if outcome==4 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_IHD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_IHD = `direct_IHD'
	scalar indirect_diff_IHD_inc = total_IHD-direct_IHD
	scalar indirect_IHD_inc = exp_med*med_IHD
	scalar prop_zcsi_IHD_inc = (indirect_IHD_inc/total_IHD)*100
	scalar diff_prop_zcsi_IHD_inc = (indirect_diff_IHD_inc/total_IHD)*100
	}
	}
	
end	
local i = 0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc{
local i=`i'+1
local x=18
local y=18
local z=18
	putexcel set RESULTS_FILE, sheet(`out') modify
	ivreg2 zcsi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==`i' , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
bootstrap indirect = indirect_`out', reps(5) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]

	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	logit `out' education $max_covar $PCs 
	ivreg2 zcsi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==`i' , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph

bootstrap prop_med = prop_zcsi_`out', reps(5) nodrop : boot  // indirect effect
        matrix results = r(table)
        local proportion = results[1,1]
        local proportion_LCI = results[5,1]
        local proportion_UCI = results[6,1]
	
	putexcel X`y'=`proportion' Y`y'=`proportion_LCI' Z`y'= `proportion_UCI'

		metan or_education lci_education uci_education if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
		logit `out' education $max_covar $PCs 

bootstrap diff_prop_med = diff_prop_zcsi_`out', reps(5) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	putexcel U`z'=`proportion' V`z'=`proportion_LCI' W`z'= `proportion_UCI'

}

*******************************************
* BMI MR mediation on the Risk diff scale *
*******************************************

capture program drop boot
program boot, rclass
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		
	if (`out' == CVD_inc)  {
	ivreg2 CVD_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_CVD = _b[eduyears_scaled]
	ivreg2 CVD_inc (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
	scalar med_CVD = _b[zbmi]
	scalar direct_CVD = _b[eduyears_scaled]
	scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar indirect_CVD_inc = exp_med*med_CVD
	scalar prop_bmi_CVD_inc = (indirect_CVD_inc/total_CVD)*100 
	scalar diff_prop_bmi_`out' = (indirect_diff_CVD_inc/total_CVD)*100 
	}
	if (`out' == Stroke_inc) {
	ivreg2 Stroke_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_Stroke = _b[eduyears_scaled]
	ivreg2 Stroke_inc (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
	scalar med_Stroke = _b[zbmi]
	scalar direct_Stroke = _b[eduyears_scaled]
	scalar indirect_diff_Stroke_inc = total_Stroke-direct_Stroke
	scalar indirect_Stroke_inc = exp_med*med_Stroke
	scalar prop_bmi_Stroke_inc = (indirect_Stroke_inc/total_Stroke)*100 
	scalar diff_prop_bmi_`out' = (indirect_diff_Stroke_inc/total_Stroke)*100 
	}
	if (`out' == AMI_inc) {
	ivreg2 AMI_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_AMI = _b[eduyears_scaled]
	ivreg2 AMI_inc (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
	scalar med_AMI = _b[zbmi]
	scalar direct_AMI = _b[eduyears_scaled]
	scalar indirect_diff_AMI_inc = total_AMI-direct_AMI
	scalar indirect_AMI_inc = exp_med*med_AMI
	scalar prop_bmi_AMI_inc = (indirect_AMI_inc/total_AMI)*100 
	scalar diff_prop_bmi_`out' = (indirect_diff_AMI_inc/total_AMI)*100 
	}
	if (`out' == IHD_inc) {
	ivreg2 IHD_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_IHD = _b[eduyears_scaled]
	ivreg2 IHD_inc (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
	scalar med_IHD = _b[zbmi]
	scalar direct_IHD = _b[eduyears_scaled]
	scalar indirect_diff_IHD_inc = total_IHD-direct_IHD
	scalar indirect_IHD_inc = exp_med*med_IHD
	scalar prop_bmi_IHD_inc = (indirect_IHD_inc/total_IHD)*100 
	scalar diff_prop_bmi_`out' = (indirect_diff_IHD_inc/total_IHD)*100 
	}
	}

end	
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
local x=12
local y=12
local z=12
	putexcel set RESULTS_FILE, sheet(`out') modify
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
bootstrap indirect = indirect_`out', reps(5) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]
        	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	ivreg2 `out' (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 zbmi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar

bootstrap prop_med = prop_bmi_`out', reps(5) nodrop : boot  // indirect effect
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]


	putexcel X`y'=`proportion' Y`y'=`proportion_LCI' Z`y'= `proportion_UCI'
	
	ivreg2 `out' (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	ivreg2 `out' (eduyears_scaled zbmi = ea_weighted zbmi_weighted) $PCs $max_covar
bootstrap diff_prop_med = diff_prop_bmi_`out', reps(5) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	putexcel U`z'=`proportion' V`z'=`proportion_LCI' W`z'= `proportion_UCI'
}

*****************************************
*SBP MR Mediation on the risk diff scale*
*****************************************

use "indirect_bootstrap_CIs_data.dta", clear
global PCs PC01 PC02 PC03 PC04 PC05 PC06 PC07 PC08 PC09 PC10
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth

set more off

append using "sbp_results_risk_diff.dta"

rename beta_bp beta_sbp
rename LCI_bp LCI_sbp
rename UCI_bp UCI_sbp

capture program drop boot
program boot, rclass
	ivreg2 zsbp (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		
	
	if (`out' == CVD_inc)  {
	ivreg2 CVD_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_CVD = _b[eduyears_scaled]
	metan beta_sbp LCI_sbp UCI_sbp if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_CVD = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_CVD = `estimate_CVD'
	metan or_education lci_education uci_education if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_CVD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_CVD = `direct_CVD'
	scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar indirect_CVD_inc = exp_med*med_CVD
	scalar prop_sbp_CVD_inc = (indirect_CVD_inc/total_CVD)*100
	scalar diff_prop_sbp_CVD_inc = (indirect_diff_CVD_inc/total_CVD)*100
	}
	if (`out' == Stroke_inc) {
	ivreg2 Stroke_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_Stroke = _b[eduyears_scaled]
	metan beta_sbp LCI_sbp UCI_sbp if outcome==2 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_Stroke = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_Stroke = `estimate_Stroke'
	metan or_education lci_education uci_education if outcome==2 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_Stroke = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_Stroke = `direct_Stroke'
	scalar indirect_diff_Stroke_inc = total_Stroke-direct_Stroke
	scalar indirect_Stroke_inc = exp_med*med_Stroke
	scalar prop_sbp_Stroke_inc = (indirect_Stroke_inc/total_Stroke)*100
	scalar diff_prop_sbp_Stroke_inc = (indirect_diff_Stroke_inc/total_Stroke)*100
	}
	if (`out' == AMI_inc) {
	ivreg2 AMI_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_AMI = _b[eduyears_scaled]
	metan beta_sbp LCI_sbp UCI_sbp if outcome==3 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_AMI = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_AMI = `estimate_AMI'
	metan or_education lci_education uci_education if outcome==3 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_AMI = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_AMI = `direct_AMI'
	scalar indirect_diff_AMI_inc = total_AMI-direct_AMI
	scalar indirect_AMI_inc = exp_med*med_AMI
	scalar prop_sbp_AMI_inc = (indirect_AMI_inc/total_AMI)*100
	scalar diff_prop_sbp_AMI_inc = (indirect_diff_AMI_inc/total_AMI)*100
	}
	if (`out' == IHD_inc) {
	ivreg2 IHD_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_IHD = _b[eduyears_scaled]
	metan beta_sbp LCI_sbp UCI_sbp if outcome==4 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_IHD = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_IHD = `estimate_IHD'
	metan or_education lci_education uci_education if outcome==4 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_IHD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_IHD = `direct_IHD'
	scalar indirect_diff_IHD_inc = total_IHD-direct_IHD
	scalar indirect_IHD_inc = exp_med*med_IHD
	scalar prop_sbp_IHD_inc = (indirect_IHD_inc/total_IHD)*100
	scalar diff_prop_sbp_IHD_inc = (indirect_diff_IHD_inc/total_IHD)*100
	}
	}
end	
local i=0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc{
local i=`i'+1
local x=13
local y=13
local z=13
	putexcel set RESULTS_FILE, sheet(`out') modify
	ivreg2 zsbp (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	metan beta_sbp LCI_sbp UCI_sbp if outcome==`i' , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
bootstrap indirect = indirect_`out', reps(5) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]

	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	ivreg2 `out' (eduyears_scaled = ea_weighted) $max_covar $PCs 
	ivreg2 zsbp (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	metan or_education lci_education uci_education if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph

bootstrap prop_med = prop_sbp_`out', reps(5) nodrop : boot  // indirect effect
        matrix results = r(table)
        local proportion = results[1,1]
        local proportion_LCI = results[5,1]
        local proportion_UCI = results[6,1]
        
	putexcel X`y'=`proportion' Y`y'=`proportion_LCI' Z`y'= `proportion_UCI'
	
	ivreg2 `out' (eduyears_scaled = ea_weighted) $max_covar $PCs 
	metan or_education lci_education uci_education if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	
bootstrap diff_prop_med = diff_prop_sbp_`out', reps(5) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	putexcel U`z'=`proportion' V`z'=`proportion_LCI' W`z'= `proportion_UCI'

}

*****************************************
*CSI MR Mediation on the Risk Diff scale*
*****************************************

use "indirect_bootstrap_CIs_data.dta", clear

global PCs PC01 PC02 PC03 PC04 PC05 PC06 PC07 PC08 PC09 PC10
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth

set more off

append using "CSI_results_risk_diff.dta"

capture program drop boot
program boot, rclass
	ivreg2 zcsi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar exp_med = _b[eduyears_scaled]
 
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		
	
	if (`out' == CVD_inc)  {
	ivreg2 CVD_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_CVD = _b[eduyears_scaled]
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_CVD = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_CVD = `estimate_CVD'
	metan or_education lci_education uci_education if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_CVD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_CVD = `direct_CVD'
	scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar indirect_CVD_inc = exp_med*med_CVD
	scalar prop_zcsi_CVD_inc = (indirect_CVD_inc/total_CVD)*100
	scalar diff_prop_zcsi_CVD_inc = (indirect_diff_CVD_inc/total_CVD)*100
	}
	if (`out' == Stroke_inc) {
	ivreg2 Stroke_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_Stroke = _b[eduyears_scaled]
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==2 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_Stroke = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_Stroke = `estimate_Stroke'
	metan or_education lci_education uci_education if outcome==2 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_Stroke = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_Stroke = `direct_Stroke'
	scalar indirect_diff_Stroke_inc = total_Stroke-direct_Stroke
	scalar indirect_Stroke_inc = exp_med*med_Stroke
	scalar prop_zcsi_Stroke_inc = (indirect_Stroke_inc/total_Stroke)*100
	scalar diff_prop_zcsi_Stroke_inc = (indirect_diff_Stroke_inc/total_Stroke)*100
	}
	if (`out' == AMI_inc) {
	ivreg2 AMI_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_AMI = _b[eduyears_scaled]
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==3 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_AMI = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_AMI = `estimate_AMI'
	metan or_education lci_education uci_education if outcome==3 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_AMI = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_AMI = `direct_AMI'
	scalar indirect_diff_AMI_inc = total_AMI-direct_AMI
	scalar indirect_AMI_inc = exp_med*med_AMI
	scalar prop_zcsi_AMI_inc = (indirect_AMI_inc/total_AMI)*100
	scalar diff_prop_zcsi_AMI_inc = (indirect_diff_AMI_inc/total_AMI)*100
	}
	if (`out' == IHD_inc) {
	ivreg2 IHD_inc (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	scalar total_IHD = _b[eduyears_scaled]
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==4 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local estimate_IHD = r(ES)
	local lci = r(ci_low)
	local uci = r(ci_upp)	
	scalar med_IHD = `estimate_IHD'
	metan or_education lci_education uci_education if outcome==4 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_IHD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_IHD = `direct_IHD'
	scalar indirect_diff_IHD_inc = total_IHD-direct_IHD
	scalar indirect_IHD_inc = exp_med*med_IHD
	scalar prop_zcsi_IHD_inc = (indirect_IHD_inc/total_IHD)*100
	scalar diff_prop_zcsi_IHD_inc = (indirect_diff_IHD_inc/total_IHD)*100
	}
	}
end	
local i = 0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc{
local i=`i'+1
local x=14
local y=14
local z=14

	putexcel set RESULTS_FILE, sheet(`out') modify
	ivreg2 zcsi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==`i' , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
bootstrap indirect = indirect_`out', reps(5) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]

	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	ivreg2 `out' (eduyears_scaled = ea_weighted) $max_covar $PCs 
	ivreg2 zcsi (eduyears_scaled  = ea_weighted ) $PCs $max_covar
	metan beta_zcsi LCI_zcsi UCI_zcsi if outcome==`i' , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph

bootstrap prop_med = prop_zcsi_`out', reps(5) nodrop : boot  // indirect effect
        matrix results = r(table)
        local proportion = results[1,1]
        local proportion_LCI = results[5,1]
        local proportion_UCI = results[6,1]
        
	putexcel X`y'=`proportion' Y`y'=`proportion_LCI' Z`y'= `proportion_UCI'
	
	ivreg2 `out' (eduyears_scaled = ea_weighted) $max_covar $PCs 
	metan or_education lci_education uci_education if outcome==`i' , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	
bootstrap diff_prop_med = diff_prop_zcsi_`out', reps(5) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	putexcel U`z'=`proportion' V`z'=`proportion_LCI' W`z'= `proportion_UCI'

}


********************************************************
*Multiple mediators MR Mediation on the Risk Diff scale*
********************************************************

use "indirect_bootstrap_CIs_data.dta", clear

global PCs PC01 PC02 PC03 PC04 PC05 PC06 PC07 PC08 PC09 PC10
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth

append using "multiple_med_riskdiff.dta"

capture program drop boot
program boot, rclass
	
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		
	if (`out' == CVD_inc)  {
	ivreg2 CVD_inc (eduyears_scaled = ea_weighted) $maxcovar $PCs
	scalar total_CVD = _b[eduyears_scaled]
	metan or_education_controlled_for_med lci_education_controlled_for_med uci_education_controlled_for_med if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_CVD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_CVD = `direct_CVD'
	scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar diff_prop_CVD_inc = (indirect_diff_CVD_inc/total_CVD)*100
	}
	if (`out' == Stroke_inc)  {
	ivreg2 Stroke_inc (eduyears_scaled = ea_weighted) $maxcovar $PCs
	scalar total_Stroke = _b[eduyears_scaled]
	metan or_education_controlled_for_med lci_education_controlled_for_med uci_education_controlled_for_med if outcome==2 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_Stroke = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_Stroke = `direct_Stroke'
	scalar indirect_diff_Stroke_inc = total_Stroke-direct_Stroke
	scalar diff_prop_Stroke_inc = (indirect_diff_Stroke_inc/total_Stroke)*100
	}
	if (`out' == AMI_inc)  {
	ivreg2 AMI_inc (eduyears_scaled = ea_weighted) $maxcovar $PCs
	scalar total_AMI = _b[eduyears_scaled]
	metan or_education_controlled_for_med lci_education_controlled_for_med uci_education_controlled_for_med if outcome==3 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_AMI = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_AMI = `direct_AMI'
	scalar indirect_diff_AMI_inc = total_AMI-direct_AMI
	scalar diff_prop_AMI_inc = (indirect_diff_AMI_inc/total_AMI)*100
	}
	if (`out' == IHD_inc)  {
	ivreg2 IHD_inc (eduyears_scaled = ea_weighted) $maxcovar $PCs
	scalar total_IHD = _b[eduyears_scaled]
	metan or_education_controlled_for_med lci_education_controlled_for_med uci_education_controlled_for_med if outcome==4 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_IHD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_IHD = `direct_IHD'
	scalar indirect_diff_IHD_inc = total_IHD-direct_IHD
	scalar diff_prop_IHD_inc = (indirect_diff_IHD_inc/total_IHD)*100
	}
	}
end	
local i = 0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc{
local i=`i'+1
local x=15
local z=15

	putexcel set RESULTS_FILE, sheet(`out') modify
	ivreg2 `out' (eduyears_scaled = ea_weighted) $maxcovar $PCs
metan or_education_controlled_for_med lci_education_controlled_for_med uci_education_controlled_for_med if outcome==`i', lcols(sample n ) effect(RD) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
bootstrap indirect = indirect_diff_`out', reps(5) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]

	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	ivreg2 `out' (eduyears_scaled = ea_weighted) $max_covar $PCs 
metan or_education_controlled_for_med lci_education_controlled_for_med uci_education_controlled_for_med if outcome==`i', lcols(sample n ) effect(RD) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	
bootstrap diff_prop_med = diff_prop_`out', reps(5) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	putexcel U`z'=`proportion' V`z'=`proportion_LCI' W`z'= `proportion_UCI'

}

****************************************************
*Multiple mediators MR Mediation on the logOR scale*
****************************************************

use "indirect_bootstrap_CIs_data.dta", clear

global PCs PC01 PC02 PC03 PC04 PC05 PC06 PC07 PC08 PC09 PC10
global max_covar age sex northing_pob easting_pob cov_dist_lon TDI_birth

capture drop education
regress eduyears_scaled ea_weighted $PCs $max_covar
predict education, xb

append using "multiple_med_logOR.dta"


capture program drop boot
program boot, rclass
	
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc {
		
	if (`out' == CVD_inc)  {
	logit CVD_inc education $max_covar $PCs 
	scalar total_CVD = _b[education]
	metan or_education lci_education uci_education if outcome==1 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_CVD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_CVD = `direct_CVD'
	scalar indirect_diff_CVD_inc = total_CVD-direct_CVD
	scalar diff_prop_CVD_inc = (indirect_diff_CVD_inc/total_CVD)*100
	}
	if (`out' == Stroke_inc)  {
	logit Stroke_inc education $max_covar $PCs 
	scalar total_Stroke = _b[education]
	metan or_education lci_education uci_education if outcome==2 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_Stroke = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_Stroke = `direct_Stroke'
	scalar indirect_diff_Stroke_inc = total_Stroke-direct_Stroke
	scalar diff_prop_Stroke_inc = (indirect_diff_Stroke_inc/total_Stroke)*100
	}
	if (`out' == AMI_inc)  {
	logit AMI_inc education $max_covar $PCs 
	scalar total_AMI = _b[education]
	metan or_education lci_education uci_education if outcome==3 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_AMI = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_AMI = `direct_AMI'
	scalar indirect_diff_AMI_inc = total_AMI-direct_AMI
	scalar diff_prop_AMI_inc = (indirect_diff_AMI_inc/total_AMI)*100
	}
	if (`out' == IHD_inc)  {
	logit IHD_inc education $max_covar $PCs 
	scalar total_IHD = _b[education]
	metan or_education lci_education uci_education if outcome==4 , lcols(sample n ) effect(log OR) null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	local direct_IHD = r(ES)
	local direct_lci = r(ci_low)
	local direct_uci = r(ci_upp)	
	scalar direct_IHD = `direct_IHD'
	scalar indirect_diff_IHD_inc = total_IHD-direct_IHD
	scalar diff_prop_IHD_inc = (indirect_diff_IHD_inc/total_IHD)*100
	}
	}
end	
local i = 0
foreach out in CVD_inc Stroke_inc AMI_inc IHD_inc{
local i=`i'+1
local x=19
local z=19

	putexcel set RESULTS_FILE, sheet(`out') modify
	ivreg2 `out' (eduyears_scaled = ea_weighted) $maxcovar $PCs
metan or_education lci_education uci_education if outcome==`i', lcols(sample n ) effect(RD) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
bootstrap indirect = indirect_diff_`out', reps(5) nodrop : boot  // indirect effect
	matrix results = r(table)
	local indirect = results[1,1]
	local indirect_LCI = results[5,1]
	local indirect_UCI = results[6,1]

	putexcel Q`x'=`indirect' R`x'=`indirect_LCI' S`x'= `indirect_UCI'

	ivreg2 `out' (eduyears_scaled = ea_weighted) $max_covar $PCs 
metan or_education lci_education uci_education if outcome==`i', lcols(sample n ) effect(RD) nooverall null(0) xlabel(-0.0001, -0.001, -0.01, 0.0) astext(50) dp(6) nograph
	
bootstrap diff_prop_med = diff_prop_`out', reps(5) nodrop : boot
         matrix results = r(table)
         local proportion = results[1,1]
         local proportion_LCI = results[5,1]
         local proportion_UCI = results[6,1]
        
	putexcel U`z'=`proportion' V`z'=`proportion_LCI' W`z'= `proportion_UCI'

}
