*! 1.0.0 Anirban Basu Sep 22 2012
capture program drop petiv
program define petiv
	version 11.0

	syntax [varlist(min=2)] [if] [in], TRT(varname) PS(varname) CMD(string) [DEGree(integer 1) CONTRols(varlist)  DISPlay *] 

capture drop _*

gettoken (local) yvar (local) xlist : varlist	
marksample touse
tempvar _touse _mss
mark `_touse' `if' `in'
egen `_mss' = rmiss(`xlist' `yvar' `ps' `controls')

capture inspect `yvar' if (`_touse' & `_mss'==0)

if (("`cmd'"=="probit") | ("`cmd'"=="logit")) {
	if (r(N_unique) !=2) {
		noi di in err "Outcome variable is not binary"
		exit 198
	}
}
else if  (("`cmd'"=="glm") | ("`cmd'"=="pglm") | ("`cmd'"=="xtgee")) {
	if (r(N_neg) >0) {
		noi di in err "Outcome variable contains negative values"
		exit 198
	}
}
else if ("`cmd'" !="reg") {
		noi di in err "Regression command not recognized"
		exit 198
}		


if (`degree' > 4) {
	noi di in err "polynomial order greater that 4 not allowed. 4 is assumed instead."
	local degree = 4
}


capture inspect `trt' if (`_touse' & `_mss'==0)
	if (r(N_unique) !=2) {
		noi di in err "Treatment in not binary"
		exit 198
	}
	

/* Calculations of interactions and polynomials */

global xlisti " "
global j=0
global xlst "`xlist'"

qui foreach var of global xlst {
tempname p_`var'
gen `p_`var'' = `ps'*`var'
global xlisti "$xlisti  `p_`var''"
global j = $j +1
}



global k=0
global xlst2 "`xlist' `controls'"
qui foreach var of global xlst2 {
global k = $k +1
}


tempname ps2 ps3 ps4
gen `ps2' = `ps'^2
gen `ps3' = `ps'^3
gen `ps4' = `ps'^4


tempvar  ysc
qui summ `yvar' if (`_touse' & `_mss'==0), meanonly
gen `ysc' = `yvar'/r(mean)
global sc_=r(mean) 

local disp="qui"
if "`display'"!="" {
 local disp="noi"
} 


/************* BINARY OUTCOMES *******************/
/* PROBIT MODEL */
qui if (("`cmd'"=="probit")|("`cmd'"=="logit"))  {

	if (`degree'==1) {
	if ("`cmd'"=="probit") {
		`disp' probit `yvar' $xlisti $xlst2 `ps' if (`_touse' & `_mss'==0),  `options' 
	}
	else if ("`cmd'"=="logit") {
		`disp' logit `yvar' $xlisti $xlst2 `ps' if (`_touse' & `_mss'==0),  `options' 
	}

	global bps=_b[`ps']
	capture global bps2=0
	capture global bps3=0
	capture global bps4 =0
	}
	
	else if (`degree'==2) {
	if ("`cmd'"=="probit") {
		`disp' probit `yvar' $xlisti $xlst2 `ps' `ps2' if (`_touse' & `_mss'==0),  `options' 
	}
	else if ("`cmd'"=="logit") {
		`disp' logit `yvar' $xlisti $xlst2 `ps' `ps2' if (`_touse' & `_mss'==0),  `options' 
	}

	global bps=_b[`ps']
	capture global bps2=_b[`ps2']
	capture global bps3=0
	capture global bps4 =0
	}
	
	else if (`degree'==3) {
	if ("`cmd'"=="probit") {
		`disp' probit `yvar' $xlisti $xlst2 `ps' `ps2' `ps3' if (`_touse' & `_mss'==0),  `options' 
	}
	else if ("`cmd'"=="logit") {
		`disp' logit `yvar' $xlisti $xlst2 `ps' `ps2' `ps3' if (`_touse' & `_mss'==0),  `options' 
	}

	global bps=_b[`ps']
	capture global bps2=_b[`ps2']
	capture global bps3=_b[`ps3']
	capture global bps4 =0
	}
	
	else if (`degree'==4) {
	if ("`cmd'"=="probit") {
		`disp' probit `yvar' $xlisti $xlst2 `ps' `ps2' `ps3' `ps4' if (`_touse' & `_mss'==0),  `options' 
	}
	else if ("`cmd'"=="logit") {
		`disp' logit `yvar' $xlisti $xlst2 `ps' `ps2' `ps3' `ps4' if (`_touse' & `_mss'==0),  `options' 
	}	 

	global bps=_b[`ps']
	capture global bps2=_b[`ps2']
	capture global bps3=_b[`ps3']
	capture global bps4 =_b[`ps4']
	}


mat betas=e(b)
scalar c1=colsof(betas)

mat beta1_beta0=betas[1,1..$j]    /* Extract coefficients of the interactions only */ 
qui foreach var of global xlst {		/* Replace by interactions with Ps=1 */
replace `p_`var'' = `var'
}

tempname x_beta1_beta0
mat score `x_beta1_beta0'=beta1_beta0

mat beta2=betas[1,$j+1..$j+$k]    /* Extract coefficients of the main-effects only */ 
tempname maineff
mat score `maineff'=beta2
replace `maineff' = `maineff' + _b[_cons]

qui foreach var of global xlst {		/* Replace back interactions */
replace `p_`var'' = `ps'*`var'
}

tempname xb
gen `xb' =1

sort `xb'
mata: xbeta = st_data(.,("`x_beta1_beta0'"))
mata: ps = st_data(.,("`ps'"))
mata: xb =st_data(.,("`xb'"))
mata: y = st_data(.,("`yvar'"))
mata: d = st_data(.,("`trt'"))
mata: maineff = st_data(.,("`maineff'"))


	if (("`cmd'"=="probit") {
		mata: te = mtebin(xbeta, xb, ps, maineff, $bps, $bps2, $bps3, $bps4, y, d, 1, 1)
	}
	else if (("`cmd'"=="logit") {
		mata: te = mtebin(xbeta, xb, ps, maineff, $bps, $bps2, $bps3, $bps4, y, d, 1, 2)
	}
	
capture drop pet_`yvar'
mata: st_addvar(("float"), ("pet_`yvar'"))
mata: st_store(., ("pet_`yvar'"), te)
replace pet_`yvar' =. if (`_mss'!=0)
} /* end of logit or probit analyses */


/*********************** NON-NEGATIVE OUTCOMES *******************************/



qui if (("`cmd'"=="glm")|("`cmd'"=="pglm") | ("`cmd'"=="xtgee")) {

	if (`degree'==1) {
	if ("`cmd'"=="glm") {
		`disp' glm `yvar' $xlisti $xlst2 `ps' if (`_touse' & `_mss'==0),  `options' 
		capture global pwr_ = real(e(power))
		capture global dnm_ = real(e(denom))
	}
	if ("`cmd'"=="xtgee") {
		`disp' xtgee `yvar' $xlisti $xlst2 `ps' if (`_touse' & `_mss'==0),  `options' 
		capture global pwr_ = real(e(power))
		capture global dnm_ = real(e(denom))
	}
	else if ("`cmd'"=="pglm") {
		 if "`options'"=="" {
		 `disp' pglm `ysc' $xlisti $xlst2 `ps' if (`_touse' & `_mss'==0) 
		 }
		 else {
		 `disp' pglm `ysc' $xlisti $xlst2 `ps' if (`_touse' & `_mss'==0),  `options' 		 
		 }
		 capture global lambda = _b[lambda:_cons]
	}
	global bps=_b[`ps']
	capture global bps2=0
	capture global bps3=0
	capture global bps4 =0
	}
	
	else if (`degree'==2) {
	if ("`cmd'"=="glm") {
		`disp' glm `yvar' $xlisti $xlst2 `ps' `ps2' if (`_touse' & `_mss'==0),  `options'  link(`link')
		capture global pwr_ = real(e(power))
		capture global dnm_ = real(e(denom))
	}
	if ("`cmd'"=="xtgee") {
		`disp' xtgee `yvar' $xlisti $xlst2 `ps' `ps2' if (`_touse' & `_mss'==0),  `options' 
		capture global pwr_ = real(e(power))
		capture global dnm_ = real(e(denom))
	}
	else if ("`cmd'"=="pglm") {
		 if "`options'"=="" {
		 `disp' pglm `ysc' $xlisti $xlst2 `ps' `ps2' if (`_touse' & `_mss'==0) 
		 }
		 else {
		 `disp' pglm `ysc' $xlisti $xlst2 `ps' `ps2' if (`_touse' & `_mss'==0),  `options' 		 
		 }
		 capture global lambda = _b[lambda:_cons]
	}
	global bps=_b[`ps']
	capture global bps2=_b[`ps2']
	capture global bps3=0
	capture global bps4 =0
	}
	
	else if (`degree'==3) {
	if ("`cmd'"=="glm") {
		`disp' glm `yvar' $xlisti $xlst2 `ps' `ps2' `ps3' if (`_touse' & `_mss'==0),  `options'  link(`link')
		capture global pwr_ = real(e(power))
		capture global dnm_ = real(e(denom))
	}
	if ("`cmd'"=="xtgee") {
		`disp' xtgee `yvar' $xlisti $xlst2 `ps' `ps2' `ps3' if (`_touse' & `_mss'==0),  `options' 
		capture global pwr_ = e(power)
		capture global dnm_ = e(denom)
	}
	else if ("`cmd'"=="pglm") {
		 if "`options'"=="" {
		 `disp' pglm `ysc' $xlisti $xlst2 `ps' `ps2' `ps3' if (`_touse' & `_mss'==0) 
		 }
		 else {
		 `disp' pglm `ysc' $xlisti $xlst2 `ps' `ps2' `ps3' if (`_touse' & `_mss'==0),  `options'  	 
		 }
		 capture global lambda = _b[lambda:_cons]
	}
	global bps=_b[`ps']
	capture global bps2=_b[`ps2']
	capture global bps3=_b[`ps3']
	capture global bps4 =0
	}
	
	else if (`degree'==4) {
	if ("`cmd'"=="glm") {
		`disp' glm `yvar' $xlisti $xlst2 `ps' `ps2' `ps3' `ps4' if (`_touse' & `_mss'==0),  `options'  link(`link')
		capture global pwr_ = real(e(power))
		capture global dnm_ = real(e(denom))
	}
	if ("`cmd'"=="xtgee") {
		`disp' xtgee `yvar' $xlisti $xlst2 `ps' `ps2' `ps3' `ps4' if (`_touse' & `_mss'==0),  `options' 
		capture global pwr_ = real(e(power))
		capture global dnm_ = real(e(denom))
	}
	else if ("`cmd'"=="pglm") {
		 if "`options'"=="" {
		 `disp' pglm `ysc' $xlisti $xlst2 `ps' `ps2' `ps3' `ps4' if (`_touse' & `_mss'==0) 
		 }
		 else {
		 `disp' pglm `ysc' $xlisti $xlst2 `ps' `ps2' `ps3' `ps4' if (`_touse' & `_mss'==0),  `options' 		 
		 }
		 capture global lambda = _b[lambda:_cons]
	} 
	global bps=_b[`ps']
	capture global bps2=_b[`ps2']
	capture global bps3=_b[`ps3']
	capture global bps4 =_b[`ps4']
	}

mat betas=e(b)
scalar c1=colsof(betas)
	
mat beta1_beta0=betas[1,1..$j]    /* Extract coefficients of the interactions only */ 
qui foreach var of global xlst {		/* Replace by interactions with Ps=1 */
replace `p_`var'' = `var'
}

tempname x_beta1_beta0
mat score `x_beta1_beta0'=beta1_beta0

mat beta2=betas[1,$j+1..$j+$k]    /* Extract coefficients of the main-effects only */ 
tempname maineff
mat score `maineff'=beta2
replace `maineff' = `maineff' + _b[_cons]

qui foreach var of global xlst {		/* Replace back interactions */
replace `p_`var'' = `ps'*`var'
}

tempname xb
gen `xb' =1

sort `xb'
mata: xbeta = st_data(.,("`x_beta1_beta0'"))
mata: ps = st_data(.,("`ps'"))
mata: xb =st_data(.,("`xb'"))
mata: y = st_data(.,("`yvar'"))
mata: d = st_data(.,("`trt'"))
mata: maineff = st_data(.,("`maineff'"))


	if ("`cmd'"=="glm")  {
		if (e(linkt)=="Probit") {
			mata: te = mtebin(xbeta, xb, ps, maineff, $bps, $bps2, $bps3, $bps4, y, d, 1, 1)
			mata: te = te:*$dnm_
			}
		else  if  (e(linkt)=="Logit") {
				mata: te = mtebin(xbeta, xb, ps, maineff, $bps, $bps2, $bps3, $bps4, y, d, 1, 2)
			mata: te = te:*$dnm_
			}
		else  {
			mata: te = mteglm(xbeta, xb, ps, maineff, $pwr_, $bps, $bps2, $bps3, $bps4, y, d)
		}
	}
	else if ("`cmd'"=="xtgee") {
		if (e(link)=="probit") {
			mata: te = mtebin(xbeta, xb, ps, maineff, $bps, $bps2, $bps3, $bps4, y, d,1, 1)
			mata: te = te:*$dnm_
		}
		else if  (e(link)=="logit") {
			mata: te = mtebin(xbeta, xb, ps, maineff, $bps, $bps2, $bps3, $bps4, y, d,1, 2)	
			mata: te = te:*$dnm_	
		}
			else if  (substr(e(link),1,4)=="odds") {
			mata: te = mtebin(xbeta, xb, ps, maineff, $bps, $bps2, $bps3, $bps4, y, d, $pwr_, 4)	
			mata: te = te:*$dnm_	
		}
	}
	else if ("`cmd'"=="pglm") {
		mata: te = mtepglm(xbeta, xb, ps, maineff, $lambda, $sc_, $bps, $bps2, $bps3, $bps4, y, d)
	}
	
capture drop pet_`yvar'
mata: st_addvar(("float"), ("pet_`yvar'"))
mata: st_store(., ("pet_`yvar'"), te)
replace pet_`yvar' =. if (`_mss'!=0)
} /* end of glm or pglm or xtgee analyses */


/************* GENERAL OUTCOMES for LINEAR REGRESSION *******************/
/* LINEAR MODEL */
qui if ("`cmd'"=="reg") {

	if (`degree'==1) {
		`disp' reg `yvar' $xlisti $xlst2 `ps' if (`_touse' & `_mss'==0),  `options' 
	global bps=_b[`ps']
	capture global bps2=0
	capture global bps3=0
	capture global bps4 =0
	}
	
	else if (`degree'==2) {
		`disp' reg `yvar' $xlisti $xlst2 `ps' `ps2' if (`_touse' & `_mss'==0),  `options' 
	global bps=_b[`ps']
	capture global bps2=_b[`ps2']
	capture global bps3=0
	capture global bps4 =0
	}
	
	else if (`degree'==3) {
		`disp' reg `yvar' $xlisti $xlst2 `ps' `ps2' `ps3' if (`_touse' & `_mss'==0),  `options' 

	global bps=_b[`ps']
	capture global bps2=_b[`ps2']
	capture global bps3=_b[`ps3']
	capture global bps4 =0
	}
	
	else if (`degree'==4) {
		`disp' reg `yvar' $xlisti $xlst2 `ps' `ps2' `ps3' `ps4' if (`_touse' & `_mss'==0),  `options' 

	global bps=_b[`ps']
	capture global bps2=_b[`ps2']
	capture global bps3=_b[`ps3']
	capture global bps4 =_b[`ps4']
	}


mat betas=e(b)
scalar c1=colsof(betas)

mat beta1_beta0=betas[1,1..$j]    /* Extract coefficients of the interactions only */ 
qui foreach var of global xlst {		/* Replace by interactions with Ps=1 */
replace `p_`var'' = `var'
}

tempname x_beta1_beta0
mat score `x_beta1_beta0'=beta1_beta0

qui foreach var of global xlst {		/* Replace back interactions */
replace `p_`var'' = `ps'*`var'
}

tempname xb
gen `xb' =1

sort `xb'
mata: xbeta = st_data(.,("`x_beta1_beta0'"))
mata: ps = st_data(.,("`ps'"))
mata: xb =st_data(.,("`xb'"))
mata: y = st_data(.,("`yvar'"))
mata: d = st_data(.,("`trt'"))
mata: te = mtereg(xbeta, xb, ps, $bps, $bps2, $bps3, $bps4, y, d)

capture drop pet_`yvar'
mata: st_addvar(("float"), ("pet_`yvar'"))
mata: st_store(., ("pet_`yvar'"), te)
replace pet_`yvar' =. if (`_mss'!=0)
} /* end of regress analyses */


end


