/************ Question 2 ***********/
/*Model 0.0: unstructured mean and covariance matrix- considering log_time as a categorical variable*/
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits0;
class log_time group id;
model severity= group*log_time 
	 / noint s;
repeated log_time / type=un subject=id r rcorr;
run;


/*Model 0: unstructured mean and covariance matrix- considering log_time as a categorical variable*/
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits0;
class log_time group id;
model severity= group*log_time chronicity*log_time frequency*log_time age*log_time
	 / noint s;
repeated log_time / type=un subject=id r rcorr;
run;
/* 7665.1 df 21
/* Two observations:
1. Chronicity*age not significant - so delete chronicity*log_time from model
2. Treating time as categorical variable may not be required and can be treated as continuous variable
*/
data LDA.acu2;
set lda.acu2;
ltimeclss = log_time;
run;
/*Model 1*/
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits1;
class ltimeclss group id;
model severity= group age chronicity frequency group*log_time frequency*log_time age*log_time/ noint s;
repeated ltimeclss / type=un subject=ID r rcorr;
run;
/* LogL 7672.2 df 15
/* chronicity not significant, remove


/*Model 2*/
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits2;
class ltimeclss group id;
model severity= group age frequency group*log_time frequency*log_time age*log_time/ noint s;
repeated ltimeclss / type=un subject=ID r rcorr;
run;

/* LogL: 7673.2 no of parameters =14



/*Model 3*/
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time/  s;
repeated ltimeclss / type=un subject=ID r rcorr;
run;

/* LogL: 7674.5 no of parameters =13
/* this model is not sig diff from unstr as well as model 2, qnd less df SO PREFERABLE*/


/* Even though we saw the lines weren't parallel, we still try to reduce using this
/* try equal slopes for mean structure reduction on the optimum covariates model */

/*model 4 */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits4;
class ltimeclss group id;
model severity= group age frequency log_time / noint s;
repeated ltimeclss / type=un subject=ID r rcorr;
run;
/* LogL 7721.1 df 11
/*  model 3 ( with pchisq(47.9,3 df) ) was significantly different from model 2
Hence we deny the possibility of removing interactions */



/** Covariance structures--CS,TOEP, CSH,TOEPH **/
/*model 3C */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3C;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time/  s;
repeated ltimeclss / type=CS subject=ID r rcorr;
run;
/* 7721.6 df 9
/*model 3T */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3T;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time/  s;
repeated ltimeclss / type=TOEP subject=ID r rcorr;
run;
/* 7721.6 df 10 */
/*model 3CSH */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3CSH;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time/  s;
repeated ltimeclss / type=CSH subject=ID r rcorr;
run;
/* 7675.3 df 11 */
/* 2 CSH seems most appropriate as per LR test */

/*model 3TOEPH */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3TOEPH;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time/  s;
repeated ltimeclss / type=TOEPH subject=ID r rcorr;
run;
/* 7675.2 df 12 */








