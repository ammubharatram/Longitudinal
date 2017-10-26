/* import and copy dataset */
libname lda '/home/bharatramammu0/Longitudinal';
data lda.acu;
set lda.aculda;
run;
/*correlation structure */

data lda.acu2;
set lda.acu;
log_time = log(time+1);
run;


ods graphics on; 
title ’data:Correlation structure’;
proc corr data= lda.aculda2 plots(MAXPOINTS=NONE)=matrix(histogram);
var s0 s3 s12;
run;
/* the variability increases s0 to s3 and then decreases more from s3 to s12 */

ods graphics off;
ods graphics on;
title’data:Correlationstructure’;
proc corr data=lda.aculda (MAXPOINTS=NONE)=matrix(histogram);
var r0 r06 r1 r2 r3 r4 r5 r6 r7 r8 r9 r10;
run;
ods graphics off;


/************ Question 2 ***********/
/*Model 0: unstructured mean and covariance matrix- considering log_time as a categorical variable*/
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

/* equal intercepts= no group as main effect in model statement */
/*Model 2b*/
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits2;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time/ noint s;
repeated ltimeclss / type=un subject=ID r rcorr;
run;

/* LogL: 7688.1 no of parameters =12
/* this model is sig diff from unstr as well as model 2, so not consider */


/* Even though we saw the lines weren't parallel, we still try to reduce using this
/* try equal slopes for mean structure reduction on the optimum covariates model */
/*model 3 */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity= group age frequency log_time / noint s;
repeated ltimeclss / type=un subject=ID r rcorr;
run;
/* LogL 7721.1 df 11
/*  model 3 ( with pchisq(47.9,3 df) ) was significantly different from model 2
Hence we deny the possibility of removing interactions */


/** Covariance structures--CS,TOEP, CSH,TOEPH **/
/*model 2C */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity= group age frequency  group*log_time frequency*log_time age*log_time / noint s;
repeated ltimeclss / type=CS subject=ID r rcorr;
run;
/* 7720.2 df 10
/*model 2T */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity= group age frequency group*log_time frequency*log_time age*log_time / noint s;
repeated ltimeclss / type=TOEP subject=ID r rcorr;
run;
/* 7720.2 df 11 */
/*model 2CSH */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity= group age frequency group*log_time frequency*log_time age*log_time / noint s;
repeated ltimeclss / type=CSH subject=ID r rcorr;
run;
/* 7674.1 df 12 */
/* 2 CSH seems most appropriate as per LR test */

/*model 2TOEPH */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity= group age frequency group*log_time frequency*log_time age*log_time / noint s;
repeated ltimeclss / type=TOEPH subject=ID r rcorr;
run;
/* 7674 df 13 */





/*2.1 To obtain mean structure and variance structure and parsimonize them . 
*/
data lda.acu2;
set lda.acu;
ltimeclss=log_time;
run;



data lda.acu2;
set lda.acu;
timeclss=time;
run;

proc mixed data=lda.acu2 method=ml;
class timeclss group id;
model severity= group*time/ noint s;
repeated timeclss / type=un subject=ID r rcorr;
run;





proc mixed data=lda.acu2 method=ml;
class timeclss ID;
model severity= age time*age frequency time*frequency chronicity time*chronicity group time*group;
run;
/* Between subject variability is higher than within subject (residual) variability */

/* step 1 : full mean structure as chronicity is not significant, remove chronicity and remodel */
proc mixed data=lda.acu2 method=ml;
class timeclss ID;
model severity= time age frequency group group*time;
repeated timeclss / type=un subject=ID r rcorr;
run;
/* As of now, the above mean structure seems parsimonous */


/* step 2: as we see variance components close enough, toep structure could be worth trying */

proc mixed data=lda.acu2;
class timeclss ID;
model severity= time chronicity group group*time;
repeated timeclss / type=toep subject=ID r rcorr;
run;



















/* Random effects */
data test3;
input GROUP x y;

cards;

1 1 260.94
1 3 277.93
1 12 247.65
run;

/*Mixed model initial for fitted variances*/
proc mixed data=lda.acu;
class ID;
model severity= time age chronicity group group*time /solution;
random intercept /type=toep subject=ID g gcorr v vcorr;
run;


/* lets see if topelitz variance structure does better*/

proc mixed data=lda.acu;
class ID;
model severity= time age chronicity group group*time /solution;
random intercept/type=toep subject=ID g gcorr v vcorr;
run;



/*model 3  equal intercepts */

/*model 3C */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity= age frequency  group*log_time frequency*log_time age*log_time / s;
repeated ltimeclss / type=CS subject=ID r rcorr;
run;
/* 7721.6 df 9
/*model 3T */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity= age frequency group*log_time frequency*log_time age*log_time /  s;
repeated ltimeclss / type=TOEP subject=ID r rcorr;
run;
/* 7721.6 df 10 */
/* 1 increase in df for toep for same logL as CS , so chose CS*/
/*model 3CSH */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity=age frequency group*log_time frequency*log_time age*log_time /  s;
repeated ltimeclss / type=CSH subject=ID r rcorr;
run;
/* 7675.3 df 11 */
/* 2 CSH seems most appropriate as per LR test */

/*model 3TOEPH */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity= age frequency group*log_time frequency*log_time age*log_time / s;
repeated ltimeclss / type=TOEPH subject=ID r rcorr;
run;
/* 7675.2 df 12 */



