/***** Question 4 : Random effects model **** /

/*Step 1: check with previously selected mean structure and unstructured covariance matrix the residuals and compare with observe variance */

/*model 3un */
proc mixed data=lda.acu2 method=reml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity= age frequency group*log_time frequency*log_time age*log_time / s;
repeated ltimeclss / type=un subject=ID r rcorr;
run;


/* Observed variance components of the response assuming parsimonous mean structure*/
data test;
input GROUP x y;
cards;
1 0 93.6006
1 1.3862943611 166.80
1 2.5649493575 160.17
run;

/*Mixed model with random intercepts for the same mean structure*/
proc mixed data=lda.acu2 method=ml;
class ltimeclss group id;
model severity= age frequency group*log_time frequency*log_time age*log_time / solution;
random intercept /type=un subject=ID g gcorr v vcorr solution;
run;

/* the fitted variance function (based on the mixed model w.r.t logTIME)- used the values of D matrix into the linear equation below*/
data test2;
	do x=0 to 2.6 by 0.5;
		group=2;
		y=67.6864*x + 65.3624;
		output;
	end;

	/* combining the fitted and observed variances in one data set */
data test;
	set test test2;
run;

/* creating the plot*/
goptions reset=all;

proc gplot data=test;
	plot y*x=group / legend=legend1 haxis=axis1 vaxis=axis2;
	symbol1 c=red value=dot h=3 i=join w=10 l=1 mode=include;
	symbol2 c=black i=join w=10 l=1 mode=include;
	title h=2.5 '    Observed and fitted variance function';
	axis1 label=(h=2 'log(TIME+1) ') value=(h=2) order=(0 to 2.6 by 0.5) offset=(2)cm 
		minor=none w=5;
	axis2 label=(h=2 angle=90 'Var(Severity)') value=(h=2) w=5 order=(60 to 240 by 
		40) offset=(0.2)cm minor=none w=5;
	legend1 down=2 frame label=(h=1.5 'Variance:') value=(h=1.5 'Observed' 
		'Fitted') order=(1 2) mode=protect position=(inside top left) offset=(0.5cm, 
		-0.5cm);
	run;


/* Even if Tried with type= CSH, still variance param estimates only slightly better, hence will try later after random slopes and shit */




/* Random slopes included  */

/*Mixed model with random intercepts and SLOPES INCLUDED for the same mean structure*/
proc mixed data=lda.acu2 method=reml;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time / solution;
random intercept log_time /type=un subject=ID g gcorr v vcorr solution;
run;
/* G matrix not positive definite: hence removing bound for negative variance components:

/* interpretation of negative components is valid only in marginal models and not when doing hierarchial interpretation */

/*model with nobound on variance components */
proc mixed data=lda.acu2 method=reml nobound;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time /  solution;
random intercept log_time  /type=un subject=ID g gcorr v vcorr solution;
run;
/*7705.9 df 11 */
/*AIC 7713.9
/* Observed variance components of the response assuming parsimonous mean structure*/
data test;
input GROUP x y;
cards;
1 0 93.6006
1 1.3862943611 166.80
1 2.5649493575 160.17
run;


/* the fitted variance function (based on the mixed model w.r.t logTIME)- used the values of D matrix into the linear equation below*/
data test3;
	do x=0 to 3 by 0.5;
		group=2;
		y=-2.2336*(x**2)+2*16.4422*x + 33.6065+69.5782;
		output;
	end;

	/* combining the fitted and observed variances in one data set */
data test;
	set test test3;
run;

/* creating the plot*/
goptions reset=all;

proc gplot data=test;
	plot y*x=group / legend=legend1 haxis=axis1 vaxis=axis2;
	symbol1 c=red value=dot h=3 i=join w=10 l=1 mode=include;
	symbol2 c=black i=join w=10 l=1 mode=include;
	title h=2.5 '    Observed and fitted variance function';
	axis1 label=(h=2 'log(TIME+1) ') value=(h=2) order=(0 to 2.6 by 0.5) offset=(2)cm 
		minor=none w=5;
	axis2 label=(h=2 angle=90 'Var(Severity)') value=(h=2) w=5 order=(60 to 240 by 
		40) offset=(0.2)cm minor=none w=5;
	legend1 down=2 frame label=(h=1.5 'Variance:') value=(h=1.5 'Observed' 
		'Fitted') order=(1 2) mode=protect position=(inside top left) offset=(0.5cm, 
		-0.5cm);
	run;


/******* Try quadratic random time effect *******/

proc mixed data=lda.acu2 method=reml nobound;
class ltimeclss group id;
model severity= age frequency group*log_time frequency*log_time age*log_time / solution;
random intercept log_time log_time*log_time  /type=un subject=ID g gcorr v vcorr solution;
run;


/* Observed variance components of the response assuming parsimonous mean structure*/
data test;
input GROUP x y;
cards;
1 0 93.6006
1 1.3862943611 166.80
1 2.5649493575 160.17
run;


/* the fitted variance function (based on the mixed model w.r.t logTIME)- used the values of D matrix into the 4th order equation below*/
data testage2;
   DO x=0 to 3 by 0.5;
   group=2;
   y= -0.3028*(x**4)+2*2.1807*(x**3) - 11.6866*x**2 - 2*12.338*(x**2)+ 2*47.8191*x + 22.7580 + 70.8432;
   output;
   end;
	
/* combining the fitted and observed variances in one data set */
data test;
	set test testage2;
run;

/* creating the plot*/
goptions reset=all;

proc gplot data=test;
	plot y*x=group / legend=legend1 haxis=axis1 vaxis=axis2;
	symbol1 c=red value=dot h=3 i=join w=10 l=1 mode=include;
	symbol2 c=black i=join w=10 l=1 mode=include;
	title h=2.5 '    Observed and fitted variance function';
	axis1 label=(h=2 'log(TIME+1) ') value=(h=2) order=(0 to 2.6 by 0.5) offset=(2)cm 
		minor=none w=5;
	axis2 label=(h=2 angle=90 'Var(Severity)') value=(h=2) w=5 order=(60 to 240 by 
		40) offset=(0.2)cm minor=none w=5;
	legend1 down=2 frame label=(h=1.5 'Variance:') value=(h=1.5 'Observed' 
		'Fitted') order=(1 2) mode=protect position=(inside top left) offset=(0.5cm, 
		-0.5cm);
	run;

/* as fitted covariance structure of sq logtime random effects explain the observed covariance structure good enough,
we now add just simple structure for the residual variance in the repeated statement */

/* we try random slopes with uncorrelated slopes and intercepts as well (banded variance type) */

proc mixed data=lda.acu2 method=reml nobound;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time /  solution;
random intercept log_time  /type=un(1) subject=ID g gcorr v vcorr solution;
repeated ltimeclss / type=simple subject=id r rcorr;
run;
/*7727.9 df 11 (10 technically) */

/* now random slope not negative anymore */
/*plot fitted variance function and compared with observed variance again */

/* Observed variance components of the response assuming parsimonous mean structure*/
data test;
input GROUP x y;
cards;
1 0 93.6006
1 1.3862943611 166.80
1 2.5649493575 160.17
run;


/* the fitted variance function (based on the mixed model w.r.t logTIME)- used the values of D matrix into the linear equation below*/
data test3;
	do x=0 to 3 by 0.5;
		group=2;
		y=5.4903*(x**2)+ 58.5019+60.2614;
		output;
	end;

	/* combining the fitted and observed variances in one data set */
data test;
	set test test3;
run;

/* creating the plot*/
goptions reset=all;

proc gplot data=test;
	plot y*x=group / legend=legend1 haxis=axis1 vaxis=axis2;
	symbol1 c=red value=dot h=3 i=join w=10 l=1 mode=include;
	symbol2 c=black i=join w=10 l=1 mode=include;
	title h=2.5 '    Observed and fitted variance function';
	axis1 label=(h=2 'log(TIME+1) ') value=(h=2) order=(0 to 2.6 by 0.5) offset=(2)cm 
		minor=none w=5;
	axis2 label=(h=2 angle=90 'Var(Severity)') value=(h=2) w=5 order=(60 to 240 by 
		40) offset=(0.2)cm minor=none w=5;
	legend1 down=2 frame label=(h=1.5 'Variance:') value=(h=1.5 'Observed' 
		'Fitted') order=(1 2) mode=protect position=(inside top left) offset=(0.5cm, 
		-0.5cm);
	run;

/*we see that when we used uncorrelated random slopes and intercepts,
the LogL is sig diff from the unstructured model, so we dont go for it */


/*********************************

/* final model */
proc mixed data=lda.acu2 method=reml nobound;
class ltimeclss group id;
model severity= age frequency group*log_time frequency*log_time age*log_time /  solution;
random intercept log_time  /type=un subject=ID g gcorr v vcorr solution;
contrast ’Equal slopes’ group*log_time 1 -1;
repeated ltimeclss / type=simple subject=id r rcorr;
run;
/* LogL 7693.3 df 14
/* comparing this model with unstructured model using LRT indicates significant change pchisq(7665.1-7690.2,21-15) = 0.00032
/* comparing this model with model 3 using LRT indicates significant change pchisq(7673.2-7690.2,14-15) = 10^-6
/* comparing this model with model 3CSH using LRT indicates significant change pchisq(7674.1-7690.2,12-15) = 0.001

