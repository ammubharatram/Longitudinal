libname LDA '/folders/myfolders';


/* Dataset preparation */

data lda.acu2;
set lda.acu2;
by id;
prev1 = lag(frequency);
prev2= 0;
if first.id then prev1 = 0;
if last.id then prev2 = lag2(frequency);
run;

/* TRANSITIONAL MODEL */

/* First order autoregressive */
proc genmod data=lda.acu2 descending;
model frequency = group chronicity age time group*time prev1 / dist=poisson;
run;

/* Second order autoregressive */
proc genmod data=lda.acu2 descending;
model frequency = group chronicity age time group*time prev1 prev2 / dist=poisson;
run;



/* CONDITIONAL MODELS */

/* Standard GEE */

proc genmod data=lda.acu2 descending;
class id timeclass;
model frequency = age chronicity group time group*time
/ dist=poisson;
repeated subject=id / withinsubject=timeclass
type=exch covb corrw modelse;
run;

/* Linearization based model */

%glimmix(data=test, procopt=%str(method=ml empirical)
,
stmts=%str(
class idnum timeclss;
model onyresp = treatn time treatn*time / solution;
repeated timeclss / subject=idnum type=cs rcorr;
),
error=binomial,
link=logit);

proc glimmix data=lda.acu2 method=RSPL empirical;
class id;
model frequency = age chronicity group time group*time
/ dist=poisson solution;
random _residual_ / subject=id type=cs;
run;


/* GENERALIZED LINEAR MIXED MODELS */

proc glimmix data=lda.acu2 method=gauss(q=50);
class id;
model frequency = age chronicity group time group*time
/ dist=poisson solution;
random intercept / subject=id;
run;

/* inclusion of random slopes? */
