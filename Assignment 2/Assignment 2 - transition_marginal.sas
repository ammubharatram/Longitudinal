libname LDA '/folders/myfolders';


/* TRANSITIONAL MODEL */

/* Dataset preparation */

%macro dropout(data=,id=,time=,response=,out=);
%if %bquote(&data)= %then %let data=&syslast;
proc freq data=&data noprint;
tables &id /out=freqid;
tables &time / out=freqtime;
run;
proc iml;
reset noprint;
use freqid;
read all var {&id};
nsub = nrow(&id);
use freqtime;
read all var {&time};
ntime = nrow(&time);
time = &time;
use &data;
read all var {&id &time &response};
n = nrow(&response);
dropout = j(n,1,0);
ind = 1;
do while (ind <= nsub);
j=1;
if (&response[(ind-1)*ntime+j]=.) then print "First Measurement is Missing";
if (&response[(ind-1)*ntime+j]^=.) then
do;
j = ntime;
do until (j=1);
if (&response[(ind-1)*ntime+j]=.) then
do;
dropout[(ind-1)*ntime+j]=1;
j = j-1;
end;
else j = 1;
end;
end;
ind = ind+1;
end;
prev = j(n,1,1);
prev[2:n] = &response[1:n-1];
i=1;
do while (i<=n);
if &time[i]=time[1] then prev[i]=.;
i = i+1;
end;
create help var {&id &time &response dropout prev};
append;
quit;
data &out;
merge &data help;
run;
%mend;

/*One lagged observation*/
%dropout(data=lda.acu2,id=id,time=time,response=frequency,out=lda.acu2a);
data lda.acu2a;
set lda.acu2a;
prev1=prev;
drop prev;
run;

/*Two lagged observations*/
%dropout(data=lda.acu2a,id=id,time=time,response=prev1,out=lda.acu2b);
data lda.acu2b;
set lda.acu2b;
prev2=prev;
drop prev dropout;
run;

/* First order autoregressive non-stationary */
data lda.acu2c;
set lda.acu2b;
prev1a=prev1;
if time>3 then prev1a=0;
run;
proc genmod data=lda.acu2c descending;
model frequency = chronicity age time group*time chronicity*time age*time prev1 prev1a / dist=poisson;
run;

/* First order autoregressive */
proc genmod data=lda.acu2b descending;
model frequency = chronicity age time group*time chronicity*time age*time prev1 / dist=poisson;
run;
/* Reduction of mean structure */
proc genmod data=lda.acu2b descending;
model frequency = chronicity age time group*time age*time prev1 / dist=poisson;
run;

proc genmod data=lda.acu2b descending;
model frequency = chronicity time group*time chronicity*time prev1 / dist=poisson;
run;

proc genmod data=lda.acu2b descending;
model frequency = chronicity time group*time prev1 / dist=negbin;
run;

/* Second order autoregressive */
proc genmod data=lda.acu2b descending;
model frequency = chronicity time group*time prev1 prev2 / dist=poisson;
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


/* Alternating logistic regression */

proc genmod data=lda.acu2 descending;
class id timeclass;
model frequency = age chronicity group time group*time
/ dist=poisson;
repeated subject=id / withinsubject=timeclass
logor=exch covb corrw modelse;
run;



/* GENERALIZED LINEAR MIXED MODELS */

proc glimmix data=lda.acu2 method=gauss(q=50);
class id;
model frequency = age chronicity group time group*time
/ dist=poisson solution;
random intercept / subject=id;
run;

/* inclusion of random slopes? */
