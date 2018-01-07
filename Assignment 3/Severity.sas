libname LDA 'C:\Users\u0117439\Desktop\LDA';

data lda.acu2;
set lda.aculda2;
drop frequency;
timeclass=time;
run;

/* Long to wide format*/
proc transpose data=lda.acu2 out=lda.acu2b prefix=severity;
    by id age chronicity group;
    id time;
    var severity ;
run;


/* Exploration of missing data mechanisms.*/

proc sort data=lda.acu2b;
by group;
run;

/*Overview of missigness table*/
proc mi data=lda.acu2b seed=675938 simple nimpute=0;
title ’standard EM’;
em itprint outem=growthem1;
var severity3 severity12;
run;
/*Overview of missigness table per group*/
proc mi data=lda.acu2b seed=675938 simple nimpute=0;
title ’standard EM’;
em itprint outem=growthem1;
var severity3 severity12;
by group;
run;



/* Using data as it is.*/

/*LMM*/
/*method=ml*/
proc mixed data=lda.acu2 method=ml;
class timeclass group id;
model severity= age chronicity group*time age*time chronicity*time/ solution;
random intercept time  /type=un subject=ID g gcorr v vcorr;
repeated timeclass / type=un subject=ID r rcorr;
run;
/*method=reml*/
proc mixed data=lda.acu2 method=reml;
class timeclass group id;
model severity= age chronicity group*time age*time chronicity*time/ solution;
random intercept time  /type=un subject=ID g gcorr v vcorr;
repeated timeclass / type=un subject=ID r rcorr;
run;

/*GEE*/
proc genmod data=lda.acu2;
title ’data as is - GEE’;
class timeclass group id;
model severity=age chronicity time group*time age*time chronicity*time /dist=normal;
repeated subject=id / withinsubject=timeclass type=un modelse;
run;
proc glimmix data=lda.acu2;
title ’data as is - GEE - linearized version’;
nloptions maxiter=50 technique=newrap;
class timeclass group id;
model severity=age chronicity time group*time age*time chronicity*time /dist=normal;
random _residual_ / subject=id type=un;
run;
proc glimmix data=lda.acu2 empirical;
title ’data as is - GEE - linearized version - empirical’;
nloptions maxiter=50 technique=newrap;
class timeclass group id;
model severity=age chronicity time group*time age*time chronicity*time /dist=normal;
random _residual_ / subject=id type=un;
run;


/*WGEE*/

/* WGEE: macro for creating variables "dropout" and "prev" */

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


%dropout(data=lda.acu2,id=id,time=time,response=severity,out=lda.acu4);

proc genmod data=lda.acu4 descending;
class group prev timeclass;
model dropout = prev group age chronicity time / pred dist=b;
ods output obstats=pred;
ods listing exclude obstats;
run;

proc print data=pred;
run;

data pred;
set pred;
keep observation pred;
run;
data lda.acu4b;
merge pred lda.acu4;
run;
proc print data=lda.acu4b;
run;



%macro dropwgt(data=,id=,time=,pred=,dropout=,out=);
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
read all var {&id &time &pred &dropout};
n = nrow(&pred);
wi = j(n,1,1);
ind = 1;
do while (ind <= nsub);
wihlp = 1;
stay = 1;
/* first measurement */
if (&dropout[(ind-1)*ntime+2]=1)
then do;
wihlp = pred[(ind-1)*ntime+2];
stay = 0;
end;
else if (&dropout[(ind-1)*ntime+2]=0)
then wihlp = 1-pred[(ind-1)*ntime+2];
/* second to penultimate measurement */
j=2;
do while ((j <= ntime-1) & stay);
if (&dropout[(ind-1)*ntime+j+1]=1)
then do;
wihlp = wihlp*pred[(ind-1)*ntime+j+1];
stay = 0;
end;
else if (&dropout[(ind-1)*ntime+j+1]=0)
then wihlp = wihlp*(1-pred[(ind-1)*ntime+j+1]);
j = j+1;
end;
j = 1;
do while (j <= ntime);
wi[(ind-1)*ntime+j]=wihlp;
j = j+1;
end;
ind = ind+1;
end;
create help var {&id &time &pred &dropout wi};
append;
quit;
data &out;
merge &data help;
data &out;
set &out;
wi=1/wi;
run;
%mend;

%dropwgt(data=lda.acu4b,id=id,time=time,pred=pred,dropout=dropout,out=lda.acu5);
proc print data=lda.acu5;
var id time severity dropout prev pred wi;
run;

/* WGEE: models */

proc genmod data=lda.acu5;
title ’data as is - WGEE’;
scwgt wi;
class timeclass group id;
model severity=age chronicity time group*time age*time chronicity*time /dist=normal;
repeated subject=id / withinsubject=timeclass type=un modelse;
run;
proc glimmix data=lda.acu5;
title ’data as is - WGEE - linearized version’;
nloptions maxiter=50 technique=newrap;
weight wi;
class timeclass group id;
model severity=age chronicity time group*time age*time chronicity*time /dist=normal;
random _residual_ / subject=id type=un;
run;
proc glimmix data=lda.acu4 empirical;
title ’data as is - WGEE - linearized version - empirical’;
weight wi;
nloptions maxiter=50 technique=newrap;
class timeclass group id;
model severity=age chronicity time group*time age*time chronicity*time /dist=normal;
random _residual_ / subject=id type=un;
run;



/*Multiple imputation*/

/*1. MI: imputation task (MCMC method)*/
proc sort data=lda.acu2b;
by group;
run;

proc mi data=lda.acu2b seed=486048 simple out=lda.acu6 nimpute=10 round=0.1;
var age chronicity severity0 severity3 severity12;
by group;
run;

data lda.acu6b;
set lda.acu6;
array x (3) severity0 severity3 severity12;
do j=1 to 3;
severity=x(j);
time=j;
output;
end;
run;
data lda.acu6b;
set lda.acu6b;
drop severity0 severity3 severity12 j;
timeclass=time;
run;

proc sort data=lda.acu6b;
by _imputation_ id time;
run;

/*2.1. MI: analysis task: GEE*/
proc genmod data=lda.acu6b;
title ’GEE after multiple imputation’;
class timeclass group id;
by _imputation_;
model severity=age chronicity time group*time age*time chronicity*time /dist=normal covb;
repeated subject=id / withinsubject=timeclass type=un modelse;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo CovB=gmcovb;
run;

/*MI: delete redundant parameters*/
data gmpinfo;
set gmpinfo;
if parameter="Prm6" then delete;
run;
data gmparms;
set gmparms;
if level1=1 then delete;
run;

/*2.2. MI: analysis task: LMM*/
proc mixed data=lda.acu6b method=ml asycov covtest;
class timeclass group id;
by _imputation_;
model severity= age chronicity time group*time age*time chronicity*time/ s covb;
random intercept time  /type=un subject=ID g gcorr v vcorr solution;
repeated timeclass / type=un subject=ID r rcorr;
ods output solutionf = lmsolution covb = lmcovb covparms = lmcovparms asycov = lmasycov;
run;

data lmsolution;
set lmsolution;
if group=1 then delete;
data lmcovb;
set lmcovb;
if row=6 then delete;
data lmcovparms;
set lmcovparms;
if CovParm=’YEAR UN(1,1)’ then effect = ’YEARUN11’;
if CovParm=’ UN(2,1)’ then effect = ’YEARUN21’;
if CovParm=’ UN(2,2)’ then effect = ’YEARUN22’;
if CovParm=’ UN(3,1)’ then effect = ’YEARUN31’;
if CovParm=’ UN(3,2)’ then effect = ’YEARUN32’;
if CovParm=’ UN(3,3)’ then effect = ’YEARUN33’;
if CovParm=’PARENT Corr’ then effect = ’PARENTCORR’ ;
drop covparm;
data lmasycov;
set lmasycov;
Col1=CovP1;
Col2=CovP2;
Col3=CovP3;
Col4=CovP4;
Col5=CovP5;
Col6=CovP6;
Col7=CovP7;
if CovParm=’YEAR UN(1,1)’ then effect = ’YEARUN11’;
if CovParm=’ UN(2,1)’ then effect = ’YEARUN21’;
if CovParm=’ UN(2,2)’ then effect = ’YEARUN22’;
if CovParm=’ UN(3,1)’ then effect = ’YEARUN31’;
if CovParm=’ UN(3,2)’ then effect = ’YEARUN32’;
if CovParm=’ UN(3,3)’ then effect = ’YEARUN33’;
if CovParm=’PARENT Corr’ then effect = ’PARENTCORR’ ;
drop CovP1 CovP2 CovP3 CovP4 CovP5 CovP6 CovP7 covparm;
run;



/*3.1. MI: inference task: GEE*/
proc mianalyze parms=gmparms covb=gmcovb parminfo=gmpinfo wcov bcov tcov;
modeleffects age chronicity time group*time age*time chronicity*time;
run;

/*3.2. MI: inference task: LMM*/
/* Combining 10 Separate Analyses (mean structure) */
proc mianalyze parms=lmsolution covb(effectvar=rowcol)=lmcovb;
title2 ’COMBINING 5 MIXED MODEL ANALYSES (MEAN STRUCTURE)’;
modeleffects age chronicity time group*time age*time chronicity*time;
run;
/* Combining 10 Separate Analyses (covariance structure) */
proc mianalyze parms=covparms0 covb(effectvar=rowcol)=asycov0;
title2 ’COMBINING 5 MIXED MODEL ANALYSES (COVARIANCE STRUCTURE)’;
modeleffects YEARUN11 YEARUN21 YEARUN22 YEARUN31 YEARUN32 YEARUN33 PARENTCO;
run;





