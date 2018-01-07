libname LDA 'C:\Users\u0117439\Desktop\LDA';

data lda.acu2;
set lda.aculda2;
drop frequency;
logtime=log(time+1);
logtimeclass=logtime;
run;


/********************************************/
/* Exploration of missing data mechanisms.*/
/********************************************/

/* Long to wide format*/
proc transpose data=lda.acu2 out=lda.acu2b prefix=severity;
    by id age chronicity group;
    id logtime;
    var severity ;
run;

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


/********************************************/
/* Using data as it is: LMM and GEE*/
/********************************************/

/*LMM*/
/*method=ml*/
proc mixed data=lda.acu2 method=ml nobound;
class logtimeclass group id;
model severity= age chronicity group*logtime age*logtime chronicity*logtime/ solution;
random intercept logtime  /type=un subject=ID g gcorr v vcorr;
repeated logtimeclass / type=un subject=ID r rcorr;
run;
/*method=reml*/
proc mixed data=lda.acu2 method=reml nobound;
class logtimeclass group id;
model severity= age chronicity group*logtime age*logtime chronicity*logtime/ solution;
random intercept logtime  /type=un subject=ID g gcorr v vcorr;
repeated logtimeclass / type=un subject=ID r rcorr;
run;

/*GEE (not valid)*/
proc genmod data=lda.acu2;
title ’data as is - GEE’;
class logtimeclass group id;
model severity=age chronicity group*logtime age*logtime chronicity*logtime /dist=normal;
repeated subject=id / withinsubject=logtimeclass type=un modelse;
run;
proc glimmix data=lda.acu2;
title ’data as is - GEE - linearized version’;
nloptions maxiter=50 technique=newrap;
class logtimeclass group id;
model severity=age chronicity group*logtime age*logtime chronicity*logtime /dist=normal;
random _residual_ / subject=id type=un;
run;
proc glimmix data=lda.acu2 empirical;
title ’data as is - GEE - linearized version - empirical’;
nloptions maxiter=50 technique=newrap;
class logtimeclass group id;
model severity=age chronicity group*logtime age*logtime chronicity*logtime /dist=normal;
random _residual_ / subject=id type=un;
run;


/*******************************************************/
/*WGEE (not valid with intermittent missing data)*/
/*******************************************************/

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


%dropout(data=lda.acu2,id=id,time=logtime,response=severity,out=lda.acu3);

proc genmod data=lda.acu3 descending;
class group prev logtimeclass;
model dropout = prev group age chronicity logtime / pred dist=b;
ods output obstats=pred;
ods listing exclude obstats;
run;

data pred;
set pred;
keep observation pred;
run;
data lda.acu3b;
merge pred lda.acu3;
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

%dropwgt(data=lda.acu3b,id=id,time=time,pred=pred,dropout=dropout,out=lda.acu4);


/* WGEE: models */

proc genmod data=lda.acu4;
title ’data as is - WGEE’;
scwgt wi;
class logtimeclass group id;
model severity=age chronicity group*logtime age*logtime chronicity*logtime /dist=normal;
repeated subject=id / withinsubject=logtimeclass type=un modelse;
run;
proc glimmix data=lda.acu4;
title ’data as is - WGEE - linearized version’;
nloptions maxiter=50 technique=newrap;
weight wi;
class logtimeclass group id;
model severity=age chronicity group*logtime age*logtime chronicity*logtime /dist=normal;
random _residual_ / subject=id type=un;
run;
proc glimmix data=lda.acu4 empirical;
title ’data as is - WGEE - linearized version - empirical’;
weight wi;
nloptions maxiter=50 technique=newrap;
class logtimeclass group id;
model severity=age chronicity group*logtime age*logtime chronicity*logtime /dist=normal;
random _residual_ / subject=id type=un;
run;


/********************************************/
/*Multiple imputation*/
/********************************************/

/*1. MI: imputation task (MCMC method)*/
proc mi data=lda.acu2b seed=486048 simple out=lda.acu5 nimpute=10 round=0.1;
var age chronicity severity0 severity3 severity12;
by group;
run;

/*wide to long*/
data lda.acu5b;
set lda.acu5;
array x (3) severity0 severity3 severity12;
do j=1 to 3;
severity=x(j);
time=j;
output;
end;
run;
data lda.acu5b;
set lda.acu5b;
drop severity0 severity3 severity12 j;
if time eq 1 then logtime=log(1+0);
if time eq 2 then logtime=log(1+3);
if time eq 3 then logtime=log(1+12);
drop time;
logtimeclass=logtime;
if group eq 0 then group0=1;
	else group0=0;
if group eq 1 then group1=1;
	else group1=0;
run;

proc sort data=lda.acu5b;
by _imputation_;
run;

/*2.1. MI: analysis task: GEE*/
proc genmod data=lda.acu5b;
title ’GEE after multiple imputation’;
class logtimeclass group id;
by _imputation_;
model severity=age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime /dist=normal covb;
repeated subject=id / withinsubject=logtimeclass type=un modelse;
ods output GEEEmpPEst=geeparms parminfo=geepinfo CovB=geecovb;
run;

/*3.1. MI: inference task: GEE*/
proc mianalyze parms=geeparms covb=geecovb parminfo=geepinfo wcov bcov tcov;
modeleffects intercept age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime;
run;


/*2.2. MI: analysis task: LMM*/
proc mixed data=lda.acu5b method=ml asycov covtest nobound;
class logtimeclass group id;
by _imputation_;
model severity= age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime/ s covb;
random intercept logtime  /type=un subject=ID g gcorr v vcorr solution;
repeated logtimeclass / type=un subject=ID r rcorr;
ods output solutionf = lmsolution covb = lmcovb covparms = lmcovparms asycov = lmasycov;
run;

/*??*/
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



/*3.2. MI: inference task: LMM*/
/* Combining 10 Separate Analyses (mean structure) */
proc mianalyze parms=lmsolution covb(effectvar=rowcol)=lmcovb;
title2 ’COMBINING 5 MIXED MODEL ANALYSES (MEAN STRUCTURE)’;
modeleffects intercept age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime;
run;
/* Combining 10 Separate Analyses (covariance structure) */
proc mianalyze parms=covparms0 covb(effectvar=rowcol)=asycov0;
title2 ’COMBINING 5 MIXED MODEL ANALYSES (COVARIANCE STRUCTURE)’;
modeleffects YEARUN11 YEARUN21 YEARUN22 YEARUN31 YEARUN32 YEARUN33 PARENTCO;
run;



/*****************************************/
/*MI-WGEE (for intermittent missingness)*/
/******************************************/

proc mi data=lda.acu2b seed=486048 simple out=lda.acu6 nimpute=10 round=0.1;
title "Monotone multiple imputation";
mcmc impute=monotone;
var age chronicity severity0 severity3 severity12;
by group;
run;

/*wide to long*/
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
if time eq 1 then logtime=log(1+0);
if time eq 2 then logtime=log(1+3);
if time eq 3 then logtime=log(1+12);
drop time;
logtimeclass=logtime;
if group eq 0 then group0=1;
	else group0=0;
if group eq 1 then group1=1;
	else group1=0;
run;

proc sort data=lda.acu6b;
by _imputation_;
run;


/*WGEE*/
%dropout(data=lda.acu6b,id=id,time=logtime,response=severity,out=lda.acu6c);

proc gee data=lda.acu6c;
by _imputation_;
class id group logtimeclass;
model severity= age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime/ dist=normal;
repeated subject=id /within=logtimeclass corr=un ecovb;
missmodel prev age chronicity group /type=obslevel;
ods output GEEEmpPEst=wgeeparms parminfo=wgeepinfo modelinfo=wgeemodelinfo GEERCov=wgeecovb;
run;

proc mianalyze parms=wgeeparms covb=wgeecovb parminfo=wgeepinfo wcov bcov tcov;
modeleffects intercept age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime;
run;




