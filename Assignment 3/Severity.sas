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
    id time;
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
/* Using data as it is: LMM*/
/********************************************/

/*LMM*/
/*method=ml*/
proc mixed data=lda.acu2 method=ml nobound;
class logtimeclass group id;
model severity= age chronicity group*logtime age*logtime chronicity*logtime/ solution;
random intercept logtime  /type=un subject=ID g gcorr v vcorr;
run;
/*method=reml*/
proc mixed data=lda.acu2 method=reml nobound;
class logtimeclass group id;
model severity= age chronicity group*logtime age*logtime chronicity*logtime/ solution;
random intercept logtime  /type=un subject=ID g gcorr v vcorr;
run;


/********************************************/
/*Multiple imputation-GEE*/
/********************************************/

/*1. MI: imputation task (MCMC method)*/
proc mi data=lda.acu2b seed=260811 simple out=lda.acu3 nimpute=30 round=0.1 minimum=0;
var age chronicity severity0 severity3 severity12;
by group;
run;

/*wide to long*/
data lda.acu3b;
set lda.acu3;
array x (3) severity0 severity3 severity12;
do j=1 to 3;
severity=x(j);
time=j;
output;
end;
run;
data lda.acu3b;
set lda.acu3b;
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

proc sort data=lda.acu3b;
by _imputation_;
run;

/*2.1. MI: analysis task: GEE*/
proc genmod data=lda.acu3b;
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


/********************************************/
/*Multiple imputation-LMM*/
/********************************************/

/*2.2. MI: analysis task: LMM*/
proc mixed data=lda.acu3b method=reml asycov covtest nobound;
class logtimeclass group id;
by _imputation_;
model severity= age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime/ s covb;
random intercept logtime  /type=un subject=ID g gcorr v vcorr solution;
ods output solutionf = lmsolution covb = lmcovb covparms = lmcovparms asycov = lmasycov;
run;


/*3.2. MI: inference task: LMM*/
/* Combining 30 Separate Analyses (mean structure) */
proc mianalyze parms=lmsolution covb(effectvar=rowcol)=lmcovb;
title2 ’COMBINING 5 MIXED MODEL ANALYSES (MEAN STRUCTURE)’;
modeleffects intercept age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime;
run;



/*****************************************/
/*MI-WGEE (for intermittent missingness)*/
/******************************************/

proc mi data=lda.acu2b seed=486048 simple out=lda.acu4 nimpute=10 round=0.1 minimum=0;
title "Monotone multiple imputation";
mcmc impute=monotone;
var age chronicity severity0 severity3 severity12;
by group;
run;

/*wide to long*/
data lda.acu4b;
set lda.acu4;
array x (3) severity0 severity3 severity12;
do j=1 to 3;
severity=x(j);
time=j;
output;
end;
run;
data lda.acu4b;
set lda.acu4b;
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

proc sort data=lda.acu4b;
by _imputation_;
run;


/*WGEE*/

/*macro for creating variables "dropout" and "prev" */
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

%dropout(data=lda.acu4b,id=id,time=logtime,response=severity,out=lda.acu4c);

proc gee data=lda.acu4c;
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


/*************************************/
/*MODELING USING LMM*/
/*************************************/

/*full model*/
proc mixed data=lda.acu2 method=ml nobound;
class logtimeclass group id;
model severity= age chronicity group*logtime age*logtime chronicity*logtime/ solution;
random intercept logtime  /type=un subject=ID g gcorr v vcorr;
run;
/*8090.5, df=11*/

/*model random effects covariance: only intercepts*/
proc mixed data=lda.acu2 method=ml nobound;
class logtimeclass group id;
model severity= age chronicity group*logtime age*logtime chronicity*logtime/ solution;
random intercept  /type=un subject=ID g gcorr v vcorr;
run;
/*8090.9,df=9*/
/*1-pchisq(8090.9-8090.5, 11-9)=0.8187308 ok*/


/*Variance plots*/
proc mixed data=lda.acu2 method=reml;
class id group logtimeclass;
model severity= age chronicity group*logtime age*logtime chronicity*logtime/ solution;
repeated logtimeclass / subject=id type=un r rcorr;
run;

/* Observed variance components of the response assuming full mean structure*/
data test;
input GROUP x y;
cards;
1 0 260.73 
1 1.3862943611 281.16 
1 2.5649493575 247.21
run;

/* creating the plot*/
goptions reset=all;
proc gplot data=test;
plot y*x / haxis=axis1 vaxis=axis2;
	symbol1 c=red value=dot i=join w=3 l=1 mode=include;
	symbol2 c=black i=join w=10 l=1 mode=include;
	axis1 label=(h=2 'log(TIME+1) ')order=(0 to 2.6 by 0.5) offset=(2)cm 
		minor=none w=5;
	axis2 label=(h=2 angle=90 'Var(Severity)') value=(h=2) w=5 order=(60 to 280);
run;


/*reduce mean structure: delete chronicity*logtime*/
proc mixed data=lda.acu2 method=ml nobound;
class logtimeclass group id;
model severity= age chronicity group*logtime age*logtime/ solution;
random intercept  /type=un subject=ID g gcorr v vcorr;
run;
/*8090.9,df=8*/
/*1-pchisq(8090.9-8090.9, 9-8)=0.8187308 ok*/

/*delete chronicity*/
proc mixed data=lda.acu2 method=reml nobound;
class logtimeclass group id;
model severity= age group*logtime age*logtime/ solution;
random intercept  /type=un subject=ID g gcorr v vcorr;
contrast 'Equal slopes' group*logtime 1 -1;
run;
/*8090.9,df=8*/
/*1-pchisq(8093.8-8090.9, 8-7)=0.08857955 ok*/
/*1-pchisq(8093.8-8090.5, 11-7)=0.5089 ok */
