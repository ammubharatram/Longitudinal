
/********************************************/
/*Sensitivity Analysis - General*/
/********************************************/


/*NCMV*/
/*MI for data provided in previous chapter*/
/* 1st step: Convert into monotoned data */
proc mi data=lda.acu2c seed=486048 simple out=lda.acu2c_m nimpute=10 minimum=0;
title ’Monotone imputation’;
var age chronicity severity0 severity3 severity12 logfreq0 logfreq3 logfreq12;
mcmc impute=monotone;
by group;
run;


/* 2nd step*/
proc mi data=lda.acu2c_m seed=486048 simple out=lda.acu4NCMV nimpute=1 minimum=0;
title ’Model multiple imputation’;
var age chronicity severity0 severity3 severity12 logfreq0 logfreq3 logfreq12;
monotone reg;
mnar model (severity0 severity3 severity12 logfreq0 logfreq3 logfreq12 / modelobs=ncmv);
by group;
run;

/*wide to long NCMV*/
data lda.acu4NCMVb;
set lda.acu4NCMV;
array f (3) logfreq0 logfreq3 logfreq12;
array s (3) severity0 severity3 severity12;
do j=1 to 3;
logfreq=f(j);
severity=s(j);
time=j+20;
output;
end;
run;
data lda.acu4NCMVb;
set lda.acu4NCMVb;
drop logfreq0 logfreq3 logfreq12 severity0 severity3 severity12 j;
if time eq 21 then time=0;
if time eq 22 then time=3;
if time eq 23 then time=12;
timeclass=time;
frequency=exp(logfreq);
run;

proc sort data=lda.acu4NCMVb;
by _imputation_;
run;


/*ncmv: analysis task: GEE*/

proc genmod data=lda.acu4ncmvb descending;
class id timeclass;
model frequency = severity age chronicity time severity*time age*time chronicity*time group*time / covb dist=poisson ;
repeated subject=id / withinsubject=timeclass type=ar  covb modelse;
by _imputation_;
ods output GEEEmpPEst=geeparms parminfo=geepinfo CovB=geecovb;
run;

proc mianalyze parms=geeparms covb=geecovb  wcov bcov tcov;
modeleffects intercept severity age chronicity time severity*time age*time chronicity*time group*time;
run;



/*ccmv*/
/*MI for data provided in previous chapter*/
/* 1st step: Convert into monotoned data */
proc mi data=lda.acu2c seed=486048 simple out=lda.acu2c_m nimpute=10 minimum=0;
title ’Monotone imputation’;
var age chronicity severity0 severity3 severity12 logfreq0 logfreq3 logfreq12;
mcmc impute=monotone;
by group;
run;


/* 2nd step*/
proc mi data=lda.acu2c_m seed=486048 simple out=lda.acu4ccmv nimpute=1 minimum=0;
title ’Model multiple imputation’;
var age chronicity severity0 severity3 severity12 logfreq0 logfreq3 logfreq12;
monotone reg;
mnar model (severity0 severity3 severity12 logfreq0 logfreq3 logfreq12 / modelobs=ccmv);
by group;
run;

/*wide to long ccmv*/
data lda.acu4ccmvb;
set lda.acu4ccmv;
array f (3) logfreq0 logfreq3 logfreq12;
array s (3) severity0 severity3 severity12;
do j=1 to 3;
logfreq=f(j);
severity=s(j);
time=j+20;
output;
end;
run;
data lda.acu4ccmvb;
set lda.acu4ccmvb;
drop logfreq0 logfreq3 logfreq12 severity0 severity3 severity12 j;
if time eq 21 then time=0;
if time eq 22 then time=3;
if time eq 23 then time=12;
timeclass=time;
frequency=exp(logfreq);
run;

proc sort data=lda.acu4ccmvb;
by _imputation_;
run;


/*ccmv: analysis task: GEE*/

proc genmod data=lda.acu4ccmvb descending;
class id timeclass group;
model frequency = severity age chronicity time severity*time age*time chronicity*time group*time / covb dist=poisson ;
repeated subject=id / withinsubject=timeclass type=ar  covb modelse;
by _imputation_;
ods output GEEEmpPEst=geeparms parminfo=geepinfo CovB=geecovb;
run;

proc mianalyze parms=geeparms covb=geecovb parminfo=geepinfo wcov bcov tcov;
modeleffects intercept severity age chronicity time severity*time age*time chronicity*time group*time;
run;
















/*Shift*/
/********************************************/
proc mi data=lda.acu2c seed=486048 simple out=lda.acu2cshift nimpute=20  minimum=0;
title ’Shift multiple imputation’;
class group;
var age chronicity severity0 severity3 severity12 logfreq0 logfreq3 logfreq12;
fcs reg;
mnar adjust (severity3 / shift=0 adjustobs=(group=´1´));
mnar adjust (severity12 / shift=1 adjustobs=(group=´1´));
mnar adjust (logfreq3 / shift=0 adjustobs=(group=´1´));
mnar adjust (logfreq12 / shift=1 adjustobs=(group=´1´));
by group;
run;

/*wide to long shift*/
data lda.acu2cshiftb;
set lda.acu2cshift;
array f (3) logfreq0 logfreq3 logfreq12;
array s (3) severity0 severity3 severity12;
do j=1 to 3;
logfreq=f(j);
severity=s(j);
time=j;
output;
end;
run;
data lda.acu2cshiftb;
set lda.acu2cshiftb;
drop logfreq0 logfreq3 logfreq12 severity0 severity3 severity12 j;
if group eq 0 then group0=1;
	else group0=0;
if group eq 1 then group1=1;
	else group1=0;
run;

proc sort data=lda.acu2cshiftb;
by _imputation_;
run;


/*Shift: analysis task: LMM*/
proc mixed data=lda.acu2cshiftb method=reml asycov covtest nobound;
class logtimeclass group id;
by _imputation_;
model logfreq= age severity chronicity group0*logtime group1*logtime age*logtime chronicity*logtime severity*logtime/ s covb;
random intercept logtime  /type=un subject=ID g gcorr v vcorr solution;
repeated logtimeclass / type=un subject=ID r rcorr;
ods output solutionf = lmsolutionshift covb = lmcovbshift covparms = lmcovparmsshift asycov = lmasycovshift;
run;

/*Shift: inference task: LMM*/
/* Combining 20 Separate Analyses (mean structure) */
proc mianalyze parms=lmsolutionshift covb(effectvar=rowcolshift)=lmcovbshift;
title2 ’COMBINING MIXED MODEL ANALYSES (MEAN STRUCTURE)’;
modeleffects intercept age severity chronicity group0*logtime group1*logtime age*logtime chronicity*logtime severity*logtime;
run;

/*Subgroup*/
proc mi data=lda.acu2c seed=486048 simple out=lda.acu2csubgroup nimpute=20;
title ’Model multiple imputation’;
class group;
var age chronicity logfreq0 logfreq3 logfreq12;
fcs reg;
mnar model (logfreq0 / modelobs= (group=’0’));
mnar model (logfreq3 / modelobs= (group=’0’));
mnar model (logfreq12 / modelobs= (group=’0’));
run;
