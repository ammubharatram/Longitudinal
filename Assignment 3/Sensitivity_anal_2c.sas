
/********************************************/
/*Sensitivity Analysis - General*/
/********************************************/

/*Shift*/
/********************************************/
proc mi data=lda3.acu2c_mi seed=486048 simple out=lda3.acu2cshift nimpute=20 round=0.1;
title ’Shift multiple imputation’;
class group;
var age chronicity severity0 severity3 severity12 frequency0 frequency3 frequency12;
fcs reg;
mnar adjust (severity3 / shift=10 adjustobs=(group=’1’));
mnar adjust (severity12 / shift=20 adjustobs=(group=’1’));
mnar adjust (frequency3 / shift=10 adjustobs=(group=’1’));
mnar adjust (frequency12 / shift=20 adjustobs=(group=’1’));
by group;
run;

/*wide to long shift*/
data lda3.acu2cshiftb;
set lda3.acu2cshift;
array f (3) frequency0 frequency3 frequency12;
array s (3) severity0 severity3 severity12;
do j=1 to 3;
frequency=f(j);
severity=s(j);
time=j;
output;
end;
run;
data lda3.acu2cshiftb;
set lda3.acu2cshiftb;
drop frequency0 frequency3 frequency12 severity0 severity3 severity12 j;
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

proc sort data=lda3.acu2cshiftb;
by _imputation_;
run;


/*Shift: analysis task: LMM*/
proc mixed data=lda3.acu2cshiftb method=reml asycov covtest nobound;
class logtimeclass group id;
by _imputation_;
model frequency= age severity chronicity group0*logtime group1*logtime age*logtime chronicity*logtime severity*logtime/ s covb;
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
*proc mi data=lda3.acu2c seed=486048 simple out=lda3.acu2csubgroup nimpute=20;
*title ’Model multiple imputation’;
*class group;
*var age chronicity frequency0 frequency3 frequency12;
*fcs reg;
*mnar model (frequency0 / modelobs= (group=’0’));
*mnar model (frequency3 / modelobs= (group=’0’));
*mnar model (frequency12 / modelobs= (group=’0’));
*run;


/*NCMV*/
/*MI for data provided in previous chapter*/
/* 1st step: Convert into monotoned data */
proc mi data=lda3.acu2c seed=486048 simple out=lda3.acu2c_m nimpute=1;
title ’Monotone imputation’;
var age chronicity severity0 severity3 severity12 frequency0 frequency3 frequency12;
mcmc impute=monotone;
by group;
run;


/* 2nd step*/
proc mi data=lda3.acu2c_m seed=486048 simple out=lda3.acu4NCMV nimpute=1;
title ’Model multiple imputation’;
var age chronicity severity0 severity3 severity12 frequency0 frequency3 frequency12;
monotone reg;
mnar model (severity0 severity3 severity12 frequency0 frequency3 frequency12 / modelobs=ncmv);
by group;
run;

/*wide to long NCMV*/
data lda3.acu4NCMVb;
set lda3.acu4NCMV;
array f (3) frequency0 frequency3 frequency12;
array s (3) severity0 severity3 severity12;
do j=1 to 3;
frequency=f(j);
severity=s(j);
time=j;
output;
end;
run;
data lda3.acu4NCMVb;
set lda3.acu4NCMVb;
drop frequency0 frequency3 frequency12 severity0 severity3 severity12 j;
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

proc sort data=lda3.acu4NCMVb;
by _imputation_;
run;


/*NCMV: analysis task: LMM*/
proc mixed data=lda3.acu4NCMVb method=reml asycov covtest nobound;
class logtimeclass group id;
by _imputation_;
model frequency= age severity chronicity group0*logtime group1*logtime age*logtime chronicity*logtime severity*logtime/ s covb;
random intercept logtime  /type=un subject=ID g gcorr v vcorr solution;
repeated logtimeclass / type=un subject=ID r rcorr;
ods output solutionf = lmsolutionNCMV covb = lmcovbNCMV covparms = lmcovparmsNCMV asycov = lmasycovNCMV;
run;

/*NCMV: inference task: LMM*/
/* Combining 20 Separate Analyses (mean structure) */
proc mianalyze parms=lmsolutionNCMV covb(effectvar=rowcolNCMV)=lmcovbNCMV;
title2 ’COMBINING MIXED MODEL ANALYSES (MEAN STRUCTURE)’;
modeleffects intercept age severity chronicity group0*logtime group1*logtime age*logtime chronicity*logtime severity*logtime;
run;




/*CCMV*/
/*MI for data provided in previous chapter*/
proc mi data=lda3.acu2c seed=486048 simple out=lda3.acu4CCMV nimpute=1;
title ’Model multiple imputation’;
var age chronicity severity0 severity3 severity12 frequency0 frequency3 frequency12;
monotone reg;
mnar model (severity0 severity3 severity12 frequency0 frequency3 frequency12 / modelobs=ccmv);
by group;
run;

/*wide to long CCMV*/
data lda3.acu4CCMVb;
set lda3.acu4CCMV;
array f (3) frequency0 frequency3 frequency12;
array s (3) severity0 severity3 severity12;
do j=1 to 3;
frequency=f(j);
severity=s(j);
time=j;
output;
end;
run;
data lda3.acu4CCMVb;
set lda3.acu4CCMVb;
drop frequency0 frequency3 frequency12 severity0 severity3 severity12 j;
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

proc sort data=lda3.acu4CCMVb;
by _imputation_;
run;


/*CCMV: analysis task: LMM*/
proc mixed data=lda3.acu4CCMVb method=reml asycov covtest nobound;
class logtimeclass group id;
by _imputation_;
model frequency= age severity chronicity group0*logtime group1*logtime age*logtime chronicity*logtime severity*logtime/ s covb;
random intercept logtime  /type=un subject=ID g gcorr v vcorr solution;
repeated logtimeclass / type=un subject=ID r rcorr;
ods output solutionf = lmsolutionCCMV covb = lmcovbCCMV covparms = lmcovparmsCCMV asycov = lmasycovCCMV;
run;

/*CCMV: inference task: LMM*/
/* Combining 20 Separate Analyses (mean structure) */
proc mianalyze parms=lmsolutionCCMV covb(effectvar=rowcolCCMV)=lmcovbCCMV;
title2 ’COMBINING MIXED MODEL ANALYSES (MEAN STRUCTURE)’;
modeleffects intercept age severity chronicity group0*logtime group1*logtime age*logtime chronicity*logtime severity*logtime;
run;