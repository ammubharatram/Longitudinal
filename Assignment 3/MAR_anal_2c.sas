/*1. MI: imputation task (FCS method)*/
proc mi data=lda3.acu2c seed=486048 simple out=lda3.acu3 nimpute=30 round=0.01;
fcs reg(severity0=age chronicity); 
fcs reg(severity3=age chronicity severity0); 
fcs reg(severity12=age chronicity severity0 severity3); 
fcs reg(frequency0=age chronicity severity0 severity3 severity12); 
fcs reg(frequency3=age chronicity severity0 severity3 severity12 frequency0); 
fcs reg(frequency12= age chronicity severity0 severity3 severity12 frequency0 frequency3); 
var age chronicity severity0 severity3 severity12 frequency0 frequency3 frequency12;
by group;
run;

/*wide to long*/
data lda3.acu2c_mi;
set lda3.acu3;
array f (3) frequency0 frequency3 frequency12;
array s (3) severity0 severity3 severity12;
do j=1 to 3;
frequency=f(j);
severity=s(j);
time=j;
output;
end;
run;
data lda3.acu2c_mi;
set lda3.acu2c_mi;
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

proc sort data=lda3.acu3b;
by _imputation_;
run;


/* use reduced mean structure below ******/
/********************************************/
/*Multiple imputation-LMM*/
/********************************************/

/*2. MI: analysis task: LMM*/
proc mixed data=lda3.acu2c_mi method=reml asycov covtest nobound;
class logtimeclass group id;
by _imputation_;
model frequency= age severity chronicity group0*logtime group1*logtime age*logtime chronicity*logtime severity*logtime/ s covb;
random intercept logtime  /type=un subject=ID g gcorr v vcorr solution;
ods output solutionf = lmsolution covb = lmcovb covparms = lmcovparms asycov = lmasycov;
run;


/*3. MI: inference task: LMM*/
/* Combining 20 Separate Analyses (mean structure) */
proc mianalyze parms=lmsolution covb(effectvar=rowcol)=lmcovb;
title2 ’COMBINING 5 MIXED MODEL ANALYSES (MEAN STRUCTURE)’;
modeleffects intercept age severity chronicity group0*logtime group1*logtime age*logtime chronicity*logtime severity*logtime;
run;

/*Continue with sensitivity analysis 2c to check if MAR is not satisfactory */
