/*Dropout analysis MCAR vs MAR psi1*/
/********************************************/

%dropout(data=lda3.acu2,id=id,time=logtime,response=severity,out=lda3.acupsi1);

proc genmod data=lda3.acupsi1 descending;
class group;
model dropout = prev group / pred dist=b;
run;

proc genmod data=lda3.acupsi1 descending;
class group;
model dropout = prev group age chronicity group*logtime age*logtime chronicity*logtime/ pred dist=b;
run;


/********************************************/
/*Sensitivity Analysis - General*/
/********************************************/

/*Shift*/
/********************************************/
proc mi data=lda3.acu2b seed=486048 simple out=lda3.acu2bshift nimpute=20 round=0.1;
title ’Shift multiple imputation’;
class group;
var age chronicity severity0 severity3 severity12;
fcs reg;
mnar adjust (severity3 / shift=10 adjustobs=(group=’1’));
mnar adjust (severity12 / shift=20 adjustobs=(group=’1’));
by group;
run;

/*wide to long shift*/
data lda3.acu2bshiftb;
set lda3.acu2bshift;
array x (3) severity0 severity3 severity12;
do j=1 to 3;
severity=x(j);
time=j;
output;
end;
run;
data lda3.acu2bshiftb;
set lda3.acu2bshiftb;
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

proc sort data=lda3.acu2bshiftb;
by _imputation_;
run;

/*Shift analysis task: GEE*/
proc genmod data=lda3.acu2bshiftb;
title ’GEE after multiple imputation’;
class logtimeclass group id;
by _imputation_;
model severity=age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime /dist=normal covb;
repeated subject=id / withinsubject=logtimeclass type=un modelse;
ods output GEEEmpPEst=geeparmsshift parminfo=geepinfoshift CovB=geecovbshift;
run;

/*Shift inference task: GEE*/
proc mianalyze parms=geeparmsshift covb=geecovbshift parminfo=geepinfoshift wcov bcov tcov;
modeleffects intercept age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime;
run;

/*Shift: analysis task: LMM*/
proc mixed data=lda3.acu2bshiftb method=reml asycov covtest nobound;
class logtimeclass group id;
by _imputation_;
model severity= age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime/ s covb;
random intercept logtime  /type=un subject=ID g gcorr v vcorr solution;
repeated logtimeclass / type=un subject=ID r rcorr;
ods output solutionf = lmsolutionshift covb = lmcovbshift covparms = lmcovparmsshift asycov = lmasycovshift;
run;

/*Shift: inference task: LMM*/
/* Combining 20 Separate Analyses (mean structure) */
proc mianalyze parms=lmsolutionshift covb(effectvar=rowcolshift)=lmcovbshift;
title2 ’COMBINING MIXED MODEL ANALYSES (MEAN STRUCTURE)’;
modeleffects intercept age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime;
run;

/*Subgroup*/
*proc mi data=lda3.acu2b seed=486048 simple out=lda3.acu2bsubgroup nimpute=20;
*title ’Model multiple imputation’;
*class group;
*var age chronicity severity0 severity3 severity12;
*fcs reg;
*mnar model (severity0 / modelobs= (group=’0’));
*mnar model (severity3 / modelobs= (group=’0’));
*mnar model (severity12 / modelobs= (group=’0’));
*run;


/*NCMV*/
/*MI for data provided in previous chapter*/
proc mi data=lda3.acu4 seed=486048 simple out=lda3.acu4NCMV nimpute=1;
title ’Model multiple imputation’;
var age chronicity severity0 severity3 severity12;
monotone reg;
mnar model (severity0 severity3 severity12 / modelobs=ncmv);
by group;
run;

/*wide to long NCMV*/
data lda3.acu4NCMVb;
set lda3.acu4NCMV;
array x (3) severity0 severity3 severity12;
do j=1 to 3;
severity=x(j);
time=j;
output;
end;
run;
data lda3.acu4NCMVb;
set lda3.acu4NCMVb;
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

proc sort data=lda3.acu4NCMVb;
by _imputation_;
run;

/*NCMV analysis task: GEE*/
proc genmod data=lda3.acu4NCMVb;
title ’GEE after multiple imputation’;
class logtimeclass group id;
by _imputation_;
model severity=age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime /dist=normal covb;
repeated subject=id / withinsubject=logtimeclass type=un modelse;
ods output GEEEmpPEst=geeparmsNCMV parminfo=geepinfoNCMV CovB=geecovbNCMV;
run;

/*NCMV inference task: GEE*/
proc mianalyze parms=geeparmsNCMV covb=geecovbNCMV parminfo=geepinfoNCMV wcov bcov tcov;
modeleffects intercept age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime;
run;

/*NCMV: analysis task: LMM*/
proc mixed data=lda3.acu4NCMVb method=reml asycov covtest nobound;
class logtimeclass group id;
by _imputation_;
model severity= age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime/ s covb;
random intercept logtime  /type=un subject=ID g gcorr v vcorr solution;
repeated logtimeclass / type=un subject=ID r rcorr;
ods output solutionf = lmsolutionNCMV covb = lmcovbNCMV covparms = lmcovparmsNCMV asycov = lmasycovNCMV;
run;

/*NCMV: inference task: LMM*/
/* Combining 20 Separate Analyses (mean structure) */
proc mianalyze parms=lmsolutionNCMV covb(effectvar=rowcolNCMV)=lmcovbNCMV;
title2 ’COMBINING MIXED MODEL ANALYSES (MEAN STRUCTURE)’;
modeleffects intercept age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime;
run;




/*CCMV*/
/*MI for data provided in previous chapter*/
proc mi data=lda3.acu4 seed=486048 simple out=lda3.acu4CCMV nimpute=1;
title ’Model multiple imputation’;
var age chronicity severity0 severity3 severity12;
monotone reg;
mnar model (severity0 severity3 severity12 / modelobs=ccmv);
by group;
run;

/*wide to long CCMV*/
data lda3.acu4CCMVb;
set lda3.acu4CCMV;
array x (3) severity0 severity3 severity12;
do j=1 to 3;
severity=x(j);
time=j;
output;
end;
run;
data lda3.acu4CCMVb;
set lda3.acu4CCMVb;
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

proc sort data=lda3.acu4CCMVb;
by _imputation_;
run;

/*CCMV analysis task: GEE*/
proc genmod data=lda3.acu4CCMVb;
title ’GEE after multiple imputation’;
class logtimeclass group id;
by _imputation_;
model severity=age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime /dist=normal covb;
repeated subject=id / withinsubject=logtimeclass type=un modelse;
ods output GEEEmpPEst=geeparmsCCMV parminfo=geepinfoCCMV CovB=geecovbCCMV;
run;

/*CCMV inference task: GEE*/
proc mianalyze parms=geeparmsCCMV covb=geecovbCCMV parminfo=geepinfoCCMV wcov bcov tcov;
modeleffects intercept age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime;
run;

/*CCMV: analysis task: LMM*/
proc mixed data=lda3.acu4CCMVb method=reml asycov covtest nobound;
class logtimeclass group id;
by _imputation_;
model severity= age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime/ s covb;
random intercept logtime  /type=un subject=ID g gcorr v vcorr solution;
repeated logtimeclass / type=un subject=ID r rcorr;
ods output solutionf = lmsolutionCCMV covb = lmcovbCCMV covparms = lmcovparmsCCMV asycov = lmasycovCCMV;
run;

/*CCMV: inference task: LMM*/
/* Combining 20 Separate Analyses (mean structure) */
proc mianalyze parms=lmsolutionCCMV covb(effectvar=rowcolCCMV)=lmcovbCCMV;
title2 ’COMBINING MIXED MODEL ANALYSES (MEAN STRUCTURE)’;
modeleffects intercept age chronicity group0*logtime group1*logtime age*logtime chronicity*logtime;
run;



/*___________*/





/*-------------------------*/



/*Local Influence*/

%macro LocInfl(data=, idno =, varesp=, fixedef=, randomef=, outdata=);
/* Initialization */
%if %bquote(&data)= %then %let data=&syslast;
/* To calculate the number of repeated measures per subject and the time points*/
proc freq data=&data noprint;
tables &idno /out=nfrec;
run;
proc sort data=&data;
by &idno;
run;
/* Fitting the model in Proc mixed */
proc mixed data=&data method=ml asycov;
class &idno;
model &varesp = &fixedef /noint covb s;
random &randomef /type=un sub=&idno g;
ods output solutionf= fixedsol;
ods output covparms= covpar;
ods output G= Dmatr;
ods output asycov = varcov;
ods output covb= covfixed;
run;
 
proc iml;
reset noprint;
/* Matrix with the solution of the fixed effects */
use fixedsol;
read all into fixedsol;
fixedsol= fixedsol[,1];
/* Matrix with the solution for the var-cov of the random. eff. parameters */
use covpar;
read all into covpar;
n_covp = nrow(covpar);
/* Var-Cov matrix for the fixed effects */
use covfixed;
read all into fixedvar;
n_fix = nrow(fixedvar);
fixedvar = fixedvar[,2:1+n_fix];
fixedvar= -inv(fixedvar);
/* Var-Cov matrix for the random components */
use varcov;
read all into varcov;
n_vc = nrow(varcov);
varcov = varcov[,1:n_vc];
varcov= -inv(varcov);
/* Matrix with the data to obtain the design matrices for the fixed and random effects */
use &data;
labelx = {&fixedef};
labelz = {&randomef};
labely = {&varesp};
read all var labelx into fixed;
read all var labelz into random;
read all var labely into resp;
/* Matrix with the frequencies for each of the subjects */
use nfrec;
read all into nfrec;
nfrec = nfrec[,2];
/* G=D matrix */
use Dmatr;
read all into D;
n_D = nrow(D);
D = D[,3:2+n_D];
/* Procedure to calculate the first derivative */
start Deriv(fixed, random, res, D, sigma2);
/*** With respect to alpha ***/
ni = nrow(fixed);
p = ncol(fixed);
q = ncol(random);
V = random*D*random` + sigma2*I(ni);
W = inv(V);
zero = {0};
help1=repeat(zero,p,1);
help2 =-2#fixed`*W*res;
first = help1 + help2;
/*** With respect to D ***/
j = 1;
do while (j <= q);
i = 1;
do while (i <= j);
if i=j then mul=1;
else mul=0;
hulp1 = (2-(mul))#random[,i]`*W*random[,j];
hulp2 = -(2-(mul))#random[,i]`*W*res*res`*W*random[,j]; /* Check this derivatives */
first = first//(hulp1 + hulp2);
i = i + 1;
end;
j = j + 1;
end;
/*** With respect to sigma2 ***/
hulp1 = trace(W);
hulp2 = -res`*W*W*res;
first3 = hulp1 + hulp2;
first = first//first3;
first = -first/2;
return(first);
finish;
/* Procedure to calculate the second derivative */
start SecDer(fixed, random, res, D, sigma2);
/* 2nd derivative alpha and D */
ni = nrow(fixed);
p = ncol(fixed);
q = ncol(random);
V = random*D*random` + sigma2*I(ni);
W = inv(V);
zero = {0};
j=1;
do while (j<=q);
i=1;
do while (i<=j);
help1 = 0;
help2 = 2#fixed`*W*random[,i]*random[,j]`*W*res;
if i=j then mul=1;
else mul=0;
help2 = help2 + 2#(1-mul)#fixed`*W*random[,j]*random[,i]`*W*res;
ret = help1` + help2`; /* ver como arreglar esto */
fila = fila//ret;
i = i+1;
end;
j= j+1;
end;
/* the 2nd derivative alpha and sigma2 */
help1 = 0;
help2 = 2#fixed`*W*W*res;
ret = help1` + help2`;
fila = fila//ret;
return(fila);
finish;
/* Initializations for the whole loop and calculations of the derivatives */
q = ncol(random);
p = ncol(fixed);
uno = {1};
begin = 1;
sum=repeat(0,(q*(q+1)/2)+1,p);
n=nrow(nfrec);
sigma2 = covpar[n_covp,];
s=1;
do s=1 to n ;
ni = nfrec[s,];
end = begin -1 + ni;
fixedi = fixed[begin:end,];
randomi = random[begin:end,];
respi = resp[begin:end];
res = respi-fixedi*fixedsol;
begin = end +1;
first = Deriv(fixedi, randomi, res, D, sigma2);
firstT = firstT||first;
second = SecDer(fixedi, randomi, res, D, sigma2);
sum = sum - (second/2);
end;
L = (fixedvar//sum)||(sum`//varcov);
vect = repeat(sqrt(2)-1,q,1);
help = diag(vect)+1;
help = symsqr(help);
help = repeat(uno,p,1)//help//1;
firstT = firstT#help;
L = L#(help*help`);
Linv = inv(-L);
b22 =repeat(0,p+(q*(q+1)/2)+1,p+(q*(q+1)/2)+1);
b11 = b22;
help = L[(p+1):(nrow(L)), (p+1):(ncol(L))];
help = -inv(-help);
b22[(p+1):(nrow(L)), (p+1):(ncol(L))] = help;
b22 = -inv(-L) - b22;
help = L[1:p,1:p];
help = -inv(-help);
b11[1:p, 1:p] = help;
b11 = -inv(-L) -b11;
/* Part to calculate the influence measures */
zero = {0};
uno = {1};
begin = 1;
n=nrow(nfrec);
do s=1 to n ;
ni = nfrec[s,];
end = begin -1 + ni;
fixedi = fixed[begin:end,];
randomi = random[begin:end,];
respi = resp[begin:end];
res = respi-fixedi*fixedsol;
begin = end +1;
V = randomi*D*randomi` + sigma2*I(ni);
W = inv(V);
Ci = 2#firstT[,s]`*Linv*firstT[,s];
Cai = -(2#firstT[,s]`*b22*firstT[,s]);
Cbi = -(2#firstT[,s]`*b11*firstT[,s]);
onei = fixedi`*W*res;
onei = onei`*onei;
twoi = randomi`*W*randomi - randomi`*W*res*res`*W*randomi;
twoi = trace(twoi*twoi`);
threei =trace(W) - res`*W*W*res;
threei = (threei##2);
call eigen(egval,egvect,W);
WW = egvect*diag(sqrt(egval))*egvect`;
rri = WW*res;
xxi = WW*fixedi;
zzi = WW*randomi;
oneai = xxi*xxi`;
oneai = sqrt(trace(oneai*oneai`));
onebi = rri*rri`;
onebi = sqrt(trace(onebi*onebi`));
twoai = zzi*zzi`;
twoai = trace(twoai*twoai`);
twobi = I(ni) - rri*rri`;
twobi = trace(twobi*twobi`);
threeai = trace(W*W`);
C = C//Ci;
Ca = Ca//Cai;
Cb = Cb//Cbi;
one = one//onei;
two = two//twoi;
three = three//threei;
onea = onea//oneai;
oneb = oneb//onebi;
twoa = twoa//twoai;
twob = twob//twobi;
threea = threea//threeai;
end;
/*Setting the output dataset with the diagnostic measures */
out=C||Ca||Cb||onea||oneb||twoa||twob||threea;
/* onea-> norm(xixi'2);
oneb-> norm(R2);
twoa-> norm(zizi'2);
twob-> norm(I-riri'2);
threea-> norm(Vi2) */
varnames = {'C_i' 'C_i_F' 'C_i_V' 'F1' 'F2' 'V1' 'V2' 'V3'};
create &outdata from out [colname= varnames];
append from out;
close covfixed;
close nfrec;
close &data;
close varcov;
close covpar;
close fixedsol;
quit;
/* End of the IML procedure */
data temp;
set nfrec;
keep &idno count;
run;
data &outdata;
merge temp &outdata;
run;
/* Printing the results */
proc print data=&outdata;
run;
/* Plotting the results */
goptions reset=goptions device=win targetdevice=winprtm rotate=landscape;
symbol v=dot i=join;
title f=swissb h=3 'Total local influence';
proc gplot data = &outdata;
plot C_i*&idno / haxis=axis1 vaxis=axis2;
axis1 label=(f=swissb h=2 'Identification number') w=5;
axis2 label=none w=5 ;
run;
goptions reset=goptions device=win targetdevice=winprtm rotate=landscape;
symbol v=dot i=join;
title f=swissb h=3 'Local influence for fixed effects';
proc gplot data = &outdata;
plot C_i_F*&idno / haxis=axis1 vaxis=axis2;
axis1 label=(f=swissb h=2 'Identification number') w=5;
axis2 label=none w=5 ;
run;
goptions reset=goptions device=win targetdevice=winprtm rotate=landscape;
symbol v=dot i=join;
title f=swissb h=3 'Local influence for variance components';
proc gplot data = &outdata;
plot C_i_V*&idno / haxis=axis1 vaxis=axis2;
axis1 label=(f=swissb h=2 'Identification number') w=5;
axis2 label=none w=5 ;
run;
goptions reset=goptions device=win targetdevice=winprtm rotate=landscape;
symbol v=dot i=join;
title f=swissb h=3 'Covariates in mean structure';
proc gplot data = &outdata;
plot F1*&idno / haxis=axis1 vaxis=axis2;
axis1 label=(f=swissb h=2 'Identification number') w=5;
axis2 label=none w=5 ;
run;
goptions reset=goptions device=win targetdevice=winprtm rotate=landscape;
symbol v=dot i=join;
title f=swissb h=3 'Residuals for mean structure';
proc gplot data = &outdata;
plot F2*&idno / haxis=axis1 vaxis=axis2;
axis1 label=(f=swissb h=2 'Identification number') w=5;
axis2 label=none w=5 ;
run;
goptions reset=goptions device=win targetdevice=winprtm rotate=landscape;
symbol v=dot i=join;
title f=swissb h=3 'Covariates in covariance structure';
proc gplot data = &outdata;
plot V1*&idno / haxis=axis1 vaxis=axis2;
axis1 label=(f=swissb h=2 'Identification number') w=5;
axis2 label=none w=5 ;
run;
goptions reset=goptions device=win targetdevice=winprtm rotate=landscape;
symbol v=dot i=join;
title f=swissb h=3 'Residuals for covariance structure';
proc gplot data = &outdata;
plot V2*&idno / haxis=axis1 vaxis=axis2;
axis1 label=(f=swissb h=2 'Identification number') w=5;
axis2 label=none w=5 ;
run;
goptions reset=goptions device=win targetdevice=winprtm rotate=landscape;
symbol v=dot i=join;
title f=swissb h=3 'Measure of variability';
proc gplot data = &outdata;
plot V3*&idno / haxis=axis1 vaxis=axis2;
axis1 label=(f=swissb h=2 'Identification number') w=5;
axis2 label=none w=5 ;
run;
goptions reset=goptions device=win targetdevice=winprtm rotate=landscape;
symbol v=dot i=join;
title f=swissb h=3 'Number of repeated measures per subject';
proc gplot data = &outdata;
plot count*&idno / haxis=axis1 vaxis=axis2;
axis1 label=(f=swissb h=2 'Identification number') w=5;
axis2 label=none w=5 ;
run;
title ' ';
%mend;

data lda3.acu2LI;
set lda3.acu2;
int=1;
acu = (group=1);
placebo = (group=0);
aclogtime = acu*logtime;
plalogtime = placebo*logtime;
agelogtime = age*logtime;
chronlogtime = chronicity*logtime;
run;

%LocInfl(data=lda3.acu2LI,
idno= id,
varesp=severity,
fixedef= age chronicity acu placebo aclogtime plalogtime agelogtime chronlogtime,
randomef= int logtime,
outdata = lda3.acu3out);


