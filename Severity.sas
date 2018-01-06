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

data lda.acu3;
set lda.acu2b;
diff3=severity3-severity0;
diff12=severity12-severity0;
bindif3=0; if diff3 <= 0 then bindif3=1;
bindif12=0;if diff12 <= 0 then bindif12=1;
if diff3=. then bindif3=.;
if diff12=. then bindif12=.;
subject=_n_;
run;
proc sort data=lda.acu3;
by group;
run;

/*Overview of missigness table*/
proc mi data=lda.acu3 seed=675938 simple nimpute=0;
title ’standard EM’;
em itprint outem=growthem1;
var diff3 diff12;
run;
/*Overview of missigness table per group*/
proc mi data=lda.acu3 seed=675938 simple nimpute=0;
title ’standard EM’;
em itprint outem=growthem1;
var diff3 diff12;
by group;
run;


/* Using data as it is.*/

/*LMM*/
proc mixed data=lda.acu4 method=ml;
class timeclass group id;
model severity= age chronicity time group*time age*time chronicity*time/ s ddfm=kr;
random intercept time  /type=un subject=ID g gcorr v vcorr solution;
repeated timeclass / type=un subject=ID r rcorr;
run;

/*GEE*/
proc genmod data=lda.acu4;
title ’data as is - GEE’;
class timeclass group id;
model severity=age chronicity time group*time age*time chronicity*time /dist=normal;
repeated subject=id / withinsubject=timeclass type=un modelse;
run;
proc glimmix data=lda.acu4;
title ’data as is - GEE - linearized version’;
nloptions maxiter=50 technique=newrap;
class timeclass group id;
model severity=age chronicity time group*time age*time chronicity*time /dist=normal;
random _residual_ / subject=id type=un;
run;
proc glimmix data=lda.acu4 empirical;
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






