/*IMPORT DATA FROM FILE*/
/*TAKE THE LOG OF TIME*/


data LDA.ACU; 
set LDA.aculda; 
log_time=log(1+time); 
run; 


/* 1. EXPLORATORY ANALYSIS */

/* 1.0. Descriptive Statistics */
proc freq data = LDA.ACU;
tables time*group/crosslist;
run;

proc means data = LDA.ACU maxdec = 3;
var severity age chronicity frequency;
class group;
where time=0;
run;

/*drop off rate*/
data LDA.ACU;
set LDA.aculda;
log_time=log(time);
log_sev=log(severity)
censored=0;
if severity ne "." then censored=1;
run;



/*covariate scatterplot at baseline*/
ods graphics on;
title 'Correlation of covariates';
proc corr data=LDA.ACU plots (MAXPOINTS=NONE)=matrix (histogram);
var severity age chronicity frequency;
where time=0;
run;
ods graphics off;



/*explore the mean structure */
/* 1.1.a Individual Profiles Combined*/
PROC SGPLOT data=LDA.ACU noautolegend;
    title "Individual profiles";
    series x=time y=severity / group=id lineattrs=(pattern=solid thickness=2) 
        transparency=0.00 name='Series';
    xaxis grid label='TIME (months)' values=(0 3 12);
    yaxis grid label='Severity'; 
    keylegend 'Series' / title='id:' location=Outside;
RUN;

/*subset of 30 individuals*/
proc surveyselect data=LDA.ACU method=srs n=30
                  out=LDA.SampleSRS seed=1234;
samplingunit id;
run;

PROC SGPLOT data=LDA.SampleSRS noautolegend;
    title "Individual profiles";
    series x=time y=severity / group=id lineattrs=(pattern=solid thickness=2) 
        transparency=0.00 name='Series';
    xaxis grid label='TIME (months)' values=(0 3 12);
    yaxis grid label='Severity'; 
    RUN;


/* 1.1.b Individual Profiles per Treatment Group*/
PROC SORT data=LDA.ACU; by group time; 
RUN;

PROC SGPLOT data=LDA.ACU noautolegend;
    title "Individual profiles";
    by group;
    series x=time y=severity / group=id lineattrs=(pattern=solid thickness=2) 
        transparency=0.00 name='Series';
    xaxis grid label='TIME (months)' values=(0 3 12);
    yaxis grid label='Headache Severity';
    keylegend 'Series' / title='id:' location=Outside;
RUN;

/*subset of 30 individuals*/
PROC SORT data=LDA.SampleSRS; by group time; 
RUN;

PROC SGPLOT data=LDA.SampleSRS noautolegend;
    title "Individual profiles";
	by group;
    series x=time y=severity / group=id lineattrs=(pattern=solid thickness=2) 
        transparency=0.00 name='Series';
    xaxis grid label='TIME (months)' values=(0 3 12);
    yaxis grid label='Severity'; 
RUN;



/* 1.2.a. Graphic MEANS and std dev */
PROC SGPLOT data=LDA.ACU;
title 'Mean of individual profiles';
    vline time / response=severity lineattrs=(thickness=2) transparency=0.00 
        stat=mean name='Line' limitstat=stderr;
    xaxis label='TIME (months)' values=(0 3 12);
    yaxis label="Headache Severity" min=10 max=35;
RUN;

/*without modifying x axis*/
proc gplot data=LDA.ACU;
plot Severity*time / haxis=axis1 vaxis=axis2;
symbol c=red i=std1mjt w=1 mode=include;
axis1 label=(h=2 "Time (months)") value=(h=1.5) order=(0 to 12 by 1) minor=none;
axis2 label=(h=2 A=90 "Severity") value=(h=1.5) order=(10 to 35)
minor=none;
title; 
run;quit;


/*Transformation: log(time+1)*/
proc gplot data=LDA.ACU;
plot Severity*log_time / haxis=axis1 vaxis=axis2;
symbol c=red i=std1mjt w=1 mode=include;
axis1 label=(h=2 "log(1+Time)") value=(h=1.5) order=(0 to 3) minor=none;
axis2 label=(h=2 A=90 "Severity") value=(h=1.5) order=(10 to 35)
minor=none;
title; 
run;quit;


/*1.2.b. Graphic MEANS and std dev for each neuro */
proc sort data=LDA.ACU;
by group time;run;
PROC SGPLOT data=LDA.ACU;
title 'Mean of individual profiles';
by group;
    vline time / response=severity lineattrs=(thickness=2) transparency=0.00 
        stat=mean name='Line' limitstat=stderr;
    xaxis label='TIME (months)';
    yaxis label="Heqdqche Severity" min=0 max=100;
RUN;



/*1.2.c. MEANS for each neuro Combined */
                                                                                                                                        
/* Calculate the mean and standard error for each X */                                                                                  
proc means data=LDA.ACU noprint;                                                                                                           
   by group time;                                                                                                                    
   var severity;                                                                                                                            
   output out=meansout(drop=_type_ _freq_) mean=mean stderr=stderr;                                                                     
run;                                                                                                                                    
                                                                                                                                        
/* Reshape the data to contain three Y values for */                                                                                    
/* each X for use with the HILOC interpolation.   */                                                                                    
data reshape(keep=group time severity mean);                                                                                             
   set meansout;                                                                                                                        
   by group time;                                                                                                                    
                                                                                                                                        
/* Offset the X values to display two groups */                                                                                         
   if group=0 then time=time - 0.08;                                                                                               
   else if group=1 then time=time + 0.08;                                                                                          
                                                                                                                                        
   severity=mean;                                                                                                                           
   output;                                                                                                                              
                                                                                                                                        
   severity=mean - stderr;                                                                                                                  
   output;                                                                                                                              
                                                                                                                                        
   severity=mean + stderr;                                                                                                                  
   output;                                                                                                                              
run;                                                                                                                                    
                                                                                                                                        
/* Define the title */                                                                                                                  
   title1 'Means by Group';                                                        
                                                                                                                                        
/* Define the axis characteristics */                                                                                                   
   axis1 offset=(0,0) minor=none value=(t=1 ' ' t=7 ' ');                                                                                                       
   axis2 label=(angle=90) order=(10 to 35 by 5) minor=(n=1);                                                                                                                   
                                                                                                                                        
/* Define the symbol characteristics */                                                                                                 
   symbol1 interpol=hiloctj color=vibg line=1;                                                                                          
   symbol2 interpol=hiloctj color=depk line=2;                                                                                          
                                                                                                                                        
   symbol3 interpol=none color=vibg value=dot height=1.5;                                                                               
   symbol4 interpol=none color=depk value=dot height=1.5;                                                                               
                                                                                                                                        
/* Define the legend characteristics */                                                                                                 
   legend1 label=('Group:') frame;                                                                                                      
                                                                                                                                        
/* Plot the error bars using the HILOCTJ interpolation */                                                                               
/* and overlay symbols at the means. */                                                                                                 
proc gplot data=reshape;                                                                                                                
   plot severity*time=group / haxis=axis1 vaxis=axis2 legend=legend1;                                                                    
   plot2 mean*time=group / vaxis=axis2 noaxis nolegend;
   
run;          




/*the variance structure and the correlation structure. */ 

/* 1.3 icc*/
/* icc for identifying variablility within patients*/
proc mixed data = LDA.ACU;
class id group;
model severity = group time time*group/ solution;
random id;
run;


/* 1.4 Correlation and Variance Structure */
data LDA.ACU; 
	set LDA.ACU;
timeclass=time; 
run;
 
/* 1.4.a Without Covariates */ 
proc mixed data=LDA.ACU ; 
class id group timeclass; 
model severity=group time group*time/ solution; 
repeated timeclass / subject=id type=un r rcorr; 
run;

data LDA.test1;
input group x1 y1;
cards;
1 1	 261.93	
1 3	 286.55	
1 12 250.06
run;

/* 1.4.b With Covariates */ 
proc mixed data=LDA.ACU; 
class id group timeclass; 
model severity=group age chronicity frequency time group*time/ solution; 
repeated timeclass / subject=id type=un r rcorr; 
run;

data LDA.test2;
input group x1 y2;
cards;
1 1	 94.4596	
1 3	 177.74		
1 12 171.44
run;


/* 1.4.c Merging for plotting */ 
data LDA.test (keep= group x1 y1 y2);
merge LDA.test1 LDA.test2;
by x1;
run;

/*--Set output size--*/
ods graphics / reset imagemap;

/*--SGPLOT proc statement--*/
proc sgplot data=LDA.TEST noautolegend;
    /*--TITLE and FOOTNOTE--*/
    title "Observed variance with/without covariates";
	footnote j=l "1= without covariates | 2=with covariates";
    /*--Scatter plot settings--*/
    series x=x1 y=y1 / curvelabel='1' curvelabelpos=max 
        curvelabelattrs=(size=7) transparency=0.0 name='Series1';
    series x=x1 y=y2 / curvelabel='2' curvelabelpos=max 
        curvelabelattrs=(size=7) transparency=0.0 name='Series2';
    xaxis label="TIME (Months)" values=(1 5 12) grid;
    yaxis grid label="VARIANCE Severity" grid;
run;


/* 1.4.1 Correlation and Variance Structure log transformed time*/
data LDA.ACU; 
	set LDA.ACU;
timeclass=time; 
run;
 
/* 1.4.1.a Without Covariates */ 

data LDA.test4;
input group x1 y1;
cards;
1 0	 261.93	
1 1.1 286.55	
1 2.5 250.06
run;

/* 1.4.1.b With Covariates */ 
data LDA.test5;
input group x1 y2;
cards;
1 0	 94.4596	
1 1.1 177.74		
1 2.5 171.44
run;


/* 1.4.1.c Merging for plotting */ 
data LDA.test3 (keep= group x1 y1 y2);
merge LDA.test4 LDA.test5;
by x1;
run;

/*--Set output size--*/
ods graphics / reset imagemap;

/*--SGPLOT proc statement--*/
proc sgplot data=LDA.TEST3 noautolegend;
    /*--TITLE and FOOTNOTE--*/
    title "Observed variance with/without covariates";
	footnote j=l "1= without covariates | 2=with covariates";
    /*--Scatter plot settings--*/
    series x=x1 y=y1 / curvelabel='1' curvelabelpos=max 
        curvelabelattrs=(size=7) transparency=0.0 name='Series1';
    series x=x1 y=y2 / curvelabel='2' curvelabelpos=max 
        curvelabelattrs=(size=7) transparency=0.0 name='Series2';
    xaxis label="TIME (Log(time+1))" values=(0 1.1 2.5) grid;
    yaxis grid label="VARIANCE Severity" grid;
run;




/*2. Variance structure*/

/*check the mean severity at every time point*/
proc means data = LDA.ACU maxdec = 3;
var severity;
class time;
run;

/*compute the squared residuals*/
data LDA.Residuals;
set LDA.ACU ;
if time eq 0 then res=(severity-26.510)**2;
if time eq 3 then res=(severity-21.598)**2;
if time eq 12 then res=(severity-19.082)**2;
run;

/*plot squared residuals vs time*/
goptions reset=ALL;
proc gplot data=LDA.Residuals;
plot res*time / haxis=axis1 vaxis=axis2;
by group;
symbol c=red i=std1mjt w=1 mode=include;
axis1 label=(h=2 "Time (months)") value=(h=1.5) order=(0 to 12 by 1) minor=none;
axis2 label=(h=2 A=90 "Squared residuals") value=(h=1.5) order=(150 to 300)
minor=none;
title; 
run;quit;

/*Residuals after removing for systematic trends*/
proc glm data=LDA.ACU;
model severity=time group age chronicity frequency;
output out=LDA.out r=residual;
run;

data Out;
set Out;
residual2=residual**2;
run;

proc gplot data=Out;
plot residual2*time / haxis=axis1 vaxis=axis2;
symbol c=red i=std1mjt w=1 mode=include;
axis1 label=(h=2 "Time (months)") value=(h=1.5) order=(0 to 12 by 1) minor=none;
axis2 label=(h=2 A=90 "Squared residuals") value=(h=1.5) order=(50 to 200)
minor=none;
title; 
run;quit;


/*3. Correlation structure*/

proc sort data=LDA.ACU;
by id time;run;
proc transpose data=LDA.ACU out=LDA.wide prefix=time;
    by id ;
    id time;
    var severity;
run;

ods graphics on;
title 'Correlation structure';
proc corr data=LDA.wide plots (MAXPOINTS=NONE)=matrix (histogram);
var time0 time3 time12;
run;
ods graphics off;





data lda.acu2;
set lda.acu;
log_time = log(time+1);
run;


ods graphics on; 
title ’data:Correlation structure’;
proc corr data= lda.aculda2 plots(MAXPOINTS=NONE)=matrix(histogram);
var s0 s3 s12;
run;
/* the variability increases s0 to s3 and then decreases more from s3 to s12 */

ods graphics off;
ods graphics on;
title’data:Correlationstructure’;
proc corr data=lda.aculda (MAXPOINTS=NONE)=matrix(histogram);
var r0 r06 r1 r2 r3 r4 r5 r6 r7 r8 r9 r10;
run;
ods graphics off;


/************ Question 2 ***********/


/*correlation structure */

data lda.acu2;
set lda.acu;
log_time = log(time+1);
run;


ods graphics on; 
title ’data:Correlation structure’;
proc corr data= lda.acu2 plots(MAXPOINTS=NONE)=matrix(histogram);
var s0 s3 s12;
run;
/* the variability increases s0 to s3 and then decreases more from s3 to s12 */

ods graphics off;
ods graphics on;
title’data:Correlationstructure’;
proc corr data=lda.aculda (MAXPOINTS=NONE)=matrix(histogram);
var r0 r06 r1 r2 r3 r4 r5 r6 r7 r8 r9 r10;
run;
ods graphics off;




/*Model 0.0: unstructured mean and covariance matrix- considering log_time as a categorical variable*/
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits0;
class log_time group id;
model severity= group*log_time 
	 / noint s;
repeated log_time / type=un subject=id r rcorr;
run;


/*Model 0: unstructured mean and covariance matrix- considering log_time as a categorical variable*/
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits0;
class log_time group id;
model severity= group*log_time chronicity*log_time frequency*log_time age*log_time
	 / noint s;
repeated log_time / type=un subject=id r rcorr;
run;
/* 7665.1 df 21
/* Two observations:
1. Chronicity*age not significant - so delete chronicity*log_time from model
2. Treating time as categorical variable may not be required and can be treated as continuous variable
*/
data LDA.acu2;
set lda.acu2;
ltimeclss = log_time;
run;
/*Model 1*/
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits1;
class ltimeclss group id;
model severity= group age chronicity frequency group*log_time frequency*log_time age*log_time/ noint s;
repeated ltimeclss / type=un subject=ID r rcorr;
run;
/* LogL 7672.2 df 15
/* chronicity not significant, remove


/*Model 2*/
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits2;
class ltimeclss group id;
model severity= group age frequency group*log_time frequency*log_time age*log_time/ noint s;
repeated ltimeclss / type=un subject=ID r rcorr;
run;

/* LogL: 7673.2 no of parameters =14



/*Model 3*/
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time/  s;
repeated ltimeclss / type=un subject=ID r rcorr;
run;

/* LogL: 7674.5 no of parameters =13
/* this model is not sig diff from unstr as well as model 2, qnd less df SO PREFERABLE*/


/* Even though we saw the lines weren't parallel, we still try to reduce using this
/* try equal slopes for mean structure reduction on the optimum covariates model */

/*model 4 */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits4;
class ltimeclss group id;
model severity= group age frequency log_time / noint s;
repeated ltimeclss / type=un subject=ID r rcorr;
run;
/* LogL 7721.1 df 11
/*  model 3 ( with pchisq(47.9,3 df) ) was significantly different from model 2
Hence we deny the possibility of removing interactions */



/** Covariance structures--CS,TOEP, CSH,TOEPH **/
/*model 3C */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3C;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time/  s;
repeated ltimeclss / type=CS subject=ID r rcorr;
run;
/* 7721.6 df 9
/*model 3T */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3T;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time/  s;
repeated ltimeclss / type=TOEP subject=ID r rcorr;
run;
/* 7721.6 df 10 */
/*model 3CSH */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3CSH;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time/  s;
repeated ltimeclss / type=CSH subject=ID r rcorr;
run;
/* 7675.3 df 11 */
/* 2 CSH seems most appropriate as per LR test */

/*model 3TOEPH */
proc mixed data=lda.acu2 method=ml;
ods output fitstatistics=fits3TOEPH;
class ltimeclss group id;
model severity=  age frequency group*log_time frequency*log_time age*log_time/  s;
repeated ltimeclss / type=TOEPH subject=ID r rcorr;
run;
/* 7675.2 df 12 */








