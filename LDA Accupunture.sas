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
/*code from bharat have nicely distanced x axis*/
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


