/*IMPORT DATA FROM FILE*/
/*TAKE THE LOG OF TIME*/


data LDA.ACU;
set LDA.aculda;
log_time=log(time);
run;


/* 1. EXPLORATORY ANALYSIS */

/* 1.0. Descriptive Statistics */
proc freq data = LDA.ACU;
tables time*group/crosslist;
run;



proc means data = LDA.ACU maxdec = 3;
var severity age;
class group;
run;

/*  NO DROPOUT- EQUAL NO OF MEASUREMENTS PER SUBJECT AND MEASUREMENTS TAKEN AT FIXED TIME POINTS =>BALANCED DATA */

/* 1.1.a Individual Profiles Combined*/
PROC SGPLOT data=LDA.ACU noautolegend;
    title "Individual profiles";
    series x=time y=severity / group=id lineattrs=(pattern=solid thickness=2) 
        transparency=0.00 name='Series';
    xaxis grid label='TIME (months)' values=(0 3 12);
    yaxis grid label='Severity'; 
    keylegend 'Series' / title='id:' location=Outside;
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


/* 1.2.a. Graphic MEANS and std dev */
/*code from bharat have nicely distanced x axis*/
PROC SGPLOT data=LDA.ACU;
title 'Mean of individual profiles';
    vline time / response=severity lineattrs=(thickness=2) transparency=0.00 
        stat=mean name='Line' limitstat=stderr;
    xaxis label='TIME (months)' values=(0 3 12);
    yaxis label="Headache Severity" min=10 max=35;
RUN;


/*1.2.b. Graphic MEANS and std dev for each neuro */

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












