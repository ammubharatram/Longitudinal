libname LDA '/folders/myfolders';


data lda.acu2;
set lda.aculda2;
drop severity;
ltime= log(time + 1);
timeclass=time;
run;

proc format ;
value group 0="Basic medication"
	    1="Acupuncture";
run;


/*Target variable is Poisson, don't know if it makes sense to consider Gamma */

/* Quick exploratory analysis */

/* mean structure */
PROC SGPLOT data=LDA.ACU2;
title 'Mean of individual profiles';
    vline time / response=frequency lineattrs=(thickness=2) transparency=0.00 
        stat=mean name='Line' limitstat=stderr;
    xaxis label='TIME (months)' values=(0 3 12);
    yaxis label="Headache Frequency" min=0 max=25;
RUN;


/*without modifying x axis*/

proc gplot data=LDA.ACU2;
plot frequency*time / haxis=axis1 vaxis=axis2;
symbol c=red i=std1mjt w=1 mode=include;
axis1 label=(h=2 "Time (months)") value=(h=1.5) order=(0 to 12 by 1) minor=none;
axis2 label=(h=2 A=90 "Frequency") value=(h=1.5) order=(0 to 20)
minor=none;
title; 
run;quit;

/*Transformation: log(time+1)*/
proc gplot data=LDA.ACU2;
plot frequency*ltime / haxis=axis1 vaxis=axis2;
symbol c=red i=std1mjt w=1 mode=include;
axis1 label=(h=2 "log(1+Time)") value=(h=1.5) order=(0 to 3) minor=none;
axis2 label=(h=2 A=90 "Frequency") value=(h=1.5) order=(0 to 20)
minor=none;
title; 
run;quit;



                                                                                                                                        
/* Calculate the mean and standard error for each X */      
proc sort data=LDA.ACU2;
by group time;
 
proc means data=LDA.ACU2 noprint;                                                                                                           
   by group time;                                                                                                                    
   var Frequency;                                                                                                                            
   output out=meansout(drop=_type_ _freq_) mean=mean stderr=stderr;                                                                     
run;                                                                                                                                    
                                                                                                                                        
                                                                                 
data reshape(keep=group time frequency mean);                                                                                             
   set meansout;                                                                                                                        
   by group time;                                                                                                                    
                                                                                                                                        
/* Offset the X values to display two groups */                                                                                         
   if group=0 then time=time - 0.08;                                                                                               
   else if group=1 then time=time + 0.08;                                                                                          
                                                                                                                                        
   frequency=mean;                                                                                                                           
   output;                                                                                                                              
                                                                                                                                        
   frequency=mean - stderr;                                                                                                                  
   output;                                                                                                                              
                                                                                                                                        
   frequency=mean + stderr;                                                                                                                  
   output;                                                                                                                              
run;                                                                                                                                    
                                                                                                                                        
/* Define the title */                                                                                                                  
   title1 'Means by Group';                                                        
                                                                                                                                        
/* Define the axis characteristics */                                                                                                   
   axis1 offset=(0,0) minor=none value=(t=1 ' ' t=7 ' ');                                                                                                       
   axis2 label=(angle=90) order=(0 to 20 by 5) minor=(n=1);                                                                                                                   
                                                                                                                                        
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
title "Means of frequency by group";
format group group.; 
   plot frequency*time=group / haxis=axis1 vaxis=axis2 legend=legend1;                                                                    
   plot2 mean*time=group / vaxis=axis2 noaxis nolegend;
   
run;          


/*the variance structure */


/* Covariates don't help to reduce the residual variance */ 


proc mixed data=LDA.ACU2; 
class id group time; 
model frequency=group age chronicity time group*time/ solution; 
repeated time / subject=id type=un r rcorr; 
run;
data LDA.test5;input group x1 y2;
cards;
1 0  45.01
1 1.1 71.65  
1 2.5 71.47
;



/* correlation structure */ 

proc sort data=LDA.ACU2;
by id time;

proc transpose data=LDA.ACU2 out=LDA.wide prefix=time;
    by id ;
    id time;
    var frequency;
run;

ods graphics on;
title 'Correlation structure';
proc corr data=LDA.wide plots (MAXPOINTS=NONE)=matrix (histogram);
var time0 time3 time12;
run;
ods graphics off;
