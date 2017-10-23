/*IMPORT DATA FROM FILE*/

libname LDA "C:\Users\Daniel\Desktop\lda";

/*TAKE THE LOG OF TIME*/
data LDA.ACU;
set LDA.aculda;
log_time=log(1+time);
run;


/*Stage 1*/
/*We include group for the second part*/
proc reg data=LDA.ACU outest=lda.coeffs noprint;
by id group;
model severity=time;
run;


/*Stage 2*/

/*Data preparation*/
proc sort data=LDA.ACU nodupkey out=lda.nodups dupout=duplicates; 
by id;
run;
data lda.nodups;
set lda.nodups;
keep id age chronicity group frequency;run; 

/*Create dataset to be used in stage 2*/
data lda.stage2;
merge lda.nodups lda.coeffs;
by id;run;

/*Regression*/
/*For intercepts*/
proc glm data=lda.stage2;
class group;
model intercept=group age frequency chronicity /noint solution;
run;
/*For slopes*/
proc glm data=lda.stage2;
class group;
model time=group age frequency chronicity /noint solution;
run;


/*Mean structure of q2 to be included in stage 2*/
