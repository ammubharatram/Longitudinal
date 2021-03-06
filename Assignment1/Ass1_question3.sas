/*IMPORT DATA FROM FILE*/

libname LDA "C:\Users\Daniel\Desktop\lda";


/*TAKE THE LOG OF TIME*/
data LDA.ACU;
set LDA.aculda;
log_time=log(1+time);
run;

data lda.acu2;
set lda.acu;
ltimeclss=log_time;
run;


/*Stage 1*/
/*We include group for the second part*/
proc reg data=LDA.ACU outest=lda.coeffs noprint;
by id group;
model severity=log_time;
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
model intercept=group age frequency/noint solution;
run;

/*delete subjects with only one observation to work with the slopes*/
data lda.stage2c;
set lda.stage2;
if log_time eq 0 then delete;
run;
/*For slopes*/
proc glm data=lda.stage2c;
class group;
model log_time=group age frequency/noint solution;
run;

/*Plot histogram of intercepts and slopes of the subjects*/
ods graphics off;
proc univariate data=lda.stage2 noprint;
   histogram intercept log_time;
run;

/*Plot histogram of individual R2*/
proc reg data=LDA.ACU outest=lda.coeffs2 noprint;
by id group;
model severity=log_time  4/ selection= rsquare;
run;
proc univariate data=lda.coeffs2 noprint;
histogram _RSQ_;
run;
