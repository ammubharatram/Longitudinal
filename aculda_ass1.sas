/* import and copy dataset */

data lda.acu;
set lda.aculda;
run;


/************ Question 2 ***********/

/*2.1 To obtain mean structure and variance structure and parsimonize them . 
*/

data lda.acu2;
set lda.acu;
timeclss=time;
run;

proc mixed data=lda.acu2 method=ml;
class timeclss group id;
model severity= group*time/ noint s;
repeated timeclss / type=un subject=ID r rcorr;
run;





proc mixed data=lda.acu2 method=ml;
class timeclss ID;
model severity= time age frequency chronicity group group*time;
repeated timeclss / type=un subject=ID r rcorr;
run;
/* Between subject variability is higher than within subject (residual) variability */

/* step 1 : reduce mean structure as chronicity is not significant, remove chronicity and remodel */
proc mixed data=lda.acu2 method=ml;
class timeclss ID;
model severity= time age frequency group group*time;
repeated timeclss / type=un subject=ID r rcorr;
run;
/* As of now, the above mean structure seems parsimonous */


/* step 2: as we see variance components close enough, toep structure could be worth trying */

proc mixed data=lda.acu2;
class timeclss ID;
model severity= time chronicity group group*time;
repeated timeclss / type=toep subject=ID r rcorr;
run;




data test3;
input GROUP x y;

cards;

1 1 260.94
1 3 277.93
1 12 247.65
run;

/*Mixed model initial for fitted variances*/
proc mixed data=lda.acu;
class ID;
model severity= time age chronicity group group*time /solution;
random intercept /type=toep subject=ID g gcorr v vcorr;
run;


/* lets see if topelitz variance structure does better*/

proc mixed data=lda.acu;
class ID;
model severity= time age chronicity group group*time /solution;
random intercept/type=toep subject=ID g gcorr v vcorr;
run;





