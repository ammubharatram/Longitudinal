/*Reduction of mean structure for final model */

proc genmod data=lda.acu2 ;
class id group timeclass;
model frequency = severity age chronicity time severity*time age*time chronicity*time group*time / type3 dist=poisson ;
repeated subject=id / withinsubject=timeclass type=ar   modelse; run;

/* step 1 : remove age*time interaction only as chronicity*time is close to borderline */
proc genmod data=lda.acu2 ;
class id group timeclass;
model frequency = severity age chronicity time severity*time  chronicity*time group*time / type3 dist=poisson ;
repeated subject=id / withinsubject=timeclass type=ar   modelse; run;

/* step 2 : remove chronicity*time now as the p-value has increased */
proc genmod data=lda.acu2 ;
class id group timeclass;
model frequency = severity age chronicity time severity*time group*time / type3 dist=poisson ;
repeated subject=id / withinsubject=timeclass type=ar   modelse; run;

/* we can either remove maine effect chronicity or main effect for group group*time */
proc genmod data=lda.acu2 ;
class id group timeclass;
model frequency = severity age time severity*time group*time / type3 dist=poisson ;
repeated subject=id / withinsubject=timeclass type=ar   modelse; run;


proc genmod data=lda.acu2 ;
class id group timeclass;
model frequency = severity age chronicity time severity*time  / type3 dist=poisson ;
repeated subject=id / withinsubject=timeclass type=ar   modelse; run;

/* neither of them are significant when one is retained and other is removed, hence remove both both chronicity and group*time */
proc genmod data=lda.acu2 ;
class id group timeclass;
model frequency = severity age time severity*time  / type3 dist=poisson ;
repeated subject=id / withinsubject=timeclass type=ar   modelse; run;

/* full: 31287.3156 QIC */
/* reduced: -31365.2765 QIC */

/* full : -31296.5644 QICu */
/* reduced: -31374.0072 QICu */

/* step 1 : remove age*time interaction only as chronicity*time is close to borderline */
/* step 2 : remove chronicity*time now as the p-value has increased */
/* step 3 : we can either remove maine effect chronicity or main effect for group group*time */
/* step 4 : neither of them are significant when one is retained and other is removed, hence remove both both chronicity and group*time */


/* variance structure */

proc genmod data=lda.acu2 ;
class id group timeclass;
model frequency = severity age time severity*time  / type3 dist=poisson ;
repeated subject=id / withinsubject=timeclass type=unstr   modelse; run;