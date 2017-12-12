*/created library manually as conveneient for webSAS */

data lda2.aculda2;
set lda2.aculda2;
run;


/* Hierarchial Model */

/*Steps for reducing the mean structure- we mean covariates and interactions */

/*Step 1: Poisson Mixed model 1st to identify significant covariate interactions */
proc glimmix data=lda2.aculda2;
class id group ;
model frequency = time group*time chronicity age time*chronicity age*time / noint dist=poisson ;
random intercept  /type=un subject=ID g gcorr v vcorr;
run;
/*remove chronicity and chronicity*time from model


/*step 2 - without chronicity and chronicity*time  */
proc glimmix data=lda2.aculda2;
class id;
model frequency = time group*time age  age*time /noint dist=poisson;
random intercept  /type=un subject=ID g gcorr v vcorr;
run;






/* Poisson mixed model with PQL*/
proc glimmix data=lda2.aculda2 method=rspl;
class id;
model frequency = time group*time age  age*time/noint dist=poisson solution;
random intercept / subject=id g gcorr v vcorr;
run;

