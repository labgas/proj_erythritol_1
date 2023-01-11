/* AUTHOR: Lukas Van Oudenhove, LaBGAS, KU Leuven; CANlab, Dartmouth College
DATE: November 12th 2020
LOCATION: Dartmouth
FOR: erytrhitol project objective 1 
USEFUL BOOK FOR REFERENCE: http://onlinelibrary.wiley.com/book/10.1002/9781118778210 */

/*-------------------------------*/
/* IMPORT DATA FROM EXCEL FILE   */
/*-------------------------------*/

FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_05_10_2021.xlsx';

PROC IMPORT DATAFILE=REFFILE
	DBMS=XLSX
	OUT=ery_VAS;
	GETNAMES=YES;
	SHEET="Hormones_glucose_appetite";
RUN;

PROC CONTENTS DATA=ery_VAS; RUN;

/*----------------------------------------*/
/* RECODE CONDITIONS FOR FACTORIAL DESIGN */
/*----------------------------------------*/
data ery_VAS;
set ery_VAS;
attrib substance format=$10.; /* initializes substance as a categorical variable with length 10 characters */
if condition = "25g D-allulose + 450ppm lactisole" then substance = "allulose";
if condition = "25g D-alluose" then substance = "allulose";
if condition = "300mL tap water + 450ppm lactisole" then substance = "water";
if condition = "300mL tap water" then substance = "water";
if condition = "50g erythritol + 450ppm lactisole" then substance = "erythritol"; 
if condition = "50g erythritol" then substance = "erythritol";
run; 

data ery_VAS;
set ery_VAS;
attrib lactisole format=$3.;
if condition = "25g D-allulose + 450ppm lactisole" then lactisole = "yes";
if condition = "25g D-alluose" then lactisole = "no";
if condition = "300mL tap water + 450ppm lactisole" then lactisole = "yes";
if condition = "300mL tap water" then lactisole = "no";
if condition = "50g erythritol + 450ppm lactisole" then lactisole = "yes"; 
if condition = "50g erythritol" then lactisole = "no";
run; 

/*-----------------------------------------------------*/
/* CALCULATE BASELINE (average of first 2 timepoints)  */
/*-----------------------------------------------------*/
/*Fullness*/
proc sort data=ery_VAS;
by subject condition;
run;

proc univariate data=ery_VAS noprint;
var fullness_cm;
by subject condition;
where time < 0;
output out=baseline_full mean=baseline_full;
run;

proc sort data=baseline_full;
by subject condition;
run;

data ery_VAS_bl_full; 
merge ery_VAS baseline_full; 
by subject condition; 
run;

/*Hunger*/
proc sort data=ery_VAS;
by subject condition;
run;

proc univariate data=ery_VAS noprint;
var hunger_cm;
by subject condition;
where time < 0;
output out=baseline_hunger mean=baseline_hunger;
run;

proc sort data=baseline_hunger;
by subject condition;
run;

data ery_VAS_bl_hunger; 
merge ery_VAS baseline_hunger; 
by subject condition; 
run;

/*Prospective Food Consumption*/
proc sort data=ery_VAS;
by subject condition;
run;

proc univariate data=ery_VAS noprint;
var prospective_food_consumption_cm;
by subject condition;
where time < 0;
output out=baseline_pfc mean=baseline_pfc;
run;

proc sort data=baseline_pfc;
by subject condition;
run;

data ery_VAS_bl_pfc; 
merge ery_VAS baseline_pfc; 
by subject condition; 
run;

/*Satiety*/
proc sort data=ery_VAS;
by subject condition;
run;

proc univariate data=ery_VAS noprint;
var satiety_cm;
by subject condition;
where time < 0;
output out=baseline_satiety mean=baseline_satiety;
run;

proc sort data=baseline_satiety;
by subject condition;
run;

data ery_VAS_bl_satiety; 
merge ery_VAS baseline_satiety; 
by subject condition; 
run;

/*---------------------------------------------------*/
/* CHECK DISTRIBUTION OF BASELINE AND DELTA VARIABLE */
/*---------------------------------------------------*/
/*Fullness*/
proc univariate data=ery_VAS_bl_full;
where time = -10;
var baseline_full;
histogram baseline_full / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*  */

/* box-cox transformation */
data ery_VAS_bl_full;
set ery_VAS_bl_full;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_VAS_bl_full maxiter=0 nozeroconstant;
	where time = -10; 
   	model BoxCox(baseline_full/parameter=1) = identity(z);
run;
/* check lambda in output, in this case 0.5
parameter is constant to make all values positive if there are negative values, hence parameter = |minimum|, see below */

data ery_VAS_bl_full;
set ery_VAS_bl_full;
bc_baseline_full = (baseline_full**0.5 -1)/0.5;
run;
/* boxcox formula, 0.5 is lambda, google drive: lambda 0.5 --> use sqrt(variable) ??? */

/* check normality of box-cox transformed variable */
proc univariate data=ery_VAS_bl_full;
where time = -10;
var bc_baseline_full;
histogram bc_baseline_full / normal (mu=est sigma=est);
run;
/*  */

proc univariate data=ery_VAS_bl_full;
where time > 0;
var delta_fullness;
histogram delta_fullness / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*  */

/* box-cox transformation */
proc transreg data=ery_VAS_bl_full maxiter=0 nozeroconstant;
	where time > 0; 
   	model BoxCox(delta_fullness/parameter=7) = identity(z);
run;
/* check lambda in output, in this case 1
parameter is constant to make all values positive if there are negative values, hence parameter = |first integer below the minimum| 
here, the minimum value is -6.5, hence parameter = 7, rendering all values positive */

data ery_VAS_bl_full;
set ery_VAS_bl_full;
bc_delta_fullness = ((delta_fullness+7)**1 -1)/1;
run;
/* boxcox formula, 1 is lambda, */

/* check normality of box-cox transformed variable */
proc univariate data=ery_VAS_bl_full;
where time > 0;
var bc_delta_fullness;
histogram bc_delta_fullness / normal (mu=est sigma=est);
run;
/* */

/*Hunger*/
proc univariate data=ery_VAS_bl_hunger;
where time = -10;
var baseline_hunger;
histogram baseline_hunger / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*  */

/* box-cox transformation */
data ery_VAS_bl_hunger;
set ery_VAS_bl_hunger;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_VAS_bl_hunger maxiter=0 nozeroconstant;
	where time = -10; 
   	model BoxCox(baseline_hunger/parameter=1) = identity(z);
run;
/* check lambda in output, in this case 1
parameter is constant to make all values positive if there are negative values, hence parameter = |minimum|, see below */

data ery_VAS_bl_hunger;
set ery_VAS_bl_hunger;
bc_baseline_hunger = (baseline_hunger**1 -1)/1;
run;
/* boxcox formula, 0.5 is lambda, google drive: lambda 1 --> variable**1 ??? */

/* check normality of box-cox transformed variable */
proc univariate data=ery_VAS_bl_hunger;
where time = -10;
var bc_baseline_hunger;
histogram bc_baseline_hunger / normal (mu=est sigma=est);
run;
/*  */

proc univariate data=ery_VAS_bl_hunger;
where time > 0;
var delta_hunger;
histogram delta_hunger / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*  */

/* box-cox transformation */
proc transreg data=ery_VAS_bl_hunger maxiter=0 nozeroconstant;
	where time > 0; 
   	model BoxCox(delta_hunger/parameter=6) = identity(z);
run;
/* check lambda in output, in this case 0.75
parameter is constant to make all values positive if there are negative values, hence parameter = |first integer below the minimum| 
here, the minimum value is -5, hence parameter = 6, rendering all values positive */

data ery_VAS_bl_hunger;
set ery_VAS_bl_hunger;
bc_delta_hunger = ((delta_hunger+6)**0.75 -1)/0.75;
run;
/* boxcox formula, 0.75 is lambda, parameter = 6 */

/* check normality of box-cox transformed variable */
proc univariate data=ery_VAS_bl_hunger;
where time > 0;
var bc_delta_hunger;
histogram bc_delta_hunger / normal (mu=est sigma=est);
run;
/* */

/*Prospective Food Consumption*/
proc univariate data=ery_VAS_bl_pfc;
where time = -10;
var baseline_pfc;
histogram baseline_pfc / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*  */

/* box-cox transformation */
data ery_VAS_bl_pfc;
set ery_VAS_bl_pfc;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_VAS_bl_pfc maxiter=0 nozeroconstant;
	where time = -10; 
   	model BoxCox(baseline_pfc/parameter=1) = identity(z);
run;
/* check lambda in output, in this case 1.5
parameter is constant to make all values positive if there are negative values, hence parameter = |minimum|, see below */

data ery_VAS_bl_pfc;
set ery_VAS_bl_pfc;
bc_baseline_pfc = (baseline_pfc**1.5 -1)/1.5;
run;
/* boxcox formula, 1.5 is lambda*/

/* check normality of box-cox transformed variable */
proc univariate data=ery_VAS_bl_pfc;
where time = -10;
var bc_baseline_pfc;
histogram bc_baseline_pfc / normal (mu=est sigma=est);
run;
/*  */

proc univariate data=ery_VAS_bl_pfc;
where time > 0;
var delta_pfc;
histogram delta_pfc / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*  */

/* box-cox transformation */
proc transreg data=ery_VAS_bl_pfc maxiter=0 nozeroconstant;
	where time > 0; 
   	model BoxCox(delta_pfc/parameter=8) = identity(z);
run;
/* check lambda in output, in this case 1
parameter is constant to make all values positive if there are negative values, hence parameter = |first integer below the minimum| 
here, the minimum value is -7.5, hence parameter = 8, rendering all values positive */

data ery_VAS_bl_pfc;
set ery_VAS_bl_pfc;
bc_delta_pfc = ((delta_pfc+8)**1 -1)/1;
run;
/* boxcox formula, 1 is lambda, */

/* check normality of box-cox transformed variable */
proc univariate data=ery_VAS_bl_pfc;
where time > 0;
var bc_delta_pfc;
histogram bc_delta_pfc / normal (mu=est sigma=est);
run;
/* */

/*Satiety*/
proc univariate data=ery_VAS_bl_satiety;
where time = -10;
var baseline_satiety;
histogram baseline_satiety / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*  */

/* box-cox transformation */
data ery_VAS_bl_satiety;
set ery_VAS_bl_satiety;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_VAS_bl_satiety maxiter=0 nozeroconstant;
	where time = -10; 
   	model BoxCox(baseline_satiety/parameter=1) = identity(z);
run;
/* check lambda in output, in this case 0.25
parameter is constant to make all values positive if there are negative values, hence parameter = |minimum|, see below */

data ery_VAS_bl_satiety;
set ery_VAS_bl_satiety;
bc_baseline_satiety = (baseline_satiety**0.25 -1)/0.25;
run;
/* boxcox formula, 0.25 is lambda*/

/* check normality of box-cox transformed variable */
proc univariate data=ery_VAS_bl_satiety;
where time = -10;
var bc_baseline_satiety;
histogram bc_baseline_satiety / normal (mu=est sigma=est);
run;
/* outlier left */

proc univariate data=ery_VAS_bl_satiety;
where time > 0;
var delta_satiety;
histogram delta_satiety / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*  */

/* box-cox transformation */
proc transreg data=ery_VAS_bl_satiety maxiter=0 nozeroconstant;
	where time > 0; 
   	model BoxCox(delta_pfc/parameter=8) = identity(z);
run;
/* check lambda in output, in this case 1
parameter is constant to make all values positive if there are negative values, hence parameter = |first integer below the minimum| 
here, the minimum value is -7.2, hence parameter = 8, rendering all values positive */

data ery_VAS_bl_satiety;
set ery_VAS_bl_satiety;
bc_delta_satiety = ((delta_satiety+8)**1 -1)/1;
run;
/* boxcox formula, 1 is lambda, formula trans_variable = variable**1 ; /* for lambda = 1 */

/* check normality of box-cox transformed variable */

proc univariate data=ery_VAS_bl_satiety;
where time > 0;
var bc_delta_satiety;
histogram bc_delta_satiety / normal (mu=est sigma=est);
run;
/* */


/*-----------------------*/
/* MARGINAL MIXED MODELS */
/*-----------------------*/
/*Fullness*/
/* compare baselines between condition*/

proc mixed data=ery_VAS_bl_full;
where time = -10;
class subject condition;
model baseline_full = condition / ddfm=kr2 influence(iter=2) solution residual;
repeated condition / subject=subject type=un r rcorr;
lsmeans condition / diff=all adjust=tukey;
run;

/* no significant baseline differences, but this does not imply you don't have to control for baseline in the subsequent models on delta, since you do want to control for the covariance between the baseline and the subsequent change over time*/
proc mixed data=ery_VAS_bl_full;
where time > 0;
class subject time condition substance lactisole;
model delta_fullness = substance | lactisole | time baseline_full / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
run;

/* "Tor style" coding with custom hypothesis testing on the substance*lactisole interaction */
/* the joint option performs an F-test on the three separate t-tests for pairwise contrasts, note that this is exactly the same as the F-test for the main effect of lactisole in the more complex full factorial parametrization of the model above */
/* I use the classic positional syntax here, but check out the nonpositional syntax too */
/* see https://support.sas.com/resources/papers/proceedings11/351-2011.pdf */

/* we will now include only the update specific hypothesis tests as formulated by Anne & Fabienne on google drive*/
proc mixed data=ery_VAS_bl_full;
where time > 0;
class subject time condition;
model delta_fullness = condition | time baseline_full / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 0 -1 -1 1 1 divisor=2, /* tests the effect of erythritol versus water */
		 'effect of lactisole within erythritol' 0 0 0 0 -1 1,  /* tests the effect of lactisole within erythritol */
		 'effect of allulose versus water' 1 1 -1 -1 divisor=2, /* tests the effect of allulose versus water */
         'effect of lactisole within allulose' 1 -1 / adjdfe=row adjust=bon stepdown; /* tests the effect of lactisole within allulose */
run;


/*Hunger*/
/* compare baselines between conditions */
proc mixed data=ery_VAS_bl_hunger;
where time = -10;
class subject condition;
model baseline_hunger = condition / ddfm=kr2 influence(iter=2) solution residual;
repeated condition / subject=subject type=un r rcorr;
lsmeans condition / diff=all adjust=tukey;
run;

/* */

proc mixed data=ery_VAS_bl_hunger;
where time > 0;
class subject time condition substance lactisole;
model delta_hunger = substance | lactisole | time baseline_hunger / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
run;


/* "Tor style" coding with custom hypothesis testing on the substance*lactisole interaction */
/* the joint option performs an F-test on the three separate t-tests for pairwise contrasts, note that this is exactly the same as the F-test for the main effect of lactisole in the more complex full factorial parametrization of the model above */
/* I use the classic positional syntax here, but check out the nonpositional syntax too */
/* see https://support.sas.com/resources/papers/proceedings11/351-2011.pdf */

/* we will now include only the update specific hypothesis tests as formulated by Anne & Fabienne on google drive*/
proc mixed data=ery_VAS_bl_hunger;
where time > 0;
class subject time condition;
model delta_hunger = condition | time baseline_hunger / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 0 -1 -1 1 1 divisor=2, /* tests the effect of erythritol versus water */
		 'effect of lactisole within erythritol' 0 0 0 0 -1 1,  /* tests the effect of lactisole within erythritol */
		 'effect of allulose versus water' 1 1 -1 -1 divisor=2, /* tests the effect of allulose versus water */
         'effect of lactisole within allulose' 1 -1 / adjdfe=row adjust=bon stepdown; /* tests the effect of lactisole within allulose */
run;

/*Prospective Food Consumption*/
/* compare baselines between conditions */
proc mixed data=ery_VAS_bl_pfc;
where time = -10;
class subject condition;
model baseline_pfc = condition / ddfm=kr2 influence(iter=2) solution residual;
repeated condition / subject=subject type=un r rcorr;
lsmeans condition / diff=all adjust=tukey;
run;
/* */

proc mixed data=ery_VAS_bl_pfc;
where time > 0;
class subject time condition substance lactisole;
model delta_pfc = substance | lactisole | time baseline_pfc / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
run;


/* "Tor style" coding with custom hypothesis testing on the substance*lactisole interaction */
/* the joint option performs an F-test on the three separate t-tests for pairwise contrasts, note that this is exactly the same as the F-test for the main effect of lactisole in the more complex full factorial parametrization of the model above */
/* I use the classic positional syntax here, but check out the nonpositional syntax too */
/* see https://support.sas.com/resources/papers/proceedings11/351-2011.pdf */

/* we will now include only the update specific hypothesis tests as formulated by Anne & Fabienne on google drive */
proc mixed data=ery_VAS_bl_pfc;
where time > 0;
class subject time condition;
model delta_pfc = condition | time baseline_pfc / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
lsmeans condition / diff=all adjust=tukey adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 0 -1 -1 1 1 divisor=2, /* tests the effect of erythritol versus water */
		 'effect of lactisole within erythritol' 0 0 0 0 -1 1,  /* tests the effect of lactisole within erythritol */
		 'effect of allulose versus water' 1 1 -1 -1 divisor=2, /* tests the effect of allulose versus water */
         'effect of lactisole within allulose' 1 -1 / adjdfe=row adjust=bon stepdown; /* tests the effect of lactisole within allulose */
run;

/*Satiety*/
/* compare baselines between conditions */

proc mixed data=ery_VAS_bl_satiety;
where time = -10;
class subject condition;
model baseline_satiety = condition / ddfm=kr2 influence(iter=2) solution residual;
repeated condition / subject=subject type=un r rcorr;
lsmeans condition / diff=all adjust=tukey;
run;
/* */

proc mixed data=ery_VAS_bl_satiety;
where time > 0;
class subject time condition substance lactisole;
model delta_satiety = substance | lactisole | time baseline_satiety / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
run;

/* "Tor style" coding with custom hypothesis testing on the substance*lactisole interaction */
/* the joint option performs an F-test on the three separate t-tests for pairwise contrasts, note that this is exactly the same as the F-test for the main effect of lactisole in the more complex full factorial parametrization of the model above */
/* I use the classic positional syntax here, but check out the nonpositional syntax too */
/* see https://support.sas.com/resources/papers/proceedings11/351-2011.pdf */

/* we will now include only the update specific hypothesis tests as formulated by Anne & Fabienne on google drive*/
proc mixed data=ery_VAS_bl_satiety;
where time > 0;
class subject time condition;
model delta_satiety = condition | time baseline_satiety / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 0 -1 -1 1 1 divisor=2, /* tests the effect of erythritol versus water */
		 'effect of lactisole within erythritol' 0 0 0 0 -1 1,  /* tests the effect of lactisole within erythritol */
		 'effect of allulose versus water' 1 1 -1 -1 divisor=2, /* tests the effect of allulose versus water */
         'effect of lactisole within allulose' 1 -1 / adjdfe=row adjust=bon stepdown; /* tests the effect of lactisole within allulose */
run;