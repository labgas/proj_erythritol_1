/* Generierter Code (IMPORT) */
/* Quelldatei: SAS_PolyAlluLac_15_11_2022.xlsx */
/* Quellpfad: /home/u50127452/sasuser.v94/PolyAlluLac */
/* Code generiert am: 15.11.22 12:56 */

%web_drop_table(ery_ghrelin);


FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_15_11_2022.xlsx';


PROC IMPORT DATAFILE=REFFILE
	DBMS=XLSX
	OUT=ery_ghrelin;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=ery_ghrelin; RUN;


%web_open_table(ery_ghrelin);


/*----------------------------------------*/
/* RECODE CONDITIONS FOR FACTORIAL DESIGN */
/*----------------------------------------*/

data ery_ghrelin;
set ery_ghrelin;
attrib substance format=$10.; /* initializes substance as a categorical variable with length 10 characters */
if condition = "25g D-alluose" then substance = "allulose";
if condition = "300mL tap water" then substance = "water";
if condition = "50g erythritol" then substance = "erythritol";
run; 

/*-----------------------------------------------------*/
/* CALCULATE BASELINE (average of first 2 timepoints)  */
/*-----------------------------------------------------*/
proc sort data=ery_ghrelin;
by subject condition;
run;

proc univariate data=ery_ghrelin noprint;
var ghrelin_pg_mL;
by subject condition;
where time < 0;
output out=baseline_ghrelin mean=baseline_ghrelin;
run;

proc sort data=baseline_ghrelin;
by subject condition;
run;

data ery_ghrelin_bl_ghrelin; 
merge ery_ghrelin baseline_ghrelin; 
by subject condition; 
run;

/*---------------------------------------------------*/
/* CHECK DISTRIBUTION OF BASELINE AND DELTA VARIABLE */
/*---------------------------------------------------*/

proc univariate data=ery_ghrelin_bl_ghrelin normal;
where time = -10;
var baseline_ghrelin;
histogram baseline_ghrelin / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/* <0.010 */

/* box-cox transformation */
data ery_ghrelin_bl_ghrelin;
set ery_ghrelin_bl_ghrelin;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_ghrelin_bl_ghrelin maxiter=0 nozeroconstant;
	where time = -10; 
   	model BoxCox(baseline_ghrelin/parameter=0) = identity(z);
run;
/* check lambda in output, in this case -0.25
parameter is constant to make all values positive if there are negative values, hence parameter = |minimum|, see below */

data ery_ghrelin_bl_ghrelin;
set ery_ghrelin_bl_ghrelin;
bc_baseline_ghrelin = (baseline_ghrelin**-0.25 -1)/-0.25;
run;
/* boxcox formula, -0.25 is lambda, -1 is fixed, add parameter to variable if needed above */

/* check normality of box-cox transformed variable */
proc univariate data=ery_ghrelin_bl_ghrelin normal;
where time = -10;
var bc_baseline_ghrelin;
histogram bc_baseline_ghrelin / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/* much better */

proc univariate data=ery_ghrelin_bl_ghrelin normal;
where time > 0;
var delta_ghrelin;
histogram delta_ghrelin / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/* looks not too bad */

/* box-cox transformation */
proc transreg data=ery_ghrelin_bl_ghrelin maxiter=0 nozeroconstant;
	where time > 0; 
   	model BoxCox(delta_ghrelin/parameter=195) = identity(z);
run;
/* check lambda in output, in this case 0.75
parameter is constant to make all values positive if there are negative values, hence parameter = |first integer below the minimum| 
here, the minimum value is -194.07089 hence parameter = 195, rendering all values positive */

data ery_ghrelin_bl_ghrelin;
set ery_ghrelin_bl_ghrelin;
bc_delta_ghrelin = ((delta_ghrelin+195)**0.75 -1)/0.75;
run;
/* boxcox formula, 0.75 is lambda, -1 is fixed, 195 is parameter above */

/* check normality of box-cox transformed variable */
proc univariate data=ery_ghrelin_bl_ghrelin normal;
where time > 0;
where delta_ghrelin > 0;
var bc_delta_ghrelin;
histogram bc_delta_ghrelin / normal (mu=est sigma=est);
run;
/* worse - stick with untransformed variable */


/*
data ery_ghrelin_bl_ghrelin;
set ery_ghrelin_bl_ghrelin;
log_delta_ghrelin = log(delta_ghrelin+195);
run;


proc univariate data=ery_ghrelin_bl_ghrelin normal;
where time > 0;
where delta_ghrelin > 0;
var log_delta_ghrelin;
histogram log_delta_ghrelin / normal (mu=est sigma=est);
run;

/* b


/*-----------------------*/
/* MARGINAL MIXED MODELS FOR Ghrelin */
/*-----------------------*/
/* compare baselines between conditions */
proc mixed data=ery_ghrelin_bl_ghrelin;
where time = -10;
class subject condition;
model baseline_ghrelin = condition / ddfm=kr2 influence(iter=2) solution residual;
repeated condition / subject=subject type=un r rcorr;
lsmeans condition / diff=all adjust=tukey;
run;
/* no significant baseline differences, but this does not imply you don't have to control for baseline in the subsequent models on delta, since you do want to control for the covariance between the baseline and the subsequent change over time */

proc mixed data=ery_ghrelin_bl_ghrelin;
where time > 0;
class subject time condition substance;
model delta_ghrelin = substance | time bc_baseline_ghrelin/ ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
run;


/*------------------------------------------------------------------------------------*/
/* Exploration of time point 30min post D-allulose and erythritol consumption */
/* Study by Sorrentino et al. 2020 found no effect on total ghrelin AUC, however, */
/* specific time points in response to erythritol compared to aspartame were significant*/
/*----------------------------------------------------------------------------------------------------*/


proc mixed data=ery_ghrelin_bl_ghrelin;
where time > 0;
class subject time condition;
model delta_ghrelin = condition | time bc_baseline_ghrelin / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 -1 1, /*tests the effect of erythritol versus water*/
		 'effect of allulose versus water' 1 -1 0 / adjdfe=row adjust=bon stepdown; /*tests the effect of allulose versus water*/
lsmestimate time*condition
		 'effect of erythritol versus water at time point 30' 0 0 0 0 -1 1, /*tests the effect of erythritol versus water at time point 30*/
		 'effect of allulose versus water at time point 30' 0 0 0 1 -1 0 / adjdfe=row adjust=bon stepdown; /*tests the effect of allulose versus water at time point 30*/
run;




