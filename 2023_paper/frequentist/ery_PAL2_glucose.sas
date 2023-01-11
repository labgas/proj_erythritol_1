/* Generierter Code (IMPORT) */
/* Quelldatei: SAS_PolyAlluLac_15_02_2022.xlsx */
/* Quellpfad: /home/u50127452/sasuser.v94/PolyAlluLac */
/* Code generiert am: 15.02.22 14:12 */

%web_drop_table(ery_gluc);


FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_15_02_2022.xlsx';

PROC IMPORT DATAFILE=REFFILE
	DBMS=XLSX
	OUT=ery_gluc;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=ery_gluc; RUN;


%web_open_table(ery_gluc);

/*----------------------------------------*/
/* RECODE CONDITIONS FOR FACTORIAL DESIGN */
/*----------------------------------------*/

data ery_gluc;
set ery_gluc;
attrib substance format=$10.; /* initializes substance as a categorical variable with length 10 characters */
if condition = "25g D-alluose" then substance = "allulose";
if condition = "300mL tap water" then substance = "water";
if condition = "50g erythritol" then substance = "erythritol";
run; 

/*-----------------------------------------------------*/
/* CALCULATE BASELINE (average of first 2 timepoints)  */
/*-----------------------------------------------------*/
proc sort data=ery_gluc;
by subject condition;
run;

proc univariate data=ery_gluc noprint;
var glucose_mmol_L;
by subject condition;
where time < 0;
output out=baseline_glucose mean=baseline_glucose;
run;

proc sort data=baseline_glucose;
by subject condition;
run;

data ery_gluc_bl_gluc; 
merge ery_gluc baseline_glucose; 
by subject condition; 
run;

/*---------------------------------------------------*/
/* CHECK DISTRIBUTION OF BASELINE AND DELTA VARIABLE */
/*---------------------------------------------------*/

proc univariate data=ery_gluc_bl_gluc normal;
where time = -10;
var baseline_glucose;
histogram baseline_glucose / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/* looks good */

/* box-cox transformation */
data ery_gluc_bl_gluc;
set ery_gluc_bl_gluc;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_gluc_bl_gluc maxiter=0 nozeroconstant;
	where time = -10; 
   	model BoxCox(baseline_glucose/parameter=0) = identity(z);
run;
/* check lambda in output, in this case -2
parameter is constant to make all values positive if there are negative values, hence parameter = |minimum|, see below */

data ery_gluc_bl_gluc;
set ery_gluc_bl_gluc;
bc_baseline_glucose = (baseline_glucose**-2 -1)/-2;
run;
/* boxcox formula, -2 is lambda, -1 is fixed, add parameter to variable if needed above */

/* check normality of box-cox transformed variable */
proc univariate data=ery_gluc_bl_gluc normal;
where time = -10;
var bc_baseline_glucose;
histogram bc_baseline_glucose / normal (mu=est sigma=est);
run;
/* doesn not make much of a difference - stick with untransformed variable */

proc univariate data=ery_gluc_bl_gluc normal;
where time > 0;
var delta_glucose;
histogram delta_glucose / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/* right outliers are driven by subject with low baseline, hence real datapoints */

/* box-cox transformation */
proc transreg data=ery_gluc_bl_gluc maxiter=0 nozeroconstant;
	where time > 0; 
   	model BoxCox(delta_glucose/parameter=2) = identity(z);
run;
/* check lambda in output, in this case 2
parameter is constant to make all values positive if there are negative values, hence parameter = |first integer below the minimum| 
here, the minimum value is -1.25, hence parameter = 2, rendering all values positive */

data ery_gluc_bl_gluc;
set ery_gluc_bl_gluc;
bc_delta_glucose = ((delta_glucose+2)**2 -1)/2;
run;
/* boxcox formula, 2 is lambda, -1 is fixed, 2 is parameter above */

/* check normality of box-cox transformed variable */
proc univariate data=ery_gluc_bl_gluc normal;
where time > 0;
where delta_glucose > 0;
var bc_delta_glucose;
histogram bc_delta_glucose / normal (mu=est sigma=est);
run;
/* does not make much of a difference - stick with untransformed variable */



/*-----------------------*/
/* MARGINAL MIXED MODELS FOR GLUCOSE */
/*-----------------------*/
/* compare baselines between conditions */
proc mixed data=ery_PAL2_bl_gluc;
where time = -10;
class subject condition;
model baseline_gluc = condition / ddfm=kr2 influence(iter=2) solution residual;
repeated condition / subject=subject type=un r rcorr;
lsmeans condition / diff=all adjust=tukey;
run;
/* no significant baseline differences, but this does not imply you don't have to control for baseline in the subsequent models on delta, since you do want to control for the covariance between the baseline and the subsequent change over time */

proc mixed data=ery_PAL2_bl_gluc;
where time > 0;
class subject time condition substance;
model delta_glucose = substance | time baseline_gluc / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
run;

proc mixed data=ery_gluc_bl_gluc;
where time > 0;
class subject time condition;
model delta_glucose = condition | time baseline_glucose / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 -1 1, /*tests the effect of erythritol versus water*/
		 'effect of allulose versus water' 1 -1 0 / adjdfe=row adjust=bon stepdown; /*tests the effect of allulose versus water*/
run;
