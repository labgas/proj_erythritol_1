/* Generierter Code (IMPORT) */
/* Quelldatei: SAS_PolyAlluLac_15_02_2022.xlsx */
/* Quellpfad: /home/u50127452/sasuser.v94/PolyAlluLac */
/* Code generiert am: 15.02.22 14:12 */

%web_drop_table(ery_insu);


FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_15_02_2022.xlsx';

PROC IMPORT DATAFILE=REFFILE
	DBMS=XLSX
	OUT=ery_insu;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=ery_insu; RUN;


%web_open_table(ery_insu);
/*----------------------------------------*/
/* RECODE CONDITIONS FOR FACTORIAL DESIGN */
/*----------------------------------------*/

data ery_insu;
set ery_insu;
attrib substance format=$10.; /* initializes substance as a categorical variable with length 10 characters */
if condition = "25g D-alluose" then substance = "allulose";
if condition = "300mL tap water" then substance = "water";
if condition = "50g erythritol" then substance = "erythritol";
run; 

/*-----------------------------------------------------*/
/* CALCULATE BASELINE (average of first 2 timepoints)  */
/*-----------------------------------------------------*/
/*Insulin*/

proc sort data=ery_insu;
by subject condition;
run;

proc univariate data=ery_insu noprint;
var insulin_mIU_L;
by subject condition;
where time < 0;
output out=baseline_insulin mean=baseline_insulin;
run;

proc sort data=baseline_insulin;
by subject condition;
run;

data ery_insu_bl_insulin; 
merge ery_insu baseline_insulin; 
by subject condition; 
run;

/*---------------------------------------------------*/
/* CHECK DISTRIBUTION OF BASELINE AND DELTA VARIABLE */
/*---------------------------------------------------*/

proc univariate data=ery_insu_bl_insulin normal;
where time = -10;
var baseline_insulin;
histogram baseline_insulin / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*outlier right*/

/* box-cox transformation */
data ery_insu_bl_insulin;
set ery_insu_bl_insulin;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_insu_bl_insulin maxiter=0 nozeroconstant;
	where time = -10; 
   	model BoxCox(baseline_insulin/parameter=0) = identity(z);
run;
/* check lambda in output, in this case 0 according to google drive code examples use:
 trans_variable = log(variable) ; /* for lambda = 0 */

data ery_insu_bl_insulin;
set ery_insu_bl_insulin;
bc_baseline_insulin = log(baseline_insulin);
run;

/* check normality of box-cox transformed variable */
proc univariate data=ery_insu_bl_insulin normal;
where time = -10;
var bc_baseline_insulin;
histogram bc_baseline_insulin / normal (mu=est sigma=est);
run;
*way better*;

proc univariate data=ery_insu_bl_insulin normal;
where time > 0;
var delta_insulin;
histogram delta_insulin / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;

/* box-cox transformation */
proc transreg data=ery_insu_bl_insulin maxiter=0 nozeroconstant;
	where time > 0; 
   	model BoxCox(delta_insulin/parameter=15) = identity(z);
run;
/* check lambda in output, in this case 1.5
parameter is constant to make all values positive if there are negative values, hence parameter = |first integer below the minimum| 
here, the minimum value is -14, hence parameter = 15, rendering all values positive */

data ery_insu_bl_insulin;
set ery_insu_bl_insulin;
bc_delta_insulin = ((delta_insulin+15)**1.5 -1)/1.5;
run;
/* boxcox formula, 1.5 is lambda, -1 is fixed, 15 is parameter above */

/* check normality of box-cox transformed variable */
proc univariate data=ery_insu_bl_insulin normal;
where time > 0;
var bc_delta_insulin;
histogram bc_delta_insulin / normal (mu=est sigma=est);
run;
*also better*;

/*-----------------------*/
/* MARGINAL MIXED MODELS */
/*-----------------------*/
/* compare baselines between conditions */

proc mixed data=ery_insu_bl_insulin;
where time = -10;
class subject condition;
model baseline_insulin = condition / ddfm=kr2 influence(iter=2) solution residual;
repeated condition / subject=subject type=un r rcorr;
lsmeans condition / diff=all adjust=tukey;
run;

proc mixed data=ery_insu_bl_insulin;
where time > 0;
class subject time condition substance;
model delta_insulin = substance | time baseline_insulin / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
run;


/* we will now include only the update specific hypothesis tests as formulated by Anne & Fabienne in C:\Users\u0027997\Google Drive\LaBGAS\proj-erythritol\Basel\Objective_1\Results\Objective_1_hypotheses_peptides_GE_glucose_insulin.docx*/

proc mixed data=ery_insu_bl_insulin;
where time > 0;
class subject time condition;
model bc_delta_insulin = condition | time baseline_insulin / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 -1 1, /*tests the effect of erythritol versus water*/
		 'effect of allulose versus water' 1 -1 0 / adjdfe=row adjust=bon stepdown; /*tests the effect of allulose versus water*/
run;


