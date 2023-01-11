/*-------------------------------*/
/* IMPORT DATA FROM EXCEL FILE   */
/*-------------------------------*/
%web_drop_table(ery_glp1);


FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_05_10_2021.xlsx';

PROC IMPORT DATAFILE=REFFILE
	DBMS=XLSX
	OUT=ery_glp1;
	GETNAMES=YES;
	SHEET="Hormones_glucose_appetite";
RUN;

PROC CONTENTS DATA=ery_glp1; RUN;


%web_open_table(ery_glp1);

/*----------------------------------------*/
/* RECODE CONDITIONS FOR FACTORIAL DESIGN */
/*----------------------------------------*/
data ery_glp1;
set ery_glp1;
attrib substance format=$10.; /* initializes substance as a categorical variable with length 10 characters */
if condition = "25g D-allulose + 450ppm lactisole" then substance = "allulose";
if condition = "25g D-alluose" then substance = "allulose";
if condition = "300mL tap water + 450ppm lactisole" then substance = "water";
if condition = "300mL tap water" then substance = "water";
if condition = "50g erythritol + 450ppm lactisole" then substance = "erythritol"; 
if condition = "50g erythritol" then substance = "erythritol";
run; 

data ery_glp1;
set ery_glp1;
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
proc sort data=ery_glp1;
by subject condition;
run;

proc univariate data=ery_glp1 noprint;
var GLP_1_pmol_L;
by subject condition;
where time < 0;
output out=baseline_glp1 mean=baseline_glp1;
run;

proc sort data=baseline_glp1;
by subject condition;
run;

data ery_glp1_bl_glp1; 
merge ery_glp1 baseline_glp1; 
by subject condition; 
run;

/*---------------------------------------------------*/
/* CHECK DISTRIBUTION OF BASELINE AND DELTA VARIABLE */
/*---------------------------------------------------*/
proc univariate data=ery_glp1_bl_glp1;
where time = -10;
var baseline_glp1;
histogram baseline_glp1 / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/* */

/* box-cox transformation */
data ery_glp1_bl_glp1;
set ery_glp1_bl_glp1;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_glp1_bl_glp1 maxiter=0 nozeroconstant;
	where time = -10; 
   	model BoxCox(baseline_glp1/parameter=0) = identity(z);
run;
/* check lambda in output, in this case -0.25
parameter is constant to make all values positive if there are negative values, hence parameter = |minimum|, see below */

data ery_glp1_bl_glp1;
set ery_glp1_bl_glp1;
bc_baseline_glp1 = (baseline_glp1**-0.25 -1)/-0.25;
run;
/* boxcox formula, -0.25 is lambda*/

/* check normality of box-cox transformed variable */
proc univariate data=ery_glp1_bl_glp1;
where time = -10;
var bc_baseline_glp1;
histogram bc_baseline_glp1 / normal (mu=est sigma=est);
run;
/*  */

proc univariate data=ery_glp1_bl_glp1;
where time > 0;
var delta_GLP_1;
histogram delta_GLP_1 / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*  */

/* box-cox transformation */
proc transreg data=ery_glp1_bl_glp1 maxiter=0 nozeroconstant;
	where time > 0; 
   	model BoxCox(delta_GLP_1/parameter=18) = identity(z);
run;
/* check lambda in output, in this case 0.75
parameter is constant to make all values positive if there are negative values, hence parameter = |first integer below the minimum| 
here, the minimum value is -17, hence parameter = 18, rendering all values positive */

data ery_glp1_bl_glp1;
set ery_glp1_bl_glp1;
bc_delta_GLP_1 = ((delta_GLP_1+18)**0.75 -1)/0.75;
run;
/* boxcox formula is lambda 0.75  */

/* check normality of box-cox transformed variable */
proc univariate data=ery_glp1_bl_glp1;
where time > 0;
var bc_delta_GLP_1;
histogram bc_delta_GLP_1 / normal (mu=est sigma=est);
run;
/*two outliers left*/


/*-----------------------*/
/* MARGINAL MIXED MODELS */
/*-----------------------*/
/* compare baselines between conditions */
proc mixed data=ery_glp1_bl_glp1;
where time = -10;
class subject condition;
model baseline_glp1 = condition / ddfm=kr2 influence(iter=2) solution residual;
repeated condition / subject=subject type=un r rcorr;
lsmeans condition / diff=all adjust=tukey;
run;
/* */

proc mixed data=ery_glp1_bl_glp1;
where time > 0;
class subject time condition substance lactisole;
model delta_GLP_1 = substance | lactisole | time baseline_glp1 / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
run;


/* we will now include only the update specific hypothesis tests as formulated by Anne & Fabienne in C:\Users\u0027997\Google Drive\LaBGAS\proj-erythritol\Basel\Objective_1\Results\Objective_1_hypotheses_peptides_GE_glucose_insulin.docx
the hypothesis is: 
	1) release of GLP-1 in response to allulose and erythritol 
	2) reduced relase of GLP-1 in response to allulose&lactisole and erythritol&lactisole */
	
proc mixed data=ery_glp1_bl_glp1;
where time > 0;
class subject time condition;
model delta_GLP_1 = condition | time baseline_glp1 / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 0 -1 -1 1 1 divisor=2, /* tests the effect of erythritol versus water */
		 'effect of lactisole within erythritol' 0 0 0 0 -1 1,  /* tests the effect of lactisole within erythritol */
		 'effect of allulose versus water' 1 1 -1 -1 divisor=2, /* tests the effect of allulose versus water */
         'effect of lactisole within allulose' 1 -1 / adjdfe=row adjust=bon stepdown; /* tests the effect of lactisole within allulose */
run;
