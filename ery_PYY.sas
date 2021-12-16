/*-------------------------------*/
/* IMPORT DATA FROM EXCEL FILE   */
/*-------------------------------*/
%web_drop_table(ery_PYY);


FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_05_10_2021.xlsx';

PROC IMPORT DATAFILE=REFFILE
	DBMS=XLSX
	OUT=ery_PYY;
	GETNAMES=YES;
	SHEET="Hormones_glucose_appetite";
RUN;

PROC CONTENTS DATA=ery_PYY; RUN;


%web_open_table(ery_PYY);

/*----------------------------------------*/
/* RECODE CONDITIONS FOR FACTORIAL DESIGN */
/*----------------------------------------*/
data ery_PYY;
set ery_PYY;
attrib substance format=$10.; /* initializes substance as a categorical variable with length 10 characters */
if condition = "25g D-allulose + 450ppm lactisole" then substance = "allulose";
if condition = "25g D-alluose" then substance = "allulose";
if condition = "300mL tap water + 450ppm lactisole" then substance = "water";
if condition = "300mL tap water" then substance = "water";
if condition = "50g erythritol + 450ppm lactisole" then substance = "erythritol"; 
if condition = "50g erythritol" then substance = "erythritol";
run; 

data ery_PYY;
set ery_PYY;
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
proc sort data=ery_PYY;
by subject condition;
run;

proc univariate data=ery_PYY noprint;
var PYY_pg_mL;
by subject condition;
where time < 0;
output out=baseline_PYY mean=baseline_PYY;
run;

proc sort data=baseline_PYY;
by subject condition;
run;

data ery_PYY_bl_PYY; 
merge ery_PYY baseline_PYY; 
by subject condition; 
run;

/*---------------------------------------------------*/
/* CHECK DISTRIBUTION OF BASELINE AND DELTA VARIABLE */
/*---------------------------------------------------*/
proc univariate data=ery_PYY_bl_PYY;
where time = -10;
var baseline_PYY;
histogram baseline_PYY / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/* */

data ery_PYY_bl_PYY;
set ery_PYY_bl_PYY;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_PYY_bl_PYY maxiter=0 nozeroconstant;
	where time = -10; 
   	model BoxCox(baseline_PYY/parameter=0) = identity(z);
run;
/* check lambda in output, in this case 0.25
parameter is constant to make all values positive if there are negative values, hence parameter = |minimum|, see below */

data ery_PYY_bl_PYY;
set ery_PYY_bl_PYY;
bc_baseline_PYY = (baseline_PYY**0.25 -1)/0.25;
run;
/* boxcox formula, 0.25 is lambda*/

/* check normality of box-cox transformed variable */
proc univariate data=ery_PYY_bl_PYY;
where time = -10;
var bc_baseline_PYY;
histogram bc_baseline_PYY / normal (mu=est sigma=est);
run;
/*  */

proc univariate data=ery_PYY_bl_PYY;
where time > 0;
var delta_PYY;
histogram delta_PYY / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*  */

/* box-cox transformation */
proc transreg data=ery_PYY_bl_PYY maxiter=0 nozeroconstant;
	where time > 0; 
   	model BoxCox(delta_PYY/parameter=111) = identity(z);
run;
/* check lambda in output, in this case 0.5
parameter is constant to make all values positive if there are negative values, hence parameter = |first integer below the minimum| 
here, the minimum value is -105, hence parameter = 111, rendering all values positive */

data ery_PYY_bl_PYY;
set ery_PYY_bl_PYY;
bc_delta_PYY = ((delta_PYY+111)**0.5 -1)/0.5;
run;
/* boxcox formula is lambda 0.5  */

/* check normality of box-cox transformed variable */
proc univariate data=ery_PYY_bl_PYY;
where time > 0;
var bc_delta_PYY;
histogram bc_delta_PYY / normal (mu=est sigma=est);
run;
/* */


/*-----------------------*/
/* MARGINAL MIXED MODELS */
/*-----------------------*/
/* compare baselines between conditions */
proc mixed data=ery_PYY_bl_PYY;
where time = -10;
class subject condition;
model baseline_PYY = condition / ddfm=kr2 influence(iter=2) solution residual;
repeated condition / subject=subject type=un r rcorr;
lsmeans condition / diff=all adjust=tukey;
run;
/* */
/*changed type=un@ar(1) to type=un@cs*/
proc mixed data=ery_PYY_bl_PYY;
where time > 0;
class subject time condition substance;
model delta_PYY = substance | time baseline_PYY / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@cs r rcorr;
run;


/* "Tor style" coding with custom hypothesis testing on the substance*lactisole interaction */
/* the joint option performs an F-test on the three separate t-tests for pairwise contrasts, note that this is exactly the same as the F-test for the main effect of lactisole in the more complex full factorial parametrization of the model above */
/* I use the classic positional syntax here, but check out the nonpositional syntax too */
/* see https://support.sas.com/resources/papers/proceedings11/351-2011.pdf */

/* we will now include only the update specific hypothesis tests as formulated by Anne & Fabienne in C:\Users\u0027997\Google Drive\LaBGAS\proj-erythritol\Basel\Objective_1\Results\Objective_1_hypotheses_peptides_GE_glucose_insulin.docx
the hypothesis is: 
	1) release of PYY in response to allulose, reduced release of PYY in response to allulose&lactisole
	2) relase of PYY in response to erythritol, reduced release of PYY in response to erythritol&lactisole */
	
proc mixed data=ery_PYY_bl_PYY;
where time > 0;
class subject time condition;
model delta_PYY = condition | time baseline_PYY / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@cs r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 0 -1 -1 1 1 divisor=2, /* tests the effect of erythritol versus water */
		 'effect of lactisole within erythritol' 0 0 0 0 -1 1,  /* tests the effect of lactisole within erythritol */
		 'effect of allulose versus water' 1 1 -1 -1 divisor=2, /* tests the effect of allulose versus water */
         'effect of lactisole within allulose' 1 -1 / adjdfe=row adjust=bon stepdown; /* tests the effect of lactisole within allulose */
run;
