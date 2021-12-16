/* AUTHOR: Lukas Van Oudenhove, LaBGAS, KU Leuven; CANlab, Dartmouth College
DATE: November 12th 2020
LOCATION: Dartmouth
FOR: erytrhitol project objective 1 
USEFUL BOOK FOR REFERENCE: http://onlinelibrary.wiley.com/book/10.1002/9781118778210 */

/*-------------------------------*/
/* IMPORT DATA FROM EXCEL FILE   */
/*-------------------------------*/
/*%let path=C:\Users\u0027997\Google Drive\LaBGAS\proj-erythritol\Basel\Objective_1\Results;
proc import
datafile="&path\SAS_PolyAlluLac_20201117.xlsx"
out=work.ery_obj1;
sheet="hormones_glucose_appetite";
run;*/

%web_drop_table(ery_ge);


FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_05_10_2021.xlsx';

PROC IMPORT DATAFILE=REFFILE
	DBMS=XLSX
	OUT=ery_ge;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=ery_ge; RUN;


%web_open_table(ery_ge);


/*----------------------------------------*/
/* RECODE CONDITIONS FOR FACTORIAL DESIGN */
/*----------------------------------------*/
data ery_ge;
set ery_ge;
attrib substance format=$10.; /* initializes substance as a categorical variable with length 10 characters */
if condition = "25g D-allulose + 450ppm lactisole" then substance = "allulose";
if condition = "25g D-alluose" then substance = "allulose";
if condition = "300mL tap water + 450ppm lactisole" then substance = "water";
if condition = "300mL tap water" then substance = "water";
if condition = "50g erythritol + 450ppm lactisole" then substance = "erythritol"; 
if condition = "50g erythritol" then substance = "erythritol";
run; 

data ery_ge;
set ery_ge;
attrib lactisole format=$3.;
if condition = "25g D-allulose + 450ppm lactisole" then lactisole = "yes";
if condition = "25g D-alluose" then lactisole = "no";
if condition = "300mL tap water + 450ppm lactisole" then lactisole = "yes";
if condition = "300mL tap water" then lactisole = "no";
if condition = "50g erythritol + 450ppm lactisole" then lactisole = "yes"; 
if condition = "50g erythritol" then lactisole = "no";
run; 

/*-----------------------------------------------------*/
/* Baseline */
/*-----------------------------------------------------*/
proc sort data=ery_ge;
by subject condition;
run;

proc univariate data=ery_ge noprint;
var GE_dosis;
by subject condition;
where time > 0;
output out=baseline_ge;
run;

proc sort data=baseline_ge;
by subject condition;
run;

data ery_ge_bl_ge; 
merge ery_ge baseline_ge; 
by subject condition; 
run;


/*---------------------------------------------------*/
/* CHECK DISTRIBUTION OF VARIABLE */
/*---------------------------------------------------*/

proc univariate data=ery_ge_bl_ge;
where time > 0;
var GE_dosis;
histogram GE_dosis / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/**/

/* box-cox transformation */
data ery_ge_bl_ge;
set ery_ge_bl_ge;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_ge_bl_ge maxiter=0 nozeroconstant;
	where time > 0; 
   	model BoxCox(GE_dosis/parameter=1) = identity(z);
run;
/* check lambda in output, in this case 0.25
parameter is constant to make all values positive if there are negative values, hence parameter = |minimum|, see below */

data ery_ge_bl_ge;
set ery_ge_bl_ge;
bc_GE_dosis = ((GE_dosis)**0.25 -1)/0.25;
run;
/* boxcox formula, 0.25 is lambda, -1 is fixed, add parameter to variable if needed above */

/* check normality of box-cox transformed variable */
proc univariate data=ery_ge_bl_ge;
where time > 0;
var bc_GE_dosis;
histogram bc_GE_dosis / normal (mu=est sigma=est);
run;
/* outlier left */


/*-----------------------*/
/* MARGINAL MIXED MODELS */
/*-----------------------*/
/* compare baselines between conditions */
proc mixed data=ery_ge_bl_ge;
where time > 0;
class subject condition;
model GE_dosis = condition / ddfm=kr2 influence(iter=2) solution residual;
repeated condition / subject=subject type=un r rcorr;
lsmeans condition / diff=all adjust=tukey;
run;

proc mixed data=ery_ge_bl_ge;
where time > 0;
class subject time condition substance lactisole;
model bc_GE_dosis = substance | lactisole | time GE_dosis / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
run;

/* we will now include only the update specific hypothesis tests as formulated by Anne & Fabienne in C:\Users\u0027997\Google Drive\LaBGAS\proj-erythritol\Basel\Objective_1\Results\Objective_1_hypotheses_peptides_GE_glucose_insulin.docx */
/* The hypothesis is: Retardation in response to Erythritol and D-Allulose, no effect of lactisole*/
/*Retardation in response to erythritol, no effect of lactisole.*/
/*Retardation in response to D-allulose, no effect of lactisole.*/
;
proc mixed data=ery_ge_bl_ge;
where time > 0 AND time < 61; 
class subject time condition;
model bc_GE_dosis = condition | time / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 0 -1 -1 1 1 divisor=2, /* tests the effect of erythritol versus water */
		 'effect of lactisole within erythritol' 0 0 0 0 -1 1,  /* tests the effect of lactisole within erythritol */
		 'effect of allulose versus water' 1 1 -1 -1 divisor=2, /* tests the effect of allulose versus water */
         'effect of lactisole within allulose' 1 -1 / adjdfe=row adjust=bon stepdown; /* tests the effect of lactisole within allulose */
lsmestimate condition*time
		'effect of erythritol versus water at time 15' 0 0 -1 -1 1 1 divisor=2;
run;

