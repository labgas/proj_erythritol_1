/* Generierter Code (IMPORT) */
/* Quelldatei: SAS_PolyAlluLac_15_02_2022.xlsx */
/* Quellpfad: /home/u50127452/sasuser.v94/PolyAlluLac */
/* Code generiert am: 15.02.22 14:12 */

%web_drop_table(ery_insu);


*FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_15_02_2022.xlsx';
FILENAME REFFILE 'G:\.shortcut-targets-by-id\1Vszbf7huRIkgZm1S6pVFaVC5A68CoNOK\proj-erythritol\Basel\Objective_1\MS_PolyAlluLac_2_Review/SAS_PolyAlluLac_15_02_2022.xlsx';

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


/* we will test the specific hypothesis tests as formulated by Anne & Fabienne in C:\Users\u0027997\Google Drive\LaBGAS\proj-erythritol\Basel\Objective_1\Results\Objective_1_hypotheses_peptides_GE_glucose_insulin.docx*/

*FREQUENTIST ANALYSIS
Marginal linear mixed model with Kronecker product of two covariance matrices;
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

*redefine equivalent model using random statement as kronecker product is not available in proc bglimm;
proc mixed data=ery_insu_bl_insulin;
where time > 0;
class subject time condition;
model bc_delta_insulin = condition | time baseline_insulin / ddfm=kr2 influence solution residual;
random int / subject=subject g gcorr;
repeated time / subject=subject(condition) group=condition type=ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 -1 1, /*tests the effect of erythritol versus water*/
		 'effect of allulose versus water' 1 -1 0 / adjdfe=row adjust=bon stepdown; /*tests the effect of allulose versus water*/
run;

proc mixed data=ery_insu_bl_insulin;
where time > 0;
class subject time condition;
model bc_delta_insulin = condition | time baseline_insulin / ddfm=kr2 influence solution residual;
random int / subject=subject(condition) g gcorr;
repeated time / subject=subject(condition) group=condition type=ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 -1 1, /*tests the effect of erythritol versus water*/
		 'effect of allulose versus water' 1 -1 0 / adjdfe=row adjust=bon stepdown; /*tests the effect of allulose versus water*/
run;

proc mixed data=ery_insu_bl_insulin;
where time > 0;
class subject time condition;
model bc_delta_insulin = condition | time baseline_insulin / ddfm=kr2 influence solution residual;
random condition / subject=subject type=un g gcorr;
repeated time / subject=subject(condition) group=condition type=ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 -1 1, /*tests the effect of erythritol versus water*/
		 'effect of allulose versus water' 1 -1 0 / adjdfe=row adjust=bon stepdown; /*tests the effect of allulose versus water*/
run;

proc mixed data=ery_insu_bl_insulin;
where time > 0;
class subject time condition;
model bc_delta_insulin = condition | time baseline_insulin / ddfm=kr2 influence solution residual;
random time / subject=condition(subject) group=condition type=ar(1) g gcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 -1 1, /*tests the effect of erythritol versus water*/
		 'effect of allulose versus water' 1 -1 0 / adjdfe=row adjust=bon stepdown; /*tests the effect of allulose versus water*/
run;

proc mixed data=ery_insu_bl_insulin;
where time > 0;
class subject time condition;
model bc_delta_insulin = condition | time baseline_insulin / ddfm=kr2 influence solution residual;
random condition / subject=subject type=un g gcorr;
random time / subject=subject group=condition type=ar(1) g gcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 -1 1, /*tests the effect of erythritol versus water*/
		 'effect of allulose versus water' 1 -1 0 / adjdfe=row adjust=bon stepdown; /*tests the effect of allulose versus water*/
run;
/* fit (BIC) nor results differ much */

*BAYESIAN ANALYSIS;
proc bglimm data=ery_insu_bl_insulin seed=59621 logpost plots=all nbi=10000 nmc=250000 nthreads=16 dic outpost=insulin_model1;
where time > 0;
class subject time(ref='180') condition;
model bc_delta_insulin = condition | time baseline_insulin / dist=normal cprior=normal(var=2) noint;
random condition / subject=subject type=un g gcorr covprior=uniform(lower=0,upper=1000);
repeated time / subject=subject(condition) group=condition type=ar(1) r rcorr;
estimate "allulose v water at t=15" condition 1 -1 0 time*condition 0 0 0 1 -1 0;
estimate "erythritol v water at t=15" condition 0 -1 1 time*condition 0 0 0 0 -1 1;
estimate "allulose v water at t=30" condition 1 -1 0 time*condition 0 0 0 0 0 0 1 -1 0;
estimate "erythritol v water at t=30" condition 0 -1 1 time*condition 0 0 0 0 0 0 0 -1 1;
estimate "allulose v water at t=45" condition 1 -1 0 time*condition 0 0 0 0 0 0 0 0 0 1 -1 0;
estimate "erythritol v water at t=45" condition 0 -1 1 time*condition 0 0 0 0 0 0 0 0 0 0 -1 1;
estimate "allulose v water at t=60" condition 1 -1 0 time*condition 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0;
estimate "erythritol v water at t=60" condition 0 -1 1 time*condition 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1;
estimate "allulose v water at t=90" condition 1 -1 0 time*condition 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0;
estimate "erythritol v water at t=90" condition 0 -1 1 time*condition 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1;
estimate "allulose v water at t=120" condition 1 -1 0 time*condition 1 -1 0;
estimate "erythritol v water at t=120" condition 0 -1 1 time*condition 0 -1 1;
estimate "allulose v water at t=180" condition 1 -1 0 time*condition 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0;
estimate "erythritol v water at t=180" condition 0 -1 1 time*condition 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1; 
estimate "allulose v water" condition 1 -1 0;
estimate "erythritol v water" condition 0 -1 1;
run;

data work.insulin_model1;
set work.insulin_model1;
rename 'condition 25g D-alluose'n = con_al 'condition 300mL tap water'n = con_wat 'condition 50g erythritol'n = con_ery 'time 120*condition 25g D-alluose'n = al_120 'time 120*condition 300mL tap wat'n = wat_120 'time 120*condition 50g erythrito'n = ery_120 'time 15*condition 25g D-alluose'n = al_15 'time 15*condition 300mL tap wate'n = wat_15 'time 15*condition 50g erythritol'n = ery_15 'time 180*condition 25g D-alluose'n = al_180 'time 180*condition 300mL tap wat'n = wat_180 'time 180*condition 50g erythrito'n = ery_180 'time 30*condition 25g D-alluose'n = al_30 'time 30*condition 300mL tap wate'n = wat_30 'time 30*condition 50g erythritol'n = ery_30 'time 45*condition 25g D-alluose'n = al_45 'time 45*condition 300mL tap wate'n = wat_45 'time 45*condition 50g erythritol'n = ery_45 'time 60*condition 25g D-alluose'n = al_60 'time 60*condition 300mL tap wate'n = wat_60 'time 60*condition 50g erythritol'n = ery_60 'time 90*condition 25g D-alluose'n = al_90 'time 90*condition 300mL tap wate'n = wat_90 'time 90*condition 50g erythritol'n = ery_90 'Residual AR_1,condition 25g D-al'n = res_ar1_al 'Residual AR_1,condition 300mL ta'n = res_ar1_wat 'Residual AR_1,condition 50g eryt'n = res_ar1_ery 'Residual Var,condition 25g D-all'n = res_var_al 'Residual Var,condition 300mL tap'n = res_var_wat 'Residual Var,condition 50g eryth'n = res_var_ery 'Random UN_1_1'n = Random_UN_1_1 'Random UN_2_1'n = Random_UN_2_1 'Random UN_2_2'n = Random_UN_2_2 'Random UN_3_1'n = Random_UN_3_1 'Random UN_3_2'n = Random_UN_3_2 'Random UN_3_3'n = Random_UN_3_3;
run;

%ESS(data=insulin_model1, var = con_al con_wat con_ery al_15 wat_15 ery_15 al_30 wat_30 ery_30 al_45 wat_45 ery_45 al_60 wat_60 ery_60 al_90 wat_90 ery_90 al_120 wat_120 ery_120 al_180 wat_180 ery_180 baseline_insulin res_ar1_al res_ar1_wat res_ar1_ery res_var_al res_var_wat res_var_ery Random_UN_1_1 Random_UN_2_1 Random_UN_2_2 Random_UN_3_1 Random_UN_3_2 Random_UN_3_3);

*cell means form without intercept to check main effect estimation;
proc bglimm data=ery_insu_bl_insulin seed=59621 logpost plots=all nbi=20000 nmc=250000 nthreads=16 dic outpost=insulin_model2;
where time > 0;
class subject time condition;
model bc_delta_insulin = condition*time baseline_insulin / dist=normal cprior=normal(var=2) noint;
random condition / subject=subject type=un g gcorr; *covprior=uniform(lower=0, upper=1000);
repeated time / subject=subject(condition) group=condition type=ar(1) r rcorr;
estimate "allulose v water" time*condition 1 -1 0 1 -1 0 1 -1 0 1 -1 0 1 -1 0 1 -1 0 1 -1 0 / divisor = 7;
estimate "erythritol v water" time*condition 0 -1 1 0 -1 1 0 -1 1 0 -1 1 0 -1 1 0 -1 1 0 -1 1 / divisor = 7;
run;

data work.insulin_model2;
set work.insulin_model2;
rename 'time 120*condition 25g D-alluose'n = al_120 'time 120*condition 300mL tap wat'n = wat_120 'time 120*condition 50g erythrito'n = ery_120 'time 15*condition 25g D-alluose'n = al_15 'time 15*condition 300mL tap wate'n = wat_15 'time 15*condition 50g erythritol'n = ery_15 'time 180*condition 25g D-alluose'n = al_180 'time 180*condition 300mL tap wat'n = wat_180 'time 180*condition 50g erythrito'n = ery_180 'time 30*condition 25g D-alluose'n = al_30 'time 30*condition 300mL tap wate'n = wat_30 'time 30*condition 50g erythritol'n = ery_30 'time 45*condition 25g D-alluose'n = al_45 'time 45*condition 300mL tap wate'n = wat_45 'time 45*condition 50g erythritol'n = ery_45 'time 60*condition 25g D-alluose'n = al_60 'time 60*condition 300mL tap wate'n = wat_60 'time 60*condition 50g erythritol'n = ery_60 'time 90*condition 25g D-alluose'n = al_90 'time 90*condition 300mL tap wate'n = wat_90 'time 90*condition 50g erythritol'n = ery_90 'Residual AR_1,condition 25g D-al'n = res_ar1_al 'Residual AR_1,condition 300mL ta'n = res_ar1_wat 'Residual AR_1,condition 50g eryt'n = res_ar1_ery 'Residual Var,condition 25g D-all'n = res_var_al 'Residual Var,condition 300mL tap'n = res_var_wat 'Residual Var,condition 50g eryth'n = res_var_ery 'Random UN_1_1'n = Random_UN_1_1 'Random UN_2_1'n = Random_UN_2_1 'Random UN_2_2'n = Random_UN_2_2 'Random UN_3_1'n = Random_UN_3_1 'Random UN_3_2'n = Random_UN_3_2 'Random UN_3_3'n = Random_UN_3_3;
run;

%ESS(data=insulin_model2, var = al_15 wat_15 ery_15 al_30 wat_30 ery_30 al_45 wat_45 ery_45 al_60 wat_60 ery_60 al_90 wat_90 ery_90 al_120 wat_120 ery_120 al_180 wat_180 ery_180 baseline_insulin res_ar1_al res_ar1_wat res_ar1_ery res_var_al res_var_wat res_var_ery Random_UN_1_1 Random_UN_2_1 Random_UN_2_2 Random_UN_3_1 Random_UN_3_2 Random_UN_3_3);
