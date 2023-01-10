/*-------------------------------*/
/* IMPORT DATA FROM EXCEL FILE   */
/*-------------------------------*/
%web_drop_table(ery_BLUA);


FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_05_10_2021.xlsx';

PROC IMPORT DATAFILE=REFFILE
	DBMS=XLSX
	OUT=ery_BLUA;
	GETNAMES=YES;
	SHEET="Blood_lipids_uric_acid";
RUN;

PROC CONTENTS DATA=ery_BLUA; RUN;


%web_open_table(ery_BLUA);


/*---------------------------------------------------*/
/* CHECK DISTRIBUTION OF ABSOLUTE VARIABLES			 */
/*---------------------------------------------------*/
/*hsCRP*/
proc univariate data=ery_BLUA normal;
where time > -11;
var hsCRP_mg_L;
histogram hsCRP_mg_L / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*  */

/* box-cox transformation */
data ery_BLUA;
set ery_BLUA;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_BLUA maxiter=0 nozeroconstant;
where time > -11; 
model BoxCox(hsCRP_mg_L/parameter=0) = identity(z);
run;
/* check lambdarun;
/* adds variable z with all zeros, needed in proc transreg */

data ery_BLUA;
set ery_BLUA;
bc_hsCRP_mg_L = (hsCRP_mg_L**-0.25 -1)/-0.25;
run;
/* boxcox formula, -0.25 is lambda

/* check normality of box-cox transformed variable */
proc univariate data=ery_BLUA;
where time > -11;
var bc_hsCRP_mg_L;
histogram bc_hsCRP_mg_L/ normal (mu=est sigma=est);
run; 


/*-----------------------*/
/* MARGINAL MIXED MODELS */
/*-----------------------*/
/*hsCRP*/
proc mixed data=ery_BLUA;
where time > -11;
class subject time condition;
model bc_hsCRP_mg_L = condition | time / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 -1 1, /*tests the effect of erythritol versus water*/
		 'effect of allulose versus water' 1 -1 0 / adjdfe=row adjust=bon stepdown; /*tests the effect of allulose versus water*/
run;

