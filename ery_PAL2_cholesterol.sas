/*-------------------------------*/
/* IMPORT DATA FROM EXCEL FILE   */
/*-------------------------------*/
%web_drop_table(ery_BLUA);


FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_14_12_2022.xlsx';

PROC IMPORT DATAFILE=REFFILE
	DBMS=XLSX
	OUT=ery_BLUA;
	GETNAMES=YES;
	SHEET="Blood_lipids_uric_acid";
RUN;

PROC CONTENTS DATA=ery_BLUA; RUN;


%web_open_table(ery_BLUA);


/*---------------------------------------------------*/
/* CHECK DISTRIBUTION OF ABSOLUTE VALUES			 */
/*---------------------------------------------------*/
/*Cholesterol*/
proc univariate data=ery_BLUA normal;
where time > -10;
var cholesterol;
histogram cholesterol / normal (mu=est sigma=est) lognormal (sigma=est theta=est zeta=est);
run;
/*sig. Shapiro and Kolmogorov, also histogram not perfect, outlier right --> transformation */

/* box-cox transformation */
data ery_BLUA;
set ery_BLUA;
z=0;
run;
/* adds variable z with all zeros, needed in proc transreg */

proc transreg data=ery_BLUA maxiter=0 nozeroconstant;
	where time > -10; 
   	model BoxCox(cholesterol/parameter=0) = identity(z);
run;
/* check lambda in output, in this case -0.5
parameter is constant to make all values positive if there are negative values, hence parameter = |minimum|, see below */

data ery_BLUA;
set ery_BLUA;
bc_cholesterol = (cholesterol**-0.5 -1)/-0.5;
run;
/* boxcox formula, -0.5 for lambda*/

/* check normality of box-cox transformed variable */
proc univariate data=ery_BLUA;
where time > -10;
var bc_cholesterol;
histogram bc_cholesterol / normal (mu=est sigma=est);
run;
/*better*/


proc mixed data=ery_BLUA;
where time > -10;
class subject time condition;
model bc_cholesterol = condition | time / ddfm=kr2 influence solution residual;
repeated condition time / subject=subject type=un@ar(1) r rcorr;
lsmeans condition / adjdfe=row;
lsmeans condition*time / slice=time slice=condition adjdfe=row;
lsmestimate condition
		 'effect of erythritol versus water' 0 -1 1, /*tests the effect of erythritol versus water*/
		 'effect of allulose versus water' 1 -1 0 / adjdfe=row adjust=bon stepdown; /*tests the effect of allulose versus water*/
run;
