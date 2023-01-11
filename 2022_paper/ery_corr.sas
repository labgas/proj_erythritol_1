/*-------------------------------*/
/* IMPORT DATA FROM EXCEL FILE   */
/*-------------------------------*/
%web_drop_table(ery_corr);


FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_10_12_2021.xlsx';

PROC IMPORT DATAFILE=REFFILE
	DBMS=XLSX
	OUT=ery_corr;
	GETNAMES=YES;
	SHEET="Associations_hormones";
RUN;

PROC CONTENTS DATA=ery_corr; RUN;


%web_open_table(ery_corr);

/*-----------------------*/
/* Multiple regression  */
/*----------------------*/
/*PFC*/
proc mixed data=ery_corr;
	class time;
	where time > 0;
	model pfc_E_W = CCK_E_W time / solution influence residual;
	repeated time / subject=subject type=arh(1) r rcorr;
	run;
quit;

proc mixed data=ery_corr;
	class time;
	where time > 0;
	model pfc_E_W = PYY_E_W time / solution influence residual;
	repeated time / subject=subject type=arh(1) r rcorr;
	run;
quit;

proc mixed data=ery_corr;
	class time;
	where time > 0;
	model pfc_E_W = GLP_E_W time / solution influence residual;
	repeated time / subject=subject type=arh(1) r rcorr;
	run;
quit;
/**/

/*Fullness*/
proc mixed data=ery_corr;
	class time;
	where time > 0;
	model fullness_E_W = CCK_E_W time / solution influence residual;
	repeated time / subject=subject type=arh(1) r rcorr;
	run;
quit;

proc mixed data=ery_corr;
	class time;
	where time > 0;
	model fullness_E_W = PYY_E_W time / solution influence residual;
	repeated time / subject=subject type=arh(1) r rcorr;
	run;
quit;

proc mixed data=ery_corr;
	class time;
	where time > 0;
	model fullness_E_W = GLP_E_W time / solution influence residual;
	repeated time / subject=subject type=arh(1) r rcorr;
	run;
quit;

/*Gastric emptying*/
proc mixed data=ery_corr;
	class time;
	where time > 0;
	model GE_E_W = PYY_E_W time / solution influence residual;
	repeated time / subject=subject type=arh(1) r rcorr;
	run;
quit;

proc mixed data=ery_corr;
	class time;
	where time > 0;
	model GE_E_W = CCK_E_W time / solution influence residual;
	repeated time / subject=subject type=arh(1) r rcorr;
	run;
quit;

proc mixed data=ery_corr;
	class time;
	where time > 0;
	model GE_E_W = GLP_E_W time /solution influence residual;
	repeated time / subject=subject type=arh(1) r rcorr;
	run;
quit;

proc mixed data=ery_corr;
	class time;
	where time > 0;
	model GE_E_W = pfc_E_W time /solution influence residual;
	repeated time / subject=subject type=arh(1) r rcorr;
	run;
quit;

proc mixed data=ery_corr;
	class time;
	where time > 0;
	model GE_E_W = fullness_E_W time /solution influence residual;
	repeated time / subject=subject type=arh(1) r rcorr;
	run;
quit;


