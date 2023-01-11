/*-------------------------------*/
/* IMPORT DATA FROM EXCEL FILE   */
/*-------------------------------*/

%web_drop_table(ery_git);


FILENAME REFFILE '/home/u50127452/sasuser.v94/PolyAlluLac/SAS_PolyAlluLac_05_10_2021.xlsx';

PROC IMPORT DATAFILE=REFFILE
	DBMS=XLSX
	OUT=ery_git;
	GETNAMES=YES;
	SHEET="GIT_symptoms";
RUN;

PROC CONTENTS DATA=ery_git; RUN;


%web_open_table(ery_git);

/*----------------------------------------*/
/* RECODE CONDITIONS FOR FACTORIAL DESIGN */
/*----------------------------------------*/

data ery_git;
set ery_git;
attrib substance format=$10.; /* initializes substance as a categorical variable with length 10 characters */
if condition = "25g D-allulose + 450ppm lactisole" then substance = "allulose";
if condition = "25g D-alluose" then substance = "allulose";
if condition = "300mL tap water + 450ppm lactisole" then substance = "water";
if condition = "300mL tap water" then substance = "water";
if condition = "50g erythritol + 450ppm lactisole" then substance = "erythritol"; 
if condition = "50g erythritol" then substance = "erythritol";
run; 

data ery_git;
set ery_git;
attrib lactisole format=$3.;
if condition = "25g D-allulose + 450ppm lactisole" then lactisole = "yes";
if condition = "25g D-alluose" then lactisole = "no";
if condition = "300mL tap water + 450ppm lactisole" then lactisole = "yes";
if condition = "300mL tap water" then lactisole = "no";
if condition = "50g erythritol + 450ppm lactisole" then lactisole = "yes"; 
if condition = "50g erythritol" then lactisole = "no";
run; 


/*---------------------------------------------------*/
/* CHECK DISTRIBUTION OF VARIABLE DIARRHEA           */
/*---------------------------------------------------*/

/*proc univariate data=ery_git;
histogram diarrhea;
run;

proc univariate data=ery_git;
var diarrhea;
histogram diarrhea;
run;*/

proc freq data=ery_git;
table diarrhea;
table diarrhea*condition*time;
run;

proc freq data=ery_git;
table abdominal_pain;
table abdominal_pain*condition*time;
run;

proc freq data=ery_git;
table nausea;
table nausea*condition*time;
run;

proc freq data=ery_git;
table vomiting;
table vomiting*condition*time;
run;

proc freq data=ery_git;
table bowel_sounds;
table bowel_sounds*condition*time;
run;

proc freq data=ery_git;
table bloating;
table bloating*condition*time;
run;

proc freq data=ery_git;
table eructation;
table eructation*condition*time;
run;

proc freq data=ery_git;
table flatulence;
table flatulence*condition*time;
run;
/*abdominal pain, nausea, vomiting, bowel sound, bloating, eructation, flatulence*/


