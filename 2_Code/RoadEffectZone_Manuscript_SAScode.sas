* Standard starting stuff;
dm "output;clear;log;clear";
ods listing;ods html close;ods graphics off;
options pageno = 1 linesize =166;

options nocenter;
proc datasets lib = work kill memtype = data;
run;

*-------------------------------------------------------------------------------------------------------;
* Load the data;
proc import out = data
datafile = "L:\TAMUK\Data_and_Analyses\Chapter_REZ\Analysis_SAS\Data\PERMANOVA_data_REZ_20231206.xlsx"
dbms = xlsx replace; sheet = "all_PD_rm"; getnames = yes;
run;

proc print data = data; 
run; quit;

*-------------------------------------------------------------------------------------------------------;
* Run the mixed effect model on each species;
** Need to brute force these because I can't figure out how to add the appropriate species label to each lsmeans output in a macro;
*** Armadillo;
proc mixed data = data; 
	class Transect Location Distance month; 
	model armadillo = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = toep r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm_out;
run;quit;

data lsm_out;
	set lsm_out;
	Length species $ 100;
	species = 'armadillo';
run; 

***Badger; 
proc mixed data = data; 
	class Transect Location Distance month; 
	model badger = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = toep r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm;
run;quit;

data lsm;
	set lsm;
	species = 'badger';
run; 


data lsm_out;
	set lsm_out lsm;
run;

***Bobcat;
proc mixed data = data; 
	class Transect Location Distance month; 
	model bobcat = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = ar(1) r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm;
run;quit;

data lsm;
	set lsm;
	species = 'bobcat';

data lsm_out;
	set lsm_out lsm;
run;

***Cottontail;
proc mixed data = data; 
	class Transect Location Distance month; 
	model cottontail = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = ar(1) r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm;
run;quit;

data lsm;
	set lsm;
	species = 'cottontail';

data lsm_out;
	set lsm_out lsm;
run;

***Coyote;
proc mixed data = data; 
	class Transect Location Distance month; 
	model coyote = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = ar(1) r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm;
run;quit;

data lsm;
	set lsm;
	species = 'coyote';

data lsm_out;
	set lsm_out lsm;
run;

***Feral Hog;
proc mixed data = data; 
	class Transect Location Distance month; 
	model feral_hog = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = ar(1) r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm;
run;quit;

data lsm;
	set lsm;
	species = 'feral_hog';

data lsm_out;
	set lsm_out lsm;
run;

***Javelina;
proc mixed data = data; 
	class Transect Location Distance month; 
	model javelina = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = ar(1) r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm;
run;quit;

data lsm;
	set lsm;
	species = 'javelina';

data lsm_out;
	set lsm_out lsm;
run;

***Nilgai;
proc mixed data = data; 
	class Transect Location Distance month; 
	model nilgai = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = ar(1) r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm;
run;quit;

data lsm;
	set lsm;
	species = 'nilgai';

data lsm_out;
	set lsm_out lsm;
run;

***Opossum;
proc mixed data = data; 
	class Transect Location Distance month; 
	model opossum = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = ar(1) r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm;
run;quit;

data lsm;
	set lsm;
	species = 'opossum';

data lsm_out;
	set lsm_out lsm;
run;

***Raccoon;
****This has to be run using a different correlation structure than the others because of infinite likelihood issues;
****ar(1), vc do not work. cs, toep, arma(1,1) work;
proc mixed data = data; 
	class Transect Location Distance month; 
	model raccoon = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = toep r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm;
run;quit;

data lsm;
	set lsm;
	species = 'raccoon';

data lsm_out;
	set lsm_out lsm;
run;

***Striped skunk;
proc mixed data = data; 
	class Transect Location Distance month; 
	model striped_skunk = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = ar(1) r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm;
run;quit;

data lsm;
	set lsm;
	species = 'striped_skunk';

data lsm_out;
	set lsm_out lsm;
run;

***White-tailed deer;
proc mixed data = data; 
	class Transect Location Distance month; 
	model white_tailed_deer = Location Distance Location*Distance month Location*month Distance*month Location*Distance*month/outpm = out ddfm = satterth;
	random Transect(Location);
	random Distance*Transect(Location);
	repeated month/subject = Transect*Distance*Location type = ar(1) r = 1, 2; 
	lsmeans Location|Distance|month;
	ods output lsmeans = lsm;
run;quit;

data lsm;
	set lsm;
	species = 'white-tailed_deer';

data lsm_out;
	set lsm_out lsm;
run;

proc print data = lsm_out;
run;quit;

*-------------------------------------------------------------------------------------------------------;
*Export the lsmeans data for use in R;
proc export data = lsm_out
outfile = "L:\TAMUK\Data_and_Analyses\Chapter_REZ\Analysis_SAS\LSM_REZ_20240104.csv" dbms = csv;
run;quit;
