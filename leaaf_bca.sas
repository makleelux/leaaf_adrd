*******************************************************************************************************

Compute LE-AAF BCA Intervals

*******************************************************************************************************

* parallelize process with multiple SAS sessions via MP connect 
(see https://blogs.sas.com/content/sgf/2021/01/13/running-sas-programs-in-parallel-using-sas-connect/);
options threads cpucount=actual sascmd="sas" mlogic mprint mlogic 
	symbolgen /* macro options */;

proc options option=cpucount;
run;

* SAS library where analytic datset and data from previous steps stored;
libname mre /* enter path */;
* SAS libraries where bootstrapped and jackknifed samples are stored;
libname bootdata /* enter path */;
libname jackdata /* enter path */;
* Assign number of bootstrap samples;
%let NumSamples=3000;
* load SAS macros;
%include /* enter path */;
* create library for saving bootstrap and jackknife analyses;
libname bootdoc /* enter path */;
libname jackdoc /* enter path */;

**** CREATE BOOTSTRAP SAMPLES ****;
data final_LEAAF;
	set mre.tidy_long(where=(final_drop=0 and usetime=1));
	study_id=hhidpn;
run;

/* Generate bootstrap samples */
proc surveyselect data=final_Leaaf seed=678 
		out=bootdata.bootsample(rename=(Replicate=SampleID)) samprate=1 
		method=urs /* resample with replacement */
		OUTHITS /* option to suppress the frequency var */
		reps=&NumSamples; /* generate NumSamples bootstrap resamples*/
	samplingunit study_id; /*Identify sampling units as groups of observations. This way all time intervals are kept per id.*/
run;

/* parallelize separatation of bootstraps */
/** process 1 **/
signon task1;
rsubmit task1 wait=no inheritlib=(bootdata);
%include /* path to leaaf_macro.sas */;
%separate(start=1, end=100);
endrsubmit;

/** process ... **/
waitfor _all_;
signoff _all_;

**** BOOTSTRAP LEAAF ANALYSES ****;
/** process 1 **/
signon task1;
rsubmit task1 wait=no inheritlib=(tor bootdata bootdoc);
%include /* path to leaaf_macro.sas */;

proc printto 
		log=/* enter path */'\bootdoc\log-file_leaafout1';
run;

proc printto 
		print=/* enter path */'bootdoc\output-file_leaafout1';
run;

%allsteps(start=1, end=100, libdata=bootdata.bootsample, lib=bootdoc);
endrsubmit;

/** process ... **/
waitfor _all_;
signoff _all_;

**** JACKNIFE ****;
**** CREATE JACKNIFE SAMPLES ****;
data first;
	set mre.tidy_long(where=(final_drop=0 and usetime=1));
	by hhidpn;

	if first.hhidpn;
run;

data original;
	set first;
run;

%global nobs;
* save this dataset to jackdoc library so it will be accessible by parallel processing steps;

data jackdata.origjack;
	/* create a new data set which contains observation */
	set original end=eof;

	/* numbers 1 to &nobs (no. obs in data set) */
	obsnum=_n_;

	if eof then
		call symput('nobs', put(obsnum, 2.));
run;

/* parallelize jackknife sampling */
/** process 1 **/
signon task1;
rsubmit task1 wait=no inheritlib=(mre tor jackdata);
%include /* path to leaaf_macro.sas */;

proc printto log=/* enter path */'\jackdoc\log-file1';
run;

proc printto print=/* enter path */'\jackdoc\output-file1';
run;

%jackdata(start=1, end=100);
endrsubmit;

/** process ... **/
waitfor _all_;
signoff _all_;

**** JACKNIFE LEAAF ANALYSES ****;
/** process 1 **/
signon task1;
rsubmit task1 wait=no inheritlib=(tor jackdata jackdoc);
%include /* path to leaaf_macro.sas */;

proc printto 
		log=/* enter path */'\jackdoc\log-file_leaafout1';
run;

proc printto 
		print=/* enter path */'\jackdoc\output-file_leaafout1';
run;

%allsteps(start=1, end=100, libdata=jackdata.jackdata, lib=jackdoc);
endrsubmit;

/** process ... **/
waitfor _all_;
signoff _all_;

***** COMBINE BOOTSTRAP AND JACKNIFE RESULTS FOR BCA *****;
** Calculate bias and acceleration factors **;
* access the original point estimate of the LE_AAF;

data original;
	set mre.Leaaf1_primaryanalysis(keep=r_comp);
	rename r_comp=LE_AAFA;
run;

* transpose this dataset;
proc transpose data=original out=out_original;
run;

data mre.original_leaaf_trans;
	set out_original(where=(_name_='LE_AAFA'));
	col9=sum(of col1-col8);
	col11=sumabs(of col1-col8);
	array cond1[8] col1-col8;
	array cond2[8] aol1-aol8;

	do i=1 to 8;
		cond2[i]=cond1[i];

		if cond2[i] < 0 then
			cond2[i]=0;
	end;
	col10=sum(of aol1-aol8);
	drop _name_ i aol1-aol8;
run;

proc print;
run;


* combine bootstrap estimates from LEAAF;
%combineLEAAF(end=3000, libdata=bootdoc.leaaf, outlib=bootdoc, name=boot);

* combine jackknife estimates from LEAAF;
%combineLEAAF(end=9944, libdata=jackdoc.leaaf, outlib=jackdoc, name=jack);

%bca_debug(libboot=bootdoc, libjack=jackdoc, cond=1, labcond=ami, numboot=3000);
%bca_debug(libboot=bootdoc, libjack=jackdoc, cond=2, labcond=ATRIAL_FIB, numboot=3000);
%bca_debug(libboot=bootdoc, libjack=jackdoc, cond=3, labcond=diabetes, numboot=3000);
%bca_debug(libboot=bootdoc, libjack=jackdoc, cond=4, labcond=chf, numboot=3000);
%bca_debug(libboot=bootdoc, libjack=jackdoc, cond=5, labcond=hypert, numboot=3000);
%bca_debug(libboot=bootdoc, libjack=jackdoc, cond=6, labcond=ischemicheart, numboot=3000);
%bca_debug(libboot=bootdoc, libjack=jackdoc, cond=7, labcond=RA_OA, numboot=3000);
%bca_debug(libboot=bootdoc, libjack=jackdoc, cond=8, labcond=STROKE_TIA, numboot=3000);
%bca_debug(libboot=bootdoc, libjack=jackdoc, cond=9, labcond=Total, numboot=3000);
%bca_debug(libboot=bootdoc, libjack=jackdoc, cond=10, labcond=Total_Pos, numboot=3000);
%bca_debug(libboot=bootdoc, libjack=jackdoc, cond=11, labcond=Total_AbsPos, numboot=3000);

* concatenate leaaf and bca interval results;
data mre.final_leaaf_bca;
	set orig_ci_col1-orig_ci_col11;
run;