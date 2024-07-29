************************
  MACROS FOR BOOTSTRAP
************************;
** separate samples from the larger bootstrap sample;

%macro separate(start, end);
	%do i=&start %to &end;

		data bootdata.bootsample_r_&i;
			set bootdata.bootsample_r(where=(sampleid=&i));
		run;

	%end;
%mend;

************************
  MACROS FOR JACKKNIFE
************************;

%macro jackdata(start, end);
	/* use macro for %do processing utility; enable parallel processing */
%do i=&start %to &end;

		/* do loop to create all samples */;

		data jackdata_&i;
			set jackdata.origjack;

			if obsnum ne &i;
			repeat=&i;
			keep hhidpn repeat;
		run;

		data jackdata.jackdata_r_&i;
			merge mre.leaaf_cms_dem_long_primary(where=(final_drop=0 and usetime=1)) 
				jackdata_&i(in=in1);
			by hhidpn;

			if in1;
		run;

	%end;
%mend;

********************
  MACROS FOR LEAAF
********************;
** Run all steps that lead to LEAAF (design, order, model) then LEAAF;
** NOTE: as currently written, the macro will only export the LEAAF results to a permanent dataset;

%macro allsteps(start, end, libdata, lib);
	%do i=&start %to &end;

		/* Create dataset to order variables */
		data order_prep_&i;
			retain ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA 
				female coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH 
				time1 - time18 event time;
			set &libdata._r_&i;
			keep ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA 
				female coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH 
				time1 - time18 event time;
		run;

		/* Create design datasets with risk factors and covariates. Retain statement is used to order the variables: Risk factors first, then covariates */
		data design_prep;
			retain ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA 
				female coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH;
			set &libdata._r_&i;
			keep ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA 
				female coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH;
		run;

		/* Create final design datasets with all observed combinations of risk factors and covariates using the noduplicates option on proc sort */
		proc sort data=design_prep noduplicates out=design_prep_&i;
			by _all_;
		run;

		/* Create vector of coefficients from modeling */
		ods output ParameterEstimates=point_Parms1;

		proc logistic data=&libdata._r_&i;
			class ami(ref='0') ATRIAL_FIB(ref='0') diabetes(ref='0') chf(ref='0') 
				hypert(ref='0') ischemicheart(ref='0') RA_OA(ref='0') STROKE_TIA(ref='0') 
				female(ref='0') coupled_cf(ref='0') b_first_age(ref='0') 
				b_networth(ref='0') b_bmi(ref='0') b_edu(ref='0') b_adl(ref='0') 
				raceB (ref='0') ethnH (ref='0') / param=ref;
			model event(event='1')=time ami ATRIAL_FIB diabetes chf hypert ischemicheart 
				RA_OA STROKE_TIA female b_edu b_networth b_bmi coupled_cf b_first_age b_adl 
				raceB ethnH;
		run;

		data outbeta_&i;
			set point_Parms1;
			keep estimate;
		run;

		/* LEAAF analyses */
		proc iml;
			/* Read in the ordering dataset */
			varnames={ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA 
				female coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH 
				time1 time2 time3 time4 time5 time6 time7 time8 time9 time10 time11 time12 
				time13 time14 time15 time16 time17 time18 event time};
			use order_prep_&i;
			read all var _num_ into originaldata[colname=varnames];
			close order_prep_&i;

			/* Read the design dataset into matrix desmtx*/
			varnames2={ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA 
				female coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH};
			use design_prep_&i;
			read all var _num_ into desmtx[colname=varnames2];
			close design_prep_&i;

			/* Enter the number of risk factors you have read into matrices.
			Risk factors only not covariates*/
			nrisk=8;

			/* Enter the number of covariates you have read into matrices.*/
			ncov=9;

			/* Enter the maximum number of time points.*/
			maxyear=18;
			nx=nrisk+ncov;
			ninter=0;
			nmandc=nrisk;
			x=originaldata[, (1:(nx+maxyear+1))];
			*x=originaldata[, (1:(nx+2))] ;
			xsub=originaldata[, (1:(nx))];
			year=originaldata[, ((nx+1):(nx+maxyear))];
			*year=originaldata[, (1:(nx+2))];
			outcomea=originaldata[, (nx+maxyear+1)];
			N=nrow(x);
			package load listutil;
			lendesmtx=nrow(desmtx);

			/* Select the columns for the risk factors from the design matrix */
			sumnonzero1=desmtx[, 1:Nrisk];
			* Sum the cols of the design matrix risk factors only;
			sumnonzeros=sumnonzero1[, +];
			* Get maximum # of risk factors;
			maxnonzeros=max(sumnonzeros);

			/* Function that creates factorial matrix of 2 ## n risk factors*/
			start FullFactorial(n);
			r=2 ## n;
			xf=j(r, n);
			pot=2 ## ((n-1):0);

			do i=1 to r;
				xf[i, ]=band(i-1, pot) > 0;
			end;
			return(xf);
			finish;
			nz2enum=listcreate(maxnonzeros);
			subsets=listcreate(maxnonzeros);
			loczero=listcreate(maxnonzeros);

			do i=1 to maxnonzeros;
				valf=0:(2##i-1);
				call ListSetItem (nz2enum, i, valf);
				call ListSetItem (subsets, i, {0});
				call ListSetItem (loczero, i, {0});
				a=FullFactorial(i);
				call ListSetItem (subsets, i, a);
				call ListSetItem (loczero, i, loc(a));

				/*Locations of ones*/
			end;

			/* Read SAS dataset of vector of coefficients in.
			SAS Order should be intercept, time, risk factors, covariates.
			Should only be one variable in dataset */
			varnames3={betacoeff};
			use outbeta_&i;
			read all var _num_ into outbeta_condIntercova[colname=varnames3];
			close outbeta_&i;
			betayeara=outbeta_condIntercova[1];
			betayrona=outbeta_condIntercova[2];
			rowbetaa=nrow(outbeta_condIntercova);
			betaa=outbeta_condIntercova[(3):rowbetaa];
			betaconda=betaa[1:nrisk];
			betaintera=betaa[(nrisk+1):(nrisk+ninter)];

			if NCov ^=0 then
				do;
					betacova=outbeta_condIntercova[(2+nrisk+ninter+1):rowbetaa];
				end;
			AAFbyYearMtxa=j(maxyear, nmandc, 0);
			NskByYearMtx=j(maxyear, lendesmtx, 0);
			N1a=j(maxyear, 1, 0);
			Nyr=j(maxyear, 1, 0);

			do yr=1 to maxyear;
				N1a[yr]=sum(outcomea[loc(year[, yr]=1), ] );
				beta0a=betayeara[1]+(betayrona[1]*yr);
				xmat=xsub[loc(year[, yr]=1), ];
				Nyr[yr]=nrow(xmat);
				xPs0va=j(NMandC, 1, 0);
				sumbetacova=j(1, 1, 0);

				/* Below looking for row of xmat that matches row i of
				desgnmatrix. If sum of logical elements are equal then it is a match */
				do r=1 to lendesmtx;
					xmat2=repeat(desmtx[r, ], Nyr[yr], 1);
					xmat3=(xmat[, ]=xmat2[, ]);
					temv=xmat3[, +];
					temv2=(temv[, ]=nx);
					NskByYearMtx[yr, r]=temv2[+, ];
					row=DesMtx[r, ];
					rowMandC=row[, 1:NMandC];
					numNonZeros=row[, +];
					* C Loop;

					if Ncov > 0 then
						do;
							rowcov=row[, (1+NMandC):NX];
							sumBetaCova=rowcov*betacova;
						end;

					/* C Loop */
					rowcond=row[, 1:NMandC];

					If Ncov=0 then
						sumbetaCova=0;

					/* First get indicators for columns */
					/* Find indices for nonzero elements in design matrix row, for meds and conditions */
					nzsloc=loc(rowMandC);
					numNzs=ncol(nzsloc);

					/* Count how many not zero */
					* D Loop;

					if numNzs > 0 then
						do;
							subs=FullFactorial(numNzs);
							numSubs=nrow(subs);
							ps0vectora=j(numsubs, 1, 0);
							total=j(numsubs, 1, 0);
							* E Loop;

							do j=1 to numsubs;
								rowsubMandC=j(1, NMandC, 0);
								mat_row=j(j, numNzs, 0);
								* F Loop;

								do k=1 to numnzs;
									mat_row[j, k]=subs[j, k]*nzsloc[1, k];
									scale=mat_row[j, k];

									if scale ^=0 then
										rowsubMandC[1, scale]=1;
								end;

								/* F Loop */
								rowsubcond=rowsubMandC[1, 1:NMandC];
								condsuma=rowsubcond*betaconda;
								sumcoeffa=((-1*beta0a) - condsuma -sumBetaCova);
								denoma=1+exp(sumcoeffa);
								ps0vectora[j]=1 / denoma;
							end;

							/* E loop */
							/* DO Loop G */
							*********************************************;
							*Need to generate a nz2enum list with 3 rows.;
							*These are cells with 0,1  0,1,2,3 0,1,2,3,4,5,6,7;
							************************************************;
							enum=(listgetitem(nz2enum, numnzs))`;

							do g=1 to numnzs;
								xk=nzsloc[g];
								p=2##(numnzs-g);
								x=(band(enum, p)=p);
								y1a=sum(Ps0vectora[loc(x[, 1]=1), ]);
								y2a=sum(Ps0vectora[loc(x[, 1]=0), ]);
								combinea=(y1a-y2a)/(numsubs/2);
								combinea=combinea/n1a[yr];
								combinea=NskByYearMtx[yr, r]*combinea;
								xPs0va[xk]=xPs0va[xk]+combinea;
							end;

							/* G Loop */
						end;

					/* D Loop */
				end;

				/* A Loop */
				AAFbyYearMtxa[yr, ]=xPs0va`;
				AAFbyYearMtx2a=AAFbyYearMtxa#Nyr;
			end;
			allaveragea=(AAFbyYearMtx2a[+, ])/N;
			sumLEAAFa=sum(allaveragea);
			create outa from allaveragea;

			/* create data set and variables */
			append from allaveragea;

			/* write all rows */
			close;

			/* Specify SAS dataset name to output LE-AAFs */
			create &lib..leaaf_r_&i var {allaveragea sumLEAAFa};
			append;
			close &lib..leaaf_r_&i;
		run;

		quit;
	%end;

	proc datasets nolist noprint;
		delete order_prep_&i design_prep design_prep_&i 
		point_Parms1 outbeta_&i outa;
	quit;

%mend;

* Combine LE_AAF estimates from bootstrap or jackknife steps ;

%macro combineLEAAF(end, libdata, outlib, name);
	proc datasets nolist noprint library=work;
		delete &name._leaaf;
	quit;

	%do i=1 %to &end;

		data create_&i;
			set &libdata._&i;
			LE_AAFA=ALLAVERAGEA*100;
			LE_AAFA=round(LE_AAFA, .01);
		run;

		proc transpose data=create_&i out=out_&i;
		run;

		data out_&i;
			set out_&i;
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
			sample=&i;
			drop i aol1-aol8;
			label col1='ami' col2='ATRIAL_FIB' col3='diabetes' col4='chf' col5='hypert' 
				col6='ischemicheart' col7='RA_OA' col8='STROKE_TIA';
		run;

		proc append base=&name._leaaf data=out_&i(where=(_name_='LE_AAFA'));
		run;

	%end;

	data &outlib..&name._leaaf_COMBINED;
		set &name._leaaf;
	run;

	proc datasets library=work nolist noprint;
		delete create_: out_:;
	quit;

%mend;

/* Produce BCA Bounds */
* reference: https://www.lexjansen.com/phuse/2005/pk/pk02.pdf;

%macro bca_debug(libboot, libjack, cond, labcond, numboot);
	proc sort data=&libboot..boot_leaaf_COMBINED;
		/* sort estimates for LE_AAF from low to high */
		by col&cond;

		/* each column has LEAAF Estimates for a condition*/
	run;

	/* Bootstrap Interval by Percentile */
	data boot_leaaf_combined_complete;
		set &libboot..boot_leaaf_COMBINED end=eof;

		if col&cond=. then
			delete;
	run;

	data ci_perc;
		set boot_leaaf_combined_complete end=eof;
		retain conf_lo conf_hi;

		/* Pull out .025 x n and .975x n values as percentile confidence limits. */
		if _n_=floor(0.025 * &numboot) then
			conf_lo=col&cond;

		if _n_=ceil(0.975 * &numboot) then
			conf_hi=col&cond;

		if eof then
			output;
		keep conf_lo conf_hi;
	run;

	proc print;
		title "Percentile CI for col&cond &labcond";
	run;

	/* set original sample single leaaf estimate dataset and
	bootstrap sample LEAAF estimates dataset together to calculate means*/
	data bootorig(drop=_name_);
		set tor.original_leaaf_trans_r(in=a) boot_leaaf_combined_complete end=eof;

		if a then
			sample=0;
	run;

	/* Calculate Means for Original and Bootstrap LEAAF estimates */
	proc means data=bootorig noprint nway;
		class sample;
		var col&cond;
		output out=bootmean /* mean for each sample */
		mean=mean_col&cond;
	run;

	/* Estimate Bias (from bootstrap samples) */
	data bootdist(keep=sample mean_col&cond) bias_col&cond(keep=bias propless);
		set bootmean end=eof;
		retain orig_leaaf;

		if sample=0 then
			orig_leaaf=mean_col&cond;

		if mean_col&cond lt orig_leaaf then
			lessthan=1;

		/* flag if bootstrap sample gives lower estimate */
		else
			lessthan=0;
		retain nless 0;

		if sample gt 0 then
			nless=nless+lessthan;

		/* count samples with flag for less than comparison with original estimate */
		if sample ne 0 then
			output bootdist;

		/* bootstrap le_aaf */
		if eof then
			do;

				/* for the last value calculate: */
				propless=nless/sample;

				/* 1. proportion of values below original estimate */
				bias=probit(propless);

				/* 2. inverse normal of that proportion */
				output bias_col&cond;

				/* 3. output only that record to new data set */
			end;
	run;

	proc print data=bias_col&cond;
	run;

	/* Estimate Acceleration (from jackknife samples) */
	/* Access jackknife samples */
	data jack_leaaf_combined_complete;
		set &libjack..jack_leaaf_COMBINED;

		if col&cond=. then
			delete;
	run;

	proc means data=jack_leaaf_combined_complete noprint nway;
		class sample;
		var col&cond;
		output out=jacksum 						/* mean for each sample */
		mean=jmean_col&cond;
	run;

	proc sql noprint;
		select mean(jmean_col&cond) 
			/* put mean of jackknifed values into macro variable */
			into :meanjack1 from jacksum;
	quit;

	data jacksum2;
		set jacksum;
		cubed=(&meanjack1 - jmean_col&cond)**3;
		/* create cubed value of difference */
		squared=(&meanjack1 - jmean_col&cond)**2;
		/* create squared value of difference */
	run;

	proc means data=jacksum2 noprint;
		output out=jacksum3 
		sum(cubed)=sumcube /* find sum of cubed values */
		sum(squared)=sumsquar; /* find sum of squared values */
	run;

	data accel_col&cond;
		set jacksum3;
		accel=sumcube / (6 * (sumsquar**1.5));

		/* plug values into equation for */
		/* the acceleration statistic */
		keep accel;
	run;

	proc print;
		var accel;
	run;

	/* obtain BCA confidence intervals */
	data ciends;
		merge accel_col&cond 
	   bias_col&cond;
		part1=(bias + probit(0.025)) / (1 - (accel*(bias + probit(0.025))));
		part2=(bias + probit(0.975)) / (1 - (accel*(bias + probit(0.975))));
		alpha1=probnorm(bias + part1);
		alpha2=probnorm(bias + part2);
		n1=alpha1*&numboot;
		n2=alpha2*&numboot;
		call symput('n1', put(floor(n1), 5.));

		/* Create macro variables with values */
		call symput('n2', put(ceil(n2), 5.));

		/* of N1 and N2 for later use */
	run;

	proc sort data=boot_leaaf_combined_complete;
		by col&cond;
	run;

	data ci_bca;
		set boot_leaaf_combined_complete end=eof;
		retain bca_conf_lo bca_conf_hi;

		if _n_=&n1 then
			bca_conf_lo=col&cond;

		/* select values for upper and lower */
		if _n_=&n2 then
			bca_conf_hi=col&cond;

		/* limits using N1 and N2 values */
		if eof then
			output;
		keep bca_conf_lo bca_conf_hi;
	run;

	proc print;
		var bca_conf_lo bca_conf_hi;
		title 'Before Merge with Original Estimate';
	run;

	data orig_ci_col&cond;
		merge ci_bca tor.original_leaaf_trans_r;
		keep col&cond bca_conf:;
	run;

	proc print;
		var col&cond bca_conf_lo bca_conf_hi;
		title "col&cond &labcond  Original Estimate and BCA interval";
	run;

	/* reformat output dataset for concatenation with other conditions */
	data orig_ci_col&cond;
		set orig_ci_col&cond;
		rename col&cond=leeaf;
		format desc $20.;
		order=&cond;
		desc="&&labcond";
	run;

%mend;