*******************************************************************************************************

Compute LE-AAF point estimates

*******************************************************************************************************

*** LE-AAF ANALYSES AND RELATED STEPS ***;
** These analyses include the following steps:
   1. order the variables for LE-AAF analyses
   2. create the design matrix for LE-AAF analyses
   3. perform pooled logistic regression to estimate coefficients for LE-AAF analyses
   4. conduct LE-AAF analyses with the same risk factors and covariates
   
   IML code for LE-AAF analyses are from the GRASP website: https://www.peppercenter.org/public/grasp.cfm

** set SAS library;
libname mre /* enter path */;

************
** STEP 1 **
************;
** order variables by type in model;

data leaaf0;
	retain ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA 
		female coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH 
		time1 - time18 event time;
	set mre.tidy_long(where=(final_drop=0 and usetime=1));
	keep ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA female 
		coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH 
		time1 - time18 event time;
run;

** check order of variables;

proc contents data=leaaf0 position;
run;

************
** STEP 2 **
************;
** create design matrix;

data design0_prep;
	retain ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA 
		female coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH;
	set mre.tidy_long(where=(final_drop=0 and usetime=1));
	keep ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA female 
		coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH;
run;

proc sort data=design0_prep noduplicates out=design0;
	by _all_;
run;

* check position of variables;

proc contents data=design0 position;
run;

************
** STEP 3 **
************;
** run pooled logistic regression & create column of coefficients;

ods output ParameterEstimates=beta0;

proc logistic data=mre.tidy_long(where=(final_drop=0 and usetime=1));
	class ami(ref='0') ATRIAL_FIB(ref='0') diabetes(ref='0') chf(ref='0') 
		hypert(ref='0') ischemicheart(ref='0') RA_OA(ref='0') STROKE_TIA(ref='0') 
		female(ref='0') coupled_cf(ref='0') b_first_age(ref='0') b_networth(ref='0') 
		b_bmi(ref='0') b_edu(ref='0') b_adl(ref='0') raceB (ref='0') ethnH (ref='0') 
		/ param=ref;
	model event(event='1')=time ami ATRIAL_FIB diabetes chf hypert ischemicheart 
		RA_OA STROKE_TIA female b_edu b_networth b_bmi coupled_cf b_first_age b_adl 
		raceB ethnH;
run;

** output results;

data beta0;
	set beta0;
	keep estimate;
run;

************
** STEP 4 **
************;
** IML code for LE-AAF analyses;

proc iml;
	varnames={ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA 
		female coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH time1 
		time2 time3 time4 time5 time6 time7 time8 time9 time10 time11 time12 time13 
		time14 time15 time16 time17 time18 event time};
	use leaaf0;
	read all var _num_ into originaldata[colname=varnames];
	close leaaf0;

	/* Read the design dataset into matrix desmtx */
	varnames2={ami ATRIAL_FIB diabetes chf hypert ischemicheart RA_OA STROKE_TIA 
		female coupled_cf b_first_age b_networth b_bmi b_edu b_adl raceB ethnH};
	use design0;
	read all var _num_ into desmtx[colname=varnames2];
	close design0;

	/* Enter the number of risk factors you have read into matrices.
	Risk factors only not covariates */
	nrisk=8;

	/* Enter the number of covariates you have read into matrices. */
	ncov=9;

	/* Enter the maximum number of time points. */
	maxyear=18;
	nx=nrisk+ncov;
	ninter=0;
	nmandc=nrisk;
	x=originaldata[, (1:(nx+maxyear+1))];
	xsub=originaldata[, (1:(nx))];
	year=originaldata[, ((nx+1):(nx+maxyear))];
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

	/* Function that creates factorial matrix of 2 ## n risk factors */
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

		/* Locations of ones */
	end;

	/* Read SAS dataset of vector of logistic regression coefficients into
	SAS Order should be intercept, time, risk factors, covariates.
	Should only be one variable in dataset */
	varnames3={betacoeff};
	use beta0;
	read all var _num_ into outbeta_condIntercova[colname=varnames3];
	close beta0;
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
			/* Find indices for nonzero elements in design matrix row */
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
	create out_LEAAF0 var {allaveragea sumLEAAFa};
	append;
	close out_LEAAF0;
run;

quit;

/* Specify SAS dataset in the set statement below */

data out_LEAAF0;
	set out_LEAAF0;
	LE_AAFA=ALLAVERAGEA*100;
	LE_AAFA=round(LE_AAFA, .01);
	
	* create labels for LE-AAF results;
	format label $20.;

	if _n_=1 then
		label="ami";
	else if _n_=2 then
		label="ATRIAL_FIB";
	else if _n_=3 then
		label="diabetes";
	else if _n_=4 then
		label="chf";
	else if _n_=5 then
		label="hypert";
	else if _n_=6 then
		label="ischemicheart";
	else if _n_=7 then
		label="RA_OA";
	else if _n_=8 then
		label="STROKE_TIA";
run;

* output LE-AAF results with limited columns;

proc sql;
	create table out_LEAAF0_ as select label
  		 , ALLAVERAGEA
		 , SUMLEAAFA
		 , LE_AAFA from out_LEAAF0;
quit;

* combine results of initial LEA-AF models to export to R;

proc sort data=out_LEAAF0_;
	by label;
run;

data leaaf1;
	merge out_LEAAF0_(keep=label sumleaafa LE_AAFA rename=(sumleaafa=r_sum 
		LE_AAFA=r_comp));
	by label;
	* add order variable;

	if label=:"ami" then
		order=1;
	else if label=:"ATRIAL_FIB" then
		order=2;
	else if label=:"diabetes" then
		order=3;
	else if label=:"chf" then
		order=4;
	else if label=:"hypert" then
		order=5;
	else if label=:"ischemicheart" then
		order=6;
	else if label=:"RA_OA" then
		order=7;
	else if label=:"STROKE_TIA" then
		order=8;
run;

proc sort data=leaaf1;
	by order;
run;

data mre.leaaf1_PrimaryAnalysis;
	set leaaf1;
run;

