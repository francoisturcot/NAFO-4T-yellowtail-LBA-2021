//age-atructured length goup model
// Jan 22 2021
// 2 M groups <25 cm, 25+cm  by age: m1 1-3; m2 5+; 4 average of m1 and m2
// 5 M blocks 1985-1990, 1991-1996, 1997-2002, 2003-2008, 2009-2014, 2015-2020
// projections at different catch levels
// recrate added from m3j11c_fixFmcmc
DATA_SECTION
  init_int nT; // last year 36
  init_int nTretro; // last year to use for fitting the likelihoods, used for retrospective analyses
  init_int pT; // last projection year 46
  init_int sage; // start and end ages in catch at age
  init_int lage; // lage is a plus group
  vector age(1,8);
  !! age.fill_seqadd(1,1);
  // Number of M groups
  init_int nSplitAge; // 3: 1-3,4,5+
  // Lower ages of M groups
  init_vector splitAge(1,nSplitAge); // 1 4 5
  // First age for avg'ing F1 ages 1-3
  init_int a1F1;
  // Last age for avg'ing F
  init_int a2F1;
  // First age for avg'ing F2 ages 5-8
  init_int a1F2;
  // Last age for avg'ing F
  init_int a2F2;
  // Number of fisheries
  init_int nFisheries;  
  // Number of selectivity time-blocks
  init_ivector nSelBlocks(1,nFisheries);
  // Selectivity time block boundaries
  init_matrix tBlock(1,nFisheries,1,nSelBlocks);  
  // 8.5/12
  init_number fracYearSurvey;
  // Baranov iterations
  init_int baranovIter;

  // PRIORS
  // recruitment
  init_number priorSD_R;   //  0.5
  // natural mortality
  init_vector prior_Minit(1,2);  //ages 1-3, 5-8
  init_vector mpriorSD(1,2);
  // survey catchability
  init_number logq_prior;     // log(0.7)= -0.36  
  init_number logq_priorSD;   // 0.15 
  // inital abundances
  //init_vector log_priorN(sage,lage+1);
  //init_vector sd_priorN(sage,lage+1); 
  // Initial values for Lgrp prop observation std errors and qs 
  init_vector LikeWeight(1,3); // 1=RVjuv, 2=RVadt, 3=CatchLgrps 
  init_number init_tauLgrp; 
  // PHASES
  init_int ph_log_AvgR; // 1
  //init_int ph_log_initN;   // 1
  init_int ph_S50;         // 1
  init_int ph_S95;         // 1 
  init_int ph_q;           // 1 
  init_int ph_sig;         // 1 
  init_int ph_recDevs;     // 2
  init_int ph_initRecDevs; // 2
  init_int ph_M;           // 3  
  init_int ph_gammaR;      // -1
  // DATA
  init_vector Ct(1,nT); // landings in tonnes
  init_matrix Cwtg(1,2,1,nT); // mean weight (kg) per individual in the fishery removals by stage (juv vs adult) 
  init_matrix Iy(1,2,1,nT); // abundance juv and Adult (1000s)  
  init_matrix Iwtg(1,2,1,nT); // mean weight (kg) per individual in the survey catch by stage (juv vs adult) 
  //Observed proportions-at-age
  init_matrix ageObsProp(1,2,1,nT);  // fishery age props, 1=juv 2=1dt, 1 44
  // observed mean weights at age by year
  init_matrix waaSep(sage,lage,1,nT);
  // observed proportions mature by age and year
  init_matrix matAge(sage,lage,1,nT);
  // Catch levels for projections
  init_vector Cpro(1,3);
  // tracking number of function evaluations
  int neval;
  init_int eof;

      //projection variables for output
        vector Mpro(sage,lage);    // average of last 5 yr
        vector Spro1(nT+1,pT);
        vector Spro2(nT+1,pT);
        vector Spro3(nT+1,pT);
        matrix Fpro(1,3,nT+1,pT); //  Ftg 
        matrix nAgpPro1(1,3,nT+1,pT);    
        matrix nAgpPro2(1,3,nT+1,pT);    
        matrix nAgpPro3(1,3,nT+1,pT);
        matrix Npro(nT+1,pT,sage,lage);
        vector PRav(sage,lage);

 LOC_CALCS
   if(eof!=54321)
   {
     cout<<" data entry error 1985 land "<< Ct(1) <<endl;
     cout<<" data entry error 1985 Adt RV "<< Iy(2,1) <<endl;
     cout<<" data entry error 1985 Ct adt prop "<< ageObsProp(2,1) <<endl;
     cout<<" data entry error eof "<< eof <<endl;
     exit(1);
   }
 END_CALCS

   //*******************************************************************/
PARAMETER_SECTION
  objective_function_value objFunc;
  // Average recruitment: initialized at 
  init_number log_AvgR(ph_log_AvgR);
  // Recruitment deviations: initial abundance 1985
  init_bounded_vector init_recDevs(sage,lage,-5.,5.,ph_initRecDevs);
  // Recruitment deviations: 1977-2014
  init_bounded_vector recDevs(2,nT,-5.,5.,ph_recDevs);
  // Initial natural mortality rates - Mgrps 1 and 3, 2 is average of 1 and 3
  init_bounded_vector log_Mjuv(1,6,-3.0,0.69,ph_M);
  init_bounded_vector log_Madt(1,6,-3.0,0.69,ph_M); //log(2) upper bound
  // Selectivity parameters 1=fishery, 2=survey
  init_bounded_matrix log_S50_gi(1,nFisheries,1,nSelBlocks,0.6,2.,ph_S50);
  init_bounded_matrix log_S95_step_gi(1,nFisheries,1,nSelBlocks,0.0,2.0,ph_S95);
  init_bounded_number lnq(-5,1.5,ph_q);
  init_bounded_vector log_sig(1,2,-3.0,1.5,ph_sig);
  // Initial fishing mortality prior to 1985
  init_bounded_number logFinit(-3.0,0.69,1);
  // Recruitment autocorrelation
  init_number logit_gamma_R(ph_gammaR);

  // Parameters: natural scale
  number avgR;
  number gamma_R;  
  number rvq;
  vector rvsig(1,2);
  number Finit;
 
  matrix Mta(1,nT,sage,lage);
  // Selectivity parameters: age-at-50%
  matrix S50_gt(1,nFisheries,1,nT);
  // Selectivity parameters: age-at-95%
  matrix S95_gt(1,nFisheries,1,nT);
  
  // Posterior quantities
  number recPrior;
  number mPrior;
  number qPrior;
  vector indexLikelihood(1,2);
  number validObs;
  number LgrpLikelihood;
  number tauSquareLgrp
  number ss;
  matrix rvpred(1,2,1,nT);
  matrix rvres(1,2,1,nT);  

  // Derived variables
  vector rw_recDevs(2,nT);  
  matrix nLgrp(1,2,1,nT);
  matrix bLgrp(1,2,1,nT);
  vector ssb(1,nT);
  matrix Sta(1,nT,sage,lage);
  matrix nAgrp(1,3,1,nT);
  matrix expNtg(1,2,1,nT);
  matrix Nta(1,nT,sage,lage);
  matrix Zta(1,nT,sage,lage); 
  vector Ft(1,nT);
  vector recRate(2,nT);  

  sdreport_vector avgF(1,nT);
  vector avgF2(1,nT);
  
  matrix  Bprime(1,2,sage,lage);
  matrix  uCgat(1,2,1,nT);  // fishery only, juv or adt, year
  
  // Selectivity by age
  //matrix sel(1,nFisheries,1,plusGroupAge);
  3darray sel_gta(1,nFisheries,1,nT,sage,lage);
  
  // REPORT_SECTION objects
  // Subgroup exploitation rates, ages juv and adult
  matrix subHt(1,2,1,nT);

GLOBALS_SECTION
    void solveBaranov( const int& t, const int& iter);
    
    // Flag to control whether header written to mcmc output file.
    int mcmcHeader = 0;
    
    // Flag to control which parameters are done in MCMC phase.
    int mcmcFlag = 1;
    
    // Uncomment this line for 32-bit compiler.
    // #include <fstream.h>
    #include <admodel.h>
    ofstream mcoutM1("mcoutM1.dat");
    ofstream mcoutM2("mcoutM2.dat");
    ofstream mcoutM3("mcoutM3.dat");
    ofstream mcoutnAgrp3("mcoutnAgrp3.dat");
    ofstream mcoutnAgrp2("mcoutnAgrp2.dat");
    ofstream mcoutnAgrp1("mcoutnAgrp1.dat");
    ofstream mcoutnAgpPro11("mcoutnAgpPro11.dat");
    ofstream mcoutnAgpPro12("mcoutnAgpPro12.dat");
    ofstream mcoutnAgpPro13("mcoutnAgpPro13.dat");
    ofstream mcoutnAgpPro21("mcoutnAgpPro21.dat");
    ofstream mcoutnAgpPro22("mcoutnAgpPro22.dat");
    ofstream mcoutnAgpPro23("mcoutnAgpPro23.dat");
    ofstream mcoutnAgpPro31("mcoutnAgpPro31.dat");
    ofstream mcoutnAgpPro32("mcoutnAgpPro32.dat");
    ofstream mcoutnAgpPro33("mcoutnAgpPro33.dat");
    ofstream mcoutFg1("mcoutFg1.dat");
    ofstream mcoutFg2("mcoutFg2.dat");
    ofstream mcoutF8("mcoutF8.dat"); 
    ofstream mcoutF6("mcoutF6.dat");
    ofstream mcoutF4("mcoutF4.dat");
    ofstream mcoutF2("mcoutF2.dat");
    ofstream mcoutFpro1("mcoutFpro1.dat");
    ofstream mcoutFpro2("mcoutFpro2.dat");
    ofstream mcoutFpro3("mcoutFpro3.dat");
    ofstream mcoutSpro1("mcoutSpro1.dat");
    ofstream mcoutSpro2("mcoutSpro2.dat");
    ofstream mcoutSpro3("mcoutSpro3.dat");
    ofstream mcoutSSB("mcoutSSB.dat");
    ofstream mcoutnSsb1("mcoutssb1.dat");
    ofstream mcoutnSsb2("mcoutssb2.dat");
    ofstream mcoutnSsb3("mcoutssb3.dat");
    ofstream mcoutnSsb4("mcoutssb4.dat");
    ofstream mcoutnSsb5("mcoutssb5.dat");
    ofstream mcoutnSsb6("mcoutssb6.dat");
    ofstream mcoutnSsb7("mcoutssb7.dat");
    ofstream mcoutnSsb8("mcoutssb8.dat");
    ofstream mcoutRVq("mcoutRVq.dat");
    ofstream mcoutSel1("mcoutSel1.dat");
    ofstream mcoutSel2("mcoutSel2.dat");
    ofstream mcoutSel3("mcoutSel3.dat");
    ofstream mcoutAgeProp("mcoutAgeProp.dat");
    ofstream mcoutR("mcoutR.dat");
    ofstream mcoutrecR("mcoutrecR.dat");
    ofstream mcoutRVpred2log("mcoutLogPredRVadt.dat");
    ofstream mcoutRVpred1log("mcoutLogPredRVjuv.dat");
    ofstream mcoutRVpred2("mcoutPredRVadt.dat");
    ofstream mcoutRVpred1("mcoutPredRVjuv.dat");

    int iseed=1337;
    random_number_generator rng(iseed);
    

TOP_OF_MAIN_SECTION
    arrmblsize=20000000;
    gradient_structure::set_CMPDIF_BUFFER_SIZE(25000000);
    gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
    
PRELIMINARY_CALCS_SECTION
    //Cwtg *= 1.e-3; // convert weight (kg) to model weight (tonnes)  NO model scale changed to 1000s for N
    //Iwtg *= 1.e-3; // convert weight (kg) to model weight (tonnes)
    neval = 0;
    tauSquareLgrp  = square(init_tauLgrp);
  
PROCEDURE_SECTION
  initModelParameters();
  //cout << "initModelParameters done..." << endl;
  popInit();
  //cout << "popInit done..." << endl;
  popDynamics();
  //cout << "popDynamics done..." << endl;
  calc_index_likelihood();
  //cout << "index_likelihood done..." << endl;
  calc_age_likelihood();
  //cout << "age_likelihood done..." << endl;
  calc_rec_prior();
  //cout << "rec_prior done..." << endl;
  calc_M_prior();
  //cout << "M_prior done..." << endl;
  calc_q_prior();
  //cout << "q_prior done..." << endl;
  // Exploitation rates at the end
  if( last_phase() )
    calc_exploitation_rates();
  //cout << "calc_exploitation_rates done..." << endl;
  objFunc  = 0.;
  objFunc  += LikeWeight(1)*indexLikelihood(1) + LikeWeight(2)*indexLikelihood(2);
  objFunc  += LikeWeight(3)*LgrpLikelihood;
  objFunc += recPrior;
  objFunc += mPrior;
  objFunc += qPrior;  

  if ( mceval_phase() )
  {
    projection(); 
    mcoutM1   << column(Mta,2) << endl;
    mcoutM2   << column(Mta,4) << endl;
    mcoutM3   << column(Mta,6) << endl;  
    mcoutnAgrp3    << nAgrp(3)    << endl; 
    mcoutnAgrp2    << nAgrp(2)    << endl;  
    mcoutnAgrp1    << nAgrp(1)    << endl; 
    mcoutnAgpPro11 << nAgpPro1(1) << endl;
    mcoutnAgpPro12 << nAgpPro1(2) << endl;
    mcoutnAgpPro13 << nAgpPro1(3) << endl;
    mcoutnAgpPro21 << nAgpPro2(1) << endl;
    mcoutnAgpPro22 << nAgpPro2(2) << endl;
    mcoutnAgpPro23 << nAgpPro2(3) << endl;
    mcoutnAgpPro31 << nAgpPro3(1) << endl;
    mcoutnAgpPro32 << nAgpPro3(2) << endl;
    mcoutnAgpPro33 << nAgpPro3(3) << endl;
    mcoutFg1   << avgF   << endl;   
    mcoutFg2   << avgF2   << endl;   
    mcoutF8   << elem_prod(column(sel_gta(1),8),Ft)   << endl;   
    mcoutF6   << elem_prod(column(sel_gta(1),6),Ft)   << endl;   
    mcoutF4   << elem_prod(column(sel_gta(1),4),Ft)   << endl;   
    mcoutF2   << elem_prod(column(sel_gta(1),2),Ft)   << endl;   
    mcoutFpro1     << Fpro(1)     << endl;
    mcoutFpro2     << Fpro(2)     << endl;
    mcoutFpro3     << Fpro(3)     << endl;
    mcoutSpro1     << Spro1     << endl;
    mcoutSpro2     << Spro2     << endl;
    mcoutSpro3     << Spro3     << endl;
    mcoutSSB   << ssb    << endl;
    mcoutnSsb1  << column(Sta,1)    << endl;
    mcoutnSsb2  << column(Sta,2)    << endl;
    mcoutnSsb3  << column(Sta,3)    << endl;
    mcoutnSsb4  << column(Sta,4)    << endl;
    mcoutnSsb5  << column(Sta,5)    << endl;
    mcoutnSsb6  << column(Sta,6)    << endl;
    mcoutnSsb7  << column(Sta,7)    << endl;
    mcoutnSsb8  << column(Sta,8)    << endl;
    mcoutRVq   << exp(lnq)*sel_gta(2)(1)(sage,lage) << endl;
    mcoutSel1   << sel_gta(1)(1)(sage,lage)/max(sel_gta(1)(1)(sage,lage)) << endl;
    mcoutSel2   << sel_gta(1)(25)(sage,lage)/max(sel_gta(1)(25)(sage,lage)) << endl;
    mcoutSel3   << sel_gta(1)(29)(sage,lage)/max(sel_gta(1)(29)(sage,lage)) << endl;
    mcoutAgeProp << uCgat(2)  << endl;   
    mcoutR   << column(Nta,1)   << endl;   
    mcoutrecR   << recRate   << endl;   
    mcoutRVpred2log << log(rvpred(2))    << endl;
    mcoutRVpred1log << log(rvpred(1))    << endl;
    mcoutRVpred2 << rvpred(2)    << endl;
    mcoutRVpred1 << rvpred(1)    << endl;

    mcmcHeader = 1;
  }


FUNCTION initModelParameters
  {
  avgR = mfexp( log_AvgR );
  Finit = mfexp( logFinit );
  gamma_R = exp( logit_gamma_R )/(1.+exp(logit_gamma_R));  
  rvq     = mfexp( lnq );
  rvsig = mfexp( log_sig );
  Mta.initialize();
  calc_Mta_Split();  
  calc_sel_gta();
  
  }

FUNCTION calc_Mta_Split
  {
    for( int t=1; t<=6; t++)
    {
      Mta(t)(sage,splitAge(2)-1) = mfexp(log_Mjuv(1));
      Mta(t)(splitAge(2)+1,lage) = mfexp(log_Madt(1));
      Mta(t,splitAge(2)) = (mfexp(log_Mjuv(1)) + mfexp(log_Madt(1)))/2.0;
    }  
    for( int t=7; t<=12; t++)
    {
      Mta(t)(sage,splitAge(2)-1) = mfexp(log_Mjuv(2));
      Mta(t)(splitAge(2)+1,lage) = mfexp(log_Madt(2));
      Mta(t,splitAge(2)) = (mfexp(log_Mjuv(2)) + mfexp(log_Madt(2)))/2.0;
    }  
    for( int t=13; t<=18; t++)
    {
      Mta(t)(sage,splitAge(2)-1) = mfexp(log_Mjuv(3));
      Mta(t)(splitAge(2)+1,lage) = mfexp(log_Madt(3));
      Mta(t,splitAge(2)) = (mfexp(log_Mjuv(3)) + mfexp(log_Madt(3)))/2.0;
    }  
    for( int t=19; t<=24; t++)
    {
      Mta(t)(sage,splitAge(2)-1) = mfexp(log_Mjuv(4));
      Mta(t)(splitAge(2)+1,lage) = mfexp(log_Madt(4));
      Mta(t,splitAge(2)) = (mfexp(log_Mjuv(4)) + mfexp(log_Madt(4)))/2.0;
    }  
    for( int t=25; t<=30; t++)
    {
      Mta(t)(sage,splitAge(2)-1) = mfexp(log_Mjuv(5));
      Mta(t)(splitAge(2)+1,lage) = mfexp(log_Madt(5));
      Mta(t,splitAge(2)) = (mfexp(log_Mjuv(5)) + mfexp(log_Madt(5)))/2.0;
    }  
    for( int t=31; t<=36; t++)
    {
      Mta(t)(sage,splitAge(2)-1) = mfexp(log_Mjuv(6));
      Mta(t)(splitAge(2)+1,lage) = mfexp(log_Madt(6));
      Mta(t,splitAge(2)) = (mfexp(log_Mjuv(6)) + mfexp(log_Madt(6)))/2.0;
    }  
  }

FUNCTION calc_sel_gta
  {
  // 
  int lt, ut;
  sel_gta.initialize();
  S50_gt.initialize();
  S95_gt.initialize();  
  for( int g=1; g<=nFisheries; g++ )
  {
    if( nSelBlocks(g) > 1 )
    {
      // time-varying selectivity in blocks.
      for( int i=1; i<=nSelBlocks(g)-1; i++ )
      {
        lt = tBlock(g,i); ut = tBlock(g,i+1)-1;
        S50_gt(g)(lt,ut) = mfexp( log_S50_gi(g,i) );
        S95_gt(g)(lt,ut) = S50_gt(g)(lt,ut) + mfexp( log_S95_step_gi(g,i) );
      }
      S50_gt(g)(ut+1,nT) = mfexp( log_S50_gi(g,nSelBlocks(g)) );
      S95_gt(g)(ut+1,nT) = S50_gt(g)(ut+1,nT) + mfexp( log_S95_step_gi(g,nSelBlocks(g)) );
    }
    else
    {
      // constant sel: i is vector 1,nT
      S50_gt(g)(1,nT) = mfexp( log_S50_gi(g,1) );
      S95_gt(g)(1,nT) = S50_gt(g) + mfexp( log_S95_step_gi(g,1) );
    }
    // time loop to fill
    for( int t=1; t<=nT; t++ )
    {
      dvariable tmp = log(19.)/( S95_gt(g)(t) - S50_gt(g)(t) );
      sel_gta(g)(t) = 1./( 1. + exp(-tmp*(age - S50_gt(g)(t) ) ) );
    }
  }
  }


FUNCTION popInit
  {
  nLgrp.initialize(); bLgrp.initialize();
  uCgat.initialize(); expNtg.initialize();
  nAgrp.initialize();
  ssb.initialize(); Sta.initialize();
  
  // Initialize abundance in first year assuming each cohort
  // had independent recruitment deviations and age-dependent M.
  // age-1 does not involve M

  Nta(1)(sage) = avgR*mfexp( init_recDevs(sage) );
  
  for( int a=sage+1; a<lage; a++ )
  {
    Nta(1)(a) = avgR*mfexp( init_recDevs(a) - sum( sel_gta(1)(1)(sage,a-1)*Finit + Mta(1)(sage,a-1) ) );
  }
  Nta(1)(lage) = avgR*mfexp( init_recDevs(lage) - sum( sel_gta(1)(1)(sage,lage-1)*Finit + Mta(1)(sage,lage-1) ) );
  Nta(1)(lage) *= 1./( 1.-mfexp(-( sel_gta(1)(1)(lage)*Finit + Mta(1)(lage) ) ) );
  
  // nLgrp(1,1) and bLgrp(1,1) : 
  nLgrp(1,1) = 0.5*Nta(1,splitAge(2));
  for( int a=sage; a<=splitAge(2)-1; a++ ) nLgrp(1,1) += Nta(1,a);
  nLgrp(2,1) = 0.5*Nta(1,splitAge(2));
  for( int a=splitAge(2)+1; a<=lage; a++ ) nLgrp(2,1) += Nta(1,a);
  for( int g=1; g<=2; g++ )  bLgrp(g,1) = nLgrp(g,1)*Iwtg(g,1); 
  
  // nAgrp(,1)  
   for(int a=1; a<=3; a++ ) nAgrp(1,1) = sum(Nta(1)(1,3));
   for(int a=4; a<=5; a++ ) nAgrp(2,1) = sum(Nta(1)(4,5));
   for(int a=6; a<=8; a++ ) nAgrp(3,1) = sum(Nta(1)(6,8));
   
   //SSB in tonnes
   for(int a=sage; a<=lage; a++) Sta(1,a) = Nta(1,a)*matAge(a,1)*waaSep(a,1);
   ssb(1) = sum(Sta(1)(sage,lage)); 

  // Solve catch equations, t=1: returns Zta(1) and Ftg(1)
  int t=1;
  Zta.initialize(); Ft.initialize();
  //Zta(1) = sel_gta(1)(t)*Finit + Mta(1);
  Zta(1) =  Mta(1);
  //cout << "Zta(1) = " << Zta(1) << endl;
  solveBaranov(t,baranovIter);
  //cout << "Zta(1) = " << Zta(1) << endl;
 
  // predicted exploitable abundance for RV (in individuals- need to divide by 1000 to compare to RV index) 
      expNtg(1,1) = 0.5*Nta(1,splitAge(2))*sel_gta(2,1,splitAge(2))*exp(-fracYearSurvey*Zta(1,splitAge(2)));
      expNtg(2,1) = 0.5*Nta(1,splitAge(2))*sel_gta(2,1,splitAge(2))*exp(-fracYearSurvey*Zta(1,splitAge(2)));
      for (int a=sage; a<=splitAge(2)-1; a++) expNtg(1,1) += Nta(1,a)*sel_gta(2,1,a)*exp(-fracYearSurvey*Zta(1,a));
      for (int a=splitAge(2)+1; a<=lage; a++) expNtg(2,1) += Nta(1,a)*sel_gta(2,1,a)*exp(-fracYearSurvey*Zta(1,a));

  // predicted age-proportions uCgat in this year's fishery
  // based on numbers-at-age
      dvariable tmp1 = 0.5*Nta(1,splitAge(2))*sel_gta(1,1,splitAge(2)); 
      tmp1 += sum(elem_prod(Nta(1)(sage,splitAge(2)-1),sel_gta(1)(1)(sage,splitAge(2)-1))); 
      dvariable tmp2 = 0.5*Nta(1,splitAge(2))*sel_gta(1,1,splitAge(2)); 
      tmp2 += sum(elem_prod(Nta(1)(splitAge(2)+1,lage),sel_gta(1)(1)(splitAge(2)+1,lage))); 
      uCgat(1,1) = tmp1/(tmp1+tmp2);
      uCgat(2,1) = tmp2/(tmp1+tmp2);

  }

FUNCTION popDynamics
  {
  int t, g;
  for( int t=2; t<=nT; t++ )
  {
    // Age-2 recruitment
      Nta(t,sage) = mfexp( gamma_R*log(Nta(t-1,sage)) + (1.-gamma_R)*log_AvgR + recDevs(t) );

    // age-3 to age-(lage)
    for( int a=sage+1; a<=lage-1; a++ )
    {
      Nta(t,a) = Nta(t-1,a-1)*exp( -Zta(t-1,a-1) );
    }
    Nta(t,lage) = Nta(t-1,lage-1)*exp( -Zta(t-1,lage-1) ) + 
                  Nta(t-1,lage)*exp( -Zta(t-1,lage) );

    // nLgrp(g,t) and bLgrp(g,t) : 
    nLgrp(1,t) = 0.5*Nta(t,splitAge(2));
    for( int a=sage; a<=splitAge(2)-1; a++ ) nLgrp(1,t) += Nta(t,a);
    nLgrp(2,t) = 0.5*Nta(t,splitAge(2));
    for( int a=splitAge(2)+1; a<=lage; a++ ) nLgrp(2,t) += Nta(t,a);
    for( int g=1; g<=2; g++ )  bLgrp(g,t) = nLgrp(g,t)*Iwtg(g,t); 

   // nAgrp(,1)  
     for(int a=1; a<=3; a++ ) nAgrp(1,t) = sum(Nta(t)(1,3));
     for(int a=4; a<=5; a++ ) nAgrp(2,t) = sum(Nta(t)(4,5));
     for(int a=6; a<=8; a++ ) nAgrp(3,t) = sum(Nta(t)(6,8));
     
   //SSB in tonnes
     for(int a=sage; a<=lage; a++) Sta(t,a) = Nta(t,a)*matAge(a,t)*waaSep(a,t);
     ssb(t) = sum(Sta(t)(sage,lage)); 
     recRate(t) = Nta(t,1)/ssb(t-1);
     

   // Solve catch equations, t: returns Zta(t) and Ftg(t)(1,nFisheries)
    solveBaranov(t,baranovIter);
    
   // predicted exploitable abundance for RV (in individuals- need to divide by 1000 to compare to RV index) 
      expNtg(1,t) = 0.5*Nta(t,splitAge(2))*sel_gta(2,t,splitAge(2))*exp(-fracYearSurvey*Zta(t,splitAge(2)));
      expNtg(2,t) = 0.5*Nta(t,splitAge(2))*sel_gta(2,t,splitAge(2))*exp(-fracYearSurvey*Zta(t,splitAge(2)));
      for (int a=sage; a<=splitAge(2)-1; a++) expNtg(1,t) += Nta(t,a)*sel_gta(2,t,a)*exp(-fracYearSurvey*Zta(t,a));
      for (int a=splitAge(2)+1; a<=lage; a++) expNtg(2,t) += Nta(t,a)*sel_gta(2,t,a)*exp(-fracYearSurvey*Zta(t,a));
   
  // predicted age-proportions uCgat in this year's fishery
  // based on numbers-at-age
      dvariable tmp1 = 0.5*Nta(t,splitAge(2))*sel_gta(1,t,splitAge(2)); 
      tmp1 += sum(elem_prod(Nta(t)(sage,splitAge(2)-1),sel_gta(1)(t)(sage,splitAge(2)-1))); 
      dvariable tmp2 = 0.5*Nta(t,splitAge(2))*sel_gta(1,t,splitAge(2)); 
      tmp2 += sum(elem_prod(Nta(t)(splitAge(2)+1,lage),sel_gta(1)(t)(splitAge(2)+1,lage))); 
      uCgat(1,t) = tmp1/(tmp1+tmp2);
      uCgat(2,t) = tmp2/(tmp1+tmp2);

  } // end t-loop

  if ( last_phase() ) neval+=1;	
  else neval=-2.;  // start neval counter at -2 at start of last phase so it equals admb screen output
  }

  
FUNCTION calc_rec_prior
  {
    recPrior = 0.;
    if( active(init_recDevs) )
        recPrior =  0.5*norm2(init_recDevs)/square(priorSD_R);
    if( active(recDevs) )
        recPrior += 0.5*norm2(recDevs)/square(priorSD_R);        
  }

  
FUNCTION calc_M_prior
  {
  mPrior = 0.;
  // Normal prior for initial M
  if( active( log_Mjuv) )
  {
      mPrior += 0.5*square(mfexp(log_Mjuv(1)) - prior_Minit(1))/square(mpriorSD(1));
      mPrior += 0.5*square(mfexp(log_Madt(1)) - prior_Minit(2))/square(mpriorSD(2));
  }
  }

FUNCTION calc_q_prior
  {
    qPrior = 0.;
    if( active(lnq) ) 
    {
        qPrior += 0.5*square(lnq - logq_prior)/square(logq_priorSD);
    }
  }  


FUNCTION calc_index_likelihood
  {
    int t; int g;
    rvpred.initialize();
    rvres.initialize();
    indexLikelihood.initialize();
    for( g=1; g<=2; g++ )
    {
      validObs.initialize();
      ss.initialize();
      for( t=1; t<=nTretro; t++)
      {
        rvpred(g,t) = expNtg(g,t)*rvq;
        if( Iy(g,t) > 0 )
        {
          rvres(g,t) = log(Iy(g,t)/rvpred(g,t));
          validObs += 1;
          ss += pow( rvres(g,t), 2);
        }  
      }
      indexLikelihood(g) = validObs*log(rvsig(g)) + (ss/(2.0*square(rvsig(g))));
    }  
  }
 
FUNCTION calc_age_likelihood
  {
  dvariable meanDiff;
  dvariable nYearsAges;
  dvariable etaSumSq;
  dvariable mnLL;
  tauSquareLgrp.initialize();

  // Calculate predicted age-proportions for each index gear.
    nYearsAges = 0.;
    etaSumSq   = 0.;
    mnLL       = 0.;
    for ( int t=1; t<=nTretro; t++ )
    {
      dvar_vector res(1,2);
      res.initialize();
      nYearsAges += 1;
      dvar_vector obsPat  = column(ageObsProp,t);
      dvar_vector predPat = ( column(uCgat,t) )/sum( column(uCgat,t) ); //denom not required but doesn't hurt
      
      // Column means for p_at residuals.
      int iRes = 0; dvariable sumRes = 0.;
      for( int a=1; a<=2; a++ )
      {
        if( obsPat(a) > 0 )
        {
          //cout << obsPat(a) << endl;
          iRes   += 1;
          res(iRes) = log( obsPat(a) ) - log( predPat(a) );
          sumRes += res(iRes);
        }
        else
        {
          //cout << obsPat(a) << endl;
          mnLL += -50.*log(1.-predPat(a));
        }
      }
      if( iRes > 0 )
      {
        meanDiff = sumRes/iRes;
        // Proportion residual function.
        etaSumSq += norm2(res(1,iRes)-meanDiff);
      }
    }
    // MLE of variance in age-proportions.
    tauSquareLgrp = etaSumSq/nYearsAges;  // or denom ((maxAge(g)-minAge(g))*nYearsAges) if nages > 2 
    LgrpLikelihood = nYearsAges*log(tauSquareLgrp) + mnLL; // or (maxAge(g)-minAge(g))*nYearsAges if nages > 2
  }



FUNCTION  void solveBaranov( int t, int nIter )
  {
  RETURN_ARRAYS_INCREMENT();
  dvariable f;
  dvariable J;
  dvar_vector Bprime(sage,lage);
  dvar_vector tmp(sage,lage);
  dvar_vector Za(sage,lage); 
  dvar_vector ZaNew(sage,lage);  

  // Initialize Z to current vector of Mta
  ZaNew.initialize(); 
  // Initial approximation of F...
  // Selected biomass
  // Cwtg is average for Lgrp - assume same for all ages within Lgrp
  for (int a=sage; a<=splitAge(2)-1; a++) Bprime(a) = Nta(t,a)*Cwtg(1,t)*sel_gta(1,t,a);
  Bprime(splitAge(2)) = 0.5*Nta(t,splitAge(2))*Cwtg(1,t)*sel_gta(1,t,splitAge(2));
  Bprime(splitAge(2)) += 0.5*Nta(t,splitAge(2))*Cwtg(2,t)*sel_gta(1,t,splitAge(2));
  for (int a=splitAge(2)+1; a<=lage; a++) Bprime(a) = Nta(t,a)*Cwtg(2,t)*sel_gta(1,t,a);
  Ft(t) = Ct(t)/sum( Bprime );
  ZaNew = Mta(t) + sel_gta(1)(t)*Ft(t);
  // refine F for fisheries (surveys not critical)
  for( int i=1; i<=nIter; i++ )
  {
    // Total mortality
    Za=ZaNew; ZaNew=Mta(t);
    // Predicted catch given F
    tmp    = elem_div( elem_prod( Bprime*Ft(t),1.-exp(-Za) ), Za );
    //cout << "predCatch = " << tmp << endl;
    // Function value: difference of pred - obs catch
    f =  sum(tmp) - Ct(t); tmp=0;
    // Jacobian
    dvar_vector tmp1 = elem_div( Bprime, Za );
    dvar_vector tmp2 = elem_prod( sel_gta(1)(t), Za )*Ft(t);
    dvar_vector tmp3 = 1. - mfexp( -Za );
    dvar_vector tmp4 = elem_prod( elem_div( sel_gta(1)(t)*Ft(t), Za ), 1.-exp(-Za)  );
    
    tmp = elem_prod( tmp1, tmp2 + tmp3 - tmp4 );
    J = sum(tmp); tmp=0;
    Ft(t) -= f/J;
    ZaNew += sel_gta(1)(t)*Ft(t); 
    //cout <<"t = " << t<< " iter = "<< i << " f = "<< f <<" J = "<< J << " f/J = " << f/J << endl;   
    //cout <<"iter = "<< i << " Ftg = "<< Ftg(t, 1) << endl;   
  }
  Zta(t) = ZaNew;
  avgF(t) = sum( elem_prod( Nta(t)(a1F1,a2F1),sel_gta(1)(t)(a1F1,a2F1)*Ft(t) ) );
  avgF(t) /= sum(Nta(t)(a1F1,a2F1));
  avgF2(t) = sum( elem_prod( Nta(t)(a1F2,a2F2),sel_gta(1)(t)(a1F2,a2F2)*Ft(t) ) );
  avgF2(t) /= sum(Nta(t)(a1F2,a2F2));
  RETURN_ARRAYS_DECREMENT();
  }


FUNCTION calc_exploitation_rates
  {
  subHt.initialize();
  dvar_matrix subC(1,2,1,nT);  // juv adult
  dvar_vector prop(1,2);
  dvar_vector propC(1,2);
  for( int t=1; t<=nT; t++ )
  {
    prop = column( ageObsProp, t );
    propC = elem_prod(prop,column(Cwtg,t))/sum(elem_prod(prop,column(Cwtg,t)));
    for (int a=1; a<=2; a++) subC(a,t) = propC(a)*Ct(t);
  }  
  subHt = elem_div( subC, bLgrp);
  // Harvest rate matrix (1,2,1,nT)
  }

FUNCTION  double getFforward( const int& t, const int& p, const int& nIter, const dvector& Np, const dvector& Cw);
  {
  double fp;
  double Jp;
  double Fest;
  dvector Bp(sage,lage);
  dvector tmp(sage,lage);
  dvector Zap(sage,lage); 
  dvector ZapNew(sage,lage);  
  
  // Initialize Z to current vector of Mta
  ZapNew.initialize(); 
  // Initial approximation of F...
  // Selected biomass
  for (int a=sage; a<=splitAge(2)-1; a++) Bp(a) = Np(a)*Cw(1)*PRav(a);
  Bp(splitAge(2)) = 0.5*Np(splitAge(2))*Cw(1)*PRav(splitAge(2));
  Bp(splitAge(2)) += 0.5*Np(splitAge(2))*Cw(2)*PRav(splitAge(2));
  for (int a=splitAge(2)+1; a<=lage; a++) Bp(a) = Np(a)*Cw(2)*PRav(a);
  Fest = Cpro(p)/sum( Bp );
  ZapNew = Mpro + PRav*Fest;
  // refine F for fisheries (surveys not critical)
  for( int i=1; i<=nIter; i++ )
  {
    // Total mortality
    Zap=ZapNew; ZapNew=Mpro;
    // Predicted catch given F
    tmp    = elem_div( elem_prod( Bp*Fest,1.-exp(-Zap) ), Zap );
    //cout << "predCatch = " << tmp << endl;
    // Function value: difference of pred - obs catch
    fp =  sum(tmp) - Cpro(p); tmp=0;
    // Jacobian
    dvector tmp1 = elem_div( Bp, Zap );
    dvector tmp2 = elem_prod( PRav, Zap )*Fest;
    dvector tmp3 = 1. - mfexp( -Zap );
    dvector tmp4 = elem_prod( elem_div( PRav*Fest, Zap ), 1.-exp(-Zap)  );
    
    tmp = elem_prod( tmp1, tmp2 + tmp3 - tmp4 );
    Jp = sum(tmp); tmp=0;
    Fest -= fp/Jp;
    ZapNew += PRav*Fest; 
    //cout <<"t = " << t<< " iter = "<< i << " f = "<< f <<" J = "<< J << " f/J = " << f/J << endl;   
    //cout <<"iter = "<< i << " Ftg = "<< Ftg(t, 1) << endl;   
  }
  return Fest;
  }


FUNCTION projection
 {
   int i; int j; int k; int p;
   dmatrix usewate(sage,lage,nT+1,pT);
   double rnyr;
   double maxSel;
   dvector NnT(sage,lage);
   dvector FnT(sage,lage);
   dvector SSB(1,nT);
   dvector rRateUse(nT+1,pT);
   dmatrix CwtUse(1,2,nT+1,pT);

  // randomly pick weight-at-age vectors and catch-weight from last 10 years
   for (j=nT+1; j<=pT; j++)
   {
     k=9.0*randu(rng);
     for (i=sage; i<=lage; i++) usewate(i,j)=waaSep(i,nT-k);
     for (i=1; i<=2; i++) CwtUse(i,j) = Cwtg(i,nT-k);
   }
  
   // use last M time block for M
   Mpro = value(Mta(nT));

  // randomly pick recruitment rate from last 10 years
   for (j=nT+1; j<=pT; j++)
   {
     k=9.0*randu(rng) + 1;
     rRateUse(j) = value(recRate(nT-k+1));
   }

  NnT = value(Nta(nT));
  SSB = value(ssb(1,nT));
  FnT = value(sel_gta(1)(nT)(sage,lage)*Ft(nT));
  
  // calculate fishery selectivity adjusted to a maximum of 1
  maxSel = max(value(sel_gta(1)(nT)(sage,lage)));
  for (i=sage; i<=lage; i++) PRav(i) = value(sel_gta(1,nT,i))/maxSel;
  
//*************************************** Projection 1, Catch = 0.0 t  ************************************

  // PROJECTION 1, Cpro = 0 t
  Fpro(1) = 0.0;
  Npro = 0.0;
  Spro1 = 0.0;
  nAgpPro1 = 0.0;
  
  // first projection year
  Npro(nT+1,sage) = rRateUse(nT+1) * SSB(nT);  
  for (i=sage+1; i<=lage-1; i++) Npro(nT+1,i) = NnT(i-1)*mfexp(-Mpro(i-1)-FnT(i-1));
  Npro(nT+1,lage) =  NnT(lage-1)*mfexp(-Mpro(lage-1)-FnT(lage-1));
  Npro(nT+1,lage) +=  NnT(lage)*mfexp(-Mpro(lage)-FnT(lage));
  for (i=sage; i<=lage; i++) Spro1(nT+1) += Npro(nT+1,i)*usewate(i,nT+1)*matAge(i,nT);
  for (i=1; i<=3; i++) nAgpPro1(1,nT+1) = sum(Npro(nT+1)(1,3));
  for (i=4; i<=5; i++) nAgpPro1(2,nT+1) = sum(Npro(nT+1)(4,5));
  for (i=6; i<=8; i++) nAgpPro1(3,nT+1) = sum(Npro(nT+1)(6,8));
  
  // remaining projection years
  for( j=nT+2; j<=pT; j++)
  {
    Npro(j,sage) = rRateUse(j) *Spro1(j-1);  
    for (i=sage+1; i<=lage-1; i++) Npro(j,i) = Npro(j-1,i-1)*mfexp(-Mpro(i-1));
    Npro(j,lage) =  Npro(j-1,lage-1)*mfexp(-Mpro(lage-1));
    Npro(j,lage) +=  Npro(j-1,lage)*mfexp(-Mpro(lage));
    for (i=sage; i<=lage; i++) Spro1(j) += Npro(j,i)*usewate(i,j)*matAge(i,nT);
    for (i=1; i<=3; i++) nAgpPro1(1,j) = sum(Npro(j)(1,3));
    for (i=4; i<=5; i++) nAgpPro1(2,j) = sum(Npro(j)(4,5));
    for (i=6; i<=8; i++) nAgpPro1(3,j) = sum(Npro(j)(6,8));
  }

 
//*************************************** END PROJECTION1 ************************************

  // PROJECTION 2, Cpro = 100 t
  Fpro(2) = 0.0;
  Npro = 0.0;
  Spro2 = 0.0;
  nAgpPro2 = 0.0;
  
  // first projection year
  Npro(nT+1,sage) = rRateUse(nT+1) *SSB(nT);  
  for (i=sage+1; i<=lage-1; i++) Npro(nT+1,i) = NnT(i-1)*mfexp(-Mpro(i-1)-FnT(i-1));
  Npro(nT+1,lage) =  NnT(lage-1)*mfexp(-Mpro(lage-1)-FnT(lage-1));
  Npro(nT+1,lage) +=  NnT(lage)*mfexp(-Mpro(lage)-FnT(lage));
  for (i=sage; i<=lage; i++) Spro2(nT+1) += Npro(nT+1,i)*usewate(i,nT+1)*matAge(i,nT);
  for (i=1; i<=3; i++) nAgpPro2(1,nT+1) = sum(Npro(nT+1)(1,3));
  for (i=4; i<=5; i++) nAgpPro2(2,nT+1) = sum(Npro(nT+1)(4,5));
  for (i=6; i<=8; i++) nAgpPro2(3,nT+1) = sum(Npro(nT+1)(6,8));
  Fpro(2,nT+1)=getFforward( nT+1, 2, baranovIter, Npro(nT+1), column(CwtUse,nT+1));
   
  // remaining projection years
  for( j=nT+2; j<=pT; j++)
  {
    Npro(j,sage) = rRateUse(j) *Spro2(j-1);  
    for (i=sage+1; i<=lage-1; i++) Npro(j,i) = Npro(j-1,i-1)*mfexp(-Mpro(i-1)-(PRav(i-1)*Fpro(2,j-1)));
    Npro(j,lage) =  Npro(j-1,lage-1)*mfexp(-Mpro(lage-1)-(PRav(lage-1)*Fpro(2,j-1)));
    Npro(j,lage) +=  Npro(j-1,lage)*mfexp(-Mpro(lage)-(PRav(lage)*Fpro(2,j-1)));
    for (i=sage; i<=lage; i++) Spro2(j) += Npro(j,i)*usewate(i,j)*matAge(i,nT);
    for (i=1; i<=3; i++) nAgpPro2(1,j) = sum(Npro(j)(1,3));
    for (i=4; i<=5; i++) nAgpPro2(2,j) = sum(Npro(j)(4,5));
    for (i=6; i<=8; i++) nAgpPro2(3,j) = sum(Npro(j)(6,8));
    Fpro(2,j)=getFforward( j, 2, baranovIter, Npro(j), column(CwtUse,j) );
  }

 
//*************************************** END PROJECTION2 ************************************

  // PROJECTION 3, Cpro = 300 t
  Fpro(3) = 0.0;
  Npro = 0.0;
  Spro3 = 0.0;
  nAgpPro3 = 0.0;
  
  // first projection year
  Npro(nT+1,sage) = rRateUse(nT+1) *SSB(nT);  
  for (i=sage+1; i<=lage-1; i++) Npro(nT+1,i) = NnT(i-1)*mfexp(-Mpro(i-1)-FnT(i-1));
  Npro(nT+1,lage) =  NnT(lage-1)*mfexp(-Mpro(lage-1)-FnT(lage-1));
  Npro(nT+1,lage) +=  NnT(lage)*mfexp(-Mpro(lage)-FnT(lage));
  for (i=sage; i<=lage; i++) Spro3(nT+1) += Npro(nT+1,i)*usewate(i,nT+1)*matAge(i,nT);
  for (i=1; i<=3; i++) nAgpPro3(1,nT+1) = sum(Npro(nT+1)(1,3));
  for (i=4; i<=5; i++) nAgpPro3(2,nT+1) = sum(Npro(nT+1)(4,5));
  for (i=6; i<=8; i++) nAgpPro3(3,nT+1) = sum(Npro(nT+1)(6,8));
  Fpro(3,nT+1)=getFforward( nT+1, 3, baranovIter, Npro(nT+1), column(CwtUse,nT+1));
   
  // remaining projection years
  for( j=nT+2; j<=pT; j++)
  {
    Npro(j,sage) = rRateUse(j) *Spro3(j-1);  
    for (i=sage+1; i<=lage-1; i++) Npro(j,i) = Npro(j-1,i-1)*mfexp(-Mpro(i-1)-(PRav(i-1)*Fpro(3,j-1)));
    Npro(j,lage) =  Npro(j-1,lage-1)*mfexp(-Mpro(lage-1)-(PRav(lage-1)*Fpro(3,j-1)));
    Npro(j,lage) +=  Npro(j-1,lage)*mfexp(-Mpro(lage)-(PRav(lage)*Fpro(3,j-1)));
    for (i=sage; i<=lage; i++) Spro3(j) += Npro(j,i)*usewate(i,j)*matAge(i,nT);
    for (i=1; i<=3; i++) nAgpPro3(1,j) = sum(Npro(j)(1,3));
    for (i=4; i<=5; i++) nAgpPro3(2,j) = sum(Npro(j)(4,5));
    for (i=6; i<=8; i++) nAgpPro3(3,j) = sum(Npro(j)(6,8));
    Fpro(3,j)=getFforward( j, 3, baranovIter, Npro(j), column(CwtUse,j));
  }

 
//*************************************** END PROJECTION3 ************************************

 }


REPORT_SECTION
  {
  // Parameter estimates
  report << "## Parameter estimates" << endl;
  report << "# avgR "<< endl;            report << avgR << endl;
  report << "# Mjuv " << endl;           report << mfexp(log_Mjuv) << endl;
  report << "# Madt " << endl;           report << mfexp(log_Madt) << endl;
  report << "# recDevs " << endl;           report << recDevs << endl;
  report << "# initRecDevs" << endl;     report << init_recDevs << endl;
  report << "# gamma_R" << endl;     report << gamma_R << endl;  
  
  // Selectivity parameters in blocks
  report << "## Selectivity parameters: S50_gi" << endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    for( int j=1; j<=nSelBlocks(g); j++ )
    {
      report << "# S50_"<< g << j << endl; 
      report << mfexp(log_S50_gi(g,j)) << endl;
    }
  }
  report << "## Selectivity parameters: S95_gi" << endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    for( int j=1; j<=nSelBlocks(g); j++ )
    {
      report << "# S95_"<< g << j << endl; 
      report << mfexp(log_S50_gi(g,j))+mfexp(log_S95_step_gi(g,j)) << endl;
    }
  }

  // Derived likelihood components
  report << "## Standard error estimates" << endl;
  report << "# tauIndex" << endl;    report << rvsig << endl;
  report << "# tauLgrp" << endl;     report << sqrt(tauSquareLgrp) << endl;
  report << "# rvq" << endl;          report << exp( lnq ) << endl;
  report << "# rvpred " << endl;          report << rvpred << endl;
  report << "# rvres " << endl;          report << rvres << endl;
  report << "# expNtg " << endl;          report << expNtg << endl;
    
  // Minimization performance
  report << endl;
  report << "## Minimization performance" << endl;
  report << "# indexLikelihood" << endl;  report << indexLikelihood << endl;
  report << "# LgrpLikelihood" << endl;    report << LgrpLikelihood << endl;
  report << "# objFun" << endl;       report << *objective_function_value::pobjfun << endl;
  report << "# maxGrad" << endl;      report << objective_function_value::gmax << endl;
  report << "# exitCode" << endl;     report << iexit << endl;
  report << "# funEvals" << endl;     report << neval << endl;

  // Derived variables
  report << endl;
  report << "## Derived variables" << endl;
  report << "# Catch" << endl;         report << Ct << endl;
  report << "# bLgrp" << endl;        report << bLgrp << endl;
  report << "# nLgrp" << endl;        report << nLgrp << endl;
  report << "# nAgrp" << endl;        report << nAgrp << endl;
  report << "# ssb" << endl;        report << ssb << endl;
  report << "# Sta" << endl;        report << Sta << endl;
  report << "# Mt2" << endl;         report << column(Mta,2) << endl;
  report << "# Mt4" << endl;         report << column(Mta,4) << endl;
  report << "# Mt5" << endl;         report << column(Mta,5) << endl;
  //report << "# subBt" << endl;      report << subBt << endl;
  report << "# subHt" << endl;      report << subHt << endl;
  
  report <<"## Fishing mortality rates"<< endl;
  report << "# avgF"<< endl;        report << avgF << endl;
  report << "# avgF2"<< endl;        report << avgF2 << endl;
  report << "# Ft" << endl;        report << Ft<< endl;
  report << "# Finit" << endl;        report << Finit<< endl;

  report << "# Nta" << endl;        report << Nta << endl;
  report <<"## Scaled RV index)"<< endl;
  report << "# RVscaled"<< endl;   report << Iy/exp(lnq) << endl;
  
  report <<"## Predicted age proportions by fishery (if obs series was used)"<< endl;
  report << "# Catprops " << endl;  report << uCgat << endl;
  //--------------------------------------------------------------------------//
  // Echo other inputs
  //--------------------------------------------------------------------------//
  report << "# firstAge" << endl;        report << sage << endl;
  report << "# plusGroupAge" << endl;    report << lage << endl;
  report << "# nSplitAge" << endl;       report << nSplitAge << endl;
  report << "# splitAge" << endl;        report << splitAge << endl;
  report << "# baranovIter" << endl;     report << baranovIter << endl;
  report << "# ph_M" << endl;  report << ph_M << endl;
  report << "# ph_S50" << endl;  report << ph_S50 << endl;
  report << "# ph_S95" << endl;  report << ph_S95 << endl;
  report << "# ph_recDevs" << endl;  report << ph_recDevs << endl;  
  report << "# a1F1" << endl;  report << a1F1 << endl;
  report << "# a2F1" << endl;  report << a2F1 << endl;
  report << "# a1F2" << endl;  report << a1F2 << endl;
  report << "# a2F2" << endl;  report << a2F2 << endl;
  
 
  // Model parameter priors.
  report << "# prior_initM" << endl;     report << prior_Minit << endl;
  report << "# priorSD_initM" << endl;   report << mpriorSD << endl;

  // Catch Inputs
  report << "# nT" << endl;              report << nT << endl;
  report << "# Land" << endl;      report << Ct << endl;
  report << "# Cwtg" << endl; report << Cwtg << endl;

  // RV Inputs 
  report << "# RVndx" << endl;        report << Iy << endl;
  report << "# Iwtg" << endl;     report << Iwtg << endl;
  report << "# idxLikeWeight" << endl;   report << LikeWeight << endl;
  report << "# fracYearSurvey" << endl;  report << fracYearSurvey << endl;

  // Age Inputs
  report << "# firstAge" << endl;        report << sage << endl;
  report << "# nSplitAge" << endl;       report << nSplitAge << endl;
  report << "# splitAge" << endl;        report << splitAge << endl;
  report << "# plusGroupAge" << endl;    report << lage << endl;
  report << "# ageObsProp " << endl;    report << ageObsProp << endl;
  report << "# init_tauLgrp " << endl;    report << init_tauLgrp << endl;

  report <<"## Selectivity function by fishery"<< endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# sel_"<<g<< endl; report << sel_gta(g) << endl;
  }  
  //--------------------------------------------------------------------------//
  // End echo of inputs                                                       //
  //--------------------------------------------------------------------------//
  }
