#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
 
 # include "admodel.h"
  //adstring model_name;
  //adstring data_file;
  /// for R report 
	#undef R_Report 
	#define R_Report(object) R_report << #object "\n" << object << endl;
	
	#undef log_input
	#define log_input(object) write_input_log << "# " #object "\n" << object << endl;
  ofstream write_input_log("input.log");
  ofstream R_report("pop_R.rep");
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <pop.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_evalout = new ofstream("evalout.prj");;
  ctrl_file.allocate("ctrl_file");
  mat_file.allocate("mat_file");
 ad_comm::change_datafile_name(ctrl_file);        // Read in phases, penalties and priors from "tem.ctl"
  model_name.allocate("model_name");
  data_file.allocate("data_file");
  styr_rec_est.allocate("styr_rec_est");
  endyr_rec_est.allocate("endyr_rec_est");
 nrecs_est = endyr_rec_est-styr_rec_est+1;
  ph_Fdev.allocate("ph_Fdev");
  ph_avg_F.allocate("ph_avg_F");
  ph_recdev.allocate("ph_recdev");
  ph_fydev.allocate("ph_fydev");
  ph_historic_F.allocate("ph_historic_F");
  ph_fish_sel.allocate("ph_fish_sel");
  ph_fish_sel_dlog.allocate("ph_fish_sel_dlog");
                   ph_fish_sel_doublelog = ph_fish_sel_dlog;        
  ph_cubic_sel.allocate("ph_cubic_sel");
  ph_srv1_sel.allocate("ph_srv1_sel");
  ph_srv2_sel.allocate("ph_srv2_sel");
  mprior.allocate("mprior");
 log_mprior = log(mprior);
  cvmprior.allocate("cvmprior");
  ph_m.allocate("ph_m");
  sigrprior.allocate("sigrprior");
  cvsigrprior.allocate("cvsigrprior");
  ph_sigr.allocate("ph_sigr");
  q_srv1prior.allocate("q_srv1prior");
 log_q_srv1prior = log(q_srv1prior);
  cvq_srv1prior.allocate("cvq_srv1prior");
  ph_q_srv1.allocate("ph_q_srv1");
  q_srv2prior.allocate("q_srv2prior");
 log_q_srv2prior = log(q_srv2prior);
  cvq_srv2prior.allocate("cvq_srv2prior");
  ph_q_srv2.allocate("ph_q_srv2");
  yr_catchwt.allocate("yr_catchwt");
  wt_ssqcatch.allocate("wt_ssqcatch");
  wt_ssqcatch2.allocate("wt_ssqcatch2");
  wt_cpue.allocate("wt_cpue");
  wt_srv1.allocate("wt_srv1");
  wt_srv2.allocate("wt_srv2");
  wt_fish_age.allocate("wt_fish_age");
  wt_srv1_age.allocate("wt_srv1_age");
  wt_fish_size.allocate("wt_fish_size");
  wt_srv1_size.allocate("wt_srv1_size");
  wt_srv2_size.allocate("wt_srv2_size");
  wt_rec_var.allocate("wt_rec_var");
  wt_fy_var.allocate("wt_fy_var");
  wt_hf_pen.allocate("wt_hf_pen");
  wt_fmort_reg.allocate("wt_fmort_reg");
  wt_avg_sel.allocate("wt_avg_sel");
  initial_LMR.allocate("initial_LMR");
  yieldratio.allocate("yieldratio");
  fishselopt.allocate("fishselopt");
  fyopt.allocate("fyopt");
  selp_in.allocate(1,3,"selp_in");
  num_yrs_sel_ch.allocate("num_yrs_sel_ch");
  yrs_sel_ch.allocate(1,num_yrs_sel_ch,"yrs_sel_ch");
  sigma_sel_ch.allocate(1,num_yrs_sel_ch,"sigma_sel_ch");
  n_sel_ch_fsh=num_yrs_sel_ch;
  R_Bk.allocate("R_Bk");
  R_Bk_Yrs.allocate(1,R_Bk+1,"R_Bk_Yrs");
  n_yr_nodes.allocate("n_yr_nodes");
  n_age_nodes.allocate("n_age_nodes");
  wt_spl_avg.allocate("wt_spl_avg");
  wt_spl_dome.allocate("wt_spl_dome");
  wt_spl_2d_ages.allocate("wt_spl_2d_ages");
  wt_spl_1d_yrs.allocate("wt_spl_1d_yrs");
  wt_spl_2d_yrs.allocate("wt_spl_2d_yrs");
  scal_yr_nodes.allocate(1,n_yr_nodes);
  scal_age_nodes.allocate(1,n_age_nodes);
 ad_comm::change_datafile_name(data_file);                 // Read data from the data file
  styr.allocate("styr");
  endyr.allocate("endyr");
  recage.allocate("recage");
  nages_D.allocate("nages_D");
  nages_M.allocate("nages_M");
  nlenbins.allocate("nlenbins");
  n_ageage_mat.allocate("n_ageage_mat");
  n_sizeage_mat.allocate("n_sizeage_mat");
  len_bin_labels.allocate(1,nlenbins,"len_bin_labels");
  age_vector.allocate(1,nages_M);
 for (int j=recage;j<=recage+nages_M-1;j++) age_vector(j-recage+1) = j;
  nyrs = endyr - styr + 1;
 styr_rec = (styr - nages_M) + 1;                          // First year of recruitment
 styr_sp  = styr_rec - recage ;                            // First year of spawning biomass  
 endyr_sp = endyr   - recage - 1;                          // endyr year of (main) spawning biomass
  yy.allocate(styr,endyr);
 yy.fill_seqadd(styr,1) ;
  aa.allocate(1,nages_M);
 aa.fill_seqadd(recage,1) ;
  aa_D.allocate(1,nages_D);
 aa_D.fill_seqadd(recage,1) ;
 ph_F50 = 4;
  spawn_fract.allocate("spawn_fract");
 spawn_fract = (spawn_fract - 1) / 12;
  wt.allocate(1,nages_M,"wt");
  obs_catch_early.allocate(styr,yr_catchwt,"obs_catch_early");
  obs_catch_later.allocate(yr_catchwt+1,endyr,"obs_catch_later");
  historic_catch.allocate("historic_catch");
  nyrs_cpue.allocate("nyrs_cpue");
  yrs_cpue.allocate(1,nyrs_cpue,"yrs_cpue");
  obs_cpue.allocate(1,nyrs_cpue,"obs_cpue");
 if (nyrs_cpue>0) mean_obs_cpue = exp(mean(log(obs_cpue))); 
  nyrs_srv1.allocate("nyrs_srv1");
  yrs_srv1.allocate(1,nyrs_srv1,"yrs_srv1");
  obs_srv1_biom.allocate(1,nyrs_srv1,"obs_srv1_biom");
  obs_srv1_se.allocate(1,nyrs_srv1,"obs_srv1_se");
  obs_srv1_lci.allocate(1,nyrs_srv1,"obs_srv1_lci");
  obs_srv1_uci.allocate(1,nyrs_srv1,"obs_srv1_uci");
  nyrs_srv2.allocate("nyrs_srv2");
  yrs_srv2.allocate(1,nyrs_srv2,"yrs_srv2");
  obs_srv2_biom.allocate(1,nyrs_srv2,"obs_srv2_biom");
  obs_srv2_se.allocate(1,nyrs_srv2,"obs_srv2_se");
  obs_srv2_lci.allocate(1,nyrs_srv2,"obs_srv2_lci");
  obs_srv2_uci.allocate(1,nyrs_srv2,"obs_srv2_uci");
  nyrs_fish_age.allocate("nyrs_fish_age");
  yrs_fish_age.allocate(1,nyrs_fish_age,"yrs_fish_age");
  nsamples_fish_age.allocate(1,nyrs_fish_age,"nsamples_fish_age");
  nhauls_fish_age.allocate(1,nyrs_fish_age,"nhauls_fish_age");
  age_age_ind_fsh.allocate(1,nyrs_fish_age,"age_age_ind_fsh");
  oac_fish.allocate(1,nyrs_fish_age,1,nages_D,"oac_fish");
  nmulti_fish_age.allocate(1,nyrs_fish_age);
  nyrs_srv1_age.allocate("nyrs_srv1_age");
  yrs_srv1_age.allocate(1,nyrs_srv1_age,"yrs_srv1_age");
  nsamples_srv1_age.allocate(1,nyrs_srv1_age,"nsamples_srv1_age");
  nhauls_srv1_age.allocate(1,nyrs_srv1_age,"nhauls_srv1_age");
  age_age_ind_srv.allocate(1,nyrs_srv1_age,"age_age_ind_srv");
  oac_srv1.allocate(1,nyrs_srv1_age,1,nages_D,"oac_srv1");
  nmulti_srv1_age.allocate(1,nyrs_srv1_age);
  nyrs_fish_size.allocate("nyrs_fish_size");
  yrs_fish_size.allocate(1,nyrs_fish_size,"yrs_fish_size");
  nsamples_fish_size.allocate(1,nyrs_fish_size,"nsamples_fish_size");
  nhauls_fish_size.allocate(1,nyrs_fish_size,"nhauls_fish_size");
  siz_age_ind_fsh.allocate(1,nyrs_fish_size,"siz_age_ind_fsh");
  osc_fish.allocate(1,nyrs_fish_size,1,nlenbins,"osc_fish");
  nmulti_fish_size.allocate(1,nyrs_fish_size);
  nyrs_srv1_size.allocate("nyrs_srv1_size");
  yrs_srv1_size.allocate(1,nyrs_srv1_size,"yrs_srv1_size");
  nsamples_srv1_size.allocate(1,nyrs_srv1_size,"nsamples_srv1_size");
  nhauls_srv1_size.allocate(1,nyrs_srv1_size,"nhauls_srv1_size");
  siz_age_ind_srv1.allocate(1,nyrs_srv1_size,"siz_age_ind_srv1");
  osc_srv1.allocate(1,nyrs_srv1_size,1,nlenbins,"osc_srv1");
  nyrs_srv2_size.allocate("nyrs_srv2_size");
  yrs_srv2_size.allocate(1,nyrs_srv2_size,"yrs_srv2_size");
  nsamples_srv2_size.allocate(1,nyrs_srv2_size,"nsamples_srv2_size");
  nhauls_srv2_size.allocate(1,nyrs_srv2_size,"nhauls_srv2_size");
  siz_age_ind_srv2.allocate(1,nyrs_srv2_size,"siz_age_ind_srv2");
  osc_srv2.allocate(1,nyrs_srv2_size,1,nlenbins,"osc_srv2");
  sizeage.allocate(1,n_sizeage_mat,1,nages_M,1,nlenbins,"sizeage");
  ageage.allocate(1,n_ageage_mat,1,nages_M,1,nages_D,"ageage");
  eof.allocate("eof");
 cout <<obs_catch_later<<endl;
 ad_comm::change_datafile_name(mat_file);                 // Read data from the data file
  nages_mat.allocate("nages_mat");
  ages_mat.allocate(1,nages_mat,"ages_mat");
  L_tot_na.allocate(1,nages_mat,"L_tot_na");
  L_mat_na.allocate(1,nages_mat,"L_mat_na");
  C_tot_na.allocate(1,nages_mat,"C_tot_na");
  C_mat_na.allocate(1,nages_mat,"C_mat_na");
   if(eof==42) cout<<"The data has been read correctly!";
   else { cout <<"You f'ed up your data file!"<<endl;exit(1); }
   if(wt_rec_var==0) {
     if (ph_sigr>0) {
       cout << "Warning, wt_rec_var is zero, so can't estimate sigr!@"<<endl;
       cout << "turning sigr off "<<endl;
       ph_sigr =-4;
       cout << "hit any key, then enter to continue"<<endl;
       char  xxx; cin >> xxx;
     }
   }
  // reset phases to take out parameters not being used in the estimation of initial numbers at age
  if(fyopt==1) ph_historic_F = -1;
  if(fyopt==2)
    {
     ph_fydev = -1;
     if(historic_catch<0.01) ph_historic_F = -1;
    }
   
   if (fishselopt==1)
     {
       ph_fish_sel_doublelog=-1;
       ph_cubic_sel=-1;
     }
   else if (fishselopt==2)    // double logisic selectivity 
     { 
      ph_fish_sel=-1;
      ph_cubic_sel=-1;
      } 
   else if (fishselopt==3)   // cubic spline selectivity
    {
      isel_npar = n_age_nodes;
      jsel_npar = n_sel_ch_fsh+1;
      ph_fish_sel_doublelog=-1;
      ph_fish_sel=-1;
    }
   else if (fishselopt==4)  // bicubic fishery selectivity 
   {
    isel_npar = n_yr_nodes;
    jsel_npar = n_age_nodes;
    scal_age_nodes.fill_seqadd(0,1./(n_age_nodes-1));
    scal_yr_nodes.fill_seqadd(0,1./(n_yr_nodes-1));
    ph_fish_sel_doublelog=-1;
    ph_fish_sel=-1;
    }
  rescaled_fish_sel.allocate(styr,endyr,1,nages_M);
  rescaled_F.allocate(styr,endyr);
}

void model_parameters::initializationfunction(void)
{
  logm.set_initial_value(log_mprior);
  log_mean_rec.set_initial_value(initial_LMR);
  sigr.set_initial_value(sigrprior);
  a50.set_initial_value(2.5);
  delta.set_initial_value(4.5);
  a502.set_initial_value(6);
  delta2.set_initial_value(1.5);
  a503.set_initial_value(2.5);
  delta3.set_initial_value(4.5);
  a50_srv1.set_initial_value(7.3);
  delta_srv1.set_initial_value(3.8);
  a50_srv2.set_initial_value(7.3);
  delta_srv2.set_initial_value(3.8);
  selp.set_initial_value(selp_in);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  Sp_Biom.allocate(styr,endyr,"Sp_Biom");
  #ifndef NO_AD_INITIALIZE
    Sp_Biom.initialize();
  #endif
  log_mean_rec.allocate(1,R_Bk,1,"log_mean_rec");
  sigr.allocate(1,R_Bk,0.3,10,ph_sigr,"sigr");
  a50.allocate(ph_fish_sel,"a50");
  delta.allocate(ph_fish_sel,"delta");
  a502.allocate(ph_fish_sel,"a502");
  delta2.allocate(ph_fish_sel,"delta2");
  a503.allocate(ph_fish_sel,"a503");
  delta3.allocate(ph_fish_sel,"delta3");
  selp.allocate(1,3,0,n_sel_ch_fsh,ph_fish_sel_doublelog,"selp");
  expa50.allocate("expa50");
  #ifndef NO_AD_INITIALIZE
  expa50.initialize();
  #endif
  expa502.allocate("expa502");
  #ifndef NO_AD_INITIALIZE
  expa502.initialize();
  #endif
  fish_sel1.allocate(1,nages_M,"fish_sel1");
  #ifndef NO_AD_INITIALIZE
    fish_sel1.initialize();
  #endif
  fish_sel2.allocate(1,nages_M,"fish_sel2");
  #ifndef NO_AD_INITIALIZE
    fish_sel2.initialize();
  #endif
  fish_sel3.allocate(1,nages_M,"fish_sel3");
  #ifndef NO_AD_INITIALIZE
    fish_sel3.initialize();
  #endif
  fish_sel4.allocate(1,nages_M,"fish_sel4");
  #ifndef NO_AD_INITIALIZE
    fish_sel4.initialize();
  #endif
  fish_sel.allocate(styr,endyr,1,nages_M,"fish_sel");
  #ifndef NO_AD_INITIALIZE
    fish_sel.initialize();
  #endif
  sel_par.allocate(1,jsel_npar,1,isel_npar,ph_cubic_sel,"sel_par");
  log_fish_sel.allocate(styr,endyr,1,nages_M,"log_fish_sel");
  #ifndef NO_AD_INITIALIZE
    log_fish_sel.initialize();
  #endif
  a50_srv1.allocate(ph_srv1_sel,"a50_srv1");
  delta_srv1.allocate(ph_srv1_sel,"delta_srv1");
  srv1_sel.allocate(1,nages_M,"srv1_sel");
  #ifndef NO_AD_INITIALIZE
    srv1_sel.initialize();
  #endif
  a50_srv2.allocate(ph_srv2_sel,"a50_srv2");
  delta_srv2.allocate(ph_srv2_sel,"delta_srv2");
  srv2_sel.allocate(1,nages_M,"srv2_sel");
  #ifndef NO_AD_INITIALIZE
    srv2_sel.initialize();
  #endif
  log_avg_F.allocate(ph_avg_F,"log_avg_F");
  log_F_devs.allocate(styr,endyr,-15.,15.,ph_Fdev,"log_F_devs");
  Fmort.allocate(styr,endyr,"Fmort");
  #ifndef NO_AD_INITIALIZE
    Fmort.initialize();
  #endif
  Z.allocate(styr,endyr,1,nages_M,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  F.allocate(styr,endyr,1,nages_M,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  S.allocate(styr,endyr,1,nages_M,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  mat_a50.allocate(0.,20.,1,"mat_a50");
  mat_delta.allocate(0.,2.,1,"mat_delta");
  L_pmat.allocate(1,nages_mat,"L_pmat");
  #ifndef NO_AD_INITIALIZE
    L_pmat.initialize();
  #endif
  C_pmat.allocate(1,nages_mat,"C_pmat");
  #ifndef NO_AD_INITIALIZE
    C_pmat.initialize();
  #endif
  Pred_pmat.allocate(1,nages_mat,"Pred_pmat");
  #ifndef NO_AD_INITIALIZE
    Pred_pmat.initialize();
  #endif
  p_mature.allocate(1,nages_M,"p_mature");
  #ifndef NO_AD_INITIALIZE
    p_mature.initialize();
  #endif
  wt_mature.allocate(1,nages_M,"wt_mature");
  #ifndef NO_AD_INITIALIZE
    wt_mature.initialize();
  #endif
  Like_L_vec.allocate(1,nages_mat,"Like_L_vec");
  #ifndef NO_AD_INITIALIZE
    Like_L_vec.initialize();
  #endif
  Like_C_vec.allocate(1,nages_mat,"Like_C_vec");
  #ifndef NO_AD_INITIALIZE
    Like_C_vec.initialize();
  #endif
  Like_L.allocate("Like_L");
  #ifndef NO_AD_INITIALIZE
  Like_L.initialize();
  #endif
  Like_C.allocate("Like_C");
  #ifndef NO_AD_INITIALIZE
  Like_C.initialize();
  #endif
  zero_pen_mat.allocate("zero_pen_mat");
  #ifndef NO_AD_INITIALIZE
  zero_pen_mat.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
 k=0.00001;
  natmortv.allocate(1,nages_M,"natmortv");
  #ifndef NO_AD_INITIALIZE
    natmortv.initialize();
  #endif
  fydev.allocate(1,nages_M-2,-10,10,ph_fydev,"fydev");
  historic_F.allocate(ph_historic_F,"historic_F");
  ehc.allocate("ehc");
  #ifndef NO_AD_INITIALIZE
  ehc.initialize();
  #endif
  log_rec_dev.allocate(styr,endyr,-10.,10.,ph_recdev,"log_rec_dev");
  natage.allocate(styr,endyr,1,nages_M,"natage");
  #ifndef NO_AD_INITIALIZE
    natage.initialize();
  #endif
  batage.allocate(styr,endyr,1,nages_M,"batage");
  #ifndef NO_AD_INITIALIZE
    batage.initialize();
  #endif
  catage.allocate(styr,endyr,1,nages_M,"catage");
  #ifndef NO_AD_INITIALIZE
    catage.initialize();
  #endif
  pred_catch_early.allocate(styr,yr_catchwt,"pred_catch_early");
  #ifndef NO_AD_INITIALIZE
    pred_catch_early.initialize();
  #endif
  pred_catch_later.allocate(yr_catchwt+1,endyr,"pred_catch_later");
  #ifndef NO_AD_INITIALIZE
    pred_catch_later.initialize();
  #endif
  log_q_srv1.allocate(ph_q_srv1,"log_q_srv1");
  log_q_srv2.allocate(ph_q_srv2,"log_q_srv2");
  cv_cpue.allocate(-1,"cv_cpue");
  logm.allocate(ph_m,"logm");
  q_cpue.allocate("q_cpue");
  #ifndef NO_AD_INITIALIZE
  q_cpue.initialize();
  #endif
  pred_cpue.allocate(1,nyrs_cpue,"pred_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_cpue.initialize();
  #endif
  pred_srv1.allocate(1,nyrs_srv1,"pred_srv1");
  #ifndef NO_AD_INITIALIZE
    pred_srv1.initialize();
  #endif
  pred_srv2.allocate(1,nyrs_srv2,"pred_srv2");
  #ifndef NO_AD_INITIALIZE
    pred_srv2.initialize();
  #endif
  eac_fish.allocate(1,nyrs_fish_age,1,nages_D,"eac_fish");
  #ifndef NO_AD_INITIALIZE
    eac_fish.initialize();
  #endif
  eac_srv1.allocate(1,nyrs_srv1_age,1,nages_D,"eac_srv1");
  #ifndef NO_AD_INITIALIZE
    eac_srv1.initialize();
  #endif
  esc_fish.allocate(1,nyrs_fish_size,1,nlenbins,"esc_fish");
  #ifndef NO_AD_INITIALIZE
    esc_fish.initialize();
  #endif
  esc_srv1.allocate(1,nyrs_srv1_size,1,nlenbins,"esc_srv1");
  #ifndef NO_AD_INITIALIZE
    esc_srv1.initialize();
  #endif
  esc_srv2.allocate(1,nyrs_srv2_size,1,nlenbins,"esc_srv2");
  #ifndef NO_AD_INITIALIZE
    esc_srv2.initialize();
  #endif
  effn_fish_age.allocate(1,nyrs_fish_age,"effn_fish_age");
  #ifndef NO_AD_INITIALIZE
    effn_fish_age.initialize();
  #endif
  sdnr_fish_age.allocate(1,nyrs_fish_age,"sdnr_fish_age");
  #ifndef NO_AD_INITIALIZE
    sdnr_fish_age.initialize();
  #endif
  effn_fish_size.allocate(1,nyrs_fish_size,"effn_fish_size");
  #ifndef NO_AD_INITIALIZE
    effn_fish_size.initialize();
  #endif
  sdnr_fish_size.allocate(1,nyrs_fish_size,"sdnr_fish_size");
  #ifndef NO_AD_INITIALIZE
    sdnr_fish_size.initialize();
  #endif
  effn_srv1_age.allocate(1,nyrs_srv1_age,"effn_srv1_age");
  #ifndef NO_AD_INITIALIZE
    effn_srv1_age.initialize();
  #endif
  sdnr_srv1_age.allocate(1,nyrs_srv1_age,"sdnr_srv1_age");
  #ifndef NO_AD_INITIALIZE
    sdnr_srv1_age.initialize();
  #endif
  effn_srv1_size.allocate(1,nyrs_srv1_size,"effn_srv1_size");
  #ifndef NO_AD_INITIALIZE
    effn_srv1_size.initialize();
  #endif
  sdnr_srv1_size.allocate(1,nyrs_srv1_size,"sdnr_srv1_size");
  #ifndef NO_AD_INITIALIZE
    sdnr_srv1_size.initialize();
  #endif
  effn_srv2_size.allocate(1,nyrs_srv2_size,"effn_srv2_size");
  #ifndef NO_AD_INITIALIZE
    effn_srv2_size.initialize();
  #endif
  sdnr_srv2_size.allocate(1,nyrs_srv2_size,"sdnr_srv2_size");
  #ifndef NO_AD_INITIALIZE
    sdnr_srv2_size.initialize();
  #endif
  tot_biom.allocate(styr,endyr,"tot_biom");
  q_srv1.allocate("q_srv1");
  q_srv2.allocate("q_srv2");
  #ifndef NO_AD_INITIALIZE
  q_srv2.initialize();
  #endif
  pred_rec.allocate(styr,endyr,"pred_rec");
  expl_rate.allocate(styr,endyr,"expl_rate");
  #ifndef NO_AD_INITIALIZE
    expl_rate.initialize();
  #endif
  avg_rec.allocate("avg_rec");
  spbiom_trend.allocate("spbiom_trend");
  Depletion.allocate("Depletion");
  #ifndef NO_AD_INITIALIZE
  Depletion.initialize();
  #endif
  spawn_biom.allocate(styr,endyr,"spawn_biom");
  natmort.allocate("natmort");
  #ifndef NO_AD_INITIALIZE
  natmort.initialize();
  #endif
  LMR.allocate("LMR");
  cigar.allocate("cigar");
  q2.allocate("q2");
  nattymort.allocate("nattymort");
  mF50.allocate(0.01,1.,ph_F50,"mF50");
  mF40.allocate(0.01,1.,ph_F50,"mF40");
  mF35.allocate(0.01,1.,ph_F50,"mF35");
  F50.allocate("F50");
  F40.allocate("F40");
  F35.allocate("F35");
  SB0.allocate("SB0");
  #ifndef NO_AD_INITIALIZE
  SB0.initialize();
  #endif
  SBF50.allocate("SBF50");
  #ifndef NO_AD_INITIALIZE
  SBF50.initialize();
  #endif
  SBF40.allocate("SBF40");
  #ifndef NO_AD_INITIALIZE
  SBF40.initialize();
  #endif
  SBF35.allocate("SBF35");
  #ifndef NO_AD_INITIALIZE
  SBF35.initialize();
  #endif
  sprpen.allocate("sprpen");
  #ifndef NO_AD_INITIALIZE
  sprpen.initialize();
  #endif
  Nspr.allocate(1,4,1,nages_M,"Nspr");
  #ifndef NO_AD_INITIALIZE
    Nspr.initialize();
  #endif
  surv_like.allocate(1,3,"surv_like");
  #ifndef NO_AD_INITIALIZE
    surv_like.initialize();
  #endif
  cpue_like.allocate("cpue_like");
  #ifndef NO_AD_INITIALIZE
  cpue_like.initialize();
  #endif
  age_like.allocate(1,6,"age_like");
  #ifndef NO_AD_INITIALIZE
    age_like.initialize();
  #endif
  offset.allocate(1,6,"offset");
  #ifndef NO_AD_INITIALIZE
    offset.initialize();
  #endif
  rec_like.allocate("rec_like");
  #ifndef NO_AD_INITIALIZE
  rec_like.initialize();
  #endif
  rec_like_blk.allocate(1,R_Bk,"rec_like_blk");
  #ifndef NO_AD_INITIALIZE
    rec_like_blk.initialize();
  #endif
  norm2_recdevs.allocate(1,R_Bk,"norm2_recdevs");
  #ifndef NO_AD_INITIALIZE
    norm2_recdevs.initialize();
  #endif
  szcnt_recdevs.allocate(1,R_Bk,"szcnt_recdevs");
  #ifndef NO_AD_INITIALIZE
    szcnt_recdevs.initialize();
  #endif
  fy_like.allocate("fy_like");
  #ifndef NO_AD_INITIALIZE
  fy_like.initialize();
  #endif
  hf_pen.allocate("hf_pen");
  #ifndef NO_AD_INITIALIZE
  hf_pen.initialize();
  #endif
  norm2_fydevs.allocate("norm2_fydevs");
  #ifndef NO_AD_INITIALIZE
  norm2_fydevs.initialize();
  #endif
  szcnt_fydevs.allocate("szcnt_fydevs");
  #ifndef NO_AD_INITIALIZE
  szcnt_fydevs.initialize();
  #endif
  ssqcatch.allocate("ssqcatch");
  #ifndef NO_AD_INITIALIZE
  ssqcatch.initialize();
  #endif
  F_mort_regularity.allocate("F_mort_regularity");
  #ifndef NO_AD_INITIALIZE
  F_mort_regularity.initialize();
  #endif
  avg_sel_penalty.allocate("avg_sel_penalty");
  #ifndef NO_AD_INITIALIZE
  avg_sel_penalty.initialize();
  #endif
  spline_pen.allocate(1,5,"spline_pen");
  #ifndef NO_AD_INITIALIZE
    spline_pen.initialize();
  #endif
  priors.allocate(1,5,"priors");
  #ifndef NO_AD_INITIALIZE
    priors.initialize();
  #endif
  sigr_prior.allocate(1,R_Bk,"sigr_prior");
  #ifndef NO_AD_INITIALIZE
    sigr_prior.initialize();
  #endif
  Like.allocate("Like");
  #ifndef NO_AD_INITIALIZE
  Like.initialize();
  #endif
  obj_fun.allocate("obj_fun");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  pred_catch.allocate(styr,endyr,"pred_catch");
  #ifndef NO_AD_INITIALIZE
    pred_catch.initialize();
  #endif
  obs_catch.allocate(styr,endyr,"obs_catch");
  #ifndef NO_AD_INITIALIZE
    obs_catch.initialize();
  #endif
  N_proj.allocate(endyr+1,endyr+15,1,nages_M,"N_proj");
  #ifndef NO_AD_INITIALIZE
    N_proj.initialize();
  #endif
  FABC_proj.allocate("FABC_proj");
  #ifndef NO_AD_INITIALIZE
  FABC_proj.initialize();
  #endif
  FABC_tot_proj.allocate(1,nages_M,"FABC_tot_proj");
  #ifndef NO_AD_INITIALIZE
    FABC_tot_proj.initialize();
  #endif
  FOFL_proj.allocate("FOFL_proj");
  #ifndef NO_AD_INITIALIZE
  FOFL_proj.initialize();
  #endif
  FOFL_tot_proj.allocate(1,nages_M,"FOFL_tot_proj");
  #ifndef NO_AD_INITIALIZE
    FOFL_tot_proj.initialize();
  #endif
  ABC.allocate("ABC");
  B40.allocate("B40");
  OFL.allocate("OFL");
  #ifndef NO_AD_INITIALIZE
  OFL.initialize();
  #endif
  Z_proj.allocate(1,nages_M,"Z_proj");
  #ifndef NO_AD_INITIALIZE
    Z_proj.initialize();
  #endif
  ZOFL_proj.allocate(1,nages_M,"ZOFL_proj");
  #ifndef NO_AD_INITIALIZE
    ZOFL_proj.initialize();
  #endif
  S_proj.allocate(1,nages_M,"S_proj");
  #ifndef NO_AD_INITIALIZE
    S_proj.initialize();
  #endif
  catage_proj.allocate(endyr+1,endyr+15,1,nages_M,"catage_proj");
  #ifndef NO_AD_INITIALIZE
    catage_proj.initialize();
  #endif
  catage_proj_OFL.allocate(endyr+1,endyr+15,1,nages_M,"catage_proj_OFL");
  #ifndef NO_AD_INITIALIZE
    catage_proj_OFL.initialize();
  #endif
  pred_catch_proj.allocate(endyr+1,endyr+15,"pred_catch_proj");
  #ifndef NO_AD_INITIALIZE
    pred_catch_proj.initialize();
  #endif
  pred_catch_proj_OFL.allocate(endyr+1,endyr+15,"pred_catch_proj_OFL");
  #ifndef NO_AD_INITIALIZE
    pred_catch_proj_OFL.initialize();
  #endif
  spawn_biom_proj.allocate(endyr+1,endyr+15,"spawn_biom_proj");
  tot_biom_proj.allocate(endyr+1,endyr+15,"tot_biom_proj");
  stdev_rec.allocate("stdev_rec");
  #ifndef NO_AD_INITIALIZE
  stdev_rec.initialize();
  #endif
  FOFL.allocate("FOFL");
  #ifndef NO_AD_INITIALIZE
  FOFL.initialize();
  #endif
  FABC.allocate("FABC");
  #ifndef NO_AD_INITIALIZE
  FABC.initialize();
  #endif
  FOFL2.allocate("FOFL2");
  #ifndef NO_AD_INITIALIZE
  FOFL2.initialize();
  #endif
  FABC2.allocate("FABC2");
  #ifndef NO_AD_INITIALIZE
  FABC2.initialize();
  #endif
}

void model_parameters::userfunction(void)
{
  obj_fun =0.0;
  ofstream& evalout= *pad_evalout;
  l=l+1;                                                       // Initiate counter for random seeds in projection
  Get_Selectivity();                                           // Call function to get selectivities
  Get_Maturity();                                              // Call function to get maturity proportions at age
  Get_Mortality_Rates();                                       // Call function to get fishing and natural mortality
  Get_First_Year();                                            // Call function to get first year numbers at age
  Get_Numbers_At_Age();                                        // Call function to get numbers at age per year
  Get_Catch_at_Age();                                          // Call function to get catch at age per year
  Get_Predicted_Values();                                      // Get predicted values for catch, survbio, age and size comps
  if (last_phase())
  {
    Get_Dependent_Vars();                                      // Solve for dependent variables like total bio, recruitment etc.
    Compute_SPR_Rates();                                       // Compute f40 etc.
    Get_Population_Projection();                               // Get 15 year population projection
  }
  Evaluate_Objective_Function();                               // Minimize objective function value
  if (mceval_phase())                                          // For outputting MCMC simulations in text format
  {
     evalout<<sigr<<" "<<q_srv1<<" "<<q_srv2<<" "<<F40<<" "<<natmort<<" "<<" "<<ABC<<" "<<obj_fun<<" "<<tot_biom<<" "<<log_rec_dev<<" "<<spawn_biom<<" "<<log_mean_rec<<" "<<spawn_biom_proj<<" "<<pred_catch_proj<<" "<<N_proj(endyr+1,1)<<" "<<N_proj(endyr+2,1)<<" "<<N_proj(endyr+3,1)<<" "<<N_proj(endyr+4,1)<<" "<<N_proj(endyr+5,1)<<" "<<N_proj(endyr+6,1)<<" "<<N_proj(endyr+7,1)<<" "<<N_proj(endyr+8,1)<<" "<<N_proj(endyr+9,1)<<" "<<N_proj(endyr+10,1)<<" "<<tot_biom_proj(endyr+1)<<" "<<tot_biom_proj(endyr+2)<<" "<<pred_srv1<<" "<<endl;
  }
}

void model_parameters::Get_Selectivity(void)
{
  ofstream& evalout= *pad_evalout;
  switch (fishselopt)
  {
    case 1:
    {
    // Fishery Selectivity
      expa50 = mfexp(a50); 
      expa502 = mfexp(a503); 
      for (j=1;j<=nages_M;j++)  {
        fish_sel1(j) = 1./(1. + mfexp(-2.944438979*(double(j)-a502)/delta2));
        fish_sel3(j)=(pow(j/expa50,expa50/(0.5*(sqrt(square(expa50)+4*square(delta))-expa50)))*mfexp((expa50-j)/(0.5*(sqrt(square(expa50)+4*square(delta))-expa50))));
        fish_sel4(j)=(pow(j/expa502,expa502/(0.5*(sqrt(square(expa502)+4*square(delta3))-expa502)))*mfexp((expa502-j)/(0.5*(sqrt(square(expa502)+4*square(delta3))-expa502)))); }
      fish_sel2 = (fish_sel1 + fish_sel3)/2; 
      fish_sel3=fish_sel3/max(fish_sel3);
      fish_sel4=fish_sel4/max(fish_sel4);
      for (iyr=styr; iyr<=endyr; iyr++)
				if (iyr <=1976)
					fish_sel(iyr) = fish_sel1;
				else if (iyr<=1995)
					fish_sel(iyr) = fish_sel2;
				else if (iyr<=2006)
					fish_sel(iyr) = fish_sel3;
				else 
					fish_sel(iyr) = fish_sel4;
			// cout <<fish_sel<<endl;exit(1);
      break;
    }
    case 2: // Double logistic
    {
      /*
      sel_p1_fsh(k)  = mfexp(logsel_p1_fsh(k));
      sel_p3_fsh(k)  = mfexp(logsel_p3_fsh(k));
      */
      int isel_ch_tmp = 0 ; // base year 
      dvariable p1 = selp(1,isel_ch_tmp);
      dvariable p2 = selp(2,isel_ch_tmp);
      dvariable p3 = selp(3,isel_ch_tmp);
      dvariable i1 = p1 + p2;
      dvariable i2 = p1 + i1 + p3;
      isel_ch_tmp++;
      for (i=styr;i<=endyr;i++)
      {
        if (i==yrs_sel_ch(isel_ch_tmp)) 
        {
          p1 = selp(1,isel_ch_tmp);
          p2 = selp(2,isel_ch_tmp);
          p3 = selp(3,isel_ch_tmp);
          i1 = p1 + p2;
          i2 = p1 + i1 + p3;
          if (isel_ch_tmp<n_sel_ch_fsh)
            isel_ch_tmp++;
        }
        fish_sel(i) = exp( ( -log(1.0 + mfexp(-2.9444389792/p1 * ( age_vector - i1) )) +
               log(1. - 1./(1.0 + mfexp(-2.9444389792/p3 * ( age_vector - i2))) ) )+0.102586589) ; // constant at end is log(0.95*0.95)
       // cout << p1 << " "<<p2<<" "<<p3<<endl<<i1<<" "<<i2<<endl<<age_vector<<endl<<fish_sel(i)<<endl;exit(1);
      }
    break;
    }
    case 3: // cubic spline, varying between blocks
    {
      int nbins;
      nbins = n_sel_ch_fsh + 1;
      ivector binstart(1,nbins);
      binstart(1) = styr;
      for (i=2;i<=nbins;i++) {binstart(i) = yrs_sel_ch(i-1);}
      int bincount;
      bincount = 1;
      for (j=1;j<nbins;j++)
      {
        for (i=binstart(j);i<binstart(j+1);i++)
        {
           log_fish_sel(i)=cubic_spline(sel_par(j));
        }
      }
      for (i=binstart(nbins);i<=endyr;i++)
        {
           log_fish_sel(i) = cubic_spline(sel_par(j));
        }
    fish_sel = mfexp(log_fish_sel);   // convert to unlogged space;     
    break;
    }
    case 4: // bicubic spline
    {
    bicubic_spline(scal_yr_nodes,scal_age_nodes,sel_par,log_fish_sel);
    fish_sel = mfexp(log_fish_sel);   // convert to unlogged space;
    break;
    }
  }
  for (j=1;j<=nages_M;j++)
    srv1_sel(j) = 1./(1. + mfexp(-2.944438979*(double(j)-a50_srv1)/delta_srv1));
  srv2_sel = srv1_sel;
}

void model_parameters::Get_Maturity(void)
{
  ofstream& evalout= *pad_evalout;
  for (int i=1;i<=nages_mat;i++){
   if(L_tot_na(i)>0){L_pmat(i) = L_mat_na(i)/L_tot_na(i);}
   if(C_tot_na(i)>0){C_pmat(i) = C_mat_na(i)/C_tot_na(i);}       
   Pred_pmat(i) = 1/(1+exp(-1.0*mat_delta*(ages_mat(i)-mat_a50)));}
  for (int i=1;i<=nages_M;i++){
  p_mature(i) = 1/(1+exp(-1.0*mat_delta*((i+1)-mat_a50)));
  wt_mature(i) = wt(i)*p_mature(i)/2;}
}

void model_parameters::Get_Mortality_Rates(void)
{
  ofstream& evalout= *pad_evalout;
  natmort = mfexp(logm);                                       // setting natural mortality to arithmetic scale
  if(ph_m>0) nattymort=natmort; else nattymort=log_mean_rec(R_Bk);
  Fmort = mfexp(log_avg_F +  log_F_devs);                      // setting fishing mortaltiy to arithmetic scale
  for (iyr=styr; iyr<=endyr; iyr++)
	  F(iyr)     = Fmort(iyr)*fish_sel(iyr);
	/*
  for (iyr=styr; iyr<=1976; iyr++)
    F(iyr) = Fmort(iyr) * fish_sel1;                           // Getting fully selected fishing mortality
  for (iyr=1977; iyr<=1995; iyr++)
    F(iyr) = Fmort(iyr) * fish_sel2;                           // Getting fully selected fishing mortality
  for (iyr=1996; iyr<=2006; iyr++)
    F(iyr) = Fmort(iyr) * fish_sel3;                           // Getting fully selected fishing mortality Z = F + natmort;
  for (iyr=2007; iyr<=endyr; iyr++)
    F(iyr) = Fmort(iyr) * fish_sel4;                           // Getting fully selected fishing mortality Z = F + natmort;
	*/
  Z = F + natmort;                                             // Fully selected total mortality
  S = mfexp(-1.0*Z);                                           // Fully selected survival
}

void model_parameters::Get_First_Year(void)
{
  ofstream& evalout= *pad_evalout;
  //int itmp;
  if(fyopt==1) {                       // Option 1 -- stohastic numbers at age in first year
  for (j=2;j<nages_M;j++) {
    //itmp = styr+1-j;
    natage(styr,j) = mfexp(log_mean_rec(1) - natmort * double(j-1)+ fydev(j-1)); 
  }
  natage(styr,nages_M) = mfexp(log_mean_rec(1) - natmort * (nages_M-1)) / (1. - exp(-natmort));
  norm2_fydevs = norm2(fydev);
  szcnt_fydevs = size_count(fydev);
  }
  if(fyopt==2) {                     //  Option 2 -- eq first year natage with historic catch
    natage(styr,1) = mfexp(log_mean_rec(1));
    for (j=2;j<=nages_M;j++)
      natage(styr,j) = natage(styr,j-1)*mfexp(-(natmort + historic_F*fish_sel(styr,j-1)));
    natage(styr,nages_M) /= 1-mfexp(-(historic_F*fish_sel(styr,nages_M)+natmort));   // Plus group for first year
    if (historic_catch > 0.) {  // estimate the historical catch
      ehc = 0;
      for (j=1;j<=nages_M;j++)
          {
          ehc += natage(styr,j)*wt(j)*(historic_F*fish_sel(styr,j))*
                (1.0-mfexp(-(historic_F*fish_sel(styr,j)+natmort)))/(historic_F*fish_sel(styr,j)+natmort);
          }
        }
  }  
}

void model_parameters::Get_Numbers_At_Age(void)
{
  ofstream& evalout= *pad_evalout;
  for (b=1;b<=R_Bk;b++) {
    for(i=R_Bk_Yrs(b);i<R_Bk_Yrs(b+1);i++){
      natage(i,1) = mfexp(log_rec_dev(i) + log_mean_rec(b));
      natage(i+1)(2,nages_M) = ++elem_prod(natage(i)(1,nages_M-1),S(i)(1,nages_M-1));       // Following year
      natage(i+1,nages_M) += natage(i,nages_M)*S(i,nages_M);
      Sp_Biom(i) = natage(i) * wt_mature;
    }
  }
  natage(endyr,1) = mfexp(log_rec_dev(endyr) + log_mean_rec(R_Bk)); 
  Sp_Biom(endyr) = elem_prod(natage(endyr),pow(S(endyr),spawn_fract)) * wt_mature;  //Right way, old way was: Sp_Biom(endyr) = natage(endyr)* wt_mature;
  for (b=1;b<R_Bk;b++) {
    norm2_recdevs(b) = norm2(log_rec_dev(R_Bk_Yrs(b),(R_Bk_Yrs(b+1)-1)));
    szcnt_recdevs(b) = size_count(log_rec_dev(R_Bk_Yrs(b),(R_Bk_Yrs(b+1)-1)));
  }
  norm2_recdevs(R_Bk) = norm2(log_rec_dev(R_Bk_Yrs(R_Bk),R_Bk_Yrs(R_Bk+1)));
  szcnt_recdevs(R_Bk) = size_count(log_rec_dev(R_Bk_Yrs(R_Bk),R_Bk_Yrs(R_Bk+1)));
}

void model_parameters::Get_Catch_at_Age(void)
{
  ofstream& evalout= *pad_evalout;
  pred_catch_early.initialize();
  pred_catch_later.initialize();
  for (iyr=styr;iyr<=yr_catchwt;iyr++) {
    catage(iyr) = elem_div(elem_prod(elem_prod(natage(iyr),F(iyr)),(1.-S(iyr))),Z(iyr));
    pred_catch_early(iyr) = catage(iyr)*wt;
  }
  for (iyr=yr_catchwt+1;iyr<=endyr;iyr++) {
    catage(iyr) = elem_div(elem_prod(elem_prod(natage(iyr),F(iyr)),(1.-S(iyr))),Z(iyr));
    pred_catch_later(iyr) = catage(iyr)*wt;
  }
}

void model_parameters::Get_Predicted_Values(void)
{
  ofstream& evalout= *pad_evalout;
  offset.initialize();
  q_srv1         = exp(log_q_srv1);                            // Survey catchability at arithmetic scale
  q_srv2         = exp(log_q_srv2);                            // Survey catchability at arithmetic scale
  for (i=1;i<=nyrs_srv1;i++)
    pred_srv1(i) = q_srv1 * (natage(yrs_srv1(i))*elem_prod(srv1_sel,wt));   // Predicted Survey biomass
  if(nyrs_srv2>0) { for (i=1;i<=nyrs_srv2;i++)
    pred_srv2(i) = q_srv2 * (natage(yrs_srv2(i))*elem_prod(srv2_sel,wt));   // Predicted Survey biomass
  }
  for (i=1;i<=nyrs_fish_age;i++) {
   eac_fish(i)  = catage(yrs_fish_age(i))/sum(catage(yrs_fish_age(i))) * ageage(age_age_ind_fsh(i));
   nmulti_fish_age(i) = sqrt(nsamples_fish_age(i));
   effn_fish_age(i) = (1-eac_fish(i))*eac_fish(i)/norm2(oac_fish(i)-eac_fish(i));
   sdnr_fish_age(i) = sdnr(eac_fish(i),oac_fish(i),double(nmulti_fish_age(i)));
  }
  for (i=1; i<=nyrs_fish_age; i++) {
   oac_fish(i)/=sum(oac_fish(i));
   offset(1) -= nmulti_fish_age(i) *((oac_fish(i) + 0.00001)*log(oac_fish(i) + 0.00001)); 
  }
  for (i=1;i<=nyrs_srv1_age;i++) {
   eac_srv1(i)  = elem_prod(srv1_sel,natage(yrs_srv1_age(i)))/(natage(yrs_srv1_age(i)) * srv1_sel)* ageage(age_age_ind_srv(i));
   nmulti_srv1_age(i) = sqrt(nsamples_srv1_age(i));
   effn_srv1_age(i) = (1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i));
   sdnr_srv1_age(i) = sdnr(eac_srv1(i),oac_srv1(i),double(nmulti_srv1_age(i)));
  }
  for (i=1; i<=nyrs_srv1_age; i++) {
   oac_srv1(i)/=sum(oac_srv1(i));
   offset(2) -= nmulti_srv1_age(i)*((oac_srv1(i) + 0.00001)*log(oac_srv1(i) + 0.00001));
  }
  for (i=1;i<=nyrs_fish_size;i++) {
   esc_fish(i)  = catage(yrs_fish_size(i))/sum(catage(yrs_fish_size(i))) * sizeage(siz_age_ind_fsh(i));
   nmulti_fish_size(i) = nhauls_fish_size(i);
   effn_fish_size(i) = (1-esc_fish(i))*esc_fish(i)/norm2(osc_fish(i)-esc_fish(i));
   sdnr_fish_size(i) = sdnr(esc_fish(i),osc_fish(i),double(nmulti_fish_size(i)));
  }
  nmulti_fish_size = nmulti_fish_size/max(nmulti_fish_size)*100;
  for (i=1; i<=nyrs_fish_size; i++) {
   osc_fish(i)/=sum(osc_fish(i));
   offset(3) -= nmulti_fish_size(i)*((osc_fish(i) + 0.00001)*log(osc_fish(i) + 0.00001));
  }
  for (i=1;i<=nyrs_srv1_size;i++) {
   esc_srv1(i)  = elem_prod(srv1_sel,natage(yrs_srv1_size(i))) /(natage(yrs_srv1_size(i)) * srv1_sel)* sizeage(siz_age_ind_srv1(i));
   effn_srv1_size(i) = (1-esc_srv1(i))*esc_srv1(i)/norm2(osc_srv1(i)-esc_srv1(i));
   //sdnr_srv1_size(i) = sdnr(esc_srv1(i),osc_srv1(i),double(effn_srv1_size(i)));
  }
  for (i=1; i<=nyrs_srv1_size; i++) {
   osc_srv1(i)/=sum(osc_srv1(i));
   offset(4) -= effn_srv1_size(i)*((osc_srv1(i) + 0.00001)*log(osc_srv1(i) + 0.00001));
  }
  if(nyrs_srv2>0) {
   for ( i=1;i<=nyrs_srv2_size;i++) {
    esc_srv2(i)  = elem_prod(srv2_sel,natage(yrs_srv2_size(i))) /(natage(yrs_srv2_size(i)) * srv2_sel)* sizeage(siz_age_ind_srv2(i));
    effn_srv2_size(i) = (1-esc_srv2(i))*esc_srv2(i)/norm2(osc_srv2(i)-esc_srv2(i));
    //sdnr_srv2_size(i) = sdnr(esc_srv2(i),osc_srv2(i),double(effn_srv2_size(i)));
   }
   for (i=1; i<=nyrs_srv2_size; i++) {
    osc_srv2(i)/=sum(osc_srv2(i));
    offset(5) -= effn_srv2_size(i)*((osc_srv2(i) + 0.00001)*log(osc_srv2(i) + 0.00001));
   }
  }
  if (nyrs_cpue>0)
  {
    int yy;
    for (i=1;i<=nyrs_cpue;i++) 
    {
      yy = yrs_cpue(i);
      pred_cpue(i) = wt*elem_div(elem_prod(natage(yy),fish_sel4),Z(yy)); 
    } 
    q_cpue = mean_obs_cpue/mfexp(mean(log(pred_cpue)));
    pred_cpue *=  q_cpue;
  } 
  pred_catch(styr,yr_catchwt) = pred_catch_early;
  pred_catch(yr_catchwt+1,endyr) = pred_catch_later;
  obs_catch(styr,yr_catchwt) = obs_catch_early;
  obs_catch(yr_catchwt+1,endyr) = obs_catch_later;
  if(ph_q_srv2>0) q2=mfexp(log_q_srv2); else q2=mfexp(log_q_srv1);
  cigar= sigr(R_Bk);
  LMR = log_mean_rec(R_Bk);
}

void model_parameters::Get_Dependent_Vars(void)
{
  ofstream& evalout= *pad_evalout;
  for (i=styr;i<=endyr;i++) {
    pred_rec(i) = natage(i,1);                                 // Setting up results based on estimated paramters
    tot_biom(i) = wt * natage(i);                              // Total biomass results
    expl_rate(i) = pred_catch(i)/tot_biom(i);                  // Setting up results based on estimated paramters
    spawn_biom(i) = Sp_Biom(i) ;                               // Spawning biomass result
  }
  avg_rec        = mean(pred_rec);
  Depletion      = spawn_biom(endyr)/spawn_biom(styr);         // 1-Depletion
  spbiom_trend   = spawn_biom(endyr)/spawn_biom(endyr-1);
}

void model_parameters::Compute_SPR_Rates(void)
{
  ofstream& evalout= *pad_evalout;
  SB0=0.;
  SBF50=0.;
  SBF40=0.;
  SBF35=0.;
  dvar_vector seltmp = fish_sel(endyr);
  // Scale F-spr rates to be on full-selected values
  // Ianelli commented this out...shouldn't need right?
  F50  = mF50*max(seltmp);
  F40  = mF40*max(seltmp);
  F35  = mF35*max(seltmp);
  /*
  */
  for (i=1;i<=4;i++)
    Nspr(i,1)=1.;
  for (j=2;j<nages_M;j++) {
    Nspr(1,j)=Nspr(1,j-1)*mfexp(-1.*natmort);
    Nspr(2,j)=Nspr(2,j-1)*mfexp(-1.*(natmort+mF50*seltmp(j-1)));
    Nspr(3,j)=Nspr(3,j-1)*mfexp(-1.*(natmort+mF40*seltmp(j-1)));
    Nspr(4,j)=Nspr(4,j-1)*mfexp(-1.*(natmort+mF35*seltmp(j-1)));
  }
  Nspr(1,nages_M)=Nspr(1,nages_M-1)*mfexp(-1.*natmort)/(1.-mfexp(-1.*natmort));
  Nspr(2,nages_M)=Nspr(2,nages_M-1)*mfexp(-1.* (natmort+mF50*seltmp(nages_M-1)))/(1.-mfexp(-1.*(natmort+mF50*seltmp(nages_M))));
  Nspr(3,nages_M)=Nspr(3,nages_M-1)*mfexp(-1.* (natmort+mF40*seltmp(nages_M-1)))/ (1.-mfexp(-1.*(natmort+mF40*seltmp(nages_M))));
  Nspr(4,nages_M)=Nspr(4,nages_M-1)*mfexp(-1.* (natmort+mF35*seltmp(nages_M-1)))/ (1.-mfexp(-1.*(natmort+mF35*seltmp(nages_M))));
  for (j=1;j<=nages_M;j++) {
   // Kill them off till (spawn_fract)
    SB0    += Nspr(1,j)*wt_mature(j)*mfexp(-spawn_fract*natmort);
    SBF50  += Nspr(2,j)*wt_mature(j)*mfexp(-spawn_fract*(natmort+mF50*seltmp(j)));
    SBF40  += Nspr(3,j)*wt_mature(j)*mfexp(-spawn_fract*(natmort+mF40*seltmp(j)));
    SBF35  += Nspr(4,j)*wt_mature(j)*mfexp(-spawn_fract*(natmort+mF35*seltmp(j)));
  }
  sprpen    = 100.*square(SBF50/SB0-0.5);
  sprpen   += 100.*square(SBF40/SB0-0.4);
  sprpen   += 100.*square(SBF35/SB0-0.35);
  B40= SBF40*mean(pred_rec(1979,endyr-recage));
}

void model_parameters::Get_Population_Projection(void)
{
  ofstream& evalout= *pad_evalout;
  int k;
  if(mceval_phase()) {
    stdev_rec = sqrt(norm2(value(log_rec_dev(1977 + recage,endyr - recage))-mean(value(log_rec_dev(1977 + recage,endyr - recage))))/(size_count(value(log_rec_dev(1977 + recage,endyr - recage))) - 1));
    k=round(value(stdev_rec) * 10000);
    N_proj(endyr+1,1) = mfexp(value(log(mean(value(pred_rec(1977+recage,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l)));
    cout<<stdev_rec<<" "<<k<<" "<<l<<" "<<endl;
  }
  else {
    N_proj(endyr+1,1)= value(mean(pred_rec(1977+recage,endyr-recage))); }
  for (j=1;j<nages_M-1;j++)
    N_proj(endyr+1,j+1) = natage(endyr,j) * S(endyr,j);
  N_proj(endyr+1,nages_M) = natage(endyr,nages_M-1) * S(endyr,nages_M-1) + natage(endyr,nages_M) * S(endyr,nages_M);
  tot_biom_proj(endyr+1) = N_proj(endyr+1) * wt;
  spawn_biom_proj(endyr+1) = elem_prod(N_proj(endyr+1),pow(mfexp(-yieldratio * FABC_tot_proj-natmort),spawn_fract)) * wt_mature;
  for (i=endyr+1;i<=endyr+15;i++) {
   // F ABC 
    if (spawn_biom_proj(i)/B40 > 1.) {
      FABC_proj = F40;
      FOFL_proj = F35; }
    else {
      FABC_proj = F40 * (spawn_biom_proj(i)/B40 - 0.05)/(1 - 0.05); 
      FOFL_proj = F35 * (spawn_biom_proj(i)/B40 - 0.05)/(1 - 0.05); }
    for (j=1;j<=nages_M;j++) {  
      FOFL_tot_proj(j) = fish_sel(endyr,j) * FOFL_proj;
      FABC_tot_proj(j) = fish_sel(endyr,j) * FABC_proj;
      Z_proj(j) = FABC_tot_proj(j) + natmort;
      ZOFL_proj(j) = FOFL_tot_proj(j) + natmort;
      S_proj(j) = mfexp(-1.0 * Z_proj(j));
    }
   // Catch 
    for (j=1;j<=nages_M;j++) { 
      catage_proj(i,j) = yieldratio * N_proj(i,j) * FABC_tot_proj(j) / Z_proj(j) * (1.-S_proj(j));
      catage_proj_OFL(i,j) = yieldratio * N_proj(i,j) * FOFL_tot_proj(j) / ZOFL_proj(j) * (1. - mfexp(-ZOFL_proj(j)));
    }
    pred_catch_proj(i) = catage_proj(i) * wt / yieldratio;
    pred_catch_proj_OFL(i) = catage_proj_OFL(i) * wt / yieldratio;
   // Next year's abundance
    if (i < endyr + 15) {
     if (mceval_phase()) {
       stdev_rec = sqrt(norm2(value(log_rec_dev(1977+recage,endyr-recage))-mean(value(log_rec_dev(1977+recage,endyr-recage))))/(size_count(value(log_rec_dev(1977+recage,endyr-recage)))-1));
       k=round(value(spawn_biom(endyr)*10000))+i;
       k=k+i;
       N_proj(i+1,1)= mfexp((log(mean(value(pred_rec(1977+recage,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l))); }
     else {
       N_proj(i+1,1)= value(mean(pred_rec(1979,endyr-recage))); }
     for (j=1; j<nages_M-1;j++) {
       N_proj(i+1,j+1) = N_proj(i,j) * mfexp(-yieldratio*FABC_tot_proj(j)-natmort); }
     N_proj(i+1,nages_M) = N_proj(i,nages_M-1) * mfexp(-yieldratio*FABC_tot_proj(nages_M-1)-natmort)+ N_proj(i,nages_M) * mfexp(-yieldratio*FABC_tot_proj(nages_M)-natmort);
     spawn_biom_proj(i+1) = elem_prod(N_proj(i+1),pow(mfexp(-yieldratio*FABC_tot_proj-natmort),spawn_fract)) * wt_mature;
     tot_biom_proj(i+1) = N_proj(i+1)*wt;
     }
    }
  if (spawn_biom_proj(endyr+1)/B40 > 1.) {
      FABC = F40;
      FOFL = F35; 
      FABC2 = F40;
      FOFL2 = F35;
  }
  else {
      FABC = F40 * (spawn_biom_proj(endyr+1)/B40 - 0.05)/(1 - 0.05); 
      FOFL = F35 * (spawn_biom_proj(endyr+1)/B40 - 0.05)/(1 - 0.05);  
      FABC2 = F40 * (spawn_biom_proj(endyr+2)/B40 - 0.05)/(1 - 0.05); 
      FOFL2 = F35 * (spawn_biom_proj(endyr+2)/B40 - 0.05)/(1 - 0.05);
  }
  OFL = pred_catch_proj_OFL(endyr+1);
  ABC = pred_catch_proj(endyr+1);
}

void model_parameters::Evaluate_Objective_Function(void)
{
  ofstream& evalout= *pad_evalout;
  Like.initialize();
  Catch_Like();                                                // Catch biomass likelihood (lognormal)
  Surv_Like();                                                 // Trawl survey biomass likelihood (lognormal)
  Size_Age_Like();                                             // Age/Size composition likelihood (multinomial)
  Maturity_Like();                                             // Maturity proportion likelihood (binomial)
  Calc_Priors();                                               // Prior penalties for estimated parameters
  Rec_Like();                                                  // Penalty function for recruitment
  if(fyopt==1)  Fy_Like();                                     // Penalty function for first year deviations
  if(fyopt==2)  Hf_pen();                                      // Penalty function for historic catch   
  F_Like();                                                    // Penalty function for fishing mortality deviations
  Like              += ssqcatch ;
  Like              += sum(surv_like);
  Like              += cpue_like;
  Like              += sum(age_like);
  obj_fun           += Like;
  obj_fun           += Like_L;
  obj_fun           += Like_C;
  obj_fun           += sum(priors);
  obj_fun           += wt_rec_var * rec_like;
  if(fyopt==1)  
      obj_fun           += wt_fy_var * fy_like;
  if(fyopt==2)  
      obj_fun           += wt_hf_pen * hf_pen;
  if(fishselopt==3 || fishselopt==4 )
      obj_fun           += sum(spline_pen);
  if(active(log_F_devs))
    obj_fun         += F_mort_regularity;
  if (current_phase()<3)
      obj_fun       += norm2(F);
  if (active(mF50)&&last_phase())
    obj_fun         += sprpen;   // To solve for the F40 etc.     
}

void model_parameters::Catch_Like(void)
{
  ofstream& evalout= *pad_evalout;
  ssqcatch.initialize();
  ssqcatch  +=  wt_ssqcatch *norm2(log(obs_catch_early+.00001)-log(pred_catch_early+.00001));
  ssqcatch  +=  wt_ssqcatch2 *norm2(log(obs_catch_later+.00001)-log(pred_catch_later+.00001));
}

void model_parameters::Surv_Like(void)
{
  ofstream& evalout= *pad_evalout;
  surv_like.initialize();
  cpue_like.initialize();
  for (i=1; i<=nyrs_srv1; i++)
   //surv_like(1) += square(obs_srv1_biom(i)-pred_srv1(i) )/ (2.*square(obs_srv1_se(i)));
   surv_like(1) += square((log(obs_srv1_biom(i))-log(pred_srv1(i)) ))/ (2.*square(obs_srv1_se(i)/obs_srv1_biom(i)));
  if(nyrs_srv2>0) { 
    for (i=1; i<=nyrs_srv2; i++)
     surv_like(2) += square((log(obs_srv2_biom(i))-log(pred_srv2(i)) ))/ (2.*square(obs_srv2_se(i)/obs_srv2_biom(i)));
  }
  surv_like(1) *= wt_srv1 ;  
  surv_like(2) *= wt_srv2 ;  
  if (nyrs_cpue>0)
    cpue_like  = norm2(log(obs_cpue)-log(pred_cpue)) / (2.*cv_cpue*cv_cpue); // likelihood for fishery cpue
}

void model_parameters::Size_Age_Like(void)
{
  ofstream& evalout= *pad_evalout;
  age_like.initialize();
  for (i=1; i <= nyrs_fish_age; i++)
    age_like(1) -= nmulti_fish_age(i) * ((oac_fish(i) + 0.00001) * log(eac_fish(i) + 0.00001));
  age_like(1)   -= offset(1);
  for (i=1; i <= nyrs_srv1_age; i++)
    age_like(2) -= nmulti_srv1_age(i) * ((oac_srv1(i) + 0.00001) * log(eac_srv1(i) + 0.00001));
  age_like(2)   -= offset(2);
  for (i=1; i <= nyrs_fish_size; i++)
    age_like(3) -= nmulti_fish_size(i) * ((osc_fish(i) + 0.00001) * log(esc_fish(i) + 0.00001));
  age_like(3)   -= offset(3);
  for (i=1; i <= nyrs_srv1_size; i++)
    age_like(4) -= effn_srv1_size(i) * ((osc_srv1(i) + 0.00001) * log(esc_srv1(i) + 0.00001));
  age_like(4)   -= offset(4);
 if(nyrs_srv2>0) {
  for (i=1; i <= nyrs_srv2_size; i++)
    age_like(5) -= effn_srv2_size(i) * ((osc_srv2(i) + 0.00001) * log(esc_srv2(i) + 0.00001));
  age_like(5)   -= offset(5);
 }
  age_like(1) *= wt_fish_age;
  age_like(2) *= wt_srv1_age;
  age_like(3) *= wt_fish_size;
  age_like(4) *= wt_srv1_size;
  age_like(5) *= wt_srv2_size;
}

void model_parameters::Maturity_Like(void)
{
  ofstream& evalout= *pad_evalout;
  Like_L.initialize();
  Like_C.initialize();
  for (int i=1;i<=nages_mat;i++){
   if(L_tot_na(i)>0){
    Like_L_vec(i) = (log(L_mat_na(i)+k)+gammln(L_mat_na(i)+k)+log((L_tot_na(i)-L_mat_na(i))+k)+gammln((L_tot_na(i)-L_mat_na(i))+k)-log(L_tot_na(i)+k)-gammln(L_tot_na(i)+k)-(L_mat_na(i)+k)*log(Pred_pmat(i)+k)-((L_tot_na(i)-L_mat_na(i))+k)*log((1-Pred_pmat(i))+k));}
   if(C_tot_na(i)>0){
    Like_C_vec(i) = (log(C_mat_na(i)+k)+gammln(C_mat_na(i)+k)+log((C_tot_na(i)-C_mat_na(i))+k)+gammln((C_tot_na(i)-C_mat_na(i))+k)-log(C_tot_na(i)+k)-gammln(C_tot_na(i)+k)-(C_mat_na(i)+k)*log(Pred_pmat(i)+k)-((C_tot_na(i)-C_mat_na(i))+k)*log((1-Pred_pmat(i))+k));}}
  zero_pen_mat = 1/(1+exp(-mat_delta*(0-mat_a50)));
  Like_L = sum(Like_L_vec);
  Like_C = sum(Like_C_vec) + 1000*zero_pen_mat;
}

void model_parameters::Calc_Priors(void)
{
  ofstream& evalout= *pad_evalout;
  priors.initialize();
  sigr_prior.initialize();
  spline_pen.initialize();
  dvar_matrix trans_log_fish_sel = trans(log_fish_sel);
  dvariable s = 0.0;
  dvar_vector df1;
  dvar_vector df2;
  if (active(sigr)){
    for(b=1;b<=R_Bk;b++){
      sigr_prior(b) = square(log(sigr(b)/sigrprior))/(2.*square(cvsigrprior));
    }
    priors(1) = sum(sigr_prior);
  }
  if (active(log_q_srv1))
    priors(2)    = square(log_q_srv1-log_q_srv1prior)/(2.*square(cvq_srv1prior));
  if (active(logm))
    priors(3)    = square(logm-log(mprior))/(2.*square(cvmprior));
  if (active(log_q_srv2))
    priors(4)    = square(log_q_srv2-log_q_srv2prior)/(2.*square(cvq_srv2prior));
  if (fishselopt==3 || fishselopt==4)  // penalize the dome-shape for fishery selectivity splines  
   {
    for (i=styr;i<=endyr;i++)
     {
       for (j=1;j<=nages_M-1;j++)
        {
         if (log_fish_sel(i,j)>log_fish_sel(i,j+1))
           {
             spline_pen(1) += wt_spl_dome*square(log_fish_sel(i,j)-log_fish_sel(i,j+1));  // penalize the dome shape
           }
         if (wt_spl_1d_yrs>0 || wt_spl_2d_yrs>0 )
           {
              df1 = first_difference(trans_log_fish_sel(j));
              spline_pen(4) +=  wt_spl_1d_yrs/(endyr-styr+1)*df1*df1;                       // the penalty for interannual variation (across years)
           }
         if (wt_spl_2d_yrs>0)
          {   
            df2 = first_difference(df1);                                     
            spline_pen(5) += wt_spl_2d_yrs/(endyr-styr+1)*df2*df2;                         // the penalty for smoothness over time
          }
         }
        if (wt_spl_1d_yrs>0 || wt_spl_2d_yrs>0 )
           {
              df1 = first_difference(trans_log_fish_sel(nages_M));  
              spline_pen(4) +=  wt_spl_1d_yrs/(endyr-styr+1)*df1*df1;
           }
        if (wt_spl_2d_yrs>0)
           {   
              df2 = first_difference(df1);                                     
              spline_pen(5) += wt_spl_2d_yrs/(endyr-styr+1)*df2*df2;                         
           }     
       s = mean(log_fish_sel(i));
       spline_pen(2) += wt_spl_avg*s*s;                                 // penalize deviations from 0 (across ages)
       dvar_vector df3 = first_difference(first_difference(log_fish_sel(i)));
       spline_pen(3) += wt_spl_2d_ages/nages_M*df3*df3;                                   // penalize 2nd deriviative (smoothness across ages)
      }
    }
}

void model_parameters::Rec_Like(void)
{
  ofstream& evalout= *pad_evalout;
  rec_like.initialize();
  rec_like_blk.initialize();
  for (b=1;b<=R_Bk;b++) {
    rec_like_blk(b) = norm2_recdevs(b)/(2*square(sigr(b))) + szcnt_recdevs(b)*log(sigr(b));
  }
  rec_like = sum(rec_like_blk);
}

void model_parameters::Fy_Like(void)
{
  ofstream& evalout= *pad_evalout;
  fy_like.initialize();
  fy_like = norm2_fydevs/(2*square(sigr(1))) + szcnt_fydevs*log(sigr(1)); 
}

void model_parameters::Hf_pen(void)
{
  ofstream& evalout= *pad_evalout;
  hf_pen.initialize();  
  hf_pen = square(ehc-historic_catch);
}

void model_parameters::F_Like(void)
{
  ofstream& evalout= *pad_evalout;
  F_mort_regularity.initialize();
  if(active(log_F_devs))
    F_mort_regularity  = wt_fmort_reg * norm2(log_F_devs);
}

dvar_vector model_parameters::cubic_spline(const dvar_vector& spline_coffs)
{
  ofstream& evalout= *pad_evalout;
  {
  RETURN_ARRAYS_INCREMENT();
  int nodes=size_count(spline_coffs);
  dvector ia(1,nodes);
  dvector fa(1,nages_M);
  ia.fill_seqadd(0,1./(nodes-1));
  fa.fill_seqadd(0,1./(nages_M-1));
  vcubic_spline_function ffa(ia,spline_coffs);
  RETURN_ARRAYS_DECREMENT();
  return(ffa(fa));
  }
}

double model_parameters::round(double r)
{
  ofstream& evalout= *pad_evalout;
    return double((r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5)); 
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  // rescale F and fishery sel so that fishery selective has max of one
  if(fishselopt==3 || fishselopt==4 )
  {
    rescaled_F = value(mfexp(log_avg_F + log_F_devs));
    for (i=styr;i<=endyr;i++)
      {
        rescaled_fish_sel(i) = value(fish_sel(i));
        rescaled_F(i) = rescaled_F(i)*max(rescaled_fish_sel(i));
        rescaled_fish_sel(i) = rescaled_fish_sel(i)/max(rescaled_fish_sel(i));
      }
  }    
   cout<<"-------------Finished: "<<current_phase()<<" "<<Like<<" "<<age_like<<endl;
  if (last_phase())
    write_proj();
  report<<"~~~~~~~ Executive Summary Material ~~~~~~~"<<endl;report<<endl;
  report<<"     Model name"<<endl;
  report<<model_name<<endl;
  report<<"     .dat file"<<endl;
  report<<data_file<<endl;
  report<<"     Number parameters estimated"<<endl;
  report<<initial_params::nvarcalc()<<endl;
  report<<"     TotalBiomass for "<<endyr+1<<endl;
  report<<N_proj(endyr+1)*wt<<endl;
  report<<"     TotalBiomass for "<<endyr+2<<endl;
  report<<N_proj(endyr+2)*wt<<endl;
  report<<"     Female_Spawning Biomass for "<<endyr+1<<endl;
  report<<spawn_biom_proj(endyr+1)<<endl;
  report<<"     Female_Spawning_Biomass for "<<endyr+2<<endl;
  report<<spawn_biom_proj(endyr+2)<<endl;
  report<<"     B_zero"     <<endl;
  report<<SB0*mean(pred_rec(1977+recage,endyr-recage))<<endl;
  report<<"     B_40"<<endl;
  report<<B40<<endl;
  report<<"     B_35"<<endl;
  report<<SBF35*mean(pred_rec(1977+recage,endyr-recage))<<endl;
  report<<"     F_40"<<endl;
  report<<F40<<endl;
  report<<"     F_35"<<endl;
  report<<F35<<endl;
  report<<"     F_ABC for "<<endyr+1<<endl;
  report<<FABC<<endl;
  report<<"     F_ABC for "<<endyr+2<<endl;
  report<<FABC2<<endl;
  report<<"     ABC for "<<endyr+1<<endl;
  report<<pred_catch_proj(endyr+1)<<endl;
  report<<"     ABC for "<<endyr+2<<endl;
  report<<pred_catch_proj(endyr+2)<<endl;
  report<<"     F_OFL for "<<endyr+1<<endl;
  report<<FOFL<<endl;
  report<<"     F_OFL for "<<endyr+2<<endl;
  report<<FOFL2<<endl;
  report<<"     OFL for "<<endyr+1<<endl;
  report<<OFL<<endl; 
  report<<"     OFL for "<<endyr+2<<endl;
  report<<pred_catch_proj_OFL(endyr+2)<<endl; 
  report<<"     Total likelihood"<<endl;
  report<<obj_fun<<endl;
  report<<"     Data likelihood"<<endl;
  report<<Like<<endl;report<<endl;
  report<<"~~~~~~~ Some more key parameter estimates ~~~~~~~"<<endl;report<<endl;
  report<<"   q_trawl   "<<endl;
  report<<q_srv1<<endl;
  report<<"   nat_mort  "<<endl;
  report<<natmort<<endl;
  report<<"  sigr   "<<endl;  
  report<<sigr<<endl;  
  report<<"   log_mean_rec"<<endl;
  report<<log_mean_rec<<endl;
  report<<"   q_alt  "<<endl;
  report<<q_srv2<<endl;report<<endl;
  report<<"~~~~~~~ Rest of model/data output ~~~~~~~"<<endl;report<<endl;
  report << "Year "<< yy <<endl;
  report << "Pred_Catch "<< pred_catch_early<<pred_catch_later <<endl;
  report << "Obs_Catch "<< obs_catch_early<<obs_catch_later <<endl;report<<endl;
  report << "Catch_at_age "<<aa <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<catage(i) <<endl; report<<endl;
  report << "Numbers "<<aa <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<natage(i) <<endl; report<<endl;
  if (nyrs_cpue>0) {
    report <<"Years_CPUE: "     <<yrs_cpue <<endl; 
    report <<"Predicted_CPUE: " <<pred_cpue<<endl; 
    report <<"observed_CPUE: "  <<obs_cpue <<endl; 
    report <<"q_CPUE: "         <<q_cpue   <<endl; 
  }
  report << "Obs_P_fish_age"<<aa_D <<endl;
  for (i=1;i<=nyrs_fish_age;i++) report << yrs_fish_age(i)<<" "<<oac_fish(i) <<endl; report<<endl;
  report << "Pred_P_fish_age"<<aa_D <<endl;
  for (i=1;i<=nyrs_fish_age;i++) report << yrs_fish_age(i)<<" "<<eac_fish(i) <<endl; report<<endl;
  report << "yrs_fish_age nmulti_fish_age effn_fish_age sdnr_fish_age" <<endl;
  for (i=1;i<=nyrs_fish_age;i++) report <<yrs_fish_age(i)<<" "<<nmulti_fish_age(i)<<" "<<effn_fish_age(i)<<" "<<sdnr_fish_age(i)<<endl; report<<endl;
  report << "Obs_P_fish_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish_size;i++) report << yrs_fish_size(i)<<" "<<osc_fish(i) <<endl; report<<endl;
  report << "Pred_P_fish_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish_size;i++) report << yrs_fish_size(i)<<" "<<esc_fish(i) <<endl; report<<endl;
  report << "yrs_fish_size nmulti_fish_size effn_fish_size sdnr_fish_size" <<endl;
  for (i=1;i<=nyrs_fish_size;i++) report <<yrs_fish_size(i)<<" "<<nmulti_fish_size(i)<<" "<<effn_fish_size(i)<<" "<<sdnr_fish_size(i)<<endl; report<<endl;
  report << "Obs_P_srv1_age"<<aa_D <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i)<<" "<<oac_srv1(i) <<endl; report<<endl;
  report << "Pred_P_srv1_age"<<aa_D <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i)<<" "<<eac_srv1(i) <<endl; report<<endl;
  report << "yrs_srv1_age nmulti_srv1_age effn_srv1_age sdnr_srv1_age" <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) report <<yrs_srv1_age(i)<<" "<<nmulti_srv1_age(i)<<" "<<effn_srv1_age(i)<<" "<<sdnr_srv1_age(i)<<endl; report<<endl;
  report << "Obs_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<osc_srv1(i) <<endl; report<<endl;
  report << "Pred_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<esc_srv1(i) <<endl; report<<endl;
  report << "Obs_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<osc_srv2(i) <<endl; report<<endl;
  report << "Pred_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<esc_srv2(i) <<endl; report<<endl;
  report << "Bottom Trawl Survey Biomass " <<endl;
  report << "Year:     " << yrs_srv1  <<endl;
  report << "Predicted:   " << pred_srv1  <<endl;
  report << "Observed:   " << obs_srv1_biom  <<endl;
  report << "Observed_SE:   " << obs_srv1_se  <<endl<< endl;
  report << "Alternative Survey Biomass " <<endl;
  report << "Year:     " << yrs_srv2  <<endl;
  report << "Predicted:   " << pred_srv2  <<endl;
  report << "Observed:   " << obs_srv2_biom  <<endl;
  report << "Observed_SE:   " << obs_srv2_se  <<endl<< endl;
  report << "Year " << yy<< endl;
  report << "SpBiom "<< spawn_biom <<endl;
  report << "Tot_biom "<< tot_biom   <<endl;
  report << "Recruitment "<< pred_rec   <<endl;
  if (fishselopt==3 || fishselopt==4){
    report << "Fully_selected_F "<< rescaled_F <<endl<<endl;;
  }
  else{
    report << "Fully_selected_F "<< Fmort*max(fish_sel4) <<endl<<endl;;
  }
  report << "Age  "<<aa<< endl;
  report << "Weight "<< wt << endl;
  report << "Maturity "<<p_mature<<endl;
  if (fishselopt==3 || fishselopt==4){
    for (i=styr;i<=endyr;i++) report <<"Fishery_Selectivity_"<< i<<" "<<rescaled_fish_sel(i) <<endl;
  }
  else{
    for (i=styr;i<=endyr;i++) report <<"Fishery_Selectivity_"<< i<<" "<<fish_sel(i) <<endl;  
  }
  //report << "Fishery_Selectivity_1967-1976" << fish_sel1  <<endl;
  //report << "Fishery_Selectivity_1977-1995" << fish_sel2  <<endl;
  //report << "Fishery_Selectivity_1996-2006" << fish_sel3  <<endl;
  //report << "Fishery_Selectivity_2007-"<<endyr << fish_sel4  <<endl;
  report << "Bottom_Trawl_Survey_Selectivity " << srv1_sel / max(srv1_sel) <<endl;
  report << "Alternative_Survey_Selectivity " << srv2_sel / max(srv2_sel) <<endl<<endl;
  report << "F35 F40 F50 "<<endl;
  report <<  F35 << " "<< F40 <<" "<<  F50 <<endl;report<<endl;
  report << "~~~~~~~ Wts and Likelihoods ~~~~~~~"<< endl;report<<endl;
  report << wt_ssqcatch <<" "<<ssqcatch     <<" " ; report << "SSQ Catch Likelihood" << endl;
  report << wt_srv1 <<" "<<surv_like(1) <<" " ; report << "Bottom Trawl Survey Likelihood" << endl;
  report << wt_fish_age <<" "<<age_like(1)  <<" " ; report << "Fishery Age Composition Likelihood"  << endl;
  report << wt_srv1_age <<" "<<age_like(2)  <<" " ; report << "Bottom Trawl Survey Age Composition Likelihood" << endl;
  report << wt_fish_size <<" "<<age_like(3)  <<" " ; report << "Fishery Size Composition Likelihood" << endl;
  report << wt_cpue <<" "<<cpue_like    <<" " ; report << "Fishery CPUE Likelihood" << endl;
  report << wt_srv2 <<" "<<surv_like(2) <<" " ; report << "Alternative Survey Abundance Index Likelihood"   << endl;
  report << wt_srv1_size <<" "<<age_like(4)  <<" " ; report << "Bottom Trawl Survey Size Composition Likelihood" << endl;
  report << wt_srv2_size <<" "<<age_like(5)  <<" " ; report << "Alternative Survey Size Composition Likelihood" << endl;
  report << wt_rec_var <<" "<<rec_like     <<" " ; report << "Recruitment Deviations Likelihood" << endl;
  report << wt_fmort_reg <<" "<<F_mort_regularity<<" " ; report << "Fishing Mortality Deviations Penalty" << endl;
  report << 1 <<" "<<priors(1)  <<" " ; report << "Priors SigmaR" <<endl;
  report << 1 <<" "<<priors(2)  <<" " ; report << "Priors q Bottom Trawl Survey" <<endl;
  report << 1 <<" "<<priors(3)  <<" " ; report << "Priors M" <<endl;
  report << 1 <<" "<<priors(4)  <<" " ; report << "Priors q Alternative Survey" <<endl;
  report << " "<<obj_fun    <<" " ; report << "Objective Function" <<endl;
  report << " "<<Like       <<" " ; report << "Data Likelihood" <<endl;report<<endl;
  report<<"~~~~~~~ Projection outputs ~~~~~~~"<<endl;report<<endl;
  report << "N_at_age projected "<<endl<<N_proj<<endl<<" spawn_bio projected"<<endl<<spawn_biom_proj<<endl;
}

double model_parameters::sdnr(const dvar_vector& pred,const dvector& obs,double m)
{
  ofstream& evalout= *pad_evalout;
  RETURN_ARRAYS_INCREMENT();
  double sdnr;
  dvector pp = value(pred);
  int ntmp = -obs.indexmin()+obs.indexmax();
  sdnr = std_dev(elem_div(obs-pp,sqrt(elem_prod(pp,(1.-pp))/m)));
  RETURN_ARRAYS_DECREMENT();
  return sdnr;
}

void model_parameters::write_proj(void)
{
  ofstream& evalout= *pad_evalout;
 ofstream newproj("proj.dat");
 newproj <<"#Species name here:"<<endl;
 newproj <<model_name+"_"+data_file<<endl;
 newproj <<"#SSL Species?"<<endl;
 newproj <<"0"<<endl;
 newproj <<"#Constant buffer of Dorn?"<<endl;
 newproj <<"0"<<endl;
 newproj <<"#Number of fisheries?"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#Number of sexes?"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#5year_Average_F(endyr-4,endyr_as_estimated_by_ADmodel)"<<endl;
 newproj << mean(Fmort(endyr-4,endyr))<<endl;
 newproj <<"#_Author_F_as_fraction_F_40%"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#ABC SPR" <<endl;
 newproj <<"0.4"<<endl;
 newproj <<"#MSY SPR" <<endl;
 newproj <<"0.35"<<endl;
 newproj <<"#_Spawn_month"<<endl;
 newproj << spawn_fract*12+1<<endl;
 newproj <<"#_Number_of_ages"<<endl;
 newproj <<nages_M<<endl;
 newproj <<"#_F_ratio(must_sum_to_one_only_one_fishery)"<<endl;
 newproj <<"1"<<endl;
 for (j=1;j<=nages_M;j++) natmortv = natmort; 
 newproj <<"#_Natural_Mortality" << aa << endl;
 newproj <<natmortv<<endl;
 newproj <<"#_Maturity_divided_by_2(projection_program_uses_to_get_female_spawning_biomass_if_divide_by_2"<<aa<<endl<<p_mature<< endl;
 newproj <<"#_Wt_at_age_spawners"<<aa<<endl<<wt<< endl;
 newproj <<"#_Wt_at_age_fishery" <<aa<<endl<<wt<< endl;
 newproj <<"#_Selectivity_fishery_scaled_to_max_at_one"<<aa<<endl<<fish_sel4/max(fish_sel4)<< endl;
 newproj <<"#_Numbers_at_age_end_year"<<aa<<endl<<natage(endyr)<< endl;
 newproj <<"#_N_recruitment_years"<<endl<<endyr-recage-1979+1<< endl;
 newproj <<"#_Recruitment_start_at_1977_yearclass=1979_for_age_2_recruits"<<yy(1979,endyr-recage)<<endl<<pred_rec(1979,endyr-recage)<< endl;
 newproj <<"#_Spawners per recruitment (starting at 1977)"<<endl<<spawn_biom(1977,endyr-recage)/1000<< endl;
 newproj.close();
}

void model_parameters::final_calcs()
{
  // R_Report(fish_sel);
  R_report<<"#Selectivity"<<endl; 
  for (i=styr;i<=endyr;i++) 
    R_report<<i<<" "<<fish_sel(i)<<endl;
    // sdreport_vector      spawn_biom(styr,endyr);                 // " " for spawning biomass vector
  R_report<<"#SSB"<<endl; 
  for (i=styr;i<=endyr;i++) 
  {
    // sdreport_vector      spawn_biom(styr,endyr);                 // " " for spawning biomass vector
    double lb=value(spawn_biom(i)/exp(2.*sqrt(log(1+square(spawn_biom.sd(i))/square(spawn_biom(i))))));
    double ub=value(spawn_biom(i)*exp(2.*sqrt(log(1+square(spawn_biom.sd(i))/square(spawn_biom(i))))));
    R_report<<i<<" "<<spawn_biom(i)<<" "<<spawn_biom.sd(i)<<" "<<lb<<" "<<ub<<endl;
  }
  R_report<<"#R"<<endl; 
  for (i=styr;i<=endyr;i++) 
  {
    // sdreport_vector      pred_rec(styr,endyr);                   // " " for predicted recruitments
    double lb=value(pred_rec(i)/exp(2.*sqrt(log(1+square(pred_rec.sd(i))/square(pred_rec(i))))));
    double ub=value(pred_rec(i)*exp(2.*sqrt(log(1+square(pred_rec.sd(i))/square(pred_rec(i))))));
    R_report<<i<<" "<<pred_rec(i)<<" "<<pred_rec.sd(i)<<" "<<lb<<" "<<ub<<endl;
  }
}

void model_parameters::set_runtime(void)
{
  dvector temp("{1.e-4 1.e-4 1.e-4 1.e-7 1.e-7 1.e-7}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{1000, 1000, 1000, 10000, 20000, 20000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_evalout;
  pad_evalout = NULL;
}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
  arrmblsize = 1000000;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
