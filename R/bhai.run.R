


#'  
#' Estimate the burden of healthcare-associated infections (bhai) from point prevalence surveys 
#' 
#' @title Burden of healthcare-associated infections
#' 
#' @param num_hai_patients Named numeric containing patients having healthcare-associated infections.
#' @param num_survey_patients Number of patients in point prevalence survey.
#' @param length_of_stay Length of stay of all patients in hospitals. This is need for the prevalence to incidence conversion with the Rhame-Sudderth formula.
#' @param loi_pps A list containing length of infections from all patients in the PPS. The length of infection of all healthcare-associated infections. In PPS this is usually approximated as the time from infection onset until the date of the survey.
#' @param hospital_discharges The number of hospital discharges.
#' @param num_hai_patients_by_stratum A list containing for each infection the number of patients in each age and gender stratum.
#' @param num_hai_patients_by_stratum_prior The prior weight (counts) for each infection, age and gender stratum. This is used for smooting the age and gender distribution when small numbers are observed.
#' @param mccabe_scores_distr The observed McCabe scores (counts) for each infection, age and gender stratum from the PPS.
#' @param mccabe_by_stratum_prior The prior weight (counts) for each infection, McCabe score, age and gender stratum. This is used for smooting the age and gender distribution when small numbers are observed.
#' @param mccabe_life_exp Named list containing remaining life expectancies for each McCabe score (NONFATAL, ULTFATAL, RAPFATAL).
#' @param pop.sampling Indicating whether the sampling of the outcome tree parameters (e.g. probability of death, disability weights) should be carried out in each stratum separately or for the whole population.
#' @param sample_distr The distribution for the prevalence sampling. Default is 'rbetamix' which uses the mid-P-Pearson-Clopper interval. Alternatives are 'rpert' - which uses the PERT distribution - and 'bcode' which uses the sampling procedures as described in Cassini et al. (2016) <doi:https://doi.org/10.1371/journal.pmed.1002150>.
#' @param infections Infections for which the burden should be calculated. Default is names(num_hai_patients).
#' @param num_survey_patients_by_stratum Number of survey patients stratified by infection, age and gender. If this parameter is provided the methodology described in Cassini et al. (2016) <doi:https://doi.org/10.1371/journal.pmed.1002150> is applied.
#' @param estimate_loi_fun Function to use for the estimation of length of infection. Default is runif_bootstrap_rear_gren.
#' 
#' @usage bhai.run(num_hai_patients, num_survey_patients, length_of_stay, loi_pps,  hospital_discharges, num_hai_patients_by_stratum, num_hai_patients_by_stratum_prior=NULL, mccabe_scores_distr, mccabe_by_stratum_prior=NULL, mccabe_life_exp, pop.sampling=TRUE, nsim=1000, sample_distr="rbetamix", infections=c("UTI", "SSI", "CDI", "PN", "BSI"), num_survey_patients_by_stratum=NULL, estimate_loi_fun=runif_bootstrap_rear_gren)
#' 
#' @return A list containing estimates for all input infections.
#' 
#' @keywords internal
#' @noRd
bhai.run = function(num_hai_patients, num_survey_patients, 
                    length_of_stay, loi_pps, 
                    hospital_discharges,
                    num_hai_patients_by_stratum, num_hai_patients_by_stratum_prior=NULL,  
                    mccabe_scores_distr, mccabe_by_stratum_prior=NULL, mccabe_life_exp,
                    pop.sampling=TRUE, nsim=1000, sample_distr="rbetamix", 
                    infections=names(num_hai_patients), 
                    num_survey_patients_by_stratum=NULL,
                    estimate_loi_fun=bootstrap_mean_gren,
                    stratified_sampling=FALSE) {
  
  ncases_by_stratum = num_hai_patients_by_stratum
 
 
  age_prior = num_hai_patients_by_stratum_prior
  if(length(age_prior) > 0) {
  #  pseudocount = Reduce("+",age_prior)
  #  pseudocount = pseudocount/sum(pseudocount)+1e-3
    for(i in 1:length(age_prior)) {
      age_prior[[i]] = age_prior[[i]]+max(c(1e-3, 1e-3*sum(age_prior[[i]])))
    }
  }
   
  
  nhai = num_hai_patients
  npatients = num_survey_patients
  
  if(length(mccabe_scores_distr) > 0) {
    for(n in names(mccabe_scores_distr)) {
      mccabe_overall = sapply(mccabe_scores_distr[[n]], sum)
      mccabe_overall = mccabe_overall/sum(mccabe_overall)
      for(i in 1:length(mccabe_scores_distr[[n]])) {
        # Add a pseudocount to make sure that no cases are lost due to empty strata
        # Empty strata will have McCabe distribution as average for infection
        pseudocount = 1e-3*mccabe_overall[i]+1e-3
        mccabe_scores_distr[[n]][[i]] = mccabe_scores_distr[[n]][[i]]+pseudocount
       }
    }
  }
  
  
  mccabe_prior = mccabe_by_stratum_prior
  if(length(mccabe_prior) > 0) {
    for(n in names(mccabe_prior)) {
      for(i in 1:length(mccabe_prior[[n]])) {
        # Add a pseudocount to make sure that no cases are lost due to empty strata
        #if(any(mccabe_prior[[n]][[i]] == 0)) {
          mccabe_prior[[n]][[i]] = mccabe_prior[[n]][[i]]+
            max(c(1e-3, 1e-3*sum(mccabe_prior[[n]][[i]])))
        #}
      }
    }  
  }

  # The following have be > 0. Otherwise cases might be 'lost' when distributed in starta and McCabe categories
  if(any(unlist(mccabe_prior) == 0)) 
    warning("Empty strata in McCabe prior! Might lead to flawed estimation!")
  if(any(unlist(mccabe_scores_distr) == 0)) 
    warning("Empty strata in McCabe distribution! Might lead to flawed estimation!")
  if(any(unlist(age_prior) == 0)) 
    warning("Empty strata in Age-gender prior! Might lead to flawed estimation!")
        
  discharges = hospital_discharges
  la = length_of_stay
  loi = loi_pps
  pop.level = pop.sampling
  p_age = NULL
  if(length(num_survey_patients_by_stratum) > 0 & stratified_sampling) {
    p_age = num_survey_patients_by_stratum#/sum(num_survey_patients_by_stratum)
  }
  else if(length(num_survey_patients_by_stratum) == 0 & stratified_sampling) {
    warning("Strified sampling cannot be performed: stratified_sampling==TRUE but num_survey_patients_by_stratum is empty in PPS object!\n")
  }
  
  allsim = list()
  for(n in infections) {
    # print(n)
    allsim[[n]] = list()
    currloi = NULL
    if(class(loi) == "list" & length(loi) == 2) {
      currloi = sort(unlist(sapply(loi, function(x) x[n])))
    }
    else {
      currloi = loi[n]
    }
    #(print(currloi))
    for(k in 1:nsim) {
      
      allsim[[n]][[k]] = switch(n,
                                HAP = sim.pn(nhai[n], npatients, la, currloi, 
                                            discharges, mccabe_life_exp,  
                                            mccabe_scores_distr[[n]], ncases_by_stratum[[n]], age_prior=age_prior[[n]], p_age=p_age, pop.level=pop.level, sample_distr=sample_distr, mccabe_prior[[n]], estimate_loi=estimate_loi_fun),
                                BSI = sim.bsi(nhai[n], npatients, la, currloi, 
                                              discharges, mccabe_life_exp,  
                                              mccabe_scores_distr[[n]], ncases_by_stratum[[n]], age_prior=age_prior[[n]], p_age=p_age, pop.level=pop.level, sample_distr=sample_distr, mccabe_prior[[n]], estimate_loi=estimate_loi_fun),
                                SSI = sim.ssi(nhai[n], npatients, la, currloi, 
                                              discharges, mccabe_life_exp,  
                                              mccabe_scores_distr[[n]], ncases_by_stratum[[n]], age_prior=age_prior[[n]], p_age=p_age, pop.level=pop.level, sample_distr=sample_distr, mccabe_prior[[n]], estimate_loi=estimate_loi_fun),
                                CDI = sim.cdi(nhai[n], npatients, la, currloi, 
                                              discharges, mccabe_life_exp,  
                                              mccabe_scores_distr[[n]], ncases_by_stratum[[n]], age_prior=age_prior[[n]], p_age=p_age, pop.level=pop.level, sample_distr=sample_distr, mccabe_prior[[n]], estimate_loi=estimate_loi_fun),
                                UTI = sim.uti(nhai[n], npatients, la, currloi, 
                                              discharges, mccabe_life_exp,  
                                              mccabe_scores_distr[[n]], ncases_by_stratum[[n]], age_prior=age_prior[[n]], p_age=p_age, pop.level=pop.level, sample_distr=sample_distr, mccabe_prior[[n]], estimate_loi=estimate_loi_fun))
    }
  }
  allsim
}
