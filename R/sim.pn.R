#' Performs simulation of the outcome tree of heatlhcare-associated pneumonia (HAP) for DALY calculation
#' 
#' @param nhai Named numeric containing patients having healthcare-associated infections.
#' @param npatients Number of patients in point prevalence survey.
#' @param la Length of stay of all patients in hospitals. This is need for the prevalence to incidence conversion with the Rhame-Sudderth formula.
#' @param loi A list containing length of infections from all patients in the PPS. The length of infection of all healthcare-associated infections. In PPS this is usually approximated as the time from infection onset until the date of the survey.
#' @param discharges The number of hospital discharges.
#' @param ncases_by_stratum A list containing for each infection the number of patients in each age and gender stratum.
#' @param age_prior The prior weight (counts) for each infection, age and gender stratum. This is used for smooting the age and gender distribution when small numbers are observed.
#' @param mccabe_scores_distr The observed McCabe scores (counts) for each infection, age and gender stratum from the PPS.
#' @param mccabe_prior The prior weight (counts) for each infection, McCabe score, age and gender stratum. This is used for smooting the age and gender distribution when small numbers are observed.
#' @param mccabe_life_exp Named list containing remaining life expectancies for each McCabe score (NONFATAL, ULTFATAL, RAPFATAL).
#' @param pop.level Indicating whether the sampling of the outcome tree parameters (e.g. probability of death, disability weights) should be carried out in each stratum separately or for the whole population.
#' @param sample_distr The distribution for the prevalence sampling. Default is 'rbetamix' which uses the mid-P-Pearson-Clopper interval. Alternatives are 'rpert' - which uses the PERT distribution - and 'bcode' which uses the sampling procedures as described in Cassini et al. (2016) <doi:https://doi.org/10.1371/journal.pmed.1002150>.
#' @param infections Infections for which the burden should be calculated. Default is names(num_hai_patients).
#' @param p_age Number of survey patients stratified by infection, age and gender. If this parameter is provided the methodology described in Cassini et al. (2016) <doi:https://doi.org/10.1371/journal.pmed.1002150> is applied.
#' @param estimate_loi_fun Function to use for the estimation of length of infection. Default is runif_bootstrap_rear_gren.
#' 
#' @keywords internal
#' @seealso 
#' Cassini et al. (2016) <doi:https://doi.org/10.1371/journal.pmed.1002150>
#' Colzani et al. (2017) <doi:https://doi.org/10.1371/journal.pone.0170662>
#'  
#' @keywords internal
#' @noRd
sim.pn = function(nhai, npatients, la, loi, discharges, mccabe_life_exp, mccabe_scores_distr, ncases_by_stratum, sample_distr="rbetamix", age_prior, p_age=NULL, pop.level, mccabe_prior, estimate_loi) {
  
  
  ncases_hai = sample.ncases(nhai, npatients, la, loi, discharges, mccabe_life_exp, mccabe_scores_distr, ncases_by_stratum, sample_distr, age_prior, p_age=p_age, mccabe_prior, estimate_loi=estimate_loi)
  
  output = list()
  output$ncases_hai = ncases_hai # sum(unlist(ncases_hai))
  
  prob_death = sample.rpert(19*2, 0.001, 0.09, 0.035, pop.level, names(ncases_hai))
  ndeath_hai = list()
  for(i in names(ncases_hai)) {
    ndeath_hai[[i]] = ncases_hai[[i]]*prob_death[[i]]
  }
  names(ndeath_hai) = names(ncases_hai)
  output$ndeath_hai = ndeath_hai # sum(unlist(ndeath_hai))
  
  yll_hai = ndeath_hai
  for(i in names(mccabe_life_exp)) {
    yll_hai[[i]] = yll_hai[[i]]*mccabe_life_exp[[i]]
  }
  output$yll_hai = yll_hai # sum(unlist(yll_hai))
  
  weight_symptomatic_infection = sample.rpert(19*2, 0.104, 0.152, 0.125, pop.level, names(ncases_hai))
  duration_symptomatic_infection = sample.runif(19*2, 0.019, 0.031, pop.level, names(ncases_hai))
  daly_symptomatic_infection = list()
  for(i in names(ncases_hai)) {
    daly_symptomatic_infection[[i]] = ncases_hai[[i]]*weight_symptomatic_infection[[i]]*duration_symptomatic_infection[[i]]
  }
  output$daly_symptomatic_infection = daly_symptomatic_infection # sum(unlist(daly_symptomatic_infection))
  
  weight_severe_sepsis_shock = sample.rpert(19*2, 0.579, 0.727, 0.655, pop.level, names(ncases_hai))
  duration_severe_sepsis_shock = sample.runif(19*2, 0.027, 0.036, pop.level, names(ncases_hai))
  severe_sepsis_shock_trans = list()
  severe_sepsis_shock = list()
  for(i in names(ncases_hai)) {
    severe_sepsis_shock_trans[[i]] = matrix(rep(0.39,19*2), ncol=2)
    severe_sepsis_shock[[i]] = ncases_hai[[i]]*severe_sepsis_shock_trans[[i]]* # Transition probability
      weight_severe_sepsis_shock[[i]]*duration_severe_sepsis_shock[[i]] # DALYs
  }
  output$severe_sepsis_shock = severe_sepsis_shock # sum(unlist(severe_sepsis_shock))
  #
  prob_post_traumatic_disorder = sample.runif(19*2, 0.13, 0.21, pop.level, names(ncases_hai))
  weight_post_traumatic_disorder = sample.rpert(19*2, 0.07, 0.108, 0.088, pop.level, names(ncases_hai))
  post_traumatic_disorder = list()
  for(i in names(ncases_hai)) {
    post_traumatic_disorder_trans = severe_sepsis_shock_trans[[i]]*prob_post_traumatic_disorder[[i]]
    post_traumatic_disorder[[i]] = ncases_hai[[i]]*post_traumatic_disorder_trans* # Transition probability
      weight_post_traumatic_disorder[[i]]
    post_traumatic_disorder[[i]] = post_traumatic_disorder[[i]]*mccabe_life_exp[[i]]
  }
  output$post_traumatic_disorder = post_traumatic_disorder # sum(unlist(post_traumatic_disorder))
  #
  prob_cognitive_impairment = sample.runif(19*2, 0.11, 0.47, pop.level, names(ncases_hai))
  weight_cognitive_impairment = sample.rpert(19*2, 0.026, 0.064, 0.043, pop.level, names(ncases_hai))
  cognitive_impairment = list()
  for(i in names(ncases_hai)) {
    cognitive_impairment_trans = severe_sepsis_shock_trans[[i]]*prob_cognitive_impairment[[i]]
    cognitive_impairment[[i]] = ncases_hai[[i]]*cognitive_impairment_trans* # Transition probability
      weight_cognitive_impairment[[i]]
    cognitive_impairment[[i]] = cognitive_impairment[[i]]*mccabe_life_exp[[i]]
    
  }
  output$cognitive_impairment = cognitive_impairment # sum(unlist(cognitive_impairment))
  #
  weight_physical_impairment = sample.runif(19*2, 0.011, 0.053, pop.level, names(ncases_hai))
  physical_impairment = list()
  for(i in names(ncases_hai)) {
    physical_impairment_trans = severe_sepsis_shock_trans[[i]]*matrix(rep(1,19*2), ncol=2)
    physical_impairment[[i]] = ncases_hai[[i]]*physical_impairment_trans* # Transition probability
      weight_physical_impairment[[i]]
    physical_impairment[[i]] = physical_impairment[[i]]*mccabe_life_exp[[i]]
  }
  output$physical_impairment = physical_impairment # sum(unlist(physical_impairment))
  #
  prob_renal_failure = sample.runif(19*2, 0.009, 0.013, pop.level, names(ncases_hai))
  weight_renal_failure = sample.runif(19*2, 0.03, 0.487, pop.level, names(ncases_hai))
  renal_failure = list()
  for(i in names(ncases_hai)) {
    renal_failure_trans = severe_sepsis_shock_trans[[i]]*prob_renal_failure[[i]]
    renal_failure[[i]] = ncases_hai[[i]]*renal_failure_trans* # Transition probability
      weight_renal_failure[[i]]
    renal_failure[[i]] = renal_failure[[i]]*mccabe_life_exp[[i]]
  }
  output$renal_failure = renal_failure # sum(unlist(renal_failure))
  output
}

