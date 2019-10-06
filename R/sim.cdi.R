#' Performs simulation of the outcome tree of C. difficile infection (CDI) for DALY calculation
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
sim.cdi = function(nhai, npatients, la, loi, discharges, mccabe_life_exp, mccabe_scores_distr, ncases_by_stratum, sample_distr="rbetamix", age_prior, p_age=NULL, pop.level, mccabe_prior, estimate_loi) {
  
  
  ncases_hai = sample.ncases(nhai, npatients, la, loi, discharges, mccabe_life_exp, mccabe_scores_distr, ncases_by_stratum, sample_distr, age_prior, p_age=p_age, mccabe_prior, estimate_loi=estimate_loi)
  
  output = list()
  output$ncases_hai = ncases_hai # sum(unlist(ncases_hai))
  
  
  prob_death = sample.runif(19*2, 0, 0.11, pop.level, names(ncases_hai))
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
  #
  prop_cases_uncomplicated = sample.runif(19*2, 0.85, 0.985, pop.level, names(ncases_hai))
  weight_symptomatic_infection_uncomplicated = sample.runif(19*2, 0.073, 0.149, pop.level, names(ncases_hai))
  duration_symptomatic_infection_uncomplicated = sample.runif(19*2, 0, 0.0219, pop.level, names(ncases_hai))
  weight_symptomatic_infection_complicated = sample.rpert(19*2, 0.202, 0.285, 0.239, pop.level, names(ncases_hai))
  duration_symptomatic_infection_complicated = sample.runif(19*2, 0, 0.0219, pop.level, names(ncases_hai))
  daly_symptomatic_infection_uncomplicated = list()
  daly_symptomatic_infection_complicated = list()
  for(i in names(ncases_hai)) {
    prop_cases_complicated = 1-prop_cases_uncomplicated[[i]]
    
    daly_symptomatic_infection_uncomplicated[[i]] =  ncases_hai[[i]]*prop_cases_uncomplicated[[i]]*
      weight_symptomatic_infection_uncomplicated[[i]]*duration_symptomatic_infection_uncomplicated[[i]]
    daly_symptomatic_infection_complicated[[i]] =  ncases_hai[[i]]*prop_cases_complicated*
      weight_symptomatic_infection_complicated[[i]]*duration_symptomatic_infection_complicated[[i]]
  }
  output$daly_symptomatic_infection_uncomplicated = daly_symptomatic_infection_uncomplicated # sum(unlist(daly_symptomatic_infection_uncomplicated))
  output$daly_symptomatic_infection_complicated = daly_symptomatic_infection_complicated # sum(unlist(daly_symptomatic_infection_complicated))
  #
  prob_post_colectomy_state = sample.runif(19*2, 0.002, 0.038, pop.level, names(ncases_hai))
  weight_post_colectomy_state = sample.rpert(19*2, 0.104, 0.155, 0.125, pop.level, names(ncases_hai))
  post_colectomy_state = list()
  for(i in names(ncases_hai)) {
    post_colectomy_state_trans = prob_post_colectomy_state[[i]]
    post_colectomy_state[[i]] = ncases_hai[[i]]*post_colectomy_state_trans* # Transition probability
      weight_post_colectomy_state[[i]]
    post_colectomy_state[[i]] = post_colectomy_state[[i]]*mccabe_life_exp[[i]]
  }
  output$post_colectomy_state = post_colectomy_state # sum(unlist(post_colectomy_state))
  
  output
}
