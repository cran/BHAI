#' Performs simulation of the outcome tree of surgical site infection (SSI) for DALY calculation
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
sim.ssi = function(nhai, npatients, la, loi, discharges, mccabe_life_exp, mccabe_scores_distr, ncases_by_stratum, sample_distr="rbetamix", age_prior, p_age=NULL, pop.level, mccabe_prior, estimate_loi) {
  
  
  ncases_hai = sample.ncases(nhai, npatients, la, loi, discharges, mccabe_life_exp, mccabe_scores_distr, ncases_by_stratum, sample_distr, age_prior, p_age=p_age, mccabe_prior, estimate_loi=estimate_loi)
  
  output = list()
  output$ncases_hai = ncases_hai # sum(unlist(ncases_hai)) # 
  
  prob_death = matrix(rep(c(rep(9e-3,14), rep(0.036,5)),2), ncol=2)
  ndeath_hai = list()
  for(i in names(ncases_hai)) {
    ndeath_hai[[i]] = ncases_hai[[i]]*prob_death
  }
  names(ndeath_hai) = names(ncases_hai)
  output$ndeath_hai = ndeath_hai # sum(unlist(ndeath_hai)) # 
  
  yll_hai = ndeath_hai
  for(i in names(mccabe_life_exp)) {
    yll_hai[[i]] = yll_hai[[i]]*mccabe_life_exp[[i]]
  }
  output$yll_hai = yll_hai # sum(unlist(yll_hai)) # 
  # 
  weight_symptomatic_infection = sample.rpert(19*2, 0.039, 0.06, 0.051, pop.level, names(ncases_hai))
  duration_symptomatic_infection = 0.096
  daly_symptomatic_infection = list()
  for(i in names(ncases_hai)) {
    
    daly_symptomatic_infection[[i]] = ncases_hai[[i]]*weight_symptomatic_infection[[i]]*duration_symptomatic_infection
  }
  output$daly_symptomatic_infection = daly_symptomatic_infection # sum(unlist(daly_symptomatic_infection)) # 
  
  output
}

