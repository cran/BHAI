


#' Simulate PPS data
#' 
#' @title Simulate PPS data
#' 
#' @param pps_data The PPS object containing the data. Parameters for simulations are extracted from this data.
#' @param num_survey_patients Numeric vector indicating sample sizes for simulations.
#'
#' @seealso \code{\linkS4class{PPS}}
#' 
#' @return A simulated PPS object.
#' 
#' @examples 
#' 
#' # Specify the number of survey patients
#' sim_survey_patients = 10000
#' # Subsample data sets from european PPS
#' sim_pps = sample.pps(eu_pps, num_survey_patients = sim_survey_patients)
#' 
#' @export
sample.pps = function(pps_data, num_survey_patients) {
  

  pps_list = list()
  
  for(curr_num_survey_patients in num_survey_patients) {
    
    # get (fixed) parameters from PPS data
    prev = pps_data@num_hai_patients/pps_data@num_survey_patients
    length_of_stay = pps_data@length_of_stay
    
    # Get number of patients for current sample size
    num_hai_patients = round(curr_num_survey_patients*prev)
    
     # Sample length of infection
    loi_pps = lapply(names(pps_data@loi_pps), 
                                 function(x) sample(pps_data@loi_pps[[x]], num_hai_patients[x], replace=TRUE))
    names(loi_pps) = names(pps_data@loi_pps)
 
    # Calculate observed age-gender distribution from PPS data 
    prob_hai_patients_by_stratum = lapply(pps_data@num_hai_patients_by_stratum, function(x) x/sum(x))
    # Distribute sampled cases into age gender categories
    num_hai_patients_by_stratum = lapply(names(prob_hai_patients_by_stratum),
                                         function(x) matrix(table(factor(sample(1:38, num_hai_patients[x], 
                                                                                prob=prob_hai_patients_by_stratum[[x]], replace=TRUE), 
                                                                         levels=1:38)), ncol=2))
    names(num_hai_patients_by_stratum) = names(prob_hai_patients_by_stratum)
 
    # Create a sample of McCabe scores  
    mccabe_scores_distr = list()
    for(curr_inf in names(pps_data@mccabe_scores_distr)) {
      # Starat with > 0 cases are sampled
      sample_this = as.vector(num_hai_patients_by_stratum[[curr_inf]])
      # McCabe score distribution for each stratum
      all_probs = apply(do.call("rbind", 
                                lapply(pps_data@mccabe_scores_distr[[curr_inf]], as.vector)), 
                        2, function(x) x/sum(x))
      
      # Create McCabe score distribution for current sample
      new_mccabe = list()
      for(i in as.character(1:3)) {
        new_mccabe[[i]] = matrix(0, nrow=19, ncol=2)
      }
      for(i in 1:ncol(all_probs)) {
        # Starta > 0 are sampled 
        if(sample_this[i] > 0) {
          tab = table(factor(sample(1:3, sample_this[i], prob=all_probs[,i], replace=TRUE), levels=1:3))
          for(j in as.character(1:3)) {
            new_mccabe[[j]][i] = tab[j]
          }
        }
      }
      names(new_mccabe) = names(pps_data@mccabe_scores_distr[[curr_inf]])
      mccabe_scores_distr[[curr_inf]] = new_mccabe
    }
    
    
    # Denominator for strata
    num_survey_patients_by_stratum = pps_data@num_survey_patients_by_stratum
    prob_survey_patients_by_stratum = lapply(num_survey_patients_by_stratum, function(x) x/sum(x))
    cases_mccabe = lapply(names(prob_survey_patients_by_stratum), function(x) Reduce("+", lapply(mccabe_scores_distr, function(y) y[[x]])))
    names(cases_mccabe) = names(prob_survey_patients_by_stratum)
    # Distribute sampled cases into age gender categories
    num_survey_patients_by_stratum = table(factor(sample(1:length(unlist(prob_survey_patients_by_stratum)), curr_num_survey_patients-sum(unlist(cases_mccabe)), 
                                                         prob=unlist(prob_survey_patients_by_stratum), 
                                                         replace=TRUE), levels=1:length(unlist(prob_survey_patients_by_stratum))))
    num_survey_patients_by_stratum = list(matrix(num_survey_patients_by_stratum[as.character(1:38)], ncol=2),
                                          matrix(num_survey_patients_by_stratum[as.character(39:76)], ncol=2),
                                          matrix(num_survey_patients_by_stratum[as.character(77:114)], ncol=2))
    names(num_survey_patients_by_stratum) = names(prob_survey_patients_by_stratum)
    for(i in names(num_survey_patients_by_stratum)) {
      num_survey_patients_by_stratum[[i]] = num_survey_patients_by_stratum[[i]]+cases_mccabe[[i]]
    }
    
    curr_pps_data = pps_data
    curr_pps_data@num_hai_patients = num_hai_patients
    curr_pps_data@num_survey_patients = curr_num_survey_patients
    curr_pps_data@loi_pps = loi_pps
    curr_pps_data@num_hai_patients_by_stratum = num_hai_patients_by_stratum
    curr_pps_data@mccabe_scores_distr = mccabe_scores_distr
    curr_pps_data@num_hai_patients_by_stratum_prior = curr_pps_data@num_hai_patients_by_stratum
    curr_pps_data@mccabe_by_stratum_prior = pps_data@mccabe_scores_distr
    curr_pps_data@num_survey_patients_by_stratum = num_survey_patients_by_stratum
    
    pps_list[[as.character(curr_num_survey_patients)]] = curr_pps_data
  }
  pps_list
}
