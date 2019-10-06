#' BHAI:
#'
#'
#' The \strong{BHAI} package
#'
#' @section BHAI functions:
#' \code{\link{bhai}}:
#'
#' @docType package
#' @name BHAI
#' @importFrom prevtoinc "calculate_I_smooth"
#' @importFrom MCMCpack "MCmultinomdirichlet"
#' @importFrom plotrix "draw.circle"
#' @importFrom methods "new"
#' @import stats
#' @import graphics
#' @import grDevices
NULL


#' Aggregated data of the ECDC PPS 2010-2011.
#'
#' @docType data
#' @keywords datasets
#' @name eu_pps
#' @usage data(eu_pps_2011)
#' @format A PPS object.
NULL

#' Aggregated data of the german PPS 2010-2011 (convenience sample).
#'
#' @docType data
#' @keywords datasets
#' @name german_pps_conv
#' @usage data(german_pps_2011_conv)
#' @format A PPS object.
NULL


#' Hospital discharges in Germany (2011)
#'
#' @docType data
#' @keywords datasets
#' @name hospital_discharges
#' @usage data(german_pps_2011_repr)
#' @format A PPS object.
NULL

#' Average length of stay of survey patients in german PPS 2011 (representative sample)
#'
#' @docType data
#' @keywords datasets
#' @name length_of_stay
#' @usage data(german_pps_2011_repr)
#' @format A PPS object.
NULL


#' A list containing length of infections from all patients in the german PPS 2011 representative sample.
#' 
#' @docType data
#' @keywords datasets
#' @name loi_pps
#' @usage data(german_pps_2011_repr)
#' @format A PPS object.
NULL


#' Named list containing remaining life expectancies for each McCabe score (NONFATAL, ULTFATAL, RAPFATAL).
#' 
#' @docType data
#' @keywords datasets
#' @name mccabe_life_exp
#' @usage data(german_pps_2011_repr)
#' @format A PPS object.
NULL


#' The observed McCabe scores (counts) for each infection, age and gender stratum from the ECDC PPS 2011-2012.
#' 
#' @docType data
#' @keywords datasets
#' @name mccabe_scores_distr
#' @usage data(german_pps_2011_repr)
#' @format A PPS object.
NULL


#' Number of cases for each infection in the german PPS 2011 (representative sample)
#' 
#' @docType data
#' @keywords datasets
#' @name num_hai_patients
#' @usage data(german_pps_2011_repr)
#' @format A PPS object.
NULL

#' Stratified number of cases for each infection in the german PPS 2011 (representative sample)
#' 
#' @docType data
#' @keywords datasets
#' @name num_hai_patients_by_stratum
#' @usage data(german_pps_2011_repr)
#' @format A PPS object.
NULL

#' Stratified number of cases for each infection in the german PPS 2011 (convenience sample). This distribution is used as a Prior for the representative sample.
#' 
#' @docType data
#' @keywords datasets
#' @name num_hai_patients_by_stratum_prior
#' @usage data(german_pps_2011_repr)
#' @format A PPS object.
NULL

#' Number of survey patients in the german PPS 2011 (representative sample).
#' 
#' @docType data
#' @keywords datasets
#' @name num_survey_patients
#' @usage data(german_pps_2011_repr)
#' @format A PPS object.
NULL

#' Population size of Germany in 2011.
#' 
#' @docType data
#' @keywords datasets
#' @name population
#' @usage data(german_pps_2011_repr)
#' @format A PPS object.
NULL


#' Simulated/subsampled data sets from european PPS
#' 
#' @docType data
#' @keywords datasets
#' @name sim_pps
#' @usage data(simulations)
#' @format A PPS object.
NULL

#' BHAI with stratified sampling was applied to simulated/subsampled data sets from european PPS
#' 
#' @docType data
#' @keywords datasets
#' @name sim_pps_stratified
#' @usage data(simulations)
#' @format A PPS object.
NULL

#' BHAI with default options was applied to simulated/subsampled data sets from european PPS
#' 
#' @docType data
#' @keywords datasets
#' @name sim_pps_bhai
#' @usage data(simulations)
#' @format A PPS object.
NULL

#' BHAI with prior was applied to simulated/subsampled data sets from european PPS
#' 
#' @docType data
#' @keywords datasets
#' @name sim_pps_bhai_prior
#' @usage data(simulations)
#' @format A PPS object.
NULL


