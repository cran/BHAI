

#'
#' This class is a generic container for PPS data sets.
#' 
#' @slot infections Character vector storing names of infections in PPS
#' @slot num_hai_patients Named numeric containing patients having healthcare-associated infections.
#' @slot num_survey_patients Number of patients in point prevalence survey.
#' @slot length_of_stay Length of stay of all patients in hospitals. This is need for the prevalence to incidence conversion with the Rhame-Sudderth formula.
#' @slot loi_pps A list containing length of infections from all patients in the PPS. In PPS this is usually calculated as the time from infection onset until the date of the survey.
#' @slot hospital_discharges The number of hospital discharges.
#' @slot num_hai_patients_by_stratum A list containing for each infection the number of patients in each age and gender stratum.
#' @slot num_hai_patients_by_stratum_prior The prior weight (counts) for each infection, age and gender stratum. This is used for smooting the age and gender distribution when small numbers are observed.
#' @slot mccabe_scores_distr The observed McCabe scores (counts) for each infection, age and gender stratum from the PPS.
#' @slot mccabe_by_stratum_prior The prior weight (counts) for each infection, McCabe score, age and gender stratum. This is used for smooting the age and gender distribution when small numbers are observed.
#' @slot mccabe_life_exp Named list containing remaining life expectancies for each McCabe score (NONFATAL, ULTFATAL, RAPFATAL).
#' @slot num_survey_patients_by_stratum Number of survey patients stratified by infection, age and gender. If this parameter is provided the methodology described in Cassini et al. (2016) <doi:https://doi.org/10.1371/journal.pmed.1002150> is applied.
#' @slot population Population size
#' @slot country Name of the country in which PPS was conducted
#' @slot bhai_options Options with which bhai was run. If bhai was not run yet, this is an empty list.
#' @slot bhai_summary Summary statistics of bhai. If bhai was not run yet, this is an empty list.
#' @exportClass PPS
.PPS <- setClass("PPS", representation(infections = "character",
                                       num_hai_patients = "numeric", 
                                       num_survey_patients = "numeric", 
                                       length_of_stay = "numeric",
                                       loi_pps = "list",
                                       hospital_discharges  = "matrix",
                                       num_hai_patients_by_stratum = "list",
                                       num_hai_patients_by_stratum_prior = "list",
                                       mccabe_scores_distr = "list",
                                       mccabe_by_stratum_prior = "list",
                                       mccabe_life_exp = "list",
                                       num_survey_patients_by_stratum = "list",
                                       population = "numeric",
                                       country = "character",
                                       bhai_options = "list",
                                       bhai_summary = "list"))


#' Checks validity of PPS parameters
#' 
#' @keywords internal
#' @noRd
checkPPSParams = function(num_hai_patients, 
                          num_survey_patients, 
                          length_of_stay,
                          loi_pps,
                          hospital_discharges ,
                          num_hai_patients_by_stratum,
                          num_hai_patients_by_stratum_prior,
                          mccabe_scores_distr,
                          mccabe_by_stratum_prior,
                          mccabe_life_exp,
                          num_survey_patients_by_stratum,
                          population,
                          country) {
  
  
  if(! inherits(num_hai_patients, what=c("numeric", "integer")))
    stop("num_hai_patients must be of a numeric or integer!")
  if(! inherits(num_survey_patients, what=c("numeric", "integer")))
    stop("num_survey_patients must be of a numeric or integer!")
  if(! inherits(length_of_stay, what=c("numeric", "integer")))
    stop("num_survey_patients must be of a numeric or integer!")
  if(length(length_of_stay) != 1)
    stop("length_of_stay must have length 1!")
  if(! inherits(loi_pps, what="list"))
    stop("loi_pps must be of a list!")
  sapply(loi_pps, function(x) 
    if(! inherits(x, what=c("numeric", "integer"))) {
      stop("Elements of loi_pps must be numeric or integer vectors!")
    })
  if(! inherits(hospital_discharges, what=c("numeric", "integer", "matrix")))
    stop("hospital_discharges must be of a numeric or integer or a matrix!")
  if(inherits(hospital_discharges, what=c("matrix"))) {
    checkMatrix(hospital_discharges)
  }
  else {
    if(length(hospital_discharges) != 1)
      stop("hospital_discharges must have length 1!")
  }
  
  
  if(! inherits(num_hai_patients_by_stratum, what="list"))
    stop("num_hai_patients_by_stratum must be of a list!")
  sapply(num_hai_patients_by_stratum, function(x) 
    checkMatrix(x, "num_hai_patients_by_stratum"))

  
  if(! inherits(num_hai_patients_by_stratum_prior, what="list"))
    stop("num_hai_patients_by_stratum_prior must be of a list!")
  if(length(num_hai_patients_by_stratum_prior) > 0) {
    if(is.null(names(num_hai_patients_by_stratum_prior)))
      stop("num_hai_patients_by_stratum_prior must be named!")
    sapply(num_hai_patients_by_stratum_prior, function(x) checkMatrix(x, "num_hai_patients_by_stratum_prior"))
    sapply(num_hai_patients_by_stratum_prior, function(x) 
      if(! all(names(x) %in% c("NONFATAL", "RAPFATAL", "ULTFATAL"))) {
        stop("num_hai_patients_by_stratum_prior must contain McCabe scores: NONFATAL, RAPFATAL, ULTFATAL")
      })
    if(! all(names(num_hai_patients_by_stratum_prior) %in% names(num_hai_patients)))
      stop("Names of num_hai_patients_by_stratum_prior do not match names of num_hai_patients!")
  }
 
  if(! inherits(mccabe_by_stratum_prior, what="list"))
    stop("mccabe_by_stratum_prior must be of a list!")
  if(length(mccabe_by_stratum_prior) > 0) {
    if(is.null(names(mccabe_by_stratum_prior)))
      stop("mccabe_by_stratum_prior must be named!")
    sapply(mccabe_by_stratum_prior, function(x) sapply(x, function(y) 
      checkMatrix(y, "mccabe_by_stratum_prior")))
    sapply(mccabe_by_stratum_prior, function(x) 
      if(! all(names(x) %in% c("NONFATAL", "RAPFATAL", "ULTFATAL"))) {
        stop("mccabe_by_stratum_prior must contain McCabe scores: NONFATAL, RAPFATAL, ULTFATAL")
      })
    if(all(names(mccabe_by_stratum_prior) %in% names(num_hai_patients)))
      stop("Names of mccabe_by_stratum_prior do not match names of num_hai_patients!")
  }
  
  if(! inherits(mccabe_scores_distr, what="list"))
  if(is.null(names(mccabe_scores_distr)))
    stop("mccabe_scores_distr must be named!")
  sapply(mccabe_scores_distr, function(x) sapply(x, function(y) 
    checkMatrix(y, "mccabe_scores_distr")))
  sapply(mccabe_scores_distr, function(x) 
    if(! all(names(x) %in% c("NONFATAL", "RAPFATAL", "ULTFATAL"))) {
      stop("mccabe_scores_distr must contain McCabe scores: NONFATAL, RAPFATAL, ULTFATAL")
    })
  if(! all(names(mccabe_scores_distr) %in% names(num_hai_patients)))
    stop("Names of mccabe_scores_distr do not match names of num_hai_patients!")
  
  if(! inherits(mccabe_life_exp, what="list"))
    if(is.null(names(mccabe_life_exp)))
      stop("mccabe_life_exp must be named!")
  sapply(mccabe_life_exp, function(x) checkMatrix(x, "mccabe_life_exp"))
  if(! all(names(mccabe_life_exp) %in% c("NONFATAL", "RAPFATAL", "ULTFATAL"))) {
    stop("mccabe_life_exp must contain McCabe scores: NONFATAL, RAPFATAL, ULTFATAL")
  }

  if(! inherits(num_survey_patients_by_stratum, what="list"))
    stop("num_survey_patients_by_stratum must be of a list!")
  if(length(num_survey_patients_by_stratum) > 0) {
    if(is.null(names(num_survey_patients_by_stratum)))
      stop("num_survey_patients_by_stratum must be named!")
    sapply(num_survey_patients_by_stratum, function(x) checkMatrix(x, "num_survey_patients_by_stratum"))
    sapply(num_survey_patients_by_stratum, function(x) 
      if(! all(names(x) %in% c("NONFATAL", "RAPFATAL", "ULTFATAL"))) {
        stop("num_survey_patients_by_stratum must contain McCabe scores: NONFATAL, RAPFATAL, ULTFATAL")
      })
    if(all(names(num_survey_patients_by_stratum) %in% names(num_hai_patients)))
      stop("Names of num_survey_patients_by_stratum do not match names of num_hai_patients!")
  }
  
  if(! inherits(population, what=c("numeric", "integer"))) 
    stop("population must be of a numeric or integer!")
  if(! inherits(country, what=c("character"))) 
    stop("country must be of a character!")
  
  if(is.null(names(num_hai_patients)))
    stop("num_hai_patients must be named!")
  if(is.null(names(loi_pps)))
    stop("loi_pps must be named!")
  if(is.null(names(loi_pps)))
    stop("mccabe_scores_distr must be named!")
  if(! all(names(loi_pps) %in% names(num_hai_patients)))
    stop("Names of loi_pps do not match names of num_hai_patients!")
  if(is.null(names(num_hai_patients_by_stratum)))
    stop("num_hai_patients_by_stratum must be named!")
  if(! all(names(num_hai_patients_by_stratum) %in% names(num_hai_patients)))
    stop("Names of num_hai_patients_by_stratum do not match names of num_hai_patients!")
  if(! all(names(mccabe_scores_distr) %in% names(num_hai_patients)))
    stop("Names of mccabe_scores_distr do not match names of num_hai_patients!")
  
}


checkMatrix = function(x, mat_name) {
  
  if(! inherits(x, what=c("matrix", "data.frame"))) {
    stop(paste("Elements of ", mat_name, " must be a matrix!", sep=""))
  }
  apply(x, 2, function(y)
    if(!inherits(y, what=c("numeric", "integer"))) {
      stop(paste(mat_name, "  elements must be numeric or integer!", sep=""))})
  if(! all(colnames(x) == c("F", "M") | colnames(x) == c("Female", "Male"))) {
    stop(paste("Colnames of ", mat_name, " must be c(\"F\", \"M\")", sep=""))
  } 
  rowN_print = "c(\"[0;0]\", \"[1;4]\", \"[5;9]\", \"[10;14]\", \"[15;19]\", \"[20;24]\", \"[25;29]\", \"[30;34]\", \"[35;39]\", \"[40;44]\", \"[45;49]\", \"[50;54]\", \"[55;59]\", \"[60;64]\", \"[65;69]\", \"[70;74]\", \"[75;79]\", \"[80;84]\", \"[85;Inf]\")"
  rowN = c("[0;0]", "[1;4]", "[5;9]", "[10;14]", "[15;19]", "[20;24]", "[25;29]", "[30;34]", "[35;39]", "[40;44]", "[45;49]", "[50;54]", "[55;59]", "[60;64]", "[65;69]", "[70;74]", "[75;79]", "[80;84]", "[85;Inf]")
  if(! all(rownames(x) == rowN)) {
    stop(paste("Colnames of ", mat_name, " must be ", rowN_print, sep=""))
  }  
}

#'
#' This function creates a PPS object.
#' @title Create a PPS object
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
#' @param num_survey_patients_by_stratum Number of survey patients stratified by infection, age and gender. If this parameter is provided the methodology described in Cassini et al. (2016) <doi:https://doi.org/10.1371/journal.pmed.1002150> is applied.
#' @param population Population size.
#' @param country Name of the country.
#' 
#' @return A PPS class object.
#' 
#' @seealso \code{\linkS4class{PPS}}
#' 
#' @examples
#' 
#' data(german_pps_2011_repr)
#' german_pps_repr = PPS(num_hai_patients = num_hai_patients,
#'     num_hai_patients_by_stratum = num_hai_patients_by_stratum,
#'     num_hai_patients_by_stratum_prior = num_hai_patients_by_stratum_prior,
#'     num_survey_patients = num_survey_patients,
#'     length_of_stay = length_of_stay,
#'     loi_pps = loi_pps,
#'     mccabe_scores_distr = mccabe_scores_distr,
#'     mccabe_life_exp = mccabe_life_exp,
#'     hospital_discharges = hospital_discharges,
#'     population = population,
#'     country="Germany (representative sample)")
#' german_pps_repr
#' 
#' @export 
PPS <- function(num_hai_patients = NULL, 
                num_survey_patients = NULL, 
                length_of_stay = NULL,
                loi_pps = NULL,
                hospital_discharges  = NULL,
                num_hai_patients_by_stratum = NULL,
                num_hai_patients_by_stratum_prior = NULL,
                mccabe_scores_distr = NULL,
                mccabe_by_stratum_prior = NULL,
                mccabe_life_exp = NULL,
                num_survey_patients_by_stratum = NULL,
                population = NULL,
                country = "") {
  
  
  if(is.null(num_hai_patients_by_stratum_prior)) {
    num_hai_patients_by_stratum_prior = list()  
  }
  if(is.null(mccabe_by_stratum_prior)) {
    mccabe_by_stratum_prior = list()  
  }
  if(is.null(num_survey_patients_by_stratum)) {
    num_survey_patients_by_stratum = list()  
  }
  
  checkPPSParams(num_hai_patients, 
                 num_survey_patients, 
                 length_of_stay,
                 loi_pps,
                 hospital_discharges ,
                 num_hai_patients_by_stratum,
                 num_hai_patients_by_stratum_prior,
                 mccabe_scores_distr,
                 mccabe_by_stratum_prior,
                 mccabe_life_exp,
                 num_survey_patients_by_stratum,
                 population,
                 country)
  
  if(class(hospital_discharges) == "numeric") {
    hospital_discharges = matrix(hospital_discharges,nrow=1,ncol=1)  
  }
  
  .PPS(infections = names(num_hai_patients),
       num_hai_patients = num_hai_patients, 
       num_survey_patients = num_survey_patients, 
       length_of_stay = length_of_stay,
       loi_pps = loi_pps,
       hospital_discharges  = hospital_discharges,
       num_hai_patients_by_stratum = num_hai_patients_by_stratum,
       num_hai_patients_by_stratum_prior = num_hai_patients_by_stratum_prior,
       mccabe_scores_distr = mccabe_scores_distr,
       mccabe_by_stratum_prior = as.list(mccabe_by_stratum_prior),
       mccabe_life_exp = mccabe_life_exp,
       num_survey_patients_by_stratum = as.list(num_survey_patients_by_stratum),
       population = population,
       country = country,
       bhai_options = list(),
       bhai_summary = list())
}







#' Prints description of PPS object
#' 
#' @keywords internal
#' @noRd
setMethod(f = "show", signature = c("PPS"), function(object) {
  if(object@country != "") {
    cat("  Country: ", object@country, "\n", sep="")
  }
  cat("  Patients in PPS: ", prettyNum(object@num_survey_patients, big.mark=","), "\n", sep = "")
  cat("  Patients with HAI: ",  prettyNum(sum(object@num_hai_patients), big.mark=","), "\n", sep = "")
  cat("  HAIs: ", paste(names(object@num_hai_patients), collapse=", "), "\n", sep="")

  if(length(object@bhai_summary) > 0) {
    mat = matrix(bhai.prettyTable(object)["ALL",c("Cases", "Deaths", "DALY")], ncol=1)
    rownames(mat) = c("  Cases:", "  Deaths:", "  DALY:")
    colnames(mat) = ""
    
    print(mat, quote=FALSE)
    if(object@bhai_options$stratified_sampling) {
      cat("Note: Results were obtained using stratified sampling (not recommended)!")
    }
  }
  
    
  
})



#' Estimation of the burden of healthcare-associated infections
#' 
#' @name bhai
#' @rdname bhai-methods
#'
#' @title Main function of the package to estimation of the burden of healthcare-associated infections
#' 
#' @param pps The PPS object containing the data.
#' @param nsim Number of Monte Carlo simulations, default: 1000.
#' @param pop.sampling Specifying whether parameters of the disease outcome trees should be sampled on population level, default: TRUE.
#' @param sample_distr Distribution used for prevalence sampling, default: "rbetamix".
#' @param estimate_loi_fun Function used for estimation of the length of infection, default: bootstrap_mean_gren (recommended!).
#' @param stratified_sampling Specifying whether stratified sampling should be done.
#' @param summarize_strata Specifying whether stratum-specific summary statistics should be computed.
#' @param use_prior Specifying whether Prior distributions should be used for computations.
#' 
#' @return A PPS class object.
#' 
#' @seealso \code{\linkS4class{PPS}}
#' 
#' @examples 
#' 
#' data(german_pps_2011_repr)
#' german_pps_repr = PPS(num_hai_patients = num_hai_patients,
#'     num_hai_patients_by_stratum = num_hai_patients_by_stratum,
#'     num_hai_patients_by_stratum_prior = num_hai_patients_by_stratum_prior,
#'     num_survey_patients = num_survey_patients,
#'     length_of_stay = length_of_stay,
#'     loi_pps = loi_pps,
#'     mccabe_scores_distr = mccabe_scores_distr,
#'     mccabe_life_exp = mccabe_life_exp,
#'     hospital_discharges = hospital_discharges,
#'     population = population,
#'     country="Germany (representative sample)")
#' german_pps_repr
#' 
#' set.seed(3)
#' # The following example is run only for illustratory reasons
#' # Note that you should never run the function with only 10 Monte-Carlo simulations in practice!
#' bhai(german_pps_repr, nsim=10)
#' 
#' @export
#' @docType methods
#' @rdname bhai-methods
setGeneric("bhai", function(pps, nsim=1000, pop.sampling=TRUE, sample_distr="rbetamix", estimate_loi_fun=bootstrap_mean_gren, stratified_sampling=FALSE, summarize_strata=TRUE, use_prior=TRUE) { standardGeneric("bhai") })


#' @rdname bhai-methods
#' @aliases bhai,PPS,ANY-method
setMethod("bhai", "PPS", function(pps, nsim=1000, pop.sampling=TRUE, sample_distr="rbetamix", estimate_loi_fun=bootstrap_mean_gren, stratified_sampling=FALSE, summarize_strata=TRUE, use_prior=TRUE)
{
  
  if(!use_prior) {
    pps@mccabe_by_stratum_prior = list()
    pps@mccabe_by_stratum_prior = list()
  }
  
  res = bhai.run(pps@num_hai_patients, pps@num_survey_patients, 
                 pps@length_of_stay, pps@loi_pps, 
                 pps@hospital_discharges,
                 pps@num_hai_patients_by_stratum, 
                 num_hai_patients_by_stratum_prior=pps@num_hai_patients_by_stratum_prior,  
                 pps@mccabe_scores_distr, mccabe_by_stratum_prior=pps@mccabe_by_stratum_prior, 
                 pps@mccabe_life_exp,
                 pop.sampling=pop.sampling, nsim=nsim, sample_distr=sample_distr,
                 infections=names(pps@num_hai_patients), 
                 num_survey_patients_by_stratum=pps@num_survey_patients_by_stratum,
                 estimate_loi_fun=estimate_loi_fun,
                 stratified_sampling=stratified_sampling)
  
  sim_options = list(nsim=nsim, pop.sampling=pop.sampling, 
                     sample_distr=sample_distr, estimate_loi_fun=estimate_loi_fun,
                     stratified_sampling=stratified_sampling)
  
  summary_bhai = summary.bhai(res, pps@population, summarize_strata)
  pps@bhai_options = sim_options
  pps@bhai_summary = summary_bhai
  pps
  
})


