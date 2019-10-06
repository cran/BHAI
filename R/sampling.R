

#'  
#' Calculate lower bound of Clopper-Pearson confidence interval for proportions
#' code is taken from binom.test()
#' 
#' @title Calculate lower bound of Clopper-Pearson confidence interval for proportions
#' 
#' @param x Number of positive cases
#' @param alpha confidence level
#' @param n sample size
#' 
#' @keywords internal
#' @noRd
p.L <- function(x, alpha, n) {
  if (x == 0) 
    0
  else qbeta(alpha, x, n - x + 1)
}

#'  
#' Calculate upper bound of Clopper-Pearson confidence interval for proportions
#' code is taken from binom.test()
#' 
#' @title Calculate lower bound of Clopper-Pearson confidence interval for proportions
#' 
#' @param x Number of positive cases
#' @param alpha confidence level
#' @param n sample size
#' 
#' @keywords internal
#' @noRd
p.U <- function(x, alpha, n) {
  if (x == n) 
    1
  else qbeta(1 - alpha, x + 1, n - x)
}



#'  
#' Sample HAI prevalence of a mixture beta distributions
#' 
#' @title  Sample HAI prevalence of a mixture beta distributions
#' 
#' @param size number of values to sample
#' @param x Number of positive cases
#' @param n sample size
#' @param n Mixture probability
#' 
#' @keywords internal
#' @noRd
rbetamix <- function(size,x,n,alpha=0.5){
  U <- runif(size)
  I <- as.numeric(U<alpha)
  y <- I*rbeta(size, x, n - x + 1)+
    (1-I)*rbeta(size, x + 1, n - x)
  return(y)
}


#'  
#' Sample HAI prevalence from a PERT distribution
#' Code is taken from https://www.riskamp.com/beta-pert
#' 
#' @title  Sample HAI prevalence from a PERT distribution
#' 
#' @param n number of values to sample
#' @param x.min Lower bound of distribution
#' @param x.min Upper bound of distribution
#' @param x.mode Mode of distribution
#' @param lambda Shape parameter
#' 
#' @keywords internal
#' @noRd
rpert <- function( n, x.min, x.max, x.mode, lambda = 4 ){
  
  if( x.min > x.max || x.mode > x.max || x.mode < x.min ) stop( "invalid parameters" );
  
  x.range <- x.max - x.min;
  if( x.range == 0 ) return( rep( x.min, n ));
  
  mu <- ( x.min + x.max + lambda * x.mode ) / ( lambda + 2 );
  
  # special case if mu == mode
  if( mu == x.mode ){
    v <- ( lambda / 2 ) + 1
  }
  else {
    v <- (( mu - x.min ) * ( 2 * x.mode - x.min - x.max )) /
      (( x.mode - mu ) * ( x.max - x.min ));
  }
  
  w <- ( v * ( x.max - mu )) / ( mu - x.min );
  return ( rbeta( n, v, w ) * x.range + x.min );
}


#'  
#' Sample disease outcome tree parameters from a uniform distribution
#' 
#' @title  Sample disease outcome tree parameters from a uniform distribution
#' 
#' @param nstrata Number of strata, i.e. size of sample to draw
#' @param lower Lower bound of uniform distribution
#' @param upper Upper bound of uniform distribution
#' @param pop.level Indicating whether outcome trees are sampled on population level 
#' @param n Number of McCabe scores
#' 
#' @keywords internal
#' @noRd
sample.runif = function(nstrata, lower, upper, pop.level, list.names=NULL, n=3) {
  curr_sample = list()
  if(pop.level) {
    curr_values = rep(runif(1, lower, upper), nstrata)
    for(i in 1:n) {
      curr_sample[[i]] = matrix(curr_values, ncol=2)
    }
  }
  else {
    for(i in 1:n) {
      curr_sample[[i]] = matrix(runif(nstrata, lower, upper), ncol=2)
    }
  }
  names(curr_sample) = list.names
  curr_sample
}

#'  
#' Sample disease outcome tree parameters from a PERT distribution
#' 
#' @title   Sample disease outcome tree parameters from a PERT distribution
#' 
#' @param nstrata Number of strata, i.e. size of sample to draw
#' @param min Lower bound of PERT distribution
#' @param max Upper bound of PERT distribution
#' @param mode Mode of PERT distribution
#' @param pop.level Indicating whether outcome trees are sampled on population level 
#' @param list.names Names of output list (McCabe scores)
#' @param n Number of McCabe scores
#' 
#' @keywords internal
#' @noRd
sample.rpert = function(nstrata, min, max, mode, pop.level, list.names=NULL, n=3) {
  curr_sample = list()
  if(pop.level) {
    curr_values = rep(rpert(1, min, max, mode), nstrata)
    for(i in 1:n) {
      curr_sample[[i]] = matrix(curr_values, ncol=2)
    }
  }
  else {
    for(i in 1:n) {
      curr_sample[[i]] = matrix(rpert(nstrata, min, max, mode), ncol=2)
    }
  }
  names(curr_sample) = list.names
  curr_sample
}



#'  
#' Sample number of cases for Monte Carlo simulation
#' 
#' @title Sample number of cases for Monte Carlo simulation
#' 
#' @param nhai Number of patients in PPS with HAI
#' @param npatients Total number of patients in PPS
#' @param la Average length of stay of patients in PPS
#' @param loi Length of infection from PPS as a numeric vector
#' @param discharges Number of yearly hospital discharges
#' @param mccabe_life_exp A list containing life expectancy according to McCabe scores
#' @param mccabe_scores_distr Number of cases by HAI and McCabe scores
#' @param ncases_by_stratum Number of cases stratified by HAI, age and gender
#' @param sample_distr The distribution for the prevalence sampling. Default is 'rbetamix' which uses the mid-P-Pearson-Clopper interval. Alternatives are 'rpert' - which uses the PERT distribution - and 'bcode' which uses the sampling procedures as described in Cassini et al. (2016) <doi:https://doi.org/10.1371/journal.pmed.1002150>.
#' @param age_prior The prior weight (counts) for each infection, age and gender stratum. This is used for smooting the age and gender distribution when small numbers are observed.
#' @param p_age Number of survey patients stratified by infection, age and gender. If this parameter is provided the methodology described in Cassini et al. (2016) <doi:https://doi.org/10.1371/journal.pmed.1002150> is applied.
#' @param mccabe_prior The prior weight (counts) for each infection, McCabe score, age and gender stratum. This is used for smooting the age and gender distribution when small numbers are observed.
#' @param estimate_loi Function to use for the estimation of length of infection. Default is runif_bootstrap_rear_gren.
#' 
#' @keywords internal
#' @noRd
sample.ncases = function(nhai, npatients, la, loi, discharges, mccabe_life_exp, mccabe_scores_distr, ncases_by_stratum, sample_distr="rbetamix", age_prior=NULL, p_age=NULL, mccabe_prior=NULL, estimate_loi) {
   
  loi = estimate_loi(loi[[1]])
  
 
  ncases_hai = NULL
  if(sample_distr == "rpert") {
    point_estimate_hai = nhai/npatients*la/loi*sum(discharges)
    point_estimate_confint = c(p.L(nhai, 0.025, npatients),p.U(nhai, 0.025, npatients)) * la/loi*sum(discharges)
    ncases_hai = rpert(1, point_estimate_confint[1], point_estimate_confint[2], point_estimate_hai)
  }
  else if(sample_distr == "rbetamix") {
    ncases_hai = rbetamix(1, nhai, npatients) * la/loi*sum(discharges)
  }
  else if(sample_distr == "bcode") {
    point_estimate_hai = nhai/npatients*la/loi*sum(discharges)
    point_estimate_confint = c(p.L(nhai, 0.025, npatients),p.U(nhai, 0.025, npatients)) * la/loi*sum(discharges)
    ncases_hai = rpert(1, point_estimate_confint[1], point_estimate_confint[2], point_estimate_hai)
    if(point_estimate_hai == 0) {
      ncases_hai = runif(1, point_estimate_hai, point_estimate_confint[2])
    }
  }
  else {
    stop("sample_distr must be one of c(\"rpert\", \"rbetamix\", \"bcode\")")
  }
  
  hai_age_distr = ncases_by_stratum
  if(length(age_prior) > 0) {
    age_gender_prior = age_prior
    hai_age_distr = matrix(MCmultinomdirichlet(as.vector(ncases_by_stratum), as.vector(age_gender_prior), mc=1), ncol=2)
    rownames(hai_age_distr) = rownames(ncases_by_stratum)
    colnames(hai_age_distr) = colnames(ncases_by_stratum)
  }
  else {
    hai_age_distr = hai_age_distr #+ (1/length(hai_age_distr))
    hai_age_distr = hai_age_distr/sum(hai_age_distr)
  }
  
  if(length(mccabe_prior) > 0) {
    mccabe_scores_distr = mccabe_scores_distr[names(mccabe_prior)]
    for(i in 1:nrow(mccabe_prior[[1]])) {
      for(j in 1:ncol(mccabe_prior[[1]])) {
        curr_sample =  as.vector(MCmultinomdirichlet(as.vector(sapply(mccabe_scores_distr, function(x) x[i,j])), 
                                                     as.vector(sapply(mccabe_prior, function(x) x[i,j])), mc=1))
        names(curr_sample) = names(mccabe_prior)
        for(k in names(mccabe_prior)) {
          mccabe_scores_distr[[k]][i,j] = curr_sample[k]
        }
      }
    }
 
    
  }
  else {
    all_mccabe = Reduce("+", mccabe_scores_distr)
    for(i in 1:length(mccabe_scores_distr)) {
      mccabe_scores_distr[[i]] = mccabe_scores_distr[[i]]/all_mccabe
      mccabe_scores_distr[[i]][is.na(mccabe_scores_distr[[i]])] = 0
    }
  }
  
  ncases_hai = lapply(mccabe_scores_distr, function(x) x*hai_age_distr*ncases_hai) 
  if(length(p_age) > 0) {
    p_age_all = sum(unlist(p_age))
    p_age_red = Reduce("+", p_age)
    p_age_mccabe = p_age
    for(k in names(mccabe_scores_distr)) {
      p_age_mccabe[[k]] = p_age[[k]]/p_age_red
      p_age[[k]] = p_age[[k]]/p_age_all
      
      curr_ncases_hai = ncases_by_stratum
      for(i in 1:nrow(ncases_by_stratum)) { # iterate over age groups
        for(j in 1:ncol(ncases_by_stratum)) { # iterate over gender
          curr_n = round(npatients * p_age[[k]][i,j])
          if(curr_n == 0) {
            curr_ncases_hai[i,j] = 0
          }
          else {
            # calculate prevalence for each stratum
            curr_ncases = round(ncases_by_stratum[i,j] * mccabe_scores_distr[[k]][i,j])
            if(curr_n < curr_ncases) {
              print(c(curr_ncases, curr_n, "NO!"))
            }
            if(sample_distr == "rpert") {
              curr_prev = (curr_ncases / curr_n) #/ p_age[i,j]
              curr_conf_int = c(p.L(curr_ncases, 0.025, curr_n), 
                                p.U(curr_ncases, 0.025, curr_n)) 
              curr_prev = curr_prev * la/loi*discharges[i,j] * p_age_mccabe[[k]][i,j]
              curr_conf_int = curr_conf_int * la/loi*discharges[i,j] * p_age_mccabe[[k]][i,j]
              
              curr_ncases_hai[i,j] = rpert(1, curr_conf_int[1], curr_conf_int[2], curr_prev) 
            }
            if(sample_distr == "bcode") {
              curr_ncases_hai[i,j] = rpert(1, curr_conf_int[1], curr_conf_int[2], curr_prev) 
              if(curr_prev == 0) {
                curr_ncases_hai[i,j] = runif(1, curr_prev, curr_conf_int[2]) 
              }
            }
            
            if(sample_distr == "rbetamix") {
             curr_ncases_hai[i,j] = rbetamix(1,curr_ncases,curr_n)* 
                la/loi*discharges[i,j] * p_age_mccabe[[k]][i,j]
            }
          }
         
        }
      }
      ncases_hai[[k]] = curr_ncases_hai
    }
  }
  names(ncases_hai) = names(mccabe_scores_distr)
  
  ncases_hai
}
