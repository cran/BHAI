
#'  
#' Estimate the length of infection using the rearrengment and grenander estimator
#' 
#' @title Estimate the length of infection
#' 
#' @param loi Numeric vector containing the length of infections (LOI_PPS).
#' 
#' @usage runif_bootstrap_rear_gren(loi)
#' 
#' @return Estimated length of infection.
#' 
#' @keywords internal
#' @noRd
runif_bootstrap_rear_gren = function(loi) {
  loi = sample(loi, length(loi), replace = TRUE)
  loi = c(calculate_I_smooth(data = data.frame(A.loi=loi),
                             method = "rear")$x.loi.hat, 
          calculate_I_smooth(data = data.frame(A.loi=loi),
                             method = "gren")$x.loi.hat)
  loi = sort(loi)
  loi = runif(1, loi[1], loi[2])
  loi
}


#'  
#' Estimate the length of infection using the mean
#' 
#' @title Estimate the length of infection
#' 
#' @param loi Numeric vector containing the length of infections (LOI_PPS).
#' 
#' @usage bootstrap_mean(loi)
#' 
#' @return Estimated length of infection.
#' 
#' @keywords internal
#' @noRd
bootstrap_mean = function(loi) {
  loi = mean(sample(loi, length(loi), replace = TRUE))
  loi
}

#'  
#' Estimate the length of infection using the mean
#' 
#' @title Estimate the length of infection
#' 
#' @param loi Numeric vector containing the length of infections (LOI_PPS).
#' 
#' @usage bootstrap_mean_gren(loi)
#' 
#' @return Estimated length of infection.
#' 
#' @keywords internal
#' @noRd
bootstrap_mean_gren = function(loi) {
  loi = sample(loi, length(loi), replace = TRUE)
  loi = c(mean(loi), 
          calculate_I_smooth(data = data.frame(A.loi=loi),
                             method = "gren")$x.loi.hat)
  
  n = length(loi)
  a = 0.01
  b = 500
  alpha <- exp(a*(n-b))/(1+exp(a*(n-b)))
  loi = sum(loi*c(1-alpha, alpha))
  loi
}

