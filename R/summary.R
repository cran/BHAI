

#' Create BHAI summary table
#'
#' @title Create summary table
#' 
#' @param pps The PPS object containing the data.
#' @param pop_norm Indicating whether statistics should be computed per 100,000 population, default: TRUE.
#' @param conf.int Specifying whether confidence intervals should be computed, default: TRUE.
#' 
#' @usage bhai.prettyTable(pps, pop_norm=FALSE, conf.int=TRUE)
#' 
#' @seealso \code{\linkS4class{PPS}}
#' 
#' @return A data.frame containing the summarised results.
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
#' result = bhai(german_pps_repr, nsim=10)
#' bhai.prettyTable(result)
#'   
#' @export
bhai.prettyTable = function(pps, pop_norm=FALSE, conf.int=TRUE) {
  summary_bhai = pps@bhai_summary
  
  pop_size = pps@population
  if(pop_norm) {
    pop_size = pop_size/100000
  }
  else {
    pop_size = 1
  }
  sum_tab = lapply(c("Cases", "Deaths", "DALY", "YLL", "YLD"), function(n) 
    sapply(summary_bhai, function(x) 
    paste(prettyNum(round(x$TOTAL[n,2]/pop_size, ifelse(pop_size==1,0,1)), big.mark=","), ifelse(conf.int," (",""), 
          ifelse(conf.int,prettyNum(round(x$TOTAL[n,1]/pop_size, ifelse(pop_size==1,0,1)), big.mark=","),""), ifelse(conf.int," - ",""), 
          ifelse(conf.int,prettyNum(round(x$TOTAL[n,3]/pop_size, ifelse(pop_size==1,0,1)), big.mark=","),""), ifelse(conf.int,")",""), sep="")))
  sum_tab = do.call("rbind", sum_tab)
  rownames(sum_tab) = c("Cases", "Deaths", "DALY", "YLL", "YLD")
  sum_tab = t(sum_tab)
  sum_tab
}

#'  
#' Internal function for bhai summary calculation
#' 
#' @param allsim Monte Carlo simulations.
#' @param pop_size Population size.
#' @param summarize_strata Specifying whether stratum-specific summary statistics should be computed.
#' 
#' @keywords internal
#' @noRd
summary.bhai = function(allsim, pop_size, summarize_strata=TRUE) {
  
  out = c(summary_by_disease(allsim, pop_size, summarize_strata),
    list("ALL"=c(summary_total(allsim),
    ifelse(summarize_strata, summary_by_stratum(allsim), list("stratum_specific_results"=NULL)))))
  names(out[["ALL"]]) = c("TOTAL", "stratum_specific_results")
  out
}

#'  
#' Internal function for bhai summary calculation for a single infection
#' 
#' @param allsim Monte Carlo simulations.
#' @param pop_size Population size.
#' @param summarize_strata Specifying whether stratum-specific summary statistics should be computed.
#' 
#' @keywords internal
#' @noRd
summary_by_disease = function(allsim, pop_size, summarize_strata=TRUE) {
  out = list()
  summary_out = list()
  for(n in names(allsim)) {
    # print(n)
    out[[n]] = list()
    currsim = allsim[[n]]
    for(curr in names(currsim[[1]])) {
      # print(curr)
      currval = list()
      for(s in 1:length(currsim)) {
        currval[[s]] = c(sapply(currsim[[s]][[curr]], sum), "TOTAL"=sum(unlist(currsim[[s]][[curr]])))
      }
      out[[n]][[curr]] = t(apply(do.call("rbind", currval), 2, quantile, prob=c(0.025, 0.5, 0.975)))
      out[[n]][[curr]] = rbind(out[[n]][[curr]], apply(out[[n]][[curr]][1:3,], 2, sum))
      rownames(out[[n]][[curr]])[5] = "CHECK-WRONG"
    }
    
    alldaly = list()
    allcases = list()
    alldeath = list()
    allyll = list()
    allyld = list()
    for(s in 1:length(currsim)) {
      alldaly[[s]] = sum(sapply(currsim[[s]][3:length(currsim[[s]])], function(y) sum(unlist(y))))
      allyll[[s]] = sum(sapply(currsim[[s]][3], function(y) sum(unlist(y))))
      allyld[[s]] = sum(sapply(currsim[[s]][4:length(currsim[[s]])], function(y) sum(unlist(y))))
      allcases[[s]] = sum(sapply(currsim[[s]][1], function(y) sum(unlist(y))))
      alldeath[[s]] = sum(sapply(currsim[[s]][2], function(y) sum(unlist(y))))
    }
    
    outmat = rbind(quantile(unlist(alldaly), prob=c(0.025, 0.5, 0.975)),
                   quantile(unlist(allyll), prob=c(0.025, 0.5, 0.975)),
                   quantile(unlist(allyld), prob=c(0.025, 0.5, 0.975)),
                   quantile(unlist(allcases), prob=c(0.025, 0.5, 0.975)),#/sum(pop_demographics[,2:3])*100000
                   quantile(unlist(alldeath), prob=c(0.025, 0.5, 0.975)))#/sum(pop_demographics[,2:3])*100000
    
    rownames(outmat) = c("DALY", "YLL", "YLD", "Cases", "Deaths")
    
    straum_specific = NULL
    if(summarize_strata) {
      ncases_by_stratum = do.call("rbind", lapply(currsim, function(y) as.vector(Reduce("+", y$ncases_hai))))
      allncases_by_stratum = apply(ncases_by_stratum, 2, quantile, prob=c(0.025, 0.5, 0.975))
      colnames(allncases_by_stratum) = rep(rownames(currsim[[1]]$ncases_hai[[1]]), 2)
      allncases_by_stratum = list("F"=t(allncases_by_stratum[,1:19]), "M"=t(allncases_by_stratum[,20:38]))
      
      ndeath_by_stratum = do.call("rbind", lapply(currsim, function(y) as.vector(Reduce("+", y$ndeath_hai))))
      allndeath_by_stratum = apply(ndeath_by_stratum, 2, quantile, prob=c(0.025, 0.5, 0.975))
      colnames(allndeath_by_stratum) = rep(rownames(currsim[[1]]$ndeath_hai[[1]]), 2)
      allndeath_by_stratum = list("F"=t(allndeath_by_stratum[,1:19]), "M"=t(allndeath_by_stratum[,20:38]))
      
      daly_strat = do.call("rbind",lapply(currsim, function(x) unlist((Reduce("+", unlist(x[3:length(x)], recursive=FALSE))))))
      daly_strat = apply(daly_strat, 2, quantile, prob=c(0.025,0.5,0.975))
      daly_strat = list(F=t(daly_strat[,1:19]), M=t(daly_strat[,20:38]))
      rownames(daly_strat[[1]]) = rownames(currsim[[1]][[1]][[1]])
      rownames(daly_strat[[2]]) = rownames(currsim[[1]][[1]][[1]])
      
      yll_strat = do.call("rbind",lapply(currsim, function(x) unlist((Reduce("+", unlist(x[3], recursive=FALSE))))))
      yll_strat = apply(yll_strat, 2, quantile, prob=c(0.025,0.5,0.975))
      yll_strat = list(F=t(yll_strat[,1:19]), M=t(yll_strat[,20:38]))
      rownames(yll_strat[[1]]) = rownames(currsim[[1]][[1]][[1]])
      rownames(yll_strat[[2]]) = rownames(currsim[[1]][[1]][[1]])
      
      yld_strat = do.call("rbind",lapply(currsim, function(x) unlist((Reduce("+", unlist(x[4:length(x)], recursive=FALSE))))))
      if(n == "SSI") {
        yld_strat = do.call("rbind",lapply(currsim, function(x) as.vector((Reduce("+", unlist(x[4:length(x)], recursive=FALSE))))))
      }
      yld_strat = apply(yld_strat, 2, quantile, prob=c(0.025,0.5,0.975))
      yld_strat = list(F=t(yld_strat[,1:19]), M=t(yld_strat[,20:38]))
      rownames(yld_strat[[1]]) = rownames(currsim[[1]][[1]][[1]])
      rownames(yld_strat[[2]]) = rownames(currsim[[1]][[1]][[1]])
      
      straum_specific = list(ncases=allncases_by_stratum, ndeath=allndeath_by_stratum, DALY=daly_strat, YLL=yll_strat, YLD=yld_strat)
    }
    
    summary_out[[n]] = list(TOTAL=outmat, stratum_specific_results=straum_specific)
    
  }
  
  summary_out
  
}

#'  
#' Internal function for bhai summary calculation for over all infections in PPS
#' 
#' @param allsim Monte Carlo simulations.
#' 
#' @keywords internal
#' @noRd
summary_total = function(allsim) {
  all_out = list("Cases"=quantile(apply(sapply(allsim, function(x) (sapply(x, function(y) sum(unlist(y$ncases_hai))))), 1, sum), prob=c(0.025, 0.5, 0.975)),
                 "Deaths"=quantile(apply(sapply(allsim, function(x) (sapply(x, function(y) sum(unlist(y$ndeath_hai))))), 1, sum), prob=c(0.025, 0.5, 0.975)),
                 "DALY"=quantile(apply(sapply(allsim, function(x) (sapply(x, function(y) sum(unlist(y[3:length(y)]))))), 1, sum), prob=c(0.025, 0.5, 0.975)),
                 "YLL"=quantile(apply(sapply(allsim, function(x) (sapply(x, function(y) sum(unlist(y[3]))))), 1, sum), prob=c(0.025, 0.5, 0.975)),
                 "YLD"=quantile(apply(sapply(allsim, function(x) (sapply(x, function(y) sum(unlist(y[4:length(y)]))))), 1, sum), prob=c(0.025, 0.5, 0.975)))
  list("TOTAL"=do.call("rbind", all_out))
  
}


#'  
#' Internal function that calculates stratum-specific summary statistics
#' 
#' @param allsim Monte Carlo simulations.
#' 
#' @keywords internal
#' @noRd
summary_by_stratum = function(allsim) {
  ncases_by_stratum = lapply(allsim, function(x) do.call("rbind", lapply(x, function(y) as.vector(Reduce("+", y$ncases_hai)))))
  allncases_by_stratum = apply(Reduce("+", ncases_by_stratum), 2, quantile, prob=c(0.025, 0.5, 0.975))
  colnames(allncases_by_stratum) = rep(rownames(allsim[[1]][[1]][[1]][[1]]), 2)
  allncases_by_stratum = list("F"=t(allncases_by_stratum[,1:19]), "M"=t(allncases_by_stratum[,20:38]))
  
  ndeath_by_stratum = lapply(allsim, function(x) do.call("rbind", lapply(x, function(y) as.vector(Reduce("+", y$ndeath_hai)))))
  allndeath_by_stratum = apply(Reduce("+", ndeath_by_stratum), 2, quantile, prob=c(0.025, 0.5, 0.975))
  colnames(allndeath_by_stratum) = rep(rownames(allsim[[1]][[1]][[1]][[1]]), 2)
  allndeath_by_stratum = list("F"=t(allndeath_by_stratum[,1:19]), "M"=t(allndeath_by_stratum[,20:38]))
  
  
  yll_by_stratum = lapply(allsim, function(x) do.call("rbind", lapply(x, function(y) unlist(Reduce("+", y$yll_hai)))))
  allyll_by_stratum = apply(Reduce("+", yll_by_stratum), 2, quantile, prob=c(0.025, 0.5, 0.975))
  colnames(allyll_by_stratum) = rep(rownames(allsim[[1]][[1]][[1]][[1]]), 2)
  allyll_by_stratum = list("F"=t(allyll_by_stratum[,1:19]), "M"=t(allyll_by_stratum[,20:38]))
  
  
  daly_by_stratum = lapply(allsim, function(x) do.call("rbind", lapply(x, function(y) unlist(Reduce("+", unlist(y[3:length(y)], recursive=FALSE))))))
  alldaly_by_stratum = apply(Reduce("+", daly_by_stratum), 2, quantile, prob=c(0.025, 0.5, 0.975))
  colnames(alldaly_by_stratum) = rep(rownames(allsim[[1]][[1]][[1]][[1]]), 2)
  alldaly_by_stratum = list("F"=t(alldaly_by_stratum[,1:19]), "M"=t(alldaly_by_stratum[,20:38]))
  
  
  list("stratum_specific_results"=
         list("ncases"=allncases_by_stratum, 
              "ndeath"=allndeath_by_stratum, 
              "YLL"=allyll_by_stratum, 
              "DALY"=alldaly_by_stratum))
}

