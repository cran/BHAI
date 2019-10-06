

#' Summary plot of number of infections, deaths and DALYs
#' 
#' @title Summary plot of number of infections, deaths and DALYs
#' 
#' @param pps The PPS object containing the data.
#' @param infections Infections to be plotted.
#' @param main Title of plot.
#' @param xlim Limits of x-axis.
#' @param ylim Limits of y-axis.
#' 
#' @usage bhai.circleplot(pps, infections=NULL, main="", xlim=NULL, ylim=NULL)
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
#' result = bhai(german_pps_repr, nsim=10)
#' bhai.circleplot(pps=result)
#' 
#' @export 
bhai.circleplot = function(pps, infections=NULL, main="", xlim=NULL, ylim=NULL) {
  
  bhai_summary = pps@bhai_summary
  
  infections = names(bhai_summary)
  infections = infections[infections != "ALL"]
  
  ncases_by_hai = sapply(bhai_summary[infections], function(x) x$TOTAL["Cases",2])#/pop_size*100000
  ndeath_by_hai = sapply(bhai_summary[infections], function(x) x$TOTAL["Deaths",2])#/pop_size*100000
  dalys_by_hai = sapply(bhai_summary[infections], function(x) x$TOTAL["DALY",2])#/pop_size*100000
  dalys_by_hai = dalys_by_hai/sum(dalys_by_hai)*mean(ncases_by_hai)
  
  xLim = range(ncases_by_hai*min(c(dalys_by_hai/ncases_by_hai, ncases_by_hai/dalys_by_hai)), ncases_by_hai+1.5*dalys_by_hai)
  yLim = range(c(ndeath_by_hai, ndeath_by_hai))
  
  xyasp <- par("pin")
  xycr <- diff(c(xLim, yLim))[c(1, 3)]
  ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
  
  yLim = range(c(ndeath_by_hai-dalys_by_hai*1.2*ymult, ndeath_by_hai+dalys_by_hai*1.2*ymult))
  
  if(!is.null(ylim)) {
    yLim = ylim
  }
  if(!is.null(xlim)) {
    xLim = xlim
  }
  
  xAxis = axisTicks(xLim, log=FALSE)
  yAxis = axisTicks(yLim, log=FALSE)
  plot(ncases_by_hai, ndeath_by_hai, pch="", 
       xlim=xLim, 
       ylim=yLim, 
       xaxt="n", yaxt="n",
       xlab="Number of cases", ylab="Number of deaths", main=main)
  axis(1, xAxis, xAxis)
  axis(2, yAxis, yAxis)
  abline(v=xAxis, col="grey")
  abline(h=yAxis, col="grey")
 # segments(x0=xAxis, y0=rep(yLim[1]-dalys_by_hai*1.2, length(xAxis)), x1=xAxis, y1=rep(yLim[2]+dalys_by_hai*1.2, length(xAxis)), col="grey")
  #segments(y0=yAxis, x0=rep(xLim[1]-dalys_by_hai*1.2, length(yAxis)), y1=yAxis, x1=rep(xLim[2]+dalys_by_hai*1.2, length(yAxis)), col="grey")
  
  
  for(i in names(ncases_by_hai)) {
    draw.circle(ncases_by_hai[i], ndeath_by_hai[i],
                radius=dalys_by_hai[i], nv=10000,
                col="lightblue", border="darkblue")
    text(ncases_by_hai[i], ndeath_by_hai[i], i)
  }
} 


#' Stratified barplot of cases, deaths and DALYs.
#' 
#' @title Stratified barplot of cases, deaths and DALYs.
#' 
#' @param pps The PPS object containing the data.
#' @param infection Infection to be plotted.
#' @param what One of c("Cases", "Deaths", "DALY")
#' @param col  Color used to fill the bars.
#' @param errors Specifying whether error bars should be plotted, default: TRUE.
#' @param lwd.errors Line width of error bars.
#' @param xlab X-axis labels.
#' @param ... Further plotting arguments
#' 
#' @usage bhai.strataplot(pps, infection, what, col=NULL, errors=TRUE, lwd.errors=2, xlab=NULL, ...)
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
#' result = bhai(german_pps_repr, nsim=10)
#' bhai.strataplot(pps=result, infection="HAP", what="Cases")
#' 
#' @export 
bhai.strataplot = function(pps, infection, what, col=NULL, errors=TRUE, lwd.errors=2, xlab=NULL, ...) {
  
  bhai_summary = pps@bhai_summary[[infection]]
  data_by_stratum2 = NULL
  if(what == "DALY") {
    data_by_stratum = bhai_summary$stratum_specific_results$DALY
    data_by_stratum2 = bhai_summary$stratum_specific_results$YLL
    if(is.null(col)) {
      col = c("blue", "red")
    }
    col1 = col[1]
    col2 = col[2]
  }
  if(what == "Cases") {
    data_by_stratum = bhai_summary$stratum_specific_results$ncases
    if(is.null(col)) {
      col = c("grey")
    }
    col1 = col[1]
    col2 = col[1]
  }
  if(what == "Deaths") {
    data_by_stratum = bhai_summary$stratum_specific_results$ndeath
    if(is.null(col)) {
      col = c("grey")
    }
    col1 = col[1]
    col2 = col[1]
  }
  
  if(is.null(xlab)) {
    xlab=switch(what,
                Cases = "No. of cases ",
                Deaths =  "No. of deaths ",
                DALY = "DALY ")
    xlab = paste(xlab, "(", infection, ")", sep="")
  }
  
  oldpar <- par(mfrow=c(1,2), mar=c(5, 2, 4, 2.4))	
  on.exit(par(oldpar))			
  
  
  barcenters = barplot(-t(data_by_stratum[["F"]][,2]), main="Female", legend=TRUE, names.arg = rep("",19),
                       args.legend=list(x="topleft"), ylab="", las=2, xlab=xlab,
                       xlim=c(-max(unlist(data_by_stratum), na.rm=TRUE), 0), horiz=TRUE, col=col1, axes=FALSE, ...)
  if(!is.null(data_by_stratum2)) {
    barplot(-t(data_by_stratum2[["F"]][,2]), col=col2, add=TRUE, horiz=TRUE, axes=FALSE, names.arg = rep("",19), ...)
  }
  if(errors) {
    segments(-t(data_by_stratum[["F"]][,1]), barcenters,
             -t(data_by_stratum[["F"]][,3]), barcenters, lwd = lwd.errors)
    segments(-t(data_by_stratum[["F"]][,1]), barcenters+0.25,
             -t(data_by_stratum[["F"]][,1]), barcenters-0.25, lwd = lwd.errors)
    segments(-t(data_by_stratum[["F"]][,3]), barcenters+0.25,
             -t(data_by_stratum[["F"]][,3]), barcenters-0.25, lwd = lwd.errors)
  }
  
  currTicks = axisTicks(c(0,max(unlist(data_by_stratum), na.rm=TRUE)), log=FALSE)
  axis(1, at=-currTicks, currTicks, cex.axis=0.8)
  mynames = gsub("-0", "", gsub("-Inf", "+", gsub(";", "-", gsub("\\[|\\]", "", rownames(data_by_stratum[["F"]])))))
  barcenters = barplot(t(data_by_stratum[["M"]][,2]), main="Male", legend=TRUE, names.arg=mynames, 
                       args.legend=list(x="topleft"), ylab="", las=2, xlab=xlab,
                       xlim=c(0,max(unlist(data_by_stratum), na.rm=TRUE)), horiz=TRUE, col=col1, cex.axis=0.75, las=2, axes=FALSE, ...)
  if(!is.null(data_by_stratum2)) {
    barplot(t(data_by_stratum2[["M"]][,2]), col=col2, add=TRUE, horiz=TRUE, axes=FALSE, names.arg = rep("",19), ...)
  }
  if(errors) {
    segments(t(data_by_stratum[["M"]][,1]), barcenters,
             t(data_by_stratum[["M"]][,3]), barcenters, lwd = lwd.errors)
    segments(t(data_by_stratum[["M"]][,1]), barcenters+0.25,
             t(data_by_stratum[["M"]][,1]), barcenters-0.25, lwd = lwd.errors)
    segments(t(data_by_stratum[["M"]][,3]), barcenters+0.25,
             t(data_by_stratum[["M"]][,3]), barcenters-0.25, lwd = lwd.errors)
  }
  
  axis(1, at=currTicks, currTicks, cex.axis=0.8)

}


#' Barplot of cases, deaths and DALYs.
#' 
#' @title Barplot of cases, deaths and DALYs.
#' 
#' @param ... Further plotting arguments
#' @param what One of c("Cases", "Deaths", "DALY")
#' @param infections If sepcified only a subset of infections in \code{bhai_summary} is plotted.
#' @param cols1  Color used to fill the bars.
#' @param cols2 Specifies colors of YLDs when plotting DALYs.
#' @param ylab Y-axis labels.
#' @param ylim Limits of y-axis.
#' @param legend_labs Labels of legend.
#' @param main Title of plot
#' @param names.inf Specifying whether names of infections should be plotted.
#' @param cex.names Font size of labels.
#' @param border The color to be used for the border of the bars, default: par("fg").
#' @param lwd.errors Line width of error bars.
#' 
#' @usage bhai.barplot(..., what, infections=NULL, cols1=NULL, cols2=NULL, ylab=NULL, ylim=NULL, 
#' legend_labs=NULL, main="", names.inf=TRUE, cex.names=1, border=par("fg"), lwd.errors=2)
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
#' result_ger = bhai(german_pps_repr, nsim=10)
#' 
#' bhai.barplot(result_ger, what="Cases")
#' 
#' @export 
bhai.barplot = function(..., what, infections=NULL, cols1=NULL, cols2=NULL, ylab=NULL, ylim=NULL, legend_labs=NULL, main="", names.inf=TRUE, cex.names=1, border=par("fg"), lwd.errors=2) {
  
  pps_objects = list(...)
  if(class(pps_objects[[1]]) == "list") {
    pps_objects = pps_objects[[1]]
  }
  
  bhai_summaries = lapply(pps_objects, function(x) x@bhai_summary)
  names(bhai_summaries) = sapply(pps_objects, function(x) x@country)
  pop_demographics_list = lapply(pps_objects, function(x) x@population)

  if(is.null(ylab)) {
    ylab=switch(what,
                Cases = "No. of cases per 100,000",
                Deaths =  "No. of deaths per 100,000",
                DALY = "DALY per 100,000")
  }
  
  if(is.null(infections)) {
    infections = names(bhai_summaries[[1]])
    infections = infections[infections != "ALL"]
  }
  
  pop_norm = sapply(pop_demographics_list, function(x) 1/sum(x)*100000)
  
  median_all = do.call("rbind", lapply(infections, function(x) 
    c(sapply(1:length(bhai_summaries), function(y) bhai_summaries[[y]][[x]]$TOTAL[what,"50%"]*pop_norm[y]), NA)))
  confint_lower = as.vector(t(do.call("rbind", lapply(infections, function(x) 
    c(sapply(1:length(bhai_summaries), function(y) bhai_summaries[[y]][[x]]$TOTAL[what,"2.5%"]*pop_norm[y]), NA)))))
  confint_upper = as.vector(t(do.call("rbind", lapply(infections, function(x) 
    c(sapply(1:length(bhai_summaries), function(y) bhai_summaries[[y]][[x]]$TOTAL[what,"97.5%"]*pop_norm[y]), NA)))))
  
  if(is.null(cols1)) {
    cols1 = colorRampPalette(c("gray40", "lightgrey"))(length(bhai_summaries))
  }
  else {
    if(length(cols1) == 1) {
      cols1 = rep(cols1, length(bhai_summaries))
    }
  }
  
  median_all_2 = NULL
  if(what == "DALY") {
    median_all_2 = do.call("rbind", lapply(infections, function(x) 
      c(sapply(1:length(bhai_summaries), function(y) bhai_summaries[[y]][[x]]$TOTAL["YLL","50%"]*pop_norm[y]), NA)))
    if(is.null(cols2)) {
      cols1 = colorRampPalette(c("red", "gold"))(length(bhai_summaries))
      cols2 = colorRampPalette(c("blue", "lightblue"))(length(bhai_summaries))
    }
    else {
      if(length(cols2) == 1) {
        cols2 = rep(cols2, length(bhai_summaries))
      }
    }
  }
  
  confint_all = list()
  for(i in 1:length(bhai_summaries)) {
    confint_all[[i]] = do.call("rbind", lapply(infections, function(x) bhai_summaries[[i]][[x]]$TOTAL[what,c("2.5%", "97.5%")]))*pop_norm[[i]]
  }
  
  if(is.null(ylim)) {
    ylim = c(0, max(unlist(confint_upper), na.rm=TRUE))*1.2
  }
  
  curr_barplot = barplot(as.vector(t(median_all)), col=c(cols1, NA),   ylab=ylab, ylim=ylim, beside=TRUE, main=main, border=border)
 
  if(names.inf) {
    axis(1, at=curr_barplot[seq(ceiling(length(bhai_summaries)/2),length(curr_barplot)-ceiling(length(bhai_summaries)/2), 
                                length=length(infections))], labels=infections, cex.axis=cex.names)
  }
  else {
    axis(1, at=curr_barplot, colnames(median_all), las=2, cex.axis=cex.names)
  }
  #axis(1, at=seq(range(curr_barplot)[1], range(curr_barplot)[2], length=5), labels=infections)
  
  yAxis = axisTicks(ylim, log=FALSE)
  abline(h=yAxis, col="grey")
  
  barplot(as.vector(t(median_all)), col=c(cols1, NA), add=TRUE, beside=TRUE, border=border)
  if(!is.null(median_all_2)) {
    barplot(as.vector(t(median_all_2)), col=c(cols2, NA), add=TRUE, border=border)
  }
  
  segments(curr_barplot, confint_lower, curr_barplot, confint_upper, lwd = lwd.errors)
  segments(curr_barplot-0.1, confint_lower, curr_barplot+0.1, confint_lower, lwd = lwd.errors)
  segments(curr_barplot-0.1, confint_upper, curr_barplot+0.1, confint_upper, lwd = lwd.errors)
  
  
  
  if(what == "DALY") {
    if(is.null(legend_labs)) {
      legend_labs = unlist(lapply(unique(names(bhai_summaries)), function(x) c(paste(x, ":", sep=""), "YLD", "YLL")))
    }
    legend("topright", legend_labs, col=sapply(1:length(bhai_summaries), function(x) c(NA, cols1[x], cols2[x])), pch=15, bty="n")
  }
  else {
    if(is.null(legend_labs)) {
      legend_labs = unique(names(bhai_summaries))
    }
    legend("topright", legend_labs, col=cols1, pch=15, bty="n")
  }
  
}







