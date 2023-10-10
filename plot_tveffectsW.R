plot_tveffectsW <- function (objTVeff,objTVeff2, objCeff, lag=FALSE,cut_time = NULL, ccol, tvcol, conf.int.const=c("TRUE", "FALSE"),
                            scale = c("logHR", "HR", "slope"), include = c("const", "tv", "both")) {
  scale <- match.arg(scale)
  include <- match.arg(include)
  tv_effect <- plot(objTVeff, "tv_effect", return_values = TRUE)[[1]]
  tv_effect2 <- plot(objTVeff2, "tv_effect", return_values = TRUE)[[1]]
  tv_effect[,3]<-(tv_effect2[,3]+tv_effect[,3])/2
  tv_effect[,2]<-(tv_effect2[,2]+tv_effect[,2])/2
  tv_effect[,4]<-(tv_effect2[,4]+tv_effect[,4])/2
  const_effect <- cbind(times = tv_effect[, 1],
                        est = objCeff$statistics$postMeans$alphas,
                        low = objCeff$statistics$CIs$alphas[1, ],
                        upp = objCeff$statistics$CIs$alphas[2, ])
  
  if (!is.null(cut_time)) {
    tv_effect <- tv_effect[tv_effect[, "times"] <= cut_time, ]
    const_effect <- const_effect[const_effect[, "times"] <= cut_time, ]
  }
  if (scale == "HR") {
    tv_effect[, -1] <- exp(tv_effect[, -1])
    const_effect[, -1] <- exp(const_effect[, -1])
  }
  matplot(tv_effect[, 1], cbind(tv_effect[, -1], const_effect[, -1]), type = "n",
          xlab = "Time (years)", 
          ylab = if (scale %in% c("HR", "slope")) expression("Value + Slope of GPI("~t[w]* ~ hat(X) *")") else expression("Value of GPI("~t[w]* ~ hat(X) *")"))
  
  if (include %in% c("const", "both")) {
    polygon(c(const_effect[, 1], rev(const_effect[, 1])), 
            c(const_effect[, 3], rev(const_effect[, 4])), 
            col = NULL, border = NA)
  }
  if (include %in% c("tv", "both")) {
    polygon(c(tv_effect[, 1], rev(tv_effect[, 1])), 
            c(tv_effect[, 3], rev(tv_effect[, 4])), 
            col = NULL,border = NA)
  }
  
  if (conf.int.const =="TRUE") {
    matlines(const_effect[, 1], const_effect[, -1], lty = c(1, 3, 3), 
             col = ccol, lwd = c(4, 4, 3))
  }
  if (conf.int.const =="FALSE") {
    matlines(const_effect[, 1], const_effect[, -c(1,3,4)], lty = c(1, 3, 3), 
             col = ccol, lwd = c(4, 4, 3))
  }
  if (include %in% c("tv", "both")) {
    matlines(tv_effect[, 1], tv_effect[, -1], lty = c(1, 3, 3), 
             col = tvcol, lwd = c(4, 3, 3))
  }
  if (lag =="TRUE") {
    leg <- legend("topleft","",bty="n")$text
    text(2+leg$x,leg$y, "Lag = 0.25years",adj=c(1,1))
  }
  
}
