#
# rapa_utils.R
#
# Utility functions for rapamycin analysis
#
# 2017-05-09  WTR
#


Get_file = function(fqp, compensate = TRUE, transform = TRUE, nice_names = TRUE, derail = TRUE) {
  ff = read.FCS(fqp)

  if (compensate) {
    ff = autocomp(ff)
  }

  if (transform) {
    fl_params = which(colnames(ff) %in% colnames(keyword(ff)$SPILL))
    sc_params = 1:6
    ff = doTransform(f = ff, cols = fl_params, method = 'biexp')
    ff = doTransform(f = ff, cols = sc_params, method = 'linear', fac = 5.4 / 262143)
  }

  if (nice_names) {
    colnames(ff)[fl_params] = parameters(ff)$desc[fl_params]
  }

  if (derail) {
    ff = derail(ff = ff, parameters = c("FSC-A", "SSC-A"))
  }

  ff
}

Gate_debris = function(ff, bandwidth = .1, gridsize = 1001, gheight = .75, show = FALSE) {
  kde = bkde(exprs(ff)[, "FSC-A"], bandwidth = bandwidth, gridsize = gridsize)
  kde$y = kde$y / max(kde$y)
  gde = fit_gaussian_base(kde = kde, height = gheight)
  gde$y = gde$y / max(gde$y)

  thresh_debris = cumm_prob_thresh(gde)

  ff_gated = Subset(ff, rectangleGate("FSC-A" = c(thresh_debris, Inf)))

  if (show) {
    pplot(ff, c("FSC-A", "SSC-A"), tx = 'linear', ty = 'linear')
    lines(kde$x, 4 * kde$y)
    lines(gde$x, 4 * gde$y, col = 'red')
    xline(thresh_debris, col = 'red', lwd = 3)
  }

  ff_gated
}
