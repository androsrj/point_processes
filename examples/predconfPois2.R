predconfPois2 <- function(region, object, level, N, m,
                         what=c("estimate", "se",
                                "confidence", "prediction"),
                         new.coef=NULL) {
  what <- match.arg(what)
  stopifnot(0 < level && level < 1)
  lam <- predict(object, window=region, new.coef=new.coef) * (N / m)
  mu.hat <- integral.im(lam)
  if(what == "estimate") return(mu.hat)
  mo <- model.images(object, W=as.owin(lam))
  ZL <- unlist(lapply(mo,
                      function(z, w) integral.im(eval.im(z * w)),
                      w = lam))
  ZL <- matrix(ZL, nrow=1)
  vc <- vcov(object, new.coef=new.coef)
  var.muhat <- as.numeric(ZL %*% vc %*% t(ZL))
  sd.muhat <- sqrt(var.muhat) / sqrt(N / m)
  if(what == "se") return(sd.muhat)
  alpha2 <- (1-level)/2
  pp <- sort(c(alpha2, 1-alpha2))
  out <- switch(what,
                confidence = mu.hat + qnorm(pp) * sd.muhat,
                prediction = qmixpois(pp, mu.hat, sd.muhat, I))
  names(out) <- paste0(signif(100 * pp, 3), "%")
  out
}
