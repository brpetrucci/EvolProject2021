for (i in 1:nrow(key)) {
  comb <- i
  
  pars <- key[comb, ]
  
  meDur <- pars[["meDur"]]
  
  reStart <- meStart + meDur
  
  reDur <- pars[["reDur"]]
  
  meInt <- pars[["meInt"]]
  
  mu_me <- lambda_bg - log(1 - meInt) / meDur
  
  reInt <- pars[["reInt"]]
  
  lambda_re <- (log(reInt) - (lambda_bg - mu_me) * meDur) / reDur + mu_bg
  
  lambda_null <- stepfun(c(reStart, reStart + reDur),
                         c(lambda_bg, lambda_re, lambda_bg))
  
  mu_null <- stepfun(c(meStart, meStart + meDur),
                     c(mu_bg, mu_me, mu_bg))
  
  expNt <- Vectorize(function(t) {
    exp(integrate(function(x) lambda_null(x) - mu_null(x), 0, t)$value)
  })
  
  expTotalNt <- function(t) {
    integrate(function(x) lambda_null(x) * expNt(x), 0, t, 
              stop.on.error = FALSE)$value + 1
  }
  
  tMax <- uniroot(function(t) expTotalNt(t) - nExp, 
                  c(meStart + meDur + reDur, 
                    meStart + meDur + reDur + 10),
                  extendInt = "yes")$root
  
  bmSigma2 <- pars[["bmSigma2"]]
  
  stQ01 <- pars[["stQ01"]]
  stQ10 <- pars[["stQSum"]] - stQ01
  
  stQ <- matrix(c(0, stQ01, stQ10, 0), 2, 2)
  
  meanLambda <- function(l) {
    exponent <- Vectorize(function(t, l) {
      integrate(function(x) l - l/5*
                  expected.trait(stQ01, stQ10, x), 0, t)$value
    })
    
    integrate(function(t) t * (l - l/5 * expected.trait(stQ01, stQ10, t)) *
                exp(-exponent(t, l)), 0, Inf)$value
  }
  
  lambda_bg_0 <- uniroot(function(l) meanLambda(l) - 1/lambda_bg, 
                         interval = c(0.1, 1),
                         extendInt = "yes")$root
  
  lambdaModCont <- pars[["lambdaModCont"]]
  
  lambda <- function(t, traits) {
    reStart <- meStart + meDur
    reEnd <- reStart + reDur
    ifelse((t < reStart) || (t > reEnd),                    
           # BG
           lambda_bg_0 - lambda_bg_0/5 * traits[2],
           # RE
           lambda_re + lambdaModCont * lambda_re * traits[1])
  }
  
  meanMu <- function(m) {
    exponent <- Vectorize(function(t, m) {
      integrate(function(x) m - m/5*
                  expected.trait(stQ01, stQ10, x), 0, t)$value
    })
    
    integrate(function(t) t * (m - m/5 * expected.trait(stQ01, stQ10, t)) *
                exp(-exponent(t, m)), 0, Inf)$value
  }
  
  mu_me_0 <- uniroot(function(m) meanMu(m) - 1/mu_me, 
                     interval = c(0.1, 1),
                     extendInt = "yes")$root
  
  muModCont <- 0.05 
  
  mu <- function(t, traits) {
    ifelse((t < meStart) || (t > meStart + meDur),
           # BG
           mu_bg + muModCont * mu_bg  * traits[1],
           # ME
           mu_me_0 - mu_me_0/5 * traits[2])
  }
  
  expNt1 <- Vectorize(function(t) {
    exp(integrate(Vectorize(function(x) lambda(x, c(0, expected.trait(stQ01, stQ10, x))) - 
                                          mu(x, c(0, expected.trait(stQ01, stQ10, x)))), 0, t)$value)
  })
  
  expTotalNt1 <- function(t) {
    integrate(Vectorize(function(x) lambda(x, c(0, expected.trait(stQ01, stQ10, x))) * expNt1(x)), 0, t, 
              stop.on.error = FALSE)$value + 1
  }
  
  print(paste0(tMax, " ", reStart + reDur, ", Diff = ", tMax - (reStart + reDur)))
  print(paste0("Comb: ", comb))
  print(paste0("Null: ", expNt(tMax), ", Trait: ", expNt1(tMax)))
  print(paste0("Null: ", expTotalNt(tMax), ", Trait: ", expTotalNt1(tMax)))
}
