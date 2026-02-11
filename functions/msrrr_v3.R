## engine: msRRR algorithm via MM + SRRR
require(glmnet)
require(rrpack)


## matrix row-wise 2_1 norm for group sparsity 
norm21 <- function(A) sum(apply(A,1,function(a) sqrt(sum(a^2))))

##' @importFrom stats dnorm dpois dbinom
logLikehood <- function(Y, MU, Sigma = 1, family){
  switch(
    family,
    "gaussian" = dnorm(Y, MU, Sigma, log = TRUE), 
    "poisson"  = dpois(Y, MU, log = TRUE),
    "binomial" = dbinom(Y, 1, MU, log = TRUE)
  )
}

## neg log-lik, to minimize
objFun <- function(Y, mu, Phi, family, familygroup) {
  nfamily <- length(family) 
  cfamily <- unique(familygroup)
  n = nrow(Y)
  temp <- vector()
  for (j in 1:nfamily) {
    if(family[[j]]$family == "gaussian"){
      Sig = t(matrix(Phi[familygroup == cfamily[j]],
                     sum(familygroup == cfamily[j]), n))
    }else{Sig = 1}
    temp[j] <- -sum(logLikehood(Y[, familygroup == cfamily[j]],
                                mu[, familygroup == cfamily[j]],
                                sqrt(Sig), family[[j]]$family), na.rm = TRUE)
  }
  sum(temp)
}

# pseudo MSE...
pmse <- function(Y, mu, Phi=NULL, family, familygroup) {
  nfamily <- length(family) 
  cfamily <- unique(familygroup)
  n = nrow(Y)
  m = sum(!is.na(Y))
  Yvar = apply(Y,2,var,na.rm=T)
  Yvar = t(matrix(Yvar, ncol(Y), n))
  temp <- vector()
  for (j in 1:nfamily) {
    if(family[[j]]$family == "gaussian"){
      # Sig = t(matrix(Phi[familygroup == cfamily[j]],
      #                sum(familygroup == cfamily[j]), n))
      Sig = Yvar[, familygroup == cfamily[j]]
      temp[j] <- sum((Y-mu)[, familygroup == cfamily[j]]^2 / Sig, na.rm = T)
    }else if(family[[j]]$family == "binomial"){
      # Sig = ((1-mu)*mu)[, familygroup == cfamily[j]]
      Sig = Yvar[, familygroup == cfamily[j]]
      temp[j] <- sum((Y-mu)[, familygroup == cfamily[j]]^2 / Sig, na.rm = T)
    }else if(family[[j]]$family == "poisson"){
      # Sig = mu[, familygroup == cfamily[j]]
      Sig = Yvar[, familygroup == cfamily[j]]
      temp[j] <- sum((Y-mu)[, familygroup == cfamily[j]]^2 / Sig, na.rm = T)
    } 
  }
  return(sum(temp)/m)
}


## msrrr fit with prespecified nrank and lambda
# Z include a intercept col
msrrr.fit <- function(Y, X, Z, family, familygroup, nrank=2, lambda, init=NULL, 
                  control=list(epsilon=1e-4, maxit=200, trace=F)){ 
                  
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  # Z = cbind(1,Z)
  pz = ncol(Z)
  Y_mis = Y
  fg = family[familygroup]
  id.mis = which(is.na(Y))
  r = nrank 
  
  ## get init values via ridge
  if(is.null(init$A0)){  # use ridge to get init working Y
    C0 = matrix(NA, pz, q)
    B0 = matrix(NA, p, q)
    for(iq in 1:q){
      # fit0 <- glm(Y[,iq] ~ 0+Z+X, family = fg[[iq]])
      idx = !is.na(Y_mis[,iq])
      # ridge
      fit0 <- glmnet(cbind(Z,X)[idx,], Y[idx,iq], family=fg[[iq]], alpha=0, lambda=0.03, 
                     intercept=F, penalty.factor=c(rep(0,pz),rep(1,p))) 
      C0[,iq] = as.numeric(fit0$beta)[1:pz]
      B0[,iq] = as.numeric(fit0$beta)[-c(1:pz)]
    }
    svdB0 = svd(B0)
    A0 = svdB0$u[,1:r,drop=F] %*% diag(svdB0$d[1:r],r,r)
    V0 = svdB0$v[,1:r,drop=F] 
    init$A0 = A0
    init$V0 = V0
    init$C0 = C0
  } else{ 
    A0 = init$A0
    V0 = init$V0
    C0 = init$C0
  } 
   
  B0 = A0 %*% t(V0) # plot(Bt, B0)
  eta = (X%*%B0 + Z%*%C0) # nature parameter matrix
  mu = matrix(NA, n, q)   # mean matrix
  phi = rep(NA, q)        # dispersion par
  for(iq in 1:q) mu[,iq] = fg[[iq]]$linkinv(eta[,iq])       
  Y[id.mis] = mu[id.mis]  # Y: imputed working outcome 
  phi = ifelse(unlist(lapply(fg, function(a) a$family))=='gaussian', colMeans((Y_mis-mu)^2,na.rm=T), 1) 
  # plot(Y_mis, mu)
  
  ## scaling predictor matrix
  Kappa <- vector()
  if (is.null(init$kappa)) {
    svdX0d1 <- svd(cbind(Z,X))$d[1]
    for (j in 1:q) {
      Kappa[j] <- switch(
        fg[[j]]$family,
        'gaussian' = svdX0d1 / phi[j],
        'binomial' = svdX0d1 / 2,
        'poisson' = svdX0d1 * quantile(Y_mis[,j],0.9,na.rm=T)
      )
    }
    kappa <- max(Kappa)
    init$kappa = kappa
  }else{kappa = init$kappa }
  # kappa = 1
  X = X/kappa
  Z = Z/kappa
  A0 = A0*kappa
  B0 = B0*kappa
  C0 = C0*kappa
  
  C = C0
  B = B0
  A = A0
  V = V0
  # mean((B/kappa - Bt)^2) # 0.00417
  dif = obj = rep(NA, control$maxit+1)
  obj[1] = objFun(Y_mis, mu, phi, family, familygroup) + lambda*norm21(A0)# 
  for (iter in 1:control$maxit) {  
    # if (control$trace) message(gettextf("iteration %d", iter), domain = NA) 
    
    R = (t(X)%*%(Y-mu)%*%diag(1/phi) + B0)
    solve.srrr = srrr(Y=R,X=diag(1,p), nrank, A0=A0, V0=V0, modstr = list(lamA=c(lambda),nlam=1),
                      control=list(epsilon=control$epsilon, maxit=10))
    A = solve.srrr$A.path[1,,] 
    V = solve.srrr$V.path[1,,] 
    B = A%*%t(V)  
    # cat(sum(B[,1]==0))
    # plot(A0, A)
    
    ## update C
    for(iq in 1:q){
      # cat('.')
      # fine to use glm here as Z is low-d
      fit0 <- glm(Y_mis[,iq] ~ 0+Z, offset=X%*%B[,iq], family=fg[[iq]])
      C[,iq] = fit0$coef 
    }
    
    ## update working Y, and phi
    eta <- X%*%B + Z%*%C
    for(iq in 1:q){
      fm = fg[[iq]]
      mu[,iq] = fm$linkinv(eta[,iq])      
    }
    Y[id.mis] = mu[id.mis]       
    phi = ifelse(unlist(lapply(fg, function(a) a$family))=='gaussian', colMeans((Y_mis-mu)^2,na.rm=T), 1)
    
    ## conv  
    # if (sum((eta - etaold)^2) < tol * sum(eta^2))  break
    # dif[iter] <- sum((B - B0)^2)/(sum(B0^2) + 1e-6) 
    obj[iter+1] = objFun(Y_mis, mu, phi, family, familygroup) + lambda*norm21(A)
    dif[iter] = obj[iter+1] / obj[iter] - 1
    B0 = B
    A0 = A
    V0 = V
    C0 = C
    # if (dif[iter] < control$epsilon)  break
    if(dif[iter]>0 & control$trace==T) warning('obj increased')
    if(abs(dif[iter]) < control$epsilon) break
  }
  
  dif = dif[1:iter]
  obj = obj[1:(iter+1)]  
  
  A = A/kappa
  B = B/kappa
  C = C/kappa
  return(list(B=B,A=A,V=V,C=C,phi=phi,mu=mu, dif=dif,obj=obj,Y_imp=Y,Kappa=Kappa,kappa=kappa,iter=iter))
}


## wrapper: select lambda with prespecified nrank, by IC or CV selection 
# Z include a intercept col
msrrr.tuning = function(Y, X, Z, family, familygroup, nrank=2, init=NULL, 
                        lamseq=NULL, nlam=50, warm=T,  # lam.max, lam.min, 
                        method='CV', cv.criteria='pMSE', foldid=NULL, nfold=5, c.BIC=2, 
                        # c('CV', 'BIC', 'BICP', 'AIC', 'GIC')
                        control=list(epsilon=1e-4, maxit=200, trace=F, conv.obj=T)){
  n = nrow(Y)
  q = ncol(Y)
  p = ncol(X)
  # Z = cbind(1,Z)
  pz = ncol(Z) 
  Y_mis = Y
  fg = family[familygroup]
  id.mis = which(is.na(Y))
  r = nrank  
  
  ## distribution families ##
  nfamily <- length(family)
  if (nfamily == 1 & is.null(familygroup)) familygroup <- rep(1, q)
  ## characters of families
  cfamily <- unique(familygroup)
  
  ## default lambda seq
  if (is.null(lamseq)) {
    cat("No lambda sequence provided, use default!\n") 
    yt <- rep(0, n)
    for (i in 1:nfamily) {   
      yt <- yt + switch(family[[i]]$family,
                        'gaussian' = apply(abs(2 * Y[, familygroup == i]), 1, function(a) sum(a, na.rm = TRUE)),
                        'binomial' = apply(abs(2 - Y[, familygroup == i]), 1, function(a) sum(a, na.rm = TRUE)),
                        'poisson'  = apply(abs(2 * Y[, familygroup == i] - 2), 1, function(a) sum(a, na.rm = TRUE)) )
    }
    lam.max <- sum(yt * apply(abs(X), 1, max)) / 1000 ## ???? 30-30000 too large?
    lamseq <- 10 ^ (seq(log10(lam.max * 0.001), log10(lam.max), len = nlam)) 
    if(control$trace==T) cat('lambda:', lamseq)
  }
  nlam = length(lamseq)
  
  ## sol path
  Apath = array(NA, c(nlam,p,nrank))
  Vpath = array(NA, c(nlam,q,nrank))
  Cpath = array(NA, c(nlam,pz,q))
  phipath = matrix(NA, nlam,q)
  mupath = array(NA, c(nlam,n,q))
  for(il in 1:nlam){
    if(control$trace==T) cat(il,'...')
    init.i = init 
    if(il!=1 & warm==T) init.i=list(A0=fit$A, V0=fit$V, C0=fit$C, kappa=fit$kappa)
    fit = msrrr.fit(Y, X, Z, family, familygroup, nrank, lamseq[il], init.i, control)
    Apath[il,,] = fit$A
    Vpath[il,,] = fit$V
    Cpath[il,,] = fit$C
    phipath[il,] = fit$phi
    mupath[il,,] = fit$mu
  }
  
  ## tuning
  if(method!='CV'){ 
    xrank = min(n, p) # sum(svd(XX)$d > 0.0001)
    dev = rep(0, nlam) # -2*logL
    df = nz = rep(NA, nlam)
    for(il in 1:nlam){ 
      for(iq in 1:q){   
        # tt= logLikehood(Y[,iq],mupath[il,,iq],phipath[il,iq],fg[[iq]]$family) 
        dev[il] <- dev[il] - 2*sum(logLikehood(Y[,iq],mupath[il,,iq],phipath[il,iq],fg[[iq]]$family), na.rm=T)
        # objFun(Y_test, mu.test, fit1$dispersion)
      }
      # d.f. 
      dfu0 = sum(Apath[il,,] != 0) 
      dfv0 = sum(Vpath[il,,] != 0) 
      df[il] = dfu0 * xrank/p + dfv0 - nrank*nrank + pz*q
      nz[il] = mean(Apath[il,,1] != 0)
    }
    ## IC
    n.obs = sum(!is.na(Y_mis))
    logqn = log(n.obs)
    BIC = (dev + df*logqn)/n.obs 
    BICP = (dev + df*logqn*c.BIC)/n.obs 
    AIC = (dev + 2*df)/n.obs 
    GIC = (dev + df*log(logqn)*log(p*q))/n.obs 
    # dfqn2 = (1 - df/(q*n))^2 
    # GCV = dev/(q*n*dfqn2) 
    ICpath = data.frame(nz, df, dev, BIC, BICP, AIC, GIC)
    lam.idx = switch(method,
                     'BIC'=which.min(BIC),
                     'BICP'=which.min(BICP),
                     'AIC'=which.min(AIC),
                     'GIC'=which.min(GIC)) 
    Tunepath = ICpath
    tunepath.opt = switch(method,
                          'BIC'=min(BIC),
                          'BICP'=min(BICP),
                          'AIC'=min(AIC),
                          'GIC'=min(GIC))
  }else{ 
    # store the deviance/pseudo MSE of the test data
    cv <- matrix(NA, nlam, nfold) 
    if (is.null(foldid)) {  ## CV by row  
      foldid <- rep(1:nfold, len = n)
      foldid <- sample(foldid, n, replace = FALSE)
    }  
    
    init.i=init 
    for (ifold in 1:nfold) { 
      Y_test <- Y[foldid==ifold,]
      X_test <- X[foldid==ifold,]
      Z_test <- Z[foldid==ifold,] 
      Y_train <- Y[foldid!=ifold,]
      X_train <- X[foldid!=ifold,]
      Z_train <- Z[foldid!=ifold,]  
      mu_test <- matrix(NA, sum(foldid==ifold), q)
      for (il in 1:nlam) {
        if(il!=1 & warm==T) init.i =list(A0=fit1$A, V0=fit1$V, C0=fit1$C, kappa=fit1$kappa)
        fit1 <- msrrr.fit(Y_train, X_train, Z_train, family, familygroup, nrank, lamseq[il], init.i, control)
        eta_test = X_test%*%fit1$B + Z_test%*%fit1$C
        for(iq in 1:q){
          fm = fg[[iq]]
          mu_test[,iq] = fm$linkinv(eta_test[,iq])      
        }  
        if(cv.criteria=='deviance') cv[il, ifold] <- 2 * objFun(Y_test, mu_test, fit1$phi, family, familygroup)
        if(cv.criteria=='pMSE') cv[il, ifold] <- pmse(Y_test, mu_test, fit1$phi, family, familygroup)
      }
    }
    cv.mean <- apply(cv, 1, mean)
    lam.idx <- which.min(cv.mean)

    Tunepath = data.frame(cv, cv.mean = cv.mean)
    tunepath.opt = min(cv.mean) 
  }
  
  lam.opt = lamseq[lam.idx] 
  ## refit
  init.i = NULL
  if(lam.idx>1 & warm==T) init.i = list(A0=Apath[lam.idx,,], V0=Vpath[lam.idx,,], C0=Cpath[lam.idx,,])
  fit = msrrr.fit(Y, X, Z, family, familygroup, nrank, lam.opt, init.i, control)
     
  out <- list(lamseq=lamseq, nrank=nrank, Apath=Apath, Vpath=Vpath, Cpath=Cpath,phipath=phipath,mupath=mupath,  
              method=method, cv.criteria=cv.criteria, foldid=foldid, nfold=nfold, c.BIC=c.BIC, # c('CV', 'BIC', 'BICP', 'AIC', 'GIC'), 
              Tunepath = Tunepath, lam.idx = lam.idx, lam.opt = lam.opt, tunepath.opt = tunepath.opt,
              fit = fit, A=fit$A, V=fit$V, B=fit$B, C=fit$C)
  return(out)
}

 
## final wrapper: select nrank  
# Z not include a intercept col
msrrr = function(Y, X, Z=NULL, family, familygroup, nrankseq=c(1:3), init=NULL, 
                  lamseq=NULL, nlam=50, warm=T,   
                  method='CV', cv.criteria='pMSE', foldid=NULL, nfold=5, c.BIC=2, 
                  # c('CV', 'BIC', 'BICP', 'AIC', 'GIC')
                  control=list(epsilon=1e-4, maxit=200, trace=F, conv.obj=T)){
  n = nrow(X)
  Z = cbind(rep(1,n), Z) # 
  # pz = ncol(Z)
  nr = length(nrankseq)
  out.allrank = list()
  
  if(method=='CV' & is.null(foldid)) {   
    foldid <- rep(1:nfold, len = n)
    foldid <- sample(foldid, n, replace = FALSE)
  }  
  for(ir in 1:nr){
    if(control$trace==T) cat(ir)
    out.allrank[[ir]] = msrrr.tuning(Y, X, Z, family, familygroup, nrankseq[ir], init, lamseq, 
                                     nlam, warm, method, cv.criteria, foldid, nfold, c.BIC, control)
  }
  names(out.allrank) = paste0('nrank_', nrankseq)
  Tunepath = lapply(out.allrank, function(a) a$Tunepath)
  tunepath.opt = lapply(out.allrank, function(a) a$tunepath.opt)
  idx = which.min(tunepath.opt)
  out = out.allrank[[idx]]
  
  out$nrankseq = nrankseq
  out$out.allrank = out.allrank
  out$Tunepath = Tunepath
  out$tunepath.opt = tunepath.opt
  out$nrank.opt = nrankseq[idx]
  out$family = family
  out$familygroup = familygroup
  return(out)
}

 

predict.msrrr <- function(object, X.new, Z.new, Y.new, family=NULL, familygroup=NULL, type='response', cv.criteria="pMSE"){
  q = ncol(object$B)
  n = nrow(X.new)
  if(is.null(family)) family = object$family
  if(is.null(familygroup)) familygroup = object$familygroup
  fg = family[familygroup]
  eta.new = X.new%*%object$B + cbind(1,Z.new)%*%object$C
  mu.new = matrix(NA, n, q)

  for(iq in 1:q){
    fm = fg[[iq]]
    mu.new[,iq] = fm$linkinv(eta.new[,iq])
    # if(type=='response') 
  }
  
  # evaluate predictive performance (pMSE) if Ynew is given
  pred.perf = NULL
  if(!is.null(Y.new)) pred.perf <- pmse(Y.new, mu.new, object$fit$phi, family, familygroup)
  if(cv.criteria=='deviance') pred.perf <- 2 * objFun(Y.new, mu.new, object$fit$phi, family, familygroup)

  return(list(fit=mu.new, pred.perf=pred.perf))
}
