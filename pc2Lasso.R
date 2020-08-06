
# algorithm definitions

# helper function for pcLasso

pcLassof <- function(x, y, w, mg, aa, ulam, theta, ngroups, thr = 1e-4,
                     maxit = 1e5, family = "gaussian", verbose = FALSE) {
  no = nrow(x)
  ni = ncol(x)
  ne = ni
  nx = ni
  ng = ngroups
  if (is.null(w)) w = rep(1,no)
  nlam=length(ulam)
  verbose=1*verbose
  mode(no)="integer"
  mode(ni)="integer"
  mode(x)="double"
  mode(y)="double"
  mode(w)="double"
  mode(theta)="double"
  mode(ng)="integer"
  mode(mg)="integer"
  mode(aa)="double"
  mode(ne)="integer"
  mode(nx)="integer"
  mode(nlam)="integer"
  mode(ulam)="double"
  mode(thr)="double"
  mode(maxit)="integer"
  mode(verbose)="integer"
  if (family == "gaussian") {
    out = .Fortran("pclasso",
                   no,ni,x,y,w,theta,ng,mg,aa,ne,nx,nlam=nlam,ulam=ulam,thr,maxit,verbose,
                   ao=double(nx*nlam),
                   ia=integer(nx),
                   kin=integer(nlam),
                   nlp=integer(1),
                   jerr=integer(1),
                   PACKAGE="pcLasso")
  }
  if (family == "binomial") {
    out = .Fortran("logpclasso",
                   no,ni,x,y,w,theta,ng,mg,aa,ne,nx,nlam=nlam,ulam=ulam,thr,maxit,verbose,
                   a0=double(nlam),
                   ao=double(nx*nlam),
                   ia=integer(nx),
                   kin=integer(nlam),
                   nlp=integer(1),
                   jerr=integer(1),
                   PACKAGE="pcLasso")
  }
  
  ao = matrix(out$ao, nrow = ni)
  
  # uncompress soln
  beta <- matrix(0, ni, out$nlam)
  for (klam in 1:out$nlam) {
    temp <- out$kin[klam]
    beta[out$ia[1:temp], klam] <- ao[1:temp, klam]
  }
  
  a0 <- NA
  if (family == "binomial") a0 <- out$a0
  
  return(list(beta=beta, a0=a0, ulam=out$ulam, nlam=nlam, nlp=out$nlp,
              jerr=out$jerr))
}

msefun <- function(yhat,y) {
  (y - yhat)^2
}

binfun <- function(yhat, y) {
  - y * log(yhat) - (1 - y) * log(1 - yhat)
}

error.bars <- function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}


print.pcLasso=function(x,digits = max(3, getOption("digits") - 3),...){
  devratio=(x$dev[1]-x$dev)/x$dev[1]
  
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Nonzero = x$nzero, `%Dev` = signif(devratio, digits), 
              Lambda = signif(x$lambda, digits)))
  
  
}


pc2Lasso <-  function (x, y, w = rep(1, length(y)), family = c("gaussian", 
                                                  "binomial"), ratio = NULL, theta = NULL, groups = vector("list", 
                                                                                                           1), lambda.min.ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04), 
          nlam = 100, lambda = NULL, standardize = F, SVD_info = NULL, 
          nv = NULL, propack = T, thr = 1e-04, maxit = 1e+05, verbose = FALSE) 
{
  this.call <- match.call()
  n <- nrow(x)
  p <- ncol(x)
  y <- as.vector(y)
  if (length(y) != n) {
    stop("length of y is not equal to number of rows of x")
  }
  family <- match.arg(family)
  if (family == "binomial" && any(!(unique(y) %in% c(0, 1)))) {
    stop("if family is binomial, y can only contain 0s and 1s")
  }
  if (length(groups) == 1) 
    groups[[1]] <- 1:p
  if (length(unique(unlist(groups))) < p) {
    stop("Some features not assigned to a group")
  }
  ngroups <- length(groups)
  sizes <- unlist(lapply(groups, length))
  overlap <- F
  origx <- x
  origmx <- colMeans(x)
  origp <- ncol(x)
  origgroups <- groups
  if (length(origgroups) > 1) {
    nc <- length(unlist(origgroups))
    if (nc > p) {
      groups <- vector("list", length(origgroups))
      overlap <- T
      x <- matrix(NA, n, nc)
      i1 <- 1
      for (k in 1:ngroups) {
        i2 <- i1 + length(origgroups[[k]]) - 1
        x[, i1:i2] <- origx[, origgroups[[k]]]
        groups[[k]] <- i1:i2
        i1 <- i2 + 1
      }
    }
    p <- ncol(x)
  }
  mx <- colMeans(x)
  x <- scale(x, mx, F)
  if (standardize) {
    x <- scale(x, T, T)
  }
  my <- NA
  if (family == "gaussian") {
    my <- mean(y)
    y <- y - my
  }
  if (is.null(SVD_info)) {
    v = d = dd = vector("list", ngroups)
    if (verbose) 
      cat("Starting SVD computation", fill = T)
    for (k in 1:ngroups) {
      nvv <- nv
      if (is.null(nv)) {
        nvv <- min(nrow(x[, groups[[k]], drop = F]), 
                   ncol(x[, groups[[k]], drop = F]))
      }
      if (nrow(x) > length(groups[[k]])) {
        print("Execution of code with the last modification")
        xtemp <- t(x[, groups[[k]]]) %*% x[, groups[[k]]]             
        eig <- eigen(xtemp)                                           
        v[[k]] <- eig$vec
        d[[k]] <- eig$val
        eigensum <- sum(d[[k]])                             
        len <- length(d[[k]])
        var_exp_1 <- d[[k]][1]/eigensum
        var_exp_2 <- d[[k]][2]/eigensum
        second <- sort(d[[k]],partial=len-1)[len-1]
        dd[[k]] <- var_exp_1*max(d[[k]]) + var_exp_2*second - 2*d[[k]] + 0.01
        # dd[[k]] <- max(d[[k]]) + second - 2*d[[k]] + 0.01
      }
      else {
        if (propack) {
          sv <- svd::propack.svd(x[, groups[[k]], drop = F], 
                                 neig = nvv)
        }
        else {
          sv <- svd(x[, groups[[k]], drop = F], nv = nvv)
        }
        # print("Execution of code with the last modification 2")
        v[[k]] <- sv$v                                                               
        d[[k]] <- sv$d^2                                                             
        eigensum <- sum(d[[k]])
        len <- length(d[[k]])
        var_exp_1 <- max(d[[k]])/eigensum
        second <- sort(d[[k]],partial=len-1)[len-1]
        var_exp_2 <- second/eigensum
        dd[[k]] <- var_exp_1*max(d[[k]]) + var_exp_2*second - 2*d[[k]] + 0.01                            # modified
        # dd[[k]] <- max(d[[k]]) + second - 2*d[[k]] + 0.01
      }
    }
    if (verbose) 
      cat("SVD completed", fill = T)
  }
  else {
    aa <- SVD_info$aa
    d <- SVD_info$d
    dd <- SVD_info$dd
  }
  if (missing(ratio) && missing(theta)) {
    stop("Provide ratio or theta")
  }
  else if (!missing(ratio) && !missing(theta) && !is.null(ratio) && 
           !is.null(theta)) {
    stop("Provide only ratio or theta, not both")
  }
  else if (is.null(theta)) {
    if (ratio < 0 || ratio > 1) {
      stop("ratio must be in [0, 1]")
    }
    thetanew <- rep(NA, ngroups)
    for (k in 1:ngroups) {
      thetanew[k] <- d[[k]][2] * (1 - ratio)/(ratio * 
                                                (d[[k]][1] - d[[k]][2]))
    }
    theta <- mean(thetanew)
  }
  else if (theta < 0) {
    stop("theta must be non-negative")
  }
  if (is.null(SVD_info)) {
    aa <- matrix(0, p, max(sizes))
    i1 <- 1
    for (k in 1:ngroups) {
      i2 <- i1 + sizes[k] - 1
      aa[i1:i2, 1:sizes[k]] <- scale(v[[k]], FALSE, 1/(dd[[k]])) %*% 
        t(v[[k]])
      i1 <- i2 + 1
    }
    SVD_info <- list()
    SVD_info$aa <- aa
    SVD_info$d <- d
    SVD_info$dd <- dd
  }
  i1 <- 1
  mg <- c(1, cumsum(sizes) + 1)
  ulam <- lambda
  if (is.null(lambda)) {
    maxlam <- max(abs(t(x) %*% (y - mean(y))))
    if (family == "binomial") {
      maxlam <- 4 * maxlam
    }
    ulam <- exp(seq(log(maxlam), log(maxlam * lambda.min.ratio), 
                    length = nlam))
  }
  nlam <- length(ulam)
  out <- pcLassof(x, y, w, mg, aa, ulam, theta, ngroups, thr = thr, 
                  maxit = maxit, family = family, verbose = verbose)
  nzero <- colSums(out$beta != 0)
  if (family == "gaussian") 
    a0 <- rep(my, nlam)
  if (family == "binomial") 
    a0 <- out$a0
  if (standardize) {
    out$beta <- out$beta * matrix(attr(x, "scaled:scale"), 
                                  nrow = nrow(out$beta), ncol = ncol(out$beta))
  }
  origbeta <- NULL
  orignzero <- NULL
  if (overlap) {
    origbeta <- matrix(0, ncol(origx), ncol(out$beta))
    for (k in 1:ngroups) {
      origbeta[origgroups[[k]], ] <- origbeta[origgroups[[k]], 
                                              ] + out$beta[groups[[k]], ]
    }
    orignzero <- colSums(origbeta != 0)
  }
  out <- list(beta = out$beta, origbeta = origbeta, a0 = a0, 
              lambda = out$ulam, nzero = nzero, orignzero = orignzero, 
              jerr = out$jerr, theta = theta, origgroups = origgroups, 
              groups = groups, SVD_info = SVD_info, mx = mx, origmx = origmx, 
              my = my, overlap = overlap, nlp = out$nlp, family = family, 
              call = this.call)
  yhat <- predict.pcLasso(out, origx)
  if (family == "gaussian") 
    dev <- colSums(apply(yhat, 2, msefun, y))
  if (family == "binomial") 
    dev <- colSums(apply(yhat, 2, binfun, y))
  out$dev <- dev
  class(out) <- "pcLasso"
  return(out)
}



cv.pc2Lasso <- function (x, y, w = rep(1, length(y)), ratio = NULL, theta = NULL, 
          groups = vector("list", 1), family = "gaussian", nfolds = 10, 
          foldid = NULL, keep = FALSE, verbose = FALSE, ...) 
{
  this.call <- match.call()
  n <- nrow(x)
  p <- ncol(x)
  if (length(groups) == 1) 
    groups[[1]] <- 1:p
  ngroups <- length(groups)
  if (!missing(foldid)) {
    nfolds <- length(unique(foldid))
  }
  else {
    foldid <- sample(rep(seq(nfolds), length = length(y)))
  }
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  fit0 <- pc2Lasso(x, y, groups = groups, ratio = ratio, theta = theta, 
                  family = family, ...)
  cat("Initial fit done- including SVD", fill = T)
  fits <- vector("list", nfolds)
  for (ii in 1:nfolds) {
    cat(c("Fold=", ii), fill = T)
    oo <- foldid == ii
    xc <- x[!oo, , drop = F]
    yy <- y[!oo]
    fits[[ii]] <- pc2Lasso(xc, yy, SVD_info = fit0$SVD_info, 
                          groups = groups, theta = fit0$theta, family = family, 
                          lambda = fit0$lambda, ...)
  }
  yhat <- matrix(NA, n, length(fit0$lambda))
  for (ii in 1:nfolds) {
    oo <- foldid == ii
    out <- predict(fits[[ii]], x[oo, , drop = F])
    yhat[oo, 1:ncol(out)] <- out
  }
  if (family == "binomial") {
    yhat <- 1/(1 + exp(-yhat))
  }
  if (family == "gaussian") {
    errfun = msefun
    name = "Mean-Squared Error"
  }
  if (family == "binomial") {
    errfun = binfun
    name = "Deviance"
  }
  ym <- array(y, dim(yhat))
  err <- errfun(yhat, ym)
  cvm <- apply(err, 2, mean, na.rm = T)
  nn <- apply(!is.na(err), 2, sum, na.rm = T)
  cvse <- sqrt(apply(err, 2, var, na.rm = T)/nn)
  cvlo <- cvm - cvse
  cvup <- cvm + cvse
  yhat.preval <- NULL
  foldid_copy <- NULL
  if (keep) {
    yhat.preval <- yhat
    foldid_copy <- foldid
  }
  imin <- which.min(cvm)
  lambda.min <- fit0$lambda[imin]
  imin.1se <- which(cvm < cvm[imin] + cvse[imin])[1]
  lambda.1se <- fit0$lambda[imin.1se]
  obj <- list(glmfit = fit0, theta = fit0$theta, lambda = fit0$lambda, 
              nzero = fit0$nzero, orignzero = fit0$orignzero, fit.preval = yhat, 
              cvm = cvm, cvse = cvse, cvlo = cvlo, cvup = cvup, lambda.min = lambda.min, 
              lambda.1se = lambda.1se, foldid = foldid_copy, name = name, 
              call = this.call)
  class(obj) <- "cv.pcLasso"
  return(obj)
}

set.seed(1)

# data loading and computations