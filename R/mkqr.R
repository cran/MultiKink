
#' Test the existence of kink effect in the multi-kink quantile regression
#'
#' This function tests the existence of a kink effect in the multi-kink quantile regression.
#'
#' @param y A numeric vector of response.
#' @param thre.x A numeric vector of scalar covariate with threshold effect.
#' @param cont.z A numeric matrix of design covariates with constant slopes.
#' @param id A numeric vector of index used for longitudinal data; can be missing or NULL for iid data.
#' @param tau A numeric scalar representing the quantile level (default is 0.5).
#' @param NB An integer specifying the resampling times (default is 200).
#' @param sparsity The error term type. Specify one from "iid" and "nid" (default is "nid").
#' @param bandwidth_type The bandwidth type. Specify one from "Hall-Sheather", "Bofinger", or "Chamberlain" (default is "Hall-Sheather").
#'
#' @return A list containing the p-value (pv), the statistic based on the original data (Tn), and the statistics by wild bootstrap (Tn.NB).
#'
#' @examples
#' # Example 1: i.i.d data type
#' library(quantreg)
#' n = 200
#' Z1 <- rexp(n,1)
#' Z2 <- rbinom(n,1,0.5)
#' Z <- cbind(Z1,Z2)
#' epsilon <- rnorm(n,0,1)
#' X <- runif(n,-2,1)
#' psi <- c(-1,0)
#' k <- length(psi)
#' PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
#' XP <- matrix(rep(X,k),nrow=n)
#' XR <- cbind(1,X,pmax((XP-PSI),0),Z)
#' bet <- c(1,-1,0,0,sqrt(3),-sqrt(3))
#' Y <- XR %*% bet + epsilon
#' obj <- KinkTest(y=Y,thre.x=X,cont.z=Z,
#'                 bandwidth_type=c("Hall-Sheather"))
#' obj$pv
#'
#' \dontrun{
#' # Example 2: longitudinal data
#' library(quantreg)
#' N = 200
#' T = 5
#' subject = rep(1:N,each=T)
#' NT = N*T
#' Z1 <- rexp(NT,1)
#' Z2 <- rbinom(NT,1,0.5)
#' Z <- cbind(Z1,Z2)
#' epsilon <- rnorm(NT,0,1)
#' X <- runif(NT,-2,1)
#' psi <- c(-1,0)
#' k <- length(psi)
#' PSI <- matrix(rep(psi,rep(NT,k)),ncol=k)
#' a <- rnorm(N,0,1)
#' A <- rep(a,each=T)
#' XP <- matrix(rep(X,k),nrow=NT)
#' XR <- cbind(1,X,pmax((XP-PSI),0),Z)
#' bet <- c(1,-1,0,0,sqrt(3),-sqrt(3))
#' Y <- XR %*% bet + A + epsilon
#' obj <- KinkTest(y=Y,thre.x=X,cont.z=Z,id=subject,
#'                 bandwidth_type=c("Hall-Sheather"))
#' obj$pv
#' }
#' @export

KinkTest <- function(y,thre.x,cont.z,id,tau=0.5,NB = 200, sparsity = "nid",
                     bandwidth_type=c("Hall-Sheather","Bofinger","Chamberlain"))
{
    # require(MASS)
    # require(quantreg)
    # require(pracma)
    X <- cbind(1,thre.x,cont.z)
    n <- length(y)
    wfun <- function(u, tau) tau-(u<0)

    density_fun <- function(y,X,tau,bandwidth_type){
        eps <- .Machine$double.eps^(2/3)

        if (bandwidth_type == "Bofinger") {
            bandwidth <- n^(-1/5) * ((9/2 * dnorm(qnorm(tau))^4)/(2 * qnorm(tau)^2 + 1)^2)^(1/5)
        } else if (bandwidth_type == "Chamberlain") {
            alpha <- 0.05
            bandwidth <- qnorm(1 - tau/2) * sqrt(alpha * (1 - alpha)/n)
        } else if (bandwidth_type == "Hall-Sheather") {
            alpha <- 0.05
            bandwidth <- n^(-1/3) * qnorm(1 - alpha/2)^(2/3) *
                ((3/2 * dnorm(qnorm(tau))^2)/(2 * qnorm(tau)^2 + 1))^(1/3)
        } else {
            stop("Not a valid bandwith method!")
        }
        # Compute the density
        # Hendricks and Koenker (1992)
        bu <- suppressWarnings(fitted.values(rq(y~X[,-1],tau+bandwidth,method="br")))
        bl <- suppressWarnings(fitted.values(rq(y~X[,-1],tau-bandwidth,method="br")))
        density <- pmax(0, 2 * bandwidth/((bu - bl) - eps))

        density
    }
    if(missing(id) || is.null(id)){
        id <- seq(1:n)
    }

    sub = model.matrix(~-1+factor(id))
    N = dim(sub)[2]
    nvec <- table(id)

    ##estimate under H0
    fit <- quantreg::rq(y~X[,-1],tau,method="br")
    res <- fit$residuals
    wt <- wfun(res, tau)

    testFun <- function(tt){
        #fit <- rq(y~X[,-1],tau,method="br")
        #res <- fit$residuals
        #wt <- wfun(res, tau)
        Rn <-  rep(NA, length(tt))
        for (kk in 1:length(tt)){

            Rn[kk] <- 1/sqrt(n)*sum(wt*(thre.x-tt[kk])*(thre.x <= tt[kk]))
        }
        Tn <- max(abs(Rn))
        return(Tn)

    }

    ## perturbed method to calculate the p-value
    testFun.resample <- function(tt){


        #########################
        ## permutation random errors
        #n<- length(y)
        u <- rnorm(n, 0, 1)-qnorm(tau,0,1)
        B <- rbinom(n,1,0.5)
        w <- (B==1)*1+(B==0)*(-1)


        #########################

        #fit <- quantreg::rq(y~X[,-1],tau,method="br")
        #res <- fit$residuals
        wu <- wfun(u, tau)

        if (sparsity == "iid"){

            Sn <- (t(X) %*% X)/n
            Rn <- rep(NA,length(tt))

            for(kk in 1:length(tt)){

                Sn.t <- apply(X*((thre.x-tt[kk])*(thre.x <= tt[kk])), 2, mean)
                Rn[kk] <- 1/sqrt(n)*sum(
                    w*wu * ((thre.x-tt[kk])*(thre.x <= tt[kk]) -
                                X %*% solve(Sn) %*% Sn.t)
                )
            }

        } else if (sparsity == "nid") {

            ker <- density_fun(y,X,tau,bandwidth_type)
            #hn <- 1.06* n^(-1/5) * sd(res)
            #ker <- Ku(res/hn)/hn#  %*% diag(ker)
            #h <- 1.06* n^(-1/5)* sd(res.rank)
            #pdf<- approxfun(density(res, kernel= "epanechnikov", bw=hn))
            #ker <- pdf(res)


            Sn <- (t(X)%*%diag(ker) %*%  X)/n

            ## under H0
            Rn <-  rep(NA, length(tt))

            for (kk in 1:length(tt)){#

                Sn.t <- apply(X*(ker*(thre.x-tt[kk])*(thre.x <= tt[kk])), 2, mean)
                Rn[kk] <- 1/sqrt(n)*sum(
                    w*wu * ((thre.x-tt[kk])*(thre.x <= tt[kk]) -
                                X %*% solve(Sn) %*% Sn.t)
                )
            }
        } else {
            stop("Not a valid error type !")
        }


        Tn <- max(abs(Rn))

        return(Tn)
    }


    ## longitudinal version
    testFun.resample.long <- function(tt){
        #########################
        ## permutation random errors
        #n<- length(y)

        u <- rnorm(N, 0, 1)
        uu <- sapply(1:N,function(j){rep(u[j],nvec[j])},simplify = FALSE )
        uu <- unlist(uu)

        ker <- density_fun(y,X,tau,bandwidth_type )
        Sn <- (t(X)%*%diag(ker) %*%  X)/n

        ## under H0
        Rn <-  rep(NA, length(tt))

        for (kk in 1:length(tt)){#

            Sn.t <- apply(X*(ker*(thre.x-tt[kk])*(thre.x <= tt[kk])), 2, mean)
            Rn[kk] <- 1/sqrt(n)*sum(
                uu*wt * ((thre.x-tt[kk])*(thre.x <= tt[kk]) -
                             X %*% solve(Sn) %*% Sn.t)
            )
        }

        Tn <- max(abs(Rn))

        return(Tn)
    }

    #######################################################
    ###  calculate the p-value by wild bootstrap

    tt <- seq(quantile(thre.x,0.1), quantile(thre.x,0.9), length = 100)

    Tn <-  testFun(tt)
    if(N == n){  #cross-section data
        cat("The score test is performed for iid data type \n")
        Tn.NB <- replicate(NB, testFun.resample(tt))
    }else if(n > N){
        cat("The score test is performed for longitudinal data type \n")
        Tn.NB <- replicate(NB, testFun.resample.long(tt))
    }

    pv <- mean(Tn.NB > Tn,  na.rm = TRUE)
    return(list(pv=pv,Tn=Tn,Tn.NB=Tn.NB))


}


#' Auxiliary parameters to control the model fitting
#'
#' This function defines auxiliary parameters that control the model fitting process.
#'
#' @param toll Positive convergence tolerance.
#' @param h Positive factor (from zero to one) modifying the increments in kink parameter updates during the iterative process.
#' @param it.max Positive integer for the maximal number of iterations.
#' @param K.max Positive integer for the maximal given number of kink points.
#' @param stop.if.error Logical indicating if the estimation algorithm should be stopped if some kink point estimators belong to the non-admissible set. Default is FALSE which suggests removing the non-admissible change points automatically.
#' @param dev0 Initial objective value or deviance. Default is NULL which implies that the initial value is unknown.
#' @param visual Logical indicating if the results of the estimation process should be printed at each iteration.
#' @param visualBoot Logical indicating if the results of estimation should be printed at each iteration in the bootstrap restarting process.
#' @param pow The powers of the pseudo covariates employed by the algorithm.
#' @param digits If specified, it means the desired number of decimal points of the kink estimators to be used during the iterative algorithm.
#' @param grid It measures how close between the two adjacent change points should be merged, default is NULL.
#' @param n.boot Positive integer indicating the times of bootstrap re-sampling in the bootstrap restarting algorithm, default is 20.
#'
#' @return A list with the arguments as components to be used by \emph{mkqr.fit} and \emph{mkqr.bea}.
#'
#' @examples
#' # Example usage
#' fit.control(K.max=8)
#'
#' @export
fit.control <-
    function(toll=1e-4,h=1,it.max=50,K.max=6,stop.if.error=TRUE,dev0=NULL,visual=FALSE,
             visualBoot=FALSE,pow=c(1,1),digits=NULL,grid=NULL,n.boot=20){
        list(toll=toll,h=h,it.max=it.max,K.max=K.max,stop.if.error=stop.if.error,
             dev0=dev0,visual=visual,n.boot=n.boot,
             visualBoot=visualBoot,pow=pow,digits=digits,grid=grid)
    }


brisq <-function(y, XREG,X, PSI, tau, opz, n.boot=20,
                 size.boot=NULL, jt=FALSE,
                 nonParam=TRUE, random=FALSE)
{
    extract.psi<-function(lista){
        dev.values<-lista[[1]]
        psi.values<-lista[[2]]
        dev.values <- as.vector(dev.values, mode = "numeric")
        dev.ok <- min(dev.values,na.rm = TRUE)
        id.dev.ok<-which.min(dev.values)
        if(is.list(psi.values))  psi.values<-matrix(unlist(psi.values),
                                                    nrow=length(dev.values), byrow=TRUE)
        if(!is.matrix(psi.values)) psi.values<-matrix(psi.values)
        psi.ok<-psi.values[id.dev.ok,]
        r<-list(SumSquares.no.gap=dev.ok, psi=psi.ok)
        r
    }
    #-------------
    visualBoot<- opz$visualBoot
    opz.boot<-opz
    opz.boot$pow=c(1.1,1.2)
    opz1<-opz
    opz1$it.max <-1
    n<-length(y)
    #x <- X[,1]


    #ris <- seg.qr.fit(y,XREG,X,Z,PSI,tau=0.5,opz=seg.control())
    o0 <- try(brisq.fit(y, XREG, X, PSI, tau, opz, return.all.sol=FALSE), silent=TRUE)

    rangeZ <- apply(X, 2, range)
    if(!is.list(o0)) {
        o0 <- brisq.fit(y, XREG, X, PSI, tau, opz, return.all.sol=TRUE)
        o0<-extract.psi(o0)
        if(!nonParam) {warning("using nonparametric boot");nonParam<-TRUE}
    }
    if(is.list(o0)){
        est.psi00<-est.psi0<-o0$psi
        ss00<-o0$SumSquares.no.gap
        if(!nonParam) fitted.ok<-fitted(o0$obj)
    } else {
        if(!nonParam) stop("the first fit failed and I cannot extract fitted values for the boot sample")
        if(random) {
            est.psi00<-est.psi0<-apply(rangeZ,2,function(r)runif(1,r[1],r[2]))
            PSI1 <- matrix(rep(est.psi0, rep(nrow(X), length(est.psi0))), ncol = length(est.psi0))
            o0<-try(brisq.fit(y, XREG,X, PSI, tau, opz1), silent=TRUE)
            ss00<-o0$SumSquares.no.gap
        } else {
            est.psi00<-est.psi0<-apply(PSI,2,mean)
            ss00<-opz$dev0
        }
    }

    all.est.psi.boot<-all.selected.psi<-all.est.psi<-matrix(, nrow=n.boot, ncol=length(est.psi0))
    all.ss<-all.selected.ss<-rep(NA, n.boot)
    if(is.null(size.boot)) size.boot<-n

    #      na<- ,,apply(...,2,function(x)mean(is.na(x)))

    X.orig<- X
    if(visualBoot) cat(0, " ", formatC(opz$dev0, 3, format = "f"),"", "(No breakpoint(s))", "\n")
    count.random<-0
    for(k in seq(n.boot)){
        PSI <- matrix(rep(est.psi0, rep(nrow(X), length(est.psi0))), ncol = length(est.psi0))
        if(jt)  X<-apply(X.orig,2,jitter)
        if(nonParam){
            id<-sample(n, size=size.boot, replace=TRUE)
            o.boot<-try(brisq.fit(y[id], XREG[id,,drop=FALSE],X[id,,drop=FALSE], PSI[id,,drop=FALSE],
                                  tau=tau, opz.boot), silent=TRUE)
        } else {
            yy<-fitted.ok+sample(residuals(o0$obj),size=n, replace=TRUE)#residuals(o0)
            o.boot<-try(brisq.fit(yy, XREG, X.orig, PSI, tau=tau, opz.boot), silent=TRUE)
        }
        if(is.list(o.boot)){
            all.est.psi.boot[k,]<-est.psi.boot<-o.boot$psi
        } else {
            est.psi.boot<-apply(rangeZ,2,function(r)runif(1,r[1],r[2]))
        }
        PSI <- matrix(rep(est.psi.boot, rep(nrow(X), length(est.psi.boot))), ncol = length(est.psi.boot))
        opz$h<-max(opz$h*.9, .2)
        opz$it.max<-opz$it.max+1
        o<-try(brisq.fit(y, XREG,X.orig,  PSI, tau=tau, opz, return.all.sol=TRUE), silent=TRUE)
        if(!is.list(o) && random){
            est.psi0<-apply(rangeZ,2,function(r)runif(1,r[1],r[2]))
            PSI1 <- matrix(rep(est.psi0, rep(nrow(X), length(est.psi0))), ncol = length(est.psi0))
            o<-try(brisq.fit(y, XREG,X,  PSI1, tau, opz1), silent=TRUE)
            count.random<-count.random+1
        }
        if(is.list(o)){
            if(!"coefficients"%in%names(o$obj)) o<-extract.psi(o)
            all.est.psi[k,]<-o$psi
            all.ss[k]<-o$SumSquares.no.gap
            if(o$SumSquares.no.gap<=ifelse(is.list(o0), o0$SumSquares.no.gap, 10^12)) o0<-o
            est.psi0<-o0$psi
            all.selected.psi[k,] <- est.psi0
            all.selected.ss[k]<-o0$SumSquares.no.gap #min(c(o$SumSquares.no.gap, o0$SumSquares.no.gap))
        }
        if(visualBoot) {
            flush.console()
            spp <- if (k < 10) "" else NULL
            cat(k, spp, "", formatC(o0$SumSquares.no.gap, 3, format = "f"), "\n")
        }
    } #end n.boot
    #browser()


    all.selected.psi<-rbind(est.psi00,all.selected.psi)
    all.selected.ss<-c(ss00, all.selected.ss)

    SS.ok<-min(all.selected.ss)
    id.accept<- ((abs(all.ss-SS.ok)/SS.ok )<= 0.05)
    psi.mean<-apply(all.est.psi[id.accept,,drop=FALSE], 2, mean)

    #      est.psi0<-psi.mean
    #      #devi ristimare il modello con psi.mean
    #      PSI1 <- matrix(rep(est.psi0, rep(nrow(Z), length(est.psi0))), ncol = length(est.psi0))
    #      o0<-try(seg.lm.fit(y, XREG, Z, PSI1, w, offs, opz1), silent=TRUE)

    ris<-list(all.selected.psi=drop(all.selected.psi),all.selected.ss=all.selected.ss,
              all.psi=all.est.psi, all.ss=all.ss)

    if(is.null(o0$obj)){
        PSI1 <- matrix(rep(est.psi0, rep(nrow(X), length(est.psi0))), ncol = length(est.psi0))
        o0<-try(brisq.fit(y, XREG, X, PSI1, tau, opz1), silent=TRUE)
    }
    if(!is.list(o0)) return(0)
    o0$boot.restart<-ris
    return(o0)
}


brisq.fit <-function(y,XREG,X,PSI,tau,opz,return.all.sol=FALSE)
{

    #-----------
    psi <- PSI[1,]
    n <- length(y)
    x <- X[,1]
    c1 <- apply((X <= PSI), 2, all)
    c2 <- apply((X >= PSI), 2, all)
    if(sum(c1 + c2) != 0 || is.na(sum(c1 + c2))) stop("psi out of the range")
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ##~~~~~~~~~~~~~
    xreg.names <- opz$xreg.names
    #mtype <- opz$mtype
    grid <- opz$grid
    digits <- opz$digits
    pow<-opz$pow
    #nomiOK<-opz$nomiOK   ##
    toll<-opz$toll
    h<-opz$h ##
    #gap<-opz$gap  ##FALSE
    stop.if.error<-opz$stop.if.error
    dev.new<-opz$dev0
    visual<-opz$visual
    #id.psi.group<-opz$id.psi.group  #
    it.max<-old.it.max<-opz$it.max
    rangeZ <- apply(X, 2, range)
    #psi<-PSI[1,]
    #names(psi)<-id.psi.group
    #H<-1
    it <- 1
    epsilon <- 10
    dev.values<-psi.values <- NULL
    id.psi.ok<-rep(TRUE, length(psi))
    #abs(epsilon) > toll
    while (it < it.max) {
        U <- pmax((X - PSI), 0)
        V <- ifelse((X > PSI), -1, 0)

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        XZ <- cbind(XREG,U,V)
        rownames(XZ) <- NULL
        colnames(XZ) <- c(xreg.names,paste("U",1:ncol(U),sep=""),paste("V",1:ncol(V),sep=""))
        ####~~~~~~~~~~~~~~~~~~~``

        obj <- rq(y~XZ, tau = tau,  method ="br")
        dev.old<-dev.new
        dev.new <- dev.new1 <- obj$rho
        dev.values[[length(dev.values) + 1]] <- dev.new1
        if (visual) {
            flush.console()
            if (it == 1)
                cat(0, " ", formatC(dev.old, 3, format = "f"),
                    "", "(No breakpoint(s))", "\n")
            spp <- if (it < 10) "" else NULL
            cat(it, spp, "", formatC(dev.new, 3, format = "f"), "",length(psi),"\n")
            #cat(paste("iter = ", it, spp," dev = ",formatC(dev.new,digits=3,format="f"), " n.psi = ",formatC(length(psi),digits=0,format="f"), sep=""), "\n")
        }
        epsilon <- (dev.new - dev.old)/(dev.old + .001)
        obj$epsilon <- epsilon
        it <- it + 1
        obj$it <- it

        beta.c <- coef(obj)[paste("XZU",1:ncol(U),sep="")]
        gamma.c <-  coef(obj)[paste("XZV", 1:ncol(V), sep = "")]

        #if (it > it.max) break
        if(it>10 && epsilon<toll) break
        #if(max(abs(gamma.c))<1e-4) break

        psi.values[[length(psi.values) + 1]] <- psi.old <- psi
        #       if(it>=old.it.max && h<1) H<-h
        psi <- psi.old + h*gamma.c/beta.c
        if(!is.null(digits)) psi<-round(psi, digits)
        PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
        #check if psi is admissible..
        a <- apply((X <= PSI), 2, all) #prima era solo <
        b <- apply((X >= PSI), 2, all) #prima era solo >
        if(stop.if.error) {
            isErr<- (sum(a + b) != 0 || is.na(sum(a + b))) || (any(diff(sort(psi))<grid))
            if(isErr) {
                if(return.all.sol) return(list(dev.values, psi.values)) else stop("(Some) estimated psi gets wrong")
            }
        } else {
            id.psi.ok<-!is.na((a+b)<=0)&(a+b)<=0
            #X <- X[,id.psi.ok,drop=FALSE]
            psi <- psi[id.psi.ok]
            PSI <- PSI[,id.psi.ok,drop=FALSE]
            X <- X[,id.psi.ok,drop=FALSE]
            #nomiOK<-nomiOK[id.psi.ok]
            #id.psi.group<-id.psi.group[id.psi.ok]
            #names(psi)<-id.psi.group

            ##~~~~~~~~~~~~~~~~~~~~~~~~~~
            if(any(diff(sort(psi))<grid)){
                psi <- sort(psi)
                PSI <- PSI[,order(psi)]
                X <- X[,order(psi)]
                id.psi.grid <- which(diff(psi)<grid)+1
                psi <- psi[-id.psi.grid]
                PSI <- matrix(PSI[,-id.psi.grid],nrow=n)
                X <- matrix(X[,-id.psi.grid],nrow=n)
                #nomiOK<-nomiOK[-id.psi.grid]
                #id.psi.group<-id.psi.group[-id.psi.grid]
                #names(psi)<-id.psi.group
            }
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~

            #if(length(psi)<=0) return(0)
            if(length(psi)==0){
                obj <- rq(y~XREG,tau=tau,method="br")
                obj$n.psi <- 0
                obj$psi <- NULL
                return(obj)
            }
        } #end else
        #obj$psi <- psi
    } #end while
    psi <- sort(psi)
    names(psi) <- paste("psi",1:length(psi),sep="")
    #psi<-unlist(tapply(psi, id.psi.group, sort))
    #names(psi)<-id.psi.group
    X <- matrix(rep(x,length(psi)),nrow=n)
    PSI <- matrix(rep(psi, rep(n, length(psi))), ncol = length(psi))
    U <- pmax((X - PSI), 0)
    #V <- ifelse((X > PSI), -1, 0)
    XZ <- cbind(XREG,U)
    rownames(XZ) <- NULL
    colnames(XZ) <- c(xreg.names,paste("U",1:ncol(U),sep=""))

    obj <- rq(y~XZ, tau = tau,  method ="br")
    SS.new <- obj$rho

    #fino a qua..
    obj<-list(obj=obj,psi=psi,psi.values=psi.values,n.psi=length(psi), SumSquares.no.gap=SS.new)
    return(obj)
}




Cnmat <- function(y,XREG,thre.x,ht,id,psi.est,bet.est,beta.est,tau,p,resid,error.type)
{

    ##main functions
    n <- length(y)
    q <- length(psi.est)+length(bet.est)

    g <- function(xreg,x,psi.est,beta.est){
        k <- length(psi.est)
        vec <- c(1,xreg,pmax((rep(x,k)-psi.est),0),(-1)*beta.est*ifelse(rep(x,k)>psi.est,1,0))
        return(rbind(vec))
    }

    odiag=function(A, at = 0)
    {
        if (is.matrix(A)) {
            y <- A[row(A) == col(A) - at]
            return(y)
        }
        else {
            len <- length(A)
            B <- matrix(0, nrow = len + abs(at), ncol = len + abs(at))
            B[row(B) == col(B) - at] <- A
            return(B)
        }
    }

    gamfit=function(xdist,xsum,y){
        return(gam::gam(y~s(xdist,4)+s(xsum,4),family=binomial))
    }

    design=model.matrix(~-1+factor(id))
    N=dim(design)[2]
    XREG.lst <- list()
    x.lst <- list()
    resid.lst <- list()
    for(i in 1:N){
        XREG.lst[[i]] <- as.matrix(subset(XREG,design[,i]==1))
        x.lst[[i]] <- subset(thre.x,design[,i]==1)
    }
    nn <- rep(0,N)
    for(l in 1:N){
        resid.lst[[l]] <- subset(resid,design[,l]==1)
        nn[l] <- length(resid.lst[[l]])
    }

    if(error.type == "Compound"){
        yy <- NULL
        nsum <- 0
        for(i in 1:N){
            for(j in (1:nn[i])[-nn[i]]){
                for(j.pr in (j+1):nn[i]){
                    yy <- c(yy,(resid.lst[[i]][j]<0)*(resid.lst[[i]][j.pr]<0))
                }
            }
            nsum <- nsum+nn[i]*(nn[i]-1)/2
        }
        pjoint <- nsum^(-1)*sum(yy)
        Jn.1 <- tau*(1-tau)/n*t(ht)%*%ht
        Jn.2 <- matrix(0,q,q)
        for(i in 1:N){
            for(j in 1:nn[i]){
                for(j.pr in (1:nn[i])[-j]){
                    Jn.2 <- Jn.2 + t(g(XREG.lst[[i]][j,],x.lst[[i]][j],psi.est,beta.est)) %*%
                        (g(XREG.lst[[i]][j.pr,],x.lst[[i]][j.pr],psi.est,beta.est))*
                        (pjoint-tau^2)
                }
            }
        }
        Jn <- Jn.1 + Jn.2/n
        if(min(eigen(Jn)$values)<1e-8){Jn=as.matrix(Matrix::nearPD(Jn)$mat)}

    } else if(error.type == "AR"){
        di=max(nn)
        if(di>1){
            pr.m=matrix(0,di,di)
            count=rep(0,di-1)
            for(i in 1:N)
                if(nn[i]>1){
                    for(j in 1:(nn[i]-1))
                        for(j.pr in (j+1):nn[i]){
                            pr.m[j,j.pr]=pr.m[j,j.pr]+(resid.lst[[i]][j]<0)*(resid.lst[[i]][j.pr]<0)
                            pr.m[j.pr,j]=pr.m[j.pr,j]+(resid.lst[[i]][j.pr]<0)*(resid.lst[[i]][j]<0)
                        }
                    count[1:(nn[i]-1)]=count[1:(nn[i]-1)]+2*seq(nn[i]-1,1)
                }
            delta=rep(0,di-1)
            for(hh in 1:(di-1)) delta[hh]=sum(odiag(pr.m,hh)+odiag(pr.m,-hh))/count[hh]
        }
        Jn.1 <- tau*(1-tau)/n*t(ht)%*%ht
        Jn.2 <- matrix(0,q,q)

        for(i in 1:N){
            for(j in 1:nn[i]){
                for(j.pr in (1:nn[i])[-j]){
                    Jn.2 <- Jn.2 + t(g(XREG.lst[[i]][j,],x.lst[[i]][j],psi.est,beta.est)) %*%
                        (g(XREG.lst[[i]][j.pr,],x.lst[[i]][j.pr],psi.est,beta.est))*
                        (delta[abs(j-j.pr)]-tau^2)
                }
            }
        }
        Jn <- Jn.1 + Jn.2/n
        if(min(eigen(Jn)$values)<1e-8){Jn=as.matrix(Matrix::nearPD(Jn)$mat)}
    } else if(error.type == "general"){

        Jn.1 <- tau*(1-tau)/n*t(ht)%*%ht
        Jn.2 <- matrix(0,q,q)
        for(i in 1:N){
            for(j in 1:nn[i]){
                for(j.pr in (1:nn[i])[-j]){
                    Jn.2=Jn.2+t(g(XREG.lst[[i]][j,],x.lst[[i]][j],psi.est,beta.est)) %*%
                        (g(XREG.lst[[i]][j.pr,],x.lst[[i]][j.pr],psi.est,beta.est))*
                        ((resid.lst[[i]][j]<0)*(resid.lst[[i]][j.pr]<0)-(tau*(resid.lst[[i]][j.pr]<0)+
                                                                             (resid.lst[[i]][j]<0)*tau)/2)
                }
            }
        }
        Jn <- Jn.1 + Jn.2/n
        if(min(eigen(Jn)$values)<1e-8){Jn=as.matrix(Matrix::nearPD(Jn)$mat)}
    }else{
        stop("Please set the correct error structure !")
    }

    return(Jn)
}



fitting.k.GEE <- function(y,XREG,X,nk,tau,k,p,type,betaint)
{
    cn <- c(0,cumsum(nk))
    n <- length(y)
    N <- length(nk)
    nsub <- N
    ya <- list()
    for(i in 1:nsub){
        ya[[i]] <- y[(cn[i]+1):cn[i+1]]
    }
    nx = 2 + 2*k + p
    Oma <- diag(1,(nx))
    Q <- 1e+6
    betadiff <- 1
    iteration <- 0
    ga <- betaint
    ganew = ga
    maxit = 500
    w =1

    mu <- dotmu <- vmu <- r1 <- r2 <- D <- list()

    while(betadiff > 1e-6 & iteration < maxit){
        iteration = iteration + 1
        ga = ganew
        arsumg =0;arsumc =0; arsumgfirstdev =0;
        Q0 = Q; Oma0 =Oma;
        psi.est = ga[(2+p+k+1):(2+p+2*k)]
        psi.slo = ga[(2+p+1):(2+k+p)]
        eta.est = ga[1:(2+p+k)]
        PSI.est <- matrix(rep(psi.est,rep(n,k)),ncol=k)
        XREG.est <- cbind(1,XREG,pmax((X-PSI.est),0))
        BB <- matrix(rep(psi.slo,n),nrow=n,ncol=k,byrow=T)
        gdev <- cbind(XREG.est,BB*ifelse(X>PSI.est,-1,0))

        for(i in 1:nsub){
            xx <- as.matrix(XREG.est[(cn[i]+1):cn[i+1],])
            ggxx <- gdev[(cn[i]+1):cn[i+1],]
            mu[[i]] <- xx %*% eta.est
            dotmu[[i]] <- t(ggxx)
            vmu[[i]] <- diag(nk[i])
            r <- (max(diag(ggxx%*%Oma%*%t(ggxx)),1e-8))^(0.5)
            r1[[i]] <- pnorm(nsub^(0.5)*(ya[[i]]-mu[[i]])/r,0,1)-(1-tau)*matrix(1,nk[i],1)
            r2[[i]] <- nsub^(0.5)*dnorm(nsub^(0.5)*(ya[[i]]-mu[[i]])/r,0,1)/r
            D[[i]] <- diag(as.vector(r2[[i]]))

            if( type == "cs"){
                m1 <- matrix(1,nk[i],nk[i])-diag(nk[i])
                gi <- rbind(dotmu[[i]]%*%r1[[i]],dotmu[[i]]%*%m1%*%r1[[i]])/nsub
                di <- -rbind(dotmu[[i]]%*%D[[i]]%*%t(dotmu[[i]]),
                             dotmu[[i]]%*%m1%*%D[[i]]%*%t(dotmu[[i]]))/nsub
            }else if(type == "ar"){
                m1 <- matrix(0,nk[i],nk[i])
                m1[1,2] <- 1
                m1[nk[i],(nk[i]-1)] = 1
                for(j in 2:(nk[i]-1)){
                    m1[j,(j-1)] <- 1
                    m1[j,(j+1)] <- 1
                }
                gi <- rbind(dotmu[[i]]%*%r1[[i]],dotmu[[i]]%*%m1%*%r1[[i]])/nsub
                di <- -rbind(dotmu[[i]]%*%D[[i]]%*%t(dotmu[[i]]),
                             dotmu[[i]]%*%m1%*%D[[i]]%*%t(dotmu[[i]]))/nsub
            }else{
                stop("Please set the correct correlation structure within subject !")
            }
            arsumg <- arsumg + gi
            arsumc <- arsumc + gi%*%t(gi)
            arsumgfirstdev <- arsumgfirstdev + di

        }

        arsumc <- arsumc * nsub#-arsumg%*%t(arsumg)
        arcinv<- pinv(arsumc)

        arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
        arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
        invarqif2dev <- pinv(arqif2dev)
        ganew <- ganew - w*invarqif2dev%*%(arqif1dev)

        Oma <- pinv(t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev)
        Q <- t(arsumg) %*% arcinv %*% arsumg

        diff1 <- max(abs(Oma-Oma0))
        betadiff = sum(abs(ganew-ga))
        if(max(abs(arqif1dev)) < 1e-6) break

        if(abs(Q-Q0)<1e-6) break

        w = w/2

    }


    psi.est = ganew[(2+p+k+1):(2+p+2*k)]
    psi.slo = ganew[(2+p+1):(2+k+p)]
    eta.est = ganew[1:(2+p+k)]
    PSI.est <- matrix(rep(psi.est,rep(n,k)),ncol=k)
    XREG.est <- cbind(1,XREG,pmax((X-PSI.est),0))

    ese <- diag(pinv(t(arsumgfirstdev)%*%arcinv%*%arsumgfirstdev)/N)^(0.5)

    Err = y - XREG.est %*% eta.est
    rho <- sum((tau-(Err<0))*Err)

    return(list(ese=ese,rho=rho,bet.est=eta.est,psi.est=psi.est))

}

#' Fit the multi-kink quantile regression conditional on a given or pre-specified number of change points.
#'
#' @param y A numeric vector of response.
#' @param thre.x A numeric vector of scalar covariate with threshold effect.
#' @param cont.z A numeric matrix of design with constant slopes.
#' @param id A numeric vector of index used for longitudinal data; can be missing or NULL for iid data.
#' @param tau The quantile level that belongs to (0,1).
#' @param k The pre-specified number of change points.
#' @param psi Numeric vector to indicate the starting values for the change points. When psi=NULL (default), k quantiles are assumed.
#' @param bandwidth_type The bandwidth type. Specify one from "Hall-Sheather", "Bofinger", and "Chamberlain". Default is "Hall-Sheather".
#' @param control A list returned by fit.control.
#' @param est.type The estimation type for the longitudinal data. Specify one from "WI", "WC", corresponding to the working independence (WI) estimator and working correlation (WC) estimator. Default is "WI".
#' @param wi.type If est.type = "WI", then set the error structure of the variance-covariance matrix estimation. Specify one from "Compound", "AR" and "general".
#' @param wc.type If est.type = "WC", then set the correlation structure within subject. Specify one from "ar" and "cs". Default is "cs".
#'
#' @return A list containing the estimated regression coefficients with intercept (bet.est), the estimated standard error of the regression coefficients (bet.se), the estimated change points (psi.est), the estimated standard errors of threshold parameters (psi.se), and the fitted quantile objective value (rho).
#'
#' @examples
#' \dontrun{
#' # Simple examples for iid data type
#' n <- 500
#' Z1 <- rexp(n, 1)
#' Z2 <- rbinom(n, 1, 0.5)
#' Z <- cbind(Z1, Z2)
#' epsilon <- rnorm(n, 0, 1)
#' X <- runif(n, -2, 1)
#' psi <- c(-1, 0)
#' k <- length(psi)
#' Y <- XR %*% bet + epsilon
#' result_iid <- mkqr.fit(Y, X, Z, tau = 0.5, k = k)
#'
#' # Simple examples for longitudinal data
#' N <- 200
#' T <- 5
#' subject <- rep(1:N, each = T)
#' NT <- N * T
#' Z1 <- rexp(NT, 1)
#' Z2 <- rbinom(NT, 1, 0.5)
#' Z <- cbind(Z1, Z2)
#' epsilon <- rnorm(NT, 0, 1)
#' X <- runif(NT, -2, 1)
#' psi <- c(-1, 0)
#' k <- length(psi)
#' PSI <- matrix(rep(psi, rep(NT, k)), ncol = k)
#' a <- rnorm(N, 0, 1)
#' A <- rep(a, each = T)
#' XP <- matrix(rep(X, k), nrow = NT)
#' XR <- cbind(1, X, pmax((XP - PSI), 0), Z)
#' bet <- c(1, -1, 3, -3, sqrt(3), -sqrt(3))
#' Y <- XR %*% bet + A + epsilon
#' tau = 0.5
#' k = 2
#'
#' # Example 1: the working independence estimator; the error structure is "general"
#' est.type = "WI";
#' wi.type = "Compound"
#' result_WI_Compound <- mkqr.fit(y = Y, thre.x = X, cont.z = Z, id = subject, tau = tau,
#'                                k = k, est.type = est.type, wi.type = wi.type)
#'
#' # Example 2: the working correlated estimator; the correlation structure is "cs"
#' est.type = "WC";
#' wc.type = "cs"
#' result_WC_cs <- mkqr.fit(y = Y, thre.x = X, cont.z = Z, id = subject, tau = tau,
#'                          k = k, est.type = est.type, wc.type = wc.type)
#'
#' }
#' @export
mkqr.fit <- function(y,thre.x,cont.z,id,tau=0.5,k,psi=NULL,
                     bandwidth_type="Hall-Sheather",
                     control=fit.control(),
                     est.type = "WI",
                     wi.type = "general",
                     wc.type = "cs"){

    if(missing(k) || is.null(psi)){
        if(missing(k) || is.null(k)){
            stop("either psi or k must be given")
        }else{
            psi <- quantile(thre.x,seq(0,1,l=(k+2)))[-c(1,k+2)]
        }
    }
    k <- length(psi)
    if(missing(id) || is.null(id)){
        id <- seq(1:n)
    }
    density_fun <- function(y,XREG,X,PSI,residuals,bandwidth_type="Hall-Sheather"){
        eps <- .Machine$double.eps^(2/3)

        if (bandwidth_type == "Bofinger") {
            bandwidth <- n^(-1/5) * ((9/2 * dnorm(stats::qnorm(tau))^4)/(2 * qnorm(tau)^2 + 1)^2)^(1/5)
        } else if (bandwidth_type == "Chamberlain") {
            alpha <- 0.05
            bandwidth <- qnorm(1 - tau/2) * sqrt(alpha * (1 - alpha)/n)
        } else if (bandwidth_type == "Hall-Sheather") {
            alpha <- 0.05
            bandwidth <- n^(-1/3) * qnorm(1 - alpha/2)^(2/3) *
                ((3/2 * stats::dnorm(qnorm(tau))^2)/(2 * qnorm(tau)^2 + 1))^(1/3)
        } else {
            stop("Not a valid bandwith method!")
        }
        # Compute the density
        bu <- suppressWarnings(fitted.values(brisq(y,XREG,X,PSI,tau+bandwidth,control)$obj))
        bl <- suppressWarnings(fitted.values(brisq(y,XREG,X,PSI,tau-bandwidth,control)$obj))
        density <- pmax(0, 2 * bandwidth/((bu - bl) - eps))
        density
    }
    sub = model.matrix(~-1+factor(id))
    N = dim(sub)[2]
    n <- length(y)
    dev0 <- control$dev0
    if(missing(cont.z) || is.null(cont.z)){
        p <- 2; XREG <- matrix(thre.x,nrow=n,ncol=1);xreg.names <- c("x")
    }else{
        XREG <- cbind(cont.z,thre.x); p <- ncol(XREG)+1
        if(p==3) xreg.names <- c("z","x") else xreg.names <- c(paste("z",1:ncol(cont.z),sep=""),"x")
    }
    coef.names <- c("cons",xreg.names,paste0("(x-psi",1:k,")+"))
    psi.names <- paste0("psi",1:k)
    if(is.null(dev0)){
        obj0 <- rq(y ~XREG,tau,method="br")
        control$dev0 <- obj0$rho
    }
    control$xreg.names <- xreg.names
    control$stop.if.error <- TRUE
    X <- matrix(rep(thre.x,k),nrow=n)
    PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
    obj <- suppressWarnings(brisq(y,XREG,X,PSI,tau,control))
    rho <- obj$obj$rho
    residuals <- obj$obj$residuals
    df <- density_fun(y,XREG,X,PSI,residuals)
    psi.est <- obj$psi
    bet.est <- obj$obj$coefficients
    U.est <- bet.est[(p+1):(p+k)]

    PSI <- matrix(rep(psi.est,rep(n,k)),ncol=k)
    BB <- matrix(rep(U.est,n),nrow=n,ncol=k,byrow=T)
    ht <- cbind(1,XREG,pmax(X-PSI,0),BB*ifelse(X>PSI,-1,0))

    if(N == n){
        cat("Current estimation is performed for iid data. \n")
        Cn <- tau*(1-tau)/n*t(ht)%*%ht
        Dn <- t(ht)%*% diag(df) %*% ht/n
        Dn.inv <- solve(Dn+1e-8)
        Sig_thet <- n^(-1) * Dn.inv %*% Cn %*% Dn.inv
        ese <- as.vector(sqrt(diag(Sig_thet)))
        bet.se <- ese[1:(p+k)]
        psi.se <- ese[(p+k+1):(p+2*k)]
    }else if(N < n){
        if(est.type == "WI"){
            cat("Current estimation is performed for longitudinal data using working independence estimator. \n")
            Cn <- Cnmat(y,XREG,thre.x,ht,id,psi.est,bet.est,U.est,tau,p,residuals,wi.type)
            Dn <- t(ht)%*% diag(df) %*% ht/n
            Dn.inv <- solve(Dn+1e-8)
            Sig_thet <- n^(-1) * Dn.inv %*% Cn %*% Dn.inv
            ese <- as.vector(sqrt(diag(Sig_thet)))
            bet.se <- ese[1:(p+k)]
            psi.se <- ese[(p+k+1):(p+2*k)]
        }else if(est.type == "WC"){
            cat("Current estimation is performed for longitudinal data using working correlated estimator. \n")
            nk = table(id)
            betaint <- c(bet.est,psi.est)
            obj.wc <- fitting.k.GEE(y,XREG,X,nk,tau,k,p-2,wc.type,betaint)
            ese <- obj.wc$ese
            bet.est <- obj.wc$bet.est
            psi.est <- obj.wc$psi.est
            bet.se <- ese[1:(p+k)]
            psi.se <- ese[(p+k+1):(p+2*k)]
            rho <- obj.wc$rho
        }

    }


    names(bet.se) <- names(bet.est) <- coef.names
    names(psi.se) <- names(psi.est) <- psi.names

    return(list(bet.est=bet.est, bet.se=bet.se,
                psi.est=psi.est,psi.se=psi.se,
                rho = rho))

}

#' Fit the multi-kink quantile regression in the absence of the number of change points.
#'
#' @param y A numeric vector of response.
#' @param thre.x A numeric vector of scalar covariate with threshold effect.
#' @param cont.z A numeric matrix of design with constant slopes.
#' @param id A numeric vector of index used for longitudinal data; can be missing or NULL for iid data.
#' @param tau The quantile level that belongs to (0,1). Default is 0.5.
#' @param Cn A positive number corresponding to different types of BIC. Default is 1.
#' @param bandwidth_type The bandwidth type. Specify one from "Hall-Sheather", "Bofinger", and "Chamberlain". Default is "Hall-Sheather".
#' @param control A list returned by fit.control.
#' @param est.type The estimation type for the longitudinal data. Specify one from "WI", "WC", corresponding to the working independence (WI) estimator and working correlation (WC) estimator. Default is "WI".
#' @param wi.type If est.type = "WI", then set the error structure of the variance-covariance matrix estimation. Specify one from "Compound", "AR", and "general".
#' @param wc.type If est.type = "WC", then set the correlation structure within subject. Specify one from "ar" and "cs". Default is "cs".
#'
#' @return A list containing the estimated number of kink points (n.psi), the fitted quantile objective value (rho), estimated regression coefficients with intercept (bet.est), the estimated standard error of the regression coefficients (bet.se), the estimated change points (psi.est), and the estimated standard errors of threshold parameters (psi.se).
#'
#' @examples
#' \dontrun{
#' # Simple examples for iid data type
#' n <- 500
#' Z1 <- rexp(n,1)
#' Z2 <- rbinom(n,1,0.5)
#' Z <- cbind(Z1,Z2)
#' epsilon <- rnorm(n,0,1)
#' X <- runif(n,-2,1)
#' psi <- c(-1,0)
#' k <- length(psi)
#' PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
#' XP <- matrix(rep(X,k),nrow=n)
#' XR <- cbind(1,X,pmax((XP-PSI),0),Z)
#' bet <- c(1,-1,3,-3,sqrt(3),-sqrt(3))
#' Y <- XR %*% bet + epsilon
#'
#' # Estimation setting
#' tau <- 0.5
#' K.max <- 5
#' control <- fit.control(K.max = K.max)
#' Cn <- 1
#' mkqr.bea(y = Y, thre.x = X, cont.z = Z, tau = tau, Cn = Cn, control = control)
#'
#' # Simple examples for longitudinal data
#' N <- 200
#' T <- 5
#' subject <- rep(1:N, each = T)
#' NT <- N * T
#' Z1 <- rexp(NT, 1)
#' Z2 <- rbinom(NT, 1, 0.5)
#' Z <- cbind(Z1, Z2)
#' epsilon <- rnorm(NT, 0, 1)
#' X <- runif(NT, -2, 1)
#' psi <- c(-1, 0)
#' k <- length(psi)
#' PSI <- matrix(rep(psi, rep(NT, k)), ncol = k)
#' a <- rnorm(N, 0, 1)
#' A <- rep(a, each = T)
#' XP <- matrix(rep(X, k), nrow = NT)
#' XR <- cbind(1, X, pmax((XP - PSI), 0), Z)
#' bet <- c(1, -1, 3, -3, sqrt(3), -sqrt(3))
#' Y <- XR %*% bet + A + epsilon
#' tau <- 0.5
#' k <- 2
#'
#' # Example 1: the working independence estimator; the error structure is "general"
#' est.type <- "WI"
#' wi.type <- "Compound"
#' tau <- 0.5
#' K.max <- 5
#' control <- fit.control(K.max = K.max)
#' Cn <- 1
#' mkqr.bea(y = Y, thre.x = X, cont.z = Z, id = subject, tau = tau, Cn = Cn,
#'          control = control, est.type = est.type, wi.type = wi.type)
#'
#' # Example 2: the working correlated estimator; the correlation structure is "cs"
#' est.type <- "WC"
#' wc.type <- "cs"
#' tau <- 0.5
#' K.max <- 5
#' control <- fit.control(K.max = K.max)
#' Cn <- 1
#' mkqr.bea(y = Y, thre.x = X, cont.z = Z, id = subject, tau = tau, Cn = Cn,
#'          control = control, est.type = est.type, wc.type = wc.type)
#' }
#'
#' @export
mkqr.bea <- function(y,thre.x,cont.z,id,tau=0.5,Cn=1,
                     bandwidth_type="Hall-Sheather",
                     control=fit.control(),
                     est.type = "WI",
                     wi.type = "general",
                     wc.type = "cs"){
    n <- length(y)
    if(missing(cont.z) || is.null(cont.z)){
        p <- 2; XREG <- matrix(thre.x,nrow=n,ncol=1);xreg.names <- c("x")
    }else{
        XREG <- cbind(cont.z,thre.x); p <- ncol(XREG)+1
        if(p==3) xreg.names <- c("z","x") else xreg.names <- c(paste("z",1:ncol(cont.z),sep=""),"x")
    }
    dev0 <- control$dev0
    if(is.null(dev0)){
        obj0 <- rq(y~XREG,tau,method="br")
        control$dev0 <- obj0$rho
    }
    if(is.null(control$grid)) control$grid <- (max(thre.x)-min(thre.x))/30
    if(missing(id) || is.null(id)){
        id <- seq(1:n)
    }
    control$xreg.names <- xreg.names
    kmax <- control$K.max
    psi0 <- quantile(thre.x,seq(0,1,l=kmax+2)[-c(1,(kmax+2))],names=FALSE)
    control$stop.if.error <- FALSE
    X <- matrix(rep(thre.x,kmax),nrow=n)
    PSI <- matrix(rep(psi0,rep(n,kmax)),ncol=kmax)
    obj <- suppressWarnings(brisq.fit(y,XREG,X,PSI,tau,control,return.all.sol=FALSE))
    if(obj$n.psi==0) stop("There is no kink effect. \n")
    psi0 <- obj$psi
    k <- obj$n.psi
    bic <- rep(NA,k)
    control$stop.if.error <- TRUE
    for(kk in 1:k){
        obj.k <- try(mkqr.fit(y=y,thre.x=thre.x,cont.z=cont.z,id=id,tau=tau,k=kk, bandwidth_type=bandwidth_type,
                              control=control,est.type=est.type, wi.type=wi.type,wc.type=wc.type),silent=TRUE)
        if(!is(obj.k,"try-error")){
            bic[kk] <- log(obj.k$rho) + (p+2*kk)*log(n)/2/n*Cn
        }
    }
    n.psi <- which.min(bic)
    oo <- mkqr.fit(y=y,thre.x=thre.x,cont.z=cont.z,id=id,tau=tau,k=n.psi, bandwidth_type=bandwidth_type,
                   control=control,est.type=est.type, wi.type=wi.type,wc.type=wc.type)
    oo$n.psi <- n.psi
    return(oo)
}




