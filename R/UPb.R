#' Extract a group of samples from a U/Pb dataset
#'
#' Returns a subset of U/Pb data matching a certain prefix
#'
#' @param x a dataset of class 'UPb'
#' @param subset logical expression indicating elements or rows to
#'     keep: missing values are taken as false.
#' @param select a vector of sample names
#' @param ... optional arguments for the generic subset function
#' @return an object of class 'UPb'
#' @examples
#' fname <- system.file("UPb.csv",package="IsoplotR")
#' UPb <- read.UPb(fname)
#' #data(UPb)
#' ABC <- subset(UPb,select=c('A','B','C','D','E','F','G','H','I'))
#' concordia.plot(ABC,wetherill=TRUE)
#' @export
subset.UPb <- function(x,subset=NULL,select=NULL,...){
    out <- x
    if (!is.null(subset)){
        i <- which(subset,arr.ind=TRUE)
    } else if (!is.null(select)){
        i <- which(select %in% names(x))
    } else {
        return(out)
    }
    out$x <- x$x[i,]
    out
}

#' Get the names of the samples contained in a U/Pb dataset
#'
#' Returns the rownames of a U/Pb dataset
#'
#' @param x an object of class 'UPb'
#' @return an object of class 'UPb'
#' @examples
#' data(UPb)
#' names(UPb)
#' @export
names.UPb <- function(x){
    rownames(x$x)
}

# Get the predicted parent daughter ratio
get.DP <- function(age,P){
    out <- list()
    out$x <- exp(lambda(P)*age)-1
    out$e <- abs(lambda.err(P)*(exp(lambda(P)*age)*age))
    out
}

#' Get the covariance matrix of a U/Pb samples
#'
#' Returns the covariance matrix of the U/Pb and Pb/Pb ratios of the ith sample
#'
#' @param X an object of class 'UPb'
#' @param i the index of the sample of interest
#' @return a 3x3 covariance matrix
#' @examples
#' data(UPb)
#' get.covmat.UPb(UPb,2)
#' @export
get.covmat.UPb <- function(X,i){
    covmat <- matrix(rep(0,9),nrow=3)
    if (X$format == 1){
        rownames(covmat) <- c('Pb207Pb206','Pb206U238','Pb207U235')
        relvar207 <- 0.5 * ((X$x[i,'errPb207Pb206']/X$x[i,'Pb207Pb206'])^2 +
                            (X$x[i,'errPb207U235']/X$x[i,'Pb207U235'])^2 -
                            (X$x[i,'errPb206U238']/X$x[i,'Pb206U238'])^2)
        relvar206 <- 0.5 * ((X$x[i,'errPb207Pb206']/X$x[i,'Pb207Pb206'])^2 -
                            (X$x[i,'errPb207U235']/X$x[i,'Pb207U235'])^2 +
                            (X$x[i,'errPb206U238']/X$x[i,'Pb206U238'])^2)
        relvar238 <- 0.5 * ((X$x[i,'errPb207U235']/X$x[i,'Pb207U235'])^2 +
                            (X$x[i,'errPb206U238']/X$x[i,'Pb206U238'])^2 -
                            (X$x[i,'errPb207Pb206']/X$x[i,'Pb207Pb206'])^2)
        covmat[1,1] <- X$x[i,'errPb207Pb206']^2
        covmat[2,2] <- X$x[i,'errPb206U238']^2
        covmat[3,3] <- X$x[i,'errPb207U235']^2
        covmat[1,2] <- -relvar206*X$x[i,'Pb207Pb206']*X$x[i,'Pb206U238']
        covmat[1,3] <- relvar207/(X$x[i,'Pb207Pb206']*X$x[i,'Pb207U235'])
        covmat[2,3] <- relvar238*X$x[i,'Pb206U238']*X$x[i,'Pb207U235']
        covmat[2,1] <- covmat[1,2]
        covmat[3,2] <- covmat[2,3]
        covmat[3,1] <- covmat[1,3]
    }
    colnames(covmat) <- rownames(covmat)
    covmat
}

#' Calculate a Pb207/U235 age
#'
#' Returns the U/Pb age for a given Pb207/U235 ratio
#'
#' @param Pb207U235 a scalar
#' @return the Pb207/U235 age [in Ma]
#' @examples
#' get.Pb207U235age(88.8)
#' @export
get.Pb207U235age <- function(Pb207U235){
    log(1+Pb207U235)/lambda('U235')
}

#' Calculate a Pb206/U238 age
#'
#' Returns the U/Pb age for a given Pb206/U238 ratio
#'
#' @param Pb206U238 a scalar
#' @return the Pb206/U238 age [in Ma]
#' @examples
#' get.Pb206U238age(1.03)
#' @export
get.Pb206U238age <- function(Pb206U238){
    log(1+Pb206U238)/lambda('U238')
}

#' Calculate a Pb207/Pb206 age
#'
#' Numerically solves the Pb207/Pb206 age equation
#' @param Pb207Pb206 the Pb207/Pb206 ratio
#' @return the Pb207/Pb206 age [in Ma]
#' @examples
#' get.Pb207Pb206age(0.625)
#' @export
get.Pb207Pb206age <- function( Pb207Pb206){
    Pb207Pb206.misfit <- function(x,y) { (get.Pb207Pb206(x)$x - y)^2 }
    out <- stats::optimize(Pb207Pb206.misfit,c(0,4600),y=Pb207Pb206)
    out$minimum
}

#' Calculate a Pb206/U238 ratio
#'
#' Returns the expected Pb206/U238 ratio for a given age
#'
#' @param age a U/Pb age [in Ma]
#' @return a Pb206/U238 ratio
#' @examples
#' get.Pb206U238(4567)
#' @export
get.Pb206U238 <- function(age){
    get.DP(age,'U238')
}

#' Calculate a Pb207/U235 ratio
#'
#' Returns the expected Pb207/U235 ratio for a given age
#'
#' @param age a U/Pb age [in Ma]
#' @return a Pb207/U235 ratio
#' @examples
#' get.Pb207U235(4567)
#' @export
get.Pb207U235 <- function(age){
    get.DP(age,'U235')
}

#' Calculate a Pb207/Pb206 ratio
#'
#' Returns the expected Pb207/Pb206 ratio for a given age
#'
#' @param age a Pb/Pb age [in Ma]
#' @return a Pb207/Pb206 ratio
#' @examples
#' get.Pb207Pb206(4567)
#' @export
get.Pb207Pb206 <- function(age){
    Pb206U238 <- get.Pb206U238(age)
    Pb207U235 <- get.Pb207U235(age)
    out <- list()
    out$x <- (I.A('U235')*Pb207U235$x)/(I.A('U238')*Pb206U238$x)
    out$e <- out$x * sqrt((Pb207U235$e/Pb207U235$x)^2 + (Pb207U235$e/Pb207U235$x)^2)
    out
}
