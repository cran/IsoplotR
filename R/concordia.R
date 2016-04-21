#' Plot concordia diagrams
#'
#' Wetherill and Tera-Wasserburg concordia diagrams
#'
#' @param X a U-Pb dataset, i.e. a matrix with the following columns:
#'     ID, Pb207Pb206, sePb207Pb206, Pb206U238, sePb206U238,
#'     Pb207U235, sePb207U235
#' @param limits age limits of the concordia line
#' @param alpha confidence cutoff for the error ellipses
#' @param wetherill boolean flag (FALSE for Tera-Wasserburg)
#' @param show.labels boolean flag (TRUE to show grain numbers)
#' @param ellipse.col background colour of the error ellipses
#' @param concordia.col colour of the concordia line
#' @importFrom grDevices rgb
#' @examples
#' data(UPb)
#' concordia.plot(UPb)
#' @export
concordia.plot <- function(X,limits=NA,alpha=0.05,wetherill=TRUE,show.labels=FALSE,
                           ellipse.col=rgb(0,1,0,0.5),concordia.col=rgb(1,0,0,0.5)){
    concordia.line(X,limits,wetherill,concordia.col)
    if (wetherill){
        vars <- c('Pb207U235','Pb206U238')
    } else {
        vars <- c('Pb206U238','Pb207Pb206')
    }
    for (sname in names(X)){
        x0 <- X$x[sname,vars[1]]
        y0 <- X$x[sname,vars[2]]
        covmat <- get.covmat.UPb(X,sname)[vars,vars]
        if (!wetherill){
            J <- matrix(c(-1/(x0*x0),0,0,1),nrow=2)
            covmat <- J %*% covmat %*% t(J)
            x0 <- 1/x0
        }
        ell <- get.ellipse(x0,y0,covmat,alpha=alpha)
        graphics::polygon(ell$x,ell$y,col=ellipse.col)
        graphics::points(x0,y0,pch=19,cex=0.25)
        if (show.labels) graphics::text(x0,y0,sname)
    }
}

# helper function for plot.concordia
concordia.line <- function(X,limits,wetherill,col){
    lims <- get.concordia.limits(X,limits,wetherill)
    nn <- 100
    range.t <- lims$max.t-lims$min.t
    m <- max(0.8*lims$min.t,lims$min.t-range.t/20)
    M <- min(1.2*lims$max.t,lims$max.t+range.t/20)
    tt <- seq(m,M,length.out=nn)
    concordia.x <- 0*rep(tt,2)
    concordia.y <- 0*rep(tt,2)
    for (i in 1:nn){
        if (wetherill){
            concordia.x[i] <- get.Pb207U235(tt[i])$x - 2*get.Pb207U235(tt[i])$e
            concordia.y[i] <- get.Pb206U238(tt[i])$x - 2*get.Pb206U238(tt[i])$e
            concordia.x[2*nn-i+1] <- get.Pb207U235(tt[i])$x + 2*get.Pb207U235(tt[i])$e
            concordia.y[2*nn-i+1] <- get.Pb206U238(tt[i])$x + 2*get.Pb206U238(tt[i])$e
        } else {
            U238Pb206 <- 1/get.Pb206U238(tt[i])$x
            errU238Pb206 <- U238Pb206 * get.Pb206U238(tt[i])$e / get.Pb206U238(tt[i])$x
            concordia.x[i] <- U238Pb206 - 2*errU238Pb206
            concordia.y[i] <- get.Pb207Pb206(tt[i])$x - 2*get.Pb207Pb206(tt[i])$e
            concordia.x[2*nn-i+1] <- U238Pb206 + 2*errU238Pb206
            concordia.y[2*nn-i+1] <- get.Pb207Pb206(tt[i])$x + 2*get.Pb207Pb206(tt[i])$e
        }
    }
    ticks <- pretty(tt)
    xt <- 0*ticks
    yt <- 0*ticks
    for (i in 1:length(ticks)){
        if (wetherill){
            xt[i] <- get.Pb207U235(ticks[i])$x
            yt[i] <- get.Pb206U238(ticks[i])$x
        } else {
            xt[i] <- 1/get.Pb206U238(ticks[i])$x
            yt[i] <- get.Pb207Pb206(ticks[i])$x
        }
    }
    if (wetherill){
        x.lab <- expression(paste(""^"207","Pb/"^"235","U"))
        y.lab <- expression(paste(""^"206","Pb/"^"238","U"))
    } else {
        x.lab <- expression(paste(""^"238","U/"^"206","Pb"))
        y.lab <- expression(paste(""^"207","Pb/"^"206","Pb"))
    }
    graphics::plot(c(lims$min.x,lims$max.x),c(lims$min.y,lims$max.y),
                   type='n',xlab=x.lab,ylab=y.lab)
    graphics::polygon(concordia.x,concordia.y,col=col)
    graphics::points(xt,yt,pch=21,bg='white')
    graphics::text(xt,yt,as.character(ticks),pos=2)
}

get.concordia.limits <- function(X,limits,wetherill){
    out <- list(min.t=0,max.t=0,min.x=0,max.x=0,min.y=0,max.y=0)
    if (!is.na(limits) && wetherill){
        out$min.t <- limits[1]
        out$max.t <- limits[2]
        out$min.x <- get.Pb207U235(out$min.t)$x
        out$max.x <- get.Pb207U235(out$max.t)$x
        out$min.y <- get.Pb206U238(out$min.t)$x
        out$max.y <- get.Pb206U238(out$max.t)$x
    } else if (!is.na(limits) && !wetherill){
        out$min.t <- limits[1]
        out$max.t <- limits[2]
        out$min.x <- 1/get.Pb206U238(out$max.t)$x
        out$max.x <- 1/get.Pb206U238(out$min.t)$x
        out$min.y <- get.Pb207Pb206(out$max.t)$x
        out$max.y <- get.Pb207Pb206(out$min.t)$x
    } else if (is.na(limits) && wetherill) {
        out$min.x <- min(X$x[,'Pb207U235']-2*X$x[,'errPb207U235'])
        out$max.x <- max(X$x[,'Pb207U235']+2*X$x[,'errPb207U235'])
        out$min.y <- min(X$x[,'Pb206U238']-2*X$x[,'errPb206U238'])
        out$max.y <- max(X$x[,'Pb206U238']+2*X$x[,'errPb206U238'])
        out$min.t <- get.Pb206U238age(out$min.y)
        out$max.t <- get.Pb207U235age(out$max.x)
    } else if (is.na(limits) && !wetherill){
        out$min.x <- 1/(max(X$x[,'Pb206U238']+2*X$x[,'errPb206U238']))
        out$max.x <- 1/min(X$x[,'Pb206U238']-2*X$x[,'errPb206U238'])
        out$min.y <- min(X$x[,'Pb207Pb206']-2*X$x[,'errPb207Pb206'])
        out$max.y <- max(X$x[,'Pb207Pb206']+2*X$x[,'errPb207Pb206'])
        out$min.t <- get.Pb206U238age(1/out$max.x)
        out$max.t <- get.Pb207Pb206age(out$max.y)
    }
    out
}
