#' @title
#' Visualise U-Th-He data on a logratio plot or ternary diagram
#'
#' @description
#' Plot U-Th(-Sm)-He data on a (log[He/Th] vs. log[U/He]) logratio
#' plot or U-Th-He ternary diagram
#'
#' @details
#' U, Th, Sm and He are \emph{compositional} data.  This means that it
#' is not so much the absolute concentrations of these elements that
#' bear the chronological information, but rather their relative
#' proportions. The space of all possible U-Th-He compositions fits
#' within the constraints of a ternary diagram or `helioplot'
#' (Vermeesch, 2008, 2010). If Sm is included as well, then this
#' expands to a three-dimensional tetrahaedral space (Vermeesch,
#' 2008).  Data that fit within these constrained spaces must be
#' subjected to a logratio transformation prior to statistical
#' analysis (Aitchison, 1986).  In the case of the U-Th-He-(Sm)-He
#' system, this is achieved by first defining two (or three) new
#' variables:
#'
#' \eqn{u \equiv \ln[U/He]}
#' \eqn{v \equiv \ln[Th/He]}
#' \eqn{(, w \equiv \ln[Sm/He] )}
#'
#' and then performing the desired statistical analysis (averaging,
#' uncertainty propagation, ...) on the transformed data. Upon
#' completion of the mathematical operations, the results can then be
#' mapped back to U-Th-(Sm)-He space using an inverse logratio
#' transformation:
#'
#' \eqn{[He] = 1/[e^{u}+e^{v}+(e^{w})+1]},
#' \eqn{[U] = e^{u}/[e^{u}+e^{v}+(e^{w})+1]}\cr
#' \eqn{[Th] = e^{v}/[e^{u}+e^{v}+(e^{w})+1]},
#' \eqn{([Sm] = e^{w}/[e^{u}+e^{v}+(e^{w})+1])}
#'
#' where \eqn{[He] + [U] + [Th] (+ [Sm]) = 1}. In the context of
#' U-Th-(Sm)-He dating, the \emph{barycentric} age (which is
#' equivalent to the 'central age' of Vermeesch, 2008) is defined as
#' the date that corresponds to the compositional mean, which is
#' equivalent to the arithmetic mean composition in logratio space.
#' \code{IsoplotR}'s \code{helioplot} function performs this
#' calculation using the same algorithm that is used to obtain the
#' weighted mean U-Pb composition for the \code{\link{concordia}} age
#' calculation. Overdispersion is treated similarly as in a regression
#' context (see \code{\link{isochron}}).  Thus, there are options to
#' augment the uncertainties with a factor \eqn{\sqrt{MSWD}} (model
#' 1); to ignore the analytical uncertainties altogether (model 2); or
#' to add a constant overdispersion term to the analytical
#' uncertainties (model 3).  The \code{helioplot} function visualises
#' U-Th-(Sm)-He data on either a ternary diagram or a bivariate
#' \eqn{\ln[Th/U]} vs. \eqn{\ln[U/He]} contour plot. These diagrams
#' provide a convenient way to simultaneously display the isotopic
#' composition of samples as well as their chronological meaning. In
#' this respect, they fulfil the same purpose as the U-Pb
#' \code{\link{concordia}} diagram and the U-series
#' \code{\link{evolution}} plot.
#'
#' @param x an object of class \code{UThHe}
#' @param logratio Boolean flag indicating whether the data should be
#'     shown on bivariate log[He/Th] vs. log[U/He] diagram, or a
#'     U-Th-He ternary diagram.
#' @param show.barycentre show the mean composition as a white
#'     ellipse?
#' @param show.numbers show the grain numbers inside the error
#'     ellipses?
#' @param oerr indicates whether the analytical uncertainties of the
#'     output are reported in the plot title as:
#' 
#' \code{1}: 1\eqn{\sigma} absolute uncertainties.
#' 
#' \code{2}: 2\eqn{\sigma} absolute uncertainties.
#' 
#' \code{3}: absolute (1-\eqn{\alpha})\% confidence intervals, where
#' \eqn{\alpha} equales the value that is stored in
#' \code{settings('alpha')}.
#' 
#' \code{4}: 1\eqn{\sigma} relative uncertainties (\eqn{\%}).
#' 
#' \code{5}: 2\eqn{\sigma} relative uncertainties (\eqn{\%}).
#'
#' \code{6}: relative (1-\eqn{\alpha})\% confidence intervals, where
#' \eqn{\alpha} equales the value that is stored in
#' \code{settings('alpha')}.
#'
#' @param contour.col two-element vector with the fill colours to be
#'     assigned to the minimum and maximum age contour
#' @param levels a vector with additional values to be displayed as
#'     different background colours within the error ellipses.
#' @param clabel label of the colour scale
#' @param ellipse.fill
#' Fill colour for the error ellipses. This can either be a single
#' colour or multiple colours to form a colour ramp. Examples:
#'
#' a single colour: \code{rgb(0,1,0,0.5)}, \code{'#FF000080'},
#' \code{'white'}, etc.;
#'
#' multiple colours: \code{c(rbg(1,0,0,0.5)},
#' \code{rgb(0,1,0,0.5))}, \code{c('#FF000080','#00FF0080')},
#' \code{c('blue','red')}, \code{c('blue','yellow','red')}, etc.;
#'
#' a colour palette: \code{rainbow(n=100)},
#' \code{topo.colors(n=100,alpha=0.5)}, etc.; or
#'
#' a reversed palette: \code{rev(topo.colors(n=100,alpha=0.5))}, etc.
#' 
#' For empty ellipses, set \code{ellipse.fill=NA}
#' @param ellipse.stroke the stroke colour for the error
#'     ellipses. Follows the same formatting guidelines as
#'     \code{ellipse.fill}
#'
#' @param sigdig number of significant digits for the barycentric age
#' @param xlim optional limits of the x-axis (log[U/He]) of the
#'     logratio plot. If \code{xlim=NA}, the axis limits are
#'     determined automatically.
#' @param ylim optional limits of the y-axis (log[Th/He]) of the
#'     logratio plot. If \code{ylim=NA}, the axis limits are
#'     determined automatically.
#' @param fact three-element vector with scaling factors of the
#'     ternary diagram if \code{fact=NA}, these will be determined
#'     automatically
#' @param model choose one of the following statistical models:
#'
#' \code{1}: weighted mean. This model assumes that the scatter between
#' the data points is solely caused by the analytical uncertainty. If
#' the assumption is correct, then the MSWD value should be
#' approximately equal to one. There are three strategies to deal with
#' the case where MSWD>1. The first of these is to assume that the
#' analytical uncertainties have been underestimated by a factor
#' \eqn{\sqrt{MSWD}}.
#'
#' \code{2}: unweighted mean. A second way to deal with over- or
#' underdispersed datasets is to simply ignore the analytical
#' uncertainties.
#'
#' \code{3}: weighted mean with overdispersion: instead of attributing
#' any overdispersion (MSWD > 1) to underestimated analytical
#' uncertainties (model 1), it can also be attributed to the presence
#' of geological uncertainty, which manifests itself as an added
#' (co)variance term.
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from the barycentric age calculation.
#' @param omit.fill fill colour that should be used for the omitted
#'     aliquots.
#' @param omit.stroke stroke colour that should be used for the
#'     omitted aliquots.
#' @param ... optional arguments to the generic \code{plot} function
#' @seealso \code{\link{radialplot}}
#' @references Aitchison, J., 1986, The statistical analysis of
#'     compositional data: London, Chapman and Hall, 416 p.
#'
#' Vermeesch, P., 2008. Three new ways to calculate average (U-Th)/He
#' ages. Chemical Geology, 249(3), pp.339-347.
#'
#' Vermeesch, P., 2010. HelioPlot, and the treatment of overdispersed
#'     (U-Th-Sm)/He data. Chemical Geology, 271(3), pp.108-111.
#' @examples
#' attach(examples)
#' helioplot(UThHe)
#' dev.new()
#' helioplot(UThHe,logratio=FALSE)
#' @export
helioplot <- function(x,logratio=TRUE,model=1,show.barycentre=TRUE,
                      show.numbers=FALSE,oerr=3,contour.col=c('white','red'),
                      levels=NULL,clabel="",ellipse.fill=c("#00FF0080","#0000FF80"),
                      ellipse.stroke='black',sigdig=2,xlim=NA,
                      ylim=NA,fact=NA,hide=NULL,omit=NULL,
                      omit.fill=NA,omit.stroke='grey',...){
    ns <- length(x)
    calcit <- (1:ns)%ni%c(hide,omit)
    plotit <- (1:ns)%ni%hide
    x2calc <- clear(x,hide,omit)
    x2plot <- clear(x,hide)
    fit <- central.UThHe(x2calc,compositional=TRUE,model=model)
    fill <- set_ellipse_colours(ns=ns,levels=levels,
                                col=ellipse.fill,hide=hide,
                                omit=omit,omit.col=omit.fill)
    stroke <- set_ellipse_colours(ns=ns,levels=levels,
                                  col=ellipse.stroke,hide=hide,
                                  omit=omit,omit.col=omit.stroke)
    if (logratio) {
        plot_logratio_contours(x2plot,contour.col=contour.col,
                               xlim=xlim,ylim=ylim)
        if (model==2){
            u <- log(x[,'U']/x[,'He'])
            v <- log(x[,'Th']/x[,'He'])
            plot_points(u,v,show.numbers=show.numbers,
                        hide=hide,omit=omit,bg=fill,col=stroke,...)
        } else {
            plot_logratio_ellipses(x,fill=fill,stroke=stroke,
                                   oerr=oerr,levels=levels,
                                   show.numbers=show.numbers,hide=hide)
        }
    } else {
        if (all(is.na(fact))) fact <- getfact(x2plot,fit)
        plot_helioplot_contours(x2plot,fact=fact,
                                contour.col=contour.col,
                                xlim=xlim,ylim=ylim)
        if (model==2){
            plot_helioplot_points(x,show.numbers=show.numbers,
                                  fact=fact,hide=hide,omit=omit,
                                  bg=fill,col=stroke,...)
        } else {
            plot_helioplot_ellipses(x,fill=fill,stroke=stroke,fact=fact,
                                    oerr=oerr,levels=levels,
                                    show.numbers=show.numbers,hide=hide)
        }
    }
    if (show.barycentre){
        plot_barycentre(fit,fact=fact,logratio=logratio,
                        oerr=oerr,doSm=doSm(x))
    }
    fit$n <- length(which(calcit))
    graphics::title(helioplot_title(fit,sigdig=sigdig,oerr=oerr))
    invisible(colourbar(z=levels[calcit],fill=ellipse.fill,
                        stroke=ellipse.stroke,clabel=clabel))
}

plot_logratio_frame <- function(lims,...){
    graphics::plot(lims[1:2],lims[3:4],type='n',bty='n',
                   xlab='log[U/He]',ylab='log[Th/He]',...)
}

plot_helioplot_frame <- function(lims,fact=c(1,1,1),fill.col=NA,...){
    graphics::plot(c(0,1),c(0,1),type='n',xaxt='n',yaxt='n',
                   xlab='',ylab='',asp=1,bty='n',...)
    corners <- xyz2xy(matrix(c(1,0,0,1,0,1,0,0,0,0,1,0),ncol=3))
    graphics::polygon(corners,col=fill.col)
    HeLabel <- paste0(fact[1],' x He')
    ULabel <- paste0(fact[3],' x U')
    ThLabel <- paste0(fact[2],' x Th')
    labels <- c(HeLabel,ULabel,ThLabel)
    graphics::text(corners[1:3,],labels=labels,pos=c(3,1,1),xpd=NA)
}

plot_logratio_ellipses <- function(x,fill,stroke,oerr=3,
                                   show.numbers=FALSE,levels=NULL,hide=NULL){
    sn <- clear(1:length(x),hide)
    for (i in sn){
        uvc <- UThHe2uv_covmat(x,i)
        x0 <- uvc$uv[1]
        y0 <- uvc$uv[2]
        ell <- ellipse(x=x0,y=y0,covmat=uvc$covmat,alpha=oerr2alpha(oerr))
        graphics::polygon(ell,col=fill[i],border=stroke[i])
        if (show.numbers) graphics::text(x0,y0,labels=i)
        else graphics::points(x0,y0,pch=19,cex=0.25)
    }
}
plot_helioplot_ellipses <- function(x,fill,stroke,fact=c(1,1,1),oerr=3,
                                    show.numbers=FALSE,levels=NULL,hide=NULL){
    sn <- clear(1:length(x),hide)
    for (i in sn){
        uvc <- UThHe2uv_covmat(x,i)
        x0 <- uvc$uv[1]
        y0 <- uvc$uv[2]
        ell <- ellipse(x=x0,y=y0,covmat=uvc$covmat,alpha=oerr2alpha(oerr))
        HeUTh0 <- uv2HeUTh(uvc$uv)
        HeUTh <- uv2HeUTh(ell)
        xyz <- renormalise(HeUTh,fact=fact)
        xy <- xyz2xy(xyz)
        graphics::polygon(xy,col=fill[i],border=stroke[i])
        xyz0 <- renormalise(HeUTh0,fact=fact)
        x0y0 <- xyz2xy(xyz0)
        if (show.numbers) graphics::text(x0y0[1],x0y0[2],i)
        else graphics::points(x0y0[1],x0y0[2],pch=19,cex=0.25)
    }
}
plot_helioplot_points <- function(x,fact=c(1,1,1),bg=NA,
                                  show.numbers=FALSE,hide=NULL,
                                  omit=NULL,...){
    xyz <- renormalise(x[,c('He','U','Th'),drop=FALSE],fact=fact)
    xy <- xyz2xy(xyz)
    plot_points(xy[,1],xy[,2],show.numbers=show.numbers,
                hide=hide,omit=omit,bg=bg,...)
}

plot_barycentre <- function(fit,fact=c(1,1,1),logratio=TRUE,
                            oerr=3,doSm=TRUE,...){
    ell <- ellipse(x=fit$uvw[1],y=fit$uvw[2],
                   covmat=fit$covmat[1:2,1:2],alpha=oerr2alpha(oerr))
    if (logratio){
        graphics::polygon(ell,col="#FFFFFFBF")
    } else {
        HeUTh <- uv2HeUTh(ell)
        xyz <- renormalise(HeUTh,fact=fact)
        xy <- xyz2xy(xyz)
        graphics::polygon(xy,col="#FFFFFFBF")
    }
}

plot_logratio_contours <- function(x,contour.col=c('white','red'),
                                   xlim=NA,ylim=NA,...){
    cntrs <- get_logratio_contours(x,xlim=xlim,ylim=ylim)
    plot_logratio_frame(cntrs$lims,...)
    tticks <- cntrs$tticks
    nt <- length(tticks)
    crp <- grDevices::colorRampPalette(contour.col)(nt)
    for (i in 1:nt){
        uv.plot <- cntrs$uv[[i]]
        uv.first <- uv.plot[1,]
        uv.last <- uv.plot[nrow(uv.plot),]
        uv.plot <- rbind(uv.plot,c(uv.last[1],cntrs$lims[3]))
        uv.plot <- rbind(uv.plot,cntrs$lims[c(1,3)])
        uv.plot <- rbind(uv.plot,c(cntrs$lims[1],uv.first[2]))
        uv.plot <- rbind(uv.plot,uv.first)
        graphics::polygon(uv.plot,col=crp[i])
        graphics::text(uv.first[1],uv.first[2],
                       labels=paste(tticks[i],'Ma'),pos=4)
    }
}

plot_helioplot_contours <- function(x,fact=c(1,1,1),
                                    contour.col=c('white','red'),
                                    xlim=NA,ylim=NA,...){
    cntrs <- get_helioplot_contours(x,fact=fact)
    plot_helioplot_frame(cntrs$lims,fact=fact,
                         fill.col=contour.col[2],...)
    tticks <- cntrs$tticks
    nt <- length(tticks)
    crp <- grDevices::colorRampPalette(contour.col)(nt+1)
    for (i in nt:1){
        xyz <- cntrs$xyz[[i]]
        xyz.first <- xyz[1,]
        xyz <- rbind(xyz,c(0,0,1),c(0,1,0),xyz.first)
        xyz <- renormalise(xyz,fact=fact)
        xy <- xyz2xy(xyz)
        graphics::polygon(xy,col=crp[i])
        label <- paste(signif(cntrs$tticks[[i]],2),'Ma')
        graphics::text(xy[1,1],xy[1,2],labels=label,pos=2)
    }
}

helioplot_title <- function(fit,sigdig=2,oerr=3,...){
    line1 <- maintit(x=fit$age[1],sx=fit$age[-1],n=fit$n,df=fit$df,
                     sigdig=sigdig,oerr=oerr,prefix="barycentric age =")
    line1line <- 1
    if (fit$model==1){
        line2 <- mswdtit(mswd=fit$mswd,p=fit$p.value,sigdig=sigdig)
        mymtext(line2,line=0,...)
    } else if (fit$model==2){
        line1line <- 0
    } else if (fit$model==3){
        line2 <- disptit(fit$disp[1],fit$disp[-1],sigdig=sigdig,oerr=oerr)
        mymtext(line2,line=0,...)
    }
    mymtext(line1,line=line1line,...)
}

get_logratio_contours <- function(x,xlim=NA,ylim=NA,res=50){
    out <- list()
    R <- iratio('U238U235')[1]
    L8 <- lambda('U238')[1]
    L5 <- lambda('U235')[1]
    L2 <- lambda('Th232')[1]
    L7 <- lambda('Sm147')[1]
    f147 <- f147Sm()[1]
    doSm <- doSm(x)
    out$lims <- get.logratioplot.limits(x)
    if (doSm){
        uvw <- UThHe2uvw(x)
        w <- mean(uvw[,'w'],na.rm=TRUE)
        UVW <- colMeans(UThHe2uvw(x),na.rm=TRUE)
        Sm <- exp(UVW[3])/(exp(UVW[1])+exp(UVW[2])+exp(UVW[3])+1)
    } else {
        uv <- UThHe2uv(x)
        w <- 0
        Sm <- 0
        cc <- 0
    }
    if (all(is.finite(xlim))) out$lims[1:2] <- xlim
    if (all(is.finite(ylim))) out$lims[3:4] <- ylim
    du <- out$lims[2]-out$lims[1]
    dv <- out$lims[4]-out$lims[3]
    tticks <- get_logratio_tticks(out$lims,Sm=Sm)
    out$uv <- list()
    out$tticks <- NULL
    nt <- 0
    for (i in 1:length(tticks)){
        tt <- tticks[i]
        aa <- 8*(exp(L8*tt)-1)*R/(1+R) + 7*(exp(L5*tt)-1)/(1+R)
        bb <- 6*(exp(L2*tt)-1)
        if (doSm) cc <- f147*(exp(L7*tt)-1)
        uv.plot <- NULL
        # evaluate the maximum v value
        pred.exp.u <- 1-bb*exp(out$lims[4])-cc*exp(w)
        if (pred.exp.u > 0){
            u4maxv <- log(pred.exp.u) - log(aa)
            if (u4maxv > out$lims[1] && u4maxv < out$lims[2])
                uv.plot <- rbind(uv.plot,c(u4maxv,out$lims[4]))
        }
        # evaluate all the whole range of u values
        for (j in 0:res){
            u <- out$lims[1]+du*j/res
            exp.v <- (1-aa*exp(u)-cc*exp(w))/bb
            if (exp.v > exp(out$lims[3]) & exp.v < exp(out$lims[4])){
                v <- log(exp.v)
                uv.plot <- rbind(uv.plot,c(u,v))
            }
        }
        # evaluate the minimum v value
        pred.exp.u <- 1-bb*exp(out$lims[3])-cc*exp(w)
        if (pred.exp.u > 0){
            u4minv <- log(pred.exp.u) - log(aa)
            if (u4minv > out$lims[1] && u4minv < out$lims[2])
                uv.plot <- rbind(uv.plot,c(u4minv,out$lims[3]))
        }
        # add to the list if any solutions were found
        if (!is.null(uv.plot)){
            nt <- nt + 1
            out$tticks[nt] <- tt
            out$uv[[nt]] <- uv.plot
        }
    }
    out
}
get_helioplot_contours <- function(x,fact=c(1,1,1),res=50){
    out <- list()
    doSm <- doSm(x)
    if (doSm){
        uvw <- UThHe2uvw(x)
        SmU <- exp(mean(uvw[,'w']-uvw[,'u'],na.rm=TRUE))
    }
    out$tticks <- get_helioplot_tticks(fact)
    out$xyz <- list()
    nt <- length(out$tticks)
    for (i in 1:nt){
        tt <- out$tticks[i]
        U <- seq(1,0,length.out=res)
        Th <- seq(0,1,length.out=res)
        if (doSm) Sm <- SmU*U
        else Sm <- 0
        He <- get_He(tt,U,Th,Sm)
        out$xyz[[i]] <- cbind(He,U,Th)
    }
    out
}

# f = the distance from the X- and Y-margins for
# the maximum and minimum age contours to be plotted
get_logratio_tticks <- function(lims,f=0.05,Sm=0){
    uv.lims <- rep(0,4)
    uv.lims[1] <- lims[1]+f*(lims[3]-lims[1])
    uv.lims[2] <- lims[2]+f*(lims[4]-lims[2])
    uv.lims[3] <- lims[3]-f*(lims[3]-lims[1])
    uv.lims[4] <- lims[4]-f*(lims[4]-lims[2])
    Umax <- exp(uv.lims[1])/(exp(uv.lims[1])+exp(uv.lims[3])+1)
    Thmax <- exp(uv.lims[3])/(exp(uv.lims[1])+exp(uv.lims[3])+1)
    Hemax <- 1/(exp(uv.lims[1])+exp(uv.lims[3])+1)
    Umin <- exp(uv.lims[2])/(exp(uv.lims[2])+exp(uv.lims[4])+1)
    Thmin <- exp(uv.lims[4])/(exp(uv.lims[2])+exp(uv.lims[4])+1)
    Hemin <- 1/(exp(uv.lims[2])+exp(uv.lims[4])+1)
    mint <- get_UThHe_age(Umin,0,Thmin,0,Hemin,0,Sm,0)[1]
    maxt <- get_UThHe_age(Umax,0,Thmax,0,Hemax,0,Sm,0)[1]
    get_tticks(mint,maxt)
}
# f = the distance from the ternary corners for
# the maximum and minimum age contours to be plotted
get_helioplot_tticks <- function(fact,f=0.05,Sm=0){
    m <- c(f,1-f,0)/fact
    M <- c(1-f,f,0)/fact
    mint <- get_UThHe_age(m[2],0,m[3],0,m[1],0,Sm,0)[1]
    maxt <- get_UThHe_age(M[2],0,M[3],0,M[1],0,Sm,0)[1]
    get_tticks(mint,maxt)
}
get_tticks <- function(mint,maxt){
    m <- log10(mint)
    M <- log10(maxt)
    grDevices::axisTicks(usr=c(m,M),log=TRUE)
}

get.logratioplot.limits <- function(x,nse=3){
    ns <- length(x)
    doSm <- doSm(x)
    minu <- Inf
    maxu <- -Inf
    minv <- Inf
    maxv <- -Inf
    for (i in 1:ns){
        d <- UThHe2uv_covmat(x,i)
        uv <- d$uv
        uv.err <- sqrt(diag(d$covmat)[c('u','v')])
        umin <- uv['u'] - nse*uv.err['u']
        umax <- uv['u'] + nse*uv.err['u']
        vmin <- uv['v'] - nse*uv.err['v']
        vmax <- uv['v'] + nse*uv.err['v']
        if (umax>maxu) maxu <- umax
        if (umin<minu) minu <- umin
        if (vmax>maxv) maxv <- vmax
        if (vmin<minv) minv <- vmin
    }
    c(minu,maxu,minv,maxv)
}

# x is an object of class UThHe
UThHe2uvw <- function(x){
    if (is.UThHe(x)){
        logHe <- log(x[,'He'])
        u <- log(x[,'U']) - logHe
        v <- log(x[,'Th']) - logHe
        w <- log(x[,'Sm']) - logHe
        out <- cbind(u,v,w)
    } else {
        out <- matrix(log(x[c('U','Th','Sm')])-log(x['He']),1,3)
    }
    colnames(out) <- c('u','v','w')
    out
}
UThHe2uv <- function(x){
    if (is.UThHe(x)){
        logHe <- log(x[,'He'])
        u <- log(x[,'U']) - logHe
        v <- log(x[,'Th']) - logHe
        out <- cbind(u,v)
    } else {
        out <- matrix(log(x[c('U','Th')])-log(x['He']),1,2)
    }
    colnames(out) <- c('u','v')
    out
}

uvw2UThHe <- function(uvw,covmat=matrix(0,3,3)){
    u <- uvw[1]
    v <- uvw[2]
    w <- uvw[3]
    D <- exp(u)+exp(v)+exp(w)+1
    U <- exp(u)/D
    Th <- exp(v)/D
    Sm <- exp(w)/D
    He <- 1/D
    J <- matrix(0,4,3)
    J[1,1] <- exp(u)*(exp(v)+exp(w)+1)/D^2
    J[1,2] <- -exp(u)*exp(v)/D^2
    J[1,3] <- -exp(u)*exp(w)/D^2
    J[2,1] <- -exp(v)*exp(u)/D^2
    J[2,2] <- exp(v)*(exp(u)+exp(w)+1)/D^2
    J[2,3] <- -exp(v)*exp(w)/D^2
    J[3,1] <- -exp(w)*exp(u)/D^2
    J[3,2] <- -exp(w)*exp(v)/D^2
    J[3,3] <- exp(w)*(exp(u)+exp(v)+1)/D^2
    J[4,1] <- -exp(u)/D^2
    J[4,2] <- -exp(v)/D^2
    J[4,3] <- -exp(w)/D^2
    E <- J %*% covmat %*% t(J)
    sU <- sqrt(E[1,1])
    sTh <- sqrt(E[2,2])
    sSm <- sqrt(E[3,3])
    sHe <- sqrt(E[4,4])
    out <- c(U,sU,Th,sTh,Sm,sSm,He,sHe)
    names(out) <- c('U','sU','Th','sTh','Sm','sSm','He','sHe')
    out
}
uv2UThHe <- function(uv,covmat=matrix(0,2,2)){
    u <- uv[1]
    v <- uv[2]
    D <- exp(u)+exp(v)+1
    U <- exp(u)/D
    Th <- exp(v)/D
    He <- 1/D
    J <- matrix(0,3,2)
    J[1,1] <- exp(u)*(exp(v)+1)/D^2
    J[1,2] <- -exp(u)*exp(v)/D^2
    J[2,1] <- -exp(v)*exp(u)/D^2
    J[2,2] <- exp(v)*(exp(u)+1)/D^2
    J[3,1] <- -exp(u)/D^2
    J[3,2] <- -exp(v)/D^2
    E <- J %*% covmat %*% t(J)
    sU <- sqrt(E[1,1])
    sTh <- sqrt(E[2,2])
    sHe <- sqrt(E[3,3])
    out <- c(U,sU,Th,sTh,He,sHe)
    names(out) <- c('U','sU','Th','sTh','He','sHe')
    out
}
uv2HeUTh <- function(uv){
    if (is.matrix(uv)){
        u <- uv[,1]
        v <- uv[,2]
    } else {
        u <- uv[1]
        v <- uv[2]
    }
    D <- exp(u)+exp(v)+1
    U <- exp(u)/D
    Th <- exp(v)/D
    He <- 1/D
    cbind(He,U,Th)
}

UThHe2uvw_covmat <- function(x,i,w=0){
    U <- x[i,'U']
    sU <- x[i,'errU']
    Th <- x[i,'Th']
    sTh <- x[i,'errTh']
    Sm <- x[i,'Sm']
    sSm <- x[i,'errSm']
    He <- x[i,'He']
    sHe <- x[i,'errHe']
    out <- list()
    out$uvw <- UThHe2uvw(x[i,])
    out$covmat <- matrix(0,3,3)
    J <- matrix(0,3,4)
    E <- matrix(0,4,4)
    diag(E) <- c(sU,sTh,sSm,sHe)^2
    J[1,1] <- 1/U   # du.dU
    J[1,4] <- -1/He # du.dHe
    J[2,2] <- 1/Th  # dv.dTh
    J[2,4] <- -1/He # dv.dHe
    J[3,3] <- 1/Sm  # dw.dSm
    J[3,4] <- -1/He # dw.dHe
    out$covmat <- (J %*% E %*% t(J)) + diag(3)*w^2
    names(out$uvw) <- c("u","v","w")
    rownames(out$covmat) <- c("u","v","w")
    colnames(out$covmat) <- c("u","v","w")
    out
}
UThHe2uv_covmat <- function(x,i,w=0){
    out <- list()
    U <- x[i,'U']
    sU <- x[i,'errU']
    Th <- x[i,'Th']
    sTh <- x[i,'errTh']
    He <- x[i,'He']
    sHe <- x[i,'errHe']
    out$uv <- UThHe2uv(x[i,])
    out$covmat <- matrix(0,2,2)
    J <- matrix(0,2,3)
    E <- matrix(0,3,3)
    diag(E) <- c(sU,sTh,sHe)^2
    J[1,1] <- 1/U   # du.dU
    J[1,3] <- -1/He # du.dHe
    J[2,2] <- 1/Th  # dv.dTh
    J[2,3] <- -1/He # dv.dHe
    out$covmat <- (J %*% E %*% t(J)) + diag(2)*w^2
    names(out$uv) <- c("u","v")
    rownames(out$covmat) <- c("u","v")
    colnames(out$covmat) <- c("u","v")
    out
}

# ternary compositions to plot coordinates
xyz2xy <- function(xyz){
    if (is.matrix(xyz)){
        n <- nrow(xyz)
        x <- xyz[,1]
        y <- xyz[,2]
        z <- xyz[,3]
    } else {
        n <- 1
        x <- xyz[1]
        y <- xyz[2]
        z <- xyz[3]
    }
    xy <- matrix(0,n,2)
    xy[,1] <- 0.5*(x+2*z)/(x+y+z)
    xy[,2] <- sin(pi/3)*x/(x+y+z)
    xy
}

renormalise <- function(xyz,fact=c(1,1,1)){
    if (is.matrix(xyz)){
        nr <- nrow(xyz)
        FACT <- matrix(rep(fact,nr),nrow=nr,byrow=TRUE)
        out <- xyz*FACT
        NORM <- rowSums(out)
        out <- out/NORM
    } else {
        out <- xyz*fact/sum(xyz*fact)
    }
    out
}

getfact <- function(x,fit){
    HeUTh <- uv2HeUTh(fit$uvw[1:2])
    fact <- signif(1/HeUTh,1)
}
