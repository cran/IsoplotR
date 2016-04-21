#' Read U-Pb data
#'
#' Create a 'UPb' object from a .csv file with U-Pb data
#'
#' @param fname file name (.csv format)
#' @param format formatting option, one of either:
#'
#' 1: 'Pb207Pb206','errPb207Pb206', 'Pb206U238','errPb206U238', 'Pb207U235','errPb207U235'
#' @param header Boolean flag indicating whether the file contains a header
#' @examples
#' # load one of the built-in .csv files:
#' fname <- system.file("UPb.csv",package="IsoplotR")
#' UPb <- read.UPb(fname)
#' concordia.plot(UPb)
#' @export
read.UPb <- function(fname,format=1,header=TRUE){
    x <- utils::read.csv(fname,header=header)
    as.UPb(x,format)
}

#' Create object of class 'UPb'
#'
#' Create a 'UPb' object from a matrix with U-Pb data
#'
#' @param x matrix with U-Pb data
#' @param format formatting option
#'
#' 1: 'Pb207Pb206','errPb207Pb206', 'Pb206U238','errPb206U238', 'Pb207U235','errPb207U235'
#'
#' @examples
#' fname <- system.file("UPb.csv",package="IsoplotR")
#' dat <- read.csv(fname,row.names=1,header=TRUE)
#' concordia.plot(as.UPb(dat))
#' @export
as.UPb <- function(x,format=1){
    out <- list()
    out$format <- format
    nc <- ncol(x)
    nr <- nrow(x)
    if (format == 1){
        if (nc == 6){
            out$x <- as.matrix(x)
            rownames(out$x) <- paste0('X',1:nr)
        } else if (nc == 7) {
            out$x <- as.matrix(x[,2:nc])
            rownames(out$x) <- x[,1]
        } else {
            e <- simpleError('Incorrect number of columns selected')
            stop(e)
        }
        colnames(out$x) <- c('Pb207Pb206','errPb207Pb206',
                             'Pb206U238','errPb206U238',
                             'Pb207U235','errPb207U235')
    }
    class(out) <- "UPb"
    out
}

#' Create a geochronological data object
#'
#' Cast a matrix of data into one of IsoplotR's input formats
#'
#' @param x matrix with U-Pb data
#' @param method one of 'U-Pb', 'Ar-Ar', 'Rb-Sr', 'Sm-Nd', 'Re-Os', 'U-Th-He'
#' 'fission tracks', 'cosmogenic nuclides' or 'other'
#' @param format formatting option [integer]
#' @examples
#' fname <- system.file("UPb.csv",package="IsoplotR")
#' dat <- read.csv(fname,row.names=1,header=TRUE)
#' UPb <- read.matrix(dat,method='U-Pb',format=1)
#' concordia.plot(UPb)
#' @export
read.matrix <- function(x,method,format){
    if (identical(method,"U-Pb")){
        out <- as.UPb(x,format)
    } else {
        out <- NULL
    }
    out
}
