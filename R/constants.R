#' Decay constants
#'
#' Returns the decay constants of radioactive istopes
#'
#' @param nuclide the nuclide name
#' @return the decay constant  [in Ma-1]
#' @examples
#' lambda('U238')
#' @export
lambda <- function(nuclide){
    if (nuclide == 'U238') return(1.55125e-4)
    if (nuclide == 'U235') return(9.8485e-4)
    if (nuclide == 'Th232') return(4.9475e-5)
}

#' Decay constant uncertainties
#'
#' Returns the standard errors of radioactive decay constants
#'
#' @param nuclide the nuclide name
#' @return the analytical uncertainty of the decay constant  [in Ma-1]
#' @examples
#' lambda.err('U238')
#' @export
lambda.err <- function(nuclide){
    if (nuclide == 'U238') return(lambda(nuclide)*0.0008)
    if (nuclide == 'U235') return(lambda(nuclide)*0.0010)
    if (nuclide == 'Th232') return(0)
}

#' Isotope abundance
#'
#' Returns the natural abundance of isotopes
#'
#' @param nuclide the nuclide name
#' @return a number between 0 (absent) and 1 (dominant)
#' @examples
#' I.A('U238')
#' @export
I.A <- function(nuclide){
    if (nuclide == 'U238') return(137.818/138.818)
    if (nuclide == 'U235') return(1/138.818)
    if (nuclide == 'U232') return(1)
}
