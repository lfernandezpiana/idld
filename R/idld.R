#### INTEGRATED DUAL LOCAL DEPTH

#' @name idld_function
#' @title Integrated Dual Local Depth
#'
#' @description Calculates the integrated dual local depth 
#'
#' @param Z data to apply depth. It should be an array of dimension (n,p,l) where l is the number of functional coordinates. 
#' Z[,,i] is a numeric matrix where each row represents a functional observation for i=1,...,l.
#' @param data_mf data on which depth is based. Same format than Z.
#' @param beta locality parameter between 0 and 1
#' @param m number of random projections
#' @param type the data type to apply the idld, "multivariate", "functional" or "multi_functional".
#' @param verbose if TRUE prints the algorithm progress.
#'
#' @return A numeric vector object that contains the depth for each point. 
#' 
#' @export
#' @importFrom dplyr "%>%"



## idld

idld = function(z, data, beta, m, type, verbose) {
  switch(type, 
         multivariate = idld_m(z, data, beta, m, verbose),
         functional = idld_f(z, data, beta, m, verbose),
         multi_functional = idld_mf(z, data, beta, m, verbose)
  )
}