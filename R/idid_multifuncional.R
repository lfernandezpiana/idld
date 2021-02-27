#### CALCULA LA PROFUNDIDAD DUAL INTEGRADA LOCAL VERSION c++

#' @name idld_mf
#' @title Functional Integrated Dual Local Depth
#'
#' @description Calculates the integrated dual local depth for functional data and multivariate functional data.
#'
#' @param Z data to apply depth. It should be an array of dimension (n,p,l) where l is the number of functional coordinates. 
#' Z[,,i] is a numeric matrix where each row represents a functional observation for i=1,...,l.
#' @param data_f data on which depth is based. Same format than Z.
#' @param beta locality parameter between 0 and 1
#' @param m number of random projections
#' @param verbose if TRUE prints the algorithm progress.
#'
#' @return A numeric vector object that contains the depth for each point. 
#' 
#' @examples
#' data(growth_domestic_product)
#' data(inflation_rate)
#' library(abind)
#' m_functional_data = abind(as.matrix(growth_domestic_product[,-1]/400), 
#'                           as.matrix(inflation_rate[,-1]), 
#'                           along = 3)
#' local_depth = idld_mf(m_functional_data, m_functional_data, 0.3, 500, TRUE)
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom sde BM

## Local Integrated Dual Depth

## Z (numeric array): data to apply depth. Z must be a numeric array of dimension (n,p,l) where  
# l is the number of fuctional coordinates. Z[,,i] is a matrix where each row represents an observation for i=1,...,l.
## data_f (numeric array): data on which the depth is based. Same format than Z.
## beta (float): locality parameter, beta range in (0,1]. For beta = 1 idld_mf returns global Integrated Dual Depth.
## m (int): number of random projections.
## verbose (bool): if TRUE prints the algorithm progress.

require(sde)

idld_mf = function(Z,data_mf,beta,m,verbose)
{
  d_mf = dim(data_mf)
  d_Z = dim(Z)
  
  ## GENERO LAS PROYECCIONES
  Q = array(0, dim=c(m, d_mf[2], d_mf[3]))  
  if (verbose==TRUE) print("Generating random projections")
  for (l in 1:d_mf[3]) {
    for (i in 1:m)
    {
      Q[i,,l] = BM(x=0, t0=0, T=1, N=(d_mf[2]-1))
    }
  }
  
  ### PROYECTO LOS DATOS
  Proyecciones = eigenMapMatMult(data_mf[,,1],t(Q[,,1]))
  W = eigenMapMatMult(Q[,,1],t(Z[,,1]))
  for (l in 2:d_mf[3]) {
    Proyecciones = Proyecciones + eigenMapMatMult(data_mf[,,l],t(Q[,,l]))
    W = W + eigenMapMatMult(Q[,,l],t(Z[,,l]))
  }                 
  ### CALCULO LAS UNIVARIDAS
  if (verbose==TRUE) print("Procesing univariated depths")
  u = numeric(d_Z[1])
  for (l in 1:d_Z[1])
    u[l] = cppLdaux(as.numeric(W[,l]),Proyecciones,beta)
  return(u)
}
