#### CALCULA LA PROFUNDIDAD DUAL INTEGRADA LOCAL VERSION c++

#' @name idld_m
#' @title Multivariate Integrated Dual Local Depth
#'
#' @description Calculates the integrated dual local depth for multivariate data
#'
#' @param Z data to apply depth. It should be a numeric matrix where each row represents an observation.
#' @param data_m data on which depth is based. It should be a numeric matrix where each row represents an observation.
#' @param beta locality parameter between 0 and 1
#' @param m number of random projections
#' @param verbose if TRUE prints the algorithm progress.
#'
#' @return A numeric vector that contains the depth for each point. 
#' 
#' @examples
#' library(mvnfast)
#' X = rmvn(90,  c(-2,-2), sigma=diag(rep(1,2)))
#' Y = rmvn(110, c(2,2), sigma=diag(rep(1,2)))
#' Z = rmvn(150, c(4,-4), sigma = rbind(c(2,0.8),c(0.8,1)))
#' W = rbind(X,Y,Z)
#' local_depth = idld_m(W, W ,0.3, 500, TRUE)
#' plot(W, pch=20)
#' local_depth_center_region = which(local_depth>quantile(local_depth,0.9))
#' points(W[local_depth_center_region,], col="blue", pch=20)
#' 
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom mvnfast rmvn


## Integrated Local Dual Depth

## Z (numeric matrix): data to apply depth. Z must be a numeric matrix where  
# each row represents an observation.
## datos (numeric matrix): data on which the depth is based. 
# datos must be a numeric matrix where each row represents an observation.
## beta (float): locality parameter, beta range in (0,1]. For beta = 1 lidd_cpp 
# returns global Integrated Dual Depth.
## m (int): number of random projections.
## verbose (bool): if TRUE prints the algorithm progress.

require(mvnfast)

norma = function(t)
{
  u = sqrt(sum(t^2))
  return(u)
}

idld_m = function(Z,data_m,beta,m,verbose=TRUE)
  {
    n = nrow(data_m) # Cantidad de muestras
    p = ncol(data_m) # Cantidad de variables
    q = nrow(Z) # Cantidad de datos a quienes hay que calcularles la profundidad
## GENERO LAS PROYECCIONES
    Q = rmvn(m,rep(0,p),diag(1,nrow=p,ncol=p))
    if (verbose==TRUE) print("Generating random projections")
    for (i in 1:m)
    {
      Q[i,] = Q[i,]/norma(Q[i,])
    } 
### PROYECTO LOS DATOS 
    Proyecciones = eigenMapMatMult(data_m,t(Q))
    #Proyecciones = datos%*%t(Q)
### CALCULO LAS UNIVARIDAS
    if (verbose==TRUE) print("Procesing univariated depths")
    u = numeric(q)
    W = eigenMapMatMult(Q,t(Z))
    for (l in 1:q)
        u[l] = cppLdaux(as.numeric(W[,l]),Proyecciones,beta)
    return(u)
  }

    


