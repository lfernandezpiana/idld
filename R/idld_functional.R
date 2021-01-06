#### CALCULA LA PROFUNDIDAD DUAL INTEGRADA LOCAL VERSION c++

#' @name idld_f
#' @title Functional Integrated Dual Local Depth
#'
#' @description Calculates the integrated dual local depth for functional data and multivariate functional data.
#'
#' @param Z data to apply depth. It should be a numeric matrix where each row represents an observation.
#' @param data_f data on which depth is based. It should be a numeric matrix where each row represents an observation.
#' @param beta locality parameter between 0 and 1
#' @param m number of random projections
#' @param verbose if TRUE prints the algorithm progress.
#'
#' @return A numeric vector object that contains the depth for each point. 
#' 
#' @examples
#' data(growth_domestic_product)
#' data(inflation_rate)
#' m_functional_data = as.matrix(cbind(growth_domestic_product[,-1]/400,inflation_rate[,-1]))
#' local_depth = idld_f(m_functional_data, m_functional_data, 0.3, 500, TRUE)
#' @export
#' @importFrom dplyr "%>%"

## Local Integrated Dual Depth

## Z (numeric matrix): data to apply depth. Z must be a numeric matrix where  
# each row represents an observation.
## datos (numeric matrix): data on which the depth is based. 
# datos must be a numeric matrix where each row represents an observation.
## beta (float): locality parameter, beta range in (0,1]. For beta = 1 lidd_cpp 
# returns global Integrated Dual Depth.
## m (int): number of random projections.
## verbose (bool): if TRUE prints the algorithm progress.

require(sde)

idld_f = function(Z,data_f,beta,m,verbose)
  {
    n = nrow(data_f) # Cantidad de muestras
    p = ncol(data_f) # Cantidad de variables
    q = nrow(Z) # Cantidad de datos a quienes hay que calcularles la profundidad
## GENERO LAS PROYECCIONES
    Q=matrix(ncol=p,nrow=m)
    if (verbose==TRUE) print("Generating random projections")
    for (i in 1:m)
    {
      Q[i,] = BM(x=0, t0=0, T=1, N=(p-1))
    } 
### PROYECTO LOS DATOS 
    Proyecciones = eigenMapMatMult(data_f,t(Q))
### CALCULO LAS UNIVARIDAS
    if (verbose==TRUE) print("Procesing univariated depths")
    u = numeric(q)
    W = eigenMapMatMult(Q,t(Z))
    for (l in 1:q)
        u[l] = cppLdaux(as.numeric(W[,l]),Proyecciones,beta)
    return(u)
  }

    


