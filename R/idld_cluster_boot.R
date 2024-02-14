#### INTEGRATED DUAL LOCAL DEPTH CLUSTER METHOD

#' @name idld_cluster_boot
#' @title IDLD Clustering with Bootstrap
#'
#' @description It is partition-based clustering technique based on local depth and distance measurement applied to data.
#' This functional selects the optimal alpha_quantile for the procedure.  
#'
#' @param Z data to apply depth. It should be an array of dimension (n,p,l) where l is the number of functional coordinates. 
#' Z[,,i] is a numeric matrix where each row represents a functional observation for i=1,...,l.
#' @param beta locality parameter between 0 and 1
#' @param m number of random projections
#' @param K number of clusters
#' @param B number of bootstrap samples
#' @param type the data type to apply the idld, "multivariate", "functional" or "multi_functional".
#' @param verbose if TRUE prints the algorithm progress.
#'
#' @return returns a list with the following components: 
#' \itemize { 
#' \item local_depth a numeric vector object that contains the depth for each point.
#' \item region vector of booleans indicating which data points is in the central region. 
#' \item clusters numeric vector with the partition.
#' }
#' 
#' @examples
#' library(funHDDC)
#' library(abind)
#' data("triangle")
#' triangle_data = abind(triangle[,1:101], triangle[,102:202], along=3)
#' d = dim(triangle_data)
#' triang_cl = idld_cluster_boot(triangle_data, 0.2, 100, 3, 20, "multi_functional", TRUE)
#' par(mfrow=c(1,2))
#' plot(triangle_data[1,,1], type="n", ylim=c(0,8))
#' for (i in 1:d[1]) lines(triangle_data[i,,1], col=triang_cl$clusters[i])
#' plot(triangle_data[1,,2], type="n", ylim=c(0,8))
#' for (i in 1:d[1]) lines(triangle_data[i,,2], col=triang_cl$clusters[i])
#' 
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom sde BM
#' @importFrom class knn1
#' @importFrom cluster pam
#' @importFrom mclust adjustedRandIndex


## idid cluster boot



idld_cluster_boot = function(Z, beta, m, K, B, type, verbose=FALSE) {
  #INPUT
  #Z: data that will be assigned local depth. Matrix each row is an observation, for bivariate functional data put one coordinate next to the other, both coordinates are assumed to have the same length.
  #datos: data set regarding which the local depth is measured. Matrix each row is an observation, for multifunctional data put one coordinate next to the other.
  #beta: locality parameter.
  #m=number or random directions where the data is projected.
  #alpha_quantile: grid of values to select de deepest region  (Ojo Lucas que vos siempre pensas en q y esto es alfa, es decir 1-q).
  #k: number of clusters.
  #B: number of boostrap samples, is a_grid has only one value put B=0.
  #verbose = if TRUE prints the progress 
  #OUTPUT
  #a.optimal= alpha optimal value
  #clusters_optimal: clustering allocation.
  alpha_quantile = seq(0.5,0.9,0.1)
  n_grid = length(alpha_quantile)
  aris = matrix(nrow=B, ncol=n_grid)
  idld_original = idld_cluster(Z, beta, m, alpha_quantile, K, type, verbose=FALSE)
  cluster_original = idld_original$clusters
  ## Generating bootstrap samples
  if (type == "multi_functional") {
    for (b in 1:B) {
      b_ind = sample(1:dim(Z)[1], dim(Z)[1], replace=TRUE) # boostrap sample indices
      Z_b = Z[b_ind,,] # bootstrap sample
      c_b = idld_cluster(Z_b, beta, m, alpha_quantile, K, type, verbose=FALSE)$clusters
      # compare only the same observations
      c_b = c_b[unique(b_ind),]
      c_o = cluster_original[b_ind,]
      c_o = c_o[unique(b_ind),]
      # Rand index
      for (j in 1:n_grid) {
        aris[b,j]=adjustedRandIndex(c_o[,j],c_b[,j])
      }
    }
  } else {
    for (b in 1:B) {
      b_ind = sample(1:dim(Z)[1], dim(Z)[1], replace=TRUE) # boostrap sample indices
      Z_b = Z[b_ind,] # bootstrap sample
      c_b = idld_cluster(Z_b, beta, m, alpha_quantile, K, type, verbose=FALSE)$clusters
      # compare only the same observations
      c_b = c_b[unique(b_ind),]
      c_o = cluster_original[b_ind,]
      c_o = c_o[unique(b_ind),]
      # Rand index
      for (j in 1:n_grid) {
        aris[b,j]=adjustedRandIndex(c_o[,j],c_b[,j])
      }
    }
  }
  mean_aris = apply(aris,2,mean)
  opt_alpha = which.max(mean_aris)
  salida = list(idld_original$depth, idld_original$region[,opt_alpha], 
                cluster_original[,opt_alpha], opt_alpha, mean_aris)
  names(salida) = c("depth","region","clusters", "alpha", "aris")
  return(salida)
}
