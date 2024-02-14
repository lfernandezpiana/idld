#### INTEGRATED DUAL LOCAL DEPTH CLUSTER METHOD

#' @name idld_cluster
#' @title IDLD Clustering
#'
#' @description It is partition-based clustering technique based on local depth and distance measurement applied to data 
#'
#' @param Z data to apply depth. It should be an array of dimension (n,p,l) where l is the number of functional coordinates. 
#' Z[,,i] is a numeric matrix where each row represents a functional observation for i=1,...,l.
#' @param data_mf data on which depth is based. Same format than Z.
#' @param beta locality parameter between 0 and 1
#' @param m number of random projections
#' @param alpha_quantile proportions of data points to include in the deepest regions. It could be a numeric vector.  
#' @param K number of clusters
#' @param type the data type to apply the idld, "multivariate", "functional" or "multi_functional".
#' @param verbose if TRUE prints the algorithm progress.
#'
#' @return returns a list with the following components:
#' \itemize{  
#' \item local_depth: A numeric vector object that contains the depth for each point.
#' \item region: A matrix containing, in each column, the data which is in the central region related to alpha_quantile selected. 
#' \item clusters a matrix containing, in each column, the data partition related to alpha_quantile selected.
#' }
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom sde BM
#' @importFrom class knn1
#' @importFrom cluster pam

require(class)
require(cluster)

idld_cluster = function(Z, beta, m, alpha_quantile, K, type, verbose=FALSE) {
  #INPUT
  #Z: data that will be assigned local depth. Matrix each row is an observation, for bivariate functional data put one coordinate next to the other, both coordinates are assumed to have the same length.
  #datos: data set regarding which the local depth is measured. Matrix each row is an observation, for multifunctional data put one coordinate next to the other.
  #beta: locality parameter.
  #m: number or random directions where the data is projected.
  #alpha_quantile: grid of values to select de deepest region  (Ojo Lucas que vos siempre pensas en q y esto es alfa, es decir 1-q).
  #k: number of clusters.
  #verbose = if TRUE prints the progress 
  #OUTPUT
  #a.optimal= alpha optimal value
  #clusters_optimal: clustering allocation.
  
  n_grid = length(alpha_quantile) 
  d_Z = dim(Z)
  cluster_orig = matrix(nrow=d_Z[1], ncol=n_grid)
  local_depth = idld(Z, Z, beta, m, type, verbose)
  region = matrix(nrow=d_Z[1], ncol=n_grid)
  
  # Clustering 
  if (type=="multi_functional") {
    Z_new = as.matrix(Z[,,1])
    for (l in 2:d_Z[3]) Z_new = cbind(Z_new, Z[,,l])
    for (q in 1:n_grid) {
      # Clustering deepest region
      region[,q] = local_depth > quantile(local_depth, alpha_quantile[q])
      d = dist(Z_new[region[,q],])
      #for (l in 2:d_Z[3]) d = d + dist(Z[region[,q],,l])
      cluster_orig[region[,q],q] = pam(d, K, cluster.only=TRUE)
      # Clustering the remaining obs 
      cluster_orig[which(region[,q]==FALSE),q] = knn1(train = Z_new[region[,q],],
                                         test = Z_new[which(region[,q]==FALSE),],
                                         cl = cluster_orig[region[,q],q])
    }
  } else {
    for (q in 1:n_grid) {
      region[,q] = local_depth > quantile(local_depth, alpha_quantile[q])
      d = dist(Z[region[,q],])
      cluster_orig[region[,q],q] = pam(d, K, cluster.only=TRUE)
      # Clustering the remaining obs
      print(length(cluster_orig[which(region[,q]==FALSE),q]))
      cluster_orig[which(region[,q]==FALSE),q] = knn1(train = Z[region[,q],],
                                        test = Z[which(region[,q]==FALSE),],
                                        cl = cluster_orig[region[,q],q])
    }
  }
  salida = list(local_depth, region, cluster_orig)
  names(salida) = c("depth", "region", "clusters")
  return(salida)
}
  
  
# idld_cluster_boot = function(Z, beta, m, K, B, type, verbose=FALSE) {
#   #INPUT
#   #Z: data that will be assigned local depth. Matrix each row is an observation, for bivariate functional data put one coordinate next to the other, both coordinates are assumed to have the same length.
#   #datos: data set regarding which the local depth is measured. Matrix each row is an observation, for multifunctional data put one coordinate next to the other.
#   #beta: locality parameter.
#   #m=number or random directions where the data is projected.
#   #alpha_quantile: grid of values to select de deepest region  (Ojo Lucas que vos siempre pensas en q y esto es alfa, es decir 1-q).
#   #k: number of clusters.
#   #B: number of boostrap samples, is a_grid has only one value put B=0.
#   #verbose = if TRUE prints the progress 
#   #OUTPUT
#   #a.optimal= alpha optimal value
#   #clusters_optimal: clustering allocation.
#   alpha_quantile = seq(0.5,0.9,0.1)
#   n_grid = length(alpha_quantile)
#   aris = matrix(nrow=B, ncol=n_grid)
#   cluster_original = idld_cluster(Z, beta, m, alpha_quantile, K, type, verbose=FALSE)$clusters
#   ## Generating bootstrap samples
#   if (type == "multi_functional") {
#     for (b in 1:B) {
#       b_ind = sample(1:dim(Z)[1], dim(Z)[1], replace=TRUE) # boostrap sample indices
#       Z_b = Z[b_ind,,] # bootstrap sample
#       c_b = idld_cluster(Z_b, beta, m, alpha_quantile, K, type, verbose=FALSE)$clusters
#       # compare only the same observations
#       c_b = c_b[unique(b_ind),]
#       c_o = cluster_original[b_ind,]
#       c_o = c_o[unique(b_ind),]
#       # Rand index
#       for (j in 1:n_grid) {
#         aris[b,j]=adjustedRandIndex(c_o[,j],c_b[,j])
#       }
#     }
#   } else {
#     for (b in 1:B) {
#       b_ind = sample(1:dim(Z)[1], dim(Z)[1], replace=TRUE) # boostrap sample indices
#       Z_b = Z[b_ind,] # bootstrap sample
#       c_b = idld_cluster(Z_b, beta, m, alpha_quantile, K, type, verbose=FALSE)$clusters
#       # compare only the same observations
#       c_b = c_b[unique(b_ind),]
#       c_o = cluster_original[b_ind,]
#       c_o = c_o[unique(b_ind),]
#       # Rand index
#       for (j in 1:n_grid) {
#         aris[b,j]=adjustedRandIndex(c_o[,j],c_b[,j])
#       }
#     }
#   }
#   mean_aris = apply(aris,2,mean)
#   opt_alpha = which.max(mean_aris)
#   salida = list(aris, mean_aris, cluster_original[,opt_alpha])
#   names(salida) = c("aris","mean_aris","clusters_optimal")
#   return(salida)
# }
  

  
