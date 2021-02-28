#### CALCULA LOS CLASTERS BASADOS EN LAS REGIONES DE MAYOR PROFUNDIDAD LOCAL
#### EL PARAMETRO ALPHA LO ELIGE EN FORMA ADAPTIVA A LOS DATOS.
require(class)
require(mclust)
require(cluster)
# Load required c++ rutines
# Load lidd cluster function


## Auxiliary function applyed to each data type

idld = function(z, data, beta, m, type, verbose) {
  switch(type, 
         multivariate = idld_m(z, data, beta, m, verbose),
         functional = idld_f(z, data, beta, m, verbose),
         multi_functional = idld_mf(z, data, beta, m, verbose)
         )
}

idld_cluster = function(Z, data, beta, m, alpha_quantile, K, type, verbose=FALSE) {
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
  
  
  
  #cluster_orig[deep_region_selector,i] = kmeans(z[deep_region_selector,],k, nstart=10)$cluster
  # Agrupo por distancia al mas cercano del nucleo
  #cluster_orig[-c(deep_region_selector),i] = knn1(train = z[deep_region_selector,],
   #                                               test = z[-c(deep_region_selector),],
    #                                              cl = cluster_orig[deep_region_selector])

matrix(c(TRUE,TRUE,FALSE,TRUE),2,2)
