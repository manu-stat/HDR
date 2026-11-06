#######################
# KDE
#######################
KDE = function(data, H = NULL){
  data = as.matrix(data); n = dim(data)[1]; d = dim(data)[2]
  if (d > 6) stop("Sample data dimension must be less than 6")
  if (is.null(H)) H = ks::Hpi(data)
  ks::kde(data, eval.points = data, H = H)$estimate
}



#######################
# KNN
#######################
KNN <- function (data, k = NULL) {
  data = as.matrix(data); n = nrow(data)[1]; d = dim(data)[2]
  if (d > 6) stop("Sample data dimension must be less than 6")
  
  if (is.null(k)) k = round(sqrt(n * 0.5))
  if (length(k) != 1 || !is.finite(k) || k != floor(k) || k < 1) stop("Parameter k must be a positive integer")
  if (k >= n) {warning("Parameter k should be smaller than the sample size n. \n The maximum value corresponding to k = n-1 neighbors is taken."); k = n-1}
  
  apply(data, 1, function(x){
    R = sweep(data, 2, x, FUN = "-")
    d = sort(sqrt(rowSums(R^2)))[-1]
    1/sum(d[seq_len(k)])
  })
}



#######################
# CDF KNN
#######################
CDF_KNN = function(data, k = 30, CDF = NULL){
  data = as.matrix(data); n = nrow(data)[1]; d = dim(data)[2]
  if (dim(data)[2] > 6) stop("Sample data sn must be tridimensional, with the three dimensions on the columns")
  
  if (is.null(k)) k = round(sqrt(n * 0.5))
  if (length(k) != 1 || !is.finite(k) || k != floor(k) || k < 1) stop("Parameter k must be a positive integer")
  if (k >= n) {warning("Parameter k should be smaller than the sample size n. \n The maximum value corresponding to k = n-1 neighbors is taken."); k = n-1}
  
   if (is.null(CDF)) CDF = apply(data, 2, ecdf)# | is.na(sum(colSums(CDF)))) CDF = apply(sn, 2, ecdf)
  
  apply(data, 1, function(x){
    R = sweep(data, 2, x, FUN = "-")
    d = sqrt(rowSums(R^2))
    
    idx = order(d)[1:(k + 1)]
    cdf = sapply(seq_len(length(CDF)), function(j) CDF[[j]](data[idx,j]))
    cdf_x = sapply(seq_len(length(CDF)), function(j) CDF[[j]](x[j]))
    cdf_r = sweep(cdf, 2, cdf_x, FUN = "-")
    cdf_d = sqrt(rowSums(cdf_r^2))
    
    1/sum(d[idx[-1]] / cdf_d[-1])
  })
}



#######################
# CDF e
#######################
CDF_e = function(data, eps = NULL){
  data = as.matrix(data); n = nrow(data)[1]; d = dim(data)[2]
  if (dim(data)[2] > 6) stop("Sample data sn must be tridimensional, with the three dimensions on the columns")
  
  if (is.null(eps) == T) eps = rep(exp(2.13 - 0.3 * log(n)), d)
  if (any(eps <= 0)) stop("Parameter eps must be positive")
  if (length(eps) != d) warning("Parameter eps should be smaller than the sample size n. \n The first element is used as the half-length of the square neighborhood.")
  
  Sign = expand.grid(rep(list(c(1, -1)), d))
  weight = apply(Sign, 1, prod)
  
  V = lapply(1:n, function(i) sweep(Sign*eps, 2, data[i,], "+"))
  
  cdf = sapply(seq_along(V), function(i) apply(V[[i]], 1, function(v) mean(rowSums(sweep(data, 2, v, FUN = "<=")) == d) ))
  colSums(cdf*weight)
}



#######################
# NP Copula
#######################
NP_Cop = function(data){
  data = as.matrix(data); n = nrow(data)[1]; d = dim(data)[2]
  if(d > 6) stop("Sample data sn must be tridimensional, with the three dimensions on the columns")
  
  fit_x = lapply(1:d, function(i) kde1d::kde1d(data[,i]))
  d_x = sapply(1:d, function(i) kde1d::dkde1d(data[,i], obj = fit_x[[i]]))
  
  U = copula::pobs(data)
  fit_cop = kdevine::kdevine(U)
  d_cop = kdevine::dkdevine(U, fit_cop)
  
  d_cop * apply(d_x, 1, prod)
}



#######################
# P Copula
#######################
fit_margin_cop = function(sn, margin = "norm", cop_family = 1:6, ...){
  sn = as.matrix(sn); n = dim(sn)[1]; d = dim(sn)[2]
  
  if(dim(sn)[2] > 6) stop("Sample data sn must be tridimensional, with the three dimensions on the columns")
  if(length(margin) != d) margin = rep(margin, d)
  
  param = lapply(1:d, function(j) if(margin[j] == "norm"){c(mean = mean(sn[,j]), sd = sqrt(mean((sn[,j] - mean(sn[,j]))^2)))} else {c(rate = 1 / mean(sn[,j]))})
  f = sapply(1:d, function(j) do.call(get(paste0("d", margin[j]), mode = "function"), c(list(x = sn[,j]), as.list(param[[j]]))))
    
  U = sapply(1:d, function(j) do.call(get(paste0("p", margin[j]), mode = "function"), c(list(q = sn[,j]), as.list(param[[j]]))))
  U = pmax(pmin(U, 1 - 1e-10), 1e-10)
  cop = VineCopula::RVineStructureSelect(U, familyset = cop_family)
  
  return(list(f_marg = f, pobs = U, copula = cop, param = param, margin = margin))
}
P_Cop = function(data, margin = "norm", ...){
  data = as.matrix(data); n = nrow(data)[1]; d = dim(data)[2]
  if(d > 6) stop("Sample data sn must be tridimensional, with the three dimensions on the columns")
  
  param = fit_margin_cop(data, margin, ...)
  
  d_cop = VineCopula::RVinePDF(param$pobs, param$copula)
  
  d_cop * apply(param$f_marg, 1, prod)
}




#######################
# CDF P Copula
#######################
CDF_P_Cop = function(data, margin = "norm", eps = NULL, ...){
  data = as.matrix(data); n = nrow(data)[1]; d = dim(data)[2]
  if (d > 6) stop("Sample data sn must be tridimensional, with the three dimensions on the columns")
  
  if (is.null(eps) == T) eps = rep(exp(2.13 - 0.3 * log(n)), d)
  if (any(eps <= 0)) stop("Parameter eps must be positive")
  if (length(eps) != d) warning("Parameter eps should be smaller than the sample size n. \n The first element is used as the half-length of the square neighborhood.")
  
  param = fit_margin_cop(data, margin, ...)
  Sign = expand.grid(rep(list(c(1, -1)), d))
  weights = apply(Sign, 1, prod)
  apply(data, 1, function(x){
    V = sweep(Sign * eps, 2, x, FUN = "+")
    cdf = sapply(1:d, function(j) do.call(get(paste0("p", param$margin[j]), mode = "function"), c(list(q = V[,j]), as.list(param$param[[j]]))))
    cdf_cop = apply(cdf, 1, function(x) VineCopula::RVineCDF(x, param$copula))
    sum(cdf_cop*weights)
  })
}