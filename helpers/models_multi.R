simulation_multi_latent = function(n0, n1, d, scenario, unif_w = 0.2, transform = FALSE){

  n = n0 + n1

  group_labels = c(rep(0, n0), rep(1, n1))
  n = n0 + n1
  data = matrix(NA, nrow = n, ncol = d)

  indices_0 = which(group_labels == 0)
  indices_1 = which(group_labels == 1)

  size_0 = length(indices_0)
  size_1 = length(indices_1)

  differences_log_dens = numeric(n)

  if(scenario == "latent_location_shift"){

    Q_full = randortho(d)

    m = 4
    U = Q_full[, 1:m]

    dim(U)

    mean_vec0 = c(-0.5,0,0,0)
    mean_vec1 = c(0.5,0,0,0)

    sigma0 = diag(c(1,1,1,1)^2)
    sigma1 = diag(c(1,1,1,1)^2)

    data_u = matrix(NA, nrow = n, ncol = m)

    sd_small = 0.1

    mean_unif = c(0,0,0,0)
    sd_unif = 2
    sigma_unif = diag(rep(sd_unif, m)^2)

    for(i in 1:n){

      prob_mix = runif(1)

      if(group_labels[i] == 0){
        if(prob_mix > unif_w){
          u = t(rmvnorm(1, mean_vec0, sigma0))
        }else{
          u = t(rmvnorm(1, mean_unif, sigma_unif))
        }
      }else{
        if(prob_mix > unif_w){
          u = t(rmvnorm(1, mean_vec1, sigma1))
        }else{
          u = t(rmvnorm(1, mean_unif, sigma_unif))
        }
      }
      data_u[i,] = u
    }

    data_before_transform = data_u %*% t(U) + matrix(rnorm(d*n,mean=0,sd=sd_small),nrow=n,ncol=d)

    # compute the means and covariances of the d-dim distribution
    mean0_vec_true = U %*% mean_vec0
    mean1_vec_true = U %*% mean_vec1
    mean_unif_vec_true = U %*% mean_unif

    sigma0_true = U %*% sigma0 %*% t(U) + sd_small^2 * diag(d)
    sigma1_true = U %*% sigma1 %*% t(U) + sd_small^2 * diag(d)
    sigma_unif_true = U %*% sigma_unif %*% t(U) + sd_small^2 * diag(d)

    # if necessary, transform to change the marginal distributions
    if(transform){
      a_beta = 0.5
      b_beta = 5
      for(j in 1:d){
        data[,j] = qbeta(pnorm(data_before_transform[,j], mean_unif_vec_true[j], sqrt(sigma_unif_true[j,j])),a_beta,b_beta)
      }
    }else{
      data = data_before_transform
    }

    log_dens_0 = log((1-unif_w) * dmvnorm(data_before_transform, mean0_vec_true, sigma0_true) +
                       unif_w * dmvnorm(data_before_transform, mean_unif_vec_true, sigma_unif_true))

    log_dens_1 = log((1-unif_w) * dmvnorm(data_before_transform, mean1_vec_true, sigma1_true) +
                       unif_w * dmvnorm(data_before_transform, mean_unif_vec_true, sigma_unif_true))

    differences_log_dens = log_dens_0 - log_dens_1
  }


  if(scenario == "latent_dispersion"){

    Q_full = randortho(d)

    m = 4
    U = Q_full[, 1:m]

    dim(U)

    mean_vec0 = c(0,0,0,0)
    mean_vec1 = c(0,0,0,0)

    sigma0 = diag(c(1,1,1,1)^2)
    sigma1 = diag(c(0.5,1,1,1)^2)

    data_u = matrix(NA, nrow = n, ncol = m)

    sd_small = 0.1

    mean_unif = c(0,0,0,0)
    sd_unif = 2
    sigma_unif = diag(rep(sd_unif, m)^2)

    for(i in 1:n){

      prob_mix = runif(1)

      if(group_labels[i] == 0){
        if(prob_mix > unif_w){
          u = t(rmvnorm(1, mean_vec0, sigma0))
        }else{
          u = t(rmvnorm(1, mean_unif, sigma_unif))
        }
      }else{
        if(prob_mix > unif_w){
          u = t(rmvnorm(1, mean_vec1, sigma1))
        }else{
          u = t(rmvnorm(1, mean_unif, sigma_unif))
        }
      }
      data_u[i,] = u
    }

    data_before_transform = data_u %*% t(U) + matrix(rnorm(d*n,mean=0,sd=sd_small),nrow=n,ncol=d)

    # compute the means and covariances of the d-dim distribution
    mean0_vec_true = U %*% mean_vec0
    mean1_vec_true = U %*% mean_vec1
    mean_unif_vec_true = U %*% mean_unif

    sigma0_true = U %*% sigma0 %*% t(U) + sd_small^2 * diag(d)
    sigma1_true = U %*% sigma1 %*% t(U) + sd_small^2 * diag(d)
    sigma_unif_true = U %*% sigma_unif %*% t(U) + sd_small^2 * diag(d)

    # if necessary, transform to change the marginal distributions
    if(transform){
      a_beta = 0.5
      b_beta = 5
      for(j in 1:d){
        data[,j] = qbeta(pnorm(data_before_transform[,j], mean_unif_vec_true[j], sqrt(sigma_unif_true[j,j])),a_beta,b_beta)
      }
    }else{
      data = data_before_transform
    }

    log_dens_0 = log((1-unif_w) * dmvnorm(data_before_transform, mean0_vec_true, sigma0_true) +
                       unif_w * dmvnorm(data_before_transform, mean_unif_vec_true, sigma_unif_true))

    log_dens_1 = log((1-unif_w) * dmvnorm(data_before_transform, mean1_vec_true, sigma1_true) +
                       unif_w * dmvnorm(data_before_transform, mean_unif_vec_true, sigma_unif_true))

    differences_log_dens = log_dens_0 - log_dens_1
  }

  out = list("data" = data,
             "data_u" = data_u,
             "group_labels" = group_labels,
             "true_log_w_obs" = differences_log_dens)

  return(out)
}
