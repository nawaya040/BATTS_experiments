# this code reproduces the 20D experiment provided in Awaya, Xu, and Ma (2025)

library(mvtnorm)
library(balancePM)
library(ada)
library(pracma)

# notice that the original densratio package contains a bug for evaluating density ratios
# the fixed version is available on https://github.com/nawaya040/densratio
# library(densratio)

# read the helper functions
source("./helpers/models_multi.R")
source("./helpers/utilities_for_experiment.R")


# set the parameter values
max_depth = 4

K_CV = 5
num_trees_min = 10
num_trees_max = 1000
num_trees_Bayes = 200

learn_rate = 0.01
n_min_obs_per_node = 1

n_bins = 32

size_burnin = 2000
size_backfitting = 1000

# function to compute the squared error
compute_sq_error = function(log_w_estimated, log_w_true_true){
  out = mean((log_w_estimated - log_w_true_true)^2)
  return(out)
}

balance_to_log_ratio = function(x){
  return(log(x)*2)
}

# input the settings
data_settings = list()

data_settings[[1]] = list("latent_location_shift", 5000, 5000, 0.2, FALSE)
data_settings[[2]] = list("latent_location_shift", 9000, 1000, 0.2, FALSE)
data_settings[[3]] = list("latent_dispersion", 5000, 5000, 0.2, FALSE)
data_settings[[4]] = list("latent_dispersion", 9000, 1000, 0.2, FALSE)

methods = c("Boosting(Gradient)", "Boosting(helinger)", "BART", "Real adaBoost", "KLIEP", "uLSIF")
n_methods = length(methods)

error_group0_vec = numeric(n_methods)
error_group1_vec = numeric(n_methods)
#error_grid_vec = numeric(n_methods)

# in the paper, we run the following program
# for index_settings = 1,2,...,4 and
# index_repeat = 1,2,...,50

index_settings = 1
index_repeat = 1

# simulate the data
data_settings_current = data_settings[[index_settings]]

scenario = data_settings_current[[1]]
n0 = data_settings_current[[2]]
n1 = data_settings_current[[3]]
unif_w = data_settings_current[[4]]
transform = data_settings_current[[5]]

set.seed(index_repeat)
out_data = simulation_multi_latent(n0, n1, d=20, scenario, unif_w, transform)

data = out_data$data
group_labels = out_data$group_labels

indices_group0 = which(group_labels == 0)
indices_group1 = which(group_labels == 1)

log_w_true_data = out_data$true_log_w_obs


########################################################################################################
# method 1: proposed boosting (gradient based)
result_estimation = estimate_balancing_weight_boosting(data = data,
                                                       group_labels = group_labels,
                                                       num_trees = num_trees_max,
                                                       K_CV = K_CV,
                                                       max_resol = max_depth,
                                                       learn_rate = learn_rate,
                                                       n_bins = n_bins,
                                                       use_gradient = T,
                                                       quiet = T
)

# result of the boosting
log_ratio_boosting_grad_data = balance_to_log_ratio(result_estimation$balance_weight_boosting_data)

# compute errors
error_group0_vec[1] = compute_sq_error(log_ratio_boosting_grad_data[indices_group0], log_w_true_data[indices_group0])
error_group1_vec[1] = compute_sq_error(log_ratio_boosting_grad_data[indices_group1], log_w_true_data[indices_group1])

########################################################################################################
# method 2: proposed boosting (hellinger based)
result_estimation = estimate_balancing_weight_boosting(data = data,
                                                       group_labels = group_labels,
                                                       num_trees = num_trees_max,
                                                       K_CV = K_CV,
                                                       max_resol = max_depth,
                                                       learn_rate = learn_rate,
                                                       n_bins = n_bins,
                                                       use_gradient = F,
                                                       quiet = T
)

loss_hel = colMeans(result_estimation$loss_CV_store)
log_ratio_boosting_hel_data = balance_to_log_ratio(result_estimation$balance_weight_boosting_data)

# compute errors
error_group0_vec[2] = compute_sq_error(log_ratio_boosting_hel_data[indices_group0], log_w_true_data[indices_group0])
error_group1_vec[2] = compute_sq_error(log_ratio_boosting_hel_data[indices_group1], log_w_true_data[indices_group1])


########################################################################################################
# method 3: proposed back-fitting
result_estimation = estimate_balancing_weight_Bayes(data = data,
                                                    group_labels = group_labels,
                                                    num_trees = num_trees_Bayes,
                                                    n_bins = n_bins,
                                                    size_burnin = size_burnin,
                                                    size_backfitting = size_backfitting,
                                                    output_BART_ensembles = FALSE,
                                                    quiet = T,
                                                    update_lambda = FALSE
)

# result of the BART
log_ratio_BART_data = balance_to_log_ratio(result_estimation$balance_weight_BART_data)
log_ratio_BART_data_mean = rowMeans(log_ratio_BART_data)

error_group0_vec[3] = compute_sq_error(log_ratio_BART_data_mean[indices_group0], log_w_true_data[indices_group0])
error_group1_vec[3] = compute_sq_error(log_ratio_BART_data_mean[indices_group1], log_w_true_data[indices_group1])

########################################################################################################
# method 4: real AdaBoost

result_estimation = ada_log_weight(data,
                                   group_labels,
                                   NULL,
                                   K_CV = K_CV,
                                   num_trees_max = num_trees_max,
                                   num_trees_min = num_trees_min,
                                   maxdepth = max_depth,
                                   learn_rate = learn_rate
)

log_ratio_ada_data = result_estimation$log_w_hat_ada_data

# compute errors
error_group0_vec[4] = compute_sq_error(log_ratio_ada_data[indices_group0], log_w_true_data[indices_group0])
error_group1_vec[4] = compute_sq_error(log_ratio_ada_data[indices_group1], log_w_true_data[indices_group1])

########################################################################################################
# method 5 and 6: densratio methods
X0 = data[which(group_labels==0),]
X1 = data[which(group_labels==1),]

# perform KLIEP
result_estimation = KLIEP(X0, X1, verbose = FALSE)
log_ratio_KLIEP_data = log(result_estimation$compute_density_ratio(data))

# compute errors
error_group0_vec[5] = compute_sq_error(log_ratio_KLIEP_data[indices_group0], log_w_true_data[indices_group0])
error_group1_vec[5] = compute_sq_error(log_ratio_KLIEP_data[indices_group1], log_w_true_data[indices_group1])

# perform uLSIF
result_estimation = uLSIF(X0, X1, verbose = FALSE)
log_ratio_uLSIF_data = log(result_estimation$compute_density_ratio(data))

# compute errors
error_group0_vec[6] = compute_sq_error(log_ratio_uLSIF_data[indices_group0], log_w_true_data[indices_group0])
error_group1_vec[6] = compute_sq_error(log_ratio_uLSIF_data[indices_group1], log_w_true_data[indices_group1])

########################################################################################################
# show the symmetric MSE
(error_group0_vec + error_group1_vec) / 2
