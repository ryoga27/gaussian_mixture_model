library(mvtnorm)

set_initial_mu = function(x, K){
    n = nrow(x)
    d = ncol(x)
    mu = matrix(NA, nrow = d, ncol = K)
    x_sample = sample(x = 1:n, size = K, replace = FALSE)
    for(k in 1:K){
        mu[, k] = x[x_sample[k], ]
    }
    return(mu)
}

rmulti = function(p, n = 1){
    K = length(p)
    rv = rmultinom(n = n, size = 1, prob = p)
    z = rep(NA, length = n)
    for(i in 1:n){
        for(k in 1:K){
            if(rv[k, i] == 1){
                z[i] = k
            }
        }
    }
    return(z)
}

gibbs_sampling_z = function(x, K, mu){
    n = nrow(x)
    z = rep(NA, n)
    for(i in 1:n){
        p = rep(NA, length = K)
        for(k in 1:K){
            p[k] = dmvnorm(x = x[i, ], mean = mu[, k])
        }
        z[i] = rmulti(p = p/sum(p))
    }
    return(z)
}

calc_x_bar_k = function(x, z, K){
    n = nrow(x)
    d = ncol(x)
    n_k = rep(NA, legnth = K)
    x_bar_k = matrix(NA, nrow = d, ncol = K)
    for(k in 1:K){
        n_k[k] = length(z[z == k])
    }
    for(j in 1:d){
        for(k in 1:K){
            x_bar_k[j, k] = (1/n_k[k])*sum(x[z == k, j])
        }
    }
    args_list = list(n_k = n_k, x_bar_k = x_bar_k)
    return(args_list)
}

gibbs_sampling_mu = function(x, z, K, mu0, sigma0){
    d = ncol(x)
    list_x_bar_k = calc_x_bar_k(x = x, z = z, K = K)
    n_k = list_x_bar_k$n_k
    x_bar_k = list_x_bar_k$x_bar_k
    mu = matrix(NA, nrow = d, ncol = K)
    for(k in 1:K){
        mean_k = (n_k[k]/(n_k[k] + sigma0^2))*x_bar_k[, k] + (sigma0^2/(n_k[k] + sigma0^2)*mu0)
        Sigma_k = (n_k[k] + sigma0^2)^(-1)*diag(d)
        mu[, k] = rmvnorm(n = 1, mean = mean_k, sigma = Sigma_k)
    }
    return(mu)
}

z_map = function(z, K){
    n = nrow(z)
    z_count = matrix(NA, nrow = n, ncol = K)
    for(i in 1:n){
        for(k in 1:K){
            z_count[i, k] = length(z[i, z[i, ] == k])
        }
    }
    z_map = rep(NA, n)
    for(i in 1:n){
        for(k in 1:K){
            if(z_count[i, k] == max(z_count[i, ])){
                z_map[i] = k
            }
        }
    }
    return(z_map)
}

gaussian_mixture_model = function(
    x, K, iter_max = 1000, burn_in = 100,
    mu0 = rep(0, length = ncol(x)), sigma0 = 1, p0 = rep(1/K, length = K)
){
    n = nrow(x)
    d = ncol(x)

    mu = array(NA, dim = c(d, K, iter_max + 1))
    z = array(NA, dim = c(n, iter_max))

    mu[, , 1] = set_initial_mu(x = x, K = K)

    for(s in 1:iter_max){
        z[, s] = gibbs_sampling_z(x = x, K = K, mu = mu[, , s])
        mu[, , s + 1] = gibbs_sampling_mu(x = x, K = K, z = z[, s], mu0 = mu0, sigma0 = sigma0)
        if(s %% 100 == 0){
            cat("number of iteration is:", s, "\n")
            cat("class is:", z_map(z[, 1:s], K = K), "\n")
        }
    }

    args_list = list(
        z = z[, burn_in:iter_max],
        mu = mu[, , burn_in:iter_max],
        K = K,
        z_map  = z_map(z = z[, burn_in:iter_max], K = K)
    )
    return(args_list)
    message("the process has finished")
}
