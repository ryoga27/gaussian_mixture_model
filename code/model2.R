if(require(MCMCpack) == FALSE){
    install.packages("mvtnorm", dependencies = TRUE)
}
if(require(mvtnorm) == FALSE){
    install.packages("mvtnorm", dependencies = TRUE)
}
library(MCMCpack)
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

rmulti = function(prob, n = 1){
    K = length(prob)
    rv = rmultinom(n = n, size = 1, prob = prob)
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

gibbs_sampling_z = function(x, K, mu, p){
    n = nrow(x)
    z = rep(NA, n)
    for(i in 1:n){
        prob = rep(NA, length = K)
        for(k in 1:K){
            prob[k] = dmvnorm(x = x[i, ], mean = mu[, k])*p[k]
        }
        z[i] = rmulti(prob = prob/sum(prob))
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

gibbs_sampling_mu = function(x, z, K, tau, mu0, rho0){
    d = ncol(x)
    list_x_bar_k = calc_x_bar_k(x = x, z = z, K = K)
    n_k = list_x_bar_k$n_k
    x_bar_k = list_x_bar_k$x_bar_k
    mu = matrix(NA, nrow = d, ncol = K)
    for(k in 1:K){
        mean_k = (n_k[k]/(n_k[k] + rho0))*x_bar_k[, k] + (rho0/(n_k[k] + rho0)*mu0)
        Sigma_k = (tau*(n_k[k] + rho0))^(-1)*diag(d)
        mu[, k] = rmvnorm(n = 1, mean = mean_k, sigma = Sigma_k)
    }
    return(mu)
}

gibbs_sampling_tau = function(x, K, z, mu0, rho0, a0, b0){
    n = nrow(x)
    d = ncol(x)
    list_x_bar_k = calc_x_bar_k(x = x, z = z, K = K)
    n_k = list_x_bar_k$n_k
    x_bar_k = list_x_bar_k$x_bar_k

    a = a0 + (n*d)/2

    b = NA

    b1 = NA
    b1_ = rep(NA, length = n)
    b1__ = rep(NA, length = K)
    for(k in 1:K){
        for(i in 1:n){
            if(z[i] == k){
                b1_[i] = (1/2)*t(x[i, ] - x_bar_k[, k])%*%(x[i, ] - x_bar_k[, k])
            }
            if(z[i] != k){
                b1_[i] = 0
            }
        }
        b1__[k] = sum(b1_)
    }
    b1 = sum(b1__)

    b2 = NA
    b2_ = rep(NA, length = K)
    for(k in 1:K){
        b2_[k] = (n_k[k]*rho0)/(2*(n_k[k] + rho0))*t(x_bar_k[, k] - mu0)%*%(x_bar_k[, k] - mu0)
    }
    b2 = sum(b2_)

    b = b0 + b1 + b2
    tau = rgamma(n = 1, shape = a, scale = 1/b)
    return(tau)
}

gibbs_sampling_p = function(x, z, K, alpha){
    n = nrow(x)
    d = ncol(x)
    list_x_bar_k = calc_x_bar_k(x = x, z = z, K = K)
    n_k = list_x_bar_k$n_k
    x_bar_k = list_x_bar_k$x_bar_k
    a = alpha + n_k
    p = rdirichlet(n = 1, a = n_k)
    return(p)
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
    mu0 = rep(0, length = ncol(x)),
    tau0 = 1, p0 = rep(1/K, length = K), alpha0 = rep(1/K, length = K), a0 = 1, b0 = 1, rho0 = 1
){
    n = nrow(x)
    d = ncol(x)

    mu = array(NA, dim = c(d, K, iter_max + 1))
    z = array(NA, dim = c(n, iter_max))
    tau = rep(NA, length = iter_max + 1)
    p = array(NA, dim = c(K, iter_max + 1))

    mu[, , 1] = set_initial_mu(x = x, K = K)
    tau[1] = tau0
    p[, 1] = p0

    for(s in 1:iter_max){
        z[, s] = gibbs_sampling_z(x = x, K = K, mu = mu[, , s], p = p[, 1])
        mu[, , s + 1] = gibbs_sampling_mu(x = x, K = K, z = z[, s], mu0 = mu0, tau = tau[1], rho0 = rho0)
        tau[s + 1] = gibbs_sampling_tau(x = x, K = K, z = z[, s], mu0 = mu0, rho0 = rho0, a0 = a0, b0 = b0)
        p[, s + 1] = gibbs_sampling_p(x = x, K = K, z[, s], alpha = alpha0)
        if(s %% 100 == 0){
            cat("number of iteration is:", s, "\n")
            cat("class is:", z_map(z[, 1:s], K = K), "\n")
        }
    }

    args_list = list(
        z = z[, burn_in:iter_max],
        mu = mu[, , burn_in:iter_max + 1],
        tau = tau[burn_in:iter_max + 1],
        p = p[, burn_in:iter_max + 1],
        K = K,
        z_map  = z_map(z = z[, burn_in:iter_max], K = K)
    )
    return(args_list)
    message("the process has finished")
}
