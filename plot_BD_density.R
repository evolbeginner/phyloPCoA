#! /bin/env Rscriot

##########################################
# plotting a pdf under a BD process (Yang & Rannala 2006 MBE)


##########################################
g_function <- function(t, lambda, mu, rho, t1) {
    if (lambda == mu) {
        return((1 + rho * lambda * t1) / (t1 * (1 + rho * lambda * t)^2))
    }
    
    P_0_t <- rho * (lambda - mu) / (rho * lambda + (lambda * (1 - rho) - mu) * exp((mu - lambda) * t))
    p1_t <- (1/rho) * P_0_t^2 * exp((mu - lambda) * t)
    P_0_t1 <- rho * (lambda - mu) / (rho * lambda + (lambda * (1 - rho) - mu) * exp((mu - lambda) * t1))
    v_t1 <- 1 - (1/rho) * P_0_t1 * exp((mu - lambda) * t1)
    
    return((lambda * p1_t) / v_t1)
}

plot_g_function <- function(lambda, mu, rho, root_age, t_max = NULL) {
    t1 <- root_age
    t_max <- ifelse(is.null(t_max), t1, t_max)
    t_values <- seq(0, t_max, length.out = 100)
    
    g_values <- sapply(t_values, g_function, lambda = lambda, mu = mu, rho = rho, t1 = t1)
    
    plot(t_values, g_values, type = "l", col = "blue",
         xlab = "Time (t)", ylab = "g(t)", ylim = c(0, max(g_values)*1.5),
         main = ifelse(lambda == mu,
                       "Special Case (λ=μ)",
                       "General Case"))
    
    abline(v = t1, lty = 2, col = "gray")
    legend("topright", legend = paste0("λ=", lambda, ", μ=", mu, ", ρ=", rho, ", t1=", t1))
}


##########################################
plot_g_function(lambda = 5, mu = 5, rho = 0.0001, root_age = 7.5)
plot_g_function(lambda = 5, mu = 3, rho = 0.0001, root_age = 10)


