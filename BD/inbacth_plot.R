#!/bin/env Rscript

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
    
    # Handle potential NA/Inf values
    valid_idx <- is.finite(g_values)
    if (sum(valid_idx) < 2) {
        plot(1, type = "n", xlab = "", ylab = "",
             main = paste0("λ=", lambda, ",μ=", mu, ",ρ=", rho),
             cex.main = 0.5)
        text(1, 1, "Invalid", col = "red", cex = 0.6)
        return()
    }
    
    y_max <- max(g_values[valid_idx], na.rm = TRUE)
    
    plot(t_values, g_values, type = "l", col = "blue", lwd = 0.8,
         xlab = "", ylab = "", 
         ylim = c(0, min(y_max * 1.2, y_max + 0.5)),
         main = paste0("λ=", lambda, ",μ=", mu, ",ρ=", rho),
         cex.main = 0.4, cex.axis = 0.4,
         mgp = c(1, 0.2, 0), tcl = -0.2)
    
    abline(v = t1, lty = 2, col = "gray", lwd = 0.5)
}

##########################################
# Define parameter ranges based on the image
##########################################

lambda_seq <- seq(0.05, 0.5, by = 0.05)       # λ: 0.05 to 0.5, step 0.05 (10 values)
mu_seq <- seq(0, 0.45, by = 0.05)             # μ: 0 to 0.45, step 0.05 (10 values)
age_seq <- c(1, 5, 10, 50, 100)               # ages: 1, 5, 10, 50, 100 (5 values)
f_seq <- seq(0.2, 0.8, by = 0.2)              # f (rho): 0.2 to 0.8, step 0.2 (4 values)

cat("Parameter counts:\n")
cat("  lambda:", length(lambda_seq), "values\n")
cat("  mu:", length(mu_seq), "values\n")
cat("  rho:", length(f_seq), "values\n")
cat("  ages:", length(age_seq), "values\n")

##########################################
# Generate PDF with one page per root age
# Each page contains all combinations of lambda, mu, rho
##########################################

# Calculate grid dimensions for each page
# For each age: lambda x mu x rho combinations where mu <= lambda
# With 10 lambda, 10 mu, 4 rho values, and mu <= lambda constraint:
# Valid (lambda, mu) pairs: 10 + 9 + 8 + ... + 1 = 55 pairs
# Total per age: 55 * 4 = 220 plots (not 315)

# If you want exactly 315, we may need to adjust parameters
# Let's first see what we get with current parameters

# Create parameter grid for one age (to count)
test_grid <- expand.grid(lambda = lambda_seq, mu = mu_seq, rho = f_seq)
test_grid <- test_grid[test_grid$mu <= test_grid$lambda, ]
cat("\nPlots per age (with mu <= lambda):", nrow(test_grid), "\n")

# Adjust grid layout based on actual count
plots_per_age <- nrow(test_grid)

# Calculate optimal grid dimensions (approximately square)
ncol_plots <- ceiling(sqrt(plots_per_age * 1.5))  # wider than tall
nrow_plots <- ceiling(plots_per_age / ncol_plots)

cat("Grid layout:", nrow_plots, "rows x", ncol_plots, "cols =", nrow_plots * ncol_plots, "slots\n")

# Open PDF device (large page for many plots)
pdf("BD_density_plots_by_age.pdf", width = 24, height = 20)

for (age in age_seq) {
    cat("\nProcessing root age =", age, "...\n")
    
    # Generate parameter grid for this age
    param_grid <- expand.grid(
        lambda = lambda_seq,
        mu = mu_seq,
        rho = f_seq
    )
    
    # Filter: keep only valid combinations where mu <= lambda
    param_grid <- param_grid[param_grid$mu <= param_grid$lambda, ]
    
    # Sort by rho, then lambda, then mu for organized display
    param_grid <- param_grid[order(param_grid$rho, param_grid$lambda, param_grid$mu), ]
    
    # Set up the layout for this page
    par(mfrow = c(nrow_plots, ncol_plots), 
        mar = c(1.5, 1.5, 1.5, 0.5),  # margins: bottom, left, top, right
        oma = c(2, 2, 4, 1))          # outer margins
    
    # Plot all combinations for this age
    for (i in 1:nrow(param_grid)) {
        plot_g_function(
            lambda = param_grid$lambda[i],
            mu = param_grid$mu[i],
            rho = param_grid$rho[i],
            root_age = age
        )
    }
    
    # Fill remaining slots with empty plots
    remaining <- (nrow_plots * ncol_plots) - nrow(param_grid)
    if (remaining > 0) {
        for (j in 1:remaining) {
            plot.new()
        }
    }
    
    # Add page title with root age
    mtext(paste("BD Process Density Plots - Root Age =", age), 
          outer = TRUE, cex = 1.8, font = 2, line = 1)
    mtext("x-axis: Time (t), y-axis: g(t)", 
          outer = TRUE, cex = 1.0, line = -0.5)
    
    cat("  Completed", nrow(param_grid), "plots\n")
}

dev.off()

cat("\n========================================\n")
cat("PDF saved as 'BD_density_plots_by_age.pdf'\n")
cat("Total pages:", length(age_seq), "\n")
cat("Plots per page:", plots_per_age, "\n")
cat("========================================\n")


