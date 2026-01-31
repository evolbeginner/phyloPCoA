#!/bin/env Rscript

##########################################
# Classify BD density functions by shape
##########################################

g <- function(t, lambda, mu, rho, t1) {
    if (lambda == mu) {
        return((1 + rho * lambda * t1) / (t1 * (1 + rho * lambda * t)^2))
    }
    
    P_0_t <- rho * (lambda - mu) / (rho * lambda + (lambda * (1 - rho) - mu) * exp((mu - lambda) * t))
    p1_t <- (1/rho) * P_0_t^2 * exp((mu - lambda) * t)
    P_0_t1 <- rho * (lambda - mu) / (rho * lambda + (lambda * (1 - rho) - mu) * exp((mu - lambda) * t1))
    v_t1 <- 1 - (1/rho) * P_0_t1 * exp((mu - lambda) * t1)
    
    return((lambda * p1_t) / v_t1)
}

##########################################
# Extract shape features from a density
##########################################

extract_features <- function(lambda, mu, rho, t1 = 10, n_points = 500) {
    t_values <- seq(0.001, t1, length.out = n_points)
    g_values <- sapply(t_values, g, lambda = lambda, mu = mu, rho = rho, t1 = t1)
    
    # Handle invalid cases
    if (any(!is.finite(g_values))) {
        return(list(valid = FALSE))
    }
    
    # Basic statistics
    g_max <- max(g_values)
    g_min <- min(g_values)
    g_mean <- mean(g_values)
    g_at_0 <- g_values[1]
    g_at_t1 <- g_values[n_points]
    
    # Find mode (peak location)
    mode_idx <- which.max(g_values)
    mode_t <- t_values[mode_idx]
    mode_value <- g_values[mode_idx]
    
    # Monotonicity
    diffs <- diff(g_values)
    is_increasing <- all(diffs >= -1e-10)
    is_decreasing <- all(diffs <= 1e-10)
    is_monotonic <- is_increasing || is_decreasing
    
    # Skewness (where is the mass concentrated?)
    # Left-skewed: mass toward t1, Right-skewed: mass toward 0
    cumsum_g <- cumsum(g_values) / sum(g_values)
    median_idx <- which(cumsum_g >= 0.5)[1]
    median_t <- t_values[median_idx]
    relative_median <- median_t / t1  # 0 to 1, lower = mass near 0
    
    # Curvature at different regions
    n_third <- floor(n_points / 3)
    slope_early <- mean(diffs[1:n_third])
    slope_mid <- mean(diffs[(n_third+1):(2*n_third)])
    slope_late <- mean(diffs[(2*n_third+1):(n_points-1)])
    
    # Concentration (how peaked is the distribution?)
    # Using coefficient of variation
    cv <- sd(g_values) / g_mean
    
    # Ratio of endpoints
    endpoint_ratio <- g_at_t1 / g_at_0
    
    return(list(
        valid = TRUE,
        lambda = lambda,
        mu = mu,
        rho = rho,
        turnover = mu / lambda,  # extinction/speciation ratio
        net_div = lambda - mu,   # net diversification
        g_max = g_max,
        g_min = g_min,
        g_mean = g_mean,
        g_at_0 = g_at_0,
        g_at_t1 = g_at_t1,
        mode_t = mode_t,
        relative_mode = mode_t / t1,
        is_monotonic = is_monotonic,
        is_increasing = is_increasing,
        is_decreasing = is_decreasing,
        relative_median = relative_median,
        slope_early = slope_early,
        slope_mid = slope_mid,
        slope_late = slope_late,
        cv = cv,
        endpoint_ratio = endpoint_ratio
    ))
}

##########################################
# Classify based on extracted features
##########################################

classify_shape <- function(features) {
    if (!features$valid) return("Invalid")
    
    # Classification logic based on shape characteristics
    
    # 1. Monotonically decreasing (exponential-like decay from t=0)
    if (features$is_decreasing && features$relative_mode < 0.1) {
        if (features$cv > 0.5) {
            return("Steep Decay")
        } else {
            return("Gradual Decay")
        }
    }
    
    # 2. Monotonically increasing (mass concentrated near t1)
    if (features$is_increasing && features$relative_mode > 0.9) {
        if (features$cv > 0.5) {
            return("Steep Rise")
        } else {
            return("Gradual Rise")
        }
    }
    
    # 3. Unimodal with interior peak
    if (!features$is_monotonic && features$relative_mode > 0.1 && features$relative_mode < 0.9) {
        if (features$relative_mode < 0.4) {
            return("Early Peak")
        } else if (features$relative_mode > 0.6) {
            return("Late Peak")
        } else {
            return("Central Peak")
        }
    }
    
    # 4. Nearly uniform
    if (features$cv < 0.15) {
        return("Near Uniform")
    }
    
    # 5. Boundary cases
    if (features$relative_mode <= 0.1) {
        return("Left-Boundary Mode")
    }
    if (features$relative_mode >= 0.9) {
        return("Right-Boundary Mode")
    }
    
    return("Other")
}

##########################################
# Main analysis
##########################################

# Parameter ranges
lambda_seq <- seq(0.05, 0.5, by = 0.05)
mu_seq <- seq(0, 0.45, by = 0.05)
rho_seq <- seq(0.2, 0.8, by = 0.2)
t1 <- 10  # fixed root age

# Generate all valid parameter combinations
param_grid <- expand.grid(lambda = lambda_seq, mu = mu_seq, rho = rho_seq)
param_grid <- param_grid[param_grid$mu <= param_grid$lambda, ]

cat("Total parameter combinations:", nrow(param_grid), "\n\n")

# Extract features and classify each combination
results <- list()
for (i in 1:nrow(param_grid)) {
    features <- extract_features(
        lambda = param_grid$lambda[i],
        mu = param_grid$mu[i],
        rho = param_grid$rho[i],
        t1 = t1
    )
    features$class <- classify_shape(features)
    results[[i]] <- features
}

# Convert to data frame
results_df <- do.call(rbind, lapply(results, function(x) {
    if (!x$valid) {
        return(data.frame(
            lambda = x$lambda, mu = x$mu, rho = x$rho,
            turnover = NA, net_div = NA, class = "Invalid",
            relative_mode = NA, relative_median = NA, cv = NA
        ))
    }
    data.frame(
        lambda = x$lambda, mu = x$mu, rho = x$rho,
        turnover = x$turnover, net_div = x$net_div, class = x$class,
        relative_mode = x$relative_mode, relative_median = x$relative_median, cv = x$cv
    )
}))

##########################################
# Summary of classification
##########################################

cat("========================================\n")
cat("CLASSIFICATION SUMMARY (root_age = 10)\n")
cat("========================================\n\n")

class_table <- table(results_df$class)
print(class_table)

cat("\n\nPercentage breakdown:\n")
print(round(prop.table(class_table) * 100, 1))

##########################################
# Detailed breakdown by class
##########################################

cat("\n\n========================================\n")
cat("DETAILED CLASS DESCRIPTIONS\n")
cat("========================================\n")

classes <- unique(results_df$class)
for (cls in classes) {
    subset_df <- results_df[results_df$class == cls, ]
    cat("\n---", cls, "---\n")
    cat("Count:", nrow(subset_df), "\n")
    cat("Lambda range:", min(subset_df$lambda), "-", max(subset_df$lambda), "\n")
    cat("Mu range:", min(subset_df$mu), "-", max(subset_df$mu), "\n")
    cat("Rho range:", min(subset_df$rho), "-", max(subset_df$rho), "\n")
    cat("Turnover (mu/lambda) range:", 
        round(min(subset_df$turnover, na.rm=T), 2), "-", 
        round(max(subset_df$turnover, na.rm=T), 2), "\n")
    cat("Net diversification range:", 
        round(min(subset_df$net_div, na.rm=T), 3), "-", 
        round(max(subset_df$net_div, na.rm=T), 3), "\n")
}

##########################################
# Key parameter relationships
##########################################

cat("\n\n========================================\n")
cat("KEY PARAMETER RELATIONSHIPS\n")
cat("========================================\n")

cat("\nClass by Turnover Rate (mu/lambda):\n")
results_df$turnover_bin <- cut(results_df$turnover, 
                                breaks = c(-0.01, 0.25, 0.5, 0.75, 1.01),
                                labels = c("Low (0-0.25)", "Med-Low (0.25-0.5)", 
                                          "Med-High (0.5-0.75)", "High (0.75-1)"))
print(table(results_df$turnover_bin, results_df$class))

cat("\nClass by Sampling Fraction (rho):\n")
print(table(results_df$rho, results_df$class))

##########################################
# Generate visualization PDF
##########################################

pdf("BD_density_classification.pdf", width = 16, height = 12)

# Plot 1: All densities colored by class
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

# Define colors for each class
class_colors <- c(
    "Steep Decay" = "#E41A1C",
    "Gradual Decay" = "#FF7F00",
    "Steep Rise" = "#377EB8",
    "Gradual Rise" = "#4DAF4A",
    "Early Peak" = "#984EA3",
    "Central Peak" = "#FFFF33",
    "Late Peak" = "#A65628",
    "Near Uniform" = "#F781BF",
    "Left-Boundary Mode" = "#999999",
    "Right-Boundary Mode" = "#66C2A5",
    "Other" = "#000000",
    "Invalid" = "#FFFFFF"
)

t_values <- seq(0.001, t1, length.out = 200)

# Initialize plot
plot(NULL, xlim = c(0, t1), ylim = c(0, 2),
     xlab = "Time (t)", ylab = "g(t)",
     main = "BD Density Functions Colored by Shape Class (root_age = 10)",
     cex.main = 1.5, cex.lab = 1.2)

# Plot each density
for (i in 1:nrow(param_grid)) {
    g_values <- sapply(t_values, g, 
                       lambda = param_grid$lambda[i],
                       mu = param_grid$mu[i],
                       rho = param_grid$rho[i],
                       t1 = t1)
    if (all(is.finite(g_values))) {
        cls <- results_df$class[i]
        lines(t_values, g_values, col = adjustcolor(class_colors[cls], alpha = 0.4), lwd = 0.8)
    }
}

# Add legend
legend("topright", legend = names(class_colors)[names(class_colors) %in% unique(results_df$class)],
       col = class_colors[names(class_colors) %in% unique(results_df$class)],
       lwd = 2, cex = 0.8, bg = "white")

# Plot 2: Representative examples from each class
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))

for (cls in unique(results_df$class)) {
    subset_df <- results_df[results_df$class == cls, ]
    if (nrow(subset_df) == 0) next
    
    # Pick middle example
    idx <- ceiling(nrow(subset_df) / 2)
    example <- subset_df[idx, ]
    
    g_values <- sapply(t_values, g, 
                       lambda = example$lambda,
                       mu = example$mu,
                       rho = example$rho,
                       t1 = t1)
    
    plot(t_values, g_values, type = "l", col = class_colors[cls], lwd = 2,
         xlab = "Time", ylab = "g(t)", ylim=c(0, max(g_values)*1.1),
         main = paste0(cls, "\n(λ=", example$lambda, ", μ=", example$mu, ", ρ=", example$rho, ")"),
         cex.main = 0.9)
    abline(v = t1, lty = 2, col = "gray")
}

# Plot 3: Heatmap of classes by lambda and mu (faceted by rho)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

for (r in rho_seq) {
    subset_df <- results_df[results_df$rho == r, ]
    
    # Create matrix for plotting
    class_matrix <- matrix(NA, nrow = length(lambda_seq), ncol = length(mu_seq))
    rownames(class_matrix) <- lambda_seq
    colnames(class_matrix) <- mu_seq
    
    class_numeric <- as.numeric(factor(subset_df$class))
    
    for (j in 1:nrow(subset_df)) {
        li <- which(lambda_seq == subset_df$lambda[j])
        mi <- which(mu_seq == subset_df$mu[j])
        class_matrix[li, mi] <- class_numeric[j]
    }
    
    image(lambda_seq, mu_seq, class_matrix,
          xlab = "Lambda", ylab = "Mu",
          main = paste("Shape Classes (rho =", r, ")"),
          col = rainbow(length(unique(results_df$class))))
    abline(0, 1, lty = 2, col = "black", lwd = 2)
    text(0.35, 0.15, "mu <= lambda\nregion", cex = 0.8)
}

dev.off()

cat("\n\nVisualization saved to 'BD_density_classification.pdf'\n")

##########################################
# Export classification results
##########################################

write.csv(results_df, "BD_density_classification.csv", row.names = FALSE)
cat("Classification data saved to 'BD_density_classification.csv'\n")
