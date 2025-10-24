min_sigma <- 0; max_sigma <- 1

#set.seed(123)

# ---- Packages ----
pkgs <- c("ape", "TreeSim", "mvMORPH", "ggplot2", "vegan", "patchwork")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Parameters ----
tnum       <- 100   # number of extant taxa (tips)
lambda     <- 5   # birth rate
mu         <- 4.5   # death rate
rho        <- 0.1  # sampling fraction
age        <- 10   # crown age

n_microbes <- 20
char_names <- c(paste0("microbe_", seq_len(n_microbes)), "body_size")
k          <- length(char_names)

# ---- Helpers ----
clr <- function(P) {
    L <- log(P)
    sweep(L, 1, rowMeans(L), FUN = "-")
}

make_group_circles <- function(df, group_col = "group", xcol = "x", ycol = "y",
                               n_points = 360, expand = 1.05) {
    groups <- unique(df[[group_col]])
    circles <- lapply(groups, function(g) {
        sub <- df[df[[group_col]] == g, , drop = FALSE]
        cx <- mean(sub[[xcol]]); cy <- mean(sub[[ycol]])
        r <- max(sqrt((sub[[xcol]] - cx)^2 + (sub[[ycol]] - cy)^2))
        theta <- seq(0, 2 * pi, length.out = n_points)
        data.frame(group = g,
                   x = cx + expand * r * cos(theta),
                   y = cy + expand * r * sin(theta))
    })
    do.call(rbind, circles)
}

safe_chol <- function(M, jitter = 1e-10, max_tries = 8) {
    for (i in 0:max_tries) {
        diag_add <- if (i == 0) 0 else jitter * 2^(i - 1)
        Mj <- M
        if (diag_add > 0) diag(Mj) <- diag(Mj) + diag_add
        out <- try(chol(Mj), silent = TRUE)
        if (!inherits(out, "try-error")) return(out)
    }
    stop("Cholesky factorization failed; matrix may be singular.")
}

# Cholesky-based phylogenetic transform with GLS root values
phylo_transform_L <- function(Y, C, scale_by_first_diag = TRUE) {
    stopifnot(all(rownames(Y) %in% rownames(C)))
    C <- C[rownames(Y), rownames(Y), drop = FALSE]
    
    U <- safe_chol(C)   # upper triangular, t(U) %*% U = C
    L <- solve(t(U))    # C^{-1/2}
    
    if (isTRUE(scale_by_first_diag)) {
        L <- L * sqrt(diag(C)[1])
    }
    
    C_inv <- chol2inv(U)
    ones  <- matrix(1, nrow = nrow(C), ncol = 1)
    
    root_values <- apply(Y, 2, function(col) {
        numerator   <- as.numeric(t(col) %*% C_inv %*% ones)
        denominator <- as.numeric(t(ones) %*% C_inv %*% ones)
        numerator / denominator
    })
    
    X <- L %*% t(t(Y) - root_values)
    X <- t(t(X) + root_values)
    
    eig_decomp <- eigen(C)
    L <- eig_decomp$vectors %*% diag(1/sqrt(eig_decomp$values)) %*% t(eig_decomp$vectors)
    L <- L * sqrt(diag(C)[1])
    X <- L %*% t((t(Y) - root_values))
    X <- t(t(X) + root_values)
    #X <- L %*% Y
    
    rownames(X) <- rownames(Y); colnames(X) <- colnames(Y)
    X
}

# PCoA via ape::pcoa with percent variance for axes 1 and 2
run_pcoa2 <- function(d) {
    pc <- ape::pcoa(d)
    eig <- pc$values$Eigenvalues
    pos_sum <- sum(eig[eig > 0])
    ax1_var <- if (length(eig) >= 1 && eig[1] > 0 && pos_sum > 0) eig[1] / pos_sum else 0
    ax2_var <- if (length(eig) >= 2 && eig[2] > 0 && pos_sum > 0) eig[2] / pos_sum else 0
    list(scores = pc$vectors[, 1:2, drop = FALSE],
         var = c(ax1 = ax1_var, ax2 = ax2_var))
}

# ---- 1) Simulate a time-tree phylogeny ----
sim_bd_tree_retry <- function(tnum, lambda, mu, rho, age, max_tries = 300) {
    for (i in seq_len(max_tries)) {
        trees <- try(TreeSim::sim.bd.taxa.age(tnum, 1, lambda, mu, frac = rho, age = age),
                     silent = TRUE)
        if (!inherits(trees, "try-error") && length(trees) == 1 && inherits(trees[[1]], "phylo")) {
            tr <- trees[[1]]
            if (ape::Ntip(tr) == tnum) return(tr)
        }
    }
    stop("Failed to simulate a time-tree with the requested parameters after retries.")
}

tree <- sim_bd_tree_retry(tnum, lambda, mu, rho, age)
tree$tip.label <- paste0("host_", seq_len(ape::Ntip(tree)))
tip_order <- tree$tip.label

# Host BM covariance (shared path lengths)
C <- ape::vcv(tree)
C <- C[tip_order, tip_order, drop = FALSE]

# ---- 2) Among-character variance–covariance matrix Sigma (k x k) ----
loading_microbe <- 0.05
loading_body    <- 13.95
B  <- c(rep(loading_microbe, n_microbes), loading_body)
U2 <- c(rep(1 - loading_microbe^2, n_microbes), 1 - loading_body^2)
Sigma <- outer(B, B) + diag(U2)

Sigma <- crossprod(matrix(runif((n_microbes+1)*(n_microbes+1), min_sigma, max_sigma),n_microbes+1))

#Sigma[(n_microbes+1),1:(n_microbes+1)] = runif(n_microbes+1,0.8,1); Sigma[1:(n_microbes+1),(n_microbes+1)] = Sigma[(n_microbes+1),1:(n_microbes+1)] #SW
#Sigma[(n_microbes+1),1:(n_microbes+1)] = runif(n_microbes+1,0.9,1); Sigma[1:(n_microbes+1),(n_microbes+1)] = Sigma[(n_microbes+1),1:(n_microbes+1)]; Sigma[n_microbes+1,n_microbes+1]=1 #SW
#Sigma <- ifelse(row(Sigma) == col(Sigma), Sigma, runif(1,0,1))
#print(Sigma)

dimnames(Sigma) <- list(char_names, char_names)
#invisible(safe_chol(Sigma))
theta <- setNames(rep(0, k), char_names)
theta <- setNames(runif(k,0,1), char_names) #SW

# ---- 3) Joint simulation with mvSIM (BM1) ----
param <- list(sigma = Sigma, theta = theta)
sim_raw <- mvMORPH::mvSIM(tree, nsim = 1, model = "BM1", param = param)
sim_mat <- as.matrix(sim_raw)

# ---- 4) Robust alignment by names ----
if (is.null(rownames(sim_mat))) stop("mvSIM returned no rownames; cannot align tips.")
missing_tips <- setdiff(tip_order, rownames(sim_mat))
if (length(missing_tips)) stop(sprintf("Missing tips in mvSIM output: %s", paste(missing_tips, collapse = ", ")))
row_idx <- match(tip_order, rownames(sim_mat))

if (is.null(colnames(sim_mat))) {
    if (ncol(sim_mat) == length(char_names)) {
        colnames(sim_mat) <- char_names
    } else {
        stop("mvSIM returned no colnames and dimensions do not match expected traits.")
    }
}
missing_traits <- setdiff(char_names, colnames(sim_mat))
if (length(missing_traits)) stop(sprintf("Missing traits in mvSIM output: %s", paste(missing_traits, collapse = ", ")))
col_idx <- match(char_names, colnames(sim_mat))

sim_mat <- sim_mat[row_idx, col_idx, drop = FALSE]
stopifnot(identical(rownames(sim_mat), tip_order))
stopifnot(identical(colnames(sim_mat), char_names))

# ---- 5) Split traits: microbial latents -> compositions; body size (log) ----
latent_log_abund <- sim_mat[, paste0("microbe_", seq_len(n_microbes)), drop = FALSE]
body_size_log    <- sim_mat[, "body_size"]     # simulated as logody size)
body_size_raw    <- exp(body_size_log)
#print(body_size_log)
#print(sim_mat)

exp_latent <- exp(latent_log_abund)
microbe_prop <- exp_latent / rowSums(exp_latent)
stopifnot(all(abs(rowSums(microbe_prop) - 1) < 1e-10))

# ---- 6) CLR and groupings ----
clr_mat <- clr(microbe_prop)  # hosts x microbes
clr_mat <- latent_log_abund #SW
rownames(clr_mat) <- tip_order

# RAW groups (used for all panels per your update)
QUANTILE <- 0.5
group_raw <- factor(
    ifelse(body_size_raw > quantile(body_size_raw, QUANTILE), "Large (raw)", "Small (raw)"),
    levels = c("Small (raw)", "Large (raw)")
)


# ---- 7) PCoAs ----
# 7.1 Aitchison (non-phylogenetic)
dist_aitchison <- dist(latent_log_abund)
pc1 <- run_pcoa2(dist_aitchison)
scores1 <- data.frame(
    host  = rownames(clr_mat),
    Axis1 = pc1$scores[, 1],
    Axis2 = pc1$scores[, 2],
    group = group_raw
)
circles1 <- make_group_circles(
    data.frame(group = scores1$group, x = scores1$Axis1, y = scores1$Axis2),
    group_col = "group", xcol = "x", ycol = "y", expand = 1.05
)

# 7.2 Aitchison (phylogeny-whitened) — color and circle by RAW values
Z_clr <- phylo_transform_L(clr_mat, C, scale_by_first_diag = TRUE)
dist_phylo_aitch <- dist(Z_clr)
pc2 <- run_pcoa2(dist_phylo_aitch)
scores2 <- data.frame(
    host  = rownames(Z_clr),
    Axis1 = pc2$scores[, 1],
    Axis2 = pc2$scores[, 2],
    group = group_raw
)
circles2 <- make_group_circles(
    data.frame(group = scores2$group, x = scores2$Axis1, y = scores2$Axis2),
    group_col = "group", xcol = "x", ycol = "y", expand = 1.05
)

# 7.3 Bray–Curtis (relative abundance)
dist_bray <- vegan::vegdist(microbe_prop, method = "bray")
pc3 <- run_pcoa2(dist_bray)
scores3 <- data.frame(
    host  = rownames(microbe_prop),
    Axis1 = pc3$scores[, 1],
    Axis2 = pc3$scores[, 2],
    group = group_raw
)
circles3 <- make_group_circles(
    data.frame(group = scores3$group, x = scores3$Axis1, y = scores3$Axis2),
    group_col = "group", xcol = "x", ycol = "y", expand = 1.05
)

# ---- 8) Plots (separate; each uses its own limits) ----
col_raw <- c("Small (raw)" = "#377eb8", "Large (raw)" = "#e41a1c")

p1 <- ggplot(scores1, aes(x = Axis1, y = Axis2, color = group)) +
    geom_point(size = 2.6, alpha = 0.9) +
    geom_path(data = circles1, aes(x = x, y = y, color = group),
              size = 1, show.legend = FALSE, inherit.aes = FALSE) +
    coord_equal() +
    scale_color_manual(values = col_raw) +
    labs(title = "Aitchison (non-phylogenetic)",
         x = sprintf("PCoA1 (%.1f%%)", 100*pc1$var["ax1"]),
         y = sprintf("PCoA2 (%.1f%%)", 100*pc1$var["ax2"]),
         color = "Body size (raw)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold"))

p2 <- ggplot(scores2, aes(x = Axis1, y = Axis2, color = group)) +
    geom_point(size = 2.6, alpha = 0.9) +
    geom_path(data = circles2, aes(x = x, y = y, color = group),
              size = 1, show.legend = FALSE, inherit.aes = FALSE) +
    coord_equal() +
    scale_color_manual(values = col_raw) +
    labs(title = "Aitchison (phylogeny-whitened)",
         x = sprintf("PCoA1 (%.1f%%)", 100*pc2$var["ax1"]),
         y = sprintf("PCoA2 (%.1f%%)", 100*pc2$var["ax2"]),
         color = "Body size (raw)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold"))

# Combine the two Aitchison panels side-by-side
p12 <- p1 | p2

# Bray–Curtis as its own plot with its own limits
p3 <- ggplot(scores3, aes(x = Axis1, y = Axis2, color = group)) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_path(data = circles3, aes(x = x, y = y, color = group),
              size = 1, show.legend = FALSE, inherit.aes = FALSE) +
    coord_equal() +
    scale_color_manual(values = col_raw) +
    labs(title = "Bray–Curtis (relative abundance)",
         x = sprintf("PCoA1 (%.1f%%)", 100*pc3$var["ax1"]),
         y = sprintf("PCoA2 (%.1f%%)", 100*pc3$var["ax2"]),
         color = "Body size (raw)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold"))

# ---- 9) Draw ----
print(p12)
p3; plot(tree); cov2cor(Sigma)


