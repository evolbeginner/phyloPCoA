#! /bin/env Rscript

# all rights reserved
# Sishuo's idea on accounting phylogeny for proportion data

###############################################################
#set.seed(2)
suppressWarnings(
suppressPackageStartupMessages({
    library(getopt)
    library(dplyr)
    library(tidyr)

    library(TreeSim)
    library(mvMORPH)
    library(phytools)
    library(phangorn)

    #library(adephylo)
    library(castor)
    library(vegan)
    library(labdsv) #pco
    library(compositions) #clr

    library(Matrix)
    library(MASS)
    library(cluster) #silhouette()
})
)


##################################
calculate_pcoa <- function(X, method = "bray", scaled=FALSE, grp_list) {
    if(scaled == TRUE){X <- scale(X)}
    distance <- vegdist(X, method = method)
    pcoa_res <- cmdscale(distance, eig = T)
    return(pcoa_res)
}


check_clustering <- function(pcoa_res, pcoa_name, grp_list, outfile){
    lda_input <- as.data.frame(pcoa_res$points)
    group_map <- unlist(lapply(names(grp_list), function(g) {
        setNames(rep(g, length(grp_list[[g]]$labels)), grp_list[[g]]$labels)
    }))
    sample_names <- names(group_map)
    grp_factor <- factor(group_map[sample_names])

    # do LDA
    lda_model <- lda(grp_factor ~ ., data = lda_input)
    lda_pred <- predict(lda_model)
    conf_mat <- table(True = grp_factor, Pred = lda_pred$class)
    acc <- sum(diag(conf_mat))/sum(conf_mat)

    dat <- cbind(as.data.frame(pcoa_res$points), group = grp_factor)
    k <- 2 #ncol(dat)-1
    fit <- manova(as.matrix(dat[, 1:k]) ~ group, data = dat)

    # anova
    #print(summary(aov(V1 ~ group, data = cbind(as.data.frame(pcoa_res$points[, 1:2]), group = grp_factor))))

    # fisher's discriminant ratio
    coords <- as.data.frame(pcoa_res$points[, 1:2])  # we only need first two PCs
    fdr_value <- fisher_ratio(coords, grp_factor)

    # BDI Davies-Bouldin index
    db_index <- davies_bouldin(coords, grp_factor)

    write(paste(pcoa_name, acc, fdr_value, db_index), file=outfile, append=TRUE, sep="\t")
    #cat(pcoa_name, acc, fdr_value, db_index, "\n", file = outfile, append = TRUE, sep = "\t")
}


fisher_ratio <- function(coords, group) {
  # compute global Fisher's discriminant ratio
  grand_center <- colMeans(coords)
  
  # Between-group scatter
  SB <- sum(sapply(levels(group), function(g) {
    n <- sum(group == g)
    diff <- colMeans(coords[group == g, ]) - grand_center
    n * sum(diff^2)
  }))
  
  # Within-group scatter
  SW <- sum(sapply(levels(group), function(g) {
    Xg <- coords[group == g, , drop = FALSE]
    sum(rowSums((Xg - colMeans(Xg))^2))
  }))
  
  ratio <- SB / SW
  return(ratio)
}


generate_metadata <- function(metadata_file){
    metadata <- read.table(metadata_file, header = TRUE, comment.char = "")
    # Randomly select one sample_id for each species
    selected_samples <- metadata %>%
        group_by(species) %>%
        slice(1) #always the 1st, reproducible
    return(selected_samples)
}


davies_bouldin <- function(coords, group) {
  groups <- levels(group)
  centroids <- sapply(groups, function(g) colMeans(coords[group == g, ]))
  centroids <- t(centroids)
  
  # Compute within-cluster scatter (average distance to centroid)
  S <- sapply(groups, function(g) {
    mean(dist(rbind(centroids[g,], coords[group == g, ]))[1:nrow(coords[group == g, ])])
  })
  
  # Compute inter-centroid distances
  M <- as.matrix(dist(centroids))
  
  # Compute pairwise Rij and then DB
  R <- matrix(0, length(groups), length(groups))
  for (i in seq_along(groups)) {
    for (j in seq_along(groups)) {
      if (i != j) {
        R[i, j] <- (S[i] + S[j]) / M[i, j]
      } else {
        R[i, j] <- NA
      }
    }
  }
  D_i <- apply(R, 1, max, na.rm = TRUE)
  DB <- mean(D_i)
  return(DB)
}


generate_metadata2 <- function(metadata_file, abundance){
    metadata <- read.table(metadata_file, header = TRUE, comment.char = "")
    abundance$taxon <- rownames(abundance)
    rownames(abundance) <- NULL
    abundance_long <- abundance %>%
        pivot_longer(cols = -taxon, names_to = "sample_id", values_to = "abundance")
    merged_data <- inner_join(abundance_long, metadata, by = "sample_id")
    merged_data <- merged_data %>%
        group_by(sample_id) %>%
        mutate(abundance = abundance / sum(abundance)) %>%
        ungroup()

    average_abundance <- merged_data %>%
      group_by(species, taxon) %>%
      summarize(average_abundance = mean(abundance), .groups = 'drop')
      #summarize(average_abundance = first(abundance), .groups = 'drop')
      #summarize(average_abundance = exp(mean(log(abundance+1e-5), na.rm = TRUE)), .groups = 'drop')

    # Reshape the average_abundance data frame to wide format
    average_abundance_wide <- average_abundance %>%
      pivot_wider(names_from = species, values_from = average_abundance)

    average_abundance_wide = as.data.frame(average_abundance_wide)
    rownames(average_abundance_wide) <- average_abundance_wide[,1]
    average_abundance_wide <- average_abundance_wide[,-1]

    return(average_abundance_wide)
}


simulate_covariance3 <- function(n=5, rate = 1, exponent = 4) {
  # Step 1: Generate random orthogonal basis (once)
  Q <- qr.Q(qr(matrix(rnorm(n^2), n)))
  
  # Step 2: Baseline eigenvalues
  base_eigs <- rexp(n, rate = 1)
  base_eigs <- base_eigs / max(base_eigs)  # scale m1 = 1
  
  # Step 3: For each exponent, build scaled covariance matrix
    eigs_a <- base_eigs^exponent       # apply exponent
    eigs_a <- eigs_a / mean(eigs_a)  # normalize total variance if desired
    Rho <- cov2cor(Q %*% diag(eigs_a) %*% t(Q))

    sigma_sq <- sqrt(rexp(n, rate = rate))
    #print(sigma_sq)
    Sigma <- diag(sigma_sq) %*% Rho %*% diag(sigma_sq)
  
  return(list(Rho=Rho, Sigma=Sigma))
}


# XQ
# bnum: bac taxonomic unit # in the microbiota
# tnum: tip (host) species
do_sim <- function(tnum, bnum, lambda, mu, rho, age, exponent){
    # BD model
    trees <- sim.bd.taxa.age(tnum, 1, lambda, mu, rho, age)
    tree <- trees[[1]]
    write.tree(tree, file=file.path(outdir, 'sim.tree'))

    # diff root_values
    mean_a=0; sd_a=1; nsim <- bnum + 1 #the last is the phenotype trait
    #abundance <- matrix(0, nrow = tnum, ncol = nsim); 

    # generate matrix R (randomly)
    #scale_R <- 1/distRoot(tree)[1]
    scale_R <- 1/get_all_distances_to_root(tree)[1]
    # XQ, R controls sigma in mvSIM() (covariance matrix)
    #R <- crossprod(matrix(runif(bnum*bnum),bnum)) * scale_R
    #sim_cov_res <- simulate_covariance2(bnum, rate=1/scale_R, weak_range=c(-0.01,0.01), strong_range=c(0.8,1), split_quantile=0.2)

    # Uyeda2015
    sim_cov_res <- simulate_covariance3(nsim, rate=1/scale_R, exponent=exponent)
    Sigma <- sim_cov_res$Sigma
    Rho <- sim_cov_res$Rho

    rs <- rnorm(ncol(Sigma), mean_a, sd_a) # diff root_values, XQ
    abundance <- exp(mvSIM(tree, model="BM1", nsim=1, param=list(sigma=Sigma, theta=rs)))
    log_abundance <- abundance

    above_names <- get_binary(abundance, rs) #last column above root value

    # remove the last column (binary trait, not microbial abundance)
    num_cols <- ncol(abundance)
    abundance <- abundance[, 1:(num_cols - 1)]

    a_values <- rnorm(nsim, mean = mean_a, sd = sd_a)
    C <- vcv(tree)
    prop <- t(apply(abundance, 1, function(x) x/sum(x)))

    return(list(abundance=abundance, prop=prop, C=C, 
        above_names=above_names, tree=tree, 
        Sigma = Sigma, Rho = Rho,
        above_names=above_names)
    )
}


get_binary <- function(abundance, rs){
    root_value <- rs[length(rs)]
    n_trait <- dim(abundance)[2]
    a <- log(abundance[, n_trait]) - root_value
    v <- rep(0, length(a))
    v[a>0] <- 1
    names(v) <- rownames(abundance)
    return(names(v[v>0]))
}


sim_wrapper <- function(tnum=10, bnum=8, filter_P=1, lambda=1, mu=1, rho=0.001, age=1, exponent=0){
    do_sim_res <- do_sim(tnum, bnum, lambda, mu, rho, age, exponent=exponent)
    abundance <- do_sim_res$abundance
    abun <- do_sim_res$abun
    prop <- do_sim_res$prop
    C <- do_sim_res$C
    tree <- do_sim_res$tree
    above_names <- do_sim_res$above_names; cat("above zero:\t", above_names, "\n")

    #return(list(abun=abun, prop=prop, C=C, tree=tree, above_names=above_names))
    return(do_sim_res)
}


check_BM <- function(abundance, tree, filter_P, mode){
    n <- dim(abundance)[2]
    phylo_sig <- data.frame(P=rep(1,n), K=rep(1,n))
    for(i in 1:dim(abundance)[2]){
        trait <- log(abundance[,i]) - log(sum(abundance[,1]))
        trait_data <- data.frame( species=rownames(abundance), trait=trait )
        trait_data$species <- factor(trait_data$species, levels = tree$tip.label)
        trait_values <- trait_data$trait[match(tree$tip.label, trait_data$species)]
        names(trait_values) <- tree$tip.label

        phylo_signal <- suppressMessages(
            phylosig_res <- phylosig(tree, trait_values, method = "lambda", test = TRUE)
        )
        fitBM <- fitContinuous(tree, trait_values, model = "BM")

        if(mode == 'BM'){
            lrt_p_value <- pchisq(2 * (phylosig_res$logL - fitBM$opt$lnL), 
                df=1, lower.tail = FALSE)
            phylo_sig$P[i] <- lrt_p_value
        } else if(mode == 'lambda'){
            phylo_sig$P[i] <- phylo_signal$P
            #phylo_sig$K[i] <- phylo_signal$K
        }
    }
    #names(phylo_sig) <- colnames(abundance)
    return(phylo_sig)
}


get_phylo_groups <- function(tree){
    #get the two descedants of the root
    root_node <- Ntip(tree) + 1   # the root node number
    children <- tree$edge[tree$edge[,1] == root_node, 2]
    child1_tips <- sort(tree$tip.label[unlist(Descendants(tree, children[1], type = "tips"))])
    child2_tips <- sort(tree$tip.label[unlist(Descendants(tree, children[2], type = "tips"))])
    group1 <- list(labels = child1_tips, col = 'green')
    group2 <- list(labels = child2_tips, col = 'purple')
    return(list(group1=group1, group2=group2))
}


simulate_discrete_trait <- function(tree, rate = 0.1) {
  # Ensure strictly bifurcating tree (defensive)
  if (!is.binary(tree)) {
    tree <- multi2di(tree)
  }

  # Binary CTMC rate matrix
  Q <- matrix(c(-rate, rate,
                rate, -rate), nrow = 2, byrow = TRUE)
  rownames(Q) <- colnames(Q) <- c("0", "1")

  # root.value is index in states, so use 1 (state "0")
  trait_states <- ape::rTraitDisc(
    phy = tree,
    model = Q,
    states = c("0", "1"),
    root.value = 1,
    ancestor = FALSE
  )

  names(trait_states[trait_states == "1"])
}


######################################
create_single_plot <- function(matrix,
                               method = "bray",
                               color = "black",
                               title = "PCoA Plot",
                               group_samples = NULL,
                               outlier_k = 0) {
  ## ---- checks ----
  if (!is.matrix(matrix) && !is.data.frame(matrix)) {
    stop("Input must be a matrix or data.frame.")
  }

  ## ---- PCoA ----
  diss <- vegan::vegdist(matrix, method = method)
  pts  <- cmdscale(diss, k = 2, eig = TRUE)

  coords <- as.data.frame(pts$points)
  colnames(coords) <- c("PC1", "PC2")
  coords$Sample <- rownames(matrix)
  coords$Color  <- color

  ## ---- % variance explained ----
  eigvals <- pts$eig
  prop_explained <- eigvals / sum(eigvals)
  pc1_lab <- paste0("PC1 (", round(prop_explained[1] * 100, 1), "%)")
  pc2_lab <- paste0("PC2 (", round(prop_explained[2] * 100, 1), "%)")

  ## ---- group colors ----
  if (!is.null(group_samples) && is.list(group_samples)) {
    for (grp in group_samples) {
      coords$Color[coords$Sample %in% grp$labels] <- grp$col
    }
  }

  ## ---- detect outliers in PCoA space ----
  centroid <- colMeans(coords[, c("PC1", "PC2")])
  coords$dist <- sqrt(
    (coords$PC1 - centroid["PC1"])^2 +
    (coords$PC2 - centroid["PC2"])^2
  )

  thr <- quantile(coords$dist, 0.75) +
         outlier_k * IQR(coords$dist)

  coords$outlier <- coords$dist > thr
    cat("Proportion of outliers:\t", sum(coords$outlier)/length(coords$Sample), "\n")

  ## ---- plot limits based on non-outliers ----
  xrange <- range(coords$PC1[!coords$outlier])
  yrange <- range(coords$PC2[!coords$outlier])

  ## ---- main plot (outliers removed) ----
  plot(coords$PC1[!coords$outlier],
       coords$PC2[!coords$outlier],
       pch = 19,
       col = coords$Color[!coords$outlier],
       xlab = pc1_lab,
       ylab = pc2_lab,
       main = title,
       asp = 1,
       xlim = xrange * 1.2,
       ylim = yrange * 1.2)

  ## ---- sample labels (non-outliers only) ----
  cex_of_host <- 0.8 * 10 / nrow(matrix)
  text(coords$PC1[!coords$outlier],
       coords$PC2[!coords$outlier],
       labels = coords$Sample[!coords$outlier],
       pos = 4,
       cex = cex_of_host)

  ## ---- optional: show outliers ----
  if (any(coords$outlier)) {
    points(coords$PC1[coords$outlier],
           coords$PC2[coords$outlier],
           pch = 4,
           col = "red",
           lwd = 2)
  }

  ## ---- group ellipses (non-outliers only) ----
  if (!is.null(group_samples)) {
    for (grp in group_samples) {
      sel <- coords$Sample %in% grp$labels & !coords$outlier
      if (sum(sel) >= 3) {
        vegan::ordiellipse(coords[sel, c("PC1", "PC2")],
                           rep(1, sum(sel)),
                           kind = "sd",
                           conf = 0.95,
                           draw = "polygon",
                           lwd = 1.2,
                           col = grp$col)
      }
    }
  }

  ## ---- return info for reproducibility ----
  invisible(coords)
}


plot_graphs <- function(m1, m2, m3, outfile, grp_list, outlier_k) {
    par(mfrow = c(2, 2))
    create_single_plot(m1, method = "bray", "green", title = "PCoA BC-distance", group_samples = grp_list, outlier_k)
    create_single_plot(m2, method = "euclidean", "cyan", title = "PCoA: Euclidean (CLR, standard)", group_samples = grp_list, outlier_k)
    create_single_plot(m3, method = "euclidean", "red", title = "PCoA: Euclidean \n(CLR, phylo decorrelated)", group_samples = grp_list, outlier_k)

    tip_colors <- rep(NA, length(tree$tip.label))
    #tip_colors <- ifelse(tree$tip.label %in% above_names, "orange", "blue")
    for (g in grp_list) {
        tip_colors[tree$tip.label %in% g$labels] <- g$col
    }
    plot(tree, tip.color = tip_colors, cex = 0.3, label.offset = 0.01)
}


calculate_correl_with_Rho <- function(Rho, pcoas, outfile) {
	eigvals <- eigen(Rho)$values
	sorted_eigvecs <- eigen(Rho)$vectors[, order(eigvals, decreasing = TRUE)]
	sorted_norm_Rho_eigens <- sort(eigvals, decreasing = TRUE) / sum(eigvals)
	len <- length(sorted_norm_Rho_eigens)
    len2 <- min(3, len)
    len3 <- min(3, len)

    title <- paste("pcoa_name", "Correlation_coefficient", paste("RMSD_top", len2, sep=''), paste("MRB_top", len3, sep=''), sep = "\t")
    write(title, file=outfile)

	for (i in seq_along(pcoas)) {
        pcoa <- pcoas[[i]]
        pcoa_norm <- pcoa$eig / sum(pcoa$eig)
        pcoa_name <- paste("pcoa_", i, sep="")
        # correl eigen value
        correl_eigval <- cor(sorted_norm_Rho_eigens[1:len], pcoa_norm[1:len])

        # correl eigen vector
        #correl_eigvec <- mean(sapply(1:len, function(i) { cor(sorted_eigvecs[, i], pcoa$points[, i]) }))

        # RMSD
        rmsd <- sqrt(sum((sorted_norm_Rho_eigens[1:len2] - pcoa_norm[1:len2])^2) / len2)

        # mean relative bias (see Gascuel's NC paper MRB)
        len3 <- min(3, len)
        mrb <- mean((pcoa_norm[1:len3] - sorted_norm_Rho_eigens[1:len3]) / sorted_norm_Rho_eigens[1:len3])
        output <- paste(pcoa_name, correl_eigval, rmsd, mrb, sep="\t")
        write(output, file=outfile, append=TRUE)
    }
}


write_adonis_results <- function(pcoa_res, grp_list, outdir, prefix) {
  adonis_dir <- file.path(outdir, "adonis")
  dir.create(adonis_dir, recursive = TRUE, showWarnings = FALSE)

  group_map <- unlist(lapply(names(grp_list), function(g) {
    setNames(rep(g, length(grp_list[[g]]$labels)), grp_list[[g]]$labels)
  }))

  sample_names <- rownames(pcoa_res$points)
  grp_factor <- factor(group_map[sample_names])
  if (any(is.na(grp_factor))) {
    stop("adonis: group labels are missing for some PCoA samples.")
  }

  coords_all <- as.data.frame(pcoa_res$points)
  coords_2 <- coords_all[, seq_len(min(2, ncol(coords_all))), drop = FALSE]

  run_adonis <- function(coords) {
    form <- stats::as.formula("coords ~ grp_factor")
    fit <- vegan::adonis2(form, data = list(coords = coords, grp_factor = grp_factor))
    fit_df <- as.data.frame(fit)
    fit_df$term <- rownames(fit_df)
    fit_df <- fit_df[fit_df$term != "Total" & fit_df$term != "Residual", c("term", "R2", "Pr(>F)")]
    colnames(fit_df) <- c("term", "R2", "P_value")
    fit_df
  }

  all_axes <- run_adonis(coords_all)
  first_two <- run_adonis(coords_2)

  write.table(
    all_axes,
    file = file.path(adonis_dir, paste0(prefix, "_all_axes.tbl")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    first_two,
    file = file.path(adonis_dir, paste0(prefix, "_pc1_pc2.tbl")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

write_adonis_results_for_groups <- function(pcoa_res, outdir, prefix, grp_list_phylo = NULL, grp_list_trait = NULL) {
  if (!is.null(grp_list_trait)) {
    write_adonis_results(pcoa_res, grp_list_trait, outdir, paste0(prefix, "_by_trait"))
  }
  if (!is.null(grp_list_phylo)) {
    write_adonis_results(pcoa_res, grp_list_phylo, outdir, paste0(prefix, "_by_phylo"))
  }
}


write_group_list <- function(grp_list, outdir, filename = "grp_list_by_phylo.tbl") {
  adonis_dir <- file.path(outdir, "adonis")
  dir.create(adonis_dir, recursive = TRUE, showWarnings = FALSE)

  group_df <- do.call(rbind, lapply(names(grp_list), function(g) {
    data.frame(
      group = g,
      sample = grp_list[[g]]$labels,
      color = grp_list[[g]]$col,
      stringsAsFactors = FALSE
    )
  }))

  write.table(
    group_df,
    file = file.path(adonis_dir, filename),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}


##################################
read_data <- function(C){
    # reorder
    #abundance <- abundance[rownames(abundance) %in% colnames(C), , drop = FALSE]
    abundance <- abundance[match(colnames(C), rownames(abundance)), ]
    # delete all zero columns
    abundance <- abundance[, apply(abundance, 2, function(col) any(col != 0))]
    abundance <- abundance+1e-5 # to avoid abundance of zero

    Z <- 1:nrow(C)
    abundance <- abundance[Z,]
    # phylo covariance matrix, vcv()
    C <- C[Z,Z]

    abun <- abundance #old name: Y
    prop <- abundance #old name: Q
    return(list(C=C, abundance=abundance, prop=prop))
}



############ Beta #################
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

build_beta_discrete_grid <- function(K = 8, eps = 1e-6) {
  edges <- seq(0, 1, length.out = K + 1)
  mids  <- (edges[-1] + edges[-(K + 1)]) / 2
  mids  <- pmin(pmax(mids, eps), 1 - eps)
  list(edges = edges, mids = mids)
}

beta_bin_weights <- function(a, b, edges, min_w = 1e-12) {
  w <- pbeta(edges[-1], shape1 = a, shape2 = b) -
       pbeta(edges[-length(edges)], shape1 = a, shape2 = b)
  w[w < min_w] <- min_w
  w / sum(w)
}

# returns matrix: rows=features, cols=lambda classes
precompute_loglik_matrix <- function(Y, C, lam_grid) {
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)
  p <- ncol(Y); K <- length(lam_grid)
  LL <- matrix(-Inf, nrow = p, ncol = K)
  colnames(LL) <- paste0("class_", seq_len(K))
  rownames(LL) <- colnames(Y)

  for (j in seq_len(p)) {
    y <- Y[, j]
    if (any(!is.finite(y)) || sd(y) < 1e-12) next
    for (k in seq_len(K)) {
      LL[j, k] <- logLik_pagel_lam_one_trait(y, C, lam_grid[k])
    }
  }
  LL
}


estimate_pagel_lam_hierarchical <- function(Y, C, K = 8, verbose = TRUE) {
  grid <- build_beta_discrete_grid(K = K)
  lam_grid <- grid$mids
  edges <- grid$edges

  LL <- precompute_loglik_matrix(Y, C, lam_grid)  # p x K
  ok_rows <- apply(LL, 1, function(z) any(is.finite(z)))
  if (!any(ok_rows)) stop("No valid traits for hierarchical pagel_lam fit.")

  LL_use <- LL[ok_rows, , drop = FALSE]

  # Optimize over log(a), log(b) to enforce positivity
  nll <- function(theta) {
    a <- exp(theta[1]); b <- exp(theta[2])
    w <- beta_bin_weights(a, b, edges)
    logw <- log(w)

    s <- 0
    for (j in seq_len(nrow(LL_use))) {
      z <- LL_use[j, ] + logw
      s <- s + logsumexp(z)
    }
    -s
  }

  fit <- optim(
    par = c(log(1), log(1)),
    fn = nll,
    method = "L-BFGS-B",
    lower = c(log(1e-3), log(1e-3)),
    upper = c(log(200), log(200))
  )

  a_hat <- exp(fit$par[1]); b_hat <- exp(fit$par[2])
  w_hat <- beta_bin_weights(a_hat, b_hat, edges)

  # posterior per feature
  p <- nrow(LL)
  K <- ncol(LL)
  post <- matrix(NA_real_, nrow = p, ncol = K)
  rownames(post) <- rownames(LL)
  colnames(post) <- paste0("class_", seq_len(K))

  post_mean <- rep(NA_real_, p)
  post_map  <- rep(NA_real_, p)
  names(post_mean) <- rownames(LL)
  names(post_map)  <- rownames(LL)

  logw <- log(w_hat)
  for (j in seq_len(p)) {
    z <- LL[j, ] + logw
    if (!any(is.finite(z))) next
    lse <- logsumexp(z)
    pj <- exp(z - lse)
    post[j, ] <- pj
    post_mean[j] <- sum(pj * lam_grid)
    post_map[j]  <- lam_grid[which.max(pj)]
  }

  if (verbose) {
    cat("Hierarchical pagel_lam fit (K=", K, "): a=", round(a_hat, 6),
        ", b=", round(b_hat, 6), "\n", sep = "")
  }

  list(
    a = a_hat, b = b_hat,
    lam_grid = lam_grid,
    bin_weights = w_hat,
    loglik_matrix = LL,
    posterior = post,
    pagel_lam_post_mean = post_mean,
    pagel_lam_post_map = post_map,
    converged = (fit$convergence == 0),
    optim = fit
  )
}


########### MLE pagel #############
make_C_pagel_lam <- function(C, pagel_lam) {
  C_lam <- C * pagel_lam
  diag(C_lam) <- diag(C)  # Pagel's lambda: only off-diagonals scaled
  C_lam <- as.matrix(Matrix::nearPD(C_lam)$mat)  # numerical safety
  return(C_lam)
}

logLik_pagel_lam_one_trait <- function(y, C, pagel_lam) {
  C_lam <- make_C_pagel_lam(C, pagel_lam)

  invC <- tryCatch(solve(C_lam), error = function(e) NULL)
  if (is.null(invC)) return(-Inf)

  n <- length(y)
  ones <- rep(1, n)

  # GLS root/mean estimate
  mu_hat <- as.numeric((t(ones) %*% invC %*% y) / (t(ones) %*% invC %*% ones))
  r <- y - mu_hat

  quad <- as.numeric(t(r) %*% invC %*% r)
  if (!is.finite(quad) || quad <= 0) return(-Inf)

  sigma2_hat <- quad / n
  logdet <- as.numeric(determinant(C_lam, logarithm = TRUE)$modulus)

  ll <- -0.5 * (n * log(2 * pi) + n * log(sigma2_hat) + logdet + n)
  return(ll)
}

estimate_pagel_lam_mle <- function(Y, C, lower = 1e-6, upper = 1) {
  if (is.vector(Y)) Y <- matrix(Y, ncol = 1)

  obj <- function(pagel_lam) {
    ll_sum <- 0
    for (j in seq_len(ncol(Y))) {
      y <- Y[, j]
      if (any(!is.finite(y)) || sd(y) < 1e-12) next
      ll <- logLik_pagel_lam_one_trait(y, C, pagel_lam)
      if (!is.finite(ll)) return(1e100)
      ll_sum <- ll_sum + ll
    }
    -ll_sum
  }

  fit <- optimize(obj, interval = c(lower, upper))
  return(fit$minimum)
}


do_transformation <- function(transform,
                              C,
                              log_prop_geomean,
                              use_pagel_lam = TRUE,
                              pagel_lam_mode = c("global", "per_feature", "hierarchical", "none"),
                              pagel_lam_lower = 1e-6,
                              pagel_lam_upper = 1,
                              hierarchical_K = 8,
                              verbose = TRUE) {
  pagel_lam_mode <- match.arg(pagel_lam_mode)

  transform_one <- function(y, C_use, transform) {
    C_inv <- solve(C_use)
    ones <- rep(1, nrow(C_use))

    root_value <- as.numeric((t(y) %*% C_inv %*% ones) / (t(ones) %*% C_inv %*% ones))

    if (grepl("chol", transform, ignore.case = TRUE)) {
      L_inv <- solve(t(chol(C_use)))
      x <- as.numeric(L_inv %*% (y - root_value))
      x <- x + root_value
    } else if (grepl("garland", transform, ignore.case = TRUE)) {
      eig <- eigen(C_use)
      L <- eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors)
      x <- as.numeric(L %*% (y - root_value))
      x <- x + root_value
    } else {
      stop("Unknown transform: use 'chol' or 'garland'")
    }
    return(list(x = x, root = root_value))
  }

  X <- matrix(NA_real_, nrow = nrow(log_prop_geomean), ncol = ncol(log_prop_geomean))
  rownames(X) <- rownames(log_prop_geomean)
  colnames(X) <- colnames(log_prop_geomean)

  pagel_lam_est <- rep(NA_real_, ncol(log_prop_geomean))
  names(pagel_lam_est) <- colnames(log_prop_geomean)

  if (!use_pagel_lam || pagel_lam_mode == "none") {
    for (j in seq_len(ncol(log_prop_geomean))) {
      y <- log_prop_geomean[, j]
      tr <- transform_one(y, C, transform)
      X[, j] <- tr$x
      pagel_lam_est[j] <- 1
    }

  } else if (pagel_lam_mode == "global") {
    pagel_lam <- estimate_pagel_lam_mle(log_prop_geomean, C, pagel_lam_lower, pagel_lam_upper)
    C_lam <- make_C_pagel_lam(C, pagel_lam)

    if (verbose) cat("Global pagel_lam MLE:", round(pagel_lam, 6), "\n")

    for (j in seq_len(ncol(log_prop_geomean))) {
      y <- log_prop_geomean[, j]
      tr <- transform_one(y, C_lam, transform)
      X[, j] <- tr$x
      pagel_lam_est[j] <- pagel_lam
    }

  } else if (pagel_lam_mode == "per_feature") {
    for (j in seq_len(ncol(log_prop_geomean))) {
      y <- log_prop_geomean[, j]
      pagel_lam_j <- estimate_pagel_lam_mle(y, C, pagel_lam_lower, pagel_lam_upper)
      C_lam_j <- make_C_pagel_lam(C, pagel_lam_j)

      tr <- transform_one(y, C_lam_j, transform)
      X[, j] <- tr$x
      pagel_lam_est[j] <- pagel_lam_j
    }

    if (verbose) {
      cat("Per-feature pagel_lam estimated.\n")
      cat("pagel_lam summary: min=", round(min(pagel_lam_est, na.rm = TRUE), 6),
          " median=", round(median(pagel_lam_est, na.rm = TRUE), 6),
          " max=", round(max(pagel_lam_est, na.rm = TRUE), 6), "\n", sep = "")
    }
  } else if (pagel_lam_mode == "hierarchical") {
    hfit <- estimate_pagel_lam_hierarchical(
      Y = log_prop_geomean, C = C, K = hierarchical_K, verbose = verbose
    )

    # choose posterior mean (smooth shrinkage)
    pagel_lam_est <- hfit$pagel_lam_post_mean

    for (j in seq_len(ncol(log_prop_geomean))) {
      y <- log_prop_geomean[, j]
      lam_j <- pagel_lam_est[j]
      if (!is.finite(lam_j)) lam_j <- 1
      C_lam_j <- make_C_pagel_lam(C, lam_j)
      tr <- transform_one(y, C_lam_j, transform)
      X[, j] <- tr$x
    }

    return(list(
      X = X,
      pagel_lam = pagel_lam_est,
      pagel_lam_map = hfit$pagel_lam_post_map,
      pagel_lam_hyper_a = hfit$a,
      pagel_lam_hyper_b = hfit$b,
      pagel_lam_grid = hfit$lam_grid,
      pagel_lam_bin_weights = hfit$bin_weights,
      pagel_lam_posterior = hfit$posterior,
      hierarchical_converged = hfit$converged
    ))
  }

  return(list(
    X = X,
    pagel_lam = pagel_lam_est
  ))
}


##################################
get_grp_info <- function(grp_infile, species_exclude){
    grp_info <- read.table(grp_infile, fill = TRUE, stringsAsFactors = FALSE, comment.char = "")
    grp_info <- grp_info[! rownames(grp_info) %in% species_exclude, , drop = FALSE]
    return(grp_info)
}


make_grp_list <- function(df, trait_col, colors = NULL) {
  if (!is.data.frame(df)) {
    stop("df must be a data.frame")
  }
  if (is.null(rownames(df))) {
    stop("Row names are required and will be used as labels")
  }
  if (!trait_col %in% colnames(df)) {
    stop("Trait column not found: ", trait_col)
  }

  labels <- rownames(df)
  traits <- df[[trait_col]]
  split_groups <- split(labels, traits, drop = TRUE)

  if (is.null(colors)) {
    colors <- grDevices::rainbow(length(split_groups))
  }

  grp_list <- Map(
    function(lbls, col) {
      list(labels = lbls, col = col)
    },
    split_groups,
    colors
  )

  return(grp_list)
}


##################################
# some params to change, XQ
lambda <- 5
mu <- 5
rho <- 0.001
age <- 7.5
is_one <- FALSE
is_sim <- FALSE
is_inter <- FALSE

# BM means the default is pagel's lambda = 1 (standard BM) and those with p-value in a LRT >=0.05 is accepted; "lambda" means that the default is lambda = 0, and only p-value <0.05 is accepted.
# i will modify this later to enable pagel's lambda, which is better for empirical data.
# Now it works only under simulations where log(prop) perfectly follows standard BM.
filter_mode <- 'BM'
filter_P <- 1
is_check <- FALSE

# control whether features (microbial abundance) are correlated (0) or not (big value), see Rho
exponent <- 0

# host tree tip num
tnum <- 10
# bac (microbiota) taxa num
bnum <- 8

transform <- 'garland'
dist_method <- 'euclidean'
is_standardize <- FALSE
outlier_k <- 100
species_exclude <- c()

outdir <- NULL
is_force <- FALSE

grp_info <- NULL


##################################
spec = matrix(c(
    'tree', 't', 2, "character",
    'abundance', 'a', 2, "character",
    'metadata', 'm', 2, "character",
    'one', 'O', 0, "logical",
    'all', 'A', 0, "logical",
    #'transform', 'T', 2, "character",

    'sim', 's', 0, "logical",
    'check', 'c', 2, "character",
    'filter_P', 'p', 2, "double",

    'grp', 'g', 2, 'character',
    'feature', 'f', 2, 'character',

    'exponent', 'e', 2, "double",
    'tnum', 'T', 2, "integer",
    'bnum', 'B', 2, "integer",
    'bd', 'b', 2, "character",
    'dist', 'd', 2, 'character',
    'inter', 'i', 0, "logical",
    'standardize', 'S', 0, "logical",
    'outlier', '', 2, 'double',
    'species_exclude', '_', 2, 'character',
    'sim_discrete_trait', 'D', 0, 'logical',
    'binary_rate', 'r', 2, 'double',

    'help' , 'h', 0, "logical",
    'outdir', 'o', 1, "character",
    'force', 'NA', 0, 'logical'
), byrow=TRUE, ncol=4)

opt <- getopt(spec)


##################################
if (!is.null(opt$one)) {
    is_one <- TRUE
}
if (!is.null(opt$all)) {
    is_one <- FALSE
}

if(!is.null(opt$sim)){
    is_sim = TRUE
}

if(!is.null(opt$inter)){
    is_inter = TRUE
}

if (!is.null(opt$help)) {
    cat(getopt(spec, usage = T))
    q(status = 1)
}

if(! is.null(opt$tree)){
    tree <- read.tree(opt$tree)
    C <- vcv(tree)
    rownames(C) <- colnames(C) <- colnames(vcv(tree))
}

if(! is.null(opt$abundance)){
    abundance <- read.table(opt$abundance, header=TRUE, comment.char = "")
}

if (!is.null(opt$metadata)) {
    if (is_one){
        selected_samples <- generate_metadata(opt$metadata)
        overlapping_items <- intersect(colnames(abundance), selected_samples$sample_id)
        new_df <- abundance[, colnames(abundance) %in% overlapping_items]
        name_mapping <- setNames(selected_samples$species, selected_samples$sample_id)
        colnames(new_df) <- name_mapping[colnames(new_df)]
        abundance <- t(new_df)
    } else {
        abundance <- generate_metadata2(opt$metadata, abundance)
        abundance <- t(abundance)
    }
}

if(! is.null(opt$transform)){
    transform <- opt$transform
}

if(! is.null(opt$check)){
    check <- opt$check
    checks <- strsplit(check, ",")[[1]]
    filter_mode <- checks[1]
    filter_P <- checks[2]
    is_check = TRUE
}
# probably disabled later
if(! is.null(opt$filter_P)){
    filter_P <- opt$filter_P
}

if(! is.null(opt$exponent)){
    exponent <- opt$exponent
}

if(! is.null(opt$tnum)){
    tnum <- opt$tnum
}
if(! is.null(opt$bnum)){
    bnum <- opt$bnum
}

if(! is.null(opt$bd)){
    bd_param <- opt$bd
    params <- as.numeric(strsplit(bd_param, ",")[[1]])
    lambda <- params[1]
    mu     <- params[2]
    rho    <- params[3]
    age    <- params[4]
}

if(! is.null(opt$dist)){
    dist_method <- opt$dist
}

if(! is.null(opt$standardize)){
    is_standardize <- TRUE
}

if (!is.null(opt$outlier)) {
    outlier_k <- opt$outlier
}
if (!is.null(opt$species_exclude)) {
    species_exclude <- as.vector(read.table(opt$species_exclude)[, 1])
}

if(! is.null(opt$grp)){
    grp_info <- get_grp_info(opt$grp, species_exclude)
}
if(! is.null(opt$feature)){
    feature <- opt$feature
}

sim_discrete_trait <- FALSE
if (!is.null(opt$sim_discrete_trait)) {
    sim_discrete_trait <- TRUE
}
if (!is.null(opt$binary_rate)) {
    binary_rate <- opt$binary_rate
}

if(! is.null(opt$force)){
    is_force <- TRUE
}

if(! is.null(opt$outdir)){
    outdir <- opt$outdir
    if (dir.exists(outdir)) {
      if (is_force) {
        message("removing outdir ", outdir, ' ......')
        unlink(outdir, recursive = TRUE, force = TRUE)
      } else {
        stop("outdir ", outdir, " already exists. Use --force. Exiting ......")
      }
    }
    dir.create(outdir, recursive = TRUE)
}


##################################
#---------- start here ----------#
##################################
if(! is_sim){
    common <- intersect(colnames(C), rownames(abundance))
    if(length(setdiff(tree$tip.label, common)) > 0){
        print(colnames(C))
        cat("\n")
        print(rownames(abundance))
        tree <- drop.tip(tree, setdiff(tree$tip.label, common))
    }
    if(length(species_exclude) > 0){
        tree <- drop.tip(tree, species_exclude)
    }
    C <- vcv(tree)
    #C <- C[, common, drop = FALSE]
    #C <- C[common, , drop = FALSE]
    abundance <- abundance[common, , drop = FALSE]
    data_res <- read_data(C)
    prop <- data_res$prop
    prop <- prop[, colMeans(prop)>0.01, drop=FALSE]
    abundance <- data_res$abundance
} else{
    #col: bac taxa (bnum), row: host species (tnum)
    sim_res <- do_sim(tnum, bnum, lambda, mu, rho, age, exponent)
    C <- sim_res$C #phylo_cov
    tree <- sim_res$tree

    prop <- sim_res$prop
    above_names <- sim_res$above_names

    abundance <- sim_res$abundance
    Sigma <- sim_res$Sigma
    Rho <- sim_res$Rho
}

trait_group_names <- above_names
if (is_sim && isTRUE(sim_discrete_trait)) {
    trait_group_names <- simulate_discrete_trait(tree, rate = binary_rate)
    cat("using simulated discrete trait with", length(trait_group_names), "state-1 tips at rate ", binary_rate, "\n", sep = "")
}


if(is_check){
    phylo_sig <- check_BM(prop, tree, filter_P, filter_mode)
    if(filter_mode == 'lambda'){
        selected_cols <- which(phylo_sig$P < filter_P)
    } else if(filter_mode == 'BM'){
        selected_cols <- which(phylo_sig$P > filter_P)
    }
    if(length(selected_cols)<3){
        stop("check BM: passed features < 3; try a larger number for --bnum; Exiting ......")
    }
    cat("checking BM with the mode", filter_mode, "\n")
    cat(length(selected_cols), "out of", ncol(prop), "are kept")
    cat("\n\n")
    abundance <- abundance[, selected_cols, drop = FALSE]
    prop <- prop[, selected_cols, drop = FALSE]

    if(is_sim && "Sigma" %in% names(sim_res)){
        #abundance <- abundance[, selected_cols, drop = FALSE]
        Sigma <- Sigma[selected_cols, selected_cols]
        Rho <- Rho[selected_cols, selected_cols]
    }

    trait_group_names <- intersect(trait_group_names, rownames(prop))
}


###################################################
# reorder
#Y <- Y[match(colnames(C), rownames(Y)), ]

log_prop <- log(prop+1e-6)
row_geomean_log_prop <- rowMeans(log_prop)
log_prop_geomean <- log_prop

for(i in 1:dim(prop)[1]) # iterate host
{
    log_prop_geomean[i,] <- log(prop[i,]) - row_geomean_log_prop[i]
}


##################################
# P is the transformed matrix
trans_res <- do_transformation(
  transform = transform,
  C = C,
  log_prop_geomean = log_prop_geomean,
  use_pagel_lam = TRUE,
  pagel_lam_mode = "hierarchical" # "per_feature"
)

P <- trans_res$X
rownames(P) <- rownames(C)
if (is.null(colnames(P))) {
  colnames(P) <- colnames(log_prop_geomean)
}
if (is.null(colnames(P))) {
  colnames(P) <- paste0("feature_", seq_len(ncol(P)))
}

pagel_lam_tbl <- data.frame(
  feature = colnames(P),
  pagel_lam_mle = trans_res$pagel_lam
)

write.table(
  pagel_lam_tbl,
  file = file.path(outdir, "pagel_lam_estimates.tbl"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


rounded_to <- 3
write.table(round(prop, rounded_to), file = file.path(outdir, 'prop.tbl'), sep = "\t", quote = FALSE, row.names=TRUE)
write.table(round(log_prop, rounded_to), file = file.path(outdir, 'log_prop.tbl'), sep = "\t", quote = FALSE, row.names=TRUE)
write.table(round(log_prop_geomean, rounded_to), file = file.path(outdir, 'log_prop_geomean.tbl'), sep = "\t", quote = FALSE, row.names=TRUE)
write.table(round(P, rounded_to), file = file.path(outdir, 'P.tbl'), sep = "\t", quote = FALSE, row.names=TRUE)

if(is_sim){
    write.table(round(Sigma, rounded_to), file = file.path(outdir, 'Sigma.tbl'), sep = "\t", quote = FALSE, row.names=TRUE)
    write.table(round(Rho, rounded_to), file = file.path(outdir, 'Rho.tbl'), sep = "\t", quote = FALSE, row.names=TRUE)
}


##################################
# output file pdf
outfile <- file.path(outdir, "pcoa.pdf")
pdf(outfile)

#rownames(P) <- rownames(log_prop_geomean)
# normal pcoa, not phylo corrected
if(!is_sim){
    if (! is.null(grp_info)){
        grp_list <- make_grp_list(grp_info, feature)
    } else{
        trait_group_names <- c()
    }
} else{
    grp_list <- list(
        group1=list('labels'=trait_group_names, 'col'='orange'), 
        group2=list('labels'=setdiff(tree$tip.label,trait_group_names),'col'='blue')
    )
}

# BC
pcoa_0 <- calculate_pcoa(prop, 'bray', is_standardize, grp_list)
# standard PCA
pcoa_1 <- calculate_pcoa(log_prop_geomean, dist_method, is_standardize, grp_list)
# phylo corrected
if(! is_inter){
    pcoa_2 <- calculate_pcoa(P, dist_method, is_standardize, grp_list)
    plot_graphs(prop, log_prop_geomean, P, outfile, grp_list, outlier_k)
}


##################################
# generate the groups that are determined by the two descendant lineages of the root
grp_list_by_phylo <- get_phylo_groups(tree)
plot_graphs(prop, log_prop_geomean, P, outfile, grp_list_by_phylo, outlier_k)


##################################
pcoas <- list(pcoa_0, pcoa_1, pcoa_2)

write_group_list(grp_list_by_phylo, outdir)
write_adonis_results_for_groups(pcoa_0, outdir, "pcoa_0", grp_list_by_phylo, grp_list)
write_adonis_results_for_groups(pcoa_1, outdir, "pcoa_1", grp_list_by_phylo, grp_list)
if (!is_inter) {
    write_adonis_results_for_groups(pcoa_2, outdir, "pcoa_2", grp_list_by_phylo, grp_list)
}

compare_outdir <- file.path(outdir, "compare")
dir.create(compare_outdir, recursive = TRUE)

# output eigenvalues and the indices of separating groups by major PCs
for (i in seq_along(pcoas)){
    rounded_to <- 3
    pcoa <- pcoas[[i]]
    explained_file <- file.path(compare_outdir, 'explained.tbl')
    # Compute the rounded, normalized eigenvalues
    vals <- round((pcoa$eig / sum(pcoa$eig))[1:min(bnum, tnum)], rounded_to)
    # Write one record per line, separated by spaces, append to file
    write( vals, file = explained_file, ncolumns = length(vals), append = TRUE, sep = ' ')

    pcoa_name <- paste("pcoa_", i, sep='')

    # grp_by_trait
    title <- paste("pcoa_name", "LDA_acc", "fdr_value", "DBI", sep="\t")
    determined_by_trait_outfile <- file.path(compare_outdir, "determined_by_trait.tbl")
    if(i == 1){
        write(title, file=determined_by_trait_outfile, append=TRUE, sep="\t")
    }
    check_clustering(pcoa, pcoa_name=pcoa_name, grp_list=grp_list, outfile=determined_by_trait_outfile)

    # grp_by_phylo
    determined_by_phylo_outfile <- file.path(compare_outdir, "determined_by_phylo.tbl")
    if(i == 1){
        write(title, file=determined_by_phylo_outfile, append=TRUE, sep="\t")
    }
    check_clustering(pcoa, pcoa_name=pcoa_name, grp_list_by_phylo, outfile=determined_by_phylo_outfile)
}

dev.off()

if(!is_sim){
    q()
}

# compare the PC to the true R matrix (trait_cov: Sigma)
calculate_correl_with_Rho(Rho, pcoas, file.path(compare_outdir, "compared_to_R_matrix.tbl"))



##################################
files <- list.files(file.path(outdir, "compare/"), full.names = TRUE)
for (f in files) {
    cat(f,"\n")
    cat(readLines(f), sep = "\n")
    cat("\n")
}
