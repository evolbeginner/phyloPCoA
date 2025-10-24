#! /bin/env Rscript

#Sishuo's idea on accounting phylogeny for proportion data


###############################################################
# Update history
# 2025-10-21
#   greatly improved!
#   now prop, abun, log_prop_geomean to replace Q, Y, P
# 2025-07-18
#   some wrapping of the code to make it easier to read
#   -T 30 -B 8 -p 0.2
# 2025-07-16_2
#   some more updates for XQ
# 2025-07-16_1
#   some small update for XQ
# 2024-11-21
#   pcoa_plot
#   check_BM
# 2024-11-14
#   --all, --one
# 2024-10-29
#   pcoa sim
# 2024-10-26
#   Expectation transformed_X <- t(t(X)-root_values)
#   Cholesky decomposition corrected!
# 2024-10-25:
#   Expectation after phylo_var transform updatd (Sishuo)
#   pb() to generate traits
#   TreeSim to generate random tree, where lambda and mu to ctrl for "phylogenetic relateness"
#   X <- L %*% YY instead of t(L) for Cholesky decomposition
#   C^(-0.5) %*% X may be also fine
# 2024-09:
#   Initiated (Youhua Chen)


###############################################################
#set.seed(2)
suppressPackageStartupMessages({
    library(getopt)
    library(dplyr)
    library(tidyr)

    library(TreeSim)
    library(mvMORPH)
    library(phytools)
    library(adephylo)
    library(vegan)
    library(labdsv) #pco
    library(compositions) #clr
})


##################################
calculate_pcoa <- function(X, method = "bray", scaled=FALSE) {
    if(scaled == TRUE){X <- scale(X)}
    distance <- vegdist(X, method = method)
    pcoa_results <- cmdscale(distance, eig = T, k = 2)
}


generate_metadata <- function(metadata_file){
    metadata <- read.table(metadata_file, header = T)
    # Randomly select one sample_id for each species
    selected_samples <- metadata %>%
        group_by(species) %>%
        slice(1) #always the 1st, reproducible
        #sample_n(1)  # Sample one row from each group
    return(selected_samples)
}

generate_metadata2 <- function(metadata_file, abundance){
    metadata <- read.table(metadata_file, header = T)
    abundance$taxon <- rownames(abundance)
    rownames(abundance) <- NULL
    abundance_long <- abundance %>%
        pivot_longer(cols = -taxon, names_to = "sample_id", values_to = "abundance")
    merged_data <- inner_join(abundance_long, metadata, by = "sample_id")
    merged_data <- merged_data %>%
        group_by(sample_id) %>%
        mutate(abundance = abundance / sum(abundance)) %>%
        ungroup()
    #write.table(merged_data, "prop_abundance.txt")

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


# XQ
# bnum: bac taxonomic unit # in the microbiota
# tnum: tip (host) species
do_sim <- function(tnum, bnum, lambda, mu, rho, age){
    # BD model
    trees <- sim.bd.taxa.age(tnum, 1, lambda, mu, rho, age)
    tree <- trees[[1]]
    write.tree(tree, file=file.path(outdir, 'sim.tree'))

    # diff root_values
    mean_a=0; sd_a=1; nsim=bnum
    abundance <- matrix(0, nrow = tnum, ncol = nsim); 

    # generate matrix R (randomly)
    scale_R <- 1/distRoot(tree)[1]
    # XQ, R controls sigma in mvSIM() (covariance matrix)
    R <- crossprod(matrix(runif(bnum*bnum),bnum)) * scale_R
    # last more correlation
    #R[1:5,ncol(R)] <- 1/rbeta(5, 1, 1) * R[1:5,ncol(R)]
    #R[nrow(R),1:5] <- 1/rbeta(5, 10, 1) * R[nrow(R),1:5]

    rs <- rnorm(ncol(R), mean_a, sd_a) # diff root_values, XQ
    abundance <- exp(mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R, theta=rs)))
    above_names <- get_binary(abundance, rs)

    a_values <- rnorm(nsim, mean = mean_a, sd = sd_a)

    C <- vcv(tree)

    abun <- abundance
    prop <- t(apply(abundance, 1, function(x) x/sum(x)))

    return(list(abundance=abundance, abun=abun, prop=prop, C=C, above_names=above_names, tree=tree))
}


get_binary <- function(abundance, rs){
    a <- log(abundance[,dim(abundance)[2]]) - rs[length(rs)]
    v <- rep(0, length(a))
    v[a>0] <- 1
    names(v) <- rownames(abundance)
    return(names(v[v>0]))
}


sim_wrapper <- function(tnum=10, bnum=8, filter_BM_P=1, lambda=1, mu=1, rho=0.001, age=1){
    do_sim_res <- do_sim(tnum, bnum, lambda, mu, rho, age)
    abundance <- do_sim_res$abundance
    abun <- do_sim_res$abun
    prop <- do_sim_res$prop
    C <- do_sim_res$C
    tree <- do_sim_res$tree
    above_names <- do_sim_res$above_names; cat("above zero:\t", above_names, "\n")

    # XQ: set P-value for brownian motion checking (e.g., "-p 0.1")
    if(filter_BM_P < 1){
        phylo_sig <- check_BM(abundance, tree)
        selected_cols <- which(phylo_sig$abun < filter_BM_P)
        abundance <- abundance[, selected_cols, drop = FALSE]
        abun <- abun[, selected_cols, drop = FALSE]
        prop <- prop[, selected_cols, drop = FALSE]
    }

    return(list(abun=abun, prop=prop, C=C, tree=tree))
}


check_BM <- function(abundance, tree){
    n <- dim(abundance)[2]
    phylo_sig <- data.frame(P=rep(1,n), K=rep(1,n))
    for(i in 1:dim(abundance)[2]){
        trait <- log(abundance[,i])-log(sum(abundance[,1]))
        trait_data <- data.frame( species=rownames(abundance), trait=trait )
        trait_data$species <- factor(trait_data$species, levels = tree$tip.label)
        trait_values <- trait_data$trait[match(tree$tip.label, trait_data$species)]
        phylo_signal <- phylosig(tree, trait_values, method = "K", test = TRUE)
        phylo_sig$P[i] <- phylo_signal$P
        phylo_sig$K[i] <- phylo_signal$K
    }
    #names(phylo_sig) <- colnames(abundance)
    print(phylo_sig)
    return(phylo_sig)
}


plot_graphs <- function(m1, m2, m3, outfile) {
    pdf(outfile)
    par(mfrow = c(2, 2))
    create_single_plot <- function(matrix, method, color, title) {
        # Compute distance
        if (method == "euclidean") {
            diss <- vegdist(matrix, method = method)
        } else {
            diss <- vegdist(matrix, method = method)
        }
        
        # Perform PCoA
        pts <- cmdscale(diss, k = 2)
        
        plot(pts[, 1], pts[, 2],
             pch = 3,
             xlab = "PC1",
             ylab = "PC2",
             main = title,
             col = color)
        text(pts[, 1], pts[, 2],
             labels = if (!is.null(rownames(pts))) rownames(pts) else 1:nrow(pts),
             pos = 4,
             cex = 0.35)
    }
    
    # Create the three plots
    create_single_plot(m1, "bray", "green", "PCoA: Bray-Curtis (m1)")
    create_single_plot(m2, "euclidean", "blue", "PCoA: Euclidean (m2)")
    create_single_plot(m3, "euclidean", "red", "PCoA: Euclidean (m3)")
    dev.off()

    # Reset the plotting layout
    par(mfrow = c(1, 1))
}


##################################
read_data <- function(C){
    # reorder
    abundance <- abundance[match(colnames(C), rownames(abundance)), ]
    # delete all zero columns
    abundance <- abundance[, apply(abundance, 2, function(col) any(col != 0))]
    abundance <- abundance+1e-5 # to avoid abundance of zero

    Z <- 1:length(tree$tip.label)
    abundance <- abundance[Z,]
    # phylo covariance matrix, vcv()
    C <- C[Z,Z]

    abun <- abundance #old name: Y
    prop <- abundance #old name: Q
    return(list(C=C, abun=abun, prop=prop))
}


do_transformation <- function(transform, C, log_prop_geomean){
    # calculate the inverse of C in order to calculate the root values
    C_inv <- solve(C)
    # calculate the phylo mean (root value) for each node
    ones <- rep(1, dim(C)[1])
    root_values <- apply(log_prop_geomean, 2, function(col) {
        numerator <- col %*% C_inv %*% ones
        numerator / (t(ones) %*% C_inv %*% ones)
    })

    if (grepl("chol", transform, ignore.case = T)){
        #de-correlation by Cholesky decomposition
        L <- solve(t(chol(C)))
        L <- L * sqrt(diag(C)[1])
        X <- L %*% t((t(log_prop_geomean) - root_values))
        X <- t(t(X) + root_values)
    } else if(grepl("garland", transform, ignore.case = T)){
        # using C^(-0.5)
        eig_decomp <- eigen(C)
        L <- eig_decomp$vectors %*% diag(1/sqrt(eig_decomp$values)) %*% t(eig_decomp$vectors)
        L <- L * sqrt(diag(C)[1])
        X <- L %*% t((t(log_prop_geomean) - root_values))
        X <- t(t(X) + root_values)
    }
    return(X)
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
filter_BM_P <- 1

# host tree tip num
tnum <- 10
# bac (microbiota) taxa num
bnum <- 8

transform <- "garland"
dist_method <- 'euclidean'
is_standardize <- FALSE

outdir <- NULL
is_force <- FALSE


##################################
spec = matrix(c(
    'tree', 't', 2, "character",
    'abundance', 'a', 2, "character",
    'metadata', 'm', 2, "character",
    'one', 'O', 0, "logical",
    'all', 'A', 0, "logical",
    #'transform', 'T', 2, "character",
    'sim', 's', 0, "logical",
    'filter_BM_P', 'p', 2, "double",
    'tnum', 'T', 2, "integer",
    'bnum', 'B', 2, "integer",
    'bd', 'b', 2, "character",
    'dist', 'd', 2, 'character',
    'inter', 'i', 0, "logical",
    'help' , 'h', 0, "logical",
    'outdir', 'o', 1, "character",
    'force', 'f', 0, 'logical'
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
    abundance <- read.table(opt$abundance, header=T)
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

if(! is.null(opt$filter_BM_P)){
    filter_BM_P <- opt$filter_BM_P
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
    prop <- read_data(C)$prop
} else{
    sim_wrapper_result <- sim_wrapper(tnum, bnum, filter_BM_P, lambda, mu, rho, age)
    #col: bac taxa (bnum), row: host species (tnum)
    abun <- sim_wrapper_result$abun
    prop <- sim_wrapper_result$prop
    C <- sim_wrapper_result$C #phylo_cov
    tree <- sim_wrapper_result$tree
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
P <- do_transformation(transform, C, log_prop_geomean)
rownames(P) <- rownames(C)
cat("\n\n")


rounded_to <- 4
write.table(round(prop, rounded_to), file = file.path(outdir, 'prop.tbl'), sep = "\t", quote = FALSE)
write.table(round(log_prop, rounded_to), file = file.path(outdir, 'log_prop.tbl'), sep = "\t", quote = FALSE)
write.table(round(log_prop_geomean, rounded_to), file = file.path(outdir, 'log_prop_geomean.tbl'), sep = "\t", quote = FALSE)
write.table(round(P, rounded_to), file = file.path(outdir, 'P.tbl'), sep = "\t", quote = FALSE)


##################################
#rownames(P) <- rownames(log_prop_geomean)
rounded_to <- 3

cat("Original composition", "\n")
column_sorted_ind <- order( colMeans(prop, na.rm=TRUE), decreasing=TRUE )
range <- 1:10
cat("prop: before transformation (raw)", "\n")
column_sorted_ind_2 <- column_sorted_ind[range][column_sorted_ind[range] != bnum]
round(prop[, column_sorted_ind_2], rounded_to)

cat("\n")
cat("log_prop_geomean: minus geomean followed by transformation", "\n")
round(P[, column_sorted_ind_2], rounded_to)
cat("\n\n\n")

# normal pcoa, not phylo corrected
pcoa_1 <- calculate_pcoa(log_prop_geomean, dist_method, is_standardize)

# phylo corrected
if(! is_inter){
    pcoa_2 <- calculate_pcoa(P, dist_method, is_standardize)
    plot_graphs(prop, log_prop_geomean, P, file.path(outdir, "pcoa.pdf"))
}

cat(round(pcoa_1$eig/sum(pcoa_1$eig), 3), "\n", round(pcoa_2$eig/sum(pcoa_2$eig), 3), "\n")


