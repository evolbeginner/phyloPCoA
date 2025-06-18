#! /bin/env Rscript

#Sishuo's idea on accounting phylogeny for proportion data

###############################################################
#set.seed(2)
# Update history
# 2025-6-12
#   some bugs fixed for do_sim() and filter_brown()
# 2024-11-21
#   pcoa_plot
#   filter_brown
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
library(getopt)
library(dplyr)
library(tidyr)

library(TreeSim)
library(mvMORPH)
library(phytools)
library(adephylo)
library(vegan)
library(labdsv) #pco


##################################
calculate_pcoa <- function(X, method = "bray") {
    distance <- vegdist(X, method = method)
    pcoa_results <- cmdscale(distance, eig = T, k = 2)
}


generate_metadata <- function(metadata_file){
    metadata <- read.table(metadata_file, header = T)
    # Randomly select one sample_id for each species
    selected_samples <- metadata %>%
        group_by(species) %>%
        sample_n(1)  # Sample one row from each group
    return(selected_samples)
}

generate_metadata2 <- function(metadata_file, abundance){
    metadata <- read.table(metadata_file, header = T)
    abundance$taxon <- rownames(abundance)
    rownames(abundance) <- NULL

    abundance_long <- abundance %>%
        pivot_longer(cols = -taxon, names_to = "sample_id", values_to = "abundance")

    merged_data <- inner_join(abundance_long, metadata, by = "sample_id")

    average_abundance <- merged_data %>%
      group_by(species, taxon) %>%
      summarize(average_abundance = mean(abundance), .groups = 'drop')
      #summarize(average_abundance = exp(mean(log(abundance+1e-5), na.rm = TRUE)), .groups = 'drop')

    # Reshape the average_abundance data frame to wide format
    average_abundance_wide <- average_abundance %>%
      pivot_wider(names_from = species, values_from = average_abundance)

    average_abundance_wide = as.data.frame(average_abundance_wide)
    rownames(average_abundance_wide) <- average_abundance_wide[,1]
    average_abundance_wide <- average_abundance_wide[,-1]

    return(average_abundance_wide)
}


do_sim <- function(snum=10, hnum=50){
    # Do simulation
    #trees <- sim.bd.taxa.age(snum, 1, lambda=1, mu=0, rho=1, age=1); tree <- trees[[1]]
    trees <- sim.bd.taxa.age(snum, 1, lambda, mu, rho, age)
    tree <- trees[[1]]
    write.tree(tree, file='sim.tree')

    # sim trait of mu = 1, and sig2=1
    #abundance <- exp(fastBM(tree, sig2 = 1, a = 1, nsim=hnum))
    # diff root_values
    mean_a=0; sd_a=1; nsim=hnum
    abundance <- matrix(0, nrow = snum, ncol = nsim); 

    # generate matrix R (randomly)
    scale_R <- 1/distRoot(tree)[1]
    R <- crossprod(matrix(runif(hnum*hnum),hnum)) * scale_R
    # last more correlation
    #R[1:5,ncol(R)] <- 1/rbeta(5, 1, 1) * R[1:5,ncol(R)]
    #R[nrow(R),1:5] <- 1/rbeta(5, 10, 1) * R[nrow(R),1:5]

    rs <- rnorm(ncol(R), mean_a, sd_a) # diff root_values
    abundance <- exp(mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R, theta=rs)))

    above_names <- get_binary(abundance, rs)

    a_values <- rnorm(nsim, mean = mean_a, sd = sd_a)
    #for (i in 1:nsim){ abundance[,i] <- exp(fastBM(tree, sig2 = 1, a = a_values[i], nsim=1)) }

    # generate the phylo covariance matrix, C
    C=vcv(tree)

    # conversion
    Y <- abundance
    Q <- abundance

    return(list(abundance=abundance, Y=Y, Q=Q, C=C, above_names=above_names, tree=tree))
}


get_binary <- function(abundance, rs){
    a <- log(abundance[,dim(abundance)[2]]) - rs[length(rs)]
    v <- rep(0, length(a))
    v[a>0] <- 1
    names(v) <- rownames(abundance)
    return(names(v[v>0]))
}


filter_brown <- function(abudance, tree){
    n <- dim(abundance)[2]
    phylo_sig <- data.frame(P=rep(1,n), K=rep(1,n))
    for(i in 1:dim(abundance)[2]){
        trait <- log(abundance[,i])-log(sum(abundance[,1])+10^(10))
        trait_data <- data.frame( species=rownames(abundance), trait=trait )
        trait_data$species <- factor(trait_data$species, levels = tree$tip.label)
        trait_values <- trait_data$trait[match(tree$tip.label, trait_data$species)]
        phylo_signal <- phylosig(tree, trait_values, method = "K", test = TRUE)
        phylo_sig$P[i] <- phylo_signal$P
        phylo_sig$K[i] <- phylo_signal$K
    }
    #names(phylo_sig) <- colnames(abundance)
    return(phylo_sig)
}


plot_graph <- function(m1, m2, abundance){
    cex <- 0.35
    pcoa.counts1 <- m1  # Replace with your first data matrix
    pcoa.counts2 <- m2  # Replace with your second data matrix
    pcoa.bc1 <- pco(dsvdis(pcoa.counts1, index = "bray/curtis"))
    pcoa.bc2 <- pco(dsvdis(pcoa.counts2, index = "bray/curtis"))
    x_range <- range(c(pcoa.bc1$points[, 1], pcoa.bc2$points[, 1])) * 1.1
    y_range <- range(c(pcoa.bc1$points[, 2], pcoa.bc2$points[, 2])) * 1.1

    #under_zero = as.vector(rownames(abundance[abundance[,ind] <  ave, ]))
    #above_zero = as.vector(rownames(abundance[abundance[,ind] >= ave, ]))
    #print(above_zero)

    # Create the initial plot for the first dataset
    plot(pcoa.bc1$points[, 1], pcoa.bc1$points[, 2],
         pch = 3, xlab = "PC1", ylab = "PC2", main = "PCoA: Bray-Curtis", col = "green",
         xlim = x_range, ylim = y_range)

    # Add sample names as labels for the first dataset
    text(pcoa.bc1$points[, 1], pcoa.bc1$points[, 2],
         labels = rownames(pcoa.bc1$points),
         pos = 4, cex = cex)

    # Add points for the second dataset
    points(pcoa.bc2$points[, 1], pcoa.bc2$points[, 2], 
           pch = 3, col = "blue")

    # Add sample names as labels for the second dataset
    text(pcoa.bc2$points[, 1], pcoa.bc2$points[, 2], 
         #labels = 1:dim(pcoa.counts2)[1],  # Adjust labels to avoid overlap
         labels = rownames(pcoa.bc2$points),
         pos = 4, cex = cex)
}

##################################
lambda <- 1
mu <- 1
rho <- 0.001
age <- 10
is_one <- F
is_sim <- F
is_inter <- F

transform <- "garland"


##################################
spec = matrix(c(
    'tree', 't', 2, "character",
    'abundance', 'a', 2, "character",
    'metadata', 'm', 2, "character",
    'one', 'o', 0, "logical",
    'all', 'A', 0, "logical",
    'transform', 'T', 2, "character",
    'sim', 's', 0, "logical",
    'inter', 'i', 0, "logical",
    'help' , 'h', 0, "logical"
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


##################################
# reorder
abundance <- abundance[match(colnames(C), rownames(abundance)), ]

# delete all zero columns
abundance <- abundance[, apply(abundance, 2, function(col) any(col != 0))]
abundance <- abundance+1e-5

# Sim
#abundance <- exp(fastBM(tree, sig2 = 1e-3, a = 0, nsim=10))
#print(abundance)
#phylo_sig <- filter_brown(abundance, tree)
#print(phylo_sig$P)
#abundance <- abundance[, which(phylo_sig$P < 0.2), drop=FALSE]
#print(abundance)


##################################
#column_sorted_ind <- order(abs(phylo_sig$K-1), decreasing=F ) # phylo_sig P-value from low to high
#column_sorted_ind <- order( colMeans(abundance, na.rm=T), decreasing=T )
#column_sorted_ind <- order(apply(abundance, 2, max, na.rm = T), decreasing = T)
#abundance = abundance[, column_sorted_ind[1:dim(abundance)[2]]]
#abundance = abundance[, column_sorted_ind[1:5]]

Z <- 1:length(tree$tip.label)
abundance <- abundance[Z,]
C <- C[Z,Z]


Y <- abundance
Q <- abundance

#print(abundance); q()


##################################
#species-site or species-OTU table
if(is_sim){
    do_sim_res <- do_sim()
    abundance <- do_sim_res$abundance
    Q <- do_sim_res$Q
    Y <- do_sim_res$Y
    C <- do_sim_res$C
    tree <- do_sim_res$tree
    above_names <- do_sim_res$above_names; cat("above zero:\t", above_names, "\n")
}
phylo_sig <- filter_brown(abundance, tree)
print(phylo_sig)


###################################################
# reorder
#print(match(colnames(C), rownames(Y))); q()
#Y <- Y[match(colnames(C), rownames(Y)), ]

# use the last col as the ref trait (denominator): dim(Y)[2]
ref_index <- dim(Y)[2]

for(i in 1:dim(Y)[1])
{
    Y[i,]=log((Y[i,]+.01)/(Y[i,ref_index]+.01)) # avoid zero in log-transform
    Q[i,]=Q[i,]/sum(Q[i,])
}
YY=Y[,-ref_index]

#cat("Original divided (log)", "\n")


##################################
# preparation
C_inv=solve(C)


##################################
# calculate the phylo mean (root value) for each node
#ones <- rep(1, length(tree$tip.label))
ones <- rep(1, dim(C)[1])
root_values <- apply(YY, 2, function(col) {
    numerator <- col %*% C_inv %*% ones
    numerator / (t(ones) %*% C_inv %*% ones)
})


##################################
if (grepl("chol", transform, ignore.case = T)){
    #de-correlation by Cholesky decomposition
    L <- solve(t(chol(C)))
    L <- L * sqrt(diag(C)[1])
    X <- L %*% t((t(YY) - root_values))
    X <- t(t(X) + root_values)
} else if(grepl("garland", transform, ignore.case = T)){
    # using C^(-0.5)
    eig_decomp <- eigen(C)
    L <- eig_decomp$vectors %*% diag(1/sqrt(eig_decomp$values)) %*% t(eig_decomp$vectors)
    L <- L * sqrt(diag(C)[1])
    X <- L %*% t((t(YY) - root_values))
    X <- t(t(X) + root_values)
}

rownames(X) <- rownames(C)

cat("\n")
cat("\n")


##################################
#back transformation
P=vector()
for(i in 1:dim(X)[1]){
    R=1/(sum(exp(X[i,]))+1)
    X_added_one <- append(exp(X[i,]), 1, after = ref_index - 1) # add one row at each time
    #print(X_added_one); q()
    #Y_row <- R*c(exp(X[i,]),1) # add one row at each time
    Y_row <- R*X_added_one
    P=rbind(P, Y_row)
}
rownames(P) <- rownames(Q)

rounded_to <- 3

cat("Original compo Q", "\n")
column_sorted_ind <- order( colMeans(Q, na.rm=T), decreasing=T )
#column_sorted_ind <- order(apply(Q, 2, max, na.rm = TRUE), decreasing = TRUE)
range <- 1:10
print(column_sorted_ind[range])
round(Q[, column_sorted_ind[range]], rounded_to)

cat("\n")
cat("Transformed comp P", "\n")
round(P[, column_sorted_ind[range]], rounded_to)

#cat("\n", "YY", "\n")
#round(abundance[, column_sorted_ind[1:10]], rounded_to)

cat("\n\n\n")
pcoa_Q <- calculate_pcoa(Q)

if(is_inter){
    pcoa_P <- calculate_pcoa(Y, "euclidean")
    plot_graph(Q, Y)
} else{
    pcoa_P <- calculate_pcoa(P)
    plot_graph(Q, P, abundance)
}

cat(round(pcoa_Q$eig/sum(pcoa_Q$eig), 3), "\n", round(pcoa_P$eig/sum(pcoa_P$eig), 3), "\n")


q()

##################################
#great, looks similar below:
round(P[1,],5)
round(Q[1,],5)
round(P[1,],5)
round(Q[1,],5)
##################################
#PCoA(P), PCoA(Q)
#
