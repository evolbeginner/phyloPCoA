#! /bin/env Rscript

#Sishuo's idea on accounting phylogeny for proportion data

###############################################################
# Update history
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
library(phytools)
library(vegan)


##################################
calculate_pcoa <- function(X) {
    bc_distance <- vegdist(X, method = "bray")
    pcoa_results <- cmdscale(bc_distance, eig = T, k = 2)
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

    # Reshape the average_abundance data frame to wide format
    average_abundance_wide <- average_abundance %>%
      pivot_wider(names_from = species, values_from = average_abundance)

    average_abundance_wide = as.data.frame(average_abundance_wide)
    rownames(average_abundance_wide) <- average_abundance_wide[,1]
    average_abundance_wide <- average_abundance_wide[,-1]

    return(average_abundance_wide)
}


##################################
lambda <- 1
mu <- 1
rho <- 1
age <- 10
is_one <- T


##################################
spec = matrix(c(
    'tree', 't', 2, "character",
    'abundance', 'a', 2, "character",
    'metadata', 'm', 2, "character",
    'one', 'o', 0, "logical",
    'all', 'A', 0, "logical",
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

if (!is.null(opt$help)) {
    cat(getopt(spec, usage = T))
    q(status = 1)
}

if(! is.null(opt$tree)){
    tre <- read.tree(opt$tree)
    C <- vcv(tre)
    rownames(C) <- colnames(C) <- colnames(vcv(tre))
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

column_sorted_ind <- order( colMeans(abundance, na.rm=T), decreasing=T )
column_sorted_ind <- order(apply(abundance, 2, max, na.rm = T), decreasing = T)
#abundance = abundance[, -column_sorted_ind[1:500]]

# delete all zero columns
abundance <- abundance[, apply(abundance, 2, function(col) any(col != 0))]

# Sim
# abundance <- exp(fastBM(tre, sig2 = 0.01, a = 0, nsim=dim(abundance)[2]))

# reorder
abundance <- abundance[match(colnames(C), rownames(abundance)), ]
Y <- abundance
Q <- abundance
#head(abundance); q()


##################################
#species-site or species-OTU table
if(F){

snum=16
hnum=1000

##################################
#hosts species phylogeny
# using ape rtree()
#tre=rtree(snum)
# using TreeSim

# Do simulation
#trees <- sim.bd.taxa.age(snum, 1, lambda=1, mu=0, rho=1, age=1); tre <- trees[[1]]
trees <- sim.bd.taxa.age(snum, 1, lambda, mu, rho, age); tre <- trees[[1]]
write.tree(tre, file='sim.tree')

# sim trait of mu = 1, and sig2=1
#abundance <- exp(fastBM(tre, sig2 = 1, a = 1, nsim=hnum))
# diff root_values
mean_a=0; sd_a=5; nsim=hnum
a_values <- rnorm(nsim, mean = mean_a, sd = sd_a)
abundance <- matrix(0, nrow = snum, ncol = nsim); 
for (i in 1:nsim) {
    abundance[,i] <- exp(fastBM(tre, sig2 = 1, a = a_values[i], nsim=1))
}
#for(i in 1:10){ v <- rnorm(snum, 20, 2); abundance[,i] <- v }

# generate the phylo covariance matrix, C
C=vcv(tre)

# conversion
Y=abundance
Q=abundance

}


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

cat("Original divided (log)", "\n")


##################################
# preparation
C_inv=solve(C)


##################################
# calculate the phylo mean (root value) for each node
ones <- rep(1, length(tre$tip.label))
root_values <- apply(YY, 2, function(col) {
    numerator <- col %*% C_inv %*% ones
    numerator / (t(ones) %*% C_inv %*% ones)
})


##################################
#de-correlation by Cholesky decomposition
L <- solve(t(chol(C)))
L <- L * sqrt(diag(C)[1])
X <- L %*% t((t(YY) - root_values))
X <- t(t(X) + root_values)


####
####
# using C^(-0.5)
if(F){
    eig_decomp <- eigen(C)
    L <- eig_decomp$vectors %*% diag(1/sqrt(eig_decomp$values)) %*% t(eig_decomp$vectors)
    L <- L * sqrt(diag(C)[1])
    X <- L %*% t((t(YY) - root_values))
    X <- t(t(X) + root_values)
}

rownames(X) <- rownames(C)


cat("\n")

#print(L)
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

#print(Q);

rounded_to <- 6

cat("Original compo Q", "\n")
column_sorted_ind <- order( colMeans(Q, na.rm=T), decreasing=T )
#column_sorted_ind <- order(apply(Q, 2, max, na.rm = TRUE), decreasing = TRUE)
print(column_sorted_ind[1:30])
round(Q[, column_sorted_ind[1:10]], rounded_to)

cat("Transformed comp P", "\n")
round(P[, column_sorted_ind[1:10]], rounded_to)

cat("YY", "\n")
round(abundance[, column_sorted_ind[1:10]], rounded_to)

cat("\n\n\n")
pcoa_Q <- calculate_pcoa(Q)
pcoa_P <- calculate_pcoa(P)
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
