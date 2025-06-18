#! /bin/env Rscript

#Sishuo's idea on accounting phylogeny for proportion data

###############################################################
# Update history
# 2024-10-25
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
library(TreeSim)
library(phytools)


##################################
spec = matrix(c(
    'tree', 't', 2, "character",
    'abundance', 'a', 2, "character",
    'metadata', 'm', 2, "character",
    'help' , 'h', 0, "logical"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

if (!is.null(opt$help)) {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

if(! is.null(opt$tree)){
    tre <- read.tree(opt$tree)
    # get phylo-Covariance matrix
    C=vcv(tre)
    C <- diag(diag(C))
    rownames(C) <- colnames(C) <- colnames(vcv(tre))
}

if(! is.null(opt$abundance)){
    abundance <- read.table(opt$abundance, header=T)
}

if (!is.null(opt$metadata)) {
    metadata <- read.table(opt$metadata, header = TRUE)

    # Randomly select one sample_id for each species
    selected_samples <- metadata %>%
        group_by(species) %>%
        sample_n(1)  # Sample one row from each group
}
selected_samples <- as.data.frame(selected_samples)

overlapping_items <- intersect(colnames(abundance), selected_samples$sample_id)

# for only overlapped sample_id
new_df <- abundance[, colnames(abundance) %in% overlapping_items]
name_mapping <- setNames(selected_samples$species, selected_samples$sample_id)
#print(name_mapping)
colnames(new_df) <- name_mapping[colnames(new_df)]
abundance <- t(new_df)

column_means_sorted_ind <- order( colMeans(abundance, na.rm=T), decreasing=T )
abundance = abundance[, column_means_sorted_ind[1:10]]

# reorder
abundance <- abundance[match(colnames(C), rownames(abundance)), ]
Y <- abundance
Q <- abundance


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
trees <- sim.bd.taxa.age(snum, 1, lambda=1, mu=0, age=1); tre <- trees[[1]]
write.tree(tre, file='sim.tree')

# sim trait of mu = 10, and sig2=1
abundance <- exp(fastBM(tre, sig2 = 10, a = 1, nsim=hnum))
# diff root_values
mean_a=0; sd_a=0.1; nsim=hnum
a_values <- rnorm(nsim, mean = mean_a, sd = sd_a)
abundance <- matrix(0, nrow = snum, ncol = nsim); 
for (i in 1:nsim) {
    abundance[,i] <- exp(fastBM(tre, sig2 = 1, a = a_values[i], nsim=1))
}
for(i in 1:10){
    v <- rep(0,snum); v <- append(v, rnorm(1, 900, 30), after = runif(1,1,snum-1))
    v <- rnorm(snum, 90, 30)
    abundance[,i] <- v
}

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
L <- solve(t(chol(C))) #Shen
#diag(L) <- 1
#print(L); q()
X <- L %*% t((t(YY) - root_values))
X <- t(t(X) + root_values)
#head(YY)
#X <- L %*% YY


####
####
# using C^(-0.5)
if(T){
    eig_decomp <- eigen(C)
    L <- eig_decomp$vectors %*% diag(1/sqrt(eig_decomp$values)) %*% t(eig_decomp$vectors)
    X <- L %*% t((t(YY) - root_values))
    X <- t(t(X) + root_values)
}

rownames(X) <- rownames(C)
head(X)#; q()

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

cat("Original compo Q", "\n")
column_sorted_ind <- order( colMeans(Q, na.rm=T), decreasing=T )
column_sorted_ind <- order(apply(Q, 2, max, na.rm = TRUE), decreasing = TRUE)
print(column_sorted_ind[1:30])
round(Q[, column_sorted_ind[1:20]], 2)

cat("Transformed comp P", "\n")
round(P[, column_sorted_ind[1:20]], 2)

cat("YY", "\n")
round(abundance[, column_sorted_ind[1:20]], 2)


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
