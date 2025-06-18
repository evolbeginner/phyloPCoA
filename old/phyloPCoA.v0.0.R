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
    tree <- read.tree(opt$tree)
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
abundance <- abundance[, colnames(abundance) %in% overlapping_items]

C=vcv(tree)


##################################
#species-site or species-OTU table
snum=6
hnum=500


##################################
#hosts species phylogeny
# using ape rtree()
#tre=rtree(snum)
# using TreeSim
trees <- sim.bd.taxa.age(snum, 1, lambda=1, mu=0, age=1); tre <- trees[[1]]
write.tree(tre, file='sim.tree')

# sim trait of mu = 10, and sig2=1
ssmat <- exp(fastBM(tre, sig2 = 1, a = 1, nsim=hnum))
if(F){ # diff root_values
    mean_a=1; sd_a=1; nsim=hnum
    a_values <- rnorm(nsim, mean = mean_a, sd = sd_a)
    ssmat <- matrix(0, nrow = snum, ncol = nsim); 
    for (i in 1:nsim) {
        ssmat[,i] <- exp(fastBM(tre, sig2 = 1, a = a_values[i], nsim=1))
    }
}

# generate the phylo covariance matrix, C
C=vcv(tre)

# conversion
Y=ssmat
Q=ssmat

# use the last col as the ref trait (denominator): dim(Y)[2]
ref_index <- dim(Y)[2]
#ref_index <- 3


for(i in 1:dim(Y)[1])
{
    Y[i,]=log((Y[i,]+.01)/(Y[i,ref_index]+.01)) # avoid zero in log-transform
    Q[i,]=Q[i,]/sum(Q[i,])
}
YY=Y[,-ref_index]


cat("Original divided (log)", "\n")
#print(round(YY,3))


##################################
# preparation
C_inv=solve(C)


##################################
# calculate the phylo mean (root value) for each node
ones <- rep(1, snum)
root_values <- apply(YY, 2, function(col) {
    numerator <- col %*% C_inv %*% ones
    numerator / (t(ones) %*% C_inv %*% ones)
})


##################################
#de-correlation by Cholesky decomposition
L <- solve(t(chol(C))) #Shen
X <- L %*% t((t(YY) - root_values))
X <- t(t(X) + root_values)


####
####
# using C^(-0.5)
if(T){
    eig_decomp <- eigen(C)
    L <- eig_decomp$vectors %*% diag(1/sqrt(eig_decomp$values)) %*% t(eig_decomp$vectors)
    X <- L %*% t((t(YY) - root_values))
    X <- t(t(X) + root_values)
}


cat("X (transformed)", "\n")
print(round(X,3))
cat("\n")

print(L)
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

#print(Q);

cat("Original compo Q", "\n")
head(round(Q,5)[,1:20])

cat("Transformed comp P", "\n")
head(round(P,5)[,1:20])

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
