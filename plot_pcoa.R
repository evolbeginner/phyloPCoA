pcoa.counts1 <- P  # Replace with your first data matrix
pcoa.counts2 <- Q  # Replace with your second data matrix

# Perform PCoA using Bray-Curtis distance for the first dataset
pcoa.bc1 <- pco(dsvdis(pcoa.counts1, index = "bray/curtis"))

# Perform PCoA using Bray-Curtis distance for the second dataset
pcoa.bc2 <- pco(dsvdis(pcoa.counts2, index = "bray/curtis"))

# Create the initial plot for the first dataset
plot(pcoa.bc1$points[, 1], pcoa.bc1$points[, 2],
     pch = 3, xlab = "PC1", ylab = "PC2", main = "PCoA: Bray-Curtis", col = "green")

# Add sample names as labels for the first dataset
text(pcoa.bc1$points[, 1], pcoa.bc1$points[, 2], 
     labels = 1:dim(pcoa.counts1)[1], 
     pos = 4, cex = 0.8)

# Add points for the second dataset
points(pcoa.bc2$points[, 1], pcoa.bc2$points[, 2], 
       pch = 3, col = "blue")

# Add sample names as labels for the second dataset
text(pcoa.bc2$points[, 1], pcoa.bc2$points[, 2], 
     labels = 1:dim(pcoa.counts2)[1],  # Adjust labels to avoid overlap
     pos = 4, cex = 0.8)
