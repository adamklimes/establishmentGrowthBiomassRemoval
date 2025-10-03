## ordination
library(vegan)

## data_loading
datAux <- read.csv("data/estab_growth_exp_2023.csv", na.strings = c("NA", "MV", "D"))

# derived variables
datAux$BB_AB <- datAux$BB / datAux$AB
datAux$CN <- datAux$C_roots / datAux$N_roots
datAux$AB_removed <- datAux$AB / datAux$Removed_AB

# data_preparation
dat <- datAux
dat$SLA <- log(datAux$SLA)
dat$LDMC <- log(datAux$LDMC)
dat$SRL <- log(datAux$SRL)
dat$RDMC <- log(datAux$RDMC)
dat$Root_diameter <- log(datAux$Root_diameter)
dat$BB <- log(datAux$BB)
dat$AB <- log(datAux$AB)
dat$BB_AB <- log(datAux$BB_AB)
dat$CN <- log(datAux$CN)
dat$AB_removed <- log(datAux$AB_removed)
dat$Removed_AB <- log(datAux$Removed_AB)
dat$Injury <- (datAux$Injury == "Y") + 0
dat$Flowering <- (datAux$Flowering == "Y") + 0
traits <- c("SLA","LDMC","SRL","RDMC","Root_diameter","BB","AB","BB_AB","N_roots","C_roots","N_leaves","CN", "Flowering")
dat$Age_days[dat$Injury == "1"] <- 116 # coding injured plants as a different group from uninjured

# ordination
datspF <- apply(dat[, traits[!(traits %in% c("BB_AB", "CN"))]], 2, tapply, list(dat$Spec, dat$Age_days), mean, na.rm = TRUE, simplify = FALSE)
datsp <- lapply(datspF, function(x) x[, -1]) # 16 NA values in C_roots (14 in day 30)
inputNA <- function(x) {
  x[is.na(x)] <- x[which(is.na(x)) + 1]
  x[is.na(x)] <- x[which(is.na(x)) + 1]
  x
}
datspI <- lapply(datsp, function(x) (apply(x, 1, inputNA)))
datspIst <- lapply(datspI, function(x) (x - mean(x)) / sd(x))
datspA <- do.call(cbind, lapply(datspIst, as.vector))
rownames(datspA) <- paste0(rep(colnames(datspIst$SLA), each = 9), rownames(datspIst$SLA))
clon <- (tapply(dat$Clonality, dat$Spec, unique) == "clon") + 0
plotChull <- function(x, col = "black", lty = 1) {
  polygon(x[chull(x), 1], x[chull(x), 2], border = col, col = rgb(t(col2rgb(col)/255), alpha = 0.2), lty = lty, lwd = 2)
}

pca <- rda(datspA)

cols <- c("black", "red")
plotPCA <- function(pca, mai = c(0.9,0.8,0.1,0.1), main = "", pos = c(3,2,4,4)){
  datg <- apply(scores(pca)$sites, 2, tapply, rep(clon, each = 9) + 1:9*2, mean)
  expVar <- round(pca$CA$eig / sum(pca$CA$eig), 3) * 100

  par(mai = mai)
  biplot(pca, display = "sp", type = "n", xlab = paste0("PCA1 (", expVar[1], "%)"), ylab = paste0("PCA2 (", expVar[2], "%)"), main = main)
  arrows(0, 0, x1 = scores(pca)$species[, 1], y1 = scores(pca)$species[, 2], length = 0.1, col = "blue")
  plotChull(scores(pca)$sites[1:20*9-8, ][clon == 0, ], col = cols[1])
  plotChull(scores(pca)$sites[1:20*9-8, ][clon == 1, ], col = cols[2])
  plotChull(scores(pca)$sites[1:20*9-1, ][clon == 0, ], col = cols[1])
  plotChull(scores(pca)$sites[1:20*9-1, ][clon == 1, ], col = cols[2])
  plotChull(scores(pca)$sites[1:20*9, ][clon == 0, ], lty = 2, col = cols[1])
  plotChull(scores(pca)$sites[1:20*9, ][clon == 1, ], col = cols[2], lty = 2)

  lines(datg[1:8*2, ], col = cols[2], lwd = 3)
  points(datg[16, 1], datg[16, 2], col = cols[2], pch = 16, cex = 1)
  text(datg[2, 1], datg[2, 2], "30 days", pos = pos[2], cex = 0.8, col = cols[2])
  text(datg[16, 1], datg[16, 2], "115 days", pos = pos[4], cex = 0.8, col = cols[2])
  lines(datg[1:8*2-1, ], lwd = 3, col = cols[1])
  points(datg[15, 1], datg[15, 2], pch = 16, cex = 1, col = cols[1])
  text(datg[1, 1], datg[1, 2], "30 days", pos = pos[1], cex = 0.8, col = cols[1])
  text(datg[15, 1], datg[15, 2], "115 days", pos = pos[3], cex = 0.8, col = cols[1])
  points(datg[18, 1], datg[18, 2], col = cols[2], pch = 8, cex = 1.3)
  text(datg[18, 1], datg[18, 2], "Removal", pos = 2, cex = 0.8, col = cols[2])
  points(datg[17, 1], datg[17, 2], col = cols[1], pch = 8, cex = 1.3)
  text(datg[17, 1], datg[17, 2], "Removal", pos = 2, cex = 0.8, col = cols[1])
  legend("topleft", lwd = 3, col = cols, legend = c("Non-clonal", "Clonal"), bty = "n")
}


# fig 1
# png("figures/Fig1_ordination.png", width = 480*13, height = 480*9, res = 720)
plotPCA(pca)
text(scores(pca)$species[, 1], scores(pca)$species[, 2], c("SLA", "LDMC", "SRL", "RDMC", "Root diameter", "Belowground biomass", "Aboveground biomass", "N roots", "C roots", "N leaves", "Flowering"), col = "blue", cex = 0.8, pos = c(1,3,2,4,4,3,1,1,4,2,4))

# PCA without inputed values
pca2 <- rda(datspA[, colnames(datspA) != "C_roots"])

# png("figures/FigS1_compareOrd.png", width = 480*28, height = 480*10, res = 72*13)
par(mfrow = c(1,2))
plotPCA(pca, mai = c(0.9,0.8,0.3,0.1), main = "All traits")
text(scores(pca)$species[, 1], scores(pca)$species[, 2], c("SLA", "LDMC", "SRL", "RDMC", "Root diameter", "Belowground biomass", "Aboveground biomass", "N roots", "C roots", "N leaves", "Flowering"), col = "blue", cex = 0.8, pos = c(1,3,2,4,4,3,1,1,4,2,4))
plotPCA(pca2, mai = c(0.9,0.8,0.3,0.1), main = "Without C roots", pos = c(2,3,1,4))
text(scores(pca2)$species[, 1], scores(pca2)$species[, 2], c("SLA", "LDMC", "SRL", "RDMC", "Root diameter", "Belowground biomass", "Aboveground biomass", "N roots", "N leaves", "Flowering"), col = "blue", cex = 0.8, pos = c(1,1,2,4,4,3,1,2,2,4))
#_
