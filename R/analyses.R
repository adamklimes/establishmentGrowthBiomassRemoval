# analyses
library(rstan)
rstan_options(auto_write = TRUE)

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
dat$Flowering <- (datAux$Flowering == "Y") + 0
dat$Injury <- (datAux$Injury == "Y") + 0
st <- function(x) (x - mean(x)) / sd(x)
traits <- c("SLA","LDMC","SRL","RDMC","Root_diameter","BB","AB","BB_AB","N_roots","C_roots","N_leaves","CN", "Flowering")
dat$Age_days[dat$Injury == "1"] <- 116 # coding injured plants as a different group from uninjured

prepDatStan <- function(trait, dat){
  sel <- !is.na(dat[, trait])
  timeCode <- (1:10)[factor(dat$Age_days)][sel]
  resp <- dat[sel, trait]
  datStan <- list(
    N = nrow(dat[sel, ]),
    NspecTime = length(unique(dat$Spec)) * length(unique(dat$Age_days)),
    NclonTime = length(unique(dat$Clonality)) * length(unique(dat$Age_days)),
    resp = if (all(dat[sel, trait] %in% 0:1)) resp else st(resp),
    clonTime = timeCode + (dat$Clonality == "clon")[sel] * 10,
    specTime = timeCode + (as.numeric(factor(dat$Spec))[sel] - 1) * 10
  )
  list(dat = datStan, traitMean = mean(resp), traitSD = sd(resp))
}

# analyses
runStan <- function(trait, dat, seed = 121) {
  datStan <- prepDatStan(trait, dat)
  fileStan <- if (all(datStan$dat$resp %in% 0:1)) "R/mod_bernoulli.stan" else "R/mod_univariate.stan"
  mod <- stan(file = fileStan, data = datStan$dat,
              iter = 8000, seed = seed, pars = c("randSpecTimeStd", "randSpecTime"), include = FALSE)
  list(mod = mod, datStan = datStan)
}
mods <- lapply(traits, runStan, dat)
names(mods) <- traits
# save(mods, file = "data/analyses/mods.RData")

# results
load(file = "data/analyses/mods.RData")
sumFn <- function(x) c(mean = mean(x), quantile(x, c(0.025,0.25,0.75,0.975)))
age <- sort(unique(dat$Age_days))
age <- age[age != 116]
addDiffMark <- function(diffB, ageM, col, yshift = 0){
  if (sign(diffB[2]) == sign(diffB[5])){
    ps <- mean(c(ageM, age[which(age == ageM) - 1]))
    xr <- diff(par("usr")[1:2])
    yr <- diff(par("usr")[3:4])
    if (sign(diffB[2]) > 0) arrows(ps - xr/50, par("usr")[4] - yr/7 + yshift * yr, x1 = ps + xr/50, y1 = par("usr")[4] - yr/20 + yshift * yr, length = 0.1, col = col, lwd = 2)
    if (sign(diffB[2]) < 0) arrows(ps - xr/50, par("usr")[4] - yr/20 + yshift * yr, x1 = ps + xr/50, y1 = par("usr")[4] - yr/7 + yshift * yr, length = 0.1, col = col, lwd = 2)
  }
}
invLogit <- function(x) exp(x) / (1 + exp(x))
plotTrait <- function(trait, mods, cols, ylab, lab, fn = log, axisX = FALSE){
  parB <- extract(mods[[trait]]$mod, pars = "B")$B * mods[[trait]]$datStan$traitSD + mods[[trait]]$datStan$traitMean
  if (all(mods[[trait]]$datStan$dat$resp %in% 0:1)) parB <- invLogit(extract(mods[[trait]]$mod, pars = "B")$B)
  parBsum <- data.frame(apply(parB, 2, sumFn))
  diffB <- data.frame(apply(t(apply(parB, 1, diff)), 2, sumFn))

  plot(range(age), range(parBsum), pch = 16, type = "n", xlab = "",
    ylab = "", axes = FALSE)
  box(bty = "l")
  if (axisX) {
    axis(1, labels = "Age (days)", at = mean(par("usr")[1:2]), tick = FALSE, line = 1.5)
    axis(1, labels = age, at = age)
  }
  axis(2, labels = ylab, at = mean(par("usr")[3:4]), tick = FALSE, line = 2.5)
  axis(2, labels = lab, at = fn(lab), las = 2)
  lines(age-0.5, parBsum["mean", 1:9], col = cols[1], lwd = 2)
  lines(age+0.5, parBsum["mean", 11:19], col = cols[2], lwd = 2)
  points(age-0.5, parBsum["mean", 1:9], pch = 16, col = cols[1])
  invisible(Map(function(y, x) lines(rep(x, 2), y[c(2, 5)], col = cols[1]), parBsum[, 1:9], age-0.5))
  invisible(Map(function(y, x) lines(rep(x, 2), y[c(3, 4)], lwd = 3, col = cols[1]), parBsum[, 1:9], age-0.5))
  points(age+0.5, parBsum["mean", 11:19], pch = 16, col = cols[2])
  invisible(Map(function(y, x) lines(rep(x, 2), y[c(2, 5)], col = cols[2]), parBsum[, 11:19], age+0.5))
  invisible(Map(function(y, x) lines(rep(x, 2), y[c(3, 4)], lwd = 3, col = cols[2]), parBsum[, 11:19], age+0.5))
#  invisible(Map(addDiffMark, diffB[, 1:8], age[2:9], col = cols[1]))
#  invisible(Map(addDiffMark, diffB[, 10:17], age[2:9], col = cols[2], yshift = -1/20))
  points(115-0.5, parBsum["mean", 10], pch = 8, col = cols[1], cex = 1.5)
  lines(c(0, 114.5), rep(parBsum["mean", 10], 2), col = cols[1], lty = 2)
  points(115+0.5, parBsum["mean", 20], pch = 8, col = cols[2], cex = 1.5)
  lines(c(0, 115.5), rep(parBsum["mean", 20], 2), col = cols[2], lty = 2)
  lines(c(114.5,114.5), parBsum[c(2, 5), 10], col = cols[1])
  lines(c(115.5,115.5), parBsum[c(2, 5), 20], col = cols[2])
  lines(c(114.5,114.5), parBsum[3:4, 10], col = cols[1], lwd = 3)
  lines(c(115.5,115.5), parBsum[3:4, 20], col = cols[2], lwd = 3)
}
cols <- c("black", "red")
#cols <- c("orange", "purple")

# png("figures/Fig2_traitsPaired.png", width = 480*13, height = 480*8, res = 72*7.5)
par(mfrow = c(4,4), mai = c(0.5,0.6,0.1,0.1))
plotTrait("AB", mods, cols, "Aboveground biomass (g)", c(0.001,0.01,0.1,1,10))
plotTrait("SLA", mods, cols, bquote(SLA*" "*(cm^2/g)), c(250,300,400,500,700))

plotTrait("LDMC", mods, cols, "LDMC", c(0.1,0.15,0.2))
plotTrait("N_leaves", mods, cols, "N in leaves (%)", fn = function(x) x, 4:6)

plotTrait("BB", mods, cols, "Belowground biomass (g)", c(0.001,0.01,0.1,1))
plotTrait("SRL", mods, cols, "SRL (m/g)", c(200,300,400,500,600))
plotTrait("RDMC", mods, cols, "RDMC", c(0.05,0.07,0.1))
plotTrait("N_roots", mods, cols, "N in roots (%)", fn = function(x) x, 2:5, axisX = TRUE)

plotTrait("C_roots", mods, cols, "C in roots (%)", fn = function(x) x, c(35,40,45,50))
plotTrait("CN", mods, cols, "C:N ratio in roots", c(6,10,15,25))
plotTrait("Root_diameter", mods, cols, "Root diameter (mm)", c(0.2,0.22,0.25,0.28), axisX = TRUE)
plot(1,1, type = "n", axes = FALSE, ann = FALSE)

plotTrait("BB_AB", mods, cols, "Root:shoot ratio", c(0.1,0.2,0.5,0.8), axisX = TRUE)
plotTrait("Flowering", mods, cols, "Flowering (%)", fn = function(x) x/100, c(0,20,50,80), axisX = TRUE)
plot(1,1, type = "n", axes = FALSE, ann = FALSE)
legend("bottomleft", pch = c(16, 16, NA, 8, 8), lwd = c(2, 2, NA, NA, NA), col = c(cols, NA, cols), legend = c("Non-clonal", "Clonal", NA, "Removal treatment non-clonal", "Removal treatment clonal"), bty = "n", xpd = NA, cex = 1.5)
plot(1,1, type = "n", axes = FALSE, ann = FALSE)


 # overall comparision clon vs. non-clon
compClon <- function(mod, tr){
  parB <- extract(mod$mod, pars = "B")$B #* mod$datStan$traitSD + mod$datStan$traitMean
  #if (all(mod$datStan$dat$resp %in% 0:1)) parB <- invLogit(extract(mod$mod, pars = "B")$B)
  clonEff <- rowMeans(parB[, 11:19] - parB[, 1:9])
  #if (tr) clonEff <- exp(clonEff)
  sumFn(clonEff)
}
clonEff <- do.call(rbind, Map(compClon, mods, c(1,1,1,1,1,1,1,1,0,0,0,1,0)))
do.call(rbind, apply(round(clonEff, 2), 1, function(x) paste0(x[1]," [", x[2],",", x[5], "]"), simplify = FALSE))[c(7,6,1,3,2,4,11,9,10,12,5,8,13), , drop = FALSE]

# Table 2
calcEffs <- function(mod){
  parB <- extract(mod$mod, pars = "B")$B #* mod$datStan$traitSD + mod$datStan$traitMean
  #if (all(mod$datStan$dat$resp %in% 0:1)) parB <- invLogit(extract(mod$mod, pars = "B")$B)
  effs <- round(rbind(effClon = sumFn(parB[, 19] - parB[, 9]),
    effRem = sumFn((parB[, 10] + parB[, 20] - parB[, 9] - parB[, 19]) / 2),
    effClonRem = sumFn(parB[, 20] - parB[, 19] - parB[, 10] + parB[, 9])), 2)
  setNames(paste0(effs[, "mean"], " [", effs[, "2.5%"], ", ", effs[, "97.5%"], "]"), c("effClon", "effRem", "effClonRem"))
}
do.call(rbind, lapply(mods, calcEffs))[c(7,6,1,3,2,4,11,9,10,12,5,8,13), ]

# Removed_AB and AB_removed
prepDatStanSupp <- function(trait, dat){
  sel <- !is.na(dat[, trait])
  resp <- dat[sel, trait]
  datStan <- list(
    N = nrow(dat[sel, ]),
    Nspec = length(unique(dat$Spec)) * length(unique(dat$Age_days)),
    resp = st(resp),
    clon = (dat$Clonality == "clon")[sel] + 1,
    spec = as.numeric(factor(dat$Spec))[sel]
  )
  list(dat = datStan, traitMean = mean(resp), traitSD = sd(resp))
}
runStanSupp <- function(trait, dat, seed = 121) {
  datStan <- prepDatStanSupp(trait, dat)
  fileStan <- "R/mod_removedB.stan"
  mod <- stan(file = fileStan, data = datStan$dat,
              iter = 8000, seed = seed, pars = c("randSpecStd", "randSpec"), include = FALSE)
  list(mod = mod, datStan = datStan)
}
modsSupp <- lapply(c("Removed_AB", "AB_removed"), runStanSupp, dat)
names(modsSupp) <- c("Removed_AB", "AB_removed")
# save(modsSupp, file = "data/analyses/modsSupp.RData")

compClonSupp <- function(mod){
  parB <- extract(mod$mod, pars = "B")$B
  sumFn(parB[, 2] - parB[, 1])
}
effClonSupp <- round(do.call(rbind, lapply(modsSupp, compClonSupp)), 2)
paste0(effClonSupp[, "mean"], " [", effClonSupp[, "2.5%"], ", ", effClonSupp[, "97.5%"], "]")
#_
