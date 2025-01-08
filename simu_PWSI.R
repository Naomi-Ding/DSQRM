######### Estimation for Pointwise Single-index
args <- commandArgs(TRUE)
# print(args)
n <- eval(parse(text = args[[1]])) # 100 or 200
p <- eval(parse(text = args[[2]])) # p = 2
m <- eval(parse(text = args[[3]])) # 100 or 200
N <- eval(parse(text = args[[4]])) # 500 or 1000
nsimu <- eval(parse(text = args[[5]]))
print(paste0("n=", n, ", p=", p, ", m=", m, ", N=", N, ", nsimu=", nsimu))

if (!require(np)) {
  install.packages("np")
}
library("np")
if (!require(R.matlab)) {
  install.packages("R.matlab")
}
library("R.matlab")
if (!require(fdadensity)) {
  install.packages("fdadensity")
}
library("fdadensity")
if (!require(pracma)) {
  install.packages("pracma")
}
library("pracma")

data_file <- sprintf("simu_results/simu_n%d_N%d_m%d_p%d_nsimu%d.mat", n, N, m, p, nsimu)
data <- readMat(data_file)
ISE_f <- data$ISE.f
simu.x <- data$simu.x
simu.ally <- data$simu.ally
tmp_f <- colMeans((data$all.fhat - data$truef)^2, dims = 2)
ISE_f <- c(mean(tmp_f), std(tmp_f))

## Estimation
betaest <- array(dim = c(p, m, nsimu))
gest <- array(dim = c(n, m, nsimu))
gest_inv <- array(dim = c(n, m, nsimu))
start_time <- Sys.time()
for (nn in 1:nsimu) {
  print(nn)
  x <- simu.x[, , nn]
  for (s in 1:m) {
    y <- simu.ally[, s, nn]
    bw <- npindexbw(xdat = x, ydat = y)
    model <- npindex(bws = bw, gradients = TRUE)
    beta <- model$beta
    beta <- beta / norm(beta, "2")
    betaest[, s, nn] <- beta
    gest[, s, nn] <- model$mean
  }
  ## get the inversed link function
  gest_inv[, , nn] <- MakeDENsample(gest[, , nn])$DEN
  print(Sys.time() - start_time)
}

## Calculate mean & std of ISE
beta0 <- data$beta0
tmp <- (betaest - array(data = repmat(beta0, 1, nsimu), dim = c(p, m, nsimu)))^2
tmp <- rowMeans(tmp, dims = 2)
ISE_beta_PWSI <- cbind(rowMeans(tmp), apply(tmp, 1, sd))

simu.g <- data$simu.g
tmp <- colMeans(colMeans((gest - simu.g)^2))
ISE_g_PWSI <- c(mean(tmp), sd(tmp))

g_inv <- data$g.inv
tmp <- colMeans(colMeans((gest_inv - g_inv)^2))
ISE_g_inv_PWSI <- c(mean(tmp), sd(tmp))

T_PWSI <- c(
  n, N, ISE_f, m, "PWSI", ISE_beta_PWSI[1, ], ISE_beta_PWSI[2, ],
  ISE_g_PWSI, ISE_g_inv_PWSI
)
tmp_t <- data.frame(t(T_PWSI))
varnames <- c(
  "n", "N", "ISE_f_mean", "ISE_f_std", "m", "tau", "ISE_beta1_mean",
  "ISE_beta1_std", "ISE_beta2_mean", "ISE_beta2_std", "ISE_g_mean", "ISE_g_std",
  "ISE_Psi^{-1}(g)_mean", "ISE_Psi^{-1}(g)_std"
)
colnames(tmp_t) <- varnames
print(tmp_t)

write.table(tmp_t,
  file = "simu_results/simu_estimation_error.csv",
  append = TRUE, row.names = FALSE, col.names = FALSE,
  sep = ",", qmethod = "double"
)
