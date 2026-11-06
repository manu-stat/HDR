#| Scenario	Proprietà da verificare	Forma attesa dell'HDR
#| 1. Normale	Caso di riferimento	Ellisse circolare
#| 2. Normale correlata	Invarianza rispetto alla correlazione	Ellisse inclinata
#| 3. Miscela bimodale	Capacità di rilevare più cluster	Due regioni distinte
#| 4. Student-t	Robustezza alle code pesanti	Regione centrale più compatta rispetto alle code
#| 5. Esponenziali	Gestione dell'asimmetria	Regione asimmetrica
#| 6. Copula Gumbel	Dipendenza non lineare e di coda	Regione non ellittica



library(mvtnorm)
library(copula)

if(exists("seed")) set.seed(seed) else set.seed(123)
if(!exists("N")) N = 500

# 1. Normal distribution: elliptical HDR
S1 = rmvnorm(N, mean = c(0, 0), sigma = diag(2))
f1 = dmvnorm(S1, mean = c(0, 0), sigma = diag(2))
margin1 = "norm"

# 2. Nonlinear dependence: Gumbel copula
cop = gumbelCopula(param = 1.8, dim = 2)
U = rCopula(N, cop)
S2 = cbind(qnorm(U[, 1]), qnorm(U[, 2]))
f2 = dCopula(U, cop) * dnorm(S2[,1]) * dnorm(S2[,2])
margin2 = "norm"


# 3. Nonlinear dependence: t copula
cop = tCopula(param = 0.5, dim = 2, df = 2)
U = rCopula(N, cop)
S3 = cbind(qnorm(U[, 1]), qnorm(U[, 2]))
f3 = dCopula(U, cop) * dnorm(S3[,1]) * dnorm(S3[,2])
margin3 = "norm"


# 4. Exponential dependence: Gumbel copula
cop = gumbelCopula(param = 1.5, dim = 2)
U = rCopula(N, cop)
S4 = cbind(qexp(U[, 1], rate = 1), qexp(U[, 2], rate = 0.5))
f4 = dCopula(U, cop) * dexp(S4[,1], rate = 1) * dexp(S4[,2], rate = 0.5)
margin4 = "exp"

scenarios = list(
  Normal = list(S1, f1, margin1),
  Gumbel = list(S2, f2, margin2),
  Tail = list(S3, f3, margin3),
  Skewed = list(S4, f4, margin4)
)



# # Visualize 
# par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
# invisible(lapply(names(scenarios), function(nm) plot(scenarios[[nm]][[1]], pch = 16, cex = 0.45, col = "#42314A", xlab = "X1", ylab = "X2", main = nm)))
# par(mfrow = c(1, 1))


