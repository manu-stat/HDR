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
if(!exists("N")) N = 1000

# 1. Normal distribution: elliptical HDR
S1 = rmvnorm(N, mean = c(0, 0, 0), sigma = diag(3))
f1 = dmvnorm(S1, mean = c(0, 0, 0), sigma = diag(3))
margin1 = c("norm", "norm", "norm")


# 2. Nonlinear dependence: Gumbel copula
cop = gumbelCopula(param = 1.8, dim = 3)
U = rCopula(N, cop)
S2 = cbind(qnorm(U[, 1]), qnorm(U[, 2]), qexp(U[, 3], rate = .8))
f2 = dCopula(U, cop) * dnorm(S2[,1]) * dnorm(S2[,2]) * dexp(S2[, 3], rate = .8)
margin2 = c("norm", "norm", "exp")


# 3. Nonlinear dependence: t copula
cop = tCopula(param = 0.3, dim = 3, df = 1)
U = rCopula(N, cop)
S3 = cbind(qnorm(U[, 1]), qnorm(U[, 2]), qexp(U[,3], 1))
f3 = dCopula(U, cop) * dnorm(S3[,1]) * dnorm(S3[,2]) * dexp(S3[,3], 1)
margin3 = c("norm", "norm", "exp")


# 4. Exponential dependence: t copula
cop = gumbelCopula(param = 1.5, dim = 3)
U = rCopula(N, cop)
S4 = cbind(qexp(U[, 1], rate = 1), qexp(U[, 2], rate = 0.5), qexp(U[, 3], rate = 0.2))
f4 = dCopula(U, cop) * dexp(S4[,1], rate = 1) * dexp(S4[,2], rate = 0.5) * dexp(S4[,3], rate = 0.2)
margin4 = c("exp", "exp", "exp")


scenarios = list(
  Normal = list(S1, f1, margin1),
  Gumbel = list(S2, f2, margin2),
  Tail = list(S3, f3, margin3),
  Skewed = list(S4, f4, margin4)
)



# # Visualize
# invisible(lapply(names(scenarios), function(nm) {
#   fig = plotly::plot_ly(data = data.frame(scenarios[[nm]][[1]]), x = ~X1, y = ~X2, z = ~X3, type = "scatter3d", mode = "markers", marker = list(size = 4, color = "#42314A"))
#   print(fig)
#   readline("Next: ENTER")
# }))