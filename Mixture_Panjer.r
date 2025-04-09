# ===================================================================
# Fitting Distribution Models: Binomial, Poisson, NB, ZMNB,
# Polya-Aeppli, Panjer, and Mixtures (Poisson, NB, Panjer)
# ===================================================================

# ----------------------
# Load required libraries
# ----------------------
library(reshape2)
library(ggplot2)
library(gamlss)
library(gamlss.dist)
library(polyaAeppli)

# ----------------------
# Define observed data
# ----------------------
C0 <- c(122628, 21686, 4014, 832, 224, 68, 17, 7, 7)  # Frequency of claims
k <- 0:(length(C0) - 1)                               # Number of claims
limitek <- length(k)                                  # Max index for k
Prob <- C0 / sum(C0)                                  # Empirical probability
Esp <- sum(k * Prob)                                  # Empirical mean
Var <- sum(Prob * (k - Esp)^2)                        # Empirical variance
Indice <- Var / Esp                                   # Dispersion index
n_total <- sum(C0)                                    # Total observations

# ----------------------
# Binomial distribution
# ----------------------
pest <- Esp / n_total
binomial <- dbinom(k, 14, pest) * n_total
binomial_rounded <- round(binomial, 2)

# ----------------------
# Poisson distribution
# ----------------------
Poisson <- dpois(k, Esp) * n_total
Poisson_rounded <- round(Poisson, 2)

# ----------------------
# Negative Binomial distribution
# ----------------------
r_est <- Esp^2 / (Var - Esp)
p_est <- Esp / Var
BN <- dnbinom(k, r_est, p_est) * n_total
BN_rounded <- round(BN, 2)

# ----------------------
# Zero-Modified Negative Binomial (ZMNB) distribution
# ----------------------
data_zmnb <- data.frame(k = k, n_k = C0)

# Fit the model using gamlss
zm_model <- gamlss(
  k ~ 1, 
  family = ZANBI(), 
  weights = C0,
  data = data_zmnb,
  control = gamlss.control(n.cyc = 350, c.crit = 0.01, trace = FALSE)
)

# Extract parameters
mu <- exp(coef(zm_model, "mu"))
sigma <- exp(coef(zm_model, "sigma"))
nu <- plogis(coef(zm_model, "nu"))

# Define custom PMF for ZMNB
pmf_zmnb <- function(x) {
  ifelse(x == 0, 
         nu + (1 - nu) * dnbinom(0, mu = mu, size = 1/sigma),
         (1 - nu) * dnbinom(x, mu = mu, size = 1/sigma))
}

ZM_NB_dist <- round(n_total * pmf_zmnb(k), 2)

# ----------------------
# Polya-Aeppli distribution
# ----------------------

# Define negative log-likelihood for grouped data
neg_loglik_grouped <- function(params, k, n_k) {
  lambda <- params[1]
  prob <- params[2]
  if (lambda <= 0 || prob <= 0 || prob >= 1) return(Inf)
  -sum(n_k * dPolyaAeppli(k, lambda = lambda, prob = prob, log = TRUE))
}

# Estimate parameters using optimization
fit <- optim(
  par = c(1, 0.5), 
  fn = neg_loglik_grouped,
  k = k, n_k = C0,
  method = "L-BFGS-B",
  lower = c(1e-8, 1e-8), 
  upper = c(Inf, 1 - 1e-8)
)

lambda_hat <- fit$par[1]
prob_hat <- fit$par[2]
theta_hat <- 1 - prob_hat
pmf_polya <- dPolyaAeppli(k, lambda_hat, prob_hat)
polyaAeppli_dist <- round(pmf_polya * n_total, 2)

# ----------------------
# Panjer distribution
# ----------------------

qm <- 1 - Esp / Var
rmo <- Esp^2 / (Var - Esp)
alpha <- qm
beta <- (rmo - 1) * qm

Panjer <- gamma(rmo + k) / (gamma(rmo) * factorial(k)) * alpha^k * (1 - alpha)^rmo * n_total
Panjer_rounded <- round(Panjer, 2)

# ----------------------
# Mixture Models: Factorial Moments Method
# ----------------------

# Function to compute the first 4 factorial moments
MoFactorial <- function(Prob, limitek) {
  sapply(1:4, function(r) {
    sum(sapply(r:(limitek - 1), function(l) gamma(l + 1) / gamma(l + 1 - r) * Prob[l + 1]))
  })
}

m <- MoFactorial(Prob, limitek)

# Fit Poisson mixture distribution
theta <- (m[3] - m[1] * m[2]) / (m[2] - m[1]^2)
gamma <- (m[1] * m[3] - m[2]^2) / (m[2] - m[1]^2)
lambda1 <- theta / 2 + sqrt((theta / 2)^2 - gamma)
lambda2 <- theta / 2 - sqrt((theta / 2)^2 - gamma)
omega <- (m[1] - lambda2) / (lambda1 - lambda2)
Poismix <- omega * dpois(k, lambda1) + (1 - omega) * dpois(k, lambda2)
Poismix_rounded <- round(Poismix * n_total, 2)

# ----------------------
# Chi-Squared Goodness-of-Fit Statistic
# ----------------------

Chi <- function(C0, limitek, Probest, p) {
  suma <- 0
  indicetemp <- 0
  for (i in 1:limitek) {
    if (i + 1 <= limitek && Probest[i + 1] < 2) {
      indicetemp <- i
      break
    } else {
      suma <- suma + (C0[i] - Probest[i])^2 / Probest[i]
    }
  }
  suma1 <- if (indicetemp > 0) {
    sumac <- sum(C0[indicetemp:limitek])
    sumap <- sum(Probest[indicetemp:limitek])
    (sumap - sumac)^2 / sumap
  } else 0
  
  chi_stat <- suma + suma1
  df <- (if (indicetemp == 0) limitek else indicetemp) - 1 - p
  c(chi_stat, df)
}

# Compute chi-squared statistics
Chipois <- Chi(C0, limitek, Poisson, 1)
Chibinom <- Chi(C0, limitek, binomial, 1)
ChiBN <- Chi(C0, limitek, BN, 2)
ChiBN_ZM <- Chi(C0, limitek, ZM_NB_dist, 2)
ChiPolya <- Chi(C0, limitek, polyaAeppli_dist, 2)

# ----------------------
# Print Results
# ----------------------
cat("Chi-squared statistics results:\n\n")
cat("Poisson:", Chipois, "\n")
cat("Binomial:", Chibinom, "\n")
cat("Negative Binomial:", ChiBN, "\n")
cat("Zero-Modified NB:", ChiBN_ZM, "\n")
cat("Polya-Aeppli:", ChiPolya, "\n")


# ------------------------
#  Create the plot using ggplot2
# ------------------------

library(ggplot2)

ggplot(plot_data_melted, aes(x = k, y = Frequency,
                             color = Distribution,
                             linetype = Linetype)) +
  geom_line(linewidth = 0.8) +
  
  # Set custom line types for mixture vs. non-mixture
  scale_linetype_manual(values = c("Mixture" = "solid", "Non-Mixture" = "dashed"),
                        name = "Distribution Type") +
  
  # Use a custom color palette for better visual distinction
  scale_color_manual(values = c(
    "Portafolio" = "black",
    "Binomial" = "red",
    "Poisson" = "blue",
    "ZeroModifiedNB" = "purple",
    "PolyaAeppli" = "darkred",
    "Panjer" = "brown",
    "PoissonMixture" = "orange",
    "PanjerMX" = "darkgreen"
  ), name = "Distribution") +
  
  # Customize axis limits and labels
  scale_x_continuous(limits = c(min(k), max(k) + 2)) +
  labs(
    title = "Comparison of Fitted Distributions (Excluding Zero Claims)",
    x = "Number of Claims (k)",
    y = "Frequency"
  ) +
  
  # Apply minimal theme and adjust styles
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),     # Center the plot title
    legend.position = "bottom",                 # Move legend to the bottom
    axis.text = element_text(size = 12),        # Font size for axis ticks
    axis.title = element_text(size = 14),       # Font size for axis titles
    text = element_text(size = 12)              # General text size
  )

