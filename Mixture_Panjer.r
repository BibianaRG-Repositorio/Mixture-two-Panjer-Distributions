library(reshape2)
library(ggplot2)

# Final code

# Define the portfolio C0 with the provided data.
# You can change the portfolio to analyze different datasets.
C0 <- c(122628, 21686, 4014, 832, 224, 68, 17, 7, 7)

# Create a vector 'k' representing the events (0, 1, 2, ...).
k <- c(0:(length(C0) - 1))

# Calculate the length of vector 'k'.
limitek <- length(k)

# Calculate the probabilities of each event.
Prob <- (C0) / sum(C0, na.rm = TRUE)

# Calculate the expectation (expected value) of the distribution.
Esp <- sum((k * C0) / sum(C0, na.rm = TRUE), na.rm = TRUE)

# Calculate the variance of the distribution.
Var <- sum(C0 / sum(C0, na.rm = TRUE) * (k - Esp)^2, na.rm = TRUE)

# Calculate the dispersion index (variance/expectation).
Indice <- Var / Esp

# Binomial Distribution
# Calculate the binomial distribution for each event in 'k'.
binomial <- dbinom(k, 14, Esp / sum(C0, na.rm = TRUE)) * sum(C0, na.rm = TRUE)
# Round the binomial distribution values to 2 decimal places.
binomial_rounded <- round(binomial, 2)
# Calculate the estimated probability 'pest'.
pest <- round(Esp / sum(C0, na.rm = TRUE), 6)
# Calculate the sum of the binomial distribution values.
sum_binomial <- sum(binomial)

# Poisson Distribution
# Calculate the Poisson distribution for each event in 'k'.
Poisson <- dpois(k, Esp) * sum(C0, na.rm = TRUE)
# Round the Poisson distribution values to 2 decimal places.
Poisson_rounded <- round(Poisson, 2)
# Calculate the sum of the Poisson distribution values.
sum_Poisson <- sum(Poisson)

# Negative Binomial Distribution
# Calculate the negative binomial distribution for each event in 'k'.
BN <- dnbinom(k, Esp^2 / (Var - Esp), Esp / Var) * sum(C0, na.rm = TRUE)
# Round the negative binomial distribution values to 2 decimal places.
BN_rounded <- round(BN, 2)
# Calculate the estimated parameters 'r_est' and 'p_est'.
r_est <- Esp^2 / (Var - Esp)
p_est <- Esp / Var
# Calculate the sum of the negative binomial distribution values.
sum_BN <- sum(BN)

# Panjer Distribution
# Calculate the parameters 'qm', 'rmo', 'alpha', and 'beta'.
qm <- 1 - Esp / Var
rmo <- Esp^2 / (Var - Esp)
alpha <- qm
beta <- (rmo - 1) * qm

# Calculate the Panjer distribution for each event in 'k'.
Panjer <- gamma(rmo + k) / (gamma(rmo) * factorial(k)) * alpha^k * (1 - alpha)^rmo * sum(C0, na.rm = TRUE)
# Round the Panjer distribution values to 2 decimal places.
Panjer_rounded <- round(Panjer, 2)
# Calculate the sum of the Panjer distribution values.
sum_Panjer <- sum(Panjer)

# Mixture of distributions

# Necessary functions for estimations

# Function to calculate sample factorial moments.
# Prob is the relative frequency of an event (0,1,2,...)
# limitek is the length of the absolute frequency vector.
# Returns the first 4 factorial moments.
MoFactorial <- function(Prob, limitek) {
  momf <- vector()
  for (r in 1:4) {
    suma <- 0
    for (l in r:(limitek - 1)) {
      a <- gamma(l + 1) / gamma(l + 1 - r) * Prob[l + 1]
      suma <- suma + a
      momf[r] <- suma
    }
  }
  return(momf)
}

# Poisson Mixture
# Calculate sample factorial moments.
m <- MoFactorial(Prob, limitek)
# Calculate parameter 'theta'.
theta <- (m[3] - m[1] * m[2]) / (m[2] - m[1]^2)
# Calculate parameter 'gamma'.
gamma <- (m[1] * m[3] - m[2]^2) / (m[2] - m[1]^2)
# Calculate parameters 'lambda1' and 'lambda2'.
lambda1 <- theta / 2 + sqrt((theta / 2)^2 - gamma)
lambda2 <- theta / 2 - sqrt((theta / 2)^2 - gamma)
# Calculate parameter 'omega'.
omega <- (m[1] - lambda2) / (lambda1 - lambda2)

# Calculate Poisson distribution.
Pois <- c(dpois(k, m[1]) * sum(C0))
# Calculate the Poisson mixture distribution.
Poismix <- omega * c(dpois(k, lambda1) * sum(C0)) + (1 - omega) * c(dpois(k, lambda2) * sum(C0))
# Round the Poisson mixture values to 2 decimal places.
Poismix_rounded <- round(Poismix, 2)

# Negative Binomial Mixture
# Function to calculate Negative Binomial Distribution
# m is the vector of sample moments
# Esp is the sample expectation
# C0 is the portfolio
# k is the vector of events 0,1,2 ...
# Returns the distribution that fits the data
MixBN <- function(m, Esp, C0, k) {
  x <- (7 * m[2]^3 + 4 * m[3]^2 - 3 * m[2] * m[4] + 4 * Esp^2 * m[4] - 12 * Esp * m[2] * m[3])
  y <- (16 * m[2]^3 + 3 * m[3]^2 - 2 * m[2] * m[4] + 5 * Esp^2 * m[4] - 22 * Esp * m[2] * m[3])
  z <- 2 * (6 * m[2]^3 + Esp^2 * m[4] - 6 * Esp * m[2] * m[3])
  w <- (m[2]^3 + m[3]^2 - m[2] * m[4] + Esp^2 * m[4] - 2 * Esp * m[2] * m[3])
  
  pol <- c((z / w), (y / w), (x / w), 1)
  raices <- polyroot(pol)
  
  for (i in 1:3) {
    if (Im(raices[i]) < 0.001) {
      raizmayor <- raices[i]
      break
    }
  }
  
  for (i in 1:3) {
    if (Im(raices[i]) < 0.001 && Re(raices[i]) > Re(raizmayor)) {
      raizmayor <- raices[i]
    }
  }
  r <- Re(raizmayor)
  
  theta <- (r * m[3] - (r + 2) * Esp * m[2]) / ((r + 2) * (r * m[2] - (r + 1) * Esp^2))
  game <- ((r + 1) * Esp * m[3] - (r + 2) * m[2]^2) / ((r + 1) * (r + 2) * (r * m[2] - (r + 1) * Esp^2))
  
  pp1 <- theta / 2 + sqrt((theta / 2)^2 - game)
  p1 <- 1 / (pp1 + 1)
  pp2 <- theta / 2 - sqrt((theta / 2)^2 - game)
  p2 <- 1 / (pp2 + 1)
  
  omega2 <- ((Esp * p2 - r * (1 - p2)) * p1) / (r * (p2 - p1))
  BNMix <- c((gamma(r + k) / (gamma(r) * factorial(k))) * (omega2 * p1^r * (1 - p1)^k + (1 - omega2) * p2^r * (1 - p2)^k) * sum(C0, na.rm = TRUE))
}

# Calculate the negative binomial mixture distribution.
MixturaBN_result <- MixBN(m, Esp, C0, k)
# Round the negative binomial mixture values to 2 decimal places.
MixturaBN_rounded <- round(MixturaBN_result, 2)

# Panjer Mixture r1
# Function to calculate the positive roots of a cubic polynomial
# m is the vector of sample factorial moments
# Returns the estimated value of r
raices_validas <- function(m) {
  x <- (7 * m[2]^3 + 4 * m[3]^2 - 3 * m[2] * m[4] + 4 * m[1]^2 * m[4] - 12 * m[1] * m[2] * m[3])
  y <- (16 * m[2]^3 + 3 * m[3]^2 - 2 * m[2] * m[4] + 5 * m[1]^2 * m[4] - 22 * m[1] * m[2] * m[3])
  z <- 2 * (6 * m[2]^3 + m[1]^2 * m[4] - 6 * m[1] * m[2] * m[3])
  w <- (m[2]^3 + m[3]^2 - m[2] * m[4] + m[1]^2 * m[4] - 2 * m[1] * m[2] * m[3])
  
  pol <- c((z / w), (y / w), (x / w), 1)
  raices <- polyroot(pol)
  
  raices_validas <- vector()
  j <- 1
  for (i in 1:3) {
    if (Im(raices[i]) < 0.001 & Re(raices[i]) >= 0) {
      j <- j + 1
      raices_validas[j] <- raices[i]
    }
  }
  r <- Re(raices_validas)
  r <- r[!is.na(r)]
}

# Calculate the valid roots of the polynomial.
r <- raices_validas(m)

# Function to estimate parameters alpha and omega
# r is the vector of positive roots of a cubic polynomial
# m is the vector of sample factorial moments
# Returns the parameters of the mixture
Parametros <- function(r, m) {
  Matriz_parametros <- matrix(nrow = length(r), ncol = 6)
  for (i in 1:length(r)) {
    theta <- (r[i] * m[3] - (r[i] + 2) * m[1] * m[2]) / ((r[i] + 2) * (r[i] * m[2] - (r[i] + 1) * m[1]^2))
    game <- ((r[i] + 1) * m[1] * m[3] - (r[i] + 2) * m[2]^2) / ((r[i] + 1) * (r[i] + 2) * (r[i] * m[2] - (r[i] + 1) * m[1]^2))
    pp1 <- theta / 2 + sqrt((theta / 2)^2 - game)
    alpha1 <- pp1 / (pp1 + 1)
    pp2 <- theta / 2 - sqrt((theta / 2)^2 - game)
    alpha2 <- pp2 / (pp2 + 1)
    beta1 <- alpha1 * (r[i] - 1)
    beta2 <- alpha2 * (r[i] - 1)
    omega <- ((m[1] * (1 - alpha2) - r[i] * alpha2)) * (1 - alpha1) / (r[i] * (alpha1 - alpha2))
    vector_parametros <- c(r[i], alpha1, beta1, alpha2, beta2, omega)
    Matriz_parametros[i, ] <- vector_parametros
  }
  Matriz_parametros
}

# Estimate the parameters of the Panjer mixture.
parest <- Parametros(r, m)

# Function to calculate the Panjer mixture distribution
MixPanjer <- function(r, alpha1, alpha2, omega) {
  PjMix <- c((gamma(r + k) / (gamma(r) * factorial(k))) * (omega * alpha1^k * (1 - alpha1)^r + (1 - omega) * alpha2^k * (1 - alpha2)^r) * sum(C0))
  PjMix
}

# Obtain the two distributions.
PanjerMX_list <- list()
for (i in 1:nrow(parest)) {
  sol <- MixPanjer(parest[i, 1], parest[i, 2], parest[i, 4], parest[i, 6])
  PanjerMX_list[[i]] <- sol
}

# store the last Panjer mixed distribution in PanjerMX
PanjerMX <- PanjerMX_list[[length(PanjerMX_list)]]
# Round the PanjerMX values to 2 decimal places.
PanjerMX_rounded <- round(PanjerMX, 2)



# Zero-Modified Negative Binomial (ZM-NB)
# Function to calculate the ZM-NB distribution
ZM_NB <- function(C0, k) {
  p0_obs <- C0[1] / sum(C0)  # Observed probability of zero
  if (p0_obs == 1) {
    zm_nb <- ifelse(k == 0, sum(C0), 0)
  } else {
    # Estimate parameters using method of moments (simplified)
    mean_pos <- sum(k * C0[-1]) / sum(C0[-1])
    var_pos <- sum((k - mean_pos)^2 * C0[-1]) / sum(C0[-1])
    
    if (var_pos > mean_pos) {
      r_est <- mean_pos^2 / (var_pos - mean_pos)
      p_est <- mean_pos / var_pos
      p0_nb <- (p_est)^r_est
      omega <- (C0[1] - sum(C0) * p0_nb) / (sum(C0) * (1 - p0_nb))
      
      if (omega < 0 || omega > 1) {
        zm_nb <- rep(NA, length(k)) # Indicate invalid parameters
      } else {
        zm_nb <- ifelse(k == 0, sum(C0) * omega, (1 - omega) * sum(C0) * dnbinom(k[k > 0], r_est, p_est))
        zm_nb[1] <- sum(C0) * omega
      }
    } else {
      zm_nb <- rep(NA, length(k)) # Indicate invalid parameters
    }
  }
  return(zm_nb)
}

ZM_NB_dist <- ZM_NB(C0, k)
ZM_NB_rounded <- round(ZM_NB_dist, 2)

# Pólya-Aeppli Distribution
# Function to calculate the Pólya-Aeppli distribution
Polya_Aeppli <- function(C0, k) {
  mean_val <- sum(k * C0) / sum(C0)
  var_val <- sum((k - mean_val)^2 * C0) / sum(C0)
  
  if (var_val > mean_val) {
    lambda_est <- mean_val
    beta_est <- (var_val - mean_val) / mean_val
    
    polya_aeppli <- vector("numeric", length(k))
    polya_aeppli[1] <- exp(-lambda_est)
    for (i in 1:(length(k) - 1)) {
      polya_aeppli[i + 1] <- (lambda_est / factorial(i)) * sum(exp(-lambda_est) * beta_est^(0:(i)) * (1 - beta_est)^(i - (0:(i))) * choose(i, 0:(i)))
    }
  } else {
    polya_aeppli <- rep(NA, length(k)) # Indicate invalid parameters
  }
  return(polya_aeppli * sum(C0))
}

Polya_Aeppli_dist <- Polya_Aeppli(C0, k)
Polya_Aeppli_rounded <- round(Polya_Aeppli_dist, 2)




# Function to calculate the Chi-squared statistic for each distribution
# C0 the portfolio data
# limitek
# Probest vector of estimated probabilities
# p number of parameters to estimate
# Returns the statistic value
Chi <- function(C0, limitek, Probest, p) {
  indicetemp <- 0
  suma <- 0
  for (i in 1:limitek) {
    if (i + 1 <= limitek & Probest[i + 1] < 2) {
      indicetemp <- i
      break
    } else {
      suma <- suma + (C0[i] - Probest[i])^2 / Probest[i]
    }
  }
  suma1 <- 0
  if (indicetemp > 0) {
    sumac <- 0
    sumap <- 0
    for (i in indicetemp:limitek) {
      sumap <- sumap + Probest[i]
      sumac <- sumac + C0[i]
    }
    suma1 <- (sumap - sumac)^2 / sumap
  }
  if (indicetemp == 0) {
    indicetemp <- limitek
  }
  
  chi_stat <- suma + suma1
  degrees_of_freedom <- (indicetemp - 1 - p)
  c(chi_stat, degrees_of_freedom)
}

# Calculate the Chi-squared statistic for each distribution.
Chipois <- Chi(C0, limitek, Poisson, 1)
Chibinom <- Chi(C0, limitek, binomial, 1)
ChiBN <- Chi(C0, limitek, BN, 2)
ChiBN_ZM <- Chi(C0, limitek,ZM_NB_dist , 2)
ChiPolya <- Chi(C0, limitek, Polya_Aeppli_dist, 2)
ChiPanjer <- Chi(C0, limitek, Panjer, 2)
Chimixpois <- Chi(C0, limitek, Poismix, 1)
ChimixBN <- Chi(C0, limitek, MixturaBN_result, 1)
ChimixPanjer <- Chi(C0, limitek, PanjerMX, 1)

# Print the results
cat("Results:\n\n")

cat("Binomial Distribution:\n")
print(binomial_rounded)
cat("Chi-square", Chibinom, "\n\n")

cat("Poisson Distribution:\n")
print(Poisson_rounded)
cat("Chi-square", Chipois, "\n\n")

cat("Negative Binomial Distribution:\n")
print(BN_rounded)
cat("Chi-square", ChiBN, "\n\n")

cat("Zero-Modified Negative Binomial (ZM-NB):\n")
print(ZM_NB_rounded)
cat("Chi-square", ChiBN_ZM, "\n\n")

cat("Pólya-Aeppli Distributionn:\n")
print(Polya_Aeppli_rounded)
cat("Chi-square", ChiPolya, "\n\n")


cat("Panjer Distribution:\n")
print(Panjer_rounded)
cat("Chi-square", ChiPanjer, "\n\n")

cat("Poisson Mixture Distribution:\n")
print(Poismix_rounded)
cat("Chi-square", Chimixpois, "\n\n")

cat("Negative Binomial Mixture Distribution:\n")
print(MixturaBN_rounded)
cat("Chi-square", ChimixBN, "\n\n")

cat("Panjer Mixture Distribution:\n")
print(PanjerMX_rounded)
cat("Chi-square", ChimixPanjer, "\n\n")

# Create a data frame with all the distributions, including C0
plot_data <- data.frame(
  k = k,
  Portafolio = C0,    # Corrected variable name to C0
  Binomial = binomial_rounded,
  Poisson = Poisson_rounded,
  ZeroModifiedNB = ZM_NB_rounded,
  PolyaAeppli = Polya_Aeppli_rounded,
  Panjer = Panjer_rounded,
  PoissonMixture = Poismix_rounded,
  PanjerMX = PanjerMX_rounded    # Corrected consistency
)

# Filter out the row where k is 0
plot_data_filtered <- plot_data[plot_data$k != 0,]

# Melt the data frame to long format
plot_data_melted <- melt(plot_data_filtered, id.vars = "k",
                         variable.name = "Distribution", value.name = "Frequency")    # Changed value.name to Frequency

# Add a 'linetype' column to distinguish mixture vs. non-mixture
plot_data_melted$Linetype <- ifelse(
  plot_data_melted$Distribution %in% c("PoissonMixture", "PanjerMX"),    # Corrected consistency
  "Mixture",
  "Non-Mixture"
)

# 2. Create the Plot

#   - Use ggplot for flexibility
#   - Plot 'k' on the x-axis, 'Frequency' on the y-axis (corrected)
#   - Use different colors for each distribution
#   - Use 'linetype' to control line style
#   - Use a more distinct color palette
#   - Show legends, but don't duplicate them on the lines

library(ggplot2)

ggplot(plot_data_melted, aes(x = k, y = Frequency, color = Distribution, linetype = Linetype)) +    # Corrected y-axis label
  geom_line(linewidth = 0.8) +
  scale_linetype_manual(values = c("Mixture" = "solid", "Non-Mixture" = "dashed"),
                        name = "Distribution Type") + # Add a name to the linetype legend
  # Use a more distinct color palette
  scale_color_manual(values = c(
    "Portafolio" = "black",    # Corrected variable name to C0
    "Binomial" = "red",
    "Poisson" = "blue",
    "ZeroModifiedNB" = "purple",
    "PolyaAeppli" = "darkred",
    "Panjer" = "brown",
    "PoissonMixture" = "orange",
    "PanjerMX" = "darkgreen"    # Corrected consistency
  ), name = "Distribution") + # Add a name to the color legend
  scale_x_continuous(limits = c(min(k), max(k) + 2)) +  # Extend x-axis for labels
  labs(
    title = "Comparison of Fitted Distributions (Excluding Zero)",
    x = "Number of Claims (k)",
    y = "Frequency"    # Corrected y-axis label
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),    # Center title
    legend.position = "bottom",    # Position the legend at the bottom
    axis.text = element_text(size = 12),   # Adjust axis text size
    axis.title = element_text(size = 14),  # Adjust axis title size
    text = element_text(size = 12)        # Adjust general text size
  )
