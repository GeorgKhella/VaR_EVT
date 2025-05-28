load("/Users/georgkhella/Downloads/Data/dataset_03.Rdata")

# Visualizza la struttura dell'oggetto dataset
str(dataset)

# Se è una lista, possiamo accedere agli elementi così:
prices <- dataset$prices
weights <- dataset$weights
dates.calib <- dataset$dates.calib
dates.test <- dataset$dates.test

start_date <- dates.calib$Begin
end_date <- dates.calib$End

# Filtra i prezzi nel periodo di calibrazione
library(dplyr)
prices_calib <- prices %>%
  filter(Date >= start_date & Date <= end_date)

# Calcola i log-prezzi e le log-loss giornaliere
log_prices <- log(prices_calib[, c("GE", "KO", "MMM")])
log_losses <- dplyr::lag(log_prices) - log_prices  # ℓ_t = log(P_{t-1}) - log(P_t)

# Calcola la log-loss del portafoglio con i pesi fissi
w <- as.numeric(weights[1, ])
log_loss_portfolio <- as.numeric(as.matrix(log_losses) %*% w)

# Rimuovi NA (prima riga)
log_loss_portfolio <- na.omit(log_loss_portfolio)
dates <- prices_calib$Date[-1]

# Serie temporale
log_loss_ts <- ts(log_loss_portfolio, frequency = 252)

# Statistiche descrittive
summary(log_loss_ts)
sd(log_loss_ts)
library(e1071)
skewness(log_loss_ts)
kurtosis(log_loss_ts)

# Plot delle log-loss
plot(dates, log_loss_ts, type = "l", xlab = "Date", ylab = "Log-loss")

# Correlogrammi
acf(
  log_loss_ts,
  main = "",            # nessun titolo
  lwd = 2,
  col = "darkblue",
  xlab = "Lag",
  ylab = "ACF"
)

acf(abs(log_loss_ts), main = "")

# Test Ljung–Box (per autocorrelazione e clustering)
Box.test(log_loss_ts, lag = 10, type = "Ljung-Box")
Box.test(abs(log_loss_ts), lag = 10, type = "Ljung-Box")




# Carica il pacchetto
library(rugarch)

# Funzione per stimare modelli AR(k)-GARCH(1,1)
fit_ar_garch_model <- function(arOrder, series) {
  spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(arOrder, 0), include.mean = TRUE),
    distribution.model = "norm"
  )
  fit <- ugarchfit(spec = spec, data = series)
  return(fit)
}

# Stima i modelli con AR(0), AR(1), AR(2)
fit_ar0 <- fit_ar_garch_model(0, log_loss_ts)
fit_ar1 <- fit_ar_garch_model(1, log_loss_ts)
fit_ar2 <- fit_ar_garch_model(2, log_loss_ts)

# Tabella AIC e BIC
ic_table <- data.frame(
  Model = c("AR(0)-GARCH(1,1)", "AR(1)-GARCH(1,1)", "AR(2)-GARCH(1,1)"),
  AIC = c(infocriteria(fit_ar0)[1], infocriteria(fit_ar1)[1], infocriteria(fit_ar2)[1]),
  BIC = c(infocriteria(fit_ar0)[2], infocriteria(fit_ar1)[2], infocriteria(fit_ar2)[2])
)
print(ic_table)

# Optional: stampa riassunti dettagliati
show(fit_ar0)
show(fit_ar1)
show(fit_ar2)



# Residui standardizzati
res_std <- residuals(fit_ar0, standardize = TRUE)
res_std_vec <- as.numeric(res_std)

# Plot della serie dei residui
plot(dates, res_std, type = "l", ylab = "Residual", xlab = "Time")

# ACF dei residui
acf(res_std, main = "ACF of standardized residuals")

# ACF dei residui al quadrato
acf(res_std^2, main = "ACF of squared standardized residuals")

# Test Ljung–Box su residui standardizzati
Box.test(res_std, lag = 10, type = "Ljung-Box")

# Test Ljung–Box su residui standardizzati al quadrato
Box.test(res_std^2, lag = 10, type = "Ljung-Box")

# QQ plot (normalità dei residui)
qqnorm(res_std, main = "")
qqline(res_std, col = "blue")

# Test di normalità (opzionale)
shapiro.test(res_std_vec)

# Carica il pacchetto
library(tseries)

# Assicurati che i residui siano un vettore numerico
res_std_vec <- as.numeric(residuals(fit_ar0, standardize = TRUE))

# Esegui il test di Jarque–Bera
jarque.bera.test(res_std_vec)




# STEP EVT: Model innovation tails with Generalized Pareto Distribution (GPD)

# Carica il pacchetto necessario
library(evir)

# Residui standardizzati dal modello AR(0)–GARCH(1,1)
z <- res_std_vec

# -----------------------------
# 1. Plot diagnostici per scelta soglia u
# -----------------------------

# Mean Residual Life Plot (Mean Excess Plot)
mrlplot(z, main = "")

# Plot della stabilità dei parametri della GPD
thresholds <- quantile(z, probs = seq(0.90, 0.99, by = 0.01))
shape_estimates <- sapply(thresholds, function(u) {
  fit <- gpd(z, threshold = u)
  fit$par.ests["xi"]
})
plot(thresholds, shape_estimates, type = "b",
     main = "", xlab = "Threshold u", ylab = "Estimated xi")

# -----------------------------
# 2. Fit della GPD sopra la soglia scelta
# -----------------------------

# Scegli la soglia (puoi modificarla in base ai grafici precedenti)
u <- quantile(z, 0.95)

# Fit GPD ai residui z sopra la soglia u
# Carica il pacchetto corretto
library(extRemes)

# Fit GPD con extRemes
gpd_fit <- fevd(z, threshold = u, type = "GP", method = "MLE")

# Plot separati
plot(gpd_fit, type = "qq", main = "")        # QQ-plot classico
plot(gpd_fit, type = "qq2")       # QQ-plot eccedenze
plot(gpd_fit, type = "hist", main = "")      # Istogramma con fit
plot(gpd_fit, type = "density")   # Densità stimata

# Estrai parametri stimati
beta_hat <- gpd_fit$results$par["scale"]
xi_hat   <- gpd_fit$results$par["shape"]


# -----------------------------
# 3. Calcolo dei quantili delle innovazioni
# -----------------------------

# Probabilità di superamento della soglia
n <- length(z)
nu <- sum(z > u)
Fu_bar <- nu / n  # coda empirica oltre u

# Funzione quantile per GPD (VaR a livello alpha)
quantile_gpd <- function(alpha, u, beta, xi, Fu_bar) {
  u + (beta / xi) * (((1 - alpha) / Fu_bar)^(-xi) - 1)
}

# Calcola quantili per innovazioni a livello 0.95 e 0.99
v_0.95 <- quantile_gpd(0.95, u, beta_hat, xi_hat, Fu_bar)
v_0.99 <- quantile_gpd(0.99, u, beta_hat, xi_hat, Fu_bar)

# Output dei quantili stimati
cat("Quantile v_0.95 (innovations):", v_0.95, "\n")
cat("Quantile v_0.99 (innovations):", v_0.99, "\n")

# -----------------------------
# 4. Calcolo del Value-at-Risk in-sample
# -----------------------------

# Estrai la stima della volatilità condizionata σ_t dal modello GARCH
sigma_t <- sigma(fit_ar0)  # vettore della stessa lunghezza di z

# Calcola il VaR giornaliero in-sample (negativo perché perdita)
VaR_0.95 <- -sigma_t * v_0.95
VaR_0.99 <- -sigma_t * v_0.99


# Plot della serie di log-loss vs soglie VaR
plot(dates, log_loss_ts, type = "l", col = "black", lwd = 1.5,
     ylab = "Log-loss", xlab = "Date", main = "In-sample VaR")
lines(dates, VaR_0.95, col = "red", lwd = 1.2, lty = 2)
lines(dates, VaR_0.99, col = "blue", lwd = 1.2, lty = 2)
legend("topright", legend = c("Log-loss", "VaR 95%", "VaR 99%"),
       col = c("black", "red", "blue"), lty = c(1, 2, 2), lwd = 1.2)

# -----------------------------
# 5. Verifica in-sample (backtesting)
# -----------------------------

# Conta quante volte le perdite superano la soglia VaR
exceed_95 <- sum(log_loss_ts < VaR_0.95)
exceed_99 <- sum(log_loss_ts < VaR_0.99)

# Output
cat("Exceedances (VaR 95%):", exceed_95, "su", length(log_loss_ts), "giorni\n")
cat("Exceedances (VaR 99%):", exceed_99, "su", length(log_loss_ts), "giorni\n")


######## Step A6: VaR prediction for testing period ########

# Librerie necessarie
library(dplyr)
library(xts)
library(rugarch)
library(extRemes)

# Caricamento dataset
load("/Users/georgkhella/Downloads/Data/dataset_03.Rdata")

# Estrai oggetti dal dataset
prices <- dataset$prices
weights <- dataset$weights
dates.calib <- dataset$dates.calib
dates.test <- dataset$dates.test

# Date
start_date <- dates.calib$Begin
end_date   <- dates.calib$End
start_test <- dates.test$Begin
end_test   <- dates.test$End

# Calcola log-loss dell’intero periodo
log_prices_all <- log(prices[, c("GE", "KO", "MMM")])
log_losses_all <- dplyr::lag(log_prices_all) - log_prices_all
w <- as.numeric(weights[1, ])
log_loss_portfolio_all <- na.omit(as.numeric(as.matrix(log_losses_all) %*% w))
dates_all <- prices$Date[-1]  # rimuove il primo NA

# Serie temporale xts
log_loss_xts <- xts(log_loss_portfolio_all, order.by = dates_all)

# Lunghezza periodo rolling
calib_len <- sum(dates_all >= start_date & dates_all <= end_date)
test_len  <- sum(dates_all > end_date & dates_all <= end_test)
start_idx <- calib_len

# Inizializza vettori
VaR_95_roll <- numeric(test_len)
VaR_99_roll <- numeric(test_len)
actuals     <- numeric(test_len)

# Barra di avanzamento
pb <- txtProgressBar(min = 0, max = test_len, style = 3)

# Rolling forecast loop
for (i in 1:test_len) {
  setTxtProgressBar(pb, i)  # aggiorna barra
  
  # Finestra mobile di calibrazione
  window_data <- log_loss_xts[(i):(i + calib_len - 1)]
  
  # Specifica e stima GARCH(1,1)
  spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "norm"
  )
  fit <- ugarchfit(spec, data = window_data, solver = "hybrid")
  
  # Residui standardizzati
  z <- as.numeric(residuals(fit, standardize = TRUE))
  
  # EVT: GPD su coda superiore
  u <- quantile(z, 0.95)
  gpd_fit <- fevd(z, threshold = u, type = "GP", method = "MLE")
  beta_hat <- gpd_fit$results$par["scale"]
  xi_hat   <- gpd_fit$results$par["shape"]
  Fu_bar   <- mean(z > u)
  
  # Quantili delle innovazioni
  q_evt <- function(alpha) {
    u + (beta_hat / xi_hat) * (((1 - alpha) / Fu_bar)^(-xi_hat) - 1)
  }
  q_95 <- q_evt(0.95)
  q_99 <- q_evt(0.99)
  
  # Forecast GARCH t+1
  forecast <- ugarchforecast(fit, n.ahead = 1)
  mu_hat    <- fitted(forecast)
  sigma_hat <- sigma(forecast)
  
  # VaR semiparametrico
  VaR_95_roll[i] <- -(mu_hat + sigma_hat * q_95)
  VaR_99_roll[i] <- -(mu_hat + sigma_hat * q_99)
  
  # Osservazione reale
  actuals[i] <- as.numeric(log_loss_xts[i + calib_len])
}

close(pb)  # chiudi barra

# Costruisci dataframe finale
dates_test_seq <- dates_all[(start_idx + 1):(start_idx + test_len)]
VaR_out_sample <- data.frame(
  Date     = dates_test_seq,
  Actual   = actuals,
  VaR_0.95 = VaR_95_roll,
  VaR_0.99 = VaR_99_roll
)

# Salva per A7
save(VaR_out_sample, file = "VaR_out_sample.RData")



##PLOT##

# Carica i risultati se non già caricati
load("VaR_out_sample.RData")

library(ggplot2)

# Crea un data.frame lungo per ggplot
library(tidyr)
VaR_plot_df <- VaR_out_sample %>%
  pivot_longer(cols = c("VaR_0.95", "VaR_0.99"),
               names_to = "VaR_Level", values_to = "VaR") %>%
  mutate(Violation = Actual < VaR)

# Plot
ggplot(VaR_plot_df, aes(x = Date)) +
  geom_line(aes(y = Actual), color = "black", size = 0.7) +
  geom_line(aes(y = VaR, color = VaR_Level), linetype = "dashed") +
  geom_point(data = subset(VaR_plot_df, Violation),
             aes(y = Actual), color = "red", size = 1.5, shape = 4) +
  scale_color_manual(values = c("VaR_0.95" = "blue", "VaR_0.99" = "darkgreen"),
                     labels = c("VaR 95%", "VaR 99%")) +
  labs(title = "",
       y = "Log-loss", x = "Date", color = "VaR Level") +
  theme_minimal()

######## A7 ########

load("VaR_out_sample.RData")

# Numero totale di osservazioni out-of-sample
n <- nrow(VaR_out_sample)

# Coverage analysis
violations_95 <- VaR_out_sample$Actual < VaR_out_sample$VaR_0.95
violations_99 <- VaR_out_sample$Actual < VaR_out_sample$VaR_0.99

count_viol_95 <- sum(violations_95)
count_viol_99 <- sum(violations_99)

expected_viol_95 <- round(0.05 * n, 2)
expected_viol_99 <- round(0.01 * n, 2)

cat("=== Coverage Analysis ===\n")
cat("Level: 95%\n")
cat("Observed violations:", count_viol_95, "out of", n, "\n")
cat("Expected violations:", expected_viol_95, "\n\n")

cat("Level: 99%\n")
cat("Observed violations:", count_viol_99, "out of", n, "\n")
cat("Expected violations:", expected_viol_99, "\n\n")

# Binomial tests
cat("=== Binomial Test ===\n")
cat("VaR 95%:\n")
binom_95 <- binom.test(count_viol_95, n, p = 0.05, alternative = "two.sided")
print(binom_95)

cat("\nVaR 99%:\n")
binom_99 <- binom.test(count_viol_99, n, p = 0.01, alternative = "two.sided")
print(binom_99)



