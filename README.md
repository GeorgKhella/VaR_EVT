# Value-at-Risk Estimation using GARCH and Extreme Value Theory (EVT)

This project implements a semi-parametric model in R to estimate and forecast Value-at-Risk (VaR) for a low-dimensional financial portfolio. The approach integrates GARCH modeling for volatility dynamics and Extreme Value Theory (EVT) to capture tail risks.

## 📈 Objective

To build a risk estimation pipeline that captures both:
- **Volatility clustering** using AR(0)-GARCH(1,1)
- **Tail behavior** via Generalized Pareto Distribution (GPD)

## 🗃️ Dataset

- **Assets:** General Electric (GE), Coca-Cola (KO), and 3M (MMM)
- **Period:** January 3, 2000 – December 31, 2024
- **Frequency:** Daily
- **Weights:** GE 55%, KO 25%, MMM 20%

Portfolio returns are aggregated using log-loss formulation.

## 🔧 Methodology

### 1. Exploratory Analysis
- Time series visualization
- Ljung–Box test for autocorrelation
- Summary statistics (mean, skewness, kurtosis)

### 2. Volatility Modeling
- Fitting AR-GARCH models using `rugarch`
- Model selection via AIC/BIC
- Residual diagnostics and normality checks

### 3. Tail Risk Modeling with EVT
- POT method using `evd` or `evir`
- Threshold selection with mean excess and shape stability plots
- GPD fitting for excesses over the 95th percentile

### 4. VaR Forecasting
- In-sample and out-of-sample rolling forecasts
- Recursive model re-fitting with rolling window
- Semi-parametric quantile estimation

### 5. Backtesting
- Binomial test for VaR violation frequency
- Evaluation of predictive coverage at 95% and 99% confidence levels

## 📁 Project Structure

