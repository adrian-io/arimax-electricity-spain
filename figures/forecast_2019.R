# Recreates figures/forecast_2019.png: out-of-sample forecasts for 2019 from
# model 1, SARIMA(2,0,1)(0,1,1)12, without (arima.Rmd) and with (arimax.Rmd)
# calendar effects and outlier correction. Run from the repository root:
#   Rscript figures/forecast_2019.R

source("CalendarEffects.r")
source("atipics2.R")

serie <- window(ts(read.table("ConsumElec.dat"), start = 1990, freq = 12))
ultim <- c(2018, 12)
serie2 <- window(serie, end = ultim)

# SARIMA on log demand, fit through Dec 2018, intervals back-transformed
lnserie2 <- log(serie2)
mod <- arima(lnserie2, order = c(2, 0, 1), seasonal = list(order = c(0, 1, 1), period = 12), include.mean = TRUE)
p <- predict(mod, n.ahead = 12)
pr <- ts(c(tail(lnserie2, 1), p$pred), start = ultim, freq = 12)
se <- ts(c(0, p$se), start = ultim, freq = 12)
sa <- list(pr = exp(pr), tl = exp(pr - 1.96 * se), tu = exp(pr + 1.96 * se))

# ARIMAX: trading days + Easter, outliers detected and removed, refit through Dec 2018
data <- c(start(serie)[1], start(serie)[2], length(serie))
wTradDays <- Wtrad(data)
wEast <- Weaster(data)
modEC <- arima(serie, order = c(2, 0, 1), seasonal = list(order = c(0, 1, 1), period = 12),
               xreg = data.frame(wTradDays, wEast), include.mean = TRUE)
atip <- outdetec(modEC, dif = c(1, 12), crit = 2.8, LS = TRUE)
serie.lin <- lineal(serie, atip$atip)
modEC2 <- arima(window(serie.lin, end = ultim), order = c(2, 0, 1), seasonal = list(order = c(0, 1, 1), period = 12),
                xreg = data.frame(window(wTradDays, end = ultim), window(wEast, end = ultim)), include.mean = TRUE)
p <- predict(modEC2, n.ahead = 12, newxreg = window(cbind(wTradDays, wEast), start = c(ultim[1] + 1, 1)))
wLS <- sum(atip$atip[atip$atip$type_detected == "LS" & atip$atip$Obs <= length(serie) - 12, 3])  # add level shifts back
pr <- ts(c(tail(serie2, 1), p$pred + wLS), start = ultim, freq = 12)
se <- ts(c(0, p$se), start = ultim, freq = 12)
ax <- list(pr = pr, tl = pr - 1.96 * se, tu = pr + 1.96 * se)

ink <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"
grid_col <- "#e1e0d9"; axis_col <- "#c3c2b7"; surf <- "#fcfcfb"
blue <- "#2a78d6"; orange <- "#eb6834"
wash <- function(col) adjustcolor(col, alpha.f = 0.14)

obs  <- window(serie, start = c(2016, 1))
obs19 <- window(serie, start = c(2019, 1))
ylim <- range(obs, sa$tl, sa$tu, ax$tl, ax$tu)
yt <- pretty(ylim, n = 5); ylim <- range(yt)
lab_k <- function(v) paste0(v / 1000, "K")

metrics <- function(fc) {
  pr <- window(fc$pr, start = c(2019, 1)); tl <- window(fc$tl, start = c(2019, 1)); tu <- window(fc$tu, start = c(2019, 1))
  c(mape = 100 * mean(abs(obs19 - pr) / obs19), mae = mean(abs(obs19 - pr)), inside = sum(obs19 >= tl & obs19 <= tu))
}

panel <- function(fc, col, ttl) {
  m <- metrics(fc)
  plot(NA, xlim = c(2016, 2020), ylim = ylim, axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  abline(h = yt, col = grid_col, lwd = 0.45)
  abline(v = 2019 - 1/24, col = axis_col, lwd = 0.45)
  text(2019 - 1/24 + 0.04, ylim[2] - 0.03 * diff(ylim), "forecast", adj = c(0, 1), col = muted, cex = 0.72)
  x <- as.numeric(time(fc$pr))
  polygon(c(x, rev(x)), c(fc$tl, rev(fc$tu)), col = wash(col), border = NA)
  lines(as.numeric(time(obs)), obs, col = ink, lwd = 0.9)
  lines(x, fc$pr, col = col, lwd = 0.9)
  points(as.numeric(time(obs19)), obs19, pch = 21, bg = ink, col = surf, cex = 0.95, lwd = 1.1)
  axis(2, at = yt, labels = lab_k(yt), las = 1, col = NA, col.ticks = NA, col.axis = muted, cex.axis = 0.78, line = -0.6)
  axis(1, at = 2016:2020, labels = FALSE, col = axis_col, col.ticks = axis_col, lwd = 0.45, tck = -0.015)
  axis(1, at = 2016:2019 + 0.5, labels = 2016:2019, tick = FALSE, col.axis = muted, cex.axis = 0.78, line = -0.7)
  mtext(ttl, side = 3, line = 1.35, adj = 0, cex = 0.86, col = ink, font = 2)
  mtext(sprintf("MAPE %.2f%%,  MAE %.0f GWh,  %d of 12 months inside the 95%% interval", m["mape"], m["mae"], m["inside"]),
        side = 3, line = 0.35, adj = 0, cex = 0.74, col = ink2)
}

png("figures/forecast_2019.png", width = 2400, height = 1080, res = 220, bg = surf)
layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE), heights = c(0.2, 1))
par(bg = surf, family = "sans")

par(mar = c(0, 2.6, 0.5, 1)); plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
text(0, 0.9, "Monthly electricity demand in Spain (GWh): out-of-sample forecasts for 2019, models fit on 1990 to 2018", adj = c(0, 1), cex = 0.95, col = ink, font = 2, xpd = NA)
key_y <- 0.28; kx <- 0
key_line <- function(x, col, label, dot = FALSE) {
  segments(x, key_y, x + 0.025, key_y, col = col, lwd = 0.9, xpd = NA)
  if (dot) points(x + 0.0125, key_y, pch = 21, bg = col, col = surf, cex = 0.95, lwd = 1.1, xpd = NA)
  text(x + 0.033, key_y, label, adj = c(0, 0.5), cex = 0.78, col = ink2, xpd = NA)
}
key_line(0.00, ink, "Actual demand", dot = TRUE)
key_line(0.17, blue, "SARIMA forecast")
key_line(0.35, orange, "ARIMAX forecast")
rect(0.535, key_y - 0.09, 0.5475, key_y + 0.09, col = wash(blue), border = NA, xpd = NA)
rect(0.5475, key_y - 0.09, 0.56, key_y + 0.09, col = wash(orange), border = NA, xpd = NA)
text(0.568, key_y, "95% prediction interval", adj = c(0, 0.5), cex = 0.78, col = ink2, xpd = NA)

par(mar = c(2.2, 2.6, 3.0, 1.2))
panel(sa, blue, expression(bold("SARIMA(2,0,1)(0,1,1)"[12] * " on log demand")))
panel(ax, orange, expression(bold("ARIMAX: + trading days, Easter, outlier correction")))
dev.off()
