

# creating the function for Orthogonal Regression

ortho_reg <- function(y, x, alpha)
  {
  
  mean_x <- mean(x, na.rm = TRUE)
  
  mean_y <- mean(y, na.rm = TRUE)
  
  n <- length(x)
  
  Sxx <- sum((x - mean_x) ^ 2)
  
  Syy <- sum((y - mean_y) ^ 2)
  
  Sxy <- sum((y - mean_y) * (x - mean_x))
  
  beta_1 <- (Syy - Sxx + sqrt((Syy - Sxx) ^ 2 + 4 * (Sxy ^ 2 ))) / ( 2 * Sxy)
  
  beta_0 <- mean_y - beta_1 * mean_x
  
  eigen_result <- eigen(rbind(c(Sxx, Sxy), c(Sxy, Syy)))
  
  l1 <- eigen_result$value[1]
  
  l2 <- eigen_result$value[2]
  
  if(l1 < l2)
  {
  l3 <- l2
  l2 <- l1
  l1 <- l3
  }
  
  r_2 <- ((l1 - l2) ^ 2 ) / ((l1 + l2) ^ 2)
  
  o <- atan(beta_1)
  
  if(alpha > 1) alpha <- alpha / 100
  
  t <- asin(sqrt(dchisq(0.05, 1) / ((n - 1) * (l1 / l2 + l2 / l1 - 2))))
  
  lci_beta_1 <- tan(o - t)
  
  uci_beta_1 <- tan(o + t)
  
  library(boot)
  
  r <- prcomp(~x + y)
  
  slope <- r$rotation[2, 1] / r$rotation[1, 1]
  
  intercept <- r$center[2] - slope * r$center[1]
  
  # boostrap analysis
  
  boot_data <- data.frame(x = x, y = y)
  
  stat_intercept <- function(data, indices)
  {
    r <- prcomp(~x + y, data = data, subset = indices)
    
    slope <- r$rotation[2, 1] / r$rotation[1, 1]
    
    intercept <- r$center[2] - slope * r$center[1]
    
    return(intercept)
    
  }
  
  reps_intercept <- boot(boot_data, stat_intercept, R = 1000)
  
  ci_intercept <- boot.ci(reps_intercept, type = c("perc", "bca"))
  
  lci_b0_boot <- ci_intercept$percent[4]
  
  uci_b0_boot <- ci_intercept$percent[5]
  
  stat_slope <- function(data, indices)
  {
    r <- prcomp(~x + y, data = data, subset = indices)
    
    slope <- r$rotation[2, 1] / r$rotation[1, 1]
    
    intercept <- r$center[2] - slope * r$center[1]
    
    return(slope)
    
  }
  
  reps_slope <- boot(boot_data, stat_slope, R = 1000)
  
  ci_slope <- boot.ci(reps_slope, type = c("perc", "bca"))
  
  lci_b1_boot <- ci_slope$percent[4]
  
  uci_b1_boot <- ci_slope$percent[5]
  
  # plotting 
  
  par(mfrow = c(1, 2))
  
  y_pred <- beta_0 + beta_1 * x
  
  orto_res <- y_pred - y
  
  plot(y_pred, 
       y,
       pch = 21,
       col = "grey",
       bg = "grey",
       cex = 2,
       xlab = "Predicted Values", 
       ylab = "Observed Values",
       main = "Predicted vs Observed")
  
  reg_pred_obs <- lm(y ~ y_pred)
  
  lines(y_pred, reg_pred_obs$coefficients[1] + reg_pred_obs$coefficients[2] * y_pred, 
        col = "red", 
        lty = 2)
  
  abline(a = 0, 
         b = 1,
         col = "blue")
  
  legend("topleft", 
         legend = c("Adjst. Line", "Equality Line"),
         col = c("red", "blue"), lty = 1:2, cex = 0.8)
  
  # plotting residuals
  
  hist(orto_res, 
       main = "Orthogonal Residuals",
       xlab = "Residuals")
  
  par(mfrow = c(1, 1))
  
  return(list("Intercept" = beta_0, 
              "Slope" = beta_1, 
              "Slope C.I." = c(lci_beta_1, uci_beta_1), 
              "R-squared" = r_2, 
              "Intercept C.I. Boostrap" = c(lci_b0_boot, uci_b0_boot), 
              "Slope C.I. Boostrap" = c(lci_b1_boot, uci_b1_boot),
              "Predicted" = y_pred,
              "Residuals" = orto_res))
  }


