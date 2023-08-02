# Plots fitted dlmm results

dlmm_plot <- function(pred, x1 = "", x2 = "", type = "marginal") {

  if (x1 %in% names(pred$main) & x2 == "") {
    
    dat <- data.frame(x = pred$lag_vals[[x1]], y = pred[[type]][[x1]]$fit,
                      ymin = pred[[type]][[x1]]$low, ymax = pred[[type]][[x1]]$high)
    
    p <- ggplot(dat, aes(x = x, y = exp(y), ymin = exp(ymin), ymax = exp(ymax))) +
      geom_ribbon(fill = "grey") +
      geom_line() +
      labs(x = "Lag", y = "Effect")
    
  } else if (x1 %in% pred$x_names && x2 %in% pred$x_names) {
    
    if (!(paste0(x1, "-", x2) %in% names(pred$interaction))) {
      xtmp <- x2
      x2 <- x1
      x1 <- xtmp
    }
    
    n <- paste0(x1, "-", x2)
    dat <- data.frame(x = rep(pred$lag_vals[[x1]], each = length(pred$lag_vals[[x2]])),
                      y = rep(pred$lag_vals[[x2]], length(pred$lag_vals[[x1]])),
                      z = c(pred$interaction[[paste0(x1, "-", x2)]]$fit))
    
    p <- ggplot(dat, aes(x = x, y = y, z = z, fill = z)) +
      geom_tile() +
      scale_fill_viridis() +
      labs(x = x1, y = x2)
    
  }

  p <- p +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  
  return(p)
  
}
