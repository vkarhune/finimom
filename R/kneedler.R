kneedler <- function(y, h2lead, Sens = 1){
  # y <- h2_svd[seq_len(ceiling(p_eff_keepinds1))]

  x <- seq_along(y)

  # normalise
  minx <- 1
  maxx <- length(y)
  miny <- min(y)
  maxy <- max(y)

  xn <- (x - minx) / (maxx - minx)
  yn <- (y - miny) / (maxy - miny)

  # plot(xn, yn)

  xd <- xn
  yd <- yn - xn

  # xl <- xd
  # only those points that are above the lead variant h2
  whichyl <- which(yd > c(0, yd[-maxx]) & yd > c(yd[-1], 0) & y > h2lead)
  # whichyl <- whichyl[y[whichyl] > h2lead]
  
  yl <- yd[whichyl]
  xl <- xd[whichyl]

  tlmx <- yl - Sens * mean(diff(xn))

  #plot(xd, yd)
  #abline(v = xd[whichyl], lty = 2, col = "grey")
  #abline(h = tlmx, lty = 2, col = "grey")


  cond <- FALSE
  ind <- 1
  while(!(cond)){
    tstart <- whichyl[ind]
    tstop <- ifelse(ind != length(tlmx), whichyl[ind + 1] - 1, length(y))
    
    if( any(yd[seq(tstart, tstop)] < tlmx[ind]) ){
      cond <- TRUE
      knee_ind <- tstart
    } else {
      if(ind == length(tlmx)){
        knee_ind <- NA
        cond <- TRUE
      } else {
        ind <- ind + 1
      }
    }
  }

  out <- list("knee_x" = knee_ind, "knee_y" = unique(y[knee_ind]))

  return(out)

}

