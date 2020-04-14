## ----setup, include=TRUE----------------------------------------------------------------------------------------------
knitr::opts_chunk$set(cache=TRUE, fig.width=7, fig.height=5.5)
options(width=120)

## ---------------------------------------------------------------------------------------------------------------------
knitr::knit_exit()
SlicedIndex <-
  function(cal.pts,
           spdatos = NULL,
           from = NULL,
           to = NULL,
           base = NULL,
           frm = NULL,
           vars = NULL,
           bw = 1000,
           area.radius = 5 * bw,
           date = NULL,
           slice.by = NULL)
  {
    require(zoo)
    if (is.null(date))
      stop("'date' not specified with no default.")
    if (is.null(slice.by))
      stop("'slice.by' not specified with no default.")
    #
    # Convert observations to selected time interval
    #
    spdatos@data[, date] <- sapply(spdatos@data[, date], FUN = slice.by)
    spdatos <- spdatos[(spdatos@data[, date] >= from) &
                         (spdatos@data[, date] <= to), ]
    #
    # Pick only required variables, and cases complete for such variables
    #
    spdatos <- spdatos[complete.cases(spdatos@data), ]
    slices  <- sort(unique(spdatos@data[, date]))
    #
    # For each calibration point, take a subset of observations within
    # area.radius distance (projected coordinates are assumed)
    #
    xy      <- coordinates(cal.pts)
    preds   <- as.data.frame(matrix(0, length(slices), nrow(cal.pts)))
    colnames(preds) <- rownames(cal.pts)
    preds   <- zoo(preds, order.by = slices)
    
    for (i in nrow(cal.pts)) {
      sel <- rep(FALSE, nrow(spdatos@data))
      sel <- sel | (colSums((t(coordinates(spdatos)) - xy[i, ]) ^ 2) <
                      area.radius ^ 2)
      for (j in seq_along(slices)) {
        subsel <- sel & (spdatos@data[, date] == slices[j])
        cat("Processing slice:  ", as.character(slices[j]), "\n")
        cat(sum(subsel), "\n")
        preds[j, ] <- gwr(
          formula = frm,
          bandwidth = bw,
          data = spdatos[subsel, ],
          fit.points = cal.pts,
          predict = TRUE
        )$SDF@data$"(Intercept)"
      }
    }
    return(preds)
  }

