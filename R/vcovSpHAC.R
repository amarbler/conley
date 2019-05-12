## call the new version vcovSpHAC
vcovSpHAC <- function(reg,
                      unit = NULL,
                      time = NULL,
                      lat,
                      lon,
                      kernel = c("bartlett","uniform"),
                      dist_fn = c("haversine", "spherical", "chord", "flatearth"),
                      dist_cutoff,
                      lag_cutoff = 0,
                      verbose = FALSE,
                      balanced_pnl = FALSE,
                      ncores = NA,
                      maxobsmem = 50000L) {

  ## check the options on kernel and dist_fn
  kernel <- match.arg(kernel)
  dist_fn <- match.arg(dist_fn)

  noFEs <- length(unit) == 0

  ## assign the number of cores to use, default to all
  if(is.na(ncores)) ncores <- RcppParallel::defaultNumThreads()
  RcppParallel::setThreadOptions(numThreads = ncores)

  Fac2Num <- function(x) {
    as.numeric(as.character(x))
  }

  if (class(reg) == "felm") {
  #  if(length(reg$fe)==0) {
   if(noFEs == TRUE) {
       unit <- "fe1"
      time <- "fe2"
    }
    Xvars <- rownames(reg$coefficients)
    # covers xsectional case now
    dt = data.table::data.table(
      reg$cY,
      reg$cX,
      #fe1 = Fac2Num(reg$fe[[1]]),
      #fe2 = Fac2Num(reg$fe[[2]]),
      fe1 = ifelse( rep(!noFEs,length(reg$cY)) , Fac2Num(reg$fe[[1]]), 1:length(reg$cY) ),
      fe2 = ifelse( rep(!noFEs,length(reg$cY)), Fac2Num(reg$fe[[2]]), rep(1,length(reg$cY)) ),
      expand.model.felm(model = reg, extras = c(lat,lon), na.expand=T)
    #    coord1 = Fac2Num(reg$clustervar[[1]]),
     #  coord2 = Fac2Num(reg$clustervar[[2]])
    )
    # covers xsectional case now
    data.table::setnames(dt, c("fe1", "fe2"),
                        c( ifelse( !noFEs, names(reg$fe)[1], unit),
                           ifelse( !noFEs, names(reg$fe)[2], time)
                         ))
   # data.table::setnames(dt,
    #                     c("fe1", "fe2", "coord1", "coord2"),
     #                    c(names(reg$fe)[1:2], names(reg$clustervar)))
    dt = dt[, e := as.numeric(reg$residuals)]

  } else {
    message("Model class not recognized.")
    return(NULL)
  }

  n <- nrow(dt)
  k <- length(Xvars)

  # Renaming variables:
  orig_names <- c(unit, time, lat, lon)
  new_names <- c("unit", "time", "lat", "lon")
  data.table::setnames(dt, orig_names, new_names)

  # Empty Matrix:
  XeeX <- matrix(nrow = k, ncol = k, 0)

  #================================================================
  # Correct for spatial correlation:
  timeUnique <- unique(dt[, time])
  Ntime <- length(timeUnique)
  data.table::setkey(dt, time)

  if (verbose) message("Starting to loop over time periods...")

  if (balanced_pnl) {

    sub_dt <- dt[time == timeUnique[1]]
    lat <- sub_dt[, lat]
    lon <- sub_dt[, lon]
    rm(sub_dt)

    if (verbose) message("Computing Distance Matrix using RcppParallel...")
    d <- DistMatPar(cbind(lat, lon), cutoff = dist_cutoff, kernel, dist_fn)
    rm(list = c("lat", "lon"))

  }

  # this computes the sandwich as many time as there a time periods, once in the case of an xsection
  XeeXhs <- lapply(timeUnique, function(t) {
    iterateObs(dt, Xvars,
               sub_index = t,
               type = "spatial",
               cutoff = dist_cutoff,
               balanced_pnl = balanced_pnl, d=d,
               kernel = kernel,
               dist_fn = dist_fn,
               verbose = verbose,
               ncores = ncores,
               maxobsmem = maxobsmem)
  })


  if(balanced_pnl){ rm(d) }
  # First Reduce:
  XeeX <- Reduce("+",  XeeXhs)

  X <- as.matrix(dt[, eval(Xvars), with = FALSE])
  invXX <- solve(t(X) %*% X) * n

  if(Ntime > 1) {
    #================================================================
    # Correct for serial correlation:
    panelUnique <- unique(dt[, unit])
    Npanel <- length(panelUnique)
    data.table::setkey(dt, unit)

    if (verbose) {
      message("Starting to loop over units...")
    }

    XeeXhs <- lapply(panelUnique, function(t) {
      iterateObs(dt, Xvars,
                 sub_index = t,
                 type = "serial",
                 cutoff = lag_cutoff,
                 balanced_pnl = balanced_pnl,
                 kernel = kernel,
                 dist_fn = dist_fn,
                 verbose = verbose,
                 ncores = ncores,
                 maxobsmem = maxobsmem)
    })

    XeeX_serial <- Reduce("+",  XeeXhs)
    XeeX <- XeeX + XeeX_serial
    #================================================================
  }

  V_spatial_HAC <- invXX %*% (XeeX / n) %*% invXX / n
  V_spatial_HAC <- (V_spatial_HAC + t(V_spatial_HAC)) / 2
  return(V_spatial_HAC)
}

# now call the different options from within iterateObs()
iterateObs <- function(dt, Xvars, sub_index, type, cutoff, balanced_pnl, d, verbose, kernel, dist_fn, ncores, maxobsmem) {
  k <- length(Xvars)
  if (type == "spatial" & balanced_pnl) {
    sub_dt <- dt[time == sub_index]
    n1 <- nrow(sub_dt)
    if (n1 > 1000 &
        verbose) {
      message(paste("Balanced panel or cross-section. Starting on time period:", sub_index))
    }

    X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE])
    e <- sub_dt[, e]

    if (ncores > 1) {
      if (verbose) message("Computing spatial sandwich with RcppParallel...")
      XeeXhs <- Bal_XeeXhC_Par(d, X, e)
    } else {
      if (verbose) message("Computing spatial sandwich serially...")
      XeeXhs <- Bal_XeeXhC(d, X, e, n1, k)
    }

  } else if (type == "spatial" & !balanced_pnl) {
    sub_dt <- dt[time == sub_index]
    n1 <- nrow(sub_dt)
    if (n1 > 1000 &
        verbose) {
      message(paste("Starting on sub index:", sub_index))
    }

    X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE])
    e <- sub_dt[, e]
    lat <- sub_dt[, lat]
    lon <- sub_dt[, lon]

    # If n1 >= maxobsmem, then avoiding construction of distance matrix.
    # This requires more operations, but is MUCH less memory intensive.
    if (n1 < maxobsmem) {
      if (ncores > 1) {
        if (verbose) message("Computing distance matrix and spatial sandwich with RcppParallel...")
        XeeXhs <- XeeXhC_Par(cbind(lat, lon), cutoff, X, e, kernel, dist_fn)
      } else {
        if (verbose) message("Computing distance matrix and spatial sandwich serially...")
        XeeXhs <- XeeXhC(cbind(lat, lon), cutoff, X, e, n1, k, kernel, dist_fn)
      }
    } else {

      if (ncores > 1) {
        if (verbose) message("Computing distance vectors and spatial sandwich with RcppParallel...")
        XeeXhs <- XeeXhC_Lg_Par(cbind(lat, lon), cutoff, X, e, kernel = kernel, dist_fn = dist_fn)
      } else {
        if (verbose) message("Computing distance vectors and spatial sandwich serially...")
        XeeXhs <- XeeXhC_Lg(cbind(lat, lon), cutoff, X, e, n1, k,  kernel = kernel, dist_fn = dist_fn)
      }
    }

  } else if (type == "serial") {
    sub_dt <- dt[unit == sub_index]
    n1 <- nrow(sub_dt)
    if (n1 > 1000 &
        verbose) {
      message(paste("Starting on sub index:", sub_index))
    }

    X <- as.matrix(sub_dt[, eval(Xvars), with = FALSE])
    e <- sub_dt[, e]
    times <- sub_dt[, time]

    XeeXhs <- TimeDist(times, cutoff, X, e, n1, k)
  }

  XeeXhs
}

expand.model.felm <- function(model, extras, envir = environment(formula(model)),
                              na.expand = FALSE) {

  topaste <- c(names(model$fe), names(model$clustervar), extras)
  fescluext <- parse(text = paste("~", paste(topaste,
                                             collapse = "+")))[[1L]]

  data <- eval(model$call$data, envir)
  ff <- foo ~ bar + baz

  ff[[2L]] <- parse(text = paste("~", model$lhs))[[1L]][[2L]]
  ff[[3L]][[2L]] <- formula(model)[[2L]]
  ff[[3L]][[3L]] <- fescluext[[2L]]

  if (!na.expand) {
    naa <- model$call$na.action
    subset <- model$call$subset
    rval <- eval(call("model.frame", ff, data = data, subset = subset,
                      na.action = naa), envir)[, extras]
  }
  else {
    subset <- model$call$subset
    rval <- eval(call("model.frame", ff, data = data, subset = subset,
                      na.action = I), envir)
    oldmf <- model.frame(model)
    keep <- match(rownames(oldmf), rownames(rval))
    rval <- rval[keep, extras]
    class(rval) <- "data.frame"
  }
  return(rval)
}
