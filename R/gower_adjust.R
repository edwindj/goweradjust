#' @export
gower_adjust <- function(data, rules, weight=NULL, ...
            , cl = NULL
            , n = 10
            , Ncpus = getOption("Ncpus", 1)
            , timeout=60
  ){
    mip <- miprules(rules)

    vars <- mip$._vars
    nr_rows <- nrow(data)
    nr_cols <- ncol(data)
    names_cols <- colnames(data)

    missing_vars <- vars[!vars %in% names(data)]

    if (length(missing_vars)){
      # if they are part of the environment remove...
      mv_in_env <- sapply(missing_vars, exists)
      vars <- setdiff(vars, missing_vars[mv_in_env])
      mip$._vars <<- vars

      missing_vars <- missing_vars[!mv_in_env]
    }

    if (length(missing_vars)){
      stop('Missing column(s): '
           , paste0("'", missing_vars, "'", collapse = ", ")
           , ". Add them to your data and rerun."
           , call. = FALSE
      )
    }

    numvars_mip <- mip$._vars_num
    numvars_data <- names(data)[sapply(data, is.numeric)]
    numvars <- numvars_mip[numvars_mip %in% numvars_data]

    # checking size of numeric columns
    too_big <- sapply(data[numvars], function(v){
      any(abs(v) > 1e7)
    })

    if (isTRUE(any(too_big))){
      nv <- numvars[too_big]
      data[nv] <- sapply(data[nv], function(v){
        is.na(v) <- abs(v) > 1e7
        v
      })

      warning("Large values detected in: "
              , paste0("'", names(too_big)[too_big] ,"'", collapse = ", ")
              , ". Values > abs(1e7) were set to NA. "
              , "\nThis might be indication that these column(s) should be rescaled."
              , "\n(the problem because otherwise numerically unstable)"
              , call. = FALSE
      )
    }

    categorical_as_integer <- numvars_data[ numvars_data %in% vars
                                            & !numvars_data %in% numvars_mip
    ]
    if (length(categorical_as_integer)){
      stop('Categorical columns coded as integers: '
           , paste0("'", categorical_as_integer, "'", collapse = ", ")
           , ". Change them to `factor` and rerun."
           , call. = FALSE
      )
    }

    if (length(weight) == 0){
      weight <- matrix(1, nrow=nr_rows, ncol=nr_cols)
      colnames(weight) <- names_cols
    } else {
      if (is.null(dim(weight))){
        if (length(weight) == ncol(data)){
          if (!is.null(names(weight))){
            weight <- weight[names_cols]
          }
          # use recycling to fill a weight matrix and transpose...
          weight <- t(matrix(weight, nrow=nr_cols, ncol=nr_rows))
          colnames(weight) <- names_cols
        }
      }
      stopifnot(dim(weight) == dim(data))
      if (is.null(colnames(weight))){
        colnames(weight) <- names_cols
      } else {
        weight <- weight[, names_cols, drop=FALSE]
      }
      stopifnot(names(weight) == names_cols)
    }

    N <- nr_rows
    rows <- seq_len(N)


    # filter for records that are valid..., that reduces the processing
    # time considerably
    cf <- validate::confront(data, rules)
    agg <- aggregate(cf, by = "record")
    invalid <- (agg$nfail + agg$nNA) > 0

    n_invalid <- sum(invalid)

    #TODO add ref data !!!

    solve_record <- function(r, progress = invisible){
      starttime <- Sys.time()
      el <- tryCatch({
        values <- as.list(data[r,,drop=FALSE])
        mip$set_values( values = values
                        , weight[r,]
        )
        mip$execute(timeout=timeout, ...)
      }, error = function(e){
        list(solution=FALSE, s = NA, adapt = logical())
      })

      if (!isTRUE(el$solution)){
        # test for numerical instability?
        # retry because of numerical instability
        mip$set_values( values = values
                        , weight[r,]
        )
        # could be a scaling issue. Drop geometric scaling
        warning("Dropping geometric `scaling` for record ",r, " (?lpSolveAPI::lp.control.options)"
                , call. = FALSE
        )
        args <- list(...)
        args$timeout = timeout
        args$scaling = c("range", "equilibrate","integers" )
        el <- do.call(mip$execute, args)
        if (!isTRUE(el$solution)){
          dump_path <- file.path(
            tempdir(),
            paste0("no_solution_record_", r, ".mps")
          )
          mip$write_lp( dump_path,type="mps")
          warning( "dumping lp problem for record ", r
                   , " in '", dump_path, "'"
                   , call. = FALSE
          )
        }
      }

      # remove lp object, too memory hungry...
      el$lp <- NULL
      el$duration <- Sys.time() - starttime
      el$row <- r
      progress()
      el
    }
    show_progress <- interactive() && n_invalid > 0

    Ncpus <- min(Ncpus, n_invalid)
    setup_cluster <- Ncpus > 1 && is.null(cl) && n_invalid > 1

    if (setup_cluster){
      cl <- parallel::makeCluster(Ncpus, setup_strategy = "sequential")
      on.exit(parallel::stopCluster(cl))
    }

    sols <- if (is.null(cl)){ # sequential
      progress <- invisible
      if (show_progress) {
        pb <- txtProgressBar(min=0, max = n_invalid, style=3)
        e <- new.env(parent = emptyenv())
        e$val <- 0
        e$lasttime <- Sys.time()

        progress <- function(){
          e$val <- e$val + 1
          if (Sys.time() - e$lasttime > 1){
            setTxtProgressBar(pb, e$val)
            e$lasttime <- Sys.time()
          }
        }
      }
      sols <- lapply(rows[invalid], solve_record, progress = progress)
      if (show_progress){
        setTxtProgressBar(pb, n_invalid)
        close(pb)
      }
      sols
    } else if (is.numeric(cl)){
      message("starting parallel: ", cl, " cores (non-Windows)")
      parallel::mclapply( rows[invalid]
                          , solve_record
                          , mc.cores = cl
      )
    } else if (inherits(cl, "cluster")){
      message("Setting up ", class(cl)[1], " job with ", length(cl)," nodes")
      parallel::clusterEvalQ(cl, library(errorlocate))
      parallel::clusterExport( cl
                               , c("data", "weights", "mip")
                               , envir = environment()
      )
      message("Starting cluster job...")
      # parallel::parLapplyLB(cl, rows[invalid], solve_record)
      parallel::parLapplyLB(cl, rows[invalid], solve_record)
    }
    # collect info during processing
    status <- integer(N)
    duration <- numeric(N)
    solution <- logical(N)
    solution[] <- TRUE

    # adapt <- t(res)
    adapt <- matrix( FALSE
                     , nrow = nr_rows
                     , ncol = nr_cols
                     , dimnames = list(NULL, names_cols)
    )

    for (sol in sols){
      r <- sol$row
      idx <- match(names(sol$adapt), names_cols)
      adapt[r,idx] <-  sol$adapt
      duration[r] <- sol$duration
      solution[r] <- sol$solution
      status[r] <- sol$s
    }


    idx <- which(colnames(adapt) %in% colnames(weight))

    wpr <- weight
    wpr[!adapt | is.na(adapt)] <- 0

    weight_per_record <-  rowSums(wpr, na.rm=T)
    if (any(!solution)){
      warning("For some records the procedure was unsuccesful, "
              , "please check the '$solution' and '$status' of the errorlocations.\n"
              , "Records: "
              , paste0(which(!solution), collapse = " ,")
              , call. = FALSE
      )
    }

    is.na(adapt) <- is.na(data)
    create_errorlocation(
      values = adapt,
      weight = weight_per_record,
      duration = duration,
      status = status,
      solution = solution
    )
  s
}



