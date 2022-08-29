# start_parallel / stop_parallel function ------------------------------
# https://github.com/jepusto/Pusto/blob/main/R/parallel-setup.R
start_parallel <- function (cores, source_obj = NULL, packages = NULL, 
                            setup = "plan") {
  
  if (requireNamespace("snow", quietly = TRUE)) {
    
    if (!is.null(snow::getMPIcluster())) {
      
      # TACC setup
      
      cl <- snow::getMPIcluster()
      
      if (setup == "plan") {
        future::plan(future::cluster, workers = cl)
      }
      if (setup == "register") {
        if (!requireNamespace("doSNOW", quietly = TRUE)) {
          stop("The doSNOW package is required for registering the cluster. Please install it.", call. = FALSE)
        }     
        doSNOW::registerDoSNOW(cl)
      } 
      if (!is.null(source_obj)) {
        snow::clusterExport(cl, source_obj)
      }
      if (!is.null(packages)) {
        library_calls <- lapply(packages, function(lib) call("library", lib))
        snow::clusterExport(cl, "library_calls", envir = environment())
        snow::clusterEvalQ(cl, lapply(library_calls, eval))
      }
      return(cl)
    }
  }     
  
  if (missing(cores)) cores <- parallel::detectCores() - 1
  
  if (!is.na(pmatch("Windows", Sys.getenv("OS")))) {
    
    
    # Windows setup
    if (requireNamespace("multidplyr", quietly = TRUE)) {
      cl <- multidplyr::create_cluster(cores = cores)
    } else {
      cl <- parallel::makePSOCKcluster(cores)
      cat("Don't forget to use stop_parallel() to close the cluster.")
    }
    
    if (setup == "plan") {
      future::plan(future::cluster, workers = cl)
    }
    
    if (setup == "register") {
      if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("The doParallel package is required for registering the cluster. Please install it.", call. = FALSE)
      }     
      doParallel::registerDoParallel(cl)
    }
    
    if (!is.null(source_obj)) {
      parallel::clusterExport(cl, source_obj)
    }
    if (!is.null(packages)) {
      library_calls <- lapply(packages, function(lib) call("library", lib))
      parallel::clusterExport(cl, "library_calls", envir = environment())
      parallel::clusterEvalQ(cl, lapply(library_calls, eval))
    }
    return(cl)
    
  } else {
    
    # Mac setup
    
    if (setup == "plan") {
      future::plan(future::multiprocess)
    }
    if (setup == "register") {
      if (!requireNamespace("doParallel", quietly = TRUE)) {
        stop("The doParallel package is required for registering the cluster. Please install it.", call. = FALSE)
      }
      doParallel::registerDoParallel(cores = cores)
    }
    return(NULL)
  }
}

stop_parallel <- function(cluster = NULL) parallel::stopCluster(cluster)
