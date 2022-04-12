# Code adapted from https://github.com/tanlabcode/SCRABBLE:
# 
# Copyright 2019 Tao Peng, 2021 Matt Karikomi

#' Runs SCRABBLE
#'
#' This package imputes drop-out data by optimizing an objective function that consists of three terms.
#' The first term ensures that imputed values for genes with nonzero expression remain as close to their
#' original values as possible, thus minimizing unwanted bias towards expressed genes. The second term ensures
#' the rank of the imputed data matrix to be as small as possible. The rationale is that we only expect a
#' limited number of distinct cell types in the samples. The third term operates on the bulk RNA-Seq data.
#' It ensures consistency between the average gene expression of the aggregated imputed data and the
#' average gene expression of the bulk RNA-Seq data. We developed a convex optimization algorithm to minimize
#' the objective function.
#'
#'
#' @param data the input data list. The input
#' data is a list of two datasets, scRNAseq and bulk RNAseq.
#'
#' @param parameter the vector of parameters. The first parameter is the value of alpha in the mathematical model
#' , the second one is the value of beta in the mathematical model.
#'
#' @param nIter_outer the maximum iterations, the default is 20.
#'
#' @param error_out_threshold the threshold of the error between the current imputed matrix and the previous one.
#' Default is 1e-4.
#'
#' @param nIter_inner the maximum interations of calculating the sub-optimization problem. Default is 20.
#'
#' @param error_inner_threshold the threshold of the error between the current updated matrix and the previous one.
#' Default is 1e-4.
#'
#' @examples
#' # Set up the parameter used in SCRABBLE
#' parameter <- c(1, 1e-6, 1e-4)
#'
#' # Run SCRABLE
#' demo_data <- list(matrix(runif(1000000),1000,1000),matrix(runif(5000),1000,5))
#' param <- c(1,1e-6,1e-4)
#' result <- scrabble(demo_data,parameter = param)
#'
#' @return A data matrix with the same size of the input scRNAseq data
#'
#' @rdname SCRABBLE
#'
#' @export
#'
#'
scrabble <- function(data,
                     parameter,
                     nIter_outer = 20,
                     error_out_threshold = 1e-4,
                     nIter_inner = 20,
                     nSDCIters = 10000,
                     error_inner_threshold = 1e-4,
                     multi=FALSE,
                     saveImputedStep=FALSE,
                     useIrlba=TRUE){
  
  # Use the sparse matrix to store the matrix
  Y <- as(t(as.matrix(data[[1]])), "dgCMatrix")
  
  # Transpose the data matrix for the optimization
  zones <- (Y > 0)*1
  n_row <- nrow(Y)
  n_col <- ncol(Y)
  
  print(paste0('SCRABBLE begins the imputation of the data with ',n_col,' genes and ', n_row, ' cells'))
  # Define the parameters
  alpha <- parameter[1]
  beta <- parameter[2]
  gamma <- parameter[3]
  
  # prepare the bulk data
  if (isempty(data[[2]])){
    beta <- 0*beta
    Z <- matrix(1, nrow = 1, ncol = n_col)
  }else{
    Z <- getZ(as.matrix(data[[2]]), n_row)
  }
  
  if(!multi){
    print("detected SCRABBLE data")
    D <- matrix(1, nrow = 1, ncol = n_row)
    A <- getA(D, beta, gamma, n_row)
    B <- getB(D, Z, Y, beta)
  }else{
    print("detected mtSCRABBLE data")
    Z <- t(data[[2]]) # this is the actual bulk samples
    D <- data[[3]] # Bulk x Cells
    A <- getA(D, beta, gamma, n_row) # Cells x Cells
    B <- getB(D, Z, Y, beta) # Bulk x Cells
  }
  
  # initialize the X,Y,and Lambda for the iteration
  X <- Y
  newX <- X
  newY <- as.matrix(Y)
  Lambda <- matrix(0, nrow = n_row, ncol = n_col)
  
  # set up the thresholds for iterations and the errors
  k <- 0
  error <- 1
  
  print(paste0('Imputation initialization is finished'))
  print(paste0('... ....'))

  Start=Sys.time()
  while((k < nIter_outer) & (error > error_out_threshold)){
    init_start_time = Sys.time()

    # update the X
    X <- newX
    Y <- newY
    
    gamma_Y_B_Lambda <- gamma*Y + B - Lambda
    
    newX <- ToDense(newX)
    init_end_time = Sys.time()
    print(paste0('End outer init, i=', k,' time: ',init_end_time-init_start_time))
    
    cd_start_time = Sys.time()
    newX <- cDescent(gamma_Y_B_Lambda,
                     A,
                     zones,
                     newX,
                     nIter_inner,
                     nSDCIters,
                     error_inner_threshold,
                     k)
    cd_end_time = Sys.time()
    print(paste0('End descent, i=', k,' time: ',cd_end_time-cd_start_time))
    
    presvt_start_time = Sys.time()
    newX <- ToSparse(newX)
    S <- getS(newX, Lambda, gamma)
    tau <- alpha/gamma
    presvt_end_time = Sys.time()
    print(paste0('End pre-SVT, i=', k,' time: ',presvt_end_time-presvt_start_time))
    
    svt_start_time = Sys.time()
    result <- svt(S, lambda = tau)
    svt_end_time = Sys.time()
    print(paste0('End SVT i=', k,' time: ',svt_end_time-svt_start_time))
    
    gety_start_time = Sys.time()
    s <- result[[1]] - tau
    s[s < 0] <- 0
    newY <- getY(s, result[[2]], result[[3]])
    gety_end_time = Sys.time()
    print(paste0('End getY i=', k,' time: ',gety_end_time-gety_start_time))

    geterr_start_time = Sys.time()
    if(useIrlba){
      error <- irlba::irlba(log10(X+1) - log10(newX+1))["d"][[1]][1]/irlba::irlba(log10(X + 1))["d"][[1]][1]
    }else{
      error <- norm(as.matrix(log10(X + 1) - log10(newX + 1)), type = c("2"))/norm(as.matrix(log10(X + 1)), type = c("2"))
    }
    geterr_end_time = Sys.time()
    print(paste0('End get err i=', k,' time: ',geterr_end_time-geterr_start_time))
    
    if (k == 0){error = 1}
    k <- k + 1
    updatelambda_start_time = Sys.time()
    Lambda <- updateLambda(Lambda, newX, newY, gamma)
    updatelambda_end_time = Sys.time()
    print(paste0('End update lambda i=', k,' time: ',updatelambda_end_time-updatelambda_start_time))
    
    cat("\n error=",error,"\n")
  }
  print(paste0('Imputation is finished'))
  
  return(list("error"=error,"imputed"=recoverData(newX)))
}