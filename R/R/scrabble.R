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
#' @param nIter the maximum iterations, the default is 20.
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
#' result <- scrabble(demo_data,parameter = parameter)
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
                     nIter = 20,
                     error_out_threshold = 1e-4,
                     nIter_inner = 20,
                     error_inner_threshold = 1e-4){

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

  # browser()
  if(length(data)==2){
    # prepare the matrices for scrabble
    print(cat("\n detected SCRABBLE data"))
    D <- matrix(1, nrow = 1, ncol = n_row)
    A <- getA(D, beta, gamma, n_row)
    B <- getB(D, Z, Y, beta)
  }else{
    # prepare the matrices for multi-task
    print(cat("\n detected mtSCRABBLE data"))
    Z <- t(data[[2]]) # this is the actual bulk samples
    D <- data[[3]]
    # A <- t(data[[3]])%*%data[[3]] +gamma*diag(nrow(Y))
    A <- getA(D, beta, gamma, n_row)
    # B <- getB(data[[3]], t(Z), Y, beta)
    B <- getB(D, Z, Y, beta)
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

  while((k < nIter) & (error > error_out_threshold)){

    # update the X
    X <- newX
    Y <- newY

    gamma_Y_B_Lambda <- gamma*Y + B - Lambda

    newX <- ToDense(newX)

    newX <- cDescent(gamma_Y_B_Lambda,
                     A,
                     zones,
                     newX,
                     nIter_inner,
                     error_inner_threshold)

    newX <- ToSparse(newX)

    # STV
    S <- getS(newX, Lambda, gamma)

    tau <- alpha/gamma

    # print(cat("\n dim(X)= \n"))
    # print(dim(X))
    #
    # print(cat("\n ",class(S),"::S= \n",head(as.matrix(S)),"\n"))
    #
    # print(cat("\n min S: \n",min(S),"\n"))
    # print(cat("\n min newX: \n",min(newX),"\n"))
    # print(cat("\n min gamma_Y_B_Lambda: \n",min(gamma_Y_B_Lambda),"\n"))
    # print(cat("\n min Y: \n",min(Y),"\n"))
    # print(cat("\n min X: \n",min(X),"\n"))
    # print(cat("\n min Lambda: \n",min(Lambda),"\n"))
    #
    # print(cat("\n max S: \n",max(S),"\n"))
    # print(cat("\n max newX: \n",max(newX),"\n"))
    # print(cat("\n max gamma_Y_B_Lambda: \n",max(gamma_Y_B_Lambda),"\n"))
    # print(cat("\n max Y: \n",max(Y),"\n"))
    # print(cat("\n max X: \n",max(X),"\n"))
    # print(cat("\n max Lambda: \n",max(Lambda),"\n"))
    #
    # print(cat("\n norm S: \n",Matrix::norm(S),"\n"))
    # print(cat("\n norm newX: \n",Matrix::norm(newX),"\n"))
    # print(cat("\n norm gamma_Y_B_Lambda: \n",Matrix::norm(gamma_Y_B_Lambda),"\n"))
    # print(cat("\n norm Y: \n",Matrix::norm(Y),"\n"))
    # print(cat("\n norm X: \n",Matrix::norm(X),"\n"))
    # print(cat("\n norm Lambda: \n",Matrix::norm(Lambda),"\n"))

    result <- svt(S, lambda = tau)

    s <- result[[1]] - tau

    s[s < 0] <- 0
    newY <- getY(s, result[[2]], result[[3]])

    error <- norm(as.matrix(log10(X + 1) - log10(newX + 1)), type = c("2"))/norm(as.matrix(log10(X + 1)), type = c("2"))

    if (k == 0){error = 1}
    k <- k + 1
    Lambda <- updateLambda(Lambda, newX, newY, gamma)
    print(cat("\n error=",error,"\n"))
  }
  print(paste0('Imputation is finished'))

  # browser()
  return(recoverData(newX))

}
