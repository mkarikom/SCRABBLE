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

#' Runs scLDA
#'
#' MCEM version of mtSCRABBLE (multi-task SCRABBLE) using LDA to get the deconvolution weight for A
#'
#' @param data_dir where the sequencing data
#'
#' @param init dropout = don't preprocess the sc data, scr = run og scrabble before the first lda iteration
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
#' @rdname scLDA
#'
#' @export
#'
#'
scLDA <- function(data_dir,
                  ncores,
                  parameter,
                  init='none',
                  nIter = 20,
                  error_out_threshold = 1e-4,
                  nIter_inner = 20,
                  error_inner_threshold = 1e-4){
  # get the data
  T_fn = file.path(data_dir,"T.csv")
  C_fn = file.path(data_dir,"C.csv")
  pDataC_fn = file.path(data_dir,"pDataC.csv")

  data_T = read.csv(T_fn,row.names=1)
  data_C = read.csv(C_fn,row.names=1)
  data_pDataC = read.csv(pDataC_fn,row.names=1)

  # if(init=='scr'){
  #   cat("\n running pre-initiailization via SCRABBLE...")
  #   # Use the sparse matrix to store the matrix
  #   Y <- as(t(as.matrix(data_C)), "dgCMatrix")
  #
  #   # Transpose the data matrix for the optimization
  #   zones <- (Y > 0)*1
  #   n_row <- nrow(Y)
  #   n_col <- ncol(Y)
  #
  #   print(paste0('SCRABBLE begins the imputation of the data with ',n_col,' genes and ', n_row, ' cells'))
  #   # Define the parameters
  #   alpha <- parameter[1]
  #   beta <- parameter[2]
  #   gamma <- parameter[3]
  #
  #   Z <- getZ(as.matrix(data_T), n_row)
  #
  #   D <- matrix(1, nrow = 1, ncol = n_row)
  #   A <- getA(D, beta, gamma, n_row)
  #   B <- getB(D, Z, Y, beta)
  #
  #   # initialize the X,Y,and Lambda for the iteration
  #   X <- Y
  #   newX <- X
  #   newY <- as.matrix(Y)
  #   Lambda <- matrix(0, nrow = n_row, ncol = n_col)
  #
  #   # set up the thresholds for iterations and the errors
  #   k <- 0
  #   error <- 1
  #
  #   print(paste0('Initial scrabble initialization is finished'))
  #   print(paste0('... ....'))
  #
  #   while((k < nIter) & (error > error_out_threshold)){
  #     # update the X
  #     X <- newX
  #     Y <- newY
  #     gamma_Y_B_Lambda <- gamma*Y + B - Lambda
  #     newX <- ToDense(newX)
  #     newX <- cDescent(gamma_Y_B_Lambda,
  #                      A,
  #                      zones,
  #                      newX,
  #                      nIter_inner,
  #                      error_inner_threshold)
  #     newX <- ToSparse(newX)
  #     # STV
  #     S <- getS(newX, Lambda, gamma)
  #     tau <- alpha/gamma
  #     result <- svt(S, lambda = tau)
  #     s <- result[[1]] - tau
  #     s[s < 0] <- 0
  #     newY <- getY(s, result[[2]], result[[3]])
  #     error <- norm(as.matrix(log10(X + 1) - log10(newX + 1)), type = c("2"))/norm(as.matrix(log10(X + 1)), type = c("2"))
  #     if (k == 0){error = 1}
  #     k <- k + 1
  #     Lambda <- updateLambda(Lambda, newX, newY, gamma)
  #     print(cat("\n error=",error,"\n"))
  #   }
  #   xhat = recoverData(newX)
  #   print(paste0('Initial scrabble imputation is finished'))
  #   cat("\n finished pre-initiailization...")
  # }else if(init=='none'){
  #   cat("\n no pre-initialization of sc data, running scLDA...\n")
  # }

  ## run lda to initialize the coefficients

  # nbulk = dim(data_T)[2]
  # ncells = dim(data_C)[2]
  # genes_bulk = rownames(data_T)
  # genes_sc = rownames(data_C)
  #
  # bulk_ind_thresh = which(colSums(as.matrix(data_T) > 0) <= bulk_thresh) # bulk samples expressing too few genes
  # sc_ind_thresh = which(colSums(as.matrix(data_C) > 0) <= sc_thresh) # cells expressing too few genes
  # bulkgene_ind_thresh = length(which(rowSums(as.matrix(data_T) > 0) <= bulk_gene_thresh)) # genes expressed in too few bulk samples
  # scgene_ind_thresh = length(which(rowSums(as.matrix(data_C) > 0) <= sc_gene_thresh)) # genes expressed in too few cells
  # if(n_bulk_pos < nbulk){
  #   cat('\n zero-expression bulk samples detected, removing ',nbulk-n_bulk_pos,' bulk samples\n')
  # }
  # if(n_sc_pos < ncells){
  #   cat('\n zero-expression sc samples detected, removing ',ncells-n_sc_pos,' sc samples\n')
  # }
  # if(n_bulkgene_pos < ngenes){
  #   cat('\n zero-expression genes detected in bulk, removing ',ngenes-n_bulkgene_pos,' genes\n')
  # }
  # if(n_scgene_pos < ngenes){
  #   cat('\n zero-expression genes detected in sc, removing ',ngenes-n_scgene_pos,' genes\n')
  # }

  system2('/home/au/code/DTMwork/local/E_step.sh', args=c(data_dir,toString(ncores)))

  # ## construct the coeffient matrix
  # # # get the coefficient matrix for regression
  # batchinds = lapply(unique(sim@colData$Batch),function(x) which(sim@colData$Batch==x))
  # # the cell type percentages for each bulk sample
  # deconv = lapply(batchinds, function(x) table(sim@colData$Group[x])/length(x))
  # A = matrix(0,nbulk,ncells)
  # for(i in 1:nbulk){
  #   for(j in 1:ncells){
  #     # a[i,j] = % celltype(cell j) in bulk sample i / number of cells
  #     # such that A*X = D' where D[i,] is the average expression per cell of all genes in bulk sample i
  #     A[i,j] = deconv[[i]][sim@colData$Group[j]] / ncells
  #   }
  # }
  #
  # # browser()
  # # define the data list for the simulation data
  # # indcluding: data_true: true data
  # #          data_dropout: dropout data
  # #             data_bluk: bulk data
  # #      percentage_zeros: dropout rate
  # #                 group: the group label
  #
  # data = list()
  # data$data_bulk = data_bulk
  # data$data_A = A

  # Use the sparse matrix to store the matrix
  # Y <- as(t(as.matrix(data_C)), "dgCMatrix")
  #
  # # Transpose the data matrix for the optimization
  # zones <- (Y > 0)*1
  # n_row <- nrow(Y)
  # n_col <- ncol(Y)
  #
  # print(paste0('SCRABBLE begins the imputation of the data with ',n_col,' genes and ', n_row, ' cells'))
  # # Define the parameters
  # alpha <- parameter[1]
  # beta <- parameter[2]
  # gamma <- parameter[3]
  #
  # # prepare the matrices for multi-task
  # Z <- t(data_T) # this is the actual bulk samples
  # D <- data[[3]]
  # # A <- t(data[[3]])%*%data[[3]] +gamma*diag(nrow(Y))
  # A <- getA(D, beta, gamma, n_row)
  # # B <- getB(data[[3]], t(Z), Y, beta)
  # B <- getB(D, Z, Y, beta)
  #
  # # initialize the X,Y,and Lambda for the iteration
  # X <- Y
  # newX <- X
  # newY <- as.matrix(Y)
  # Lambda <- matrix(0, nrow = n_row, ncol = n_col)
  #
  # # set up the thresholds for iterations and the errors
  # k <- 0
  # error <- 1
  #
  # print(paste0('Imputation initialization is finished'))
  # print(paste0('... ....'))
  #
  # while((k < nIter) & (error > error_out_threshold)){
  #
  #   # update the X
  #   X <- newX
  #   Y <- newY
  #
  #   gamma_Y_B_Lambda <- gamma*Y + B - Lambda
  #
  #   newX <- ToDense(newX)
  #
  #   newX <- cDescent(gamma_Y_B_Lambda,
  #                    A,
  #                    zones,
  #                    newX,
  #                    nIter_inner,
  #                    error_inner_threshold)
  #
  #   newX <- ToSparse(newX)
  #
  #   # STV
  #   S <- getS(newX, Lambda, gamma)
  #
  #   tau <- alpha/gamma
  #
  #   result <- svt(S, lambda = tau)
  #
  #   s <- result[[1]] - tau
  #
  #   s[s < 0] <- 0
  #   newY <- getY(s, result[[2]], result[[3]])
  #
  #   error <- norm(as.matrix(log10(X + 1) - log10(newX + 1)), type = c("2"))/norm(as.matrix(log10(X + 1)), type = c("2"))
  #
  #   if (k == 0){error = 1}
  #   k <- k + 1
  #   Lambda <- updateLambda(Lambda, newX, newY, gamma)
  #   print(cat("\n error=",error,"\n"))
  # }
  # print(paste0('Imputation is finished'))
  #
  # # browser()
  # return(recoverData(newX))

}
