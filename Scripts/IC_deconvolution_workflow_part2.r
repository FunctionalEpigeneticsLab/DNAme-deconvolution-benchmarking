## Load in data
args = commandArgs(trailingOnly=TRUE)
props_path = args[1]
real_cell_proportions = read.csv(props_path, row.names = 1)
real_cell_proportions = real_cell_proportions[,colnames(read.csv('Cache/none_reference.csv', row.names = 1))]
## Load packages
source('./Scripts/LoadSoftware.r')
source('./Scripts/MarkerSelection.r')

## Prepare variables
normalizations <- c('none', 'z_score', 'min_max', 'col_z_score', 'col_min_max', 'QN', 'lognorm')
longlist = matrix(NA, nrow =1, ncol = 5)
colnames(longlist) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')


## Perform all deconvolutions per normalization method and store predicted and actual values in 'longlist' variable
for (u in normalizations){
  
  df <- read.csv(paste0('Cache/',method, '_mixture.csv'), row.names = 1)
  ref <- read.csv(paste0('Cache/', method, '_reference.csv'), row.names = 1)
  ref <- ref[rownames(df),]
  

  #EMeth-Normal-----
  deconv = 'emeth_normal'
  maximum_nu = 0
  maximum_iter = 50
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  predicted_cell_proportions <- emeth(Y = as.matrix(df[rownames(ref),]),eta = c(rep(0, dim(as.matrix(df[rownames(ref),]))[2])), mu = as.matrix(ref), aber = FALSE, V = 'c',init = "default", family = "normal", nu = maximum_nu, maxiter = maximum_iter, verbose = TRUE)$rho
  colnames(predicted_cell_proportions) <- colnames(real_cell_proportions)
  rownames(predicted_cell_proportions) <- rownames(real_cell_proportions)
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  #EMeth-Laplace-----
  deconv = 'emeth_laplace'
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  predicted_cell_proportions <- emeth(Y = as.matrix(df[rownames(ref),]),eta = c(rep(0, dim(as.matrix(df[rownames(ref),]))[2])), mu = as.matrix(ref), aber = FALSE, V = 'c',init = "default", family = "laplace", nu = maximum_nu, maxiter = maximum_iter, verbose = TRUE)$rho
  colnames(predicted_cell_proportions) <- colnames(real_cell_proportions)
  rownames(predicted_cell_proportions) <- rownames(real_cell_proportions)
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  
  #EMeth-Binom-----
  deconv = 'emeth_binom'
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  predicted_cell_proportions <- emeth(Y = as.matrix(df[rownames(ref),]),eta = c(rep(0, dim(as.matrix(df[rownames(ref),]))[2])), mu = as.matrix(ref), aber = FALSE, V = 'b',init = "default", family = "normal", nu = maximum_nu, maxiter = maximum_iter, verbose = TRUE)$rho
  colnames(predicted_cell_proportions) <- colnames(real_cell_proportions)
  rownames(predicted_cell_proportions) <- rownames(real_cell_proportions)
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  
  #NNLS------
  deconv = 'nnls'
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  for (i in 1:dim(real_cell_proportions)[1]){
    predicted_cell_proportions[i,] <- coef(nnls(as.matrix(ref), df[rownames(ref),i]))
  }
  colnames(predicted_cell_proportions) <- colnames(real_cell_proportions)
  rownames(predicted_cell_proportions) <- rownames(real_cell_proportions)
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  
  #BVLS------
  deconv = 'bvls'
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  for (i in 1:dim(real_cell_proportions)[1]){
    predicted_cell_proportions[i,] <- coef(bvls(as.matrix(ref), df[rownames(ref),i], bl = c(rep(0, dim(ref)[2])), bu = c(rep(1, dim(ref)[2]))))
  }
  colnames(predicted_cell_proportions) <- colnames(real_cell_proportions)
  rownames(predicted_cell_proportions) <- rownames(real_cell_proportions)
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  #EpiDISH------
  deconv = 'epidish'
  predicted_cell_proportions <- matrix(NA, ncol = dim(real_cell_proportions)[2], nrow = dim(real_cell_proportions)[1])
  predicted_cell_proportions <- epidish(as.matrix(df), as.matrix(ref))$estF
  colnames(predicted_cell_proportions) <- colnames(real_cell_proportions)
  rownames(predicted_cell_proportions) <- rownames(real_cell_proportions)
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  l <- c()
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  
  #Minfi------
  deconv = 'minfi'
  predicted_cell_proportions <- projectCellType(Y = as.matrix(df[rownames(ref),]), coefCellType = as.matrix(ref))
  colnames(predicted_cell_proportions) <- colnames(real_cell_proportions)
  rownames(predicted_cell_proportions) <- rownames(real_cell_proportions)
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  
  #MethylResolver------
  deconv = 'methylresolver'
  methylMix <- as.data.frame(df[rownames(ref),])
  methylSig <- as.data.frame(ref)
  regressionFormula = as.formula(paste0("methylMix[,i] ~ ",paste(colnames(ref),sep="",collapse=" + ")))
  alpha <- 0.5
  preds <- matrix(NA, ncol = length(colnames(ref)), nrow = length(colnames(df)))
  colnames(preds) <- colnames(ref)
  rownames(preds) <- colnames(df)
  j = 1
  for (i in colnames(df)){
    deconvoluteSample <- robustbase::ltsReg(regressionFormula, data = methylSig, alpha = 0.5)
    preds[j,] <- deconvoluteSample$coefficients[2:length(deconvoluteSample$coefficients)]
    j <- j+1
  }
  predicted_cell_proportions <- preds
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  #Meth_atlas------
  deconv = 'meth_atlas'
  predicted_cell_proportions <- read.csv(paste(method, '_results_meth_atl.csv', sep = ''))
  rownames(predicted_cell_proportions) <- predicted_cell_proportions$X
  predicted_cell_proportions$X <- NULL
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  #OLS------
  deconv = 'ols'
  m <- matrix(NA, nrow = dim(df)[2], ncol  = dim(ref)[2])
  colnames(m) <- colnames(ref)
  rownames(m) <- colnames(df)
  cmd = 'lm(formula = df[rownames(ref),i] ~ 0'
  for(i in colnames(ref)){
  cmd = paste0(cmd, paste0(' + ref[,"',i, '"]'))
  }
  cmd = paste0(cmd, ')')
  for (i in 1:dim(df)[2]){
    model <- eval(parse(text = cmd))
    m[i,] <- unlist(lapply(model$coefficients, function (x) x/sum(model$coefficients)))
  }
  m[m<0] <- 0.00000001
  predicted_cell_proportions <- m
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  #Ridge------
  deconv = 'ridge'
  m <- matrix(NA, nrow = dim(df)[2], ncol  = dim(ref)[2])
  colnames(m) <- colnames(ref)
  rownames(m) <- colnames(df)
  for (i in 1:dim(df)[2]){
    mod <- glmnet(x = as.matrix(ref),y=df[rownames(ref), i], alpha = 0)
    mod <- glmnet(x = as.matrix(ref),y=df[rownames(ref), i], alpha = 0, lambda = mod$lambda[which(mod$dev.ratio == max(mod$dev.ratio))])
    m[i,] <- coef(mod)[2:length(coef(mod))]
  }
  m[m<0] <- 0.00000001
  predicted_cell_proportions <- m
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  #Elastic net------
  deconv = 'elastic_net'
  m <- matrix(NA, nrow = dim(df)[2], ncol  = dim(ref)[2])
  colnames(m) <- colnames(ref)
  rownames(m) <- colnames(df)
  for (i in 1:dim(df)[2]){
    mod <- glmnet(x = as.matrix(ref),y=df[rownames(ref), i], alpha = 0.5)
    mod <- glmnet(x = as.matrix(ref),y=df[rownames(ref), i], alpha = 0.5, lambda = mod$lambda[which(mod$dev.ratio == max(mod$dev.ratio))])
    m[i,] <- coef(mod)[2:length(coef(mod))]
  }
  m[m<0] <- 0.00000001
  predicted_cell_proportions <- m
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
  
  #Lasso------
  deconv = 'lasso'
  m <- matrix(NA, nrow = dim(df)[2], ncol  = dim(ref)[2])
  colnames(m) <- colnames(ref)
  rownames(m) <- colnames(df)
  for (i in 1:dim(df)[2]){
    mod <- glmnet(x = as.matrix(ref),y=df[rownames(ref), i], alpha = 1)
    mod <- glmnet(x = as.matrix(ref),y=df[rownames(ref), i], alpha = 1, lambda = mod$lambda[which(mod$dev.ratio == max(mod$dev.ratio))])
    m[i,] <- coef(mod)[2:length(coef(mod))]
  }
  m[m<0] <- 0.00000001
  predicted_cell_proportions <- m
  
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, predicted_cell_proportions[i,])}
  x <- matrix(NA, nrow = length(unlist(l)), ncol = 5)
  x[,1] <- unlist(l)
  x[,3] <- rep(colnames(real_cell_proportions), dim(predicted_cell_proportions)[1])
  x[,4] <- rep(method, dim(predicted_cell_proportions)[1])
  x[,5] <- rep(deconv, dim(predicted_cell_proportions)[1])
  
  colnames(x) <- c('Predicted value', 'True value','Celltype', 'Normalization method', 'Deconvolution method')
  l <- c()
  for(i in 1:dim(predicted_cell_proportions)[1]){
    l<- c(l, real_cell_proportions[i,])}
  x[,2] <- unlist(l)
  longlist <- rbind(longlist, x)
}

##Write results into 'Deconvolution_results.csv' file
write.csv(longlist[-1,], 'Deconvolution_results.csv')