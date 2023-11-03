## Strip algorithms of inherent normalization
projectCellType <- function(Y, coefCellType, contrastCellType = NULL,
                            nonnegative = TRUE, lessThanOne = FALSE) {
  if (is.null(contrastCellType)) {
    Xmat <- coefCellType
  } else {
    Xmat <- tcrossprod(coefCellType, contrastCellType)
  }
  
  nCol <- dim(Xmat)[2]
  if (nCol == 2) {
    Dmat <- crossprod(Xmat)
    mixCoef <- t(
      apply(Y, 2, function(x) solve(Dmat, crossprod(Xmat, x))))
    colnames(mixCoef) <- colnames(Xmat)
    return(mixCoef)
  } else {
    nSubj <- dim(Y)[2]
    
    mixCoef <- matrix(0, nSubj, nCol)
    rownames(mixCoef) <- colnames(Y)
    colnames(mixCoef) <- colnames(Xmat)
    
    if (nonnegative) {
      if (lessThanOne) {
        Amat <- cbind(rep(-1, nCol), diag(nCol))
        b0vec <- c(-1, rep(0, nCol))
      } else {
        Amat <- diag(nCol)
        b0vec <- rep(0, nCol)
      }
      for (i in seq_len(nSubj)) {
        obs <- which(!is.na(Y[,i]))
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve.QP(
          Dmat = Dmat,
          dvec = crossprod(Xmat[obs,], Y[obs,i]),
          Amat = Amat,
          bvec = b0vec)$sol
      }
    } else {
      for (i in seq_len(nSubj)) {
        obs <- which(!is.na(Y[,i]))
        Dmat <- crossprod(Xmat[obs,])
        mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
      }
    }
    mixCoef
  }
}


FeatureSelect <- function(RGSet.Cancer, do.BMIQ = F, norm.method = "Funnorm",AllGap = 0.3, nnGap = 0.1, whitelist = NULL, verbose = T, feature.number = 50 , probeFDR = 0.01) {
  
  #Load libraries and stuff
  
  require(plyr)
  require(dplyr)
  require(ggplot2)
  require(reshape2)
  require(Biobase)
  require(minfi)
  require(BiocGenerics)
  require(limma)
  
  #Load RGsets for infiltrating cell line stuff
  
  data("RGSets_infiltrate.RData")
  
  if(verbose){ message("Beginning combination and evaluation processes") }
  
  if(!class(RGSet.Cancer)=="RGChannelSet"){
    stop("Your set of cancer cell lines needs to be an RGChannelSet")
  }
  
  #Combine and normalise
  
  Combined.RGSet <- BiocGenerics::combine(RGSet.Cancer,MicSet2)
  
  
  if(norm.method == "Funnorm") {
    
    if(verbose){ message ("Functional normalisation in process")}
    Combined <- preprocessFunnorm(Combined.RGSet, bgCorr = F, dyeCorr = F, verbose = verbose)
    
  }
  
  
  if(norm.method == "SWAN") {
    
    Combined <- preprocessSWAN(Combined.RGSet)
    
  }
  
  
  if(norm.method == "Raw") {
    
    Combined <- preprocessRaw(Combined.RGSet)
    
  }
  
  
  #Continue
  Pheno.cancer <- data.frame(Sample = colnames(RGSet.Cancer))%>%
    mutate(Group = "Cancer")
  
  Pheno <- rbind(Pheno.cancer, Ref.pheno2)
  Beta <- getBeta(Combined)
  rm(Combined.RGSet, RGSet.Cancer, Combined)
  Beta <- Beta[complete.cases(Beta),]
  
  #This module executes BMIQ normalisation
  
  if(!is.null(whitelist)) { Beta <- Beta[rownames(Beta) %in% whitelist ,]}
  
  if(do.BMIQ) {
    message("Performing BMIQ normalisation")
    Beta <- HMkit.BMIQ(vals = Beta, nfit = 10000, plot = F)
  }
  
  #From here I can work on actually defining signatures
  
  design <- model.matrix(~factor(Pheno$Group))
  colnames(design) <- levels(factor(Pheno$Group))
  GroupsVec <- levels(factor(Pheno$Group))
  
  #This model sweeps through for one vs all comparisons
  
  if(verbose) { message ("Beginning one versus all sweeps")}
  
  ModelOutputs.list <- list()
  
  for(i in 1:length(GroupsVec)) {
    
    message(paste0("Initiating"," ",i))
    Vec1 <- as.character(Pheno$Group)
    Vec1 <- ifelse(Vec1 %in% GroupsVec[[i]], "One","Other")
    design <- model.matrix(~0+factor(Vec1))
    colnames(design) <- levels(factor(Vec1))
    
    Fit <- find.mvp(values= Beta, type = "beta", design, coef = 1, contrast.matrix = makeContrasts(One-Other,levels = design),classes = factor(Vec1), TSBH = T, alpha.TSBH = 0.05, measure.type = "median")
    
    message(paste0("Finishing"," ", i))
    ModelOutputs.list[[i]] <- Fit
    
  }
  
  
  #Filter by initial criteria here
  
  if(verbose) { message("Filtering probes for further application here")}
  
  Models2 <- lapply(ModelOutputs.list, function(x) Annotate.object(object = x$tt,use.dB = "median"))
  names(Models2) <- GroupsVec
  
  Filtered.tt <- lapply(Models2, function(x) filter(x, adj.P.Val < probeFDR, dB > AllGap | dB < -AllGap))
  
  lapply(Filtered.tt,dim)
  
  #Then we filter by nnGap - this ensures there is some degree of separation
  
  Medians.List <- list()
  
  for (i in 1:length(GroupsVec)) {
    
    A <- GroupsVec[[i]]
    B <- Beta[,Pheno$Group==A]
    Medians.List[[i]] <- rowMedians(B)
    message(i)
  }
  
  Medians.List <- do.call(cbind,Medians.List)
  
  DiffList.closest <- list()
  
  if(verbose) {message("Generating matrix of typewise medians")}
  Medians.List <- Medians.List[complete.cases(Medians.List),]
  for(i in 1:ncol(Medians.List)) {
    
    M <- Medians.List
    M1 <- M[,i]
    M2 <- M[,-i]
    df <- data.frame(MinDiff= M1 - rowMin(M2),MaxDiff= M1 - rowMax(M2))
    DiffList.closest[[i]] <- df
    message(i)
    
  }
  
  
  DiffList2 <- lapply(DiffList.closest, function(x) transform(x,ID=rownames(x))%>%filter(MaxDiff > nnGap | MinDiff < -nnGap))
  
  lapply(DiffList2, dim)
  
  ## Intersect then
  
  Concise.set <- list()
  
  
  for(i in 1:length(GroupsVec)) {
    
    TopTab <- Filtered.tt[[i]]
    Ref <- DiffList2[[i]]$ID
    
    Concise.set[[i]] <- TopTab%>%filter(ID %in% Ref)
    message(i)
    
  }
  
  if(verbose) { message("Beginning final step of feature selection")}
  
  Concise.set2 <- lapply(Concise.set, function(x) x%>%.[order(.$t),])
  
  #This block of code basically sets feature upper bound - we would like symmetric feature numbers for cell types
  
  Threshold.val <- lapply(Concise.set2, nrow)%>%unlist(.)%>%min(.)
  if(feature.number > Threshold.val) { feature.number <- Threshold.val}
  Halfway.point <- feature.number/2
  message(feature.number)
  Sig <- lapply(Concise.set2,function(x) x[c(1:Halfway.point,(nrow(x)-Halfway.point):nrow(x)),])
  SigIDs <- do.call(rbind,Sig)%>%.$ID%>%unique(.)
  
  if(verbose){message("Signature Features Identified")}
  
  
  #Finally, this module selects and summarises the feature set for a run
  
  
  Beta2 <- Beta[rownames(Beta) %in% SigIDs,]
  Beta2 <- 100*Beta2
  Beta2 <- t(Beta2)
  Beta2 <- data.frame(Beta2)
  Beta2$ID <- Pheno$Group
  D2.summary <- Beta2%>%group_by(ID)%>%summarise_each(funs(mean))
  
  colnames(D2.summary)[[1]] <- "NAME"
  D2.summary <- t(D2.summary)
  D2.summary2 <- D2.summary[2:nrow(D2.summary),]
  colnames(D2.summary2) <- as.character(D2.summary[1,])
  D3.summary <- cbind(data.frame(NAME=rownames(D2.summary2), data.frame(D2.summary2)))
  names(Concise.set2) <- GroupsVec
  return(list(TopTables = Concise.set2, SignatureMeans = D3.summary))
  if(verbose){message("Completed")}
  
}

Prep.CancerType <- function(Beta, Probes, fname) {
  
  Common <- intersect(Probes, rownames(Beta))
  Beta <- Beta[match(Common,rownames(Beta)),]
  NAME <- data.frame(NAME = Common)
  Beta <- cbind(NAME, data.frame(Beta))
  write.table(Beta, file = paste0(method, datast, '_Mixture.txt', sep =''),sep = "\t", row.names = FALSE, quote = F )
}

Featureselect.V2 <- function(CellLines.matrix, Heatmap = FALSE, nCores = 4, reps.resamp = 20, export = TRUE, sigName, Stroma.matrix, Phenotype.stroma, Unlog = TRUE) {
  
  require(caret)
  require(glmnet)
  require(foreach)
  require(NMF)
  require(doParallel)
  require(matrixStats)
  
  if(!is.null(ncol(CellLines.matrix))) {  Pheno1 <- c(rep("Cancer",ncol(CellLines.matrix)))
  Pheno2 <- c(as.character(Pheno1), as.character(Phenotype.stroma)) }
  
  else { Pheno2 <- as.character(Phenotype.stroma)}
  
  Mat2 <- cbind(CellLines.matrix, Stroma.matrix)
  
  message("Beginning elastic net selection procedure")
  message("63.2 bootstrapping to be used")
  
  #Create foreach cluster for parallelisation
  
  Features.CVparam<- trainControl(method = "boot632",number = reps.resamp,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)
  
  if(nCores > 1) {  registerDoParallel(makeCluster(nCores))
    message( "Parallelisation schema set up")}
  
  
  Model <- train(x = t(Mat2), y = factor(Pheno2), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=c(0.5,1),.lambda = seq(0,0.05,by=0.01)), metric = "Kappa")
  
  message("Retrieving Nonzero Coefficients")
  Nonzeros <-  coef(Model$finalModel, s = Model$bestTune$lambda)
  Nonzeros <- lapply(Nonzeros, function(x) data.frame(ID = rownames(x), Coef = as.numeric(x[,1])))
  Nonzeros <- lapply(Nonzeros, function(x) filter(x, !Coef == 0))
  Nonzeros <- do.call(rbind, Nonzeros)
  
  #Then I need to do the whole shebang of getting the features
  
  Mat3 <- Mat2[rownames(Mat2) %in% Nonzeros$ID,]
  
  if(Heatmap) { aheatmap(Mat3, annCol = Pheno2, labCol = NA, distfun = "euclidean", annColors = list(c("black","dodgerblue3","orange","red","cyan","firebrick3","grey","brown","dodgerblue2","yellow","purple","goldenrod","darkorange")), color = colorRampPalette(c("dodgerblue3","orange","red","firebrick3"))(200))
  }
  
  #Print number of selected probes
  
  print(nrow(Mat3))
  if(Unlog) { Mat3 <- 2^Mat3}
  
  
  #Then , export object
  if(export) {
    
    DF <- as.data.frame(t(Mat3))
    DF <- split(DF, factor(Pheno2))
    Collapsed <- lapply(DF, function(x) colMedians(data.matrix(x)))
    Collapsed <- data.frame(do.call(cbind, Collapsed))
    Collapsed <- cbind(data.frame(NAME = rownames(Mat3), stringsAsFactors = F),
                       Collapsed)
    
    fN <- paste0(sigName, "_Signature.txt")
    write.table(Collapsed, file = fN, sep = "\t", row.names = FALSE, quote = FALSE )
    
  }
  
  return(list(SignatureMatrix = Mat3))
  
}

FeatureSelect.V3 <- function(Heatmap = FALSE, nCores = 4, reps.resamp = 20, export = TRUE, sigName, Stroma.matrix, Phenotype.stroma, Unlog = FALSE, BetaScale = TRUE) {
  
  #This function uses pairwise glmnetting
  #This may be rather slow but whatever - it shoul do better on the basis of the nearest cell types
  
  
  require(caret)
  require(glmnet)
  require(foreach)
  require(NMF)
  require(doParallel)
  require(matrixStats)
  
  if(!is.null(ncol(CellLines.matrix))) {  Pheno1 <- c(rep("Cancer",ncol(CellLines.matrix)))
  Pheno2 <- c(as.character(Pheno1), as.character(Phenotype.stroma)) } else { Pheno2 <- as.character(Phenotype.stroma)}
  
  Mat2 <- cbind(CellLines.matrix, Stroma.matrix)
  
  message("Setting up for pairwise feature selection")
  message("Beginning elastic net selection procedure")
  message("63.2 bootstrapping to be used")
  
  #Create foreach cluster for parallelisation
  
  Features.CVparam<- trainControl(method = "boot632",number = reps.resamp,verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=TRUE)
  
  Pairs <- data.frame(t(combn(unique(Pheno2),2)), stringsAsFactors = F)
  Pairs <- filter(Pairs, !X1 == X2)
  
  #Then I do the pairwise fitting here
  
  FitList <- list()
  
  for(i in 1:nrow(Pairs)) {
    
    I1 <- Phenotype.stroma == Pairs[i,]$X1 | Phenotype.stroma == Pairs[i,]$X2
    M1 <- Mat2[,I1]
    P1 <- as.character(Pheno2[I1])
    
    Model <- train(x = t(M1), y = factor(P1), trControl = Features.CVparam, method = "glmnet" , tuneGrid = expand.grid(.alpha=c(0.5,1),.lambda = seq(0,0.05,by=0.01)), metric = "Kappa")
    
    Nonzeros <- coef(Model$finalModel, s = Model$bestTune$lambda)
    Nonzeros <- as.matrix(Nonzeros)
    Nonzeros <- data.frame(ID = rownames(Nonzeros), Coef = as.numeric(Nonzeros[,1]))
    Nonzeros <- filter(Nonzeros, !Coef == 0)
    FitList[[i]] <- Nonzeros
    
    message(paste0("pair",i," done of ", nrow(Pairs)))
    
    
  }
  
  #Then I want to unify the nonzeros and then move on
  
  Nonzeros <- do.call(rbind, FitList)
  Nonzeros <- filter(Nonzeros, !duplicated(ID))
  
  Mat3 <- Mat2[rownames(Mat2) %in% Nonzeros$ID,]
  
  if(Heatmap) { aheatmap(Mat3, annCol = Pheno2, labCol = NA, distfun = "euclidean", annColors = list(c("black","dodgerblue3","orange","red","cyan","firebrick3","grey","brown","dodgerblue2","yellow","purple","goldenrod","darkorange")), color = colorRampPalette(c("dodgerblue3","orange","red","firebrick3"))(200))
  }
  
  #Print number of selected probes
  
  if(BetaScale) { Mat3 <- 100 * Mat3}
  print(nrow(Mat3))
  if(Unlog) { Mat3 <- 2^Mat3}
  
  
  #Then , export object
  if(export) {
    
    DF <- as.data.frame(t(Mat3))
    DF <- split(DF, factor(Pheno2))
    Collapsed <- lapply(DF, function(x) colMedians(data.matrix(x)))
    Collapsed <- data.frame(do.call(cbind, Collapsed))
    Collapsed <- cbind(data.frame(NAME = rownames(Mat3), stringsAsFactors = F),
                       Collapsed)
    
    fN <- paste0(sigName, "_Signature.txt")
    write.table(Collapsed, file = fN, sep = "\t", row.names = FALSE, quote = FALSE )
    
  }
  
  return(list(SignatureMatrix = Mat3))
  
}

FeatureSelect.V4 <- function( CellLines.matrix = NULL, Heatmap = FALSE, export = TRUE, sigName, Stroma.matrix, Phenotype.stroma, deltaBeta, FDR , MaxDMRs = 1000) {
  
  #This function uses pairwise limma
  #This may be rather slow but whatever - it shoul do better on the basis of the nearest cell types
  
  require(dplyr)
  require(caret)
  require(glmnet)
  require(foreach)
  require(NMF)
  require(doParallel)
  require(matrixStats)
  require(limma)
  
  if(!is.null(ncol(CellLines.matrix))) {  Pheno1 <- c(rep("Cancer",ncol(CellLines.matrix)))
  Pheno2 <- c(as.character(Pheno1), as.character(Phenotype.stroma)) } else { Pheno2 <- as.character(Phenotype.stroma)}
  
  Mat2 <- cbind(CellLines.matrix, Stroma.matrix)
  
  message("Setting up for pairwise feature selection")
  
  
  #Work out the linear model fit here
  
  ContrastMatrix <- design.pairs(levels(factor(Pheno2)))
  Des <- model.matrix(~0 + Pheno2)
  colnames(Des) <- rownames(ContrastMatrix)
  Fit <- lmFit(Mat2, Des)%>%
    contrasts.fit(., ContrastMatrix)%>%
    eBayes(.)
  
  
  FitList <- list()
  for(i in 1:ncol(ContrastMatrix)) {
    
    FitList[[i]] <- topTable(Fit, coef = i, number = nrow(Mat2))%>%
      mutate(ID = rownames(.))%>%
      filter(adj.P.Val < FDR)
    
    message(paste0(i, " done"))
    
    
  }
  
  
  
  #Here goes the pairwise dB thing
  
  Transformed <- data.frame(t(Mat2))
  Split <- split(Transformed, Pheno2)
  Split <- lapply(Split, function(x) colMedians(data.matrix(x)))
  Split <- do.call(cbind, Split)
  rownames(Split) <- rownames(Mat2)
  
  #Then I need to actually annotate each one of these comparison topTables
  
  dbList <- list()
  message("Getting Delta Beta estimates")
  for(i in 1:ncol(ContrastMatrix)) {
    
    dB <- with(data.frame(Split), eval(parse(text = colnames(ContrastMatrix)[[i]])))
    dB <- data.frame(dB = dB, ID = rownames(Split))
    dbList[[i]] <- dB
    message(paste0 (i, " done"))
  }
  
  
  #Filter by thresholds
  
  
  dbList <- lapply(dbList, function(x) filter(x, abs(dB) > deltaBeta))
  for(i in 1:length(FitList)) {
    
    A1 <- FitList[[i]]
    A1 <- filter(A1 , ID %in% dbList[[i]]$ID)
    A1 <- A1%>%.[rev(order(.$t)),]
    if(nrow(A1) > MaxDMRs) { A1 <-  A1[1:MaxDMRs,]                   }
    FitList[[i]] <- A1
  }
  
  Nonzeros <- lapply(FitList, function(x) dplyr::select(x,ID))
  Nonzeros <- do.call(rbind, Nonzeros)
  Nonzeros <- filter(Nonzeros, !duplicated(ID))
  
  Mat3 <- Mat2[rownames(Mat2) %in% Nonzeros$ID,]
  
  if(Heatmap) { aheatmap(Mat3, annCol = Pheno2, labCol = NA, distfun = "euclidean", annColors = list(c("black","dodgerblue3","orange","red","cyan","firebrick3","grey","brown","dodgerblue2","yellow","purple","goldenrod","darkorange")), color = colorRampPalette(c("dodgerblue3","orange","red","firebrick3"))(200))
  }
  
  #Print number of selected probes
  nrow(Mat3)
  Mat3 <- 100 * Mat3
  
  #Then , export object
  if(export) {
    
    DF <- as.data.frame(t(Mat3))
    DF <- split(DF, factor(Pheno2))
    Collapsed <- lapply(DF, function(x) colMedians(data.matrix(x)))
    Collapsed <- data.frame(do.call(cbind, Collapsed))
    Collapsed <- cbind(data.frame(NAME = rownames(Mat3), stringsAsFactors = F),
                       Collapsed)
    
    fN <- paste0(sigName, "_Signature.txt")
    write.table(Collapsed, file = fN, sep = "\t", row.names = FALSE, quote = FALSE )
    
  }
  
  return(list(SignatureMatrix = Mat3))
  
}



#This function creates the pairs for the pairwise matrices
design.pairs <-
  function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
      for (j in (i+1):n) {
        k <- k+1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
      }
    design
  }
Annotate.object <- function(object, use.dB = "mean"){
  
  data("probe_450K_VCs_af.RData")
  data("probe_450k_features_V4.rda")
  data("Epic_450k_consensus.RData")
  data("450k_whitelist.RData")
  probeData <- probe.features
  
  
  #compute dB
  
  if(use.dB == "mean") { object$dB <- object[,9] - object[,8] ; message(paste0("Computing diff"," ", colnames(object)[[9]] ," - " ,colnames(object)[[8]]))}
  if(use.dB == "median") { object$dB <- object[,11] - object[,10] ; message(paste0("Computing diff"," ", colnames(object)[[11]] ," - " ,colnames(object)[[10]]))}
  
  #data(probeData)
  
  if(class(object)=="list") {
    
    object <- lapply(object, function(x) merge(x,probeData,by = "ID"))
    
  } else {
    
    object <- merge(object, probeData, by="ID")
    
  }
  
  return(object)
  
}


find.mvp<-function(values, type=c("beta","M"), design, coef, contrast.matrix,classes, measure.type,TSBH=T,alpha.TSBH=0.05,... ){
  
  
  #Create results list
  
  Results.List <- list()
  
  #This block of code generates an object called Fit using limma
  
  message("Fitting Model")
  if(is.null(contrast.matrix)){
    
    Fit<-lmFit(values,design)%>%
      eBayes(.,trend=T)
    
    
  } else {
    
    Fit<-lmFit(values,design)%>%
      contrasts.fit(.,contrast.matrix)%>%
      eBayes(.,trend=T)
    
  }
  
  # Add to list here
  
  Results.List$Fit <- Fit
  
  
  
  if(length(coef) == 1) {
    
    tt <- topTable(Fit, coef= coef, number = nrow(values))%>%transform(.,ID=rownames(.))
    
    if(type=="beta"){
      
      if(is.null(classes)){stop("computing delta beta requires a vector of groups")}
      message("Computing beta value summary statistics ")
      
      dB.mean <-sumClass(classes=classes,values=values,measure="mean")%>%
        data.frame(.)
      
      colnames(dB.mean) <- paste0("mean",".",colnames(dB.mean))
      
      dB.median <- sumClass(class=classes, values=values, measure="median")%>%
        data.frame(.)
      
      colnames(dB.median) <- paste0("median",".", colnames(dB.median))
      
      dB.mean$ID <- rownames(dB.mean) ; dB.median$ID <- rownames(dB.median)
      
      dB.tab <- merge(dB.mean,dB.median,by="ID")
      
      
      tt<-merge(tt,dB.tab,by="ID")
      
    }
    
    
    
    if(TSBH) {  tt <- TSBH.adjust(topTable=tt, alpha = alpha.TSBH) }
    Results.List$tt <- tt
    
  } else {
    
    
    ## This block of code is designed to work with multiple specified coefficients
    message("Multiple coefficients specified, preparing limma topTables")
    
    topTables.list <- list()
    
    
    for( i in 1:length(coef)) {
      
      tt <- topTable(Fit, coef= i,number = nrow(values))%>%transform(., ID=rownames(.))
      
      if(type=="beta") {
        
        if(is.null(classes)){stop("computing delta beta requires a vector of groups")}
        message("Computing beta value summary statistics ")
        
        cl.v <- classes[[i]]
        dB.mean <-sumClass(classes=cl.v,values=values,measure="mean")%>%
          data.frame(.)
        
        colnames(dB.mean) <- paste0("mean",colnames(dB.mean))
        
        dB.median <- sumClass(classes = classes[[i]], values=values, measure="median")%>%
          data.frame(.)
        
        colnames(dB.median) <- paste0("median", colnames(dB.median))
        
        dB.tab <- merge(dB.mean,dB.median,by=0)%>%transform(.,ID=rownames(.))
        
        tt <-merge(tt,dB.tab,by="ID")
        
      }
      
      topTables.list[[i]] <- tt
      message(i)
      
    }
    
    topTables.list <- lapply(topTables.list, function(x) transform(x, ID=rownames(x)))
    
    if(TSBH) {
      for(i in 1:length(topTables.list)) { topTables.list[[i]] <- TSBH.adjust(topTable=tt[[i]],alpha=alpha.TSBH)}
    }
    
    
    
    Results.List$tt <- topTables.list
    
    
  }
  
  
  
  return(Results.List)
}



# New VMP stuff



find.vmp<-function(values, type=c("beta","M"), design, coef, contrast.matrix,TSBH=T,alpha.TSBH=0.05,... ){
  
  if(type=="beta"){ values = 2^(log2(values)-log2(1-values)); message("Converting to M values for variance modelling")}
  
  Results.List <- list()
  
  #This block of code generates an object called Fit using limma
  
  if(is.null(contrast.matrix)){
    
    Fit<-varFit(values,design)%>%
      eBayes(.)
    
    
    
  } else {
    
    Fit<-varFit(values,design)%>%
      contrasts.fit(.,contrast.matrix)%>%
      eBayes(.,trend=T)
    
    
  }
  
  Results.List$Fit <- Fit
  
  
  if(length(coef) == 1 ) {
    
    tt <- topVar(Fit,number=nrow(values),coef=coef)%>%
      transform(.,ID=rownames(.))
    
    if(TSBH) { tt <- TSBH.adjust.v(tt, alpha = alpha.TSBH)}
    
    Results.List$tt <- tt
    
  } else {
    
    tt.list <- list()
    
    for ( i in 1:length(coef)) {
      
      tt.list[[i]] <- topVar(Fit, number= nrow(values), coef=coef[[i]])%>%
        transform(., ID=rownames(.))
      
    }
    
    
    if(TSBH) {
      
      for( i in 1:length(coef)) { tt.list[[i]] <- TSBH.adjust.v(topTable= tt.list[[i]],alpha=alpha.TSBH)}
      
    }
    
    Results.List$tt <- tt.list
  }
  
  
  return(Results.List)
  
}





# This is the class summarisation function
sumClass<-function(classes,measure=c("median","mean"),values,verbose=T){
  
  levs<-levels(factor(classes))
  output<-list()
  
  for(i in 1:length(levs)){
    
    class<-classes %in% levs[[i]]
    val<-values[,class]
    
    if(measure=="median"){output[[i]]<-rowMedians(val)}
    if(measure=="mean")  {output[[i]]<-rowMeans(val)}
  }
  
  df<-do.call(cbind,output)
  colnames(df)<-levs
  rownames(df)<-rownames(values)
  return(df)
  
}

# This is the TSBH adjustment function
TSBH.adjust<-function(topTable,alpha){
  
  Fit<-topTable
  qvals<-mt.rawp2adjp(rawp=Fit$P.Value,proc="TSBH",alpha)
  q<-qvals$adjp[,2]
  index<-qvals$index
  Fit$adj.P.Val<-q[order(index)]
  return(Fit)
}

TSBH.adjust.v <- function(topTable,alpha){
  
  Fit<-topTable
  qvals<-mt.rawp2adjp(rawp=Fit$P.Value,proc="TSBH",alpha)
  q<-qvals$adjp[,2]
  index<-qvals$index
  Fit$Adj.P.Value<-q[order(index)]
  return(Fit)
}


#This implements the fDMR function

fDMR<-function(signif.only=F,mvp.tab,mvp.fdr,n.mvp,dmr.fdr,mvp.dB,dmr.dB,beta,dB.measure,classes,exclude.regions,os.p=F,resolution=c("high","low")){
  
  #Initial analysis...
  
  #Here, I aim to reduce stuff to candidate sigDMRs based on cutoff criteria
  
  #Annotate first with feature group
  message("annotating with feature groups")
  
  if(resolution=="low") {
    
    mvp.tab<-mvp.tab%>%
      mutate(feature2=as.character(feature),feature2=ifelse(feature2 %in% c("TSS1500","TSS200","1stExon"),
                                                            "promoter",ifelse(feature2 %in% "Body",
                                                                              "Body",feature2)),
             code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(adj.P.Val))
    
  } else {
    
    
    mvp.tab <- mvp.tab%>%mutate(feature2=as.character(feature), code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(adj.P.Val))
    
  }
  
  
  
  if(!is.null(exclude.regions)){
    
    mvp.tab<-filter(mvp.tab,!feature2 %in% exclude.regions)
    
  }
  
  
  mvp.cutoff<-mvp.tab%>%filter(adj.P.Val < mvp.fdr, dB > mvp.dB | dB < -mvp.dB)
  
  #Then carry out fDMR filtration to select dmrs that meet criteria for candidates
  
  mvp.count<-dplyr::count(mvp.cutoff,code)%>%
    filter(n > (n.mvp-1))
  
  
  #Then we reduce the original table to candidate dmrs for dmr calling and stuff
  message("setting up probe table for candidate dmrs")
  if(signif.only==FALSE){mvp.tab2<-mvp.tab%>%filter(code %in% mvp.count$code)}else
  {mvp.tab2<-mvp.cutoff%>%filter(code %in% mvp.count$code)}
  
  mvp.tab2<-mvp.tab2[order(mvp.tab2$code),]
  beta.tab2<-beta[match(mvp.tab2$ID,rownames(beta)),]
  
  
  if(os.p==T){
    
    #One sided P value conversion here...
    
    message("Converting p values to 1 sided values for p.value combination")
    P.os<-mapply(function(t,p) p.ts.os(t,p),mvp.tab2$t,as.character(mvp.tab2$P.Value))
    mvp.tab2$P.os<-P.os
  }  else {
    
    message("parameter choice is to not conver to one-sided p values")
    mvp.tab2$P.os<-as.numeric(as.character(mvp.tab2$P.Value))
    
  }
  
  #Then do the whole p value combination shenanigans.
  
  message("SLK correction")
  dmr.beta <- split(data.frame(beta.tab2), factor(mvp.tab2$code))
  corel <- lapply(dmr.beta, function(x) cor(t(x)))
  weights <- lapply(corel, function(x) 1/apply(x^2,1,sum))
  
  if(os.p==T){
    dmr.ind.p <- split(mpfr(mvp.tab2$P.os), factor(mvp.tab2$code))
    dmr.ind.p<-lapply(dmr.ind.p,function(x) asNumeric(x))
    dmr.qp <- lapply(dmr.ind.p, qnorm)
    dmr.qp.w <- mapply("*", dmr.qp, weights)
  } else {
    
    dmr.ind.p<- split(mvp.tab2$P.os, factor(mvp.tab2$code))
    dmr.qp <- lapply(dmr.ind.p, qnorm)
    dmr.qp.w <- mapply("*", dmr.qp, weights)
    
  }
  
  #Then generate p value distributions
  
  if(class(dmr.qp.w) == "matrix")
  {
    dmr.stat <- sum(dmr.qp.w)
  }else
  {
    dmr.stat <- lapply(dmr.qp.w, sum)
  }
  
  message("calculating DMR p values")
  dmr.sd <- lapply(weights, function(x) sqrt(sum(x^2)))
  dmr.p <- mapply(function(x,y) pnorm(x,0, sd=y), dmr.stat, dmr.sd)
  if(os.p==T) {dmr.p <- lapply(as.character(dmr.p), function(x) p.os.ts(x))}
  dmr.p <- unlist(dmr.p)
  
  #Once this is done, I can go on and build up a table
  
  code.p<-unique(factor(mvp.tab2$code))
  dmr.summary<-data.frame(code=code.p,dmr.p=dmr.p,dmr.FDR=p.adjust(dmr.p,method="BH"))
  
  message("Exporting results tables")
  mvp.summary<-merge(dmr.summary,mvp.tab2,by="code")
  
  #Then further annotate results table.
  
  dmr.summary<- filter(dmr.summary,dmr.FDR < dmr.fdr)
  mvp.summary<- filter(mvp.summary,code %in% dmr.summary$code, adj.P.Val < mvp.fdr, dB > mvp.dB | dB < -mvp.dB )
  
  dB.summary<-mvp.summary%>%
    group_by(code)%>%
    summarise(Mean.dB=mean(dB),Median.dB=median(dB))
  
  dmr.summary<-merge(dmr.summary,dB.summary,by="code",all.x=T)
  if(dB.measure=="mean"){dmr.summary<-filter(dmr.summary,Mean.dB > dmr.dB | Mean.dB < -dmr.dB)}
  if(dB.measure=="median"){dmr.summary<-filter(dmr.summary,Median.dB > dmr.dB | Median.dB < -dmr.dB)}
  
  
  mvp.summary<- filter(mvp.summary,code %in% dmr.summary$code, adj.P.Val < mvp.fdr, dB > mvp.dB | dB < -mvp.dB )
  
  
  return(list(Probes=mvp.summary,dmrs=dmr.summary))
  
}

#Coxfit.MVP - this function basically estimates cox model statistics for matrices of beta or m values.

coxfit.mvp<-function(time,status,strata.f,values,weighted=FALSE,verbose=T){
  
  values<-t(values)
  sData<-Surv(time,status)
  
  
  if(weighted==FALSE){
    
    results.list<-list()
    
    
    for(i in 1:ncol(values)){tryCatch({
      
      if(is.null(strata.f)){ cph.obj<-coxph(sData~values[,i])
      }  else {
        
        cph.obj<-coxph(sData~values[,i]+strata(strata.f),singular.ok=T)
        
      }
      
      results.list[[i]]<-coef(summary(cph.obj))
      rownames(results.list[[i]])<-colnames(values)[[i]]
      if(verbose==T){print(i)}
    },error=function(cond){})
      
    }
    
    
    l2<-do.call(rbind,results.list)%>%
      data.frame(.)
    
    colnames(l2)<-c("Coef","HR","se.Coef","Z","p.value")
    l2$ID<-rownames(l2)
    return(l2)
    
  } else {
    
    
    p.values<-c(rep(0,ncol(values)))
    coefs<-c(rep(0,ncol(values)))
    
    if(is.null(strata)){
      
      for(i in 1: ncol(values)) {
        
        fit<-coxphw(sData~values[,i])
        p.values[[i]]<-fit$prob
        coefs[[i]]<-fit$coef
      }
    } else {
      
      for(i in 1: ncol(values)){
        fit<-coxphw(sData~values[,i]+strata.f)
        p.values[[i]]<-fit$prob
        coefs[[i]]<-fit$coef
      }
    }
    
    ID.vec<-colnames(values)
    df<-data.frame(ID=ID.vec,coef=coefs,p.value=p.values)%>%
      mutate(HR=exp(coef))
    return(df)
    
  }
}

#Develop a version of fDMR for survival

surv.fDMR<-function(beta,mvp.cox,alpha=0.05,TSBH=T,mvp.fdr=0.05,dmr.fdr,p.os,signif.only=F,n.mvp=3,exclude.regions,resolution=c("high","low")){
  
  #This first bit performs TSBH correction if required
  
  if(TSBH){
    
    Fit<-mvp.cox
    qvals<-mt.rawp2adjp(rawp=Fit$p.value,proc="TSBH",alpha)
    q<-qvals$adjp[,2]
    index<-qvals$index
    Fit$adj.P.Val<-q[order(index)]
    
  } else {
    
    Fit$adj.P.Val<- p.adjust(Fit$adj.P.Val,method="BH")
    
  }
  
  mvp.tab<-Fit
  rm(Fit)
  
  #From the next bit we begin to pass variables into a P.value combination function after filtering.
  
  #Annotation and Filtration module
  
  message("annotating with feature groups")
  mvp.tab$ID<-rownames(mvp.tab)
  mvp.tab<-merge(mvp.tab,probe.features,by.x="ID",by.y=0,sort=F,all.x=T)
  
  if(resolution=="low") {
    
    mvp.tab<-mvp.tab%>%
      mutate(feature2=as.character(feature),feature2=ifelse(feature2 %in% c("TSS1500","TSS200","1stExon"),
                                                            "promoter",ifelse(feature2 %in% "Body",
                                                                              "Body",feature2)),
             code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(adj.P.Val))
    
  } else {
    
    
    mvp.tab <- mvp.tab%>%mutate(feature2=as.character(feature), code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(adj.P.Val))
    
  }
  
  
  
  if(!is.null(exclude.regions)){
    
    mvp.tab<-filter(mvp.tab,!feature2 %in% exclude.regions)
    
  }
  
  
  
  #The next bit filters by FDR and number of mvp cutoffs
  
  mvp.cutoff<-mvp.tab%>%filter(adj.P.Val < mvp.fdr)
  mvp.count<-dplyr::count(mvp.cutoff,code)%>%
    filter(n > (n.mvp-1))
  
  #Then we reduce the original table to candidate dmrs for dmr calling and stuff
  
  message("setting up probe table for candidate dmrs")
  if(signif.only==FALSE){mvp.tab2<-mvp.tab%>%filter(code %in% mvp.count$code)}else
  {mvp.tab2<-mvp.cutoff%>%filter(code %in% mvp.count$code)}
  
  mvp.tab2<-mvp.tab2[order(mvp.tab2$code),]
  beta.tab2<-beta[match(mvp.tab2$ID,rownames(beta)),]
  
  
  
  #P value conversion to one sided if needed
  
  if(os.p==T){
    
    #One sided P value conversion here...
    
    message("Converting p values to 1 sided values for p.value combination")
    P.os<-mapply(function(t,p) p.ts.os(t,p),mvp.tab2$t,as.character(mvp.tab2$P.Value))
    mvp.tab2$P.os<-P.os
  }  else {
    
    message("parameter choice is to not convert to one-sided p values")
    mvp.tab2$P.os<-as.numeric(as.character(mvp.tab2$P.Value))
    
  }
  
  #Correct p values for combination purposes
  message("Stouffer Liptak Correction")
  dmr.beta <- split(data.frame(beta.tab2), factor(mvp.tab2$code))
  corel <- lapply(dmr.beta, function(x) cor(t(x)))
  weights <- lapply(corel, function(x) 1/apply(x^2,1,sum))
  
  if(os.p==T){
    dmr.ind.p <- split(mpfr(mvp.tab2$P.os), factor(mvp.tab2$code))
    dmr.ind.p<-lapply(dmr.ind.p,function(x) asNumeric(x))
    dmr.qp <- lapply(dmr.ind.p, qnorm)
    dmr.qp.w <- mapply("*", dmr.qp, weights)
  } else {
    
    dmr.ind.p<- split(mvp.tab2$P.os, factor(mvp.tab2$code))
    dmr.qp <- lapply(dmr.ind.p, qnorm)
    dmr.qp.w <- mapply("*", dmr.qp, weights)
    
  }
  
  #Then generate p value distributions
  
  if(class(dmr.qp.w) == "matrix")
  {
    dmr.stat <- sum(dmr.qp.w)
  }else
  {
    dmr.stat <- lapply(dmr.qp.w, sum)
  }
  
  message("calculating DMR p values")
  dmr.sd <- lapply(weights, function(x) sqrt(sum(x^2)))
  dmr.p <- mapply(function(x,y) pnorm(x,0, sd=y), dmr.stat, dmr.sd)
  if(os.p==T) {dmr.p <- lapply(as.character(dmr.p), function(x) p.os.ts(x))}
  dmr.p <- unlist(dmr.p)
  
  #Once this is done, I can go on and build up a table
  
  message("Generating output objects")
  code.p<-unique(factor(mvp.tab2$code))
  dmr.summary<-data.frame(code=code.p,dmr.p=dmr.p,dmr.FDR=p.adjust(dmr.p,method="BH"))
  mvp.summary<-merge(dmr.summary,mvp.tab2,by="code")
  
  return(list(Probes=mvp.summary,dmrs=dmr.summary))
}


#Function to handle mpfr output

o.mpfr<-function(n) capture.output(n)[2]%>%substr(.,5,nchar(.))

#Functions for 1-sided / 2-sided p value conversions

p.ts.os<-function(t,p){
  p<-mpfr(p,base=10)
  if( t < 0 ) {px= 0.5* (1-p)} else
  { px = p*0.5}
  
  px<-asNumeric(px)
  return(as.character(px))
}

p.os.ts<-function(p){
  p<-mpfr(p)
  if( p <= 1/2) { p = 2*p} else { p = 2(1-p)}
  p<-asNumeric(p)
  return(p)
}


# fVMR for coordinate variability calling


fVMR<-function(signif.only=F,vmp.tab,vmp.fdr,n.vmp,vmr.fdr,Mtab,exclude.regions,os.p=F,resolution=c("high","low")){
  
  #Initial analysis...
  
  #Here, I aim to reduce stuff to candidate sigvmrs based on cutoff criteria
  
  #Annotate first with feature group
  message("annotating with feature groups")
  if(!"ID" %in% colnames(vmp.tab) ){vmp.tab$ID<-rownames(vmp.tab)}
  
  
  if(resolution=="low") {
    
    vmp.tab<-vmp.tab%>%
      mutate(feature2=as.character(feature),feature2=ifelse(feature2 %in% c("TSS1500","TSS200","1stExon"),
                                                            "promoter",ifelse(feature2 %in% "Body",
                                                                              "Body",feature2)),
             code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(Adj.P.Value))
    
  } else {
    
    
    vmp.tab <- vmp.tab%>%mutate(feature2=as.character(feature), code=paste0(gene,".",feature2),PVal=as.character(P.Value),FDR=as.character(Adj.P.Value))%>%
      dplyr::filter(!feature2 %in% exclude.regions)
    
  }
  
  
  
  
  
  
  if(!is.null(exclude.regions)){
    
    vmp.tab<-filter(vmp.tab,!feature2 %in% exclude.regions)
    
  }
  
  
  vmp.cutoff<-vmp.tab%>%filter(Adj.P.Value < vmp.fdr)
  
  #Then carry out fvmr filtration to select vmrs that meet criteria for candidates
  
  vmp.count<-dplyr::count(vmp.cutoff,code)%>%
    filter(n > (n.vmp-1))
  
  
  #Then we reduce the original table to candidate vmrs for vmr calling and stuff
  message("setting up probe table for candidate vmrs")
  if(signif.only==FALSE){vmp.tab2<-vmp.tab%>%filter(code %in% vmp.count$code)}else
  {vmp.tab2<-vmp.cutoff%>%filter(code %in% vmp.count$code)}
  
  vmp.tab2<-vmp.tab2[order(vmp.tab2$code),]
  Mtab.tab2<-Mtab[match(vmp.tab2$ID,rownames(Mtab)),]
  
  
  if(os.p==T){
    
    #One sided P value conversion here...
    
    message("Converting p values to 1 sided values for p.value combination")
    P.os<-mapply(function(t,p) p.ts.os(t,p),vmp.tab2$t,as.character(vmp.tab2$P.Value))
    vmp.tab2$P.os<-P.os
  }  else {
    
    message("parameter choice is to not convert to one-sided p values")
    vmp.tab2$P.os<-as.numeric(as.character(vmp.tab2$P.Value))
    
  }
  
  #Then do the whole p value combination shenanigans.
  
  message("Stouffer Liptak correction")
  vmr.Mtab <- split(data.frame(Mtab.tab2), factor(vmp.tab2$code))
  corel <- lapply(vmr.Mtab, function(x) cor(t(x)))
  weights <- lapply(corel, function(x) 1/apply(x^2,1,sum))
  
  if(os.p==T){
    vmr.ind.p <- split(mpfr(vmp.tab2$P.os), factor(vmp.tab2$code))
    vmr.ind.p<-lapply(vmr.ind.p,function(x) asNumeric(x))
    vmr.qp <- lapply(vmr.ind.p, qnorm)
    vmr.qp.w <- mapply("*", vmr.qp, weights)
  } else {
    
    vmp.tab2$P.os <- as.numeric(vmp.tab2$P.os)
    vmr.ind.p<- split(vmp.tab2$P.os, factor(vmp.tab2$code))
    vmr.qp <- lapply(vmr.ind.p, qnorm)
    vmr.qp.w <- mapply("*", vmr.qp, weights)
    
  }
  
  #Then generate p value distributions
  
  if(class(vmr.qp.w) == "matrix")
  {
    vmr.stat <- sum(vmr.qp.w)
  }else
  {
    vmr.stat <- lapply(vmr.qp.w, sum)
  }
  
  message("calculating vmr p values")
  vmr.sd <- lapply(weights, function(x) sqrt(sum(x^2)))
  vmr.p <- mapply(function(x,y) pnorm(x,0, sd=y), vmr.stat, vmr.sd)
  if(os.p==T) {vmr.p <- lapply(as.character(vmr.p), function(x) p.os.ts(x))}
  vmr.p <- unlist(vmr.p)
  
  #Once this is done, I can go on and build up a table
  
  code.p<-unique(factor(vmp.tab2$code))
  vmr.summary<-data.frame(code=code.p,vmr.p=vmr.p,vmr.FDR=p.adjust(vmr.p,method="BH"))
  
  message("Exporting results tables")
  vmp.summary<-merge(vmr.summary,vmp.tab2,by="code")
  
  #Then further annotate results table.
  
  vmr.summary<- filter(vmr.summary,vmr.FDR < vmr.fdr)
  vmp.summary<- filter(vmp.summary,code %in% vmr.summary$code, Adj.P.Value < vmp.fdr)
  
  var.summary<-vmp.summary%>%
    group_by(code)%>%
    summarise(Mean.logFC=mean(LogVarRatio),Median.logFC=median(LogVarRatio))
  
  vmr.summary<-merge(vmr.summary,var.summary,by="code",all.x=T)
  
  
  vmp.summary<- filter(vmp.summary,code %in% vmr.summary$code, Adj.P.Value < vmp.fdr)
  
  
  return(list(Probes=vmp.summary,vmrs=vmr.summary))
  
}


##Load in libraries
library(glmnet)
library(bvls)
library(preprocessCore)
library(quadprog)
library(stringr)
library(nnls)
library(EpiDISH)
library(minfi)
library(EMeth)
library(Matrix)
library(MethylResolver)
library(TOAST)
library(corpcor)
library(writexl)
library(limma)

