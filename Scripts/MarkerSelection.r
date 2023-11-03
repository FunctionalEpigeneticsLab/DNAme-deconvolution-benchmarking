## df is a data frame of methylation values, with rows representing loci and columns representing samples.
## celltypes is a list representing the respective celltype each sample corresponds to.
## p_val_cutoff is a float representing the P value cutoff.
## n is an integer representing the preferred number of loci per cell type.
identifyCpgs <- function(df, celltypes, p_val_cutoff, n = 100){
    
    df <- na.omit(df)
    rows <- rownames(df)
    cols <- colnames(df)
    df <- apply(df, 2, as.numeric)
    rownames(df) <- rows
    colnames(df) <- cols
    df <- na.omit(df)
    colnames(df) <- cols
    pvals <- matrix(NA, ncol = length(unique(celltypes)), nrow = dim(df)[1])
    abs_diffs <- matrix(NA, ncol = length(unique(celltypes)), nrow = dim(df)[1])
    vars <- matrix(NA, ncol = length(unique(celltypes)), nrow = dim(df)[1])
	logfch <- matrix(NA, ncol = length(unique(celltypes)), nrow = dim(df)[1])
    colnames(pvals) <- unique(celltypes)
    rownames(pvals) <- rownames(df)
	colnames(logfch) <- unique(celltypes)
    rownames(logfch) <- rownames(df)
    colnames(abs_diffs) <- unique(celltypes)
    rownames(abs_diffs) <- rownames(df)
    colnames(vars) <- unique(celltypes)
    rownames(vars) <- rownames(df)
    pval_f <- function(x){max(summary(lm(unlist(x)~tmp))$coefficients[,4])}
    logfch_f <- function(x){log2(mean(x[tmp == 'A'])) - log2(mean(x[tmp !='A']))}
    abs_diffs_r <-function(x){
      if(mean(x[tmp == 'A']) < mean(x[tmp != 'A'])){
        min(x[tmp!='A'])- mean(x[tmp == 'A']) 
      } else {
        mean(x[tmp == 'A']) - max(x[tmp!='A'])
      }}

    for(i in unique(celltypes)){
      tmp <- celltypes
      for(j in 1:length(tmp)){
        if (tmp[j] == i){
          tmp[j] <- 'A'
        }}
      pvals[,i] <- apply(df, 1, pval_f)
      abs_diffs[,i] <- apply(df, 1, abs_diffs_r)
      logfch[,i] <- apply(df, 1, logfch_f)


    }
   cols <- colnames(pvals)
   rows <- rownames(pvals)
   pvals <- apply(pvals, 2, p.adjust, method = 'hochberg')
   colnames(pvals) <- cols
   rownames(pvals) <- rows
    cpglist <- c()
    for (i in unique(celltypes)){
      tmp <- pvals[,i]
      cpglist <- c(cpglist, names(tmp[tmp<p_val_cutoff]))
    }
    cpglist <- na.omit(unique(cpglist))
    new_cpglist <- c()
    for (i in unique(celltypes)){
      tmp <- abs_diffs[cpglist,i]
      new_cpglist <- c(new_cpglist, names(tmp[order(tmp, decreasing = T)][1:n]))
    }

    nms <-new_cpglist
    return(list('cpgs'=nms, 'pvals'=pvals,'logfch'=logfch ,'abs_diffs'=abs_diffs))


}