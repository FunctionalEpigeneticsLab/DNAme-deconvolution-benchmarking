## Load in data
args = commandArgs(trailingOnly=TRUE)
ref_path = args[1]
val_path = args[2]
ref = read.csv(ref_path, row.names = 1)
rows = rownames(ref)
cols = colnames(ref)
ref = apply(ref, 2, as.numeric)
rownames(ref) = rows
colnames(ref) = cols
val = read.csv(val_path, row.names = 1)
rows = rownames(val)
cols = colnames(val)
val = apply(val, 2, as.numeric)
rownames(val) = rows
colnames(val) = cols
props = read.csv(args[3], row.names = 1)

## Load packages
source('./Scripts/LoadSoftware.r')
source('./Scripts/MarkerSelection.r')

## Selection of marker CpGs
celltypes = colnames(ref)
for(i in colnames(props)){
    celltypes[grepl(i, celltypes)] = i
}

output_selection = identifyCpgs(ref, celltypes, args[4], args[5])
mcpgs = output_selection$cpgs


## Production of normalized reference matrices
dataset = ref
normalizations <- c('none', 'z_score', 'min_max', 'col_z_score', 'col_min_max', 'QN', 'lognorm')
for (u in normalizations){
  if(u == 'none'){
    #None
    method <- 'none'
    normaliz <- function(x){
      return(x)
    }
    
  }else if (u == 'z_score'){
    #Z-score
    method = 'z_score'
    normaliz <- function(x){
      return((x-mean(na.omit(x)))/sd(na.omit(x)))
    }
    
  } else if (u == 'min_max'){
    #Min-max
    method <- 'min_max'
    normaliz <- function(x){
      return((x-min(na.omit(x)))/(max(na.omit(x))-min(na.omit(x))))
    }
    
  } else if (u == 'col_z_score'){
    #Column Z-score
    method <- 'col_z_score'
    normaliz <- function(x){
      cols <-function(y){
        return((y-mean(na.omit(y)))/sd(na.omit(y)))
        
      }
      apply(x, 2, cols)
      
    }
    
  } else if (u == 'col_min_max'){
    #Column min-max
    method <- 'col_min_max'
    
    normaliz <- function(x){
      cols <-function(y){
        return((y-min(na.omit(y)))/(max(na.omit(y)-min(na.omit(y)))))
        
      }
      apply(x, 2, cols)
      
    }
    
  } else if (u == 'QN'){
    #Quantile normalization
    
    method <- 'QN'
    normaliz <- function(x){
      x <- na.omit(x)
      rows <- rownames(x)
      cols <- colnames(x)
      x <- normalize.quantiles(x)
      rownames(x) <- rows
      colnames(x) <- cols
      return(x)
    }
    #
  }  else if (u == 'lognorm'){
    method <- 'lognorm'
    normaliz <- function(x){
      x <- log(na.omit(x))
      x[abs(x)==Inf] <- NA
      return(na.omit(x))
      
    }
  }
  tmp_dataset <- normaliz(dataset)
  ref <- matrix(NA, ncol = length(unique(celltypes)), nrow = length(mcpgs))
  colnames(ref) <- unique(celltypes)
  rownames(ref) <- mcpgs
  for(i in colnames(ref)){
    ref[,i] <- apply(tmp_dataset[mcpgs, celltypes == i], 1, mean)
  }
  write.csv(ref, paste0('Cache/', method, '_reference.csv'))
}


## Production of normalized validation matrices
dataset = val
normalizations <- c('none', 'z_score', 'min_max', 'col_z_score', 'col_min_max', 'QN', 'lognorm')
for (u in normalizations){
  if(u == 'none'){
    #None
    method <- 'none'
    normaliz <- function(x){
      return(x)
    }
    
  }else if (u == 'z_score'){
    #Z-score
    method = 'z_score'
    normaliz <- function(x){
      return((x-mean(na.omit(x)))/sd(na.omit(x)))
    }
    
  } else if (u == 'min_max'){
    #Min-max
    method <- 'min_max'
    normaliz <- function(x){
      return((x-min(na.omit(x)))/(max(na.omit(x))-min(na.omit(x))))
    }
    
  } else if (u == 'col_z_score'){
    #Column Z-score
    method <- 'col_z_score'
    normaliz <- function(x){
      cols <-function(y){
        return((y-mean(na.omit(y)))/sd(na.omit(y)))
        
      }
      apply(x, 2, cols)
      
    }
    
  } else if (u == 'col_min_max'){
    #Column min-max
    method <- 'col_min_max'
    
    normaliz <- function(x){
      cols <-function(y){
        return((y-min(na.omit(y)))/(max(na.omit(y)-min(na.omit(y)))))
        
      }
      apply(x, 2, cols)
      
    }
    
  } else if (u == 'QN'){
    #Quantile normalization
    
    method <- 'QN'
    normaliz <- function(x){
      x <- na.omit(x)
      rows <- rownames(x)
      cols <- colnames(x)
      x <- normalize.quantiles(x)
      rownames(x) <- rows
      colnames(x) <- cols
      return(x)
    }
    #
  }  else if (u == 'lognorm'){
    method <- 'lognorm'
    normaliz <- function(x){
      x <- log(na.omit(x))
      x[abs(x)==Inf] <- NA
      return(na.omit(x))
      
    }
  }
  tmp_dataset <- normaliz(dataset)
  
  
  write.csv(tmp_dataset, paste0('Cache/', method, '_reference.csv'))
}

