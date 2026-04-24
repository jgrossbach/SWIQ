#' function that prepares genotype and phenotype data for QTL-mapping with random forest
#'
#' @param genotype a binary matrix with strain (rows) X marker (columns) entries. Do not include markers with MAF=1.
#' @param maxNAs amount of samples that are allowed to have missing data for a given marker. If this limit is exceeded, the marker is excluded from the genotype.
#' @param propVar proportion of genotypical variance that is included in the population-structure covariates
#' @param phenotype A vector or matrix with numerical values. A matrix may contain NAs. 
#' @param sampleInfo A vector with integers indicating to which strain (row in the genotype) the measurement(s) correspond.
#' @param scale logical; Should centering and scaling be performed?
#'
#' @return List-object containing the genotype, the groupings of the markers, the mapping-covariates and the phenotype.
#' @export
#'
#' @examples
preMap <- function(genotype, maxNAs=floor(nrow(genotype)*0.5), propVar=0.75, phenotype, sampleInfo, scale=T){

  ###check input
  if(any(!is.matrix(genotype),!is.integer(genotype),genotype!=0&genotype!=1,na.rm = T)){
    stop("Genotype has to be a matrix containing integers (0/1)")
  }
  if(length(unique(sampleInfo))<nrow(genotype)){
    stop("Don't include genotype-data for strains for which no trait data is supplied. Details can be found in the manual.")
  }
  if(!is.integer(sampleInfo)){
    stop("sampleInfo has to be an integer vector.")
  }
  if(any(apply(genotype,2,FUN=function(m){length(unique(m))==1}))){
    stop("Don't include markers with MAF=1.")
  }
  if(any(!is.numeric(maxNAs),maxNAs>length(sampleInfo),maxNAs<0)){
    stop("maxNAs has to be a single integer>=0 and <= to the number of samples")
  }
  if(!is.numeric(propVar)){
    stop("propVar has to be numeric.")
  }
  if(propVar>=1|propVar<0){
    stop("propVar has to be at least zero and has to be smaller than one.")
  }
  if(any(!is.numeric(phenotype),(!is.vector(phenotype))&(!is.matrix(phenotype)))){
    stop("phenotype has to be numeric and either a vector or a matrix")
  }
  if(length(sampleInfo)!=ifelse(is.matrix(phenotype),ncol(phenotype),length(phenotype))){
    stop("sampleInfo has to be as long as there are samples")
  }
  if(ifelse(is.vector(phenotype),any(is.na(phenotype)),any(apply(is.na(phenotype),2,all)))){
    stop("Don't include samples without any trait data.")
  }
  
  
  ###reduce genotype while retaining information
  
  #find markers with too many NAs
  nNA <- apply(genotype,2,FUN=function(x){sum(is.na(x))})
  
  #group redundant markers into genotype groups, leave out markers with too many NAs
  genotype2group <- rep(NA,ncol(genotype))
  reversedGenotype <- abs(genotype-1)
  
  redGeno <- NULL
  nGroups <- 0
  
  for(i in 1:ncol(genotype)){
    if(nNA[i]>maxNAs){
      next
    }
    if(nGroups>0){
      id2G <- apply(redGeno,2,FUN=function(x){
        identical(x,genotype[,i])|identical(x,reversedGenotype[,i])
      })
      if(any(id2G)){
        genotype2group[i] <- which(id2G)
      }else{
        redGeno <- cbind(redGeno,genotype[,i])
        nGroups <- nGroups+1
        genotype2group[i] <- nGroups
      }
    }else{
      redGeno <- cbind(redGeno,genotype[,i])
      genotype2group <- 1
      nGroups <- 1
    }
  }
  group2genotype <- lapply(1:nGroups,FUN=function(x){which(genotype2group==x)})
  
  ###estimate the population-structure using the reduced genotype
  K<-emma.kinship(t(redGeno),use="pairwise.complete.obs")
  KPCA <- eigen(K)
  useEigen <- 1:which(cumsum(KPCA[[1]])>=sum(KPCA[[1]])*propVar)[1]
  populationStructure <- KPCA[[2]][,useEigen,drop=F]
  colnames(populationStructure) <- paste0("popStr",1:ncol(populationStructure))
  
  ###parse the genotype and the population-structure covariates according to sampleInfo
  mappingGenotype <- redGeno[sampleInfo,,drop=F]
  populationStructure <- populationStructure[sampleInfo,,drop=F]
  
  ###prepare a list with missing genotypes for imputation
  NApos <- extractNAs(mappingGenotype)
  
  ###scale the phenotype if needed
  if(scale){
    if(is.vector(phenotype)){
      phenotype <- scale(phenotype,center=T,scale=T)[,1]
    }else{
      phenotype <- t(apply(phenotype,1,FUN=function(x){
        out <- scale(x,center=T,scale=T)[,1]
        return(out)
      }))
    }
  }
  
  ###return output
  return(list(genotype=mappingGenotype,
              genotype2group=genotype2group,
              group2genotype=group2genotype,
              mappingCovariates=populationStructure,
              NAlist=NApos,
              phenotype=phenotype))
}

#' function that calculates the allele frequencies (1 or 0) for each marker (column) and identifies missing values.
#'
#' @param genotype genotype-matrix with one row per individual, one column per marker and binary marker information (1,0,NA)
#'
#' @return list with one vector for each marker. The first element is the frequency of allele 1. 
#' If there are missing values the indices of the respective samples follow.
#' @export
#'
#' @examples
extractNAs <- function(genotype){

  apply(genotype,2,FUN=function(predictor){
    out <- length(which(predictor==1))/length(which(!is.na(predictor)))
    if(any(is.na(predictor))){
      out <- c(out,which(is.na(predictor)))
    }
    return(out)
  })
}

#' function that replaces missing genotypes according to the local allele frequencies
#'
#' @param genotype genotype as a matrix with 1,0,NA. Predictors that are indicated in exclude may be quantitative.
#' @param NAlist list object with one vector per predictor. The vector has to contain the frequency of the 1-allele first. Indices of individuals with missing values follow.
#'
#' @return Matrix with the same dimensions as the input-genotype. NAs are replaced with alleles (0,1) according to the allele-frequencies in the input (NAlist).
#' @export
#'
#' @examples
replaceGenoNAs <- function(genotype,NAlist){

  sapply(1:ncol(genotype),FUN=function(x){ #per predictor
    lNA <- length(NAlist[[x]])-1
    if(lNA==0){
      return(genotype[,x]) #if no values are missing, return the original vector
    }
    genotype[NAlist[[x]][2:(lNA+1)],x] <- sample(c(1,0),size=lNA,prob = c(NAlist[[x]][1],1-NAlist[[x]][1]),replace=T) #sample NAs according to NAlist
    return(genotype[,x]) #return the extended genotype
  })
}

#' extracts a combined importance score from a randomForest object.
#'
#' @param rf A randomForest-object that contains both permutation importance (PI) and the increase in node purity (RSS). Random Forests that were generated with importance=TRUE are suited for this.
#' @param predictorsOfInterest Predictors for which the combined score should be computed. As a default all predictors are selected.
#' @param normalize Should the scores be normalized with the trait variance? This makes the importance scores comparable across different traits.
#'
#' @return The output consists of a numeric vector containing a score for each predictor of interest. The score is the product of RSS and PI. Negative RSS- or PI-scores are set to zero first. These scores are influenced by correlation between predictors and therefore contain a marker-specific bias.
#' @export
#'
#' @examples
cScore <- function(rf=NULL,predictorsOfInterest=NULL,normalize=TRUE){
  if(class(rf)!="randomForest"){
    stop("A randomForest object is required as input!")
  }
  if(is.null(predictorsOfInterest)){
    predictorsOfInterest <- 1:nrow(rf$importance)
  }else{
    if((!is.numeric(predictorsOfInterest))|(!is.vector(predictorsOfInterest))){
      stop("Predictors of interest have to be entered as a numeric vector!")
    }
  }
  if(ncol(rf$importance)<2){
    stop("A randomForest object with both permutation importance and increase in node purity (importance=TRUE) is required as input!")
  }
  if(normalize){ #divide scores by the trait variance. default
    if(is.null(rf$y)|rf$type!="regression"){
      stop("Normalization can only be performed if the randomForest object was created by supervised regression!")
    }
    tVar <- var(rf$y)
    PI <- rf$importance[predictorsOfInterest,1]/tVar
    PI[which(PI<0)] <- 0
    RSS <- rf$importance[predictorsOfInterest,2]/tVar
    RSS[which(RSS<0)] <- 0
  }else{ #don't normalize
    PI <- rf$importance[predictorsOfInterest,1]
    PI[which(PI<0)] <- 0
    RSS <- rf$importance[predictorsOfInterest,2]
    RSS[which(RSS<0)] <- 0
  }
  combinedScore <- PI*RSS
  names(combinedScore) <- rownames(rf$importance)[predictorsOfInterest]
  return(combinedScore)
}


#' function that generates a permutation-scheme
#'
#' @param permutationGroups list with an entry for each sample indicating within which group it may be permuted. Entries may contain more than one integer. If this is the case one is chosen randomly.
#'
#' @return permuted order of individuals
#' @export
#'
#' @examples
pscheme <- function(permutationGroups){

  permutationGroups <- sapply(permutationGroups,FUN=function(pG){
    if(length(pG)==1){
      return(pG)
    }else{
      return(sample(pG,size=1))
    }
  })
  out <- rep(NA,length(permutationGroups))
  for(i in unique(permutationGroups)){
    pos <- which(permutationGroups==i)
    out[pos] <- sample(pos,size=length(pos))
  }
  return(out)
}

#' wrapper function for qtl mapping with random forest
#'
#' @param mappingData A list containing all necessary information to map a trait as returned by preMap.
#' @param nTrait If the phenotype in mapping data is a matrix, nTrait specifies which of the multiple traits provided, should be mapped. Should be a single integer.
#' @param permutationGroups optional; list with an entry for each sample indicating within which group it may be permuted. Entries may contain more than one integer. If this is the case one is chosen randomly. Samples with the same integer may be permuted. If this list is not provided, all samples are assumed to be in the same permutation-group.
#' @param permute logical, should permutations be performed (if false the real trait-values are used)?
#' @param nPermutations number of permutations to be performed. If permute is false, nPermutations is irrelevant.
#' @param nforest number of randomizations for missing genotypes
#' @param ntree number of trees to grow per subforest
#' @param file optional; if the results are to be saved as a file, this argument specifies the exact paths
#' @param nCl integer; number of cores to use
#' @param clType character; type of cluster to create. Has to be a cluster-type implemented in snow.
#' @param pMat optional; matrix containing predefined permutation schemes. The number of rows has to equal the number of samples and the number of columns has to equal the number of permutations.
#'
#' @return for unpermuted traits: a numerical vector with importance scores for all predictors aside from those in exclude. for permuted traits: a matrix of scores with one column per permutation
#' @export
#'
#' @examples
rfMapper <- function(mappingData, nTrait=NULL, permutationGroups=NULL, permute, nPermutations=1000, nforest, ntree, file=NULL, nCl=1, clType="SOCK", pMat=NULL){
  if(!is.logical(permute)){
    stop("permute has to be TRUE or FALSE.")
  }
  if(any(!c("genotype","genotype2group","group2genotype","mappingCovariates","NAlist","phenotype")%in%names(mappingData))){
    stop("mappingData has to generated with preMap.")
  }
  if(ntree<1|nforest<1){
    stop("Choose parameters so at least one tree is built.")
  }
  if(!is.null(file)){
    dirString <- strsplit(x = file,split = "")[[1]]
    if(grepl("/",file)){
      dir <- paste(dirString[1:tail(which(dirString=="/"),n=1)],collapse = "")
      if(!file.exists(dir)){
        stop("Directory doesn't exist.")
      }
    }
    if(!grepl(".RData",file,ignore.case=T)){
      stop("Specify a RData file.")
    }
  }
  if(!permute|nCl<=1){
    out <- rfMapperNonPar(mappingData=mappingData,
                          nTrait=nTrait,
                          permutationGroups=permutationGroups,
                          permute=permute,
                          nPermutations=nPermutations,
                          nforest=nforest,
                          ntree=ntree,
                          file=file,
                          pMat=pMat)
    return(out)
  }else{
    out <- rfMapperPar(mappingData=mappingData,
                       nTrait=nTrait,
                       permutationGroups=permutationGroups,
                       nPermutations=nPermutations,
                       nforest=nforest,
                       ntree=ntree,
                       file=file,
                       nCl=nCl,
                       clType=clType,
                       pMat=pMat)
    return(out)
  }
}

#' Parallelized implementation of trait mapping with Random Forest
#'
#' @param mappingData A list containing all necessary information to map a trait as returned by preMap.
#' @param nTrait If the phenotype in mapping data is a matrix, nTrait specifies which of the multiple traits provided, should be mapped. Should be a single integer.
#' @param permutationGroups optional; list with an entry for each sample indicating within which group it may 
#' @param file optional; if the results are to be saved as a file, this argument specifies the exact path
#' @param permute logical, should permutations be performed (if false the real trait-values are used)?
#' @param nPermutations number of permutations to be performed. If permute is false, nPermutations is irrelevant.
#' @param nforest number of randomizations for missing genotypes
#' @param ntree number of trees to grow per subforest
#' @param pMat optional; matrix containing predefined permutation schemes. The number of rows has to equal the number of samples and the number of columns has to equal the number of permutations.
#' @param nCl integer; number of cores to use
#' @param clType character; type of cluster to create. Has to be a cluster-type implemented in snow.
#'
#' @return a matrix of scores with one column per permutation
#' @export
#'
#' @examples
rfMapperPar <- function(mappingData, nTrait=NULL, permutationGroups=NULL, file=NULL, permute=T, nPermutations=1000, nforest, ntree, pMat=NULL, nCl, clType="SOCK"){

  genotype <- mappingData$genotype
  phenotype <- mappingData$phenotype
  
  if(!is.vector(phenotype)){
    if(!is.null(nTrait)&is.numeric(nTrait)&length(nTrait)==1){
      if(nTrait%%1!=0|nTrait<1|nTrait>nrow(phenotype)){
        stop("Specify nTrait as a natural number between 1 and the number traits supplied.")
      }
      phenotype <- phenotype[nTrait,]
    }else{
      stop("Specify which trait to map.")
    }
  }
  
  NAlist <- mappingData$NAlist
  mappingCovariates <- mappingData$mappingCovariates
  exclude <- (ncol(genotype)+1):(ncol(genotype)+ncol(mappingCovariates))
  genotype <- cbind(genotype,mappingCovariates)
  if(is.null(permutationGroups)){
    permutationGroups <- as.list(rep(1,nrow(genotype)))
  }
  
  if(any(is.na(phenotype))&is.null(pMat)){ #remove individuals with missing trait-values if pMat
    if(all(is.na(phenotype))){
      stop("No data available for the specified trait.")
    }
    missing <- which(is.na(phenotype))
    phenotype <- phenotype[-missing]
    genotype <- genotype[-missing,]
    NAlist <- lapply(NAlist,FUN=function(markerVec){
      aFreq <- markerVec[1]
      markerVec <- markerVec[-1]
      if(length(markerVec)==0){
        return(aFreq)
      }
      markerVec <- sapply(markerVec,FUN=function(sampleInd){
        sampleInd-sum(sampleInd>missing)
      })
      return(c(aFreq,markerVec))
    })
    permutationGroups <- permutationGroups[-missing]
  }
  
  #make the cluster
  require(snow)
  cl <- makeCluster(nCl,type=clType)
  clusterExport(cl,
                list = c("cScore",
                         "pscheme",
                         "genotype",
                         "phenotype",
                         "permutationGroups",
                         "file",
                         "nforest",
                         "ntree",
                         "pMat",
                         "NAlist",
                         "replaceGenoNAs",
                         "exclude"),
                environment())  
  clusterEvalQ(cl,library(randomForest))
  
  
  
  out <- parSapply(cl=cl,X=1:nPermutations,FUN=function(nPerm){
    if(!is.null(pMat)){
      pScheme <- pMat[,nPerm] #use a permutation scheme from pMat
    }else{
      pScheme <- pscheme(permutationGroups) #generate a permutation scheme
    }
    permPhenotype <- phenotype[pScheme] #permute the trait vector
    permGenotype <- genotype 
    permGenotype[,exclude] <- permGenotype[pScheme,exclude] #preserve the association between population structure and trait values
    if(any(is.na(phenotype))&!is.null(pMat)){ #remove individuals with missing trait-values
      missing <- which(is.na(phenotype))
      phenotype <- phenotype[-missing]
      permGenotype <- permGenotype[-missing,]
      NAlist <- lapply(NAlist,FUN=function(markerVec){
        aFreq <- markerVec[1]
        markerVec <- markerVec[-1]
        if(length(markerVec)==0){
          return(aFreq)
        }
        markerVec <- sapply(markerVec,FUN=function(sampleInd){
          sampleInd-sum(sampleInd>missing)
        })
        return(c(aFreq,markerVec))
      })
    }
    rfs <- lapply(1:nforest,FUN=function(nF){
      permGenotype[,-exclude] <- replaceGenoNAs(genotype=permGenotype[,-exclude],NAlist=NAlist)
      rf <- randomForest(x=permGenotype, y=permPhenotype, ntree=ntree, importance=TRUE, keep.forest=F)
      return(rf)
    })
    
    bigRf <- do.call("combine",rfs) #combine the subforests that were generated with genotypes that differ for the NA positions
    if(!is.null(exclude)){
      scores <- cScore(bigRf)[-exclude] #don't keep the importance scores of the population-structure, since it is not permuted
    }else{
      scores <- cScore(bigRf)
    }
    return(scores)
  })
  stopCluster(cl)
  if(!is.null(file)){
    save(out,file=file)
  }
  return(out) #return the scores as a vector for unpermuted traits and as a matrix for permuted traits
}

#' Non-parallelized implementation of trait mapping with Random Forest
#'
#' @param mappingData A list containing all necessary information to map a trait as returned by preMap.
#' @param nTrait If the phenotype in mapping data is a matrix, nTrait specifies which of the multiple traits provided, should be mapped. Should be a single integer.
#' @param permutationGroups optional; list with an entry for each sample indicating within which group it may be permuted. Entries may contain more than one integer. If this is the case one is chosen randomly.
#' @param permute logical, should permutations be performed (if false the real trait-values are used)?
#' @param nPermutations number of permutations to be performed. If permute is false, nPermutations is irrelevant.
#' @param nforest number of randomizations for missing genotypes
#' @param ntree number of trees to grow per subforest
#' @param file optional; if the results are to be saved as a file, this argument specifies the exact path
#' @param pMat optional; matrix containing predefined permutation schemes. The number of rows has to equal the number of samples and the number of columns has to equal the number of permutations.
#'
#' @return a numerical vector with importance scores for all predictors, not including covariates
#' @export
#'
#' @examples
rfMapperNonPar <- function(mappingData, nTrait=NULL, permutationGroups=NULL, permute, nPermutations=1000, nforest, ntree, file=NULL, pMat=NULL){

  genotype <- mappingData$genotype
  phenotype <- mappingData$phenotype
  if(!is.vector(phenotype)){
    if(!is.null(nTrait)&is.numeric(nTrait)&length(nTrait)==1){
      if(nTrait%%1!=0|nTrait<1|nTrait>nrow(phenotype)){
        stop("Specify nTrait as a natural number between 1 and the number traits supplied.")
      }
      phenotype <- phenotype[nTrait,]
    }else{
      stop("Specify which trait to map.")
    }
  }
  NAlist <- mappingData$NAlist
  mappingCovariates <- mappingData$mappingCovariates
  exclude <- (ncol(genotype)+1):(ncol(genotype)+ncol(mappingCovariates))
  genotype <- cbind(genotype,mappingCovariates)
  if(is.null(permutationGroups)){
    permutationGroups <- as.list(rep(1,nrow(genotype)))
  }
  
  if(any(is.na(phenotype))&is.null(pMat)){ #remove individuals with missing trait-values if pMat is not provided
    if(all(is.na(phenotype))){
      stop("No data available for the specified trait.")
    }
    missing <- which(is.na(phenotype))
    phenotype <- phenotype[-missing]
    genotype <- genotype[-missing,]
    NAlist <- lapply(NAlist,FUN=function(markerVec){
      aFreq <- markerVec[1]
      markerVec <- markerVec[-1]
      if(length(markerVec)==0){
        return(aFreq)
      }
      markerVec <- sapply(markerVec,FUN=function(sampleInd){
        sampleInd-sum(sampleInd>missing)
      })
      return(c(aFreq,markerVec))
    })
    permutationGroups <- permutationGroups[-missing]
  }
  if(permute){ #do permutations
    out <- sapply(1:nPermutations,FUN=function(nPerm){
      if(!is.null(pMat)){
        pScheme <- pMat[,nPerm] #use a permutation scheme from pMat
      }else{
        pScheme <- pscheme(permutationGroups) #generate a permutation scheme
      }
      permPhenotype <- phenotype[pScheme] #permute the trait vector
      permGenotype <- genotype
      permGenotype[,exclude] <- permGenotype[pScheme,exclude] #preserve the association between population structure and trait values
      if(any(is.na(phenotype))&!is.null(pMat)){ #remove individuals with missing trait-values
        missing <- which(is.na(phenotype))
        phenotype <- phenotype[-missing]
        permGenotype <- permGenotype[-missing,]
        NAlist <- lapply(NAlist,FUN=function(markerVec){
          aFreq <- markerVec[1]
          markerVec <- markerVec[-1]
          if(length(markerVec)==0){
            return(aFreq)
          }
          markerVec <- sapply(markerVec,FUN=function(sampleInd){
            sampleInd-sum(sampleInd>missing)
          })
          return(c(aFreq,markerVec))
        })
      }
      rfs <- lapply(1:nforest,FUN=function(nF){
        permGenotype[,-exclude] <- replaceGenoNAs(genotype=permGenotype[,-exclude],NAlist=NAlist)
        rf <- randomForest(x=permGenotype, y=permPhenotype, ntree=ntree, importance=TRUE)
        return(rf)
      })
      bigRf <- do.call("combine",rfs) #combine the subforests that were generated with genotypes that differ for the NA positions
      scores <- cScore(bigRf)[-exclude] #don't keep the importance scores of the population-structure, since it is not permuted
    })
  }else{ #map real traits
    rfs <- lapply(1:nforest,FUN=function(nF){
      genotype[,-exclude] <- replaceGenoNAs(genotype=genotype[,-exclude],NAlist=NAlist)
      rf <- randomForest(x=genotype, y=phenotype, ntree=ntree, importance=TRUE)
      return(rf)
    })
    bigRf <- do.call("combine",rfs)
    out <- cScore(bigRf)[-exclude]
  }
  
  if(!is.null(file)){
    save(out,file=file)
  }
  return(out) #return the scores as a vector for unpermuted traits and as a matrix for permuted traits
}


#' function that joins consecutive significantly linked markers to QTL
#'
#' @param QTLmat matrix consisting of two columns (target, predictor) X one row per linked predictor
#' @param chrVec character-vector containing the chromosome on which a given marker is located
#'
#' @return list with one entry for each regulated trait. Each entry contains the target and a matrix with the start and end of the QTL-regions
#' @export
#'
#' @examples
joinConsecutive <- function(QTLmat, chrVec){

  #identify traits with at least one qtl
  targets <- sort(unique(QTLmat[,1]))
  
  #join consecutive significant predictors per trait
  QTLperTarget <- lapply(targets,FUN=function(target){
    targetMat <- QTLmat[which(QTLmat[,"target"]==target),] #subset for the trait
    if(is.vector(targetMat)){ #if there is only a single significant predictor, there is no joining to do
      outMat <-matrix(c(targetMat[2],targetMat[2]),nrow=1)
      colnames(outMat) <- c("start", "end")
      return(list(target=target,predictors=outMat))
    }
    predictors <- sort(targetMat[,2])
    col <- predictors[1]
    outMat <- NULL
    for(i in 2:length(predictors)){
      #join loci if they are on the same chromosome and are direct neighbors
      if(chrVec[col[1]]==chrVec[predictors[i]]&(predictors[i]-col[length(col)])==1){
        col <- c(col,predictors[i])
      }else{
        outMat <- rbind(outMat,c(min(col),max(col)))
        col <- predictors[i] #start with a new region
      }
    }
    outMat <- rbind(outMat,c(min(col),max(col)))
    colnames(outMat) <- c("start","end")
    #return one matrix per trait containing start and end points of all qtl
    return(list(target=target,predictors=outMat)) 
  })
  #return a list containg these matrices for all regulated traits
  return(QTLperTarget)
}

#' function that joins different QTL to a single one if they contain highly correlated markers 
#'
#' @param QTLlist output from joinConsecutive
#' @param corMat squared matrix of predictor-correlation coefficients
#' @param corThreshold threshold indicating how strong QTL-regions have to be correlated to be joined
#' @param distThreshold threshold indicating how close QTL-regions have to be to be joined
#' @param chrVec character vector containing the information on which chromosome a marker is located
#'
#' @return list with an entry for each regulated trait. Each entry contains the target and a matrix with the start and end of the QTL-regions after combining close and correlated QTL
#' @export
#'
#' @examples
joinNear <- function(QTLlist, corMat, corThreshold, distThreshold, chrVec){

  out <- lapply(QTLlist,FUN=function(targetList){
    qtl <- targetList[[2]] #matrix with QTL-regions for this trait
    if(nrow(qtl)==1){ #if there is only one region, no regions can be joined
      return(list(target=targetList[[1]],predictors=qtl))
    }
    prog <- T #progress-variable
    while(prog){ #as long as there is progress
      nextProg <- F #loop stops if there is no change
      #each qtl is tested to its right neighbor
      for(i in 1:(nrow(qtl)-1)){
        #are the QTL correlated and close enough?
        if(max(corMat[qtl[i,1]:qtl[i,2],qtl[i+1,1]:qtl[i+1,2]])>=corThreshold&(qtl[i+1,1]-qtl[i,2])<=distThreshold){
          qtl[i,] <- c(qtl[i,1],qtl[i+1,2]) #end of the combined QTL is the end of the second QTL
          qtl <- qtl[-(i+1),] #the second QTL is deleted as it is part of the combined QTL now
          qtl <- matrix(as.vector(qtl),ncol=2,byrow=F) #if there is one QTL remaining it still has the same format
          nextProg <- T
          break #break for loop because indices and matrix size are changed
        }
      }
      prog <- nextProg #prepare for the next iteration
      if(nrow(qtl)==1){ #if only one QTL remains, no further joining can be done
        prog <- F 
      }
    }
    colnames(qtl) <- c("start","end")
    #return the target and the QTL-matrix with start and end points
    return(list(target=targetList[[1]],predictors=qtl)) 
  })
  return(out)
}


#' join QTL-regions that are distant from each other but contain highly correlated markers to QTL-groups using hierarchical clustering
#'
#' @param QTLlist output of joinNear
#' @param corMat squared matrix of predictor-correlation coefficients
#' @param corThreshold threshold indicating how strong QTL-regions have to be correlated to be joined
#'
#' @return list with an entry for each QTL. Each entry contains the target and a matrix with the start and end of the QTL-regions after grouping highly correlated distant QTL.
#' @export
#'
#' @examples
joinCorrelated <- function(QTLlist, corMat, corThreshold){

  out <- lapply(QTLlist,FUN=function(targetList){
    qtl <- targetList[[2]] #matrix with QTL-regions for this trait
    if(nrow(qtl)==1){ #if there is only one region, no regions can be joined
      return(list(list(target=targetList[[1]],predictors=qtl)))
    }
    mat <- sapply(1:nrow(qtl), function(i) { #correlation matrix
      sapply(1:nrow(qtl), function(j) {
        return(max(corMat[qtl[i,1]:qtl[i,2],qtl[j,1]:qtl[j,1]], na.rm=T))
      })
    })
    colnames(mat) <- rownames(mat) <- 1:nrow(qtl)
    diag(mat) <- 0 # remove the diagonal
    mat[is.na(mat)] <- 0
    #here we convert the correlation-matrix into distances for clustering 
    dcorm <- as.dist((1 - mat)/2) 
    dcor_cut <- (1-corThreshold)/2 #convert the threshold accordingly
    dcorm_hclust <- hclust(dcorm)
    group <- cutree(hclust(dcorm), h=dcor_cut) #apply the threshold to the tree
    group <- split(as.numeric(names(group)), group)
    out <- lapply(group,FUN=function(gr){ #join regions to QTL-groups
      outMat <- matrix(qtl[gr,],nrow=length(gr),byrow = F)
      colnames(outMat) <- c("start","end")
      list(target=targetList[[1]],predictors=outMat)
    })
    return(out)
  })
  out <- do.call(c,out) #parse list so that one QTL-group is one entry regardless of the target
  names(out) <- NULL
  return(out)
}

#' function that uses finished permutations to estimate an empirical p-value
#'
#' @param path folder where the .RData files with the matrices of empirical values can be found. They have to conform to predictors (rows) X permutations (columns).
#' @param markersPerIteration total number of null distributions that are loaded at the same time. The possible number is dependent on available working memory and size of the null distributions.
#' @param scores real predictor-scores to be compared to the empirical null distributions. If scores is a matrix, it has to conform to traits (rows) X predictors (columns). Scores for the same predictor in matrix, it has to conform to traits (rows) X predictors (columns). Scores for the same predictor in different traits (rows) are compared to the same null distribution.
#' @param printProg logical, if TRUE the last finished predictor is printed
#' @param pCorrection Which type, if any, of multiple testing correction should be performed? Has to be either "none", "fdr" or "bonferroni".
#'
#' @return A matrix with empirical p-values is returned.
#' @export
#'
#' @examples
pEst <- function(path,markersPerIteration,scores,printProg=T,pCorrection="none"){

  #if scores is a vector (i.e. only one trait), it is transformed into a matrix
  if(is.vector(scores)){
    scores <- matrix(scores,nrow=1)
  }
  
  #identify data objects
  files <- list.files(path)
  files <- files[grep(files,pattern = ".RData")] 
  
  #separate the predictors into packages according to markersPerIteration
  blocks <- seq(1,ncol(scores),markersPerIteration)
  
  #get null distributions per block and compare it to the actual scores
  pMatBlock <- lapply(1:length(blocks),FUN=function(blockN){
    if(blockN<length(blocks)){
      preds <- blocks[blockN]:(blocks[blockN+1]-1)
    }else{
      preds <- blocks[blockN]:ncol(scores)
    }
    #load the null distributions and bind them into a single matrix
    nullMat <- lapply(files,FUN=function(file){
      if("out"%in%ls()){
        rm(out)
      }
      load(paste0(path,file))
      if(!"out"%in%ls()){
        return(NULL)
      }
      out <- out[preds,]
      return(out)
    })
    nullMat <- do.call(cbind, nullMat)
    
    #if no files were available or did not contain distributions stop here
    if(is.null(dim(nullMat))){
      stop("No null distributions can be found in the chosen directory.")
    }
    
    #compute the size of the null distribution to calculate the smallest possible p-value
    nullSize <- ncol(nullMat)
    
    #compare scores to the empirical null distributions
    pMat <- lapply(1:length(preds),FUN=function(predN){
      null <- nullMat[predN,]
      lnull <- length(null)
      pPred <- sapply(scores[,preds[predN]],FUN=function(scr){
        max(sum(null>=scr)/lnull,1/lnull)
      })
      return(pPred)
    })
    pMat <- do.call(cbind,pMat)
    print(preds[length(preds)])
    return(pMat)
  })
  #bind all matrices to a single one
  pMatBlock <- do.call(cbind,pMatBlock)
  
  if(pCorrection!="none"){
    pMatBlock[1:length(pMatBlock)] <- p.adjust(pMatBlock[1:length(pMatBlock)],method = pCorrection)
  }
  
  #return results
  if(nrow(pMatBlock)==1){
    return(pMatBlock[1,])
  }else{
    return(pMatBlock)
  }
}


#' function that identifies significant marker-trait associations and joins close markers to a single QTL if they are consecutive and/or highly correlated and joins distant QTL to QTL-groups if they contain highly correlated markers
#'
#' @param pmat matrix with p-values, traits (rows) X loci (columns)
#' @param sigThreshold significance threshold to applied on the p-values
#' @param corThreshold threshold for absolute correlation values above which signifcant loci are treated as linked to each other. If they are close they are joined to a single region, otherwise they are grouped together without including inbetween loci
#' @param distThreshold distance threshold up until which loci containing highly correlated markers including all markers in between are linked
#' @param genotype allele-information for all samples (rows) and markers (columns)
#' @param chrVec character-vector containing the chromosome on which a given marker is located
#'
#' @return list with QTLgroups. Each entry is a list of the target and a matrix with the QTL-regions involved.
#' @export
#'
#' @examples
QTLgrouper <- function(pmat, sigThreshold, corThreshold, distThreshold, genotype, chrVec){
  
  if(is.vector(pmat)){
    pmat <- matrix(pmat,nrow=1)
  }
  if(ncol(pmat)<ncol(genotype)){
    stop("Number of markers in pmat and genotype have to be identical.")
  }
  QTLmat <- which(pmat<=sigThreshold,arr.ind=T)
  
  #check if any significant QTL are present
  if(nrow(QTLmat)==0){
    return("no significant loci")
  }
  colnames(QTLmat) <- c("target","predictor")
  
  #square threshold and correlation matrix to account for extreme negative correlation
  corThreshold <- corThreshold^2
  corMat <- cor(genotype,use="pair")^2
  
  #join consecutive markers regulating the same trait to a single QTL
  QTLlist <- joinConsecutive(QTLmat, chrVec)
  
  #join markers to QTL if they are close and correlated
  QTLlist <- joinNear(QTLlist,corMat,corThreshold,distThreshold,chrVec)
  
  #link distant QTL that are correlated
  QTLlist <- joinCorrelated(QTLlist,corMat,corThreshold)
  
  #find the predictor with the smallest p-value in each qtl
  QTLlist <- minPV(QTLlist,pmat)
  return(QTLlist)
}

#' identify the predictor with the smallest p-value in each QTL
#'
#' @param QTLlist list-object with an entry for each QTL as returned by joinCorrelated
#' @param pmat matrix containing the p-values for each predictor(columns)/trait(row)
#'
#' @return object similar to QTLlist with added information regarding the most significant predictor in each QTL
#' @export
#'
#' @examples
minPV <- function(QTLlist,pmat){
  out <- lapply(QTLlist,FUN=function(qtl){
    target <- qtl$target
    predictors <- apply(qtl$predictors,1,FUN=function(preds){return(preds[1]:preds[2])})
    if(is.list(predictors)){
      predictors <- do.call("c",predictors)
    }
    sigPred <- predictors[which.min(pmat[target,predictors])]
    qtl <- c(qtl,sigPred)
    names(qtl)[length(qtl)] <- "mostSignificantPredictor"
    qtl <- c(qtl,pmat[target,sigPred])
    names(qtl)[length(qtl)] <- "minP"
    return(qtl)
  })
  return(out)
}

#' writes a QTL-list to a .qtl file as a table
#'
#' @param QTLlist QTL stored as a list-object as returned by QTLgrouper
#' @param traitNames optional; a vector with names for the traits which correspond to the rows in the p-value matrix used for QTLgrouper.
#' @param path Complete filepath ending in .qtl, specifying where the results should be saved.
#' @param markerPositions optional; Matrix with three columns specifying the chromosome, start and end for each marker.
#' @param digits optional; P-values are formatted to scientific convention and can be rounded according to digits.
#'
#' @return The output is written in a file, specified by path.
#' @export
#'
#' @examples
writeQTL <- function(QTLlist,traitNames=NULL,path,markerPositions=NULL,digits=NULL){
  if(identical(QTLlist,"no significant loci")){
    stop("no significant loci")
  }
  metaInfo <- lapply(QTLlist,names)
  if(length(QTLlist)>1){
    allId <- sapply(1:(length(metaInfo)-1),FUN=function(x){
      identical(metaInfo[x],metaInfo[x+1])
    })
    if(!all(allId)){
      stop("corrupted QTL-list")
    }
  }
  if(!is.null(markerPositions)){
    if(!is.matrix(markerPositions)&!is.data.frame(markerPositions)){
      stop("markerPositions has to be a matrix or data frame.")
    }
    if(ncol(markerPositions)!=3){
      stop("markerPositions has to contain three columns, specifying the chromosome, start and end for each marker.")
    }
  }
  if(!is.character(path)){
    stop("Specify an exact path.")
  }
  if(!grepl(".qtl",path)){
    stop("path has to end in .qtl")
  }
  dirString <- strsplit(x = path,split = "")[[1]]
  if(grepl("/",path)){
    dir <- paste(dirString[1:tail(which(dirString=="/"),n=1)],collapse = "")
    if(!file.exists(dir)){
      stop("Directory doesn't exist.")
    }
  }
  
  col_names <- c("QTL","trait number","trait name","first marker","last marker","chromosome","first position","last position","smallest p-value","most significant predictor")
  out <- lapply(QTLlist,FUN=function(x){
    qtlMat <- matrix(NA,ncol=length(col_names),nrow=nrow(x$predictors))
    qtlMat[,1] <- 0
    qtlMat[1,1] <- 1
    qtlMat[,2] <- x$target
    if(!is.null(traitNames)){
      qtlMat[,3] <- traitNames[x$target]
    }
    qtlMat[,4:5] <- x$predictors
    if(!is.null(markerPositions)){
      qtlMat[,6] <- markerPositions[x$predictors[,1],1]
      qtlMat[,7] <- markerPositions[x$predictors[,1],2]
      qtlMat[,8] <- markerPositions[x$predictors[,2],3]
    }
    qtlMat[,9] <- format(x$minP,scientific=T,digits=digits)
    qtlMat[,10] <- x$mostSignificantPredictor
    return(qtlMat)
  })
  out <- do.call("rbind",out)
  colnames(out) <- col_names
  out[,1] <- cumsum(as.numeric(out[,1]))
  write.table(x = out,file = path,append = F,quote = F,sep = "\t",col.names = T,row.names = F)
}

#' Reconstructs a QTLlist from a .qtl file
#'
#' @param path The location of the .qtl file that should be loaded.
#'
#' @return A list-object as it is returned by QTLgrouper.
#' @export
#'
#' @examples
readQTL <- function(path){
  if(!is.character(path)|!grepl(".qtl",path)){
    stop("Specify a .qtl file.")
  }
  if(!file.exists(path)){
    stop("Specify an existing .qtl file.")
  }
  qtlMat <- read.table(path,header=T,sep="\t",as.is = T)
  out <- lapply(1:max(as.numeric(qtlMat[,1])),FUN=function(i){
    submat <- qtlMat[as.numeric(qtlMat[,1])==i,,drop=F]
    target <- i
    predictors <- submat[,4:5,drop=F]
    predictors <- as.matrix(predictors)
    rownames(predictors) <- c()
    colnames(predictors) <- c("start","end")
    mostSignificantPredictor <- submat[1,10]
    minP <- as.numeric(submat[1,9])
    return(list(target=target,predictors=predictors,mostSignificantPredictor=mostSignificantPredictor,minP=minP))
  })
  return(out)
}

#Gratefully taken from emma 1.1.2 (url:http://mouse.cs.ucla.edu/emma/news.html,license: LGPL)
emma.kinship <- function(snps, method="additive", use="all") {
  n0 <- sum(snps==0,na.rm=TRUE)
  nh <- sum(snps==0.5,na.rm=TRUE)
  n1 <- sum(snps==1,na.rm=TRUE)
  nNA <- sum(is.na(snps))
  
  stopifnot(n0+nh+n1+nNA == length(snps))
  
  if ( method == "dominant" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( method == "recessive" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) & (snps == 0.5)] <- flags[!is.na(snps) & (snps == 0.5)]
  }
  else if ( ( method == "additive" ) && ( nh > 0 ) ) {
    dsnps <- snps
    rsnps <- snps
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]
    snps <- rbind(dsnps,rsnps)
  }
  
  if ( use == "all" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  }
  else if ( use == "complete.obs" ) {
    snps <- snps[rowSums(is.na(snps))==0,]
  }
  
  n <- ncol(snps)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1
  
  for(i in 2:n) {
    for(j in 1:(i-1)) {
      x <- snps[,i]*snps[,j] + (1-snps[,i])*(1-snps[,j])
      K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
      K[j,i] <- K[i,j]
    }
  }
  return(K)
}