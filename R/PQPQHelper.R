#' Remove Columns with More than N NAs
#'
#' This functionr finds the number of NAs in each column, and removes any columns with more than
#' n missing values.
#' @param df The matrix with possible missing values
#' @param n The upper limit of NAs allowed in a column
#' @keywords missing values

removeCols_nNAs <- function(df, n){
  # print('in removeCols nNAs')
  NaSbyCols <- colSums(is.na(df))
  df[,which(NaSbyCols <= n), drop=FALSE]
}

#' Identify Proteoforms from Peptides
#'
#' This function calls removeCols_nNAs(), clusterdata(), and peptideWithTheHighestAbundance().
#' It iessentially a wrapper function to make model peptide identification easier.
#' @param Model_Peptide_Candidates A n-by-m matrix including n peptides to be clustered.
#' @param AbundanceOnly Model_Peptide_Candidates, but without rownames or labels.
#' @param protein The name of the given protein, to be used for naming purposes.
#' @keywords cluster
#' @keywords proteoform
#' @seealso \code{\link{clusterdata}}

proteoform_ID_Four_Plus <- function(Model_Peptide_Candidates,AbundanceOnly, protein){
  # Allows for successful pairwise clustering
  fixed <- removeCols_nNAs(t(Model_Peptide_Candidates), dim(Model_Peptide_Candidates)[2] - 1)
  # cluster data; threshold defaults to 1
  species <- as.data.frame(clusterdata(t(fixed)))

  a <- sapply(species[,3],function(x) {peptideWithHighestAbundance(AbundanceOnly[x,])} )
  b <- paste(protein,species$clusterNumber)

  ModelPeptidePairs <-  data.frame("Protein" = b,
                                   "ModelPeptide" = a)
  return(ModelPeptidePairs)

}

#' Find Peptide with the Highest Abundance
#'
#' This function finds a given peptide/row that has the highest mean value.
#' @param AbundanceOnly An n-by-(m - 1) matrix including only the sample measurements for the n peptides.
#' @keywords Abundance
#' @keywords Model Peptide

peptideWithHighestAbundance <- function(AbundanceOnly){

  ###Highest Average Abundance
  MeansList <- rowMeans(AbundanceOnly,na.rm = TRUE)

  #Model Peptides have the highest average abundance
  ModelPeptide <-  names(MeansList[MeansList == max(MeansList)])
  if(length(ModelPeptide) > 1) ModelPeptide <- ModelPeptide[1]

  return(ModelPeptide)
}

#' Find All Nodes Beneath A Given Hclust Node
#'
#' This recursive function finds all nodes beneath the given node and returns a complete list.
#' @param RowNum Something - for later
#' @param hc The hierarchal clustering result (see the package 'cluster' for details)
#' @keywords cluster
#' @keywords tree
#' @keywords recursive

findAllNodesBeneath <- function(RowNum, hc){
  nodeList <- c(RowNum)
  mergelist <- hc$merge
  NumSingletons<- sum(mergelist[RowNum,]<0)

  if(NumSingletons == 1) {

    posValue <- mergelist[RowNum,mergelist[RowNum,]>0]
    nodeList <- c(nodeList, findAllNodesBeneath(posValue, hc))

  }

  if(NumSingletons == 0){

    nodeList <- c(nodeList, findAllNodesBeneath(mergelist[RowNum,1], hc))
    nodeList <- c(nodeList, findAllNodesBeneath(mergelist[RowNum,2], hc))
  }

  # if both are singletons, there are no other branches to follow
  return(nodeList)
}

#' Find Inconsistency Values Of Cluster Nodes
#'
#' This function computes the inconsistency values for each node in the hierarchal clustering
#' result from cluster::hclust(). The algorithm is the one described as being used in the Matlab
#' function 'inconsistent'.
#' @param hclusted This is the hierarchal clustering output, an object with

findInconsistency <- function(hclusted){
  mergelist <- hclusted$merge

  inconsistent_frame <- data.frame("mean" <- rep(c(0),dim(mergelist)[1]),
                                   "stdev" <- rep(c(0),dim(mergelist)[1]),
                                   "number" <- rep(c(1),dim(mergelist)[1]),
                                   "inconsistent" <- rep(c(0),dim(mergelist)[1]))
  colnames(inconsistent_frame) <- c("mean","stdev","number","inconsistent")

  # replaces cluster indices with height values
  mergelist[which(mergelist > 0)] <- hclusted$height[mergelist[which(mergelist > 0)]]
  # replace negative values (singletons) with zeros
  mergelist[which(mergelist < 0)] <- 0
  # 3 add height of the node itself
  mergelist <- cbind(mergelist,hclusted$height)
  # 4 Find sum (equivalent to matlab s(1)
  mergelist <- cbind(mergelist,rowSums(mergelist))
  # 5 Find sum of squares
  mergelist <- cbind(mergelist,rowSums(mergelist[,1:3]^2))
  # 6 Square column (for stdev calcs)
  mergelist <- cbind(mergelist, mergelist[,4]^2)
  # 7 Find all non-zero values by row (n)
  mergelist <- cbind(mergelist, apply(mergelist[,1:3],1,function(x) sum(x != 0)))

  #Insert mean values in final data frame
  inconsistent_frame$mean <- mergelist[,4]/mergelist[,7]

  #Insert stdev values in final data frame
  inconsistent_frame$stdev <- sqrt((mergelist[,5] - (mergelist[,6]/mergelist[,7]))/(mergelist[,7] -1))

  #Insert number value in final data frame
  inconsistent_frame$number <- mergelist[,7]

  #Insert inconsistent values
  inconsistent_frame$inconsistent <- (mergelist[,3] - inconsistent_frame$mean)/inconsistent_frame$stdev

  inconsistent_frame<- replace(inconsistent_frame, is.na(inconsistent_frame),0)

  return(inconsistent_frame)
}

#' Cluster Dataset Based on Inconsistency Values
#'
#' This function does ...stuff. For later: Actually fill in this bit.
#' @param dataset The peptides that need to be clustered
#' @param threshold The limit for inconsistnecy. Defaults to 1.

clusterdata <- function(dataset, threshold = 1) {
  # print('In Clusterdata')
  distMat <- as.dist(1 - cor(t(dataset),
                             use = 'pairwise.complete.obs',
                             method = 'pearson'))
  print(sum(which(is.na(distMat), arr.ind = TRUE)))
  valmisdat <- 1.1 * max(distMat, na.rm = TRUE)

  distMat[which(is.na(distMat))] <- valmisdat
  # print(distMat)
  hc <-
    hclust(distMat, method = "single")

  number_nodes <-  dim(hc$merge)[1]
  fI <- findInconsistency(hc)
  df <-
    data.frame(
      "ValuesBeneath" = c(1:number_nodes),
      "Inconsistencies" = c(1:number_nodes)
    )
  # print(hc$merge)
  for (node in number_nodes:1) {
    children <- findAllNodesBeneath(node, hc)

    mixedVal <-  unlist(lapply(children, function(x) {
      hc$merge[x, ]
    }))


    df$`ValuesBeneath`[[node]] <- list(abs(mixedVal[mixedVal < 0]))

    df$Inconsistencies[[node]] <- list(fI$inconsistent[children])

  }

  clusterIndex <-
    which(sapply(df$Inconsistencies, function(x) {
      sum(unlist(x) > threshold)
    }) == 0, arr.ind = T)

  clusterCand <- df[clusterIndex, 1]

  finalClusters <- data.frame("members" = 1:(number_nodes + 1)) # , 'height' = rep(0, number_nodes + 1))

  for (i in 1:(number_nodes + 1)) {
    # Searches the list of lists for anytime a value is equal to i and returns those lists
    withI <-
      clusterCand[unlist(lapply(clusterCand, function(x)
        sapply(x, function(x)
          any(x == i))))]
    # Choose the sublist with the most elements inside
    if (length(withI) > 0) {
      numelems <-
        lapply(withI, function(x)
          sapply(x, function(y)
            length(y)))
      biggestCluster <-
        withI[which(numelems == max(unlist(numelems)), arr.ind = T)][[1]]
      finalClusters$members[[i]] <-  biggestCluster
    }
  }

  finalClusters <-  unique(finalClusters)
  finalClusters$clusterNumber <-  1:(dim(finalClusters)[1])
  names <- as.matrix(rownames(dataset))
  stringCand <- sapply(clusterCand, lapply, unlist)

  numElem <- rep(0, dim(finalClusters)[1])
  height <- rep(0, dim(finalClusters)[1])

  for (i in 1:(dim(finalClusters)[1])) {
    finalClusters[i,3][[1]] <- list(names[unlist(finalClusters[i, 1])])
    numElem[i] <- length(unlist(finalClusters$members[i]))

    # If there's more than 1 cluster get heights
    if(dim(finalClusters)[1] > 1) {

      # If more than 1 element
      if (numElem[i] > 1) {

        stringPeptides <- paste(finalClusters$members[i][[1]][[1]])

        candMatch <- lapply(stringCand, function(x) {finalClusters$members[i][[1]][[1]] %in% x})
        candMatch <- which.max(lapply(candMatch, sum))

        height[i] <- hc$height[candMatch]
      }
      # If only 1 element
      else {
        # Find entry with -num in hc
        num = -1 * unlist(finalClusters$member[i])
        height[i] <- hc$height[which(hc$merge[,1] == num)]
      }
    } # End if more than 1 cluster

    # MATLAB testing
    #finalClustedMat[unlist(finalClusters[i, 1]), 2] <- i
  }
  finalClusters <- finalClusters[order(height, decreasing = FALSE),]
  numElem <- numElem[order(height)]
  finalClusters <- finalClusters[order(numElem, decreasing = TRUE),]
  return(finalClusters)
}

#' Find Peptides With Significant Correlation to Model Peptide
#'
#' This function uses the correlation and p-values matrices input to determine which peptides are higly
#' correlated with the model peptide.
#' @param ModelPeptide A number representing which row/column in the correlation matrix represents the
#' model peptide
#' @param CorMat The correlation matrix for the peptides, created by psych::cor().
#' @param p_values The p-values matrix, created by psych::corr.test()
#' @param threshold The minimum value allowed for p-values. Default is 0.4 based on recommendations by
#' the original paper.
#' #@source <https://www.ncbi.nlm.nih.gov/pubmed/21734112>

findSigPeptides <- function(ModelPeptide, CorMat, p_values, threshold = 0.4) {

  # Avoid being burned by FACTORS
  ModelPeptide <- as.character(ModelPeptide)

  # Find and isolate relevant correlation information
  SigPepCand <- CorMat[ModelPeptide,]
  CandPVals <- p_values[ModelPeptide, ]

  # Create named logical representing all "significant" peptides
  logipep <- (CandPVals < threshold) & (!is.na(CandPVals) & (SigPepCand > 0))

  SigPepCand <- SigPepCand[logipep]
  SigPeps <- names(SigPepCand[!(is.na(SigPepCand) | is.nan(SigPepCand))])

  # Return names of peptides marked TRUE
  return(SigPeps)
}

#' Remove Low Confidence Peptides
#'
#' This function removes columns from a matrix based on the confidence levels of each row.
#' If there are more than 300 rows, only the 300 most-confident are allowed.
#' @param data An n-by-m matrix, which will be returned with rows removed
#' @param confPeps A vector of length n, with confidence values from 0 to 100.
#' @param limit The lowest confidence level allowed. Defaults to 95.
#' @return A matrix with m columns (ie data but with rows removed)

removeLowConf <- function(data, confPeps, limit = 95) {

  # Check if more than 300 peptides
  numPeptides <- dim(data)[1]
  if(numPeptides > 300) {

    orderPeps <- data[order(-confPeps),]
    newConfLevel <- confPeps[order(-confPeps)[301]]

    # Remove all except 300 most-confident
    confPepIndex <- which(confPeps > newConfLevel)
    data <- data[confPepIndex,]
  }
  # Remove all low-confidence peptides
  rownames <- which(confPeps >= limit)
  data <- data[which(confPeps >= limit),]

  # Return filtered data
  return(data)
}

