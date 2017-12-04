#' PQPQ With Missing Values
#'
#' This function implements the new'n'improved version of PQPQ, with an as-yet-undecided method of
#' avoiding NAs in the distance matrix and allowing more peptides to remain when clustering.
#' @param pepDataFixed An object of type 'pepObject' (from pmartRqc), which has been preferably
#' normalized and filtered. The matrix 'e_data' has n rows.
#' @param confList A list of n confidence values from 0 to 100, corresponding to the peptides in 'e_data'.
#' @return A data frame with 3 columns: \itemize{
#' \item \code{Protein}: Which protein this peptide is associated with
#' \item \code{Peptide}: The row name of the valid peptide
#' \item \code{ID}: Which cluster the peptide is considered clustered to}
#' @keywords Missing values
#' @keywords PQPQ
#'
#' @seealso \code{\link{pqpq_orig}}

pqpq <- function(pepDataFixed, confList){

  ## Find all unique proteins
  ProteinsOfInterest <- as.vector(unique(pepDataFixed[[3]]$Protein))


  ## Allocate memory for later
  ModelPeptideFrame = data.frame(
    "Protein" = character(),
    "ModelPeptide" = character(),
    stringsAsFactors = FALSE)
  SigPeps <- character()
  Corresponding_ID <- numeric()
  Protein <-  character()

  ###### Iterate through the proteins ######
  i = 1
  for (protein in ProteinsOfInterest) {
    print(protein)
    N <- 1

    ## Find peptides associated with the protein
    peptideIndex <- which(pepDataFixed[[3]]$Protein == protein)

    ## Create data frames storing peptide information
    AssociatedPeptides <-
      dplyr::inner_join(pepDataFixed[[3]][peptideIndex, ],
                 pepDataFixed[[1]], by = 'Mass_Tag_ID')

    AbundanceOnly2 <-
      AssociatedPeptides[, 5:length(AssociatedPeptides)]
    rownames(AbundanceOnly2) <- AssociatedPeptides$Mass_Tag_ID
    confPeps <- confList[peptideIndex]

    ## Filter peptides with removeLowConf
    AbundanceOnly <- removeLowConf(data = AbundanceOnly2, confPeps, 95)

    if(plyr::empty(AbundanceOnly)) {
      warning("ERROR: No peptides left after filtration")
      next
    }

    ## Get correlation and p-values for filtered peptides
    CorMat <-
      cor(t(AbundanceOnly), use = "pairwise.complete.obs", method = "pearson")
    p_values <-
      psych::corr.test(
        t(AbundanceOnly),
        use = "pairwise.complete",
        method = "pearson",
        adjust = "none",
        ci = FALSE
      )[[4]]

    ## Choose model peptide candidates based on correlation and p-values
    # Remove correlations that are negative or have low p-values
    CorMat[which(CorMat < 0)] <- NA
    CorMat[which(p_values > 0.4)] <- NA

    # Remove correlations to itself and correlation values of 0
    CorMat <- CorMat - diag(dim(CorMat)[1])
    CorMat[which(CorMat == 0)] <- NA
    Col_All_NAs <-
      sapply(as.data.frame(CorMat), function(x)
        all(is.na(x)))

    # Remove all rows with all NAs
    Model_Peptide_Candidates <- AbundanceOnly[!Col_All_NAs, ]
    Model_Peptide_Candidates <-
      t(removeCols_nNAs(t(Model_Peptide_Candidates), dim(Model_Peptide_Candidates)[2] - 1))

    ## If they're all filtered, default to high-confident peptides
    if(plyr::empty(Model_Peptide_Candidates)) {
      Model_Peptide_Candidates <- AbundanceOnly
    }

    ## Check if more than 50% missing values before clustering
    percentMissing <- 100 * sum(is.na(Model_Peptide_Candidates)) /
      (dim(Model_Peptide_Candidates)[1] *dim(Model_Peptide_Candidates)[2])
    if(percentMissing > 50) {
      warning("More than 50% of model peptide candidate data is missing. Proceed?")
      # For later: Add in behavior for this action
    }
    Num_Model_Peptide_Candidates <- dim(Model_Peptide_Candidates)[1]

    ##### Cluster model peptide candidates, if required #####
    # Per Schematic Figure 1, S1:
    # if NMPC == 1; CorMat's Peptide -> Model Peptide
    # if NMPC <= 3 & NMPC > 1; highest mean intensity peptide -> model peptide
    # if NPMC >= 4; Cluster, and for each cluster, highest mean intensity peptide -> model peptide

    if (Num_Model_Peptide_Candidates >= 4) {

      ## Get model peptides from proteoform_ID_Four_Plus
      toAdd <-
        proteoform_ID_Four_Plus(Model_Peptide_Candidates, AbundanceOnly, protein)
      N <-  dim(toAdd)[1]
    }

    if (Num_Model_Peptide_Candidates >= 2 &
        Num_Model_Peptide_Candidates < 4) {
      toAdd <- c(protein,
                 peptideWithHighestAbundance(Model_Peptide_Candidates))
      toAdd <- matrix(toAdd, nrow = 1, ncol = 2)
      N <- 1
    }

    if (Num_Model_Peptide_Candidates == 1) {
      toAdd <- c(protein, rownames(Model_Peptide_Candidates))
      toAdd <- matrix(toAdd, nrow = 1, ncol = 2)
      N <- 1
    }


    # Check if more than 50% of data missing
    percentMissing <- 100 * sum(is.na(AbundanceOnly2)) /
      (dim(AbundanceOnly2)[1] *dim(AbundanceOnly2)[2])
    print(percentMissing)
    if(percentMissing > 50) {
      warning("More than 50% of peptide data is missing")
      # For later: Add in behavior for this action
    }

    ## Get correlation and p-values for all peptides based on model peptides
    CorMat <-
      cor(t(AbundanceOnly2), use = "pairwise.complete.obs", method = "pearson")

    p_values <-
      psych::corr.test(
        t(AbundanceOnly2),
        use = "pairwise.complete.obs",
        method = "pearson",
        adjust = "none",
        ci = FALSE
      )[[4]]

    ##### Remove clusters based on overlap #####
    ## If only 1 cluster, assign and end
    if (N == 1) {
      newSigPeps <- findSigPeptides(toAdd[1, 2], CorMat, p_values)
      SigPeps <-  c(SigPeps, newSigPeps)
      Corresponding_ID <-
        c(Corresponding_ID, rep(1, length(newSigPeps)))
      Protein <-
        c(Protein, rep(as.character(protein), length(newSigPeps)))
    }

    ## If more than 1, calculate overlap
    else {
      # Create matrix for comparison
      validPeptides <- matrix(0, nrow = dim(pepDataFixed$e_data)[1], ncol = N)

      # Fill with significant peptides
      for(clusterNum in 1:N) {
        newPeptides <- findSigPeptides(toAdd[clusterNum, 2], CorMat, p_values)
        validPeptides[as.integer(newPeptides), clusterNum] <-  as.integer(newPeptides)
      }

      # Go through each column
      for (k in 2:N) {
        # OverlapDegree measures how overlap between column k and l
        overlapDegree <- matrix(0, nrow = 1, ncol = k - 1)

        # Iterate through all previous column
        for (l in 1:(k - 1)) {
          num_peptides <- sum(validPeptides[,k] != 0)
          overlapDegree[l] <- 100 * length(intersect(validPeptides[,k], validPeptides[,l])[-1])/num_peptides
        }
        # If any overlap more than 50%, remove column k
        if(sum(overlapDegree > 50) != 0) {
          validPeptides[,k] <- rep(0, dim(validPeptides)[1])
        }
      } # End overlap for loop

      # Find all columns with nonzero entries
      indexClusters <- which(colSums(validPeptides) != 0)

      # Assign N accordingly
      N <- length(indexClusters)

      ## Assign results to 3 lists
      for(ModPeps in 1:N) {
        sigPepColumn <- validPeptides[,indexClusters[ModPeps]]
        newSigPeps <- which(sigPepColumn != 0)
        # Assign SigPeps, Corresponding_ID, and Protein as in the N = 1 case
        SigPeps <-  c(SigPeps, newSigPeps)
        Corresponding_ID <-
          c(Corresponding_ID, rep(ModPeps, length(newSigPeps)))
        Protein <-
          c(Protein, rep(as.character(protein), length(newSigPeps)))
      }
    }# End multiple cluster else statement

    i <- i + 1
  } # end iteration through proteins

  ## Assign results to data frame
  postisores <-
    data.frame("Protein" = Protein ,
               "Peptide" = SigPeps,
               "ID" = Corresponding_ID)

  postisores <- postisores[!duplicated(postisores[,2:3]),]

  return(postisores)
}
