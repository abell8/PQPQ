#' Create a Pep Object for PQPQ
#'
#' This function creates a pepObject from the provided filename, filters and normalizes it,
#' and then sends it to the relevant PQPQ version.
#' @param filename The file to be used to create a pepObject, in proteinPilot form.
#' @param isOrig Should the original Matlab version be used? Dafaults to FALSE.
#' @return Returns the results from \code{\link{pqpq}} or \code{\link{pqpq_orig}}

createPepObj <- function(filename, isOrig = FALSE){

  ## Get stuff for pepObject and confidence list from readProteinPilot
  rawData <- readProteinPilot(filename) # , isOrig)
  # print(str(rawData$e_meta))

  ## Create pep object
  pepData <- pmartRqc::as.pepData(rawData$e_data, rawData$f_data, rawData$e_meta, edata_cname = "Mass_Tag_ID",
                        fdata_cname = "Sample_ID", emeta_cname = "Protein")
  # print(str(pepData))
  pepDataFixed <- pmartRqc::group_designation(pepData,
                                              main_effects = "Condition", covariates = NULL)

  ## Normalize by hand ##
  # Calculate normalization factors
  colMedians <- apply(pepDataFixed$e_data[,-9], 2, median, na.rm = TRUE)
  aMean <- mean(colMedians)

  # Apply to e_data
  print(pepDataFixed$e_data[,-attr(pepDataFixed, 'cnames')$edata_cname])
  pepDataFixed$e_data[,-9] <- t(t(pepDataFixed$e_data[,-9]) / colMedians) * aMean

  ## Call pqpq or pqpq_orig, as requested
  if (isOrig){
    result <- pqpq_orig(pepDataFixed, rawData$confList)
  } else {
    result <- pqpq(pepDataFixed, rawData$confList)
  }

  return(result)
}


#' Read a Protein Pilot Excel File
#'
#' This function reads information from a file, and (if it is a Protein Pilot .xlsx file)
#' converts it into 3 data frames to form into a pep object.
#' @param filename The name/address of the file to be read
#' @return Returns a list of 3 items. \code{e_data} contains the measurements from samples,
#' \code{f_data} contains some groupings that are nonessential for this problem, and \code{e_meta}
#' contains the linking between proteins and peptides.
readProteinPilot <- function(filename) {


  # Read in protein information from ProteinPilot xlsx file
  rawData <- readxl::read_excel(filename)
  # print(dim(rawData))

  # Remove shared peptides
  AnnotationColumn <- match("Annotation", colnames(rawData))
  sharedPeptides <- which(rawData[,AnnotationColumn] == "auto - shared MS/MS")
  if (length(sharedPeptides > 0)) {
    rawData <- rawData[-sharedPeptides,]
  }

  # Create e_data from Area columns
  n <- grep("Area", colnames(rawData))
  e_data <- as.data.frame(rawData[n])

  e_data$Mass_Tag_ID <- 1:(dim(e_data)[1])

  # Find important columns in data frame
  n <- match(c("Accessions", "Conf", "Sequence"), colnames(rawData))
  pepInfo <- rawData[,n]

  # Create e_meta from pepInfo
  Protein <- pepInfo$Accessions # gsub(";.*", "", pepInfo$Accessions)
  uniqProtein <- unique(Protein)

  Ref_ID <- match(Protein, uniqProtein)
  Peptide_Sequence <- pepInfo$Sequence
  Mass_Tag_ID <- e_data$Mass_Tag_ID
  e_meta <- data.frame(Mass_Tag_ID, Protein, Ref_ID, Peptide_Sequence)


  # Create f_data
  Sample_ID <- colnames(e_data)
  Sample_ID <- Sample_ID[-(length(Sample_ID))]

  Condition <- rep("Infection", length(Sample_ID))
  f_data <- data.frame(Sample_ID, Condition)


  return(list(e_data = e_data, f_data = f_data, e_meta = e_meta, confList = pepInfo$Conf))
}


