#' Call gene-level copy number alterations
#'
#' @description This function calls gene-level copy number alterations (CNA; i.e., deletions and gains).
#' It takes as input a list of BAM files and a bed-like file with the regions of interest.
#'
#' @param bedFile Path to bed-like file containing the regions/genes to be analyzed. This file must not include header (colnames) and must have 4 columns: chrom, start, end, gene name.
#' @param bamFilesTumors Vector with the paths to the tumor BAM files.
#' @param bamFilesNormals Vector with the paths to the normal BAM files.
#' @param normalization Method use to normalize the counts.
#' Supported options:
#' \itemize{
#'   \item "median_gene": Normalizes the coverage of each target region by the mean coverage of the median coverage of each gene
#'   \item "mean_coverage": Normalizes the coverage of each target region by the mean coverage of the sample
#' }
#' @param lstGenes Vector with list of genes to be used for the analysis. If not provided, all gene names included in the bedFile will be used.
#' @param purity Path to tab-separated table with the purity of each tumor sample (optional). The table must contain at least two columns labelled "Sample" and "Purity".
#' @param selectNormalsByCosine If TRUE, a specified number of normal samples (see below "numOfNormalsByCosine") are selected based on their cosine
#' similarity to each tumor sample and used as reference.
#' @param numOfNormalsByCosine Number of normal samples to be used as a reference, which are selected based on cosine similarity.
#' @param tumorsToUseAsNormal If "bamFilesNormals" is NULL, this parameter defines the number of tumors samples that are used as a reference.
#' @param FFPE If TRUE, noise is force-assigned to "high".
#' @param smooth If TRUE, the coverage of each target region is smoothed taking the mean of the three contiguos regions.
#' @param callOnly If TRUE, steps 1 and 2 to process the bed file and calculate the raw coverage of each region/sample are skipped.
#' @param ccObj ccObj to be re-used if callOnly is TRUE.
#'
#' @import Rsamtools
#' @import progress
#' @import ggplot2
#' @import reshape2
#'
#' @return This function returns an object with the following elements:
#' \itemize{
#'   \item plot_SD_distribution
#'   \item plot_Normalized_coverage
#'   \item plot_SD_samples
#'   \item table_rawCountsNormal
#'   \item table_rawCountsTumor
#'   \item table_normalizedCountsNormal
#'   \item table_normalizedCountsTumor
#'   \item table_tumorAndSelectedNormals
#'   \item table_sampleMetrics
#'   \item table_CNcalls
#' }
#'
#' @details
#' The steps that this functions performs are:
#'  1. Process bed-like file
#'  2. Calculate mean coverage of each target region
#'  3. Normalize coverage per sample and
#'  4. Normalize coverage of each window
#'  5. Smooth coverage (if specified)
#'  6. Estimate noise and call gene-level copy number alterations
#'  7. Draw quality control plots
#'  8. Prepare outputs
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cpc <- cpCall(bedFile, bamFilesTumors, bamFilesNormals, normalization, lstGenes, purity)
#' }
#'
cpCall <- function(bedFile, bamFilesTumors, bamFilesNormals = NULL, normalization = "median_gene", lstGenes = NULL, purity = NULL,
                   selectNormalsByCosine = FALSE, numOfNormalsByCosine = 3, tumorsToUseAsNormal = 3,
                   FFPE = F, smooth = FALSE, callOnly = FALSE, ccObj = NULL){

  message("cpCall: starting...")

  if(isFALSE(callOnly)){

    #
    # process bed file
    #
    message("1. Processing bed file...")
    bed <- read.table(bedFile, header = F, sep = "\t", stringsAsFactors = F)[,1:4]
    colnames(bed) <- c("chrom", "start", "end", "gene")
    bed$size <- bed$end - bed$start + 1
    bed$which_label <- paste0(bed$chrom, ":", bed$start, "-", bed$end)
    ## sort bed file
    bed$chrom[bed$chrom == "X"] <- 23
    bed$chrom[bed$chrom == "Y"] <- 24
    bed$chrom <- as.numeric(as.character(bed$chrom))
    bed$start <- as.numeric(as.character(bed$start))
    bed <- bed[order(bed$chrom, bed$start, decreasing = F), ]
    bed$chrom[bed$chrom == 23] <- "X"
    bed$chrom[bed$chrom == "24"] <- "Y"
    ## subset based on lstGenes to study
    bed <- bed[bed$gene %in% lstGenes, ]
    ## make GRanges
    gr <- makeGRangesFromDataFrame(bed, keep.extra.columns = T)

    #
    # Compute coverage
    #
    message("2. Calculating coverage...")
    # flags and parameters
    flags <- scanBamFlag(isPaired = T, isProperPair = T, isDuplicate = F, isSecondaryAlignment = F, isSupplementaryAlignment = F)
    sb_params <- ScanBamParam(which = gr, flag = flags)
    p_param <- PileupParam(max_depth = 25000, min_base_quality=13, min_mapq=15,
                           min_nucleotide_depth=0, min_minor_allele_depth=0,
                           distinguish_strands=FALSE, distinguish_nucleotides=FALSE,
                           ignore_query_Ns=TRUE, include_deletions=TRUE, include_insertions=FALSE,
                           left_bins=NULL, query_bins=NULL, cycle_bins=NULL)
    # bam files tumor
    message("...starting with tumor bam files...")
    total <- length(bamFilesTumors)
    pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = total)
    pb$tick(0)
    countsTumor <- bed
    for(bam in bamFilesTumors){
      smple <- gsub(".bam", "", strsplit(bam, "/")[[1]][length(strsplit(bam, "/")[[1]])])
      res <- pileup(bam, scanBamParam = sb_params, pileupParam = p_param)
      meanCoverages <- round(sapply(unique(res$which_label), function(x) mean(res$count[res$which_label == x])), 2)
      names(meanCoverages) <- unique(res$which_label)
      countsTumor$newSample <- 0
      countsTumor$newSample[countsTumor$which_label %in% names(meanCoverages)] <- meanCoverages[match(countsTumor$which_label[countsTumor$which_label %in% names(meanCoverages)], names(meanCoverages))]
      colnames(countsTumor)[ncol(countsTumor)] <- smple
      pb$tick(1)
    }

    # bam files normal
    if(!is.null(bamFilesNormals)){
      message("...starting with normal bam files...")
      total <- length(bamFilesNormals)
      pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = total)
      pb$tick(0)
      countsNormal <- bed
      for(bam in bamFilesNormals){
        smple <- gsub(".bam", "", strsplit(bam, "/")[[1]][length(strsplit(bam, "/")[[1]])])
        res <- pileup(bam, scanBamParam = sb_params, pileupParam = p_param)
        meanCoverages <- round(sapply(unique(res$which_label), function(x) mean(res$count[res$which_label == x])), 2)
        names(meanCoverages) <- unique(res$which_label)
        countsNormal$newSample <- 0
        countsNormal$newSample[countsNormal$which_label %in% names(meanCoverages)] <- meanCoverages[match(countsNormal$which_label[countsNormal$which_label %in% names(meanCoverages)], names(meanCoverages))]
        colnames(countsNormal)[ncol(countsNormal)] <- smple
        pb$tick(1)
      }
    } else{ countsNormal <- NULL }

    # save raw counts
    rawCountsNormal <- countsNormal
    rawCountsTumor <- countsTumor

  } else{
    message("1-2. Steps skipped because running in callOnly mode. Re-using raw coverage from ccObj...")
    rawCountsNormal <- ccObj$table_rawCountsNormal
    rawCountsTumor <- ccObj$table_rawCountsTumor
    countsNormal <- rawCountsNormal
    countsTumor <- rawCountsTumor
  }


  #
  # normalize sample by (a) mean coverage or (b) median of mean gene coverage
  #
  message("3. Normalizing coverage of the samples...")
  for(j in 7:ncol(countsTumor)){
    if(normalization == "mean_coverage"){
      countsTumor[,j] <- countsTumor[,j] / mean(countsTumor[,j])
    } else if(normalization == "median_gene"){
      meanGenes <- NULL
      for(gene in unique(countsTumor$gene)){ meanGenes <- c(meanGenes, mean(countsTumor[countsTumor$gene == gene, j])) } # normalized by median of mean gene coverages
      countsTumor[,j] <- countsTumor[,j] / median(meanGenes)
    }
  }

  if(!is.null(bamFilesNormals)){
    for(j in 7:ncol(countsNormal)){
      if(normalization == "mean_coverage"){
        countsNormal[,j] <- countsNormal[,j] / mean(countsNormal[,j])
      } else if(normalization == "median_gene"){
        meanGenes <- NULL
        for(gene in unique(countsNormal$gene)){ meanGenes <- c(meanGenes, mean(countsNormal[countsNormal$gene == gene, j])) } # normalized by median of mean gene coverages
        countsNormal[,j] <- countsNormal[,j] / median(meanGenes)
      }
    }
  }


  #
  # normalize each window by the mean of the normals
  #
  message("4. Normalizing coverage of each window...")
  dbTumorAndReferences <- NULL

  if(!is.null(bamFilesNormals)){ # if normal bam files are available
    meanNormals <- rowMeans(countsNormal[,7:ncol(countsNormal)])
    names(meanNormals) <- countsNormal$gene

    for(j in 7:ncol(countsTumor)){
      if(isTRUE(selectNormalsByCosine)){ # if select 'n' normals by cosine similarity is TRUE
        cosines <- c()
        for(n in 7:ncol(countsNormal)){
          cosines <- c(cosines, sum(countsTumor[,j] * countsNormal[,n]) / (sqrt(sum(countsTumor[,j]^2)) * sqrt(sum(countsNormal[,n]^2))))
        }
        names(cosines) <- colnames(countsNormal)[7:ncol(countsNormal)]
        normalsToUse <- names(cosines)[cosines >= sort(cosines, decreasing = T)[numOfNormalsByCosine]]

        meanNormalsToUse <- rowMeans(countsNormal[,normalsToUse, drop=F])
        names(meanNormalsToUse) <- countsNormal$gene
        dbTumorAndReferences <- rbind(dbTumorAndReferences, c(colnames(countsTumor)[j], paste(normalsToUse, collapse = ";")))

        countsTumor[,j] <- (countsTumor[,j] / meanNormalsToUse) * 2

      } else{ # use all normals
        countsTumor[,j] <- (countsTumor[,j] / meanNormals) * 2
      }
    }

    if(isTRUE(selectNormalsByCosine)){
      colnames(dbTumorAndReferences) <- c("TumorSample", "NormalSamples")
      dbTumorAndReferences <- data.frame(dbTumorAndReferences, stringsAsFactors = F)
    }

    for(i in 7:ncol(countsNormal)){
      countsNormal[,i] <- (countsNormal[,i] / meanNormals) * 2
    }

  } else{ # if no normal samples are provided, normalize by the most "normal" samples within tumor samples available
    meanNormals <- c()
    countsNormal <- countsTumor[,1:6]
    covs_db <- NULL
    for(i in 1:nrow(countsTumor)){
      covs <- as.numeric(countsTumor[i, 7:ncol(countsTumor)])
      dist <- abs(covs - median(covs))
      closest_samples <- order(dist)[1:tumorsToUseAsNormal]
      closest_covs <- covs[closest_samples]
      meanNormals <- c(meanNormals, mean(closest_covs))
      covs_db <- rbind(covs_db, closest_covs)
    }
    for(j in 7:ncol(countsTumor)){
      countsTumor[,j] <- (countsTumor[,j] / meanNormals) * 2
    }

    colnames(covs_db) <- paste0("TumorUsedAsNormal_", 1:ncol(covs_db))
    countsNormal <- cbind(countsNormal, covs_db)
    for(j in 7:ncol(countsNormal)){
      countsNormal[,j] <- (countsNormal[,j] / meanNormals) * 2
    }

  }


  #
  # smooth normalized coverage of tumors in bins of 3 per gene
  #
  if(isTRUE(smooth)){
    message("5. Smoothing coverage of each window...")
    for(j in 7:ncol(countsTumor)){
      for(g in unique(countsTumor$gene)){
        originalVals <- countsTumor[countsTumor$gene == g, j]
        if(length(originalVals) <= 3){ next }
        newVals <- c()
        for(i in 1:length(originalVals)){
          if(i == 1){
            newVal <- mean(originalVals[1:2])
          }else if(i == length(originalVals)){
            newVal <- mean(originalVals[(length(originalVals)-1):length(originalVals)])
          }else{
            newVal <- mean(originalVals[(i-1):(i+1)])
          }
          newVals <- c(newVals, newVal)
        }
        countsTumor[countsTumor$gene == g, j] <- newVals
      }
    }
  } else{
    message("5. Step skipped because running in smooth = FALSE...")
  }


  #
  # call gene-level
  #
  message("6. Estimating noise and calling gene-level copy numbers...")
  CN <- NULL
  SAM <- NULL
  sdb <- melt(countsTumor, id.vars = c(1:6))
  sdbNorm <- melt(countsNormal, id.vars = c(1:6))

  for(smpl in as.character(unique(sdb$variable))){
    CNsample <- NULL
    meanSample <- mean(sdb$value[sdb$variable==smpl])
    sdGenes <- NULL
    for(gene in unique(sdb$gene)){ sdGenes <- c(sdGenes, sd(sdb$value[sdb$variable==smpl & sdb$gene==gene])) } # get sd of each gene
    sdSample <- median(sdGenes, na.rm = T)
    SAM <- rbind(SAM, c(smpl, meanSample, sdSample))
    for(gene in unique(sdb$gene)){
      obs <- sdb$value[sdb$variable==smpl & sdb$gene==gene]
      background <- sdbNorm$value[sdbNorm$gene==gene]
      CN <- rbind(CN, c(smpl, gene,
                        min(obs), quantile(obs, 0.20), mean(obs), median(obs), quantile(obs, 0.80), max(obs), sd(obs),
                        min(background), quantile(background, 0.20), mean(background), median(background), quantile(background, 0.80), max(background), sd(background),
                        log2(min(obs)/min(background)), log2(quantile(obs, 0.20)/quantile(background, 0.20)), log2(mean(obs)/mean(background)), log2(median(obs)/median(background)), log2(quantile(obs, 0.80)/quantile(background, 0.80)), log2(max(obs)/max(background))))
    }
  }
  colnames(SAM) <- c("Sample", "meanSample", "sdSample")
  SAM <- data.frame(SAM, stringsAsFactors = F)
  SAM$meanSample <- as.numeric(as.character(SAM$meanSample))
  SAM$sdSample <- as.numeric(as.character(SAM$sdSample))
  colnames(CN) <- c("Sample", "Gene",
                    "minObs", "q20Obs", "meanObs", "medianObs", "q80Obs", "maxObs", "sdObs",
                    "minBack", "q20Back", "meanBack", "medianBack", "q80Back", "maxBack", "sdBack",
                    "log2Min", "log2q20", "log2Mean", "log2Median", "log2q80", "log2Max")
  CN <- data.frame(CN, stringsAsFactors = F)
  CN[,3:ncol(CN)] <- apply(CN[,3:ncol(CN)], 2 ,as.numeric)

  # assign noise level
  if(isFALSE(FFPE)){
    SAM$Noise <- "low"
    SAM$Noise[SAM$sdSample > 0.2] <- "intermediate"
    SAM$Noise[SAM$sdSample > 0.3] <- "high"
    SAM$Noise <- factor(SAM$Noise, levels = c("low", "intermediate", "high"))
    if(isTRUE(smooth)){ SAM$Noise[SAM$Noise == "low"] <- "intermediate" }
  } else{
    SAM$Noise <- "high"
  }

  CN$Noise <- SAM$Noise[match(CN$Sample, SAM$Sample)]

  # assign Call
  CN$Call <- "wt"
  CN$Call[CN$log2q80 < -0.2 & CN$Noise == "low"] <- "del"
  CN$Call[CN$log2q80 < -1.1 & CN$Noise == "low"] <- "homodel"
  CN$Call[CN$log2q20 > 0.2 & CN$Noise == "low"] <- "gain"
  CN$Call[CN$log2q20 > 0.5 & CN$Noise == "low"] <- "highgain"

  CN$Call[CN$Noise == "intermediate" & CN$log2q80 < -0.3 & CN$q80Obs < 1.8] <- "del"
  CN$Call[CN$Noise == "intermediate" & CN$log2q80 < -1.2 & CN$q80Obs < 1.3] <- "homodel"
  CN$Call[CN$Noise == "intermediate" & CN$log2q20 > 0.3 & CN$q20Obs > 2.2] <- "gain"
  CN$Call[CN$Noise == "intermediate" & CN$log2q20 > 0.6 & CN$q20Obs > 2.7] <- "highgain"

  CN$Call[CN$Noise == "high" &  CN$log2Max < -0.4 & CN$maxObs < 1.3] <- "del"
  CN$Call[CN$Noise == "high" &  CN$log2Max < -1.4 & CN$maxObs < 0.8] <- "homodel"
  CN$Call[CN$Noise == "high" & CN$log2Min > 0.4 & CN$minObs > 2.7] <- "gain"
  CN$Call[CN$Noise == "high" & CN$log2Min > 0.8 & CN$minObs > 3.2] <- "highgain"

  # assign TotalCN, adjust TotalCN by purity (if available), and re-assign Call based on TotalCN
  CN$TotalCN <- CN$medianObs
  CN$TotalCN[CN$Call == "wt"] <- 2
  if(!is.null(purity)){
    purityDb <- read.table(purity, header = T, sep = "\t", stringsAsFactors = F)
    for(i in (1:nrow(CN))[CN$Call != "wt"]){
      sam = CN$Sample[i]
      cnCall = CN$Call[i]
      puritySample = purityDb$Purity[purityDb$Sample == sam]
      if(cnCall == "del"){
        CN$TotalCN[i] <- CN$TotalCN[i] - (1-puritySample)*(2-1)
      }else if(cnCall == "homodel"){
        CN$TotalCN[i] <- CN$TotalCN[i] - (1-puritySample)*2
      } else if(cnCall == "gain"){
        CN$TotalCN[i] <- CN$TotalCN[i] + (1-puritySample)*(2-1)
      } else if(cnCall == "highgain"){
        CN$TotalCN[i] <- CN$TotalCN[i] + (1-puritySample)*2
      }
    }
  }
  CN$TotalCN[CN$Call == "del" & CN$TotalCN < 1 & CN$TotalCN > 0.8] <- 1
  CN$Call[CN$Call == "del" & CN$TotalCN <= 0.8] <- "homodel"
  CN$TotalCN[CN$Call == "homodel" & CN$TotalCN < 0] <- 0
  CN$Call[CN$Call == "highgain" & CN$TotalCN < 3.2] <- "gain"
  CN$TotalCN[CN$Call == "gain" & CN$TotalCN > 3 & CN$TotalCN < 3.2] <- 3
  CN$Call[CN$Call == "gain" & CN$TotalCN >= 3.2] <- "highgain"


  #
  # noise levels and plots
  #
  message("7. Ploting QC plots...")
  # plot density
  p1 <- ggplot(SAM, aes(sdSample)) +
    geom_density(adjust = 1/10, color="steelblue") + theme_classic() +
    geom_vline(xintercept = c(0.2, 0.3, 0.55), col="indianred") +
    ggtitle("Distribution of the mean SD of each sample") +
    xlab("Std Dev (SD)") + ylab("Density")

  # plot normalized coverage
  sdb$SampleType <- "Tumor"
  sdbNorm$SampleType <- "Normal"
  mdb <- rbind(sdb, sdbNorm)
  mdb$Noise <- "low"
  mdb$Noise[mdb$variable %in% SAM$Sample[SAM$Noise=="intermediate"]] <- "intermediate"
  mdb$Noise[mdb$variable %in% SAM$Sample[SAM$Noise=="high"]] <- "high"
  mdb$Noise <- factor(mdb$Noise, levels = c("low", "intermediate", "high"))
  p2 <- ggplot(mdb, aes(variable, value, fill=Noise)) +
    geom_hline(yintercept = 2, col="black") +
    geom_boxplot(outlier.size = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Normalized coverage") + xlab(NULL) +
    facet_grid(.~SampleType, scales = "free_x", space = "free_x") +
    ggtitle("Normalized coverage") +
    scale_fill_manual(values=c("#008080", "#6A5ACD", "#FF6F61"))

  # plot
  p3 <- ggplot(CN, aes(Sample, as.numeric(sdObs), fill=Noise)) +
    geom_hline(yintercept = c(0.2,0.4), col="red4")+
    geom_boxplot(outlier.size = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("SD of each gene") + xlab(NULL) +
    ggtitle("Distribution of SD per sample") +
    scale_fill_manual(values=c("#008080", "#6A5ACD", "#FF6F61"))

  #
  # prepare output object
  #
  message("8. Preparing output object...")
  ccObjOut <- list(p1, p2, p3, rawCountsNormal, rawCountsTumor, countsNormal, countsTumor, dbTumorAndReferences, SAM, CN)
  names(ccObjOut) <- c("plot_SD_distribution", "plot_Normalized_coverage", "plot_SD_samples",
                       "table_rawCountsNormal", "table_rawCountsTumor",
                       "table_normalizedCountsNormal", "table_normalizedCountsTumor", "table_tumorAndSelectedNormals",
                       "table_sampleMetrics", "table_CNcalls")

  message("Done!")
  return(ccObjOut)
}

