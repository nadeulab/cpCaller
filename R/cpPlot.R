#' Plot gene-level copy number alterations
#'
#' @description This function plots gene-level copy number alterations (CNA; i.e., deletions and gains).
#' It takes as input the output of cpCall.
#'
#' @param ccObj ccObj to be used (output of cpCall function).
#' @param plotType Type of plot to accommodate multiple formats. Supported options:
#' \itemize{
#'   \item "simple": standard cpPlot output useful to plot a few genes
#'   \item "genome": option to plot all genes genome-wide
#' }
#' @param purity Path to tab-separated table with the purity of each tumor sample (optional). The table must contain at least two columns labelled "Sample" and "Purity".
#' @param lstGenes Vector with list of genes to be plotted. If not provided, all gene names will be plotted. Only applicable to ploType="simple".
#' @param chrom Allows to specify a chromosome to plot only those genes located in that chromosome. Only applicable to ploType="simple".
#' @param genomeVersion Specify genome version (hg19 or hg38) when plotting genome-wide CN calls.
#'
#' @import patchwork
#' @import ggplot2
#' @import reshape2
#' @import karyoploteR
#'
#' @return This function returns a list of ggplots, one plot per sample analyzed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cpp <- cpPlot(cpc, "simple", purity)
#' }
#'
cpPlot <- function(ccObj, plotType = "simple", purity = NULL, lstGenes = NULL, chrom = NULL, genomeVersion="hg38"){
  message("cpPlot: starting...")
  if(plotType == "simple"){
    message("... plotType = simple...")
    message("1. Preparing inputs and tables...")
    if(!is.null(purity)){ purityDb <- read.table(purity, header = T, sep = "\t", stringsAsFactors = F) }

    countsTumor <- ccObj$table_normalizedCountsTumor
    countsNormal <- ccObj$table_normalizedCountsNormal
    CN <- ccObj$table_CNcalls
    SAM <- ccObj$table_sampleMetrics

    pdb <- melt(countsTumor, id.vars = c(1:6))
    pdb <- rbind(pdb, melt(countsNormal, id.vars=c(1:6)))
    pdb$name_st <- paste0(pdb$gene, "_", pdb$start)
    pdb$NormOrTum[pdb$variable %in% colnames(countsTumor)[7:ncol(countsTumor)]] <- "tumor"
    pdb$NormOrTum[pdb$variable %in% colnames(countsNormal)[7:ncol(countsNormal)]] <- "normal"
    pdb$Call <- "reference"
    for(i in 1:nrow(pdb)){
      if(pdb$NormOrTum[i] == "normal"){ next }
      pdb$Call[i] <- CN$Call[CN$Sample == pdb$variable[i] & CN$Gene == pdb$gene[i]]
    }
    pdb$Call <- factor(pdb$Call, levels=c("reference", "wt", "del", "homodel", "gain", "highgain"))
    pdb <- pdb[nrow(pdb):1,]

    if (!is.null(lstGenes)) { pdb <- pdb[pdb$gene %in% lstGenes,] }
    if (!is.null(chrom)) { pdb <- pdb[pdb$chrom == chrom,] }

    pdb$value[pdb$value > 4] <- 4.02

    colorsCNA <- c("gray91", "#2E8B57", "#DC143C", "#800020", "#4682B4", "#191970")
    names(colorsCNA) <- c("reference", "wt", "del", "homodel", "gain", "highgain")

    message("2. Ploting...")
    pList <- list()
    for(i in 1:length(unique(pdb$variable[pdb$NormOrTum=="tumor"]))){
      sam <- as.character(sort(unique(pdb$variable[pdb$NormOrTum=="tumor"])))[i]
      if(is.null(ccObj$table_tumorAndSelectedNormals)){
        normalsToPlot <- unique(pdb$variable[pdb$NormOrTum == "normal"])
      }else{
        normalsToPlot <- strsplit(ccObj$table_tumorAndSelectedNormals$NormalSamples[ccObj$table_tumorAndSelectedNormals$TumorSample == sam], ";")[[1]]
      }

      p1 <- ggplot(pdb[pdb$variable == sam | (pdb$NormOrTum == "normal" & pdb$variable %in% normalsToPlot),],
                   aes(x = name_st, y = value, group=variable, color=Call)) +
        geom_hline(yintercept = 2, col="gray25", lty=2)+
        geom_hline(yintercept = c(1,3,4), col="gray50", lty=2, lwd=0.2)+
        geom_line() + geom_point() +
        scale_color_manual(values = colorsCNA[names(colorsCNA) %in% pdb$Call[pdb$variable==sam | pdb$NormOrTum=="normal"]] ) +
        theme_classic() +
        ylab("Normalized coverage") + xlab(NULL) +
        ylim(c(0,4.05)) + facet_grid(.~gene, scales = "free_x", space = "free_x") +
        theme(axis.text.x = element_blank()) +
        ggtitle(sam, subtitle = paste0("Purity=", ifelse(is.null(purity), "NA", purityDb$Purity[purityDb$Sample == sam]), " || Normalized coverage=", round(SAM$meanSample[SAM$Sample==sam], 2), ", SD=", round(SAM$sdSample[SAM$Sample==sam], 3), ", Noise=", SAM$Noise[SAM$Sample==sam]))

      pdbSam <- pdb[pdb$variable==sam, ]
      pdbSam$value <- CN$TotalCN[CN$Sample==sam][match(pdbSam$gene, CN$Gene[CN$Sample==sam])]
      pdbSam$value[pdbSam$value > 4] <- 4.02
      p2 <- ggplot(pdbSam, aes(x = name_st, y = value, group=variable, color=Call)) +
        geom_hline(yintercept = 2, col="gray25", lty=2)+
        geom_hline(yintercept = c(1,3,4), col="gray50", lty=2, lwd=0.2)+
        geom_line(lwd=2) +
        scale_color_manual(values = colorsCNA[names(colorsCNA) %in% pdbSam$Call[pdbSam$variable==sam]] ) +
        theme_classic() +
        ylab("Total CN") + xlab(NULL) +
        ylim(c(0,4.05)) + facet_grid(.~gene, scales = "free_x", space = "free_x") +
        # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        theme(axis.text.x = element_blank(), strip.text = element_blank())

      p <- p1 / p2 + plot_layout(heights = c(1.4, 1))
      pList[[i]] <- p
    }
    names(pList) <- as.character(sort(unique(pdb$variable[pdb$NormOrTum=="tumor"])))

    message("Done!")
    return(pList)

  } else if (plotType == "genome"){
    message("... plotType = genome...")
    message("1. Preparing inputs and tables...")
    countsTumor <- ccObj$table_normalizedCountsTumor
    countsNormal <- ccObj$table_normalizedCountsNormal
    CN <- ccObj$table_CNcalls
    SAM <- ccObj$table_sampleMetrics

    pdb <- reshape2::melt(countsTumor, id.vars = c(1:6))
    pdb <- rbind(pdb, reshape2::melt(countsNormal, id.vars=c(1:6)))
    pdb$name_st <- paste0(pdb$gene, "_", pdb$start)
    pdb$NormOrTum[pdb$variable %in% colnames(countsTumor)[7:ncol(countsTumor)]] <- "tumor"
    pdb$NormOrTum[pdb$variable %in% colnames(countsNormal)[7:ncol(countsNormal)]] <- "normal"

    pdb <- merge(
      pdb,
      CN[, c("Sample", "Gene", "Call")],
      by.x = c("variable", "gene"),
      by.y = c("Sample", "Gene"),
      all.x = TRUE
    )
    pdb$Call[is.na(pdb$Call)] <- "reference"
    pdb$Call <- factor(pdb$Call, levels=c("reference", "wt", "NoResult_Noise", "del", "homodel", "gain", "highgain"))
    pdb <- pdb[nrow(pdb):1,]

    colorsCNA <- c("gray91", "#2E8B57", "gray50", "#DC143C", "#800020", "#4682B4", "#191970")
    names(colorsCNA) <- c("reference", "wt", "NoResult_Noise", "del", "homodel", "gain", "highgain")

    message("2. Ploting...")
    pList <- list()
    for(i in 1:length(unique(pdb$variable[pdb$NormOrTum=="tumor"]))){
      sam <- as.character(sort(unique(pdb$variable[pdb$NormOrTum=="tumor"])))[i]
      samplepdb <- pdb[pdb$variable==sam,]
      samplepdb <- samplepdb[!duplicated(samplepdb$gene),]
      cpCaller.gr <- GRanges(seqnames = paste0("chr", samplepdb$chrom),
                              ranges = IRanges(start = samplepdb$start, end = samplepdb$end),
                              sample = samplepdb$variable,
                              CN = as.numeric(samplepdb$value),
                              type = samplepdb$Call, color = colorsCNA[samplepdb$Call])
      markers <- samplepdb[!duplicated(samplepdb[,"gene"]),c("chrom","start","gene")]
      markers$color <- colorsCNA[samplepdb$Call]
      if (!grepl("chr",markers$chrom[1])){
        markers$chrom <- paste0("chr", markers$chrom)
      }
      kp <- plotKaryotype(genome=genomeVersion, main = sam)
      kp <- kpPlotRegions(kp, cpCaller.gr, col=cpCaller.gr$color, r0 = 0, r1 = 0.3)
      kp <- kpPlotMarkers(kp, chr=markers$chrom, x=markers$start, labels=markers$gene, text.orientation = "horizontal",
                          y = 0, r0 = 0.3, label.color = markers$color, line.color = markers$color, adjust.label.position = T, cex=0.7)
      pList[[i]] <- recordPlot()
    }
    names(pList) <- as.character(sort(unique(pdb$variable[pdb$NormOrTum=="tumor"])))
    message("Done!")
    return(pList)
  }
}
