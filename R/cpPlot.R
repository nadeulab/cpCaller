#' Title
#'
#' @description
#'A short description...
#'
#'
#' @param ccObj
#' @param plotType
#' @param purity
#'
#' @import patchwork
#' @import ggplot2
#' @import reshape2
#'
#' @return
#' @export
#'
#' @examples
cpPlot <- function(ccObj, plotType = "simple", purity = NULL){
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
  }
}
