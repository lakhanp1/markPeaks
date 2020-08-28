
#' Annotate GRanges using \code{TxDB} object
#'
#' This is a generic function to annotate GRanges object with TxDB annotations. \cr
#' Regions are annnotated with following \strong{broad categories} and \emph{specific
#' types} (listed in decreasing order of preference):
#' \enumerate{
#' \item \strong{featureInPeak:} \emph{"include_tx", "include_CDS"}
#' \item \strong{nearStart:} \emph{"5UTR", "CDS_start", "tx_start"}
#' \item \strong{nearEnd:} \emph{"3UTR", "tx_end", "CDS_end"}
#' \item \strong{peakInFeature:} \emph{"exon", "intron", "inside_tx", "inside_CDS"}
#' \item \strong{upstreamTss:} \emph{"promoter", "upstream"}
#' \item \strong{intergenic:} \emph{"intergenic"}
#' }
#' Additionally, a \emph{pseudo} prefix is added to the peakAnnotation where a peak is
#' annotated to two target genes/features and one of it is more optimum than the other.
#' The less optimum target type is prefixed with \emph{pseudo}. If \code{blacklistRegions}
#' are provided, peak overlapping with the \code{blacklistRegions} are not used for
#' annotation and instead annotated as \emph{"blacklist"}. Please refer to the
#' \strong{Guidelines} section for specific information on this.
#'
#'
#' @param peaks A GRanges object with name column.
#' @param blacklistRegions A BED file or GRanges object with ChIPseq blacklist regions.
#' Peaks overlapping with these regions are not used for annotation.
#' @param removePseudo Logical: whether to remove peak targets which are marked as pseudo.
#' Default: FALSE
#'
#' @inheritParams get_txdb_transcripts_gr
#' @inheritParams region_annotations
#' @inheritParams upstream_annotations
#' @inheritParams select_optimal_targets
#'
#' @inheritSection upstream_annotations Use of arguments
#' @inheritSection select_optimal_targets Use of arguments
#'
#' @return A GenomicRanges object with peak annotation
#' @export
#'
#' @examples NA
annotate_ranges <- function(peaks, txdb, promoterLength, upstreamLimit,
                            bidirectionalDistance = 1000, bidirectionalSkew = 0.2,
                            includeFractionCut = 0.7, insideSkewToEndCut = 0.7,
                            txIds = NULL, blacklistRegions = NULL,
                            excludeType = c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA"),
                            bindingInGene = FALSE,
                            removePseudo = FALSE){


  stopifnot(is(object = peaks, class2 = "GRanges"))
  stopifnot(is(object = txdb, class2 = "TxDb"))

  if(! all(txIds %in% keys(txdb, keytype = "TXID"))){
    stop("unknown TXIDs in txIds: ",
         paste(c(head(txIds[which(!txIds %in% keys(txdb, keytype = "TXID"))]), "..."),
               collapse = " ")
    )
  }

  if(is.null(mcols(peaks)$name)){
    warning("name attribute not found in regions. Using region center...")
    mcols(peaks)$name <- paste("region",1:length(peaks), sep = "_")
  } else if(any(duplicated(mcols(peaks)$name)) ){
    stop("Duplicate values in name column")
  }

  ## rename columns for readability
  mcols(peaks)$peakChr <- as.character(seqnames(peaks))
  mcols(peaks)$peakStart <- start(peaks)
  mcols(peaks)$peakEnd <- end(peaks)


  if(is.null(mcols(peaks)$peak)){
    mcols(peaks)$peak <- as.integer(width(peaks) / 2)
  }

  mcols(peaks)$peakSummit <- GenomicRanges::start(peaks) + mcols(peaks)$peak

  ## ensure that the peak is given as offset from peakStart
  if(any(mcols(peaks)$peakSummit > end(peaks))){
    stop("Peak summit cannot be beyond peakEnd. ",
         "Make sure that peak is a 0 based offset from peakStart (eg. narrowPeak file 10th column)")
  }

  mcols(peaks)$relativeSummitPos <- round((mcols(peaks)$peakSummit - start(peaks)) / width(peaks), 3)


  blacklistPeaks <- NULL

  ## optionally, remove annotation from the peaks which fall in balcklist regions
  if(!is.null(blacklistRegions)){
    if(!is(object = blacklistRegions, class2 = "GRanges")){
      if(!file.exists(blacklistRegions)){
        stop("blacklist regions file not found", blacklistRegions)
      }
      blacklistGr <- rtracklayer::import(con = blacklistRegions, format = "bed")
    }

    olpBlacklist <- GenomicRanges::findOverlaps(
      query = blacklistGr, subject = peaks
    )

    blacklistPeaks <- peaks[olpBlacklist@to]

    ## only use peaks which do not overlap with blacklist regions for annotation purpose
    peaks <- peaks[-olpBlacklist@to]

  }

  ## new environment for global variables which takes time to generate
  ## create once and use multiple times
  # txdbEnv <- new.env(parent = emptyenv())

  assign(x = "pointBasedAnnotation", value = FALSE, envir = txdbEnv)
  if(all(width(peaks) < 10)){
    assign(x = "pointBasedAnnotation", value = TRUE, envir = txdbEnv)
  }

  pointBasedAnnotation <- get(x = "pointBasedAnnotation", envir = txdbEnv)

  transcriptsGr <- get_txdb_transcripts_gr(txdb = txdb, excludeType = excludeType,
                                           tx = txIds)
  txToGene <- get(x = "txToGene", envir = txdbEnv)

  ## 5' UTR annotation
  fiveUtrGr <- get_txdb_fiveUtr_gr(txdb = txdb, tx = transcriptsGr$tx_id)
  fiveUtrTargets <- splicing_unit_annotations(peaksGr = peaks, featuresGr = fiveUtrGr,
                                              featureType = "5UTR", txdb = txdb)

  ## 3' UTR region annotations
  threeUtrGr <- get_txdb_threeUtr_gr(txdb = txdb, tx = transcriptsGr$tx_id)
  threeUtrTargets <- splicing_unit_annotations(peaksGr = peaks, featuresGr = threeUtrGr,
                                               featureType = "3UTR", txdb = txdb)

  ## exons annotations
  exonsGr <- get_txdb_exons_gr(txdb = txdb, tx = transcriptsGr$tx_id)
  exonTargets <- splicing_unit_annotations(peaksGr = peaks, featuresGr = exonsGr,
                                           featureType = "exon", txdb = txdb)

  ## introns annotations
  intronsGr <- get_txdb_introns_gr(txdb = txdb, tx = transcriptsGr$tx_id)
  intronTargets <- splicing_unit_annotations(peaksGr = peaks, featuresGr = intronsGr,
                                             featureType = "intron", txdb = txdb)

  ## Transcript region annotations
  transcriptTargets <- region_annotations(peaksGr = peaks,
                                          featuresGr = transcriptsGr,
                                          includeFractionCut = includeFractionCut,
                                          name = "tx")

  ## annotate upstream targets: IMP to give excludeType so that rRNA, tRNA, snRNAs will be removed
  upstreamTargets <- upstream_annotations(peaksGr = peaks, featuresGr = transcriptsGr,
                                          txdb = txdb, promoterLength = promoterLength,
                                          upstreamLimit = upstreamLimit,
                                          bidirectionalDistance = bidirectionalDistance,
                                          bidirectionalSkew = bidirectionalSkew)

  ## prepare target preference list and peak category list
  ## this is internal preference list
  ## this order is IMP for: select 3UTR between 3UTR and inside_tx as it is more specific
  annotationTypes <- data.frame(
    peakAnnotation = c("include_tx", "include_CDS", "5UTR", "CDS_start", "tx_start", "3UTR",
                       "tx_end", "CDS_end", "EXON", "INTRON", "inside_tx", "inside_CDS",
                       "promoter", "upstream", "pseudo_promoter", "pseudo_upstream"),
    peakPosition = c("TSS", "TSS", "TSS", "TSS", "TSS", "TES",
                     "TES", "TES", "TSS", "TSS", "TSS", "TSS",
                     "TSS", "TSS", "TSS", "TSS"),
    preference = 1:16,
    stringsAsFactors = FALSE)

  ## later these peak categories will be used to decide which type of target to prefer
  peakCategories <- list(
    featureInPeak = c("include_tx", "include_CDS"),
    nearStart = c("5UTR", "CDS_start", "tx_start"),
    nearEnd = c("3UTR", "tx_end", "CDS_end"),
    peakInFeature = c("EXON", "INTRON", "inside_tx", "inside_CDS"),
    upstreamTss = c("promoter", "upstream", "pseudo_promoter", "pseudo_upstream")
  )

  peakCategoryDf <- map_dfr(.x = peakCategories,
                            .f = function(x){data.frame(peakAnnotation = x, stringsAsFactors = F)},
                            .id = "peakCategory")

  annotationTypes <- dplyr::left_join(x = annotationTypes, y = peakCategoryDf, by = "peakAnnotation")

  ## combine annotations
  combinedAnnotations <- NULL
  if(!is.null(fiveUtrTargets)){ combinedAnnotations <- append(combinedAnnotations, fiveUtrTargets) }
  if(!is.null(threeUtrTargets)){ combinedAnnotations <- append(combinedAnnotations, threeUtrTargets) }
  if(!is.null(exonTargets)){ combinedAnnotations <- append(combinedAnnotations, exonTargets) }
  if(!is.null(intronTargets)){ combinedAnnotations <- append(combinedAnnotations, intronTargets) }
  if(!is.null(transcriptTargets)){ combinedAnnotations <- append(combinedAnnotations, transcriptTargets) }
  if(!is.null(upstreamTargets)){ combinedAnnotations <- append(combinedAnnotations, upstreamTargets) }

  if(!is.null(combinedAnnotations)){

    allTargetsDf <- as.data.frame(combinedAnnotations, row.names = NULL) %>%
      dplyr::left_join(y = txToGene, by = c("tx_id" = "TXID")) %>%
      dplyr::left_join(y = annotationTypes, by = c("peakAnnotation" = "peakAnnotation"))

    ## extract best transcript region annotation for each peak-gene combination using the preference
    ## using geneId instead of tx_id: for gene which have multiple tx, 5UTR/3UTR/Intron can get
    ## annotated with any tx
    bestPeakGeneTargets <- allTargetsDf %>%
      dplyr::group_by(name, geneId) %>%
      dplyr::arrange(preference, .by_group = TRUE) %>%
      dplyr::slice(1L) %>%
      dplyr::ungroup()

    bestPeakGeneTargetsGr <- GenomicRanges::makeGRangesFromDataFrame(
      df = bestPeakGeneTargets,
      keep.extra.columns = T)

    #######################
    # ## for summary and debugging
    # tmpDf <- dplyr::group_by(bestPeakGeneTargets, name) %>%
    #   dplyr::summarise(n = n(),
    #                    peakAnnotation = paste(sort(peakAnnotation), collapse = ","),
    #                    peakCat = paste(peakCategory, collapse = ","),
    #                    uniqCat = paste(sort(unique(peakCategory)), collapse = ","),
    #                    gene = paste(geneId, collapse = ",")) %>%
    #   dplyr::filter(n > 1) %>%
    #   dplyr::distinct(peakAnnotation, .keep_all = T) %>%
    #   as.data.frame()
    #######################

    tempTargetGrl <- GenomicRanges::split(x = bestPeakGeneTargetsGr, f = mcols(bestPeakGeneTargetsGr)$name)
    # table(elementNROWS(tempTargetGrl))
    singleHitPeaks <- unlist(tempTargetGrl[which(elementNROWS(tempTargetGrl) == 1)], use.names = F)
    multipleHitPeaksGrl <- tempTargetGrl[which(elementNROWS(tempTargetGrl) != 1)]

    # #######################
    # optGr <- select_optimal_targets(
    #   targetGr = tempTargetGrl$An_kdmB_20h_HA_1_withCtrl_peak_1288,
    #   promoterLength = promoterLength,
    #   bindingInGene = bindingInGene,
    #   insideSkewToEndCut = insideSkewToEndCut)
    # #######################

    if(length(multipleHitPeaksGrl) > 0){
      ## for each peak, find optimum target/s
      multipleHitPeaks <- select_optimal_targets(
        targetGr = sort(unlist(multipleHitPeaksGrl, use.names = F)),
        promoterLength = promoterLength,
        upstreamLimit = upstreamLimit,
        bindingInGene = bindingInGene,
        insideSkewToEndCut = insideSkewToEndCut,
        bidirectionalSkew = bidirectionalSkew,
        bidirectionalDistance = bidirectionalDistance)

    } else{
      multipleHitPeaks <- NULL
    }


    peakTargetsGr <- c(singleHitPeaks, multipleHitPeaks, ignore.mcols=FALSE)

    ## add the unannotated peaks
    unannotatedPeaks <- peaks[which(!peaks$name %in% peakTargetsGr$name)]
    if(length(unannotatedPeaks) > 0){
      mcols(unannotatedPeaks)$peakAnnotation <- "intergenic"
      mcols(unannotatedPeaks)$peakCategory <- "intergenic"
    }

    if(!is.null(blacklistPeaks)){
      mcols(blacklistPeaks)$peakAnnotation <- "blacklist"
      mcols(blacklistPeaks)$peakCategory <- "blacklist"
    }


    peakTargetsGr <- c(peakTargetsGr, unannotatedPeaks, blacklistPeaks, ignore.mcols=FALSE)

    ## set upstream peaks with distance > upstreamLimit to upstream_intergenic
    farUpstreamPeaks <- which(
      peakTargetsGr$peakAnnotation == "upstream" & abs(peakTargetsGr$peakDist) > upstreamLimit
    )

    peakTargetsGr$peakAnnotation[farUpstreamPeaks] <- "upstream_intergenic"
    peakTargetsGr$peakCategory[farUpstreamPeaks] <- "intergenic"


    ## optionally filter peak targets which are marked as pseudo
    if(removePseudo){
      peakTargetsGr <- peakTargetsGr[!grepl(pattern = "pseudo_", x = mcols(peakTargetsGr)$peakAnnotation)]
    }

  } else{
    ## use the original peakset
    peakTargetsGr <- peaks

    mcols(peakTargetsGr)$geneId <- NA
    mcols(peakTargetsGr)$txName <- NA
    mcols(peakTargetsGr)$peakAnnotation <- NA
    mcols(peakTargetsGr)$peakCategory <- NA
    mcols(peakTargetsGr)$peakPosition <- NA
    mcols(peakTargetsGr)$peakDist <- NA
    mcols(peakTargetsGr)$summitDist <- NA
    mcols(peakTargetsGr)$bidirectional <- NA
    mcols(peakTargetsGr)$relativePeakPos <- NA
    mcols(peakTargetsGr)$targetOverlap <- NA
    mcols(peakTargetsGr)$peakOverlap <- NA

  }


  peakTargetsGr <- sort(peakTargetsGr)
  names(x = peakTargetsGr) <- NULL

  ## make sure the peakChr, peakStart and peakEnd columns are selected first
  peakFileCols <- c("peakChr", "peakStart", "peakEnd")
  peakFileCols <- union(peakFileCols, names(mcols(peaks)))

  ## select output colums and rename columns to standard column names:
  annotationGr <- as.data.frame(peakTargetsGr) %>%
    dplyr::select(
      seqnames, start, end, !!! peakFileCols,
      geneId, txName, peakAnnotation, peakCategory, peakPosition, peakDist, summitDist,
      bidirectional, relativePeakPos, targetOverlap, peakOverlap
    ) %>%
    dplyr::rename(
      txId = txName
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  ## remove unnecessary columns from query GRanges
  mcols(annotationGr)$peak <- NULL

  return(annotationGr)
}


##################################################################################
