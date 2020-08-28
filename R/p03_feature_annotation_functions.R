
#' Annotate peaks with 5UTR/exon/intron/3UTR overlap
#'
#' @param peaksGr GRanges object for peak data
#' @param featuresGr A transcript GRanges object
#' @param featureType Feature type. One of \code{c("5UTR", "exon", "intron", "3UTR")}
#' @param txdb TxDB object. Transcript information for the exon/intron overlapping with peak
#' is extracted from TxDB object and used to calculate other features. In case a peak
#' overlap with exons/introns from multiple transcripts of same gene, longest transcript is
#' selected.
#'
#' @return A modified peak GRanges object with additional columns: \code{tx_id,
#' peakAnnotation, targetOverlap, peakOverlap, peakDist, summitDist, relativePeakPos}
#' @export
#'
#' @examples NA
splicing_unit_annotations <- function(peaksGr, featuresGr, featureType, txdb){
  featureType <- match.arg(arg = toupper(featureType),
                           choices = c("5UTR", "EXON", "INTRON", "3UTR"))

  stopifnot(is(object = peaksGr, class2 = "GRanges"))
  stopifnot(is(object = featuresGr, class2 = "GRanges"))
  stopifnot(is(object = txdb, class2 = "TxDb"))

  ovlpHits <- GenomicRanges::findOverlaps(query = peaksGr, subject = featuresGr)

  if(length(ovlpHits) == 0){
    return(NULL)
  }

  peakTargets <- peaksGr[ovlpHits@from]
  mcols(peakTargets)$tx_id <- mcols(featuresGr)$tx_id[ovlpHits@to]
  mcols(peakTargets)$peakAnnotation <- featureType
  mcols(peakTargets)$peakDist <- 0

  ## targetOverlap is at transcript level
  txSubsetGrl <- GenomicFeatures::mapIdsToRanges(
    x = txdb,
    keys = list(tx_id = mcols(peakTargets)$tx_id),
    type = "tx", columns = c("gene_id")
  )

  txSubsetGr <- unlist(txSubsetGrl)

  mcols(peakTargets)$targetStart = start(txSubsetGr)
  mcols(peakTargets)$targetEnd = end(txSubsetGr)
  mcols(peakTargets)$targetStrand = strand(txSubsetGr)
  mcols(peakTargets)$txWidth = width(txSubsetGr)
  mcols(peakTargets)$gene_id = unlist(mcols(txSubsetGr)$gene_id)
  mcols(peakTargets)$targetOverlap <- round(
    width(pintersect(x = peakTargets, y = txSubsetGr)) / width(txSubsetGr), 3
  )
  mcols(peakTargets)$peakOverlap <- round(
    width(pintersect(x = peakTargets, y = txSubsetGr)) / width(peakTargets), 3
  )

  ## calculate relativePeakPos using peak summit position
  ## summitDist > (targetEnd - targetStart) : peak overlap at end and summit is after end
  ## summitDist < 0 : peak overlap at start and summit before start
  targetsDf <- as.data.frame(peakTargets) %>%
    dplyr::group_by(name, gene_id) %>%
    dplyr::arrange(desc(txWidth), .by_group = TRUE) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::select(-gene_id, -txWidth) %>%
    dplyr::mutate(
      summitDist = dplyr::case_when(
        targetStrand == "+" ~ peakSummit - targetStart,
        targetStrand == "-" ~ targetEnd - peakSummit
      )
    ) %>%
    dplyr::mutate(
      relativePeakPos = dplyr::case_when(
        summitDist >= (targetEnd - targetStart) ~ 1,
        summitDist <= 0 ~ 0,
        summitDist > 0 ~ round((summitDist) / (targetEnd - targetStart), 3),
        TRUE ~ 0
      )
    ) %>%
    dplyr::mutate(
      bidirectional = NA,
      relativeSummitPos = dplyr::if_else(
        condition = targetStrand == "-",
        true = round(1 - relativeSummitPos, 3), false = relativeSummitPos)
    )

  if(featureType == "5UTR"){
    ## ensure that for 5UTR region peak, peak does not go beyond target end
    targetsDf <- dplyr::filter(
      targetsDf,
      (targetStrand == "+" & peakEnd < targetEnd) |
        (targetStrand == "-" & peakStart > targetStart)
    )
  } else if(featureType == "3UTR"){
    ## ensure that for 3UTR region peak, peak does not start before target start
    targetsDf <- dplyr::filter(
      targetsDf,
      (targetStrand == "+" & peakStart > targetStart) |
        (targetStrand == "-" & peakEnd < targetEnd)
    )
  }

  targetsGr <- sort(makeGRangesFromDataFrame(df = targetsDf, keep.extra.columns = TRUE))
  return(targetsGr)
}


##################################################################################

#' Map peaks to given GRanges regions
#'
#' This function annotates the peaks onto a regions into e.g. \code{tx_start, tx_end,
#' include_tx, inside_tx} categories.
#' In addition, relativeSummitPos value is updated w.r.t. region for the peaks which are
#' annotated as e.g. \code{inside_tx}
#'
#' @param peaksGr GRanges object generated from narrowPeak or broadPeak file
#' @param featuresGr GRanges object for regions on which peaks needs to be mapped. E.g:
#'  \code{GenomicFeatures::cdsBy()} or \code{GenomicFeatures::transcripts()}
#' @param includeFractionCut A floating point number between [0, 1]. If a peak covers more
#' than this proportion of feature, it is annotated as, eg. \code{include_tx}. Default: 0.7
#' @param name Feature type to be used as suffix in peak type annotation. Eg. CDS, gene etc.
#'
#' @return A modified peak GRanges object with additional columns: \code{tx_id,
#' peakAnnotation, targetOverlap, peakOverlap, peakDist, summitDist, relativePeakPos}
#' @export
#'
#' @examples NA
region_annotations <- function(peaksGr, featuresGr, includeFractionCut = 0.7, name = "tx"){

  name <- match.arg(arg = name, choices = c("gene", "tx", "CDS", "region"))
  stopifnot(is(object = peaksGr, class2 = "GRanges"))
  stopifnot(is(object = featuresGr, class2 = "GRanges"))

  if(includeFractionCut < 0 || includeFractionCut > 1){
    stop("Invalid includeFractionCut value. Should be in range [0, 1]")
  }

  ## define feature types
  insideFeature <- paste("inside_", name, sep = "")
  includeFeature <- paste("include_", name, sep = "")
  overlapStart <- paste(name, "_start", sep = "")
  overlapEnd <- paste(name, "_end", sep = "")


  ovlpHits <- GenomicRanges::findOverlaps(query = peaksGr, subject = featuresGr)

  if(length(ovlpHits) == 0){
    return(NULL)
  }

  queryTargets <- peaksGr[ovlpHits@from]
  mcols(queryTargets)$tx_id <- mcols(featuresGr)$tx_id[ovlpHits@to]
  mcols(queryTargets)$peakAnnotation <- insideFeature
  mcols(queryTargets)$peakDist <- 0

  ##
  mcols(queryTargets)$targetStart = start(featuresGr[ovlpHits@to])
  mcols(queryTargets)$targetEnd = end(featuresGr[ovlpHits@to])
  mcols(queryTargets)$targetStrand = strand(featuresGr[ovlpHits@to])
  mcols(queryTargets)$targetOverlap <- round(
    width(pintersect(x = queryTargets, y = featuresGr[ovlpHits@to])) / width(featuresGr[ovlpHits@to]), 3
  )
  mcols(queryTargets)$peakOverlap <- round(
    width(pintersect(x = queryTargets, y = featuresGr[ovlpHits@to])) / width(queryTargets), 3
  )

  ## assign appropriate relativePeakPos using peak summit (shown as *)
  ##
  ##     |>=====>=====>======>======>======>|           |<=====<=====<======<======<======<|
  ##   -*----- 0                                                                     0 -----*--
  ##                                    -----*- 1     ----*--0.97
  ##                 ---*---0.45                                    ----*--- 0.55
  ##  -------------------------------------------    -------------------------------------------
  ##
  ## calculate relativePeakPos using peak summit position
  ## summitDist > (targetEnd - targetStart) : peak overlap at end and summit is after end
  ## summitDist < 0 : peak overlap at start and summit before start
  targetsDf <- as.data.frame(queryTargets) %>%
    dplyr::mutate(
      peakAnnotation = dplyr::case_when(
        peakStart <= targetStart & peakEnd >= targetEnd ~ includeFeature,
        targetOverlap >= includeFractionCut ~ includeFeature,
        targetStrand == "+" & peakStart <= targetStart & peakEnd < targetEnd ~ overlapStart,
        targetStrand == "+" & targetStart < peakStart & targetEnd <= peakEnd ~ overlapEnd,
        targetStrand == "-" & peakStart <= targetStart & peakEnd < targetEnd ~ overlapEnd,
        targetStrand == "-" & targetStart < peakStart & targetEnd <= peakEnd ~ overlapStart,
        TRUE ~ insideFeature
      )
    ) %>%
    dplyr::mutate(
      summitDist = dplyr::case_when(
        targetStrand == "+" ~ peakSummit - targetStart,
        targetStrand == "-" ~ targetEnd - peakSummit
      )
    ) %>%
    dplyr::mutate(
      relativePeakPos = dplyr::case_when(
        summitDist >= (targetEnd - targetStart) ~ 1,
        summitDist <= 0 ~ 0,
        summitDist > 0 ~ round((summitDist) / (targetEnd - targetStart), 3),
        TRUE ~ 0
      )
    ) %>%
    dplyr::mutate(
      bidirectional = NA,
      relativeSummitPos = dplyr::if_else(
        condition = targetStrand == "-",
        true = round(1 - relativeSummitPos, 3), false = relativeSummitPos)
    )

  ## convert back to GRanges
  queryTargets <- makeGRangesFromDataFrame(df = targetsDf, keep.extra.columns = TRUE)

  return(queryTargets)
}

##################################################################################

#' Set peakAnnotation to pseudo
#'
#' @param target A dataframe or GRanges object which has peakAnnotation column
#'
#' @return Same object with \code{pseudo} prefix to peakAnnotation column values
#' @export
#'
#' @examples NA
set_peakTarget_to_pseudo <- function(target){
  if(any(class(target) %in% "GRanges")){
    mcols(target)$peakAnnotation <- paste("pseudo_", mcols(target)$peakAnnotation, sep = "")
  } else if(any(class(target) %in% "data.frame")){
    target$peakAnnotation <- paste("pseudo_", target$peakAnnotation, sep = "")
  }

  return(target)
}

##################################################################################


#' Annotate upstream peaks on transcripts
#'
#' This function annotates the peaks with nearest downstream target.
#' See \strong{Use of arguments} section for more details.
#'
#' @section Use of arguments:
#' \subsection{upstreamOverlappingFraction}{
#' There will be cases when a peak is inside a gene and it is upstream of other gene.
#' \preformatted{
#' #                                                                            #
#' #           target1                     target2                              #
#' #         =====<=======<===       =====<=======<========<=======             #
#' #                                   -^--           --^-                      #
#' #                                 peak1           peak2                      #
#' #                         |<------>|                                         #
#' #                         |<------------------------->|                      #
#' #                                                                            #
#' In above cases, peak1 can be annotated as Upstream of target1. However not peak2
#' because target2 has bigger fraction in-between [target1, peak2] range
#'
#' Target gene inside gene case:
#' #                                                                            #
#' #         =====<=======<=======<=======<========<======= target2             #
#' #          target1 ==<==                                                     #
#' #                                  -^--           --^-                      #
#' #                                 peak1           peak2                      #
#' #                       |<-------->|                                         #
#' #                                                                            #
#' }
#' In above case, peak1 is inside targe2 and upstream of target1. These targets are
#' selected in \code{select_optimal_targets()} if peak lies within promoter range.
#' }
#'
#'
#' @param peaksGr GRanges object for peak data
#' @param featuresGr A transcript GRanges object
#' @param txdb Optional TxDB object. A 3UTR less region is created using TxDB and used
#' to find the overlap between peak and its downstream target. A max fractional
#' overlap of 0.2 with geneA is allowed in a case when peak overlaps with a geneA and
#' is upstream of geneB. This is useful for the peaks which are near TES of a geneA.
#' If TxDB object is not provided, featuresGr is used. Default: featuresGr is used.
#' @param upstreamOverlappingFraction See details. Default: 0.2
#'
#' @inheritParams nearest_upstream_bidirectional
#'
#' @inheritSection nearest_upstream_bidirectional Use of arguments
#'
#' @return A modified peak GRanges object with additional columns: \code{ tx_id,
#' peakAnnotation, targetOverlap, peakOverlap, peakDist, summitDist, relativePeakPos}
#' @export
#'
#' @examples NA
upstream_annotations <- function(peaksGr, featuresGr, txdb = NULL,
                                 upstreamOverlappingFraction = 0.2,
                                 promoterLength, upstreamLimit,
                                 bidirectionalDistance, bidirectionalSkew){

  stopifnot(is(object = peaksGr, class2 = "GRanges"))
  stopifnot(is(object = featuresGr, class2 = "GRanges"))

  if(exists("pointBasedAnnotation", envir=txdbEnv, inherits=FALSE)) {
    pointBasedAnnotation <- get(x = "pointBasedAnnotation", envir = txdbEnv)
  } else{
    stop("pointBasedAnnotation not found in txdbEnv")
  }

  ## precede(ignore.strand = FALSE) is sufficient to find a peak that preceds a subject
  ## but for bidirectional peaks, precede() will return only nearest feature which is
  ## preceded by peak. However, both the features (+ and - strand) should be returned
  ## Hence, ignore.strand = TRUE is used with both precede() and follow() to extract
  ## features upstream and downstream of the peak

  ## select immediate downstream feature to the peak
  peakDownFeatures <- GenomicRanges::precede(x = peaksGr, subject = featuresGr,
                                             select = "first", ignore.strand = TRUE)

  ## select the upstreamHits which have strand == "+"
  ## because peak is downstream to the feature with strand == "-"
  peakDownHits <- data.frame(
    from = 1:length(peaksGr), peakId = mcols(peaksGr)$name,
    to = peakDownFeatures,
    expectedFeatureStrand = "+",
    featureStrand = as.vector(strand(featuresGr))[peakDownFeatures],
    txName = featuresGr$tx_name[peakDownFeatures],
    stringsAsFactors = FALSE)

  ## select immediate upstream feature to the peak
  peakUpFeatures <- GenomicRanges::follow(x = peaksGr, subject = featuresGr,
                                          select = "last", ignore.strand = TRUE)

  ## select the upstreamHits which have strand == "-"
  ## because peak is downstream to the feature with strand == "+"
  peakUpHits <- data.frame(
    from = 1:length(peaksGr), peakId = mcols(peaksGr)$name,
    to = peakUpFeatures,
    expectedFeatureStrand = "-",
    featureStrand = as.vector(strand(featuresGr))[peakUpFeatures],
    txName = featuresGr$tx_name[peakUpFeatures],
    stringsAsFactors = FALSE)


  targetsAroundPeak <- dplyr::bind_rows(peakDownHits, peakUpHits) %>%
    dplyr::filter(!is.na(to))

  ## remove false downstream targets
  upstreamHits <- dplyr::filter(targetsAroundPeak, expectedFeatureStrand == featureStrand)

  ####################
  ## if only one gene is found i.e. either upstream of peak or downstream of peak,
  ## additionally select one immediate downstream gene in opposite direction.
  ## In following case, gene0 is found by follow() method which is upstream of peak1
  ## and is possible target. Gene2 can also be a possible target but precede() will
  ## select gene1 (ignore.strand = TRUE) and it will be filtered ultimately.
  ## So we select on gene in opposite direction with ignore.strand = TRUE. This way
  ## we will have two targets on each side of peak. if its case like gene1 and gene3,
  ## gene3 will be filtered using peak-gap overlap filter code below
  ##                           ====<=======<======<====== gene1
  ##                             ====>=====>======>======= gene2
  ##  ====<======<=                                    ====>======> gene3
  ##                   -----
  ##     gene0            peak1
  ##

  ## find if any other target overlap with false downstream target (gene1) which can
  ## be possible true downstream target
  falseTargets <- dplyr::filter(targetsAroundPeak, expectedFeatureStrand != featureStrand)

  if (nrow(falseTargets) > 0) {

    falseTargetsGr <- featuresGr[falseTargets$to]
    mcols(falseTargetsGr)$from  <- falseTargets$from
    mcols(falseTargetsGr)$to <- falseTargets$to
    mcols(falseTargetsGr)$peakId <- falseTargets$peakId
    mcols(falseTargetsGr)$expectedFeatureStrand <- falseTargets$expectedFeatureStrand

    falseTargetOverlaps <- GenomicRanges::findOverlaps(query = falseTargetsGr,
                                                       subject = featuresGr,
                                                       ignore.strand = TRUE)

    otherUpstream <- tibble::tibble(
      from = mcols(falseTargetsGr)$from[falseTargetOverlaps@from],
      peakId =  mcols(falseTargetsGr)$peakId[falseTargetOverlaps@from],
      to = falseTargetOverlaps@to,
      expectedFeatureStrand = mcols(falseTargetsGr)$expectedFeatureStrand[falseTargetOverlaps@from],
      featureStrand = as.vector(strand(featuresGr))[falseTargetOverlaps@to],
      txName = featuresGr$tx_name[falseTargetOverlaps@to],
      stringsAsFactors = FALSE) %>%
      dplyr::filter(expectedFeatureStrand == featureStrand)

    ## some of above otherUpstream target can overlap with peak. remove such targets
    ## as it is handled by splicing_unit_annotations() for 5' UTR overlap
    notOvlpWithPeak <- distance(x = peaksGr[otherUpstream$from], y = featuresGr[otherUpstream$to],
                                ignore.strand = TRUE) != 0

    otherUpstream <- otherUpstream[notOvlpWithPeak, ]

  } else{
    otherUpstream <- NULL
  }
  ####################


  ## combine putative upstream hits
  upstreamHits <- dplyr::bind_rows(upstreamHits, otherUpstream) %>%
    dplyr::arrange(from) %>%
    dplyr::mutate(hitId = row_number())

  ## return null if nothing found
  if(nrow(upstreamHits) == 0){
    return(NULL)
  }

  ## find the number of genes between peak and its target.
  ## only those targets are true where there is no other feature inbetween
  ## build a GRanges object of the gap region between peak and target gene

  ## generate peak-target gap GRanges
  peakTargetGapsGr <- GenomicRanges::pgap(x = peaksGr[upstreamHits$from],
                                          y = featuresGr[upstreamHits$to])

  peakTargetGapsGr <- unstrand(peakTargetGapsGr)
  names(peakTargetGapsGr) <- upstreamHits$hitId
  mcols(peakTargetGapsGr)$hitId <- upstreamHits$hitId

  ## build a subject GRanges for tx - 3UTR
  ## Such custom regions are used because genes have very long 3' UTR.
  ## A peak in UTR region of a gene can be upstream of another gene
  featureIdx <- tibble::tibble(
    tx_id = mcols(featuresGr)$tx_id, featureIndex = 1:length(featuresGr)
  )

  if(!is.null(txdb)){
    stopifnot(is(object = txdb, class2 = "TxDb"))

    fiveUtrGr <- get_txdb_fiveUtr_gr(txdb = txdb)
    threeUtrGr <- get_txdb_threeUtr_gr(txdb = txdb)

    threeUtrIdx <- tibble::tibble(
      tx_id = mcols(threeUtrGr)$tx_id,
      threeUtrIndex = 1:length(threeUtrGr)
    ) %>%
      dplyr::left_join(y = featureIdx, by = "tx_id")

    txMinus3Utrs <- GenomicRanges::psetdiff(
      x = featuresGr[threeUtrIdx$featureIndex],
      y = threeUtrGr[threeUtrIdx$threeUtrIndex],
      ignore.strand= TRUE)

    txMinusUtrs <- sort(
      c(txMinus3Utrs,
        featuresGr[setdiff(featureIdx$featureIndex, threeUtrIdx$featureIndex)],
        ignore.mcols = TRUE)
    )

    ## no need to remove 5UTR
  }

  ## find overlapping gene/features in gap GRanges.
  featuresInGap <- GenomicRanges::findOverlaps(query = peakTargetGapsGr,
                                               subject = txMinusUtrs,
                                               ignore.strand = TRUE)

  isFeatureInBetweenDf <- tibble::tibble(
    hitId = mcols(peakTargetGapsGr)$hitId[featuresInGap@from],
    gapWidth = width(peakTargetGapsGr)[featuresInGap@from],
    gapGrRow = featuresInGap@from,
    firstOverlapFeature = featuresInGap@to
  ) %>%
    dplyr::filter(!is.na(firstOverlapFeature))

  ## calculate fraction overlap: needed for the peak which is inside a gene
  isFeatureInBetweenDf$intersectWd <- width(
    GenomicRanges::pintersect(
      x = peakTargetGapsGr[isFeatureInBetweenDf$gapGrRow],
      y = txMinusUtrs[isFeatureInBetweenDf$firstOverlapFeature],
      ignore.strand = TRUE)
  )

  isFeatureInBetweenDf$ovlpFeatureWd <- width(txMinusUtrs[isFeatureInBetweenDf$firstOverlapFeature])
  isFeatureInBetweenDf$ovlpFeatureFrac <- isFeatureInBetweenDf$intersectWd / isFeatureInBetweenDf$ovlpFeatureWd
  isFeatureInBetweenDf$ovlpGapFrac <- isFeatureInBetweenDf$intersectWd / isFeatureInBetweenDf$gapWidth

  isFeatureInBetweenDf <- dplyr::group_by(isFeatureInBetweenDf, hitId) %>%
    dplyr::arrange(desc(ovlpFeatureFrac), .by_group = TRUE) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup()

  ## upstreamOverlappingFraction based filtering OR
  ## select if target is within promoterLength distance
  ## 0.2 is still very big for the large genomes such as human, mouse as genes are very long
  upstreamHitsFiltered <- dplyr::left_join(x = upstreamHits, y = isFeatureInBetweenDf, by = "hitId") %>%
    tidyr::replace_na(list(intersectWd = 0, ovlpFeatureWd = 0, ovlpFeatureFrac = 0, ovlpGapFrac = 0)) %>%
    dplyr::filter(ovlpFeatureFrac <= upstreamOverlappingFraction | gapWidth < promoterLength)

  ## if no upstream peak after filtering
  if(nrow(upstreamHitsFiltered) == 0){
    return(NULL)
  }

  ## build upstream peaks data and filter unnecessary peaks where there is/are genes between peak and target
  upstreamPeaks <- peaksGr[upstreamHitsFiltered$from]

  mcols(upstreamPeaks)$tx_id <- mcols(featuresGr)$tx_id[upstreamHitsFiltered$to]
  mcols(upstreamPeaks)$peakAnnotation <- "upstream"
  mcols(upstreamPeaks)$peakDist <- GenomicRanges::distance(
    x = upstreamPeaks,
    y = featuresGr[upstreamHitsFiltered$to]
  )

  mcols(upstreamPeaks)$peakDist <- mcols(upstreamPeaks)$peakDist * -1
  mcols(upstreamPeaks)$targetStart = start(featuresGr[upstreamHitsFiltered$to])
  mcols(upstreamPeaks)$targetEnd = end(featuresGr[upstreamHitsFiltered$to])
  mcols(upstreamPeaks)$targetStrand = strand(featuresGr[upstreamHitsFiltered$to])

  ## calculate summit distance and change relativeSummitPos based on target gene
  upstreamPeaks <- as.data.frame(upstreamPeaks, stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      targetOverlap = 0,
      peakOverlap = 0,
      summitDist = dplyr::case_when(
        targetStrand == "+" ~ peakSummit - targetStart,
        targetStrand == "-" ~ targetEnd - peakSummit
      )
    ) %>%
    dplyr::mutate(
      relativePeakPos = 0,
      bidirectional = NA,
      relativeSummitPos = dplyr::if_else(
        condition = targetStrand == "-",
        true = round(1 - relativeSummitPos, 3), false = relativeSummitPos),
      peakAnnotation = if_else(abs(peakDist) < promoterLength, "promoter", peakAnnotation)
    ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)


  ######
  ## select correct target/s in case of bi-directional peak
  tempTargetGrl <- GenomicRanges::split(x = upstreamPeaks, f = mcols(upstreamPeaks)$name)
  singalTargetPeaks <- tempTargetGrl[which(elementNROWS(tempTargetGrl) == 1)]
  dualTargetPeaks <- tempTargetGrl[which(elementNROWS(tempTargetGrl) > 1)]

  if(any(elementNROWS(dualTargetPeaks) > 2)){
    stop("More than two upstream targets targets found for peaks")
  }

  dualTargetFiltered <- NULL

  if(length(dualTargetPeaks) > 0){

    dualTargetPeaksDf <- sort(unlist(dualTargetPeaks, use.names = FALSE)) %>%
      as.data.frame(row.names = NULL, stringsAsFactors = FALSE) %>%
      dplyr::mutate(rowIdx = 1:n(),
                    target = TRUE)

    ## row index table for bidirectional peak
    pairTable <- dualTargetPeaksDf %>%
      dplyr::group_by(name) %>%
      dplyr::mutate(group = 1:n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(group = paste("t", group, sep = "")) %>%
      dplyr::select(name, group, rowIdx) %>%
      tidyr::spread(key = group, value = rowIdx)

    pseudoUpIdx <- nearest_upstream_bidirectional(
      targetDf = dualTargetPeaksDf,
      t1Idx = pairTable$t1,
      t2Idx = pairTable$t2,
      promoterLength = promoterLength,
      upstreamLimit = upstreamLimit,
      bidirectionalDistance = bidirectionalDistance,
      pointBasedAnnotation = pointBasedAnnotation)

    correctUpstream <- setdiff(x = c(pairTable$t1, pairTable$t2), y = pseudoUpIdx)

    bidirectAIdx <- purrr::map2(
      .x = pairTable$t1, .y = pairTable$t2,
      .f = function(x, y){
        if(any(c(x, y) %in% pseudoUpIdx)){
          return(NULL)
        } else{
          return(c(x, y))
        }
      }) %>%
      purrr::flatten_int()

    dualTargetPeaksDf$bidirectional[bidirectAIdx] <- "A"

    ## remove the targets which are too far based on nearest_upstream_bidirectional()
    dualTargetFiltered <- dualTargetPeaksDf[correctUpstream, ] %>%
      dplyr::select(-rowIdx, -target) %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  }


  upstreamPeaksAn <- sort(c(unlist(singalTargetPeaks, use.names = FALSE),
                            dualTargetFiltered))
  names(upstreamPeaksAn) <- NULL

  return(upstreamPeaksAn)

}


##################################################################################



#' True target for bidirectional peak
#'
#' In case of peak at bidirectional promoter, this function assigns one or both genes
#' as peak target. See \strong{Use of arguments} section for more details.
#'
#' @section Use of arguments:
#' \subsection{bidirectionalDistance and bidirectionalSkew}{
#' For peak at bidirectional promoter, its skewness from midpoint and summit poistion
#' w.r.t. midpoint is used to decide correct target/s. \cr
#' \code{nearest_upstream_bidirectional(..., bidirectionalDistance = 1000,
#' bidirectionalSkew = 0.2)}
#' \preformatted{
#' #                        *                                                   #
#' #                    |<-------more than 500bp------>|                        #
#' #         target1                    |                    target2            #
#' #  ==<=====<=====<===                |               ===>===>====>====>==    #
#' #                              --^---|- peak1                                #
#' #                                 ---|^--- peak2                             #
#' #                                  --|^------ peak3                          #
#' #                peak4  --^--        |                                       #
#' #                                    |                                       #
#' #                                  { | } central 10\% region                 #
#' #                                    |                                       #
#' #                        midpoint between two targets                        #
#' #                                                                            #
#' peak1 => target1: more than 80% of the peak lies on target1 side from TSS midpoint
#'          and peak summit is not in central 20\% region.
#' peak2 => target1, targe2: peak lies on the center
#' peak3 => Additionally, if peak summit (^) is within central 10\% region
#'          (= bidirectionalSkew/2), peak is assigned to both the target genes
#' peak4 => target1
#' If gap between two bidirectional genes < \code{bidirectionalDistance}
#' (default: 1000bp), peak is assigned to both the target genes
#' }
#' }
#'
#'
#' @param targetDf A dataframe which has bidirectional targets for each peak.
#' @param t1Idx target1 index vector
#' @param t2Idx target2 index vector. IMP: \code{length(t1Idx)} should be equal to
#' \code{length(t2Idx)}
#' @param bidirectionalSkew Maximum fraction of peak region allowed on the side of
#' false target from the midpoint of two target genes. Default: 0.2
#' @param bidirectionalDistance When a peak is present at bidirectional promoter where
#' distance between two TSS is < \code{bidirectionalDistance}, any gene within
#' \code{promoterLength} distance of peak is assigned to the peak as annotation. Default: 1000
#' @param promoterLength Promoter region length. Upstream peaks within \code{promoterLength}
#' distance of feature start are annotated as \code{promoter} region peaks.
#' @param upstreamLimit Maximum distance of peak for upstream annotation. Peak beyond
#' this distance can be considered as intergenic instead.
#' @param pointBasedAnnotation Logical: whether peak annotation is based on just
#' the summit or whole peak region.
#'
#' @return A vector of row index for pseudo targets.
#' @export
#'
#' @examples NA
nearest_upstream_bidirectional <- function(targetDf, t1Idx, t2Idx,
                                           promoterLength, upstreamLimit,
                                           bidirectionalSkew = 0.2,
                                           bidirectionalDistance = 1000,
                                           pointBasedAnnotation = FALSE){

  targetPairDf <- tibble::tibble(t1Idx = t1Idx, t2Idx = t2Idx,
                                 t1Select = TRUE, t2Select = TRUE,
                                 dir = "opposite")

  targetA <- targetDf[targetPairDf$t1Idx, ]
  targetB <- targetDf[targetPairDf$t2Idx, ]

  if(!all(targetA$name == targetB$name)){
    stop("peak name mismatch in bidirectional target pair")
  }

  targetPairDf$peakId <- targetA$name
  targetPairDf$peakSummit <- targetA$peakSummit
  targetPairDf$t1PeakDist <- targetA$peakDist
  targetPairDf$t2PeakDist <- targetB$peakDist
  targetPairDf$t1SummitDist <- targetA$summitDist
  targetPairDf$t2SummitDist <- targetB$summitDist

  peakGr <- GenomicRanges::makeGRangesFromDataFrame(df = targetA)
  targetPairDf$peakWidth <- width(peakGr)
  targetPairDf$peakFraction <- round(targetPairDf$peakWidth * bidirectionalSkew)

  sameDirTargets <- which(targetA$targetStrand == targetB$targetStrand)
  targetPairDf$dir[sameDirTargets] <- "same"


  targetAGr <- dplyr::select(targetA, seqnames, targetStart, targetEnd, targetStrand, name, rowIdx) %>%
    makeGRangesFromDataFrame(start.field = "targetStart", end.field = "targetEnd",
                             strand.field = "targetStrand",
                             keep.extra.columns = TRUE)

  targetBGr <- dplyr::select(targetB, seqnames, targetStart, targetEnd, targetStrand, name, rowIdx) %>%
    makeGRangesFromDataFrame(start.field = "targetStart", end.field = "targetEnd",
                             strand.field = "targetStrand",
                             keep.extra.columns = TRUE)

  zeroGapTargets <- which(distance(x = targetAGr, y = targetBGr, ignore.strand = TRUE) == 0)
  # if(any(distance(x = targetAGr, y = targetBGr, ignore.strand = TRUE) == 0)){
  #   stop("Distance between bidirectional target genes cannot be 0")
  # }

  ## gap between two bidirectional targets
  targetGapGr <- unstrand(pgap(x = targetAGr, y = targetBGr, ignore.strand = TRUE))
  targetPairDf$gapWidth <- width(targetGapGr)
  targetPairDf$midpoint <- start(targetGapGr) + round(width(targetGapGr)/2)

  ## midpoint of gap between bidirectional targets
  # midpoint <- resize(x = targetGapGr, width = 1, fix = "center", ignore.strand = TRUE)
  # targetPairDf$midpointDist <- distance(x = midpoint, y = peakGr, ignore.strand = TRUE)

  centralNoSummitZone <- bidirectionalSkew
  targetPairDf$summitPosLim <- (targetPairDf$gapWidth + targetPairDf$gapWidth*centralNoSummitZone)/2

  ## gapWidth condition has to be 1st always to ensure that the targets which are
  ## very close are not marked pseudo
  targetPairDf <- targetPairDf %>% dplyr::mutate(
    t1Select = dplyr::case_when(
      gapWidth <= bidirectionalDistance & abs(t1PeakDist) < promoterLength ~ TRUE,
      # pointBasedAnnotation == FALSE & abs(t1PeakDist) < promoterLength ~ TRUE,
      abs(t1PeakDist) > upstreamLimit & abs(t2PeakDist) <= promoterLength ~ FALSE,
      dir == "opposite" & pointBasedAnnotation == FALSE &
        (abs(t1PeakDist) + peakFraction) > (gapWidth/2) ~ FALSE,
      dir == "opposite" & pointBasedAnnotation &
        abs(t1SummitDist) > summitPosLim ~ FALSE,
      TRUE ~ t1Select
    ),
    t2Select = dplyr::case_when(
      gapWidth <= bidirectionalDistance & abs(t2PeakDist) < promoterLength ~ TRUE,
      # pointBasedAnnotation == FALSE & abs(t2PeakDist) < promoterLength ~ TRUE,
      abs(t2PeakDist) > upstreamLimit & abs(t1PeakDist) <= promoterLength ~ FALSE,
      dir == "opposite" & pointBasedAnnotation == FALSE &
        (abs(t2PeakDist) + peakFraction) > (gapWidth/2) ~ FALSE,
      dir == "opposite" & pointBasedAnnotation &
        abs(t2SummitDist) > summitPosLim ~ FALSE,
      TRUE ~ t2Select
    )
  )

  falseUpIdx <- c(targetPairDf$t1Idx[which(!targetPairDf$t1Select)],
                  targetPairDf$t2Idx[which(!targetPairDf$t2Select)])

  return(falseUpIdx)
}


##################################################################################


