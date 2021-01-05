
##################################################################################
# Flow of function calls:
# annotate_peaks()
# |- annotate_ranges()
#    |- splicing_unit_annotations()             ## 5' UTR annotation
#    |- splicing_unit_annotations()             ## 3' UTR region annotations
#    |- splicing_unit_annotations()             ## exons annotations
#    |- splicing_unit_annotations()             ## introns annotations
#    |- region_annotations()                    ## Transcript region annotations
#    |- upstream_annotations()                  ## upstream annotations
#       |- nearest_upstream_bidirectional()
#    |- select_optimal_targets()                ## select correct genes from multiple annotations
#       |- nearest_upstream_bidirectional()
#
#
##################################################################################



#' Annotate narrowPeak using \code{TxDB} object
#'
#' This function annotate the MACS2 called peaks with appropriate target transcript and
#' gene from \code{TxDB} object. Internally uses \code{annotate_ranges()} for Granges
#' annotation. Peaks are annnotated with the following \strong{broad categories} and
#' \emph{specific types} (listed in decreasing order of preference):
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
#' @section Guidelines:
#' Some important observations to do before annotating ChIPseq data:
#' \enumerate{
#' \item Whether the signal is sharp peak (normal TF peaks) or over broader region
#' (polII signal over gene body). Also check if binding is throughout the genome like
#' CTCF factor. See \code{bindingInGene, promoterLength} arguments for the details.
#' \item For the genes which are within peak region, what is the gene size (are
#' genes shorter in length than normal) and how far is the next downstream gene.
#' See \code{includeFractionCut} argument for the details.
#' \item Are there any TES or 3' UTR peaks and how confident are they?
#' \item Check the \code{TXTYPE} column in \code{TxDB} object and see which type of
#' features are of interest to you. Usually tRNA, rRNA are not needed.See
#' \code{excludeType} argument for the details.
#' }
#' These observations will help to decide appropriate parameters while annotating
#' the peaks using TxDB object.
#'
#'
#' @param peakFile A narroPeak or broadPeak file. If a broadPeak file, peak center is
#' used as summit as broadPeak file does not report summit
#' @param fileFormat Format of the peak file. One of "narrowPeak" (Default) or "broadPeak".
#' @param output Optionally store the annotation output to a file
#' @param summitRegion Region width around peak summit to use for annotation purpose.
#' If 0, whole peak region is used. If > 0, \code{summitRegion} bases around peak summit
#' are used in annotation.
#'
#' @inheritParams annotate_ranges
#'
#' @inheritSection annotate_ranges Use of arguments
#'
#' @return A GenomicRanges object with peak annotation
#' @export
#'
#' @examples NA
annotate_peaks <- function(peakFile, fileFormat = "narrowPeak",
                           summitRegion = 0,
                           txdb,
                           promoterLength, upstreamLimit,
                           bidirectionalDistance = 1000, bidirectionalSkew = 0.2,
                           includeFractionCut = 0.7, insideSkewToEndCut = 0.7,
                           txIds = NULL, blacklistRegions = NULL,
                           excludeType = c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA"),
                           bindingInGene = FALSE,
                           removePseudo = FALSE,
                           output = NULL){

  ## started working for peak_annotation on larger genomes
  fileFormat <- match.arg(arg = fileFormat, choices = c("narrowPeak", "broadPeak"))

  ## calculate peak related features
  peaks <- rtracklayer::import(con = peakFile, format = fileFormat)

  if(length(peaks) == 0){
    warning("no peak found in peak file ", basename(peakFile))
    return(NULL)
  }

  if(is.null(mcols(peaks)$peak)){
    mcols(peaks)$peak <- as.integer(width(peaks) / 2)
  }

  mcols(peaks)$peakRegion <- paste(
    as.character(seqnames(peaks)), ":", start(peaks), "-", end(peaks), sep = ""
  )

  ## update the peak region used for the annotation.
  if(summitRegion > 0){
    ## start(peaks) + peaks$peak
    peaks <- GenomicRanges::resize(
      x = GenomicRanges::shift(x = peaks, shift = peaks$peak - summitRegion),
      width = summitRegion, fix = "start"
    )

    ## update the summit
    mcols(peaks)$peak <- as.integer(width(peaks) / 2)
  }

  ## use main function annotate_ranges() for annotation
  peakTargetsGr <- annotate_ranges(
    peaks = peaks, txdb = txdb,
    promoterLength = promoterLength, upstreamLimit = upstreamLimit,
    txIds = txIds, blacklistRegions = blacklistRegions, excludeType = excludeType,
    bidirectionalDistance = bidirectionalDistance,
    includeFractionCut = includeFractionCut, bindingInGene = bindingInGene,
    insideSkewToEndCut = insideSkewToEndCut, removePseudo = removePseudo
  )

  ## rename narrowPeak/broadPeak specific columns and remove unnecessary columns
  peakTargetsGr <- as.data.frame(peakTargetsGr) %>%
    dplyr::rename(
      peakEnrichment = signalValue,
      peakPval = pValue,
      peakQval = qValue,
      peakId = name
    ) %>%
    dplyr::select(-score) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  ## optionally store the data
  if(!is.null(output)){
    readr::write_tsv(x = as.data.frame(mcols(peakTargetsGr)), file = output)
  }

  return(peakTargetsGr)
}


##################################################################################
