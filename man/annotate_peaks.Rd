% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p02_annotate_peaks.R
\name{annotate_peaks}
\alias{annotate_peaks}
\title{Annotate narrowPeak using \code{TxDB} object}
\usage{
annotate_peaks(
  peakFile,
  fileFormat = "narrowPeak",
  summitRegion = 0,
  txdb,
  promoterLength,
  upstreamLimit,
  bidirectionalDistance = 1000,
  bidirectionalSkew = 0.2,
  includeFractionCut = 0.7,
  insideSkewToEndCut = 0.7,
  txIds = NULL,
  blacklistRegions = NULL,
  excludeType = c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA"),
  bindingInGene = FALSE,
  removePseudo = FALSE,
  output = NULL
)
}
\arguments{
\item{peakFile}{A narroPeak or broadPeak file. If a broadPeak file, peak center is
used as summit as broadPeak file does not report summit}

\item{fileFormat}{Format of the peak file. One of "narrowPeak" (Default) or "broadPeak".}

\item{summitRegion}{Region width around peak summit to use for annotation purpose.
If 0, whole peak region is used. If > 0, \code{summitRegion} bases around peak summit
are used in annotation.}

\item{txdb}{\code{TxDB} object which will be used for annotation}

\item{promoterLength}{Promoter region length. Upstream peaks within \code{promoterLength}
distance of feature start are annotated as \code{promoter} region peaks.}

\item{upstreamLimit}{Maximum distance of peak for upstream annotation. Peak beyond
this distance can be considered as intergenic instead.}

\item{bidirectionalDistance}{When a peak is present at bidirectional promoter where
distance between two TSS is < \code{bidirectionalDistance}, any gene within
\code{promoterLength} distance of peak is assigned to the peak as annotation. Default: 1000}

\item{bidirectionalSkew}{Maximum fraction of peak region allowed on the side of
false target from the midpoint of two target genes. Default: 0.2}

\item{includeFractionCut}{A floating point number between [0, 1]. If a peak covers more
than this proportion of feature, it is annotated as, eg. \code{include_tx}. Default: 0.7}

\item{insideSkewToEndCut}{A floating point number in range [0, 1]. If a peak is
present inside feature/gene and the relative summit position is >
\code{insideSkewToEndCut}, it is closer to the end of the feature. Default: 0.7}

\item{txIds}{A vector of transcript IDs to be used specifically in the annotation
process instead of full transcript set. These should be internal tx_ids from \code{TxDB}
object. This is useful feature to exclude tRNA, rRNA transcripts while annotating
the regions. Default: NULL}

\item{blacklistRegions}{A BED file or GRanges object with ChIPseq blacklist regions.
Peaks overlapping with these regions are not used for annotation.}

\item{excludeType}{Types of transcripts to exclude from annotation. Should be a
character vector. Default: \code{c("tRNA", "rRNA", "snRNA", "snoRNA", "ncRNA")}.
This argument work only when \code{TxDB} object has \code{TXTYPE} column with
appropriate transcripy type values.}

\item{bindingInGene}{Logical: whether the ChIPseq TF binds in gene body. This is
useful for polII ChIPseq data. Default: FALSE}

\item{removePseudo}{Logical: whether to remove peak targets which are marked as pseudo.
Default: FALSE}

\item{output}{Optionally store the annotation output to a file}
}
\value{
A GenomicRanges object with peak annotation
}
\description{
This function annotate the MACS2 called peaks with appropriate target transcript and
gene from \code{TxDB} object. Internally uses \code{annotate_ranges()} for Granges
annotation. Peaks are annnotated with the following \strong{broad categories} and
\emph{specific types} (listed in decreasing order of preference):
\enumerate{
\item \strong{featureInPeak:} \emph{"include_tx", "include_CDS"}
\item \strong{nearStart:} \emph{"5UTR", "CDS_start", "tx_start"}
\item \strong{nearEnd:} \emph{"3UTR", "tx_end", "CDS_end"}
\item \strong{peakInFeature:} \emph{"exon", "intron", "inside_tx", "inside_CDS"}
\item \strong{upstreamTss:} \emph{"promoter", "upstream"}
\item \strong{intergenic:} \emph{"intergenic"}
}
Additionally, a \emph{pseudo} prefix is added to the peakAnnotation where a peak is
annotated to two target genes/features and one of it is more optimum than the other.
The less optimum target type is prefixed with \emph{pseudo}. If \code{blacklistRegions}
are provided, peak overlapping with the \code{blacklistRegions} are not used for
annotation and instead annotated as \emph{"blacklist"}. Please refer to the
\strong{Guidelines} section for specific information on this.
}
\section{Guidelines}{

Some important observations to do before annotating ChIPseq data:
\enumerate{
\item Whether the signal is sharp peak (normal TF peaks) or over broader region
(polII signal over gene body). Also check if binding is throughout the genome like
CTCF factor. See \code{bindingInGene, promoterLength} arguments for the details.
\item For the genes which are within peak region, what is the gene size (are
genes shorter in length than normal) and how far is the next downstream gene.
See \code{includeFractionCut} argument for the details.
\item Are there any TES or 3' UTR peaks and how confident are they?
\item Check the \code{TXTYPE} column in \code{TxDB} object and see which type of
features are of interest to you. Usually tRNA, rRNA are not needed.See
\code{excludeType} argument for the details.
}
These observations will help to decide appropriate parameters while annotating
the peaks using TxDB object.
}

\section{Use of arguments}{

\subsection{upstreamOverlappingFraction}{
There will be cases when a peak is inside a gene and it is upstream of other gene.
\preformatted{
#                                                                            #
#           target1                     target2                              #
#         =====<=======<===       =====<=======<========<=======             #
#                                   -^--           --^-                      #
#                                 peak1           peak2                      #
#                         |<------>|                                         #
#                         |<------------------------->|                      #
#                                                                            #
In above cases, peak1 can be annotated as Upstream of target1. However not peak2
because target2 has bigger fraction in-between [target1, peak2] range

Target gene inside gene case:
#                                                                            #
#         =====<=======<=======<=======<========<======= target2             #
#          target1 ==<==                                                     #
#                                  -^--           --^-                      #
#                                 peak1           peak2                      #
#                       |<-------->|                                         #
#                                                                            #
}
In above case, peak1 is inside targe2 and upstream of target1. These targets are
selected in \code{select_optimal_targets()} if peak lies within promoter range.
}


\subsection{bidirectionalDistance and bidirectionalSkew}{
For peak at bidirectional promoter, its skewness from midpoint and summit poistion
w.r.t. midpoint is used to decide correct target/s. \cr
\code{nearest_upstream_bidirectional(..., bidirectionalDistance = 1000,
bidirectionalSkew = 0.2)}
\preformatted{
#                        *                                                   #
#                    |<-------more than 500bp------>|                        #
#         target1                    |                    target2            #
#  ==<=====<=====<===                |               ===>===>====>====>==    #
#                              --^---|- peak1                                #
#                                 ---|^--- peak2                             #
#                                  --|^------ peak3                          #
#                peak4  --^--        |                                       #
#                                    |                                       #
#                                  { | } central 10\% region                 #
#                                    |                                       #
#                        midpoint between two targets                        #
#                                                                            #
peak1 => target1: more than 80% of the peak lies on target1 side from TSS midpoint
         and peak summit is not in central 20\% region.
peak2 => target1, targe2: peak lies on the center
peak3 => Additionally, if peak summit (^) is within central 10\% region
         (= bidirectionalSkew/2), peak is assigned to both the target genes
peak4 => target1
If gap between two bidirectional genes < \code{bidirectionalDistance}
(default: 1000bp), peak is assigned to both the target genes
}
}


\subsection{insideSkewToEndCut}{
Use of \code{select_optimal_targets(..., insideSkewToEndCut = 0.7, promoterLength = 500)}
\preformatted{
#                                                                            #
#                  d<500bp                                                   #
#      target1    |<---->|    target2           |<--500-->|      target3     #
#                        0  0.25  0.5  0.75  1                               #
#    ====<====<===       =====>=====>=====>=~~~~           ====>====>===     #
#                         --^--            --^--                             #
#                         peak1            peak2                             #
#                        |<----0.7---->|                                     #
#                        **                                                  #
#                                                                            #
In first example, peak1 is inside target2 and it is near the start of target2.
Even though peak1-target1 distance is < 500bp, target1 is marked as pseudo here.
In another example, peak2 is inside target2 but it is near the end at 3'UTR.
Relative position of the peak2 is >0.7 in target2. peak2 is also upstream of
target3 and within (upstreamLimit)bp. So correct annotations are:
peak1: target2
peak2: target2, target3
}
}
}

\examples{
NA
}
