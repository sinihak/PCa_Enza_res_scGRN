#' Load custom gff/gtf annotation
#'
#' @param annotation_table_path a path to gtf file
#' @param include_chromosomes remove chromosomal ranges from the returned object
#' @return GenomicRanges object
#' @import rtracklayer
#' @export
#' @examples
#' getAnnotationTable()
#'
getAnnotationTable <- function(annotation_table_path = "/scratch/svc_td_compbio/tools/Gencode_v43/gencode.v43.annotation.gtf",
                                annotation_style = "UCSC", assembly = "hg38",include_chromosomes=TRUE){
  require(rtracklayer)

  # import tables
  ann <- import.gff(annotation_table_path)
  # fix chromosome notation to chr*
  seqlevelsStyle(ann) <- annotation_style

  # set assembly
  genome(ann) <- assembly
  # ann <- sortSeqLevels(ann)

  #ann <- ann[ann$gene_name %in% rnaseq_genenames, ]
  if (!include_chromosomes){
    ann <- ann[ann$type != "chromosome", ]
  }

  return(ann)
}

