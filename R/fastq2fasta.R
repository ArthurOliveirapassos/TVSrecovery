#' Convert FASTQ file to FASTA
#'
#' This function reads a FASTQ file and generates a FASTA file with the same
#' IDs and sequences. The output filename can be defined by the user or
#' automatically generated from the input filename.
#'
#' @param file Path to the input FASTQ file.
#' @param output_name Name of the output FASTA file.
#'   If NULL or empty, it will be automatically generated as `"file_converted.fasta"`.
#'
#' @return Invisibly returns the name of the generated file.
#'
#' @examples
#' \dontrun{
#' fastq2fasta("reads.fastq")
#' fastq2fasta("reads.fastq", output_name = "data.fasta")
#' }
#'
#' @importFrom ShortRead readFastq sread id
#' @importFrom Biostrings writeXStringSet
#' @importFrom tools file_path_sans_ext
#' @export
fastq2fasta <- function(file, output_name = NULL) {
  fq <- ShortRead::readFastq(file)
  seqs <- ShortRead::sread(fq)
  names(seqs) <- as.character(ShortRead::id(fq))

  if (is.null(output_name) || output_name == "") {
    base <- tools::file_path_sans_ext(basename(file))
    output_name <- paste0(base, "_converted.fasta")
  }

  Biostrings::writeXStringSet(seqs, filepath = output_name)

  message("Conversion complete! File generated: ", output_name)
  return(invisible(output_name))
}
