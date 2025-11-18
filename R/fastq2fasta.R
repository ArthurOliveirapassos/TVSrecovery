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
