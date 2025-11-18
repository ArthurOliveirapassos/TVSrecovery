#' Converter arquivo FASTQ para FASTA
#'
#' Esta função lê um arquivo FASTQ e gera um arquivo FASTA com os mesmos IDs
#' e sequências. O nome de saída pode ser definido pelo usuário ou gerado
#' automaticamente a partir do nome do arquivo de entrada.
#'
#' @param file Caminho para o arquivo FASTQ de entrada.
#' @param output_name Nome do arquivo FASTA de saída.
#'   Se NULL ou vazio, será gerado automaticamente como `"arquivo_converted.fasta"`.
#'
#' @return Invisivelmente retorna o nome do arquivo gerado.
#'
#' @examples
#' \dontrun{
#' fastq2fasta("reads.fastq")
#' fastq2fasta("reads.fastq", output_name = "dados.fasta")
#' }
#' @import Biostrings
#' @importFrom ShortRead readFastq sread id
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
