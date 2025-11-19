#' Identify sequences containing canonical telomeric repeats
#'
#' This function searches for canonical telomeric repeats in FASTA files
#' and returns only the sequences that contain a minimum number of consecutive
#' repeat units. It also automatically considers the reverse, complement,
#' and reverse-complement versions of the repeat.
#'
#' @param fasta_file Path to the FASTA file.
#' @param telomere_repeat Canonical telomeric repeat sequence
#'   (example: `"TTAGGG"`).
#' @param min_repeats Minimum number of consecutive repeats required
#'   for the sequence to be considered as containing telomere repeats.
#' @param output_fasta Optional path to the output FASTA file.
#'
#' @return A `DNAStringSet` object containing only the sequences that meet the criteria.
#'
#' @examples
#' \dontrun{
#' hits <- Telomeresfinder("reads.fasta", telomere_repeat = "TTAGGG", min_repeats = 3)
#' }
#'
#' @importFrom Biostrings readDNAStringSet matchPattern width reverseComplement
#' @importFrom dplyr filter
#' @export
Telomeresfinder <- function(fasta_file,
                            telomere_repeat = NULL,
                            min_repeats = NULL,
                            output_fasta = "telomere_hits.fasta") {

  if (missing(fasta_file)) {
    stop("Argumento 'fasta_file' ausente. Forneça um arquivo FASTA.")
  }

  if (is.null(telomere_repeat)) {
    telomere_repeat <- readline(prompt = "Informe a sequência repetitiva: ")
  }

  if (is.null(min_repeats)) {
    min_repeats <- as.numeric(readline(prompt = "Informe o número mínimo de repetições consecutivas: "))
  }

  # Lê o FASTA
  fasta <- Biostrings::readDNAStringSet(fasta_file)

  # Cria o bloco consecutivo (exato)
  repeat_block <- paste0(rep(telomere_repeat, min_repeats), collapse = "")

  # Gera padrões equivalentes
  dna_pat <- Biostrings::DNAString(repeat_block)
  reverse_pat <- Biostrings::reverse(dna_pat)
  complement_pat <- Biostrings::complement(dna_pat)
  revcomp_pat <- Biostrings::reverseComplement(dna_pat)

  # Conta ocorrências exatas de cada padrão
  hits_direct <- Biostrings::vcountPattern(dna_pat, fasta, fixed = TRUE)
  hits_rev <- Biostrings::vcountPattern(reverse_pat, fasta, fixed = TRUE)
  hits_comp <- Biostrings::vcountPattern(complement_pat, fasta, fixed = TRUE)
  hits_revcomp <- Biostrings::vcountPattern(revcomp_pat, fasta, fixed = TRUE)

  # Mantém apenas sequências que possuem PELO MENOS UM bloco exato
  keep <- (hits_direct > 0 |
             hits_rev > 0 |
             hits_comp > 0 |
             hits_revcomp > 0)

  matched <- fasta[keep]

  # Salvar resultado
  if (length(matched) > 0) {
    Biostrings::writeXStringSet(matched, filepath = output_fasta)

    message(
      length(matched),
      " sequências contendo pelo menos 1 bloco com ",
      min_repeats, " repetições consecutivas de '", telomere_repeat,
      "' (ou reverso/complementar) foram salvas em: ",
      output_fasta
    )
  } else {
    message(
      "Nenhuma sequência contendo um bloco de ",
      min_repeats, " repetições consecutivas de '", telomere_repeat,
      "' foi encontrada."
    )
  }

  return(matched)
}



