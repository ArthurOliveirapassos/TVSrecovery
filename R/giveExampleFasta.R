#' Gerar um arquivo FASTA de exemplo com sequências teloméricas e aleatórias
#'
#' Esta função cria um arquivo FASTA contendo um conjunto misto de sequências:
#' parte delas contém repetições teloméricas canônicas (com possibilidade de
#' variantes), enquanto o restante é composto por sequências totalmente aleatórias.
#'
#' O objetivo é fornecer um dataset sintético útil para testes, demonstrações
#' e exemplos no próprio pacote, como nos chunks do README ou vinhetas.
#'
#' @param output Nome do arquivo FASTA de saída.
#'   Padrão: `"example_sequences.fasta"`.
#' @param n_seqs Número total de sequências a serem geradas.
#'   Padrão: 100.
#' @param telomere_repeat Repetição telomérica canônica a ser utilizada.
#'   Padrão: `"TTAGGG"`.
#' @param min_repeats Número mínimo de repetições consecutivas por sequência.
#' @param max_repeats Número máximo de repetições consecutivas por sequência.
#' @param variant_prob Probabilidade de gerar uma variante dentro das repetições.
#'
#' @details
#' A função cria aleatoriamente:
#'
#' * Entre **40 e 60 sequências teloméricas** contendo repetições canônicas
#'   e variantes introduzidas aleatoriamente.
#' * O restante das sequências (`n_seqs - teloméricas`) será **totalmente aleatório**.
#'
#' A função também adiciona flancos aleatórios em volta do bloco telomérico.
#'
#' Variantes são geradas substituindo uma única base dentro da repetição canônica.
#'
#' @return Um objeto `DNAStringSet` contendo todas as sequências geradas
#'   (também salvas em arquivo FASTA).
#'
#' @examples
#' \dontrun{
#' # Gerar FASTA de exemplo no diretório atual
#' example_data <- giveExampleFasta(
#'   output = "example.fasta",
#'   n_seqs = 50,
#'   telomere_repeat = "TTAGGG",
#'   variant_prob = 0.2
#' )
#'
#' # Ver as primeiras sequências
#' example_data[1:5]
#' }
#'
#' @import Biostrings
#' @export
giveExampleFasta <- function(
    output = "example_sequences.fasta",
    n_seqs = 100,
    telomere_repeat = "TTAGGG",
    min_repeats = 3,
    max_repeats = 8,
    variant_prob = 0.3  # probabilidade de criar variantes
) {

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Pacote Biostrings necessário. Instale com: BiocManager::install('Biostrings')")
  }

  # ---- Número de sequências com telômeros ----
  set.seed(123)
  n_telomeric <- sample(40:60, 1)         # entre 40 e 60 teloméricas
  n_random <- n_seqs - n_telomeric

  # ---- Função auxiliar para criar variantes ----
  create_variant <- function(rep_seq) {
    bases <- c("A", "C", "G", "T")

    # posição a modificar
    pos <- sample(1:nchar(rep_seq), 1)

    original <- substr(rep_seq, pos, pos)

    # troca por uma base diferente
    new_base <- sample(bases[bases != original], 1)

    substr(rep_seq, pos, pos) <- new_base
    return(rep_seq)
  }

  # ---- Função para gerar sequências contendo telômeros ----
  generate_telomeric_seq <- function() {
    # número de repetições consecutivas
    repeats <- sample(min_repeats:max_repeats, 1)

    # constrói bloco
    rep_unit <- telomere_repeat
    block <- character(repeats)

    for (i in seq_len(repeats)) {
      if (runif(1) < variant_prob) {
        block[i] <- create_variant(rep_unit)
      } else {
        block[i] <- rep_unit
      }
    }

    block <- paste0(block, collapse = "")

    # adiciona flancos aleatórios
    flank_left <- paste0(sample(c("A","C","G","T"), sample(20:80,1), replace = TRUE), collapse = "")
    flank_right <- paste0(sample(c("A","C","G","T"), sample(20:80,1), replace = TRUE), collapse = "")

    paste0(flank_left, block, flank_right)
  }

  # ---- Função para gerar sequências totalmente aleatórias ----
  generate_random_seq <- function() {
    len <- sample(80:200, 1)
    paste0(sample(c("A","C","G","T"), len, replace = TRUE), collapse = "")
  }

  # ---- Gerar todas as sequências ----
  seqs_telomeric <- replicate(n_telomeric, generate_telomeric_seq())
  seqs_random <- replicate(n_random, generate_random_seq())

  all_seqs <- c(seqs_telomeric, seqs_random)
  names(all_seqs) <- paste0("seq_", seq_len(n_seqs))

  # ---- Salvar FASTA ----
  dna <- Biostrings::DNAStringSet(all_seqs)
  names(dna) <- names(all_seqs)

  Biostrings::writeXStringSet(dna, output)

  message("Arquivo gerado: ", output)
  message(n_telomeric, " sequências teloméricas")
  message(n_random, " sequências aleatórias")

  return(dna)
}
