#' Detect TVSs (Telomeric Variant Sequences) in FASTA sequences
#'
#' The `detect_TVS()` function identifies canonical telomeric repeats and
#' sequence variations occurring between consecutive repeats in FASTA files.
#' It automatically processes both the forward strand and the reverse-complement
#' strand, ensuring that all potential telomeric regions are evaluated.
#'
#' The algorithm operates in three main steps:
#' 1. Identifies all **canonical** occurrences of the provided telomeric repeat.
#' 2. For each pair of consecutive repeats, extracts the intermediate sequence.
#' 3. If this sequence has between `min_var_len` and `max_var_len` bases
#'    and is not identical to the canonical repeat, it is recorded as a TVS (variant).
#'
#' The output combines both canonical repeats and variants, enabling further
#' visualization through `plot_TVS()`.
#'
#' @param arquivo_fasta Path to the FASTA file containing the sequences.
#' @param repeticao_canon Canonical repeat sequence (e.g., `"TTAGGG"`).
#' @param min_var_len Minimum allowed size for a TVS (default: 2).
#' @param max_var_len Maximum allowed size for a TVS (default: 8).
#'
#' @return A `data.frame` containing canonical repeats and variants, with columns:
#' \describe{
#'   \item{seqnames}{Name of the sequence where the repeat/TVS was found}
#'   \item{start}{Start position}
#'   \item{end}{End position}
#'   \item{width}{Length of the repeat/variant}
#'   \item{type}{Type: `"Canonical"` or `"Variant"`}
#'   \item{seq_match}{Matched sequence}
#'   \item{variant}{TRUE for variants, FALSE for canonical repeats}
#' }
#' If no repeats or TVSs are found, returns `NULL`.
#'
#' @examples
#' \dontrun{
#' result <- detect_TVS(
#'   arquivo_fasta = "telomere_hits.fasta",
#'   repeticao_canon = "TTAGGG",
#'   min_var_len = 2,
#'   max_var_len = 8
#' )
#' }
#'
#' @importFrom Biostrings readDNAStringSet reverseComplement matchPattern
#' @importFrom dplyr filter mutate select count %>%
#' @importFrom stats runif
#'
#' @export
detect_TVS <- function(arquivo_fasta, repeticao_canon, min_var_len = 2, max_var_len = 8) {

  # Leitura e preparo das sequências
  sequencias_fasta <- Biostrings::readDNAStringSet(filepath = arquivo_fasta)
  todas_sequencias <- c(sequencias_fasta, Biostrings::reverseComplement(sequencias_fasta))
  names(todas_sequencias) <- c(names(sequencias_fasta), paste0(names(sequencias_fasta), "_RC"))

  len_canon <- nchar(repeticao_canon)
  todas_repeticoes <- list()

  # Itera sobre todas as sequências (originais e complementares reversas)
  for (i in seq_along(todas_sequencias)) {
    seq_id <- names(todas_sequencias)[i]
    seq_dna <- todas_sequencias[[i]]
    seq_string <- as.character(seq_dna)

    # 1. Encontrar TODAS as ocorrências canônicas
    ocorrencias_canonicas <- Biostrings::matchPattern(pattern = repeticao_canon, subject = seq_dna) %>%
      as.data.frame() %>%
      dplyr::filter(width == len_canon) %>%
      # CORREÇÃO: Adicionar explicitamente a coluna seqnames
      dplyr::mutate(seqnames = seq_id, type = "Canonica", seq_match = repeticao_canon)

    if (nrow(ocorrencias_canonicas) < 2) {
      next
    }

    ocorrencias_canonicas <- ocorrencias_canonicas[order(ocorrencias_canonicas$start), ]

    # 2. Encontrar VARIAÇÕES entre as repetições canônicas consecutivas
    for (k in 1:(nrow(ocorrencias_canonicas) - 1)) {
      start_prox_rep <- ocorrencias_canonicas$start[k+1]
      end_rep_atual <- ocorrencias_canonicas$end[k]

      start_variante <- end_rep_atual + 1
      end_variante <- start_prox_rep - 1

      if (start_variante <= end_variante) {
        variante_seq <- substr(seq_string, start = start_variante, stop = end_variante)
        variante_len <- nchar(variante_seq)

        # O critério de variante é: tamanho correto E a sequência NÃO é a canônica
        if (variante_len >= min_var_len && variante_len <= max_var_len && variante_seq != repeticao_canon) {

          # Adicionar a variante encontrada
          todas_repeticoes[[length(todas_repeticoes) + 1]] <- data.frame(
            seqnames = seq_id,
            start = start_variante,
            end = end_variante,
            width = variante_len,
            type = "Variante",
            seq_match = variante_seq
          )
        }
      }
    }

    # 3. Adicionar as repetições canônicas originais para o cálculo da proporção
    # Agora a seleção funciona, pois a coluna seqnames foi adicionada no Passo 1
    todas_repeticoes[[length(todas_repeticoes) + 1]] <- ocorrencias_canonicas %>%
      dplyr::select(seqnames, start, end, width, type, seq_match)
  }

  if (length(todas_repeticoes) > 0) {
    # Combina todas as repetições (canônicas e variantes)
    final_df <- do.call(rbind, todas_repeticoes) %>%
      # Cria a coluna booleana 'variant' para uso no plot_TVS
      dplyr::mutate(variant = type == "Variante")

    return(final_df)
  } else {
    message("Nenhuma repetição ou variante telomérica encontrada.")
    return(NULL)
  }
}
