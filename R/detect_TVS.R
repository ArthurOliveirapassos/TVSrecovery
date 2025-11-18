#' Detectar TVSs (Telomeric Variant Sequences) em sequências FASTA
#'
#' A função `detect_TVS()` identifica repetições teloméricas canônicas e
#' variações entre repetições consecutivas em arquivos FASTA.
#' Ela processa automaticamente a fita direta e a fita reverso-complementar,
#' garantindo que todas as regiões teloméricas potenciais sejam avaliadas.
#'
#' O algoritmo funciona em três etapas principais:
#' 1. Identifica todas as ocorrências **canônicas** da repetição informada.
#' 2. Para cada par de repetições consecutivas, extrai a sequência intermediária.
#' 3. Se essa sequência tiver entre `min_var_len` e `max_var_len` bases
#'    e não for idêntica à repetição canônica, é registrada como TVS (variante).
#'
#' O output combina tanto repetições canônicas quanto variantes, permitindo
#' posterior visualização através de `plot_TVS()`.
#'
#' @param arquivo_fasta Caminho para o arquivo FASTA contendo as sequências.
#' @param repeticao_canon Sequência repetitiva canônica (ex: `"TTAGGG"`).
#' @param min_var_len Tamanho mínimo permitido para uma TVS (default: 2).
#' @param max_var_len Tamanho máximo permitido para uma TVS (default: 8).
#'
#' @return Um `data.frame` contendo repetições canônicas e variantes, com colunas:
#' \describe{
#'   \item{seqnames}{Nome da sequência onde a repetição/TBS foi encontrada}
#'   \item{start}{Posição inicial}
#'   \item{end}{Posição final}
#'   \item{width}{Tamanho da repetição/variante}
#'   \item{type}{Tipo: `"Canonica"` ou `"Variante"`}
#'   \item{seq_match}{Sequência encontrada}
#'   \item{variant}{TRUE para variantes, FALSE para canônicas}
#' }
#' Caso nenhuma repetição ou TVS seja encontrada, retorna `NULL`.
#'
#' @examples
#' \dontrun{
#' resultado <- detect_TVS(
#'   arquivo_fasta = "telomere_hits.fasta",
#'   repeticao_canon = "TTAGGG",
#'   min_var_len = 2,
#'   max_var_len = 8
#' )
#' }
#'
#'
#' @import Biostrings
#' @import dplyr
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
