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
#' @param fasta_file Path to the FASTA file containing the sequences.
#' @param canonical_repeat Canonical repeat sequence (e.g., `"TTAGGG"`).
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
detect_TVS <- function(fasta_file, canonical_repeat, min_var_len = 2, max_var_len = 8) {

  # Reading and preparing sequences
  fasta_sequences <- Biostrings::readDNAStringSet(filepath = fasta_file)
  all_sequences <- c(fasta_sequences, Biostrings::reverseComplement(fasta_sequences))
  names(all_sequences) <- c(names(fasta_sequences), paste0(names(fasta_sequences), "_RC"))

  len_canon <- nchar(canonical_repeat)
  all_repeats <- list()

  # Iterate over all sequences (original and reverse-complement)
  for (i in seq_along(all_sequences)) {
    seq_id <- names(all_sequences)[i]
    seq_dna <- all_sequences[[i]]
    seq_string <- as.character(seq_dna)

    # 1. Find ALL canonical occurrences
    canonical_hits <- Biostrings::matchPattern(pattern = canonical_repeat, subject = seq_dna) %>%
      as.data.frame() %>%
      dplyr::filter(width == len_canon) %>%
      # FIX: explicitly add seqnames column
      dplyr::mutate(seqnames = seq_id, type = "Canonical", seq_match = canonical_repeat)

    if (nrow(canonical_hits) < 2) {
      next
    }

    canonical_hits <- canonical_hits[order(canonical_hits$start), ]

    # 2. Find VARIANTS between consecutive canonical repeats
    for (k in 1:(nrow(canonical_hits) - 1)) {
      next_rep_start <- canonical_hits$start[k + 1]
      current_rep_end <- canonical_hits$end[k]

      variant_start <- current_rep_end + 1
      variant_end <- next_rep_start - 1

      if (variant_start <= variant_end) {
        variant_seq <- substr(seq_string, start = variant_start, stop = variant_end)
        variant_len <- nchar(variant_seq)

        # Variant criterion: correct length AND sequence is NOT the canonical repeat
        if (variant_len >= min_var_len && variant_len <= max_var_len && variant_seq != canonical_repeat) {

          # Add variant found
          all_repeats[[length(all_repeats) + 1]] <- data.frame(
            seqnames = seq_id,
            start = variant_start,
            end = variant_end,
            width = variant_len,
            type = "Variant",
            seq_match = variant_seq
          )
        }
      }
    }

    # 3. Add canonical repeats for plot proportion calculations
    all_repeats[[length(all_repeats) + 1]] <- canonical_hits %>%
      dplyr::select(seqnames, start, end, width, type, seq_match)
  }

  if (length(all_repeats) > 0) {
    # Combine canonical repeats and variants
    final_df <- do.call(rbind, all_repeats) %>%
      # Create the logical column 'variant' (used in plot_TVS)
      dplyr::mutate(variant = type == "Variant")

    return(final_df)
  } else {
    message("No telomeric repeats or variants were found.")
    return(NULL)
  }
}
