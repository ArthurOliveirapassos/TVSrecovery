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
    stop("Argument 'fasta_file' is missing. Please provide a FASTA file.")
  }

  if (is.null(telomere_repeat)) {
    telomere_repeat <- readline(prompt = "Enter the repetitive sequence: ")
  }

  if (is.null(min_repeats)) {
    min_repeats <- as.numeric(readline(prompt = "Enter the minimum number of consecutive repeats: "))
  }

  # Read FASTA
  fasta <- Biostrings::readDNAStringSet(fasta_file)

  # Create exact consecutive block
  repeat_block <- paste0(rep(telomere_repeat, min_repeats), collapse = "")

  # Generate equivalent patterns
  dna_pat <- Biostrings::DNAString(repeat_block)
  reverse_pat <- Biostrings::reverse(dna_pat)
  complement_pat <- Biostrings::complement(dna_pat)
  revcomp_pat <- Biostrings::reverseComplement(dna_pat)

  # Count exact occurrences of each pattern
  hits_direct <- Biostrings::vcountPattern(dna_pat, fasta, fixed = TRUE)
  hits_rev <- Biostrings::vcountPattern(reverse_pat, fasta, fixed = TRUE)
  hits_comp <- Biostrings::vcountPattern(complement_pat, fasta, fixed = TRUE)
  hits_revcomp <- Biostrings::vcountPattern(revcomp_pat, fasta, fixed = TRUE)

  # Keep only sequences that contain AT LEAST ONE exact block
  keep <- (hits_direct > 0 |
             hits_rev > 0 |
             hits_comp > 0 |
             hits_revcomp > 0)

  matched <- fasta[keep]

  # Save results
  if (length(matched) > 0) {
    Biostrings::writeXStringSet(matched, filepath = output_fasta)

    message(
      length(matched),
      " sequences containing at least one block with ",
      min_repeats, " consecutive repeats of '", telomere_repeat,
      "' (or reverse/complement) were saved to: ",
      output_fasta
    )
  } else {
    message(
      "No sequence containing a block of ",
      min_repeats, " consecutive repeats of '", telomere_repeat,
      "' was found."
    )
  }

  return(matched)
}



