#' Plotar variantes teloméricas e proporção canônica vs variante
#'
#' Esta função recebe o data frame gerado por `detect_TVS()` e produz dois
#' gráficos úteis para visualizar as repetições teloméricas:
#'
#' **1. Gráfico de barras** mostrando as variantes mais frequentes
#' **2. Gráfico de pizza** mostrando a proporção entre repetições canônicas e variantes
#'
#' @param variant_df Um data frame gerado por `detect_TVS()`, contendo
#'   repetições canônicas e variantes. É necessário que o objeto possua as
#'   colunas: `seq_match`, `type` e `variant`.
#'
#' @return Uma lista contendo dois objetos `ggplot2`:
#' \describe{
#'   \item{bar_plot}{Gráfico de barras com variantes recorrentes (`n > 1`).}
#'   \item{pie_plot}{Gráfico de pizza mostrando a proporção entre repetições do tipo canônica e variantes.}
#' }
#'
#' Caso o data frame esteja vazio, a função retorna `NULL` e emite um aviso.
#'
#' @examples
#' \dontrun{
#' resultado <- detect_TVS("exemplo.fasta", repeticao_canon = "TTAGGG")
#' graficos <- plot_TVS(resultado)
#'
#' # Exibir gráficos
#' graficos$bar_plot
#' graficos$pie_plot
#' }
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stats reorder
#' @export
plot_TVS <- function(variant_df) {

  if (is.null(variant_df) || nrow(variant_df) == 0) {
    warning("O data frame está vazio. Não é possível gerar gráficos.")
    return(NULL)
  }

  # -------------------------------
  # GRÁFICO DE BARRAS DAS VARIANTES
  # -------------------------------
  # Filtra apenas as variantes (variant == TRUE)
  bar_plot <- variant_df %>%
    dplyr::filter(variant == TRUE) %>%
    dplyr::count(seq_match, sort = TRUE) %>%
    dplyr::filter(n > 1) %>%  # mostra apenas variantes recorrentes
    ggplot2::ggplot(ggplot2::aes(x = stats::reorder(seq_match, n), y = n)) +
    ggplot2::geom_col(fill = "#0072B2") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Variantes teloméricas mais frequentes (> 1 ocorrência)",
      x = "Variante",
      y = "Frequência"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  # -------------------------------
  # GRÁFICO DE PIZZA (CANÔNICA vs VARIANTE)
  # -------------------------------
  pie_plot <- variant_df %>%
    # A coluna 'type' já existe na saída da função adaptada
    dplyr::count(type) %>%
    # Garante que as porcentagens sejam calculadas
    dplyr::mutate(percent = n / sum(n) * 100,
                  label = paste0(type, "\n", round(percent, 1), "%")) %>%
    ggplot2::ggplot(ggplot2::aes(x = "", y = n, fill = type)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::coord_polar(theta = "y") +
    # Adiciona rótulos percentuais
    ggplot2::geom_text(ggplot2::aes(y = n / 2 + c(0, cumsum(n)[-length(n)]),
                                    label = label),
                       color = "black", size = 4) +
    ggplot2::labs(
      title = "Proporção entre repetições canônicas e variantes",
      fill = "Tipo de Repetição"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  # retorno das duas figuras
  return(list(
    bar_plot = bar_plot,
    pie_plot = pie_plot
  ))
}
