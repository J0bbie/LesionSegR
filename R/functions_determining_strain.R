#' @title  Determine origin based on mutual exclusivity of H1 / H2 mutant reads.
#' @description Annotate strain of variants based on mutual presence of H1 or H2 reads.
#'
#' @param count_strain1 (numeric): H1 counts of variants.
#' @param count_strain2 (numeric): H2 counts of variants.
#'
#' @examples
#' 1 + 1
#'
#' @return (character): Vector of H1, H2, or UA status.
#'
#' @importFrom dplyr %>%
determine_origin_mutual <- function(count_strain1, count_strain2) {
    tibble::tibble(count_strain1, count_strain2) %>%
        dplyr::mutate(
            assignment = dplyr::case_when(
                .default = "UA",
                count_strain1 <= 1 & count_strain2 <= 1 ~ "UA",
                (count_strain1 > 0 & count_strain2 > 0) ~ "UA",
                count_strain1 > count_strain2 ~ "H1",
                count_strain2 > count_strain1 ~ "H2"
            )
        ) %>%
        dplyr::pull(assignment)
}
