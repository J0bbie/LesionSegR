#' @title  Determine origin based on mutual exclusivity of G1 / G2 mutant reads.
#' @description Annotate strain of variants based on mutual presence of G1 or G2 reads.
#'
#' @param count_strain1 (numeric): G1 counts of variants.
#' @param count_strain2 (numeric): G2 counts of variants.
#'
#' @examples
#' 1 + 1
#'
#' @return (character): Vector of G1, G2, or UA status.
#'
#' @importFrom dplyr %>%
determine_origin_mutual <- function(count_strain1, count_strain2) {
    tibble::tibble(count_strain1, count_strain2) %>%
        dplyr::mutate(
            assignment = dplyr::case_when(
                .default = "UA",
                count_strain1 <= 1 & count_strain2 <= 1 ~ "UA",
                (count_strain1 > 0 & count_strain2 > 0) ~ "UA",
                count_strain1 > count_strain2 ~ "G1",
                count_strain2 > count_strain1 ~ "G2"
            )
        ) %>%
        dplyr::pull(assignment)
}
