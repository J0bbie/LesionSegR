#' @title Determine mutational burden based on GRCm39.
#' @param x (tibble): Tibble containing the mutations
#' @return (tibble): Tibble with information about the mutational burden.
#' @importFrom dplyr %>%
#' @export
determine_mutational_burden <- function(x) {
    # Sanity check. ----
    checkmate::checkTibble(x)
    
    # Calc. TMB stats. ----
    tmb <- tibble::tibble(
        # Number of mappable ATCG in reference genome (GRCm39).
        TMB = nrow(x) / (2649938115 / 1E6),
        totalMutations = nrow(x),
        totalMutationsCoding = x %>% dplyr::filter(HGVSp != "") %>% dplyr::n_distinct(),
        totalH1 = x %>% dplyr::filter(origin_mutant == "H1") %>% dplyr::n_distinct(),
        totalH2 = x %>% dplyr::filter(origin_mutant == "H2") %>% dplyr::n_distinct(),
        totalUA = x %>% dplyr::filter(origin_mutant == "UA") %>% dplyr::n_distinct(),
        sample = unique(x$sample)
    )
    
    # Return statement. ----
    return(tmb)
}
