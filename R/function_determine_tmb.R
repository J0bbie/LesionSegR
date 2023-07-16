#' @title Determine mutational burden based on GRCm39.
#' @param x (SimpleVRangesList): VRangesList containing the mutations
#' @return (tibble): Tibble with information about the mutational burden.
#' @importFrom dplyr %>%
#' @export
determine_mutational_burden <- function(x) {

    checkmate::checkClass(x, "SimpleVRangesList")

    dplyr::bind_rows(future.apply::future_lapply(x, function(y){

        tibble::tibble(
            # Number of mappable ATCG in reference genome (GRCm39).
            TMB = length(y) / (2649938115 / 1E6),
            totalMutations = length(y),
            totalMutationsCoding = length(y[y$HGVSp != "", ]),
            totalG1 = sum(y$origin_mutant == "G1"),
            totalG2 = sum(y$origin_mutant == "G2"),
            totalUA = sum(y$origin_mutant == "UA"),
            sample = unique(y$sample)
        )
    }))

}
