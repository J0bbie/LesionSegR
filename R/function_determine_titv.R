#' @title Determine Ti/Tv ratio.
#' @param m_sbs (matrix): Mutational SBS matrix.
#' @return (tibble): Tibble with information about the Ti/Tv.
#' @importFrom dplyr %>%
#' @export
determine_ti_tv <- function(m_sbs) {

    # Input validation ----

    checkmate::assertMatrix(m_sbs)
    futile.logger::flog.info(glue::glue("Calculating Ti/Tv ratios for {dplyr::n_distinct(base::colnames(m_sbs))} sample(s)."))

    # Determine Ti/Tv mutations and ratio ----

    data_mutcontext <- m_sbs %>%
        reshape2::melt() %>%
        dplyr::select(mut_context = Var1, sample = Var2, value) %>%
        dplyr::mutate(
            mutational_context = base::gsub("\\[.*]", "", mut_context),
            mutational_type = base::gsub("].*", "", base::gsub(".*\\[", "", mut_context)),
            # Determine CpG mutations.
            mutational_type = base::ifelse(base::grepl("C>T", mutational_type) & base::grepl("G$", mutational_context), base::paste(mutational_type, "(CpG)"), mutational_type),
            # Determine Ti / Tv mutations.
            mutational_type = dplyr::if_else(mutational_type %in% c("C>T", "C>T (CpG)", "T>C", "G>A", "A>T"),
                                             base::paste(mutational_type, "Ti", sep = "\n"), base::paste(mutational_type, "Tv", sep = "\n")
            )
        ) %>%
        # Count total mut. contexts per sample.
        dplyr::group_by(sample, mutational_type) %>%
        dplyr::summarise(total_value = base::sum(value)) %>%
        dplyr::ungroup() %>%
        # Determine Ti/Tv ratio per sample.
        dplyr::mutate(type_ti_tv = ifelse(grepl("\nTv", mutational_type), "Transversion", "Transition")) %>%
        dplyr::group_by(sample, type_ti_tv) %>%
        dplyr::mutate(totalgroup_sample = sum(total_value)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(TiTvRatio = unique(totalgroup_sample[type_ti_tv == "Transition"]) / unique(totalgroup_sample[type_ti_tv == "Transversion"])) %>%
        dplyr::ungroup()


    # Return statement --------------------------------------------------------

    return(data_mutcontext)
}
