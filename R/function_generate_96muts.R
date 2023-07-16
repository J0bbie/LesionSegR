#' @title Generate 96-context matrix of somatic variants.
#' @param x (SimpleVRangesList): VRangesList containing the mutations of all samples.
#' @return (list): List of 96-context matrices (SBS, InDel, DBS).
#' @importFrom dplyr %>%
#' @export
generate_mutmatrices_96 <- function(x) {

    # Check input. ----

    checkmate::checkClass(x, "SimpleVRangesList")


    # Sub function - Convert mutations to correct GRanges for input into MutationalPatterns. ----
    convert_muts <- function(x, dbs = FALSE) {
        S4Vectors::mcols(x) <- S4Vectors::DataFrame(sample = x$sample)

        # Add REF and ALT as column.
        if (dbs) {
            x$REF <- Biostrings::DNAStringSet(VariantAnnotation::ref(x))
            x$ALT <- Biostrings::DNAStringSetList(base::lapply(VariantAnnotation::alt(x), Biostrings::DNAStringSet))
        } else {
            x$REF <- VariantAnnotation::ref(x)
            x$ALT <- VariantAnnotation::alt(x)
        }

        # Convert to GRangeslist, split per sample.
        x <- GenomicRanges::GRanges(x)
        x <- GenomicRanges::GRangesList(base::split(x, x$sample))

        GenomeInfoDb::genome(x) <- "mm39"

        # Return.
        return(x)
    }

    data_mutmatrices <- list()
    x <- base::unlist(x)

    ## SBS. ---
    data_mutmatrices$sbs <- x[x$mutType == "SNV", ] %>%
        convert_muts() %>%
        MutationalPatterns::mut_matrix(ref_genome = "BSgenome.Mmusculus.UCSC.mm39")

    ## InDel. ---
    data_mutmatrices$indel <- x[x$mutType == "InDel", ] %>%
        convert_muts() %>%
        MutationalPatterns::get_indel_context(., ref_genome = "BSgenome.Mmusculus.UCSC.mm39") %>%
        MutationalPatterns::count_indel_contexts(.)

    ## DBS. ----
    data_mutmatrices$dbs <- x[x$mutType == "MNV" & base::nchar(VariantAnnotation::ref(x)) == 2 & base::nchar(VariantAnnotation::alt(x)) == 2] %>%
        convert_muts(., dbs = TRUE) %>%
        MutationalPatterns::get_dbs_context(.) %>%
        MutationalPatterns::count_dbs_contexts(.)

    # Return.
    return(data_mutmatrices)
}
