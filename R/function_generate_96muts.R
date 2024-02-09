#' @title Generate 96-context matrix of somatic variants.
#' @param x (VRanges): VRanges containing the mutations of the samples.
#' @return (list): List of 96-context tibbles (SBS and InDel).
#' @importFrom dplyr %>%
#' @export
generate_mutmatrices_96 <- function(x) {
    
    # Check input. ----
    
    checkmate::checkTibble(x)
    
    # Sub function - Convert mutations to correct GRanges for input into MutationalPatterns. ----
    convert_muts <- function(x) {
        xx <- GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = T)
        xx$REF <- xx$ref
        xx$ALT <- xx$alt
        
        GenomeInfoDb::genome(xx) <- "mm39"
        
        # Return.
        return(xx)
    }
    
    data_mutmatrices <- list()
    
    ## SBS. ---
    data_mutmatrices$sbs <- x %>% dplyr::filter(mutType == 'SNV') %>%
        convert_muts() %>%
        MutationalPatterns::mut_matrix(ref_genome = "BSgenome.Mmusculus.UCSC.mm39")
    
    colnames(data_mutmatrices$sbs) <- unique(x$sample)
    data_mutmatrices$sbs <- tibble::as_tibble(data_mutmatrices$sbs, rownames = 'context')
    
    ## InDel. ---
    if(nrow(x[x$mutType == "InDel", ]) != 0){
        data_mutmatrices$indel <- x[x$mutType == "InDel", ] %>%
            convert_muts() %>%
            MutationalPatterns::get_indel_context(., ref_genome = "BSgenome.Mmusculus.UCSC.mm39") %>%
            MutationalPatterns::count_indel_contexts(.)
        
        colnames(data_mutmatrices$indel) <- unique(x$sample)
        data_mutmatrices$indel <- tibble::as_tibble(data_mutmatrices$indel, rownames = 'context')
    }
    
    # Return.
    return(data_mutmatrices)
}
