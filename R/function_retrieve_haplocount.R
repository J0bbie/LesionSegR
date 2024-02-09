#' @title Count haplotype-specific reads for a given GRanges.
#' @param granges (GRanges): GRanges object.
#' @param path_bam (character): Path to BAM file.
#' @param threads (integer): No. of threads to parallize.
#' @return (tibble): Tibble object with haplotype-specific read counts.
#' @importFrom dplyr %>%
#' @import plyranges
#' @export
retrieve_haplocount <- function(granges, path_bam, threads = 10) {
    
    # Input validation ----
    
    checkmate::assertClass(granges, "GRanges")
    checkmate::assertFile(path_bam)
    checkmate::assertNumber(threads)
    
    futile.logger::flog.info(glue::glue("Counting haplotype-specific reads for {length(granges)} region(s) for {basename(path_bam)}."))
    
    # Check if granges has 'id' column.
    if (is.null(S4Vectors::mcols(granges)$id)) {
        stop("GRanges object must have 'id' column.")
    }
    
    # Generate the which_label column.
    granges$which_label <- sprintf("%s:%s-%s", GenomicRanges::seqnames(granges), GenomicRanges::start(granges), GenomicRanges::end(granges))
    
    # Retrieve reads from BAM. ----
    
    # Only select reads with vW flag is 1 and primary reads only.
    flag <- Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE, isDuplicate = FALSE)
    
    
    # Retrieve all alignment (per region).
    count_bam <- dplyr::bind_rows(pbapply::pblapply(base::seq_along(granges), function(i){
        
        filter <- Rsamtools::ScanBamParam(flag = flag, tag = c('HP', 'vW'), which = granges[i])
        data_bam <- GenomicAlignments::readGAlignments(path_bam, param = filter, with.which_label = TRUE)
        
        # Only keep reads which cover at least 10bp overlap.
        data_bam <- IRanges::subsetByOverlaps(data_bam, granges, minoverlap = 10)
        
        # Summarize counts.
        count_bam_region <- GenomicRanges::GRanges(data_bam) %>% 
            plyranges::group_by(which_label) %>% 
            plyranges::summarize(
                # Total reads with WASP-Good or no WASP tag.
                totalCounts = sum(vW == 1 | is.na(vW)),
                # Total H1 + WASP-Good
                H1 = sum(HP[vW == 1 & !is.na(vW)] == 1, na.rm = T),
                # Total H2 + WASP-Good
                H2 = sum(HP[vW == 1 & !is.na(vW)] == 2, na.rm = T),
                # Total AU + WASP-Good
                UA = sum(is.na(HP)[vW == 1 | is.na(vW)], na.rm = T),
            ) %>% 
            tibble::as_tibble()
        
        return(count_bam_region)
        
    }, cl = threads))
    
    # Add sample name.
    count_bam$sample <- base::factor(gsub("_.*", "", base::basename(path_bam)))
    
    # Add the id column.
    count_bam$id <- granges[BiocGenerics::match(count_bam$which_label, granges$which_label)]$id
    count_bam$seqnames <- as.factor(GenomicRanges::seqnames(granges[BiocGenerics::match(count_bam$which_label, granges$which_label)]))
    count_bam$start <- GenomicRanges::start(granges[BiocGenerics::match(count_bam$which_label, granges$which_label)])
    count_bam$end <- GenomicRanges::end(granges[BiocGenerics::match(count_bam$which_label, granges$which_label)])
    
    # Convert to tibble.
    count_bam <- count_bam %>% dplyr::distinct(
        seqnames, start, end, which_label, totalCounts, H1, H2, UA, sample, id
    )
    
    # Return statement ----
    
    return(count_bam)
}
