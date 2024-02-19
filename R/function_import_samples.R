#' @title Import WGS and WTS data from the snakemake-lesionsegration workflow.
#' @param metadata (tibble): Metadata of the samples.
#' @param workflow_dir (character): Path to the workflow output directory.
#' @param gtf (tibble): Tibble to imported/filtered GTF.
#' @return (list): List of 96-context matrices (SBS, InDel, DBS).
#' @importFrom dplyr %>%
#' @export
import_samples <- function(metadata, workflow_dir, gtf = NULL){
    
    # Sanity check. ----
    checkmate::assertTibble(metadata)
    checkmate::assertAccess(workflow_dir, access = "r")
    checkmate::assertTibble(gtf, null.ok = TRUE)
    
    futile.logger::flog.info(glue::glue("Retrieving data from {dplyr::n_distinct(metadata$sample_name)} samples (WGS + WTS)"))
    
    # Loop over the samples. ----
    p <- progressr::progressor(along = metadata$sample_name)
    data_samples <- future.apply::future_lapply(unique(metadata$sample_name), function(x){
        
        # Determine the paths. ----
        
        # WGS-specific files.
        current_sample_WGS <- metadata %>% 
            dplyr::filter(sample_name == x, sequencing_type == 'WGS') %>% 
            dplyr::mutate(
                path_vcf = file.path(workflow_dir, glue::glue("{experiment_name}/variants/{sequencing_name}_{strain1}_{strain2}_VEP_haplocounted.vcf.gz")),
                path_flagstats = file.path(workflow_dir, glue::glue("{experiment_name}/alignment/{sequencing_type}/{sequencing_name}_{strain1}_{strain2}_sortedByCoord_markDup_haplotagged.bam.flagstat")),
                path_haplotag = file.path(workflow_dir, glue::glue("logs/{experiment_name}/haplotyping/haplotag_{sequencing_type}_{sequencing_name}_{strain1}_{strain2}.log"))
            )
        
        # WTS-specific files.
        current_sample_WTS <- metadata %>% 
            dplyr::filter(sample_name == x, sequencing_type == 'WTS') %>% 
            dplyr::mutate(
                path_flagstats = file.path(workflow_dir, glue::glue("{experiment_name}/alignment/{sequencing_type}/{sequencing_name}_{strain1}_{strain2}_sortedByCoord_markDup_haplotagged.bam.flagstat")),
                path_haplotag = file.path(workflow_dir, glue::glue("logs/{experiment_name}/haplotyping/haplotag_{sequencing_type}_{sequencing_name}_{strain1}_{strain2}.log")),
                path_featurecounts = file.path(workflow_dir, glue::glue("{experiment_name}/counting/{sequencing_type}/{sequencing_name}_{strain1}_{strain2}_counts.txt"))
            )
        
        if(nrow(current_sample_WGS) != 0 & nrow(current_sample_WTS) == 0) futile.logger::flog.info(glue::glue("{x}: WGS-only"))
        if(nrow(current_sample_WGS) == 0 & nrow(current_sample_WTS) != 0) futile.logger::flog.info(glue::glue("{x}: WTS-only"))
        if(nrow(current_sample_WGS) != 0 & nrow(current_sample_WTS) != 0) futile.logger::flog.info(glue::glue("{x}: WGS + WTS."))
        
        # Import data. ----
        current_sample <- list()
        
        # Import flagstats. ----
        current_sample$flagstats <- LesionSegR::read_flagstats(c(current_sample_WGS$path_flagstats, current_sample_WTS$path_flagstats))
        
        # Import haplotyping stats. ----
        current_sample$haplotyping <- LesionSegR::read_haplotag_log(c(current_sample_WGS$path_haplotag, current_sample_WTS$path_haplotag))
        
        # Process somatic data.  ----
        # Normals have no somatic variants.
        if(!is.na(current_sample_WGS$matched_group)){
            # Import somatic variants. ----
            current_sample$somaticvariants <- LesionSegR::import_vcf(current_sample_WGS$path_vcf)
            
            if(!is.null(current_sample$somaticvariants)){
                # Calc. TMB.
                current_sample$tumorburden <- LesionSegR::determine_mutational_burden(current_sample$somaticvariants)
                
                # Generate the 96-context matrices. ----
                current_sample$mutmatrices <- LesionSegR::generate_mutmatrices_96(current_sample$somaticvariants)
            }
        }
        
        # Import gene-wise counts (Total / H1 / H2 / UA) ----
        if(nrow(current_sample_WTS) != 0){
            current_sample$featurecounts <- data.table::fread(current_sample_WTS$path_featurecounts, col.names = c('gene_id', 'chr', 'start', 'end', 'strand', 'length', 'totalCounts', 'H1', 'H2', 'UA')) %>% 
                dplyr::distinct(gene_id, totalCounts, H1, H2, UA) %>% 
                dplyr::mutate(sample = factor(current_sample_WTS$sample_name))
            
            # Subset on pre-determined genes.
            if(!is.null(gtf)){
                current_sample$featurecounts <- current_sample$featurecounts %>% dplyr::inner_join(gtf, by = 'gene_id')
            }
            
            # Pivot longer.
            current_sample$featurecounts <- current_sample$featurecounts %>% 
                tidyr::pivot_longer(., cols = c('totalCounts', 'H1', 'H2', 'UA')) %>% 
                dplyr::mutate(
                    gene_id = factor(gene_id),
                    seqnames = factor(seqnames),
                    gene_type = factor(gene_type),
                    gene_name = factor(gene_name),
                    name = factor(name)
                )
        }
        
        # Return.
        p(sprintf("%s", x))
        return(current_sample)
        
    })
    
    # Combine samples.
    data_combined <- list()
    data_combined$flagstats <- dplyr::bind_rows(lapply(data_samples, function(x) x$flagstats))
    data_combined$haplotyping <- dplyr::bind_rows(lapply(data_samples, function(x) x$haplotyping))
    data_combined$somaticvariants <- dplyr::bind_rows(lapply(data_samples, function(x) x$somaticvariants))
    data_combined$tumorburden <- dplyr::bind_rows(lapply(data_samples, function(x) x$tumorburden))
    data_combined$mutmatrices_sbs <- do.call(cbind, lapply(data_samples, function(x) x$mutmatrices$sbs))
    data_combined$mutmatrices_indel <- do.call(cbind, lapply(data_samples, function(x) x$mutmatrices$indel))
    
    # Combine Featurecounts.
    data_combined$featurecounts <- tryCatch({ plyr::join_all(plyr::compact(lapply(data_samples, function(x) x$featurecounts)), by='context', type='left') }, error = function(e) { NULL })
    
    # Return statement. ----
    return(data_combined)
}