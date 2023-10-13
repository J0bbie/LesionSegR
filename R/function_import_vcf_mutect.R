#' @title Importing VCF files.
#' @description Imports VCF files from the snakemake-lesionsegregation workflow.
#' Applies several filtering steps and determines the likely strain-origin of the somatic variants.
#' @param x (character): Path to one or multiple VCF file(s).
#' @examples
#' 1 + 1
#' @return (VRanges): VRanges object containing the (filtered) somatic variants
#' and strain-specific tagging.
#' @importFrom dplyr %>%
#' @export
import_vcf_mutect <- function(x) {
    # Input validation ---
    
    checkmate::assertFile(x)
    
    # Import and clean-up VCF. ---
    data_somaticvariants <- future.apply::future_lapply(x, function(path_vcf) {
        futile.logger::flog.info(glue::glue("Importing VCF: {path_vcf}"))
        
        # Read in the VCF file.
        sample_vcf <- VariantAnnotation::readVcf(path_vcf, genome = "GRCm39")
        
        # Expand multi-allelic hits to separate rows.
        sample_vcf <- VariantAnnotation::expand(sample_vcf)
        
        # Filter VCF for PASS-only.
        filter_pass <- base::length(sample_vcf)
        sample_vcf <- sample_vcf[sample_vcf@fixed$FILTER == "PASS", ]
        
        futile.logger::flog.info(glue::glue("Filtering PASS-only: retaining  of {filter_pass} somatic variants."))
        
        # Filter on DP.
        filter_dp <- base::length(sample_vcf)
        sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$DP[, 2] >= 10, ]
        
        futile.logger::flog.info(glue::glue("Filtering on DP >= 10: retaining {base::length(sample_vcf)} of {filter_dp} somatic variants."))
        
        # Filter on AD.
        filter_ad <- base::length(sample_vcf)
        sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$AD[, , 2][,2] >= 5, ]
        
        futile.logger::flog.info(glue::glue("Filtering on alt. DP >= 5: retaining {base::length(sample_vcf)} of {filter_ad} somatic variants."))
        
        # Filter on VAF.
        filter_vaf <- base::length(sample_vcf)
        sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$AF[,2] >= 0.025, ]
        
        futile.logger::flog.info(glue::glue("Filtering on VAF >= 0.025: retaining {base::length(sample_vcf)} of {filter_vaf} somatic variants."))
        
        # Generate a VRanges containing all relevant information from the sample_vcf.
        sample_vranges <- VariantAnnotation::VRanges(
            seqname = GenomicRanges::seqnames(sample_vcf),
            ranges = IRanges::IRanges(GenomicRanges::start(sample_vcf), GenomicRanges::end(sample_vcf)),
            ref = VariantAnnotation::ref(sample_vcf),
            alt = VariantAnnotation::alt(sample_vcf),
            refDepth = VariantAnnotation::geno(sample_vcf)$AD[, , 1][,2],
            altDepth = VariantAnnotation::geno(sample_vcf)$AD[, , 2][,2],
            totalDepth = VariantAnnotation::geno(sample_vcf)$DP[, 2],
            sampleNames = gsub("-.*", "", basename(path_vcf)),
            
            # Strain-specific counts.
            UA_Ref = VariantAnnotation::geno(sample_vcf)$UA[, , 1][,2],
            H1_Ref = VariantAnnotation::geno(sample_vcf)$H1[, , 1][,2],
            H2_Ref = VariantAnnotation::geno(sample_vcf)$H2[, , 1][,2],
            UA_Alt = VariantAnnotation::geno(sample_vcf)$UA[, , 2][,2],
            H1_Alt = VariantAnnotation::geno(sample_vcf)$H1[, , 2][,2],
            H2_Alt = VariantAnnotation::geno(sample_vcf)$H2[, , 2][,2],
            
            # Scores.
            VAF = VariantAnnotation::geno(sample_vcf)$AF[,2]
        )
        
        GenomeInfoDb::genome(sample_vranges) <- "mm39"
        
        # Determine mutational type.
        sample_vranges$mutType <- "Other"
        sample_vranges$mutType <- ifelse(VariantAnnotation::isSNV(sample_vranges), "SNV", sample_vranges$mutType)
        sample_vranges$mutType <- ifelse(VariantAnnotation::isIndel(sample_vranges), "InDel", sample_vranges$mutType)
        sample_vranges$mutType <- ifelse(!VariantAnnotation::isSNV(sample_vranges) & VariantAnnotation::isSubstitution(sample_vranges), "MNV", sample_vranges$mutType)
        sample_vranges$mutType <- ifelse(VariantAnnotation::isDelins(sample_vranges), "DelIn", sample_vranges$mutType)
        
        # Determine allelic origin of mutant allele.
        sample_vranges$origin_mutant <- determine_origin_mutual(sample_vranges$H1_Alt, sample_vranges$H2_Alt)
        
        # Return clean-up VRanges.
        return(sample_vranges)
    })
    
    # Convert to VRangesList. ----
    base::names(data_somaticvariants) <- base::vapply(data_somaticvariants, function(x) levels(Biobase::sampleNames(x)), "")
    data_somaticvariants <- VariantAnnotation::VRangesList(data_somaticvariants)
    
    return(data_somaticvariants)
}