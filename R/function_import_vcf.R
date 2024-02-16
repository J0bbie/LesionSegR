#' @title Importing VCF files.
#' @description Imports VCF files from the snakemake-lesionsegregation workflow.
#' Applies several filtering steps and determines the likely strain-origin of the somatic variants.
#' @param x (character): Path to one VCF file(s).
#' @examples
#' 1 + 1
#' @return (tibble): Tibble containing the (filtered) somatic variants.
#' and strain-specific tagging.
#' @importFrom dplyr %>%
#' @export
import_vcf <- function(path_vcf) {
    # Input validation ---
    
    checkmate::assertFile(path_vcf, access = "r")
    
    # Import and clean-up VCF. ---
    
    futile.logger::flog.info(glue::glue("Importing VCF: {path_vcf}"))
    
    # Read in the VCF file.
    sample_vcf <- VariantAnnotation::readVcf(path_vcf, genome = "GRCm39", param = VariantAnnotation::ScanVcfParam(samples = 'TUMOR'))
    
    # Expand multi-allelic hits to separate rows.
    sample_vcf <- VariantAnnotation::expand(sample_vcf)
    
    # Filter VCF for PASS-only.
    filter_pass <- base::length(sample_vcf)
    sample_vcf <- sample_vcf[sample_vcf@fixed$FILTER == "PASS", ]
    
    futile.logger::flog.info(glue::glue("Filtering PASS-only: retaining {base::length(sample_vcf)} of {filter_pass} somatic variants."))
    
    # Filter on DP.
    filter_dp <- base::length(sample_vcf)
    sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$DP[, 1] >= 10, ]
    
    futile.logger::flog.info(glue::glue("Filtering on DP >= 10: retaining {base::length(sample_vcf)} of {filter_dp} somatic variants."))
    
    # Filter on VAF.
    filter_vaf <- base::length(sample_vcf)
    sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$VAF >= 0.025, ]
    futile.logger::flog.info(glue::glue("Filtering on VAF >= 0.025: retaining {base::length(sample_vcf)} of {filter_vaf} somatic variants."))
    
    # Add reference and alt counts.
    suppressWarnings(VariantAnnotation::info(sample_vcf)$altDepth <- round(unlist(VariantAnnotation::geno(sample_vcf)$VAF) * VariantAnnotation::geno(sample_vcf)$DP[, 1]))
    suppressWarnings(VariantAnnotation::info(sample_vcf)$refDepth <- VariantAnnotation::geno(sample_vcf)$DP[, 1] - VariantAnnotation::info(sample_vcf)$altDepth)
    
    filter_alt <- base::length(sample_vcf)
    sample_vcf <- sample_vcf[VariantAnnotation::info(sample_vcf)$altDepth >= 5, ]
    futile.logger::flog.info(glue::glue("Filtering on alt. count >= 5: retaining {base::length(sample_vcf)} of {filter_alt} somatic variants."))
    
    # Clean-up VEP annotations. ----
    futile.logger::flog.info("Converting and cleaning annotations")
    if (is.null(VariantAnnotation::info(sample_vcf)$CSQ)) stop(sprintf("This file does not contain a CSQ (annotation / VEP) column:\t{path_vcf}"))
    
    # Convert annotation to tibble.
    split_csq <- tibble::as_tibble(unlist(VariantAnnotation::info(sample_vcf)$CSQ)) %>%
        tidyr::separate("value",
                        into =
                            c(
                                "Allele",
                                "Consequence",
                                "IMPACT",
                                "SYMBOL",
                                "Gene",
                                "Feature_type",
                                "Feature",
                                "BIOTYPE",
                                "EXON",
                                "INTRON",
                                "HGVSc",
                                "HGVSp",
                                "cDNA_position",
                                "CDS_position",
                                "Protein_position",
                                "Amino_acids",
                                "Codons",
                                "Existing_variation",
                                "DISTANCE",
                                "STRAND",
                                "FLAGS",
                                "SYMBOL_SOURCE",
                                "HGNC_ID",
                                "TSL",
                                "HGVS_OFFSET"
                            ),
                        sep = "\\|"
        )
    
    # Add the split columns to the VCF.
    base::suppressWarnings(VariantAnnotation::info(sample_vcf) <- S4Vectors::cbind(VariantAnnotation::info(sample_vcf), split_csq))
    
    # Generate a tibble containing all relevant information from the sample_vcf.
    if(length(sample_vcf) > 0){
        sample_vranges <- VariantAnnotation::VRanges(
            seqname = GenomicRanges::seqnames(sample_vcf),
            ranges = IRanges::IRanges(GenomicRanges::start(sample_vcf), GenomicRanges::end(sample_vcf)),
            ref = VariantAnnotation::ref(sample_vcf),
            alt = VariantAnnotation::alt(sample_vcf),
            refDepth = VariantAnnotation::info(sample_vcf)$refDepth,
            altDepth = VariantAnnotation::info(sample_vcf)$altDepth,
            totalDepth = VariantAnnotation::geno(sample_vcf)$DP[, 1],
            sample = base::gsub("_.*", "", basename(path_vcf)),
            
            # Strain-specific counts.
            UA_Ref = VariantAnnotation::geno(sample_vcf)$UA[, , 1],
            H1_Ref = VariantAnnotation::geno(sample_vcf)$H1[, , 1],
            H2_Ref = VariantAnnotation::geno(sample_vcf)$H2[, , 1],
            UA_Alt = VariantAnnotation::geno(sample_vcf)$UA[, , 2],
            H1_Alt = VariantAnnotation::geno(sample_vcf)$H1[, , 2],
            H2_Alt = VariantAnnotation::geno(sample_vcf)$H2[, , 2],
            
            # Scores.
            VAF = VariantAnnotation::geno(sample_vcf)$VAF,
            
            # Annotations.
            Consequence = factor(gsub("&.*", "", VariantAnnotation::info(sample_vcf)$Consequence)),
            ConsequenceAll = factor(VariantAnnotation::info(sample_vcf)$Consequence),
            SYMBOL = VariantAnnotation::info(sample_vcf)$SYMBOL,
            ENSEMBL = VariantAnnotation::info(sample_vcf)$Gene,
            BIOTYPE = VariantAnnotation::info(sample_vcf)$BIOTYPE,
            EXON = VariantAnnotation::info(sample_vcf)$EXON,
            INTRON = VariantAnnotation::info(sample_vcf)$INTRON,
            HGVSc = VariantAnnotation::info(sample_vcf)$HGVSc,
            HGVSp = VariantAnnotation::info(sample_vcf)$HGVSp,
            IMPACT = VariantAnnotation::info(sample_vcf)$IMPACT
        )
        
        
        GenomeInfoDb::genome(sample_vranges) <- "mm39"
        
        # Determine mutational type.
        sample_vranges$mutType <- "Other"
        sample_vranges$mutType <- ifelse(VariantAnnotation::isSNV(sample_vranges), "SNV", sample_vranges$mutType)
        sample_vranges$mutType <- ifelse(VariantAnnotation::isIndel(sample_vranges), "InDel", sample_vranges$mutType)
        sample_vranges$mutType <- ifelse(!VariantAnnotation::isSNV(sample_vranges) & VariantAnnotation::isSubstitution(sample_vranges), "MNV", sample_vranges$mutType)
        sample_vranges$mutType <- ifelse(VariantAnnotation::isDelins(sample_vranges), "DelIn", sample_vranges$mutType)
        
        # Convert several columns to factor to save size.
        sample_vranges$mutType <- base::factor(sample_vranges$mutType)
        sample_vranges$BIOTYPE <- base::factor(sample_vranges$BIOTYPE)
        
        # Determine allelic origin of mutant allele.
        sample_vranges$origin_mutant <- determine_origin_mutual(sample_vranges$H1_Alt, sample_vranges$H2_Alt)
        
        # Return clean-up tibble.
        sample_tibble <- tibble::as_tibble(sample_vranges)
        sample_tibble$sampleNames <- NULL
        
        return(sample_tibble)
    }else{
        return(NULL)
    }
}
