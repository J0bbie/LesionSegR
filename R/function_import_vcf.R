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
import_vcf <- function(x) {
    # Input validation ---

    checkmate::assertFile(x)

    # Import and clean-up VCF. ---
    data_somaticvariants <- future.apply::future_lapply(x, function(path_vcf){

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
        sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$DP[, 1] >= 10, ]

        futile.logger::flog.info(glue::glue("Filtering on DP >= 10: retaining {base::length(sample_vcf)} of {filter_dp} somatic variants."))

        # Filter on AD.
        filter_ad <- base::length(sample_vcf)
        sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$AD[, , 2] >= 5, ]

        futile.logger::flog.info(glue::glue("Filtering on alt. DP >= 5: retaining {base::length(sample_vcf)} of {filter_ad} somatic variants."))

        # Filter on VAF.
        filter_vaf <- base::length(sample_vcf)
        sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$AF[, 1] >= 0.025, ]

        futile.logger::flog.info(glue::glue("Filtering on VAF >= 0.025: retaining {base::length(sample_vcf)} of {filter_vaf} somatic variants."))

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
        VariantAnnotation::info(sample_vcf) <- base::suppressWarnings(S4Vectors::cbind(VariantAnnotation::info(sample_vcf), split_csq))

        # Generate a VRanges containing all relevant information from the sample_vcf.
        sample_vranges <- VariantAnnotation::VRanges(
            seqname = GenomicRanges::seqnames(sample_vcf),
            ranges = IRanges::IRanges(GenomicRanges::start(sample_vcf), GenomicRanges::end(sample_vcf)),
            ref = VariantAnnotation::ref(sample_vcf),
            alt = VariantAnnotation::alt(sample_vcf),
            refDepth = VariantAnnotation::geno(sample_vcf)$AD[, , 1],
            altDepth = VariantAnnotation::geno(sample_vcf)$AD[, , 2],
            totalDepth = VariantAnnotation::geno(sample_vcf)$DP[, 1],
            sampleNames = VariantAnnotation::header(sample_vcf)@samples,

            # Strain-specific counts.
            UA_Ref = VariantAnnotation::geno(sample_vcf)$UA[, , 1],
            G1_Ref = VariantAnnotation::geno(sample_vcf)$G1[, , 1],
            G2_Ref = VariantAnnotation::geno(sample_vcf)$G2[, , 1],
            UA_Alt = VariantAnnotation::geno(sample_vcf)$UA[, , 2],
            G1_Alt = VariantAnnotation::geno(sample_vcf)$G1[, , 2],
            G2_Alt = VariantAnnotation::geno(sample_vcf)$G2[, , 2],

            # Scores.
            VAF = VariantAnnotation::geno(sample_vcf)$AF[, 1],

            # Annotations.
            Consequence = gsub("&.*", "", VariantAnnotation::info(sample_vcf)$Consequence),
            ConsequenceAll = VariantAnnotation::info(sample_vcf)$Consequence,
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

        # Add sample name.
        sample_vranges$sample <- base::basename(path_vcf) %>% stringr::str_remove("_withStrainCounts\\.vcf\\.gz$")

        # Determine allelic origin of mutant allele.
        sample_vranges$origin_mutant <- determine_origin_mutual(sample_vranges$G1_Alt, sample_vranges$G2_Alt)

        # Return clean-up VRanges.
        return(sample_vranges)
    })

    # Convert to VRangesList. ----
    base::names(data_somaticvariants) <- base::vapply(data_somaticvariants, function(x) unique(x$sample), "")
    data_somaticvariants <- VariantAnnotation::VRangesList(data_somaticvariants)

    return(data_somaticvariants)

}
