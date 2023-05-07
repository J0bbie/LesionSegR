library(dplyr)
library(ParallelLogger)

. <- NULL

# Import and clean-up VCF. ----

import_vcf <- function(path_vcf) {

    sprintf("\tImporting VCF: %s", path_vcf) %>% ParallelLogger::logInfo()

    # Read in the VCF file.
    sample_vcf <- VariantAnnotation::readVcf(path_vcf, genome = "GRCm39")

    # Expand multi-allelic hits to separate rows.
    sample_vcf <- VariantAnnotation::expand(sample_vcf)

    # Filter VCF for PASS-only.
    filter_pass <- base::length(sample_vcf)
    sample_vcf <- sample_vcf[sample_vcf@fixed$FILTER == "PASS", ]

    sprintf("\tFiltering PASS-only: retaining %s of %s somatic variants.", base::length(sample_vcf), filter_pass) %>% ParallelLogger::logInfo()

    # Filter on DP.
    filter_dp <- base::length(sample_vcf)
    sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$DP[, 1] >= 5, ]

    sprintf("\tFiltering on DP >= 5: retaining %s of %s somatic variants.", base::length(sample_vcf), filter_dp) %>% ParallelLogger::logInfo()

    # Clean-up VEP annotations. ----
    sprintf("\tConverting and cleaning annotations") %>% ParallelLogger::logTrace()
    if (is.null(VariantAnnotation::info(sample_vcf)$CSQ)) stop(sprintf("This file does not contain a CSQ (annotation / VEP) column:\t%s", path_vcf))

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
        UA = VariantAnnotation::geno(sample_vcf)$UA[, , 2],
        G1 = VariantAnnotation::geno(sample_vcf)$G1[, , 2],
        G2 = VariantAnnotation::geno(sample_vcf)$G2[, , 2],

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
    sample_vranges$BIOTYPE <- factor(sample_vranges$BIOTYPE)

    # Add sample name.
    sample_vranges$sampleName <- base::basename(path_vcf) %>% stringr::str_remove("_withStrainCounts\\.vcf\\.gz$")

    # Return clean-up VRanges.
    return(sample_vranges)
}


# Determine allelic origin. ----

# Determine origin based on X2 test for given 50/50 probabilities.
determine_origin_chisq <- function(count_strain1, count_strain2) {
    if (count_strain1 == 0 && count_strain2 == 0) {
        return("AU")
    }

    p <- chisq.test(data.frame(G1 = count_strain1, G2 = count_strain2), p = c(0.5, 0.5))$p.value
    if (p < 0.1) {
        return(ifelse(count_strain1 > count_strain2, "G1", "G2"))
    }

    return("AU")
}

# Determine origin based on ratio of G1 / G2.
determine_origin_ratio <- function(count_strain1, count_strain2) {
    if (count_strain1 < 2 && count_strain2 < 2) {
        return("AU")
    }

    if (count_strain1 / count_strain2 >= 2) {
        return("G1")
    }
    if (count_strain2 / count_strain1 >= 2) {
        return("G2")
    }

    return("AU")
}

# Determine genome-wide TMB
calculate_tmb <- function(x) {
    # Number of mappable ATCG in reference genome (GRCm39).
    tibble::tibble(
        TMB = length(x) / (2649938115 / 1E6),
        sampleName = unique(x$sampleName)
    )
}

# Perform mutatational signature analysis (mm10).
determine_knownsignatures <- function(x) {

    # Retrieve the COSMIC v3.2 signatures. ---

    sprintf("\tRetrieving COSMIC (v3.2) signature matrices (mm10).") %>% ParallelLogger::logInfo()

    sigs_cosmic <- list()
    sigs_cosmic$sbs <- MutationalPatterns::get_known_signatures(muttype = "snv", source = "COSMIC_v3.2", genome = "mm10", incl_poss_artifacts = FALSE)
    sigs_cosmic$indel <- MutationalPatterns::get_known_signatures(muttype = "indel", source = "COSMIC_v3.2", genome = "GRCh37", incl_poss_artifacts = FALSE)
    sigs_cosmic$dbs <- MutationalPatterns::get_known_signatures(muttype = "dbs", source = "COSMIC_v3.2", genome = "mm10", incl_poss_artifacts = FALSE)

    # Convert mutations to correct GRanges for input into MutationalPatterns.
    convert_muts <- function(x, dbs = FALSE) {
        S4Vectors::mcols(x) <- S4Vectors::DataFrame(sample = x$sampleName)

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

    # Calculate rel. contribution.
    clean_sigs <- function(x) {
        # Convert to tibble.
        tibble::as_tibble(base::sweep(x, 2, base::colSums(x), "/") * 100, rownames = "Signature") %>%
            # Melt.
            tidyr::pivot_longer(cols = !dplyr::contains("Signature"), names_to = "sampleId", values_to = "relContribution")
    }

    # Retrieve mutational motifs. ---

    sprintf("\tConverting GRangesLists into mutational matrices and performing mutational signature analysis.") %>% ParallelLogger::logInfo()

    data_mutmatrices <- list()
    data_fittedsigs <- list()

    ## SBS ---
    data_mutmatrices$sbs <- convert_muts(x[x$mutType == "SNV", ]) %>%
        MutationalPatterns::mut_matrix(., ref_genome = "BSgenome.Mmusculus.UCSC.mm39")

    data_fittedsigs$sbs <- MutationalPatterns::fit_to_signatures_strict(mut_matrix = data_mutmatrices$sbs, signatures = sigs_cosmic$sbs)
    data_fittedsigs$sbs$relativeContribution <- clean_sigs(data_fittedsigs$sbs$fit_res$contribution)
    data_fittedsigs$sbs$mutMatrix <- data_mutmatrices$sbs

    ## InDel ---
    data_mutmatrices$indel <- convert_muts(x[x$mutType == "InDel", ]) %>%
        MutationalPatterns::get_indel_context(., ref_genome = "BSgenome.Mmusculus.UCSC.mm39") %>%
        MutationalPatterns::count_indel_contexts(.)

    data_fittedsigs$indel <- MutationalPatterns::fit_to_signatures_strict(mut_matrix = data_mutmatrices$indel, signatures = sigs_cosmic$indel)
    data_fittedsigs$indel$relativeContribution <- clean_sigs(data_fittedsigs$indel$fit_res$contribution)
    data_fittedsigs$indel$mutMatrix <- data_mutmatrices$indel

    ## DBS ---
    if (sum(x$mutType == "MNV") != 0) {
        data_mutmatrices$dbs <- convert_muts(x[x$mutType == "MNV" & base::nchar(VariantAnnotation::ref(x)) == 2 & base::nchar(VariantAnnotation::alt(x)) == 2], dbs = TRUE) %>%
            MutationalPatterns::get_dbs_context(.) %>%
            MutationalPatterns::count_dbs_contexts(.)

        data_fittedsigs$dbs <- MutationalPatterns::fit_to_signatures_strict(mut_matrix = data_mutmatrices$dbs, signatures = sigs_cosmic$dbs)

        data_fittedsigs$dbs$mutMatrix <- data_mutmatrices$dbs
        data_fittedsigs$dbs$relativeContribution <- clean_sigs(data_fittedsigs$dbs$fit_res$contribution)
    } else {
        data_mutmatrices$dbs <- NULL
    }

    # Return cleaned signatures. ----
    return(data_fittedsigs)
}
