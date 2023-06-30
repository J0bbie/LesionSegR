library(dplyr)

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
    sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$DP[, 1] >= 10, ]

    sprintf("\tFiltering on DP >= 10: retaining %s of %s somatic variants.", base::length(sample_vcf), filter_dp) %>% ParallelLogger::logInfo()

    # Filter on AD.
    filter_ad <- base::length(sample_vcf)
    sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$AD[, , 2] >= 5, ]

    sprintf("\tFiltering on alt. DP >= 5: retaining %s of %s somatic variants.", base::length(sample_vcf), filter_ad) %>% ParallelLogger::logInfo()

    # Filter on VAF.
    filter_vaf <- base::length(sample_vcf)
    sample_vcf <- sample_vcf[VariantAnnotation::geno(sample_vcf)$AF[, 1] >= 0.025, ]

    sprintf("\tFiltering on VAF >= 0.025: retaining %s of %s somatic variants.", base::length(sample_vcf), filter_vaf) %>% ParallelLogger::logInfo()

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
    sample_vranges$BIOTYPE <- factor(sample_vranges$BIOTYPE)

    # Add sample name.
    sample_vranges$sample <- base::basename(path_vcf) %>% stringr::str_remove("_withStrainCounts\\.vcf\\.gz$")

    # Determine allelic origin of mutant allele.
    sample_vranges$origin_mutant <- determine_origin_mutual(sample_vranges$G1_Alt, sample_vranges$G2_Alt)

    # Return clean-up VRanges.
    return(sample_vranges)
}


# Determine allelic origin. ----

# Determine origin based on mutual exclusivity of G1 / G2 mutant reads.
determine_origin_mutual <- function(count_strain1, count_strain2) {
    tibble::tibble(count_strain1, count_strain2) %>%
        dplyr::mutate(
            assignment = dplyr::case_when(
                .default = "UA",
                count_strain1 <= 1 & count_strain2 <= 1 ~ "UA",
                (count_strain1 > 0 & count_strain2 > 0) ~ "UA",
                count_strain1 > count_strain2 ~ "G1",
                count_strain2 > count_strain1 ~ "G2"
            )
        ) %>%
        dplyr::pull(.data$assignment)
}

# Determine origin based on mutual abundance of G1 / G2 mutant reads.
determine_origin_ratio <- function(count_strain1, count_strain2) {
    tibble::tibble(count_strain1, count_strain2) %>%
        dplyr::mutate(
            assignment = dplyr::case_when(
                .default = "UA",
                count_strain1 <= 2 & count_strain2 <= 2 ~ "UA",
                count_strain1 / count_strain2 >= 2 ~ "G1",
                count_strain2 / count_strain1 >= 2 ~ "G2"
            )
        ) %>%
        dplyr::pull(.data$assignment)
}


# Determine origin based on mixture model using at least 3 reads.
determine_origin_mixture <- function(x) {
    # Subset to variants with at least 2 reads in either strains.
    x <- tibble::as_tibble(x) %>%
        dplyr::mutate(id = dplyr::row_number())

    x_filter <- x %>%
        dplyr::filter(.$G1 >= 2 | .$G2 >= 2)

    # Generate a random normal distribution.
    pseudo <- sample(seq(.5, 1, by = 0.01), replace = TRUE, size = length(x_filter))

    # Calculate the difference between G1 - G2 (with added random pseudocounts)
    x_filter <- x_filter %>%
        dplyr::mutate(difference_straincount = (log2(.$G1 + pseudo) - log2(.$G2 + pseudo)))

    # Run mixture model.
    mixture_model <- mixtools::normalmixEM(x_filter$difference_straincount, lambda = .5, mu = c(-3.5, 3.5), sigma = c(1.2, 1.2))

    # Determine 99% threshold.
    threshold_strain1 <- max(x_filter$difference_straincount[mixture_model$posterior[, "comp.1"] >= 0.99])
    threshold_strain2 <- min(x_filter$difference_straincount[mixture_model$posterior[, "comp.2"] >= 0.99])

    x_filter <- x_filter %>%
        dplyr::mutate(
            origin = ifelse(
                .$difference_straincount <= threshold_strain1 | .$difference_straincount >= threshold_strain2,
                ifelse(
                    .$difference_straincount <= threshold_strain1, "G1", ifelse(.$difference_straincount >= threshold_strain2, "G2", "AU")
                ), "AU"
            )
        )

    x <- dplyr::left_join(x, x_filter[, c("id", "origin")], by = "id")

    return(x$origin)
}

# Import and clean flagstats. ----

read_flagstats <- function(x) {
    readr::read_csv(x, col_names = FALSE, show_col_types = FALSE) %>%
        dplyr::mutate(
            X1 = base::trimws(X1),
            X2 = as.integer(base::trimws(gsub(":.*", "", X2))),
            sample = basename(x) %>% stringr::str_replace_all(pattern = "_.*", replacement = "")
        ) %>%
        dplyr::select(variable = X1, value = X2, sample)
}

# Import and clean SNPsplit reporting yamls. ----
read_snpsplit <- function(x) {
    tibble::as_tibble(yaml::read_yaml(x)$Tagging) %>%
        dplyr::mutate(
            sample = basename(x) %>% stringr::str_replace_all(pattern = "_.*", replacement = "")
        )
}

# Determine genome-wide TMB ----

determine_mutational_burden <- function(x) {
    # Number of mappable ATCG in reference genome (GRCm39).
    tibble::tibble(
        TMB = length(x) / (2649938115 / 1E6),
        totalMutations = length(x),
        totalMutationsCoding = length(x[x$HGVSp != "", ]),
        totalG1 = sum(x$origin_mutant == "G1"),
        totalG2 = sum(x$origin_mutant == "G2"),
        totalUA = sum(x$origin_mutant == "UA"),
        sample = unique(x$sample)
    )
}


# Function - Generate 96-context matrix. ----
generate_mutmatrices_96 <- function(x) {
    # Convert mutations to correct GRanges for input into MutationalPatterns.
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
        MutationalPatterns::mut_matrix(., ref_genome = "BSgenome.Mmusculus.UCSC.mm39")

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
