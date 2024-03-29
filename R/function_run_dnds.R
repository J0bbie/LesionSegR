#' @title Perform dN/dS analysis.
#' @details The database can be generated using the following code:
#'
#' ENSEMBLv110_GRCm39 <- readr::read_tsv("/omics/groups/OE0538/internal/projects/sharedData/GRCm39/annotation/GRCm39_dndscv_mart_export.txt") %>%
#'     dplyr::filter(!grepl("_", `Chromosome/scaffold name`)) %>%
#'     dplyr::filter(`CDS Length` %% 3 == 0) %>%
#'     dplyr::filter(!is.na(`Genomic coding start`)) %>%
#     dplyr::filter(!grepl("^AC[0-9]", `Gene name`)) %>%
#'     dplyr::filter(!grepl("^AP00|^RP11-", `Gene name`)) %>%
#'     dplyr::filter(!grepl("\\.", `Gene name`)) %>%
#'    dplyr::select(
#'         "gene.id" = `Gene stable ID`,
#'         "gene.name" = `Gene name`,
#'         "cds.id" = `Protein stable ID`,
#'         "chr" = `Chromosome/scaffold name`,
#'         "chr.coding.start" = `Genomic coding start`,
#'         "chr.coding.end" = `Genomic coding end`,
#'         "cds.start" = `CDS start`,
#'         "cds.end" = `CDS end`,
#'         "length" = `CDS Length`,
#'         "strand" = Strand
#'     ) %>%
#' dplyr::mutate(chr = paste0("chr", chr))
#' write.table(ENSEMBLv110_GRCm39, file = "~/test/mart_export_filtered.txt", row.names = FALSE, quote = FALSE)
#' pathCDS = '~/test/mart_export_filtered.txt'
#' pathFasta = '/omics/groups/OE0538/internal/projects/sharedData/GRCm39/genome/GRCm39.primary_assembly.genome.fa'
#' dndscv::buildref(cdsfile = pathCDS, genomefile = pathFasta, outfile = '~/test/refCDS_ENSEMBLv109_GRCm39.rda', excludechrs='chrMT', useids = TRUE)
#' @param data_muts (tibble): tibble containing the mutations of all samples.
#' @param path_db (character): Path to the dN/dS database.
#' @return (list): Tibble with
#' @importFrom dplyr %>%
#' @export
run_dnds <- function(data_muts, path_db = "/omics/groups/OE0538/internal/projects/sharedData/GRCm39/annotation/refCDS_ENSEMBLv110_GRCm39.rda") {
    # Input validation --------------------------------------------------------

    checkmate::assertTibble(data_muts)

    futile.logger::flog.info(glue::glue("Performing dN/dS analysis on {dplyr::n_distinct(data_muts$sample)} unique samples.\nThis can take some minutes."))

    # Perform dN/dS -----------------------------------------------------------

    # Convert mutations to data.frame and remove chr prefix.
    data_muts_df <- data.frame(
        sampleID = data_muts$sample,
        chr = data_muts$seqnames,
        pos = as.integer(data_muts$start),
        ref = data_muts$ref,
        mut = data_muts$alt
    )

    # Perform dN/dS algorithm.
    output_dnds <- dndscv::dndscv(data_muts_df, refdb = path_db, outp = 3)

    # Remove large (unused) annotation database.
    output_dnds$annotmuts <- NULL


    # Combine final results ---------------------------------------------------

    output_dnds$finalOutput <- output_dnds$sel_cv %>%
        dplyr::left_join(output_dnds$sel_loc %>% dplyr::select(gene_name, qall_loc, qmis_loc), by = c("gene_name" = "gene_name")) %>%
        dplyr::filter(qglobal_cv <= 0.1 | qallsubs_cv <= 0.1 | qtrunc_cv <= 0.1 | qmis_cv <= 0.1 | qall_loc <= 0.1 | qmis_loc <= 0.1) %>%
        dplyr::mutate(SYMBOL = gsub(".*:", "", gene_name), ENSEMBL = gsub(":.*", "", gene_name), gene_name = NULL)

    # Return statement --------------------------------------------------------

    return(output_dnds)
}
