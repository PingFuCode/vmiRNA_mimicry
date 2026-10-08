# Figure 4A: HIV-1 shared-target and virus-host interaction networks
# Adapted from:
# D:/0work/0wholetransctiptome/2sRNAminic/13Infect_mimic_correlation/Time/Network.R

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggVennDiagram)

target_file <- paste0(
  "D:/0work/0wholetransctiptome/2sRNAminic/",
  "3mRNATarget/0006mer_target_with_hsaID.csv"
)
virus_ppi_file <- paste0(
  "D:/0work/0wholetransctiptome/2sRNAminic/",
  "16StringDB/virus_human_ppi_all.txt"
)
human_ppi_file <- paste0(
  "D:/0work/0wholetransctiptome/2sRNAminic/",
  "16StringDB/string_interactions_short.tsv"
)
output_dir <- "D:/0work/paper/0wtite/miRNA/V7_submit/code"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tar <- readr::read_csv(target_file, show_col_types = FALSE)
virus_ppi_columns <- c(
  "GeneIDA", "GeneIDB", "UniprotACA", "UniprotACB", "GeneSymbolA",
  "GeneSymbolB", "Sources", "Types", "Methods", "Score", "TaxID",
  "OrganismName", "SpeciesTaxID", "SpeciesName", "FamilyTaxID",
  "FamilyName", "Unused"
)
# The source rows contain one trailing field that is absent from the header.
# Reading all 17 fields explicitly prevents column-shift/parsing warnings.
virus_ppi <- readr::read_tsv(
  virus_ppi_file,
  col_names = virus_ppi_columns,
  skip = 1,
  show_col_types = FALSE
)
human_ppi <- readr::read_tsv(human_ppi_file, show_col_types = FALSE)

clean_gene <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x[is.na(x) | x == "" | x == "-"] <- NA_character_
  x
}

split_target_genes <- function(x) {
  clean_gene(unlist(stringr::str_split(x[!is.na(x)], ";\\s*"))) |>
    stats::na.omit() |>
    unique()
}

# -----------------------------------------------------------------------------
# A. Large HIV-1 miRNA-target network
# -----------------------------------------------------------------------------

tar_hiv1 <- tar %>%
  dplyr::filter(abbreviation == "HIV-1")

large_edges <- tar_hiv1 %>%
  dplyr::mutate(Intersection_Genes = stringr::str_split(Intersection_Genes, ";")) %>%
  tidyr::unnest(Intersection_Genes) %>%
  dplyr::select(from = miRNAID, to = Intersection_Genes) %>%
  dplyr::mutate(
    edge_source = "Virus",
    miRNA_category = "human-mimicry vmiRNAs"
  ) %>%
  dplyr::bind_rows(
    tar_hiv1 %>%
      dplyr::mutate(Intersection_Genes = stringr::str_split(Intersection_Genes, ";")) %>%
      tidyr::unnest(Intersection_Genes) %>%
      dplyr::select(from = miRNA, to = Intersection_Genes) %>%
      dplyr::mutate(
        edge_source = "Human",
        miRNA_category = "virus-like hmiRNAs"
      )
  ) %>%
  dplyr::mutate(
    from = clean_gene(from),
    to = clean_gene(to)
  ) %>%
  dplyr::filter(!is.na(from), !is.na(to))

large_graph <- igraph::graph_from_data_frame(large_edges, directed = TRUE)
large_names <- igraph::V(large_graph)$name
virus_like_hmirnas <- unique(clean_gene(tar_hiv1$miRNA))
human_mimicry_vmirnas <- unique(clean_gene(tar_hiv1$miRNAID))

large_node_type <- ifelse(
  large_names %in% c(virus_like_hmirnas, human_mimicry_vmirnas),
  "miRNA",
  "gene"
)
large_degree_out <- igraph::degree(large_graph, mode = "out")
large_degree_in <- igraph::degree(large_graph, mode = "in")
large_degree <- ifelse(
  large_node_type == "miRNA",
  large_degree_out,
  large_degree_in
)

gene_index <- which(large_node_type == "gene")
hub_genes <- large_names[gene_index][
  order(large_degree[gene_index], decreasing = TRUE)
][seq_len(min(10, length(gene_index)))]

black_genes <- c(
  "HLA-C", "TADA2A", "GRIN2B", "OCLN", "PPIB", "GRIN2D", "HERC2", "JUN",
  "ACTG1", "MLLT1", "KAT6B", "PRPF8", "RPS16", "FASN", "PRRC2A", "TRIB1",
  "HNRNPL", "DHX36", "TAP2", "AATF", "CLSTN3", "YWHAE", "PRPF31", "HLA-B",
  "RPN2", "SUPT4H1", "BRD2", "SIRT2", "IGFBP4", "HYOU1", "PRKDC", "ACACA",
  "HNRNPH1", "RBM14", "KCNH3", "DAXX", "PSMB3", "TNF", "ARFGEF1", "LSM2",
  "PLIN3", "FDFT1", "ERP44", "RPS2", "SRC", "SSRP1", "CRIPT", "CTSB", "B2M",
  "CCL11", "RNH1", "CTTN", "PTK6", "TAPBP", "EEF1D", "ATP5F1B", "FLOT1",
  "SURF4", "E2F4", "GTF2H2", "NAP1L4", "XRCC2", "DCAF11"
)
black_genes_only <- setdiff(black_genes, hub_genes)

large_nodes <- tibble::tibble(
  name = large_names,
  node_type = large_node_type,
  degree = as.numeric(large_degree),
  node_category = dplyr::case_when(
    name %in% hub_genes ~ "hub",
    name %in% black_genes_only ~ "black_genes",
    name %in% virus_like_hmirnas ~ "virus-like hmiRNAs",
    name %in% human_mimicry_vmirnas ~ "human-mimicry vmiRNAs",
    TRUE ~ "gene"
  ),
  node_shape = ifelse(node_type == "miRNA", 24L, 21L),
  plot_size = dplyr::case_when(
    node_type == "miRNA" ~ 6,
    name %in% hub_genes ~ degree * 1.5,
    TRUE ~ degree
  ),
  is_hub = name %in% hub_genes,
  is_gene_of_interest = name %in% black_genes
)

readr::write_csv(
  large_edges,
  file.path(output_dir, "Figure4A_large_network_edges.csv")
)
readr::write_csv(
  large_nodes,
  file.path(output_dir, "Figure4A_large_network_nodes.csv")
)

large_tg <- tidygraph::tbl_graph(
  nodes = large_nodes,
  edges = large_edges,
  directed = TRUE,
  node_key = "name"
)

set.seed(123)
large_layout <- ggraph::create_layout(
  large_tg,
  layout = "fr",
  niter = 3000,
  weights = rep(3, nrow(large_edges))
)

# Pull the outer 15% of nodes towards the network body while preserving their
# angular positions. This reduces the scattered appearance of peripheral nodes.
layout_center_x <- stats::median(large_layout$x)
layout_center_y <- stats::median(large_layout$y)
layout_dx <- large_layout$x - layout_center_x
layout_dy <- large_layout$y - layout_center_y
layout_radius <- sqrt(layout_dx^2 + layout_dy^2)
layout_cutoff <- as.numeric(stats::quantile(layout_radius, 0.75, na.rm = TRUE))
compressed_radius <- ifelse(
  layout_radius > layout_cutoff,
  layout_cutoff + 0.20 * (layout_radius - layout_cutoff),
  layout_radius
)
layout_scale <- ifelse(layout_radius == 0, 1, compressed_radius / layout_radius)
large_layout$x <- layout_center_x + layout_dx * layout_scale
large_layout$y <- layout_center_y + layout_dy * layout_scale

# Move disconnected satellite components closer to the largest component.
component_membership <- igraph::components(
  igraph::as.igraph(large_tg)
)$membership
largest_component <- as.integer(
  names(which.max(table(component_membership)))
)
main_index <- which(component_membership == largest_component)
main_center_x <- mean(large_layout$x[main_index])
main_center_y <- mean(large_layout$y[main_index])

for (component_i in setdiff(unique(component_membership), largest_component)) {
  component_index <- which(component_membership == component_i)
  component_center_x <- mean(large_layout$x[component_index])
  component_center_y <- mean(large_layout$y[component_index])

  large_layout$x[component_index] <-
    main_center_x +
    0.55 * (component_center_x - main_center_x) +
    0.80 * (large_layout$x[component_index] - component_center_x)
  large_layout$y[component_index] <-
    main_center_y +
    0.55 * (component_center_y - main_center_y) +
    0.80 * (large_layout$y[component_index] - component_center_y)
}

large_network_plot <- ggraph::ggraph(large_layout) +
  ggraph::geom_edge_link(
    ggplot2::aes(edge_linetype = edge_source),
    edge_colour = "grey75",
    edge_alpha = 0.55,
    edge_width = 0.4
  ) +
  ggraph::geom_node_point(
    ggplot2::aes(
      size = plot_size,
      shape = factor(node_shape),
      fill = node_category
    ),
    colour = "transparent",
    stroke = 0,
    alpha = 0.9
  ) +
  ggraph::geom_node_text(
    ggplot2::aes(label = ifelse(is_hub, name, "")),
    size = 3,
    repel = TRUE,
    fontface = "bold",
    colour = "black",
    bg.color = "white",
    bg.r = 0.1,
    max.overlaps = Inf,
    segment.size = 0.2,
    segment.color = "grey60"
  ) +
  ggraph::scale_edge_linetype_manual(
    values = c("Human" = "dashed", "Virus" = "solid"),
    breaks = c("Human", "Virus"),
    labels = c("Human miRNA", "Viral miRNA"),
    name = "Edge source"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "gene" = "#4682B4",
      "hub" = "#E41A1C",
      "black_genes" = "#000000",
      "human-mimicry vmiRNAs" = "#FFE494",
      "virus-like hmiRNAs" = "#4D4D4D"
    ),
    breaks = c(
      "gene", "hub", "black_genes",
      "human-mimicry vmiRNAs", "virus-like hmiRNAs"
    ),
    labels = c(
      "Other genes", "Hub genes (top 10 by degree)", "Genes of interest",
      "Human-mimicry vmiRNAs", "Virus-like hmiRNAs"
    ),
    name = "Node category"
  ) +
  ggplot2::scale_shape_manual(values = c("21" = 21, "24" = 24), guide = "none") +
  ggplot2::scale_size_continuous(range = c(2, 8), guide = "none") +
  ggplot2::theme_void() +
  ggplot2::theme(
    text = ggplot2::element_text(colour = "black"),
    legend.position = "right",
    plot.margin = ggplot2::margin(10, 10, 10, 10)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      order = 1,
      override.aes = list(
        shape = c(21, 21, 21, 24, 24),
        size = 5,
        colour = "transparent"
      )
    ),
    edge_linetype = ggplot2::guide_legend(order = 2)
  )

ggplot2::ggsave(
  file.path(output_dir, "Figure4A_large_network.pdf"),
  large_network_plot,
  width = 9,
  height = 6,
  device = grDevices::cairo_pdf
)

# -----------------------------------------------------------------------------
# B. Small HLA-C-centred network
# -----------------------------------------------------------------------------

small_hub_genes <- "HLA-C"

virus_to_hub <- virus_ppi %>%
  dplyr::filter(
    SpeciesName == "Human immunodeficiency virus 1",
    GeneSymbolA %in% small_hub_genes
  ) %>%
  dplyr::transmute(
    from = clean_gene(GeneSymbolB),
    to = clean_gene(GeneSymbolA),
    edge_type = "virus-hub"
  )

small_mirna_edges <- tar_hiv1 %>%
  dplyr::mutate(Intersection_Genes = stringr::str_split(Intersection_Genes, ";")) %>%
  tidyr::unnest(Intersection_Genes) %>%
  dplyr::filter(Intersection_Genes %in% small_hub_genes) %>%
  dplyr::bind_rows(
    dplyr::mutate(., from = miRNAID, to = Intersection_Genes, source = "virus"),
    dplyr::mutate(., from = miRNA, to = Intersection_Genes, source = "human")
  ) %>%
  dplyr::select(from, to, source) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    from = clean_gene(from),
    to = clean_gene(to),
    edge_type = dplyr::recode(
      source,
      virus = "virus-miRNA-target",
      human = "human-miRNA-target"
    )
  ) %>%
  dplyr::select(-source)

# The threshold is retained from the original small-network code.
small_human_ppi <- human_ppi %>%
  dplyr::filter(combined_score > 1, !is.na(combined_score)) %>%
  dplyr::filter(node1 %in% small_hub_genes | node2 %in% small_hub_genes) %>%
  dplyr::transmute(
    from = clean_gene(node1),
    to = clean_gene(node2),
    edge_type = as.character(ifelse(
      node1 %in% small_hub_genes & node2 %in% small_hub_genes,
      "hub-hub",
      "hub-human"
    ))
  ) %>%
  dplyr::distinct()

small_edges <- dplyr::bind_rows(
  virus_to_hub,
  small_human_ppi,
  small_mirna_edges
) %>%
  dplyr::filter(!is.na(from), !is.na(to), from != to) %>%
  dplyr::distinct()

small_mirna_nodes <- unique(small_mirna_edges$from)
small_nodes <- tibble::tibble(
  name = unique(c(small_edges$from, small_edges$to))
) %>%
  dplyr::mutate(
    node_category = dplyr::case_when(
      name %in% small_hub_genes ~ "hub",
      name %in% virus_to_hub$from ~ "viral protein",
      name %in% small_mirna_nodes & name %in% tar_hiv1$miRNAID ~ "viral miRNA",
      name %in% small_mirna_nodes & name %in% tar_hiv1$miRNA ~ "human miRNA",
      TRUE ~ "human protein"
    )
  )

if (nrow(small_edges) > 0) {
  small_graph <- igraph::graph_from_data_frame(
    small_edges,
    directed = FALSE,
    vertices = small_nodes
  )
  small_nodes$degree <- as.numeric(
    igraph::degree(small_graph)[match(small_nodes$name, names(igraph::degree(small_graph)))]
  )
} else {
  small_nodes$degree <- numeric(nrow(small_nodes))
}

readr::write_csv(
  small_edges,
  file.path(output_dir, "Figure4A_small_network_edges.csv")
)
readr::write_csv(
  small_nodes,
  file.path(output_dir, "Figure4A_small_network_nodes.csv")
)

# -----------------------------------------------------------------------------
# C. Venn source data and Fisher's exact tests for every virus
# -----------------------------------------------------------------------------

# Species names that changed between the miRNA table and the virus-host PPI file.
species_alias <- list(
  "HHV-6" = c("Human betaherpesvirus 6A", "Human betaherpesvirus 6B"),
  "HPyV-1" = "Human polyomavirus 1",
  "HPyV-2" = "Human polyomavirus 2",
  "MCPyV" = "Human polyomavirus 5",
  "TTV" = "Torque teno virus 1"
)

virus_key <- tar %>%
  dplyr::distinct(abbreviation, virus_species) %>%
  dplyr::arrange(abbreviation)

get_ppi_species <- function(abbreviation, virus_species) {
  aliases <- species_alias[[abbreviation]]
  if (is.null(aliases)) virus_species else aliases
}

all_target_genes <- split_target_genes(tar$Intersection_Genes)
all_interacting_human_genes <- unique(clean_gene(virus_ppi$GeneSymbolA))

# A common, explicit background is used for all viruses.
background_genes <- unique(stats::na.omit(c(
  clean_gene(human_ppi$node1),
  clean_gene(human_ppi$node2),
  all_target_genes,
  all_interacting_human_genes
)))
background_size <- length(background_genes)
background_definition <- paste(
  "Union of human STRING-network genes, all shared-target genes,",
  "and all virus-interacting human genes"
)

venn_membership_list <- vector("list", nrow(virus_key))
venn_summary_list <- vector("list", nrow(virus_key))
fisher_list <- vector("list", nrow(virus_key))

for (i in seq_len(nrow(virus_key))) {
  abbreviation_i <- virus_key$abbreviation[i]
  virus_species_i <- virus_key$virus_species[i]
  ppi_species_i <- get_ppi_species(abbreviation_i, virus_species_i)

  target_i <- tar %>%
    dplyr::filter(abbreviation == abbreviation_i) %>%
    dplyr::pull(Intersection_Genes) %>%
    split_target_genes()

  interacting_i <- virus_ppi %>%
    dplyr::filter(SpeciesName %in% ppi_species_i) %>%
    dplyr::pull(GeneSymbolA) %>%
    clean_gene() %>%
    stats::na.omit() %>%
    unique()

  membership_i <- tibble::tibble(
    abbreviation = abbreviation_i,
    virus_species = virus_species_i,
    matched_ppi_species = paste(ppi_species_i, collapse = "; "),
    gene = base::union(target_i, interacting_i)
  ) %>%
    dplyr::mutate(
      is_shared_target_gene = gene %in% target_i,
      is_virus_interacting_human_gene = gene %in% interacting_i,
      venn_region = dplyr::case_when(
        is_shared_target_gene & is_virus_interacting_human_gene ~ "Intersection",
        is_shared_target_gene ~ "Shared target only",
        TRUE ~ "Interacting human gene only"
      )
    )
  venn_membership_list[[i]] <- membership_i

  a <- length(intersect(target_i, interacting_i))
  b <- length(setdiff(target_i, interacting_i))
  c_count <- length(setdiff(interacting_i, target_i))
  d <- background_size - a - b - c_count

  venn_summary_list[[i]] <- tibble::tibble(
    abbreviation = abbreviation_i,
    virus_species = virus_species_i,
    matched_ppi_species = paste(ppi_species_i, collapse = "; "),
    shared_target_genes = length(target_i),
    interacting_human_genes = length(interacting_i),
    intersection_genes = a,
    shared_target_only = b,
    interacting_human_gene_only = c_count,
    union_genes = length(base::union(target_i, interacting_i)),
    jaccard_index = ifelse(
      length(base::union(target_i, interacting_i)) == 0,
      NA_real_,
      a / length(base::union(target_i, interacting_i))
    )
  )

  if (length(target_i) > 0 && length(interacting_i) > 0 && d >= 0) {
    contingency <- matrix(c(a, b, c_count, d), nrow = 2, byrow = TRUE)
    dimnames(contingency) <- list(
      shared_target = c("Yes", "No"),
      virus_interacting = c("Yes", "No")
    )
    fisher_i <- stats::fisher.test(contingency, alternative = "greater")

    fisher_list[[i]] <- tibble::tibble(
      abbreviation = abbreviation_i,
      virus_species = virus_species_i,
      matched_ppi_species = paste(ppi_species_i, collapse = "; "),
      target_yes_interaction_yes = a,
      target_yes_interaction_no = b,
      target_no_interaction_yes = c_count,
      target_no_interaction_no = d,
      odds_ratio = unname(fisher_i$estimate),
      conf_low = fisher_i$conf.int[1],
      conf_high = fisher_i$conf.int[2],
      p_value = fisher_i$p.value,
      test = "One-sided Fisher's exact test (greater)",
      status = "Tested"
    )
  } else {
    fisher_list[[i]] <- tibble::tibble(
      abbreviation = abbreviation_i,
      virus_species = virus_species_i,
      matched_ppi_species = paste(ppi_species_i, collapse = "; "),
      target_yes_interaction_yes = a,
      target_yes_interaction_no = b,
      target_no_interaction_yes = c_count,
      target_no_interaction_no = d,
      odds_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_,
      test = "One-sided Fisher's exact test (greater)",
      status = ifelse(
        length(interacting_i) == 0,
        "Not tested: no matched virus-host PPI genes",
        "Not tested: empty target set or invalid table"
      )
    )
  }
}

venn_membership_all <- dplyr::bind_rows(venn_membership_list)
venn_summary_all <- dplyr::bind_rows(venn_summary_list)
fisher_results <- dplyr::bind_rows(fisher_list) %>%
  dplyr::mutate(
    p_adjust_BH = stats::p.adjust(p_value, method = "BH"),
    significant_raw_p_0.05 = !is.na(p_value) & p_value < 0.05,
    significant_BH_0.05 = !is.na(p_adjust_BH) & p_adjust_BH < 0.05,
    background_size = background_size,
    background_definition = background_definition
  ) %>%
  dplyr::arrange(p_value)

readr::write_csv(
  venn_membership_all,
  file.path(output_dir, "Figure4A_all_viruses_venn_gene_membership.csv")
)
readr::write_csv(
  venn_summary_all,
  file.path(output_dir, "Figure4A_all_viruses_venn_summary.csv")
)
readr::write_csv(
  dplyr::filter(venn_membership_all, abbreviation == "HIV-1"),
  file.path(output_dir, "Figure4A_HIV1_venn_gene_membership.csv")
)
readr::write_csv(
  dplyr::filter(venn_summary_all, abbreviation == "HIV-1"),
  file.path(output_dir, "Figure4A_HIV1_venn_summary.csv")
)
readr::write_csv(
  fisher_results,
  file.path(output_dir, "Figure4A_all_viruses_fisher_results.csv")
)
readr::write_lines(
  sort(background_genes),
  file.path(output_dir, "Figure4A_fisher_background_genes.txt")
)

hiv1_targets <- tar %>%
  dplyr::filter(abbreviation == "HIV-1") %>%
  dplyr::pull(Intersection_Genes) %>%
  split_target_genes()
hiv1_interacting <- virus_ppi %>%
  dplyr::filter(SpeciesName == "Human immunodeficiency virus 1") %>%
  dplyr::pull(GeneSymbolA) %>%
  clean_gene() %>%
  stats::na.omit() %>%
  unique()

hiv1_venn_plot <- ggVennDiagram::ggVennDiagram(
  list(
    "HIV-1 shared target genes" = hiv1_targets,
    "HIV-1 interacting human genes" = hiv1_interacting
  ),
  label_alpha = 0,
  edge_size = 0.5,
  label = "count",
  label_size = 5,
  set_size = 4.5
) +
  ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
  ggplot2::theme(
    legend.position = "none",
    text = ggplot2::element_text(colour = "black")
  )

ggplot2::ggsave(
  file.path(output_dir, "Figure4A_HIV1_venn.pdf"),
  hiv1_venn_plot,
  width = 6,
  height = 5,
  device = grDevices::cairo_pdf
)

cat("Figure 4A data and statistical results were written to:\n", output_dir, "\n")
cat("Large-network hub genes:", paste(hub_genes, collapse = ", "), "\n")
cat("Fisher background size:", background_size, "genes\n")

