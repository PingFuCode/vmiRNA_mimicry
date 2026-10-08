# Figure 3D: expression of virus-like hmiRNAs and human-mimicry vmiRNAs

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

data_dir <- "D:/0work/0wholetransctiptome/2sRNAminic/13Infect_mimic_correlation/Time"
output_dir <- "D:/0work/paper/0wtite/miRNA/V7_submit/code"
target_map_file <- "D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/0006mer_target_with_hsaID.csv"

# Retain the original data processing and RPM > 1 filter.
hiv_data <- read.csv(file.path(data_dir, "HIV-1_mimic_pairs_data.csv")) %>%
  mutate(source = "HIV-1") %>%
  select(time_point, vRPM, hRPM, source, pair_id)

hhv_data <- read.csv(file.path(data_dir, "HHV-5_mimic_pairs_data.csv")) %>%
  mutate(source = "HHV-5") %>%
  select(time_point, vRPM, hRPM, source, pair_id)

sinv_data <- read.csv(file.path(data_dir, "SINV_mimic_pairs_data.csv")) %>%
  mutate(source = "SINV") %>%
  select(time_point, vRPM, hRPM, source, pair_id)

combined_long <- bind_rows(hiv_data, hhv_data, sinv_data) %>%
  pivot_longer(
    cols = c(vRPM, hRPM),
    names_to = "RNA_type",
    values_to = "RPM"
  ) %>%
  filter(RPM > 1) %>%
  mutate(
    RNA_type = factor(
      RNA_type,
      levels = c("hRPM", "vRPM"),
      labels = c("virus-like hmiRNAs", "human-mimicry vmiRNAs")
    )
  )

# Map the viral sequence identifier in pair_id to the standard miRNAID.
virus_id_map <- read.csv(target_map_file, check.names = FALSE) %>%
  select(Virus_ID, miRNAID, miRNA) %>%
  distinct()

# Export the exact filtered plotting data in one row per miRNA pair.
figure3d_source_data <- combined_long %>%
  separate(
    pair_id,
    into = c("Virus_ID", "virus-like hmiRNAs"),
    sep = "_vs_",
    extra = "merge",
    fill = "right"
  ) %>%
  left_join(
    virus_id_map,
    by = c("Virus_ID", "virus-like hmiRNAs" = "miRNA")
  ) %>%
  mutate(
    source = factor(source, levels = c("HIV-1", "HHV-5", "SINV")),
    time_order = as.numeric(gsub("hpi", "", time_point)),
    RPM_type = recode(
      as.character(RNA_type),
      "virus-like hmiRNAs" = "virus-like hmiRNA RPM",
      "human-mimicry vmiRNAs" = "human-mimicry vmiRNA RPM"
    )
  ) %>%
  select(
    source,
    time_point,
    `virus-like hmiRNAs`,
    `human-mimicry vmiRNAs` = miRNAID,
    RPM_type,
    RPM,
    time_order
  ) %>%
  pivot_wider(
    names_from = RPM_type,
    values_from = RPM
  ) %>%
  mutate(
    `virus-like hmiRNA RPM` = ifelse(
      is.na(`virus-like hmiRNA RPM`),
      NA_character_,
      sprintf("%.3f", `virus-like hmiRNA RPM`)
    ),
    `human-mimicry vmiRNA RPM` = ifelse(
      is.na(`human-mimicry vmiRNA RPM`),
      NA_character_,
      sprintf("%.3f", `human-mimicry vmiRNA RPM`)
    )
  ) %>%
  arrange(
    source,
    time_order,
    `virus-like hmiRNAs`,
    `human-mimicry vmiRNAs`
  ) %>%
  select(
    source,
    time_point,
    `virus-like hmiRNAs`,
    `human-mimicry vmiRNAs`,
    `virus-like hmiRNA RPM`,
    `human-mimicry vmiRNA RPM`
  )

stopifnot(!any(is.na(figure3d_source_data$`human-mimicry vmiRNAs`)))

write.csv(
  figure3d_source_data,
  file.path(output_dir, "Figrue 3d.csv"),
  row.names = FALSE
)

create_boxplot <- function(data_source) {
  plot_data <- combined_long %>%
    filter(source == data_source)

  time_order <- unique(plot_data$time_point)
  time_order <- time_order[
    order(as.numeric(gsub("hpi", "", time_order)))
  ]
  plot_data$time_point <- factor(plot_data$time_point, levels = time_order)

  ggplot(
    plot_data,
    aes(x = time_point, y = RPM, fill = RNA_type, group = interaction(time_point, RNA_type))
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7,
      width = 0.65,
      color = "black",
      linewidth = 0.15,
      position = position_dodge(width = 0.8)
    ) +
    geom_point(
      color = "black",
      shape = 16,
      size = 1.3,
      alpha = 0.65,
      position = position_jitterdodge(
        jitter.width = 0.06,
        jitter.height = 0,
        dodge.width = 0.8,
        seed = 123
      )
    ) +
    facet_wrap(~source) +
    scale_fill_manual(
      values = c(
        "virus-like hmiRNAs" = "#4E79A7",
        "human-mimicry vmiRNAs" = "#BEBEBE"
      ),
      breaks = c("virus-like hmiRNAs", "human-mimicry vmiRNAs")
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    labs(
      x = "",
      y = "RPM",
      fill = NULL
    ) +
    theme_classic() +
    theme(
      text = element_text(color = "black"),
      strip.background = element_rect(
        color = "black",
        fill = "white",
        linewidth = 0.15
      ),
      strip.text = element_text(
        color = "black",
        face = "bold",
        size = 10,
        margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
      ),
      axis.title = element_text(color = "black", size = 10),
      axis.text = element_text(color = "black", size = 10),
      axis.line = element_line(color = "black", linewidth = 0.15),
      axis.ticks = element_line(color = "black", linewidth = 0.15),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.15
      ),
      legend.text = element_text(color = "black", size = 10),
      legend.position = "bottom",
      plot.margin = margin(5, 5, 5, 5, "pt")
    )
}

p_hiv1 <- create_boxplot("HIV-1")
p_hhv5 <- create_boxplot("HHV-5")
p_sinv <- create_boxplot("SINV")

figure_3d <- (p_hiv1 + p_hhv5 + p_sinv) +
  plot_layout(widths = c(2, 3, 2), guides = "collect") &
  theme(legend.position = "bottom")

print(figure_3d)

ggsave(
  file.path(output_dir, "Figure3D_boxplot.pdf"),
  plot = figure_3d,
  width = 7,
  height = 2.7,
  units = "in"
)
