### Complete plotting code 
library(tidyverse)
library(ppcor)       
library(ggplot2)
library(scales)
library(readr)
library(dplyr)
library(cowplot)

### Files ###
orfE = read_csv("outputFTembryoORF_comp.csv", col_names = TRUE)
orfN = read_csv("outputFTneuralORF_comp.csv", col_names = TRUE)
orfE <- orfE %>%
  mutate(Source = case_when(
    str_detect(Type, "dORF") ~ "dORF",
    TRUE ~ "uORF"
  ))
orfN <- orfN %>%
  mutate(Source = case_when(
    str_detect(Type, "dORF") ~ "dORF",
    TRUE ~ "uORF"
  ))

#### Set Theme ####
theme_custom <- theme_minimal() +
  theme(
    panel.grid = element_blank(),                # no grid lines
    panel.border = element_blank(),              # no border box
    axis.line.x = element_line(color = "black"), # black x-axis line
    axis.line.y = element_line(color = "black"), # black y-axis line
    axis.ticks = element_blank(),                # remove ticks
    axis.text = element_text(color = "black", size = 7),  # combine here
    axis.title = element_text(size = 8),
    legend.title=element_text(size=8), 
    legend.text=element_text(size=7),
    legend.key.size = unit(0.2, "lines"),
    legend.key.height = unit(0.2, "lines"),
    legend.key.width = unit(0.2, "lines"),
    legend.spacing.y = unit(0.1, "lines")
  )

theme_custom_x0 <- theme_minimal() +
  theme(
    panel.grid = element_blank(),                 # no grid lines
    panel.border = element_blank(),               # no border box
    axis.line.x = element_blank(),                # REMOVE bottom x-axis line
    axis.line.y = element_line(color = "black"),  # KEEP solid black y-axis
    axis.ticks = element_blank(),                 # no ticks
    axis.text = element_text(color = "black", size = 7),  # combine here
    axis.title = element_text(size = 8),
    legend.title=element_text(size=8), 
    legend.text=element_text(size=7),
    legend.key.size = unit(0.2, "lines"),
    legend.key.height = unit(0.2, "lines"),
    legend.key.width = unit(0.2, "lines"),
    legend.spacing.y = unit(0.1, "lines")
  ) 

### Characterization ###
# Start Codon

plot_codon <- function(orfs, uorf_fill, dorf_fill) {
  # Group by Source and Start Codon
  codons <- orfs %>%
    group_by(Source, Start_Codon) %>%
    summarise(Freq = n(), .groups = 'drop')
  
  # Set preferred order of start codons
  codon_order <- c("ATG", "CTG", "ACG", "GTG", "ATT", "TTG")
  source_levels <- c("uORF", "dORF")
  codons <- codons %>%
    complete(Source = source_levels, Start_Codon = codon_order, fill = list(Freq = 0)) %>%
    mutate(Start_Codon = factor(Start_Codon, levels = codon_order))
  
  # Plot start codon frequencies
  p_codon = ggplot(codons, aes(x = Start_Codon, y = Freq, fill = Source)) +
    geom_bar(stat = 'identity', position = "dodge") +
    labs(x = "Start Codon", y = "Count", fill = NULL) +
    theme_custom+
    guides(
      fill = guide_legend(title.theme = element_text(size = 8),
                          label.theme = element_text(size = 7),
                          override.aes = list(size = 3)),
      alpha = guide_legend(title.theme = element_text(size = 8),
                           label.theme = element_text(size = 7),
                           override.aes = list(size = 3))
    ) +
    theme(legend.position = "none") +
    scale_fill_manual(values = c(dorf_fill, uorf_fill))
  
  return(p_codon)
}

codonE = plot_codon(orfE, "turquoise", "steelblue")
codonN = plot_codon(orfN, "tomato", "brown")

### AA Length ###
plot_length <- function(orfs, uorf_fill, dorf_fill) {
  orfs$AALength = ((orfs$True_Stop - orfs$Start)/3)
  orfs$Length_Capped <- ifelse(orfs$AALength >= 200, 200, orfs$AALength)
  orfs$Source <- factor(orfs$Source, levels = c("uORF", "dORF"))
  p_len = ggplot(orfs, aes(x = Length_Capped, fill = Source, color = Source)) +
    geom_density(adjust = 0.9, aes(y = after_stat(count)*5),
                 linewidth = 0.6, alpha = 0.7) +
    labs(
      x = "Amino Acid Length",
      y = "Count",
      fill = NULL,
      color = NULL
    ) +
    theme_custom +
    guides(
      fill = guide_legend(title.theme = element_text(size = 8),
                          label.theme = element_text(size = 7),
                          override.aes = list(size = 3)),
      alpha = guide_legend(title.theme = element_text(size = 8),
                           label.theme = element_text(size = 7),
                           override.aes = list(size = 3))
    ) +
    theme(
      legend.position = c(0.98, 0.98),  # near top-right inside plot area
      legend.justification = c("right", "top")
    ) +
    scale_fill_manual(values = c(uorf_fill, dorf_fill)) +
    scale_color_manual(values = c(uorf_fill, dorf_fill)) +
    scale_x_continuous(breaks = c(seq(0, 200, by = 20)), labels = c(seq(0, 180, by = 20), "200+"))
  return(p_len)
}

lenE = plot_length(orfE, "turquoise", "steelblue")
lenN = plot_length(orfN, "tomato", "brown")

### AA Composition ###
get_aa_freqs_dist <- function(df, orf_label) {
  df <- df %>% mutate(Name = paste0(Gene, "_", Start))
  
  # CDS AA frequencies
  cds_freqs <- df %>%
    mutate(Split = strsplit(CDS_AASequence, "")) %>%
    unnest(Split) %>%
    filter(Split != "*") %>%
    group_by(Name, Split) %>%
    summarise(CDS_Count = n(), .groups = "drop") %>%
    group_by(Name) %>%
    mutate(CDS_Total = sum(CDS_Count)) %>%
    ungroup() %>%
    mutate(Freq = CDS_Count / CDS_Total,
           Group = "CDS") %>%
    dplyr::select(Name, Amino_Acid = Split, Freq, Group)
  
  # ORF AA frequencies
  orf_freqs <- df %>%
    mutate(Split = strsplit(AASequence, "")) %>%
    unnest(Split) %>%
    filter(Split != "*") %>%
    group_by(Name, Split) %>%
    summarise(ORF_Count = n(), .groups = "drop") %>%
    group_by(Name) %>%
    mutate(ORF_Total = sum(ORF_Count)) %>%
    ungroup() %>%
    mutate(Freq = ORF_Count / ORF_Total,
           Group = orf_label) %>%
    dplyr::select(Name, Amino_Acid = Split, Freq, Group)
  
  # Merge — ensure all CDS AAs are present in ORF (fill ORF missing as 0)
  all_freqs <- full_join(cds_freqs, orf_freqs, by = c("Name", "Amino_Acid")) %>%
    mutate(
      Freq.x = replace_na(Freq.x, 0),  # CDS freq
      Freq.y = replace_na(Freq.y, 0)   # ORF freq
    ) %>%
    pivot_longer(
      cols = c(Freq.x, Freq.y),
      names_to = "source",
      values_to = "Freq"
    ) %>%
    mutate(Group = ifelse(source == "Freq.x", "CDS", "ORF")) %>%
    dplyr::select(Name, Amino_Acid, Freq, Group)
  
  return(all_freqs)
}

plot_AAcomp <- function(orfs) {
  uorf_data <- orfs %>% filter(!str_detect(Type, "dORF"))
  dorf_data <- orfs %>% filter(str_detect(Type, "dORF"))
  aa_freqs_u <- get_aa_freqs_dist(uorf_data, "uORF")
  aa_freqs_d <- get_aa_freqs_dist(dorf_data, "dORF")
  
  wilcox_u <- aa_freqs_u %>%
    group_by(Amino_Acid) %>%
    summarise(
      p_value = tryCatch(
        wilcox.test(Freq ~ Group)$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    mutate(ORF_Type = "uORF")
  
  wilcox_d <- aa_freqs_d %>%
    group_by(Amino_Acid) %>%
    summarise(
      p_value = tryCatch(
        wilcox.test(Freq ~ Group)$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    mutate(ORF_Type = "dORF")
  
  # Combine Wilcoxon p-values
  wilcox_all <- bind_rows(wilcox_u, wilcox_d) %>%
    mutate(p_adj = p.adjust(p_value, method = "bonferroni"))
  
  # Compute medians
  medians_u <- aa_freqs_u %>%
    group_by(Amino_Acid, Group) %>%
    summarise(Median = median(Freq, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = Median) %>%
    mutate(
      Delta = (ORF - CDS)/CDS,
      Dataset = "uORF"
    )
  medians_d <- aa_freqs_d %>%
    group_by(Amino_Acid, Group) %>%
    summarise(Median = median(Freq, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = Median) %>%
    mutate(
      Delta = (ORF - CDS)/CDS,
      Dataset = "dORF"
    )
  
  # Combine
  medians_all <- bind_rows(medians_u, medians_d)
  
  diff_summary <- medians_all %>%
    left_join(wilcox_all, by = c("Amino_Acid", "Dataset" = "ORF_Type")) %>%
    mutate(Signif = ifelse(!is.na(p_adj) & p_adj < 0.05, "*", ""))
  
  
  aa_order <- c("D", "E", "A", "G", "I", "L", "M", "V", "C", "N", "P", "Q", "S", "T", "F", "W", "Y", "H", "K", "R")
  diff_summary$Amino_Acid <- factor(diff_summary$Amino_Acid, levels = aa_order)
  
  diff_summary$Type <- NA
  diff_summary$Type[diff_summary$Amino_Acid %in% c("A", "G", "I", "L", "M", "V")] <- "Nonpolar,\naliphatic"
  diff_summary$Type[diff_summary$Amino_Acid %in% c("C", "N", "P", "Q", "S", "T")] <- "Polar,\nneutral"
  diff_summary$Type[diff_summary$Amino_Acid %in% c("F", "W", "Y")] <- "Nonpolar,\naromatic"
  diff_summary$Type[diff_summary$Amino_Acid %in% c("H", "K", "R")] <- "Positive"
  diff_summary$Type[diff_summary$Amino_Acid %in% c("D", "E")] <- "Negative"
  
  # Plot with no legend
  p_1 <- ggplot(diff_summary, aes(x = Amino_Acid, y = Delta*100, fill = Type, alpha = Dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 1)) +
    geom_text(aes(label = Signif, vjust = ifelse(Delta > 0, 0.3, 1.2)),
              position = position_dodge(width = 1),
              size = 4) +
    scale_fill_manual(
      values = c(
        "Nonpolar,\naliphatic" = "#1f77b4",
        "Polar,\nneutral" = "#2ca02c",
        "Nonpolar,\naromatic" = "#ff7f0e",
        "Positive" = "#d62728",
        "Negative" = "#9467bd"
      ),
      breaks = c(
        "Negative",
        "Nonpolar,\naliphatic",
        "Nonpolar,\naromatic",
        "Polar,\nneutral",
        "Positive"
      )
    ) +
    scale_alpha_manual(values = c("uORF" = 1, "dORF" = 0.5)) +
    labs(
      x = "Amino Acid",
      y = "ORF-CDS % change in median",
      fill = "Residue",
      alpha = "Type"
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    theme_custom_x0 +
    theme(
      legend.position = "none"
    ) 
  
  # Plot with legend
  p_2 <- ggplot(diff_summary, aes(x = Amino_Acid, y = Delta*100, fill = Type, alpha = Dataset)) +
    geom_bar(stat = "identity", position = position_dodge(width = 1)) +
    geom_text(aes(label = Signif, vjust = ifelse(Delta > 0, 0.3, 1.2)),
              position = position_dodge(width = 1),
              size = 4) +
    scale_fill_manual(
      values = c(
        "Nonpolar,\naliphatic" = "#1f77b4",
        "Polar,\nneutral" = "#2ca02c",
        "Nonpolar,\naromatic" = "#ff7f0e",
        "Positive" = "#d62728",
        "Negative" = "#9467bd"
      ),
      breaks = c(
        "Negative",
        "Nonpolar,\naliphatic",
        "Nonpolar,\naromatic",
        "Polar,\nneutral",
        "Positive"
      )
    ) +
    scale_alpha_manual(values = c("uORF" = 1, "dORF" = 0.5)) +
    labs(
      x = "Amino Acid",
      y = "ORF-CDS % change in median",
      fill = "Residue",
      alpha = "Type"
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    theme_custom_x0 +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.spacing.y = unit(0.15, "lines"),
      legend.text = element_text(size = 5, margin = margin(l = 1)),  
      legend.key.size = unit(0.15, "lines")
    ) +
    guides(
      fill = "none",
      alpha = guide_legend(
        title.theme = element_text(size = 7),
        label.theme = element_text(size = 5.5),
        override.aes = list(shape = 15, size = 4, fill = c("black", "black"), color = NA),  # square boxes, no stroke
        byrow = TRUE
      )
    )
  return(list(p1 = p_1, p2 = p_2, corr_res = diff_summary))
}

aacompE = plot_AAcomp(orfE)$p2
aacompN = plot_AAcomp(orfN)$p1


# # #### Arrangement ####
# # Fig 3: start codon, length, AA comp for neural and embryo
rowN = plot_grid(lenN, codonN, labels = c("A", "B"), label_size = 10, ncol=1)
figN = plot_grid(rowN, aacompN, labels = c("", "C"), label_size = 10, rel_widths = c(1, 1.3))
rowE = plot_grid(lenE, codonE, labels = c("D", "E"), label_size = 10, ncol=1)
figE = plot_grid(rowE, aacompE, labels = c("", "F"), label_size = 10, rel_widths = c(1, 1.3))
fig3 = plot_grid(figN, figE, labels = c("", ""), label_size = 10, ncol = 1)
ggsave("fig3FIN.pdf", fig3, width = 6.3, height = 5.5, units = "in")


### RiboTish AA Comp ###
orfRT = read_csv("neuralRiboTish_comp.csv")

orfRT$Type <- ifelse(orfRT$TisType == "5'UTR", "uORF", "dORF")
orfRT = orfRT %>% rename(
  Gene = Symbol,
  AASequence = AASeq,
  CDS_AASequence = CDS.AA.Sequence
)

aacompRT = plot_AAcomp(orfRT)$p2
rowRT = plot_grid("", aacompRT)
fig3b = plot_grid(fig3, rowRT, ncol = 1, rel_heights = c(2,1))
ggsave("fig3RT.pdf", fig3b, width = 6.3, height = 8.25, units = "in")
