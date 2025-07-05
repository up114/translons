library(tidyverse)
library(ppcor)       
library(ggplot2)
library(scales)
library(readr)
library(dplyr)
library(cowplot)
library(viridis)
library(DescTools) 
library("ggExtra")

# Load files
orfE = read_csv("outputFTembryoORF_comp.csv", col_names = TRUE)
orfN = read_csv("outputFTneuralORF_comp.csv", col_names = TRUE)
dataE = read_csv("outputRNASeqORFembryo.csv", col_names = TRUE)
dataN = read_csv("outputRNASeqORFneural.csv", col_names = TRUE)

# Add Type
embryo_data <- left_join(dataE, dplyr::select(orfE, Gene, Start, Stop, Type), 
                         by = c("Gene", "Start", "Stop"))
neural_data <- left_join(dataN, dplyr::select(orfN, Gene, Start, Stop, Type), 
                         by = c("Gene", "Start", "Stop"))

# Label and combine datasets 
neural_data <- neural_data %>% mutate(Context = "Neural")
embryo_data <- embryo_data %>% mutate(Context = "Embryo")
combined_data <- bind_rows(neural_data, embryo_data)
reads_CPM <- data.frame(reads_CPM = ((combined_data$`ORF Reads` / combined_data$`Experiment Reads`) * 1000000))
combined_data <- cbind(combined_data, reads_CPM)
gene_CPM <- data.frame(gene_CPM = ((combined_data$`CDS Reads` / combined_data$`Experiment Reads`) * 1000000))
combined_data <- cbind(combined_data, gene_CPM)
combined_data <- combined_data %>% filter(!is.na(reads_CPM) & !is.na(gene_CPM))

# Combine and simplify cell lines
combined_data$Cell_Line[combined_data$Cell_Line %in% c("Liver cells", "liver tissue", "liver")] <- "Liver"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("ES cell line A3-1", "RF8 mouse embryonic stem cell line", "Embryonic stem cells", "embryonic stem cells", "A3-1", "v6.5", "R1 mouse embryonic stem cells", "ESC line CGR8", "E14", "E14Tg2a")] <- "ESC"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("brain neural stem cells", "neural tube")] <- "NSC"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("inguinal white fat", "subcutaneous white fat", "epididymal white fat", "interscapular brown fat")] <- "Fat"

combined_data$Cell_Line[combined_data$Cell_Line %in% c("Bone marrow derived primary dendritic cells", "Bone marrow derived dendritic cells", "Bone marrow derived regulatory dendritic cells", "Flt3L-DC")] <- "Bone marrow dendritic"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("3T3", "CD4 T-cells", "resting state T cells")] <- "T cells"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("Forebrain tissue (WT E14.5 mouse embryos)", "Forebrain tissue (Elp3cKO E14.5 mouse embryos)", "hippocampal", "cortical tissue", "brain cortex", "hippocampal cortex", "Cerebellum", "Fetal cortex", "dentate gyrus", "hemisphere")] <- "brain"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("Bone marrow derived macrophage", "RAW264", "peritoneal macrophage cells", "RAW264.7")] <- "Macrophage"
combined_data[combined_data == "RAW264.7"] <- "RAW264"

combined_data$Cell_Line[combined_data$Cell_Line %in% c("Mouse back skins", "primary epidermis", "Epidermal basal cells")] <- "Skin"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("Quadriceps muscle", "tibialis anterior tissue", "gastrocnemius tissue", "forelimbs")] <- "Skeletal muscle"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("ES cell derived neurons", "Primary cortical neurons", "Neurons (DIV 8) derived from CGR8 ES cells", "Cortical neuron", "Striatal cells")] <- "Neuron"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("EL4", "T-ALL")] <- "Leukemia"

combined_data$Cell_Line[combined_data$Cell_Line %in% c("Splenic B cells", "follicular B cells", "primary splenic B cells", "resting state B cells")] <- "B cells"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("1 cell", "2 cell", "4 cell", "morula", "blastocyst")] <- "Preimplantation embryo"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("duodenum", "Ileum")] <- "Small intestine"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("KH2", "P19")] <- "Tetracarcinoma"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("Spontaneously Immortalized Mouse Keratinocyte Cult", "Primary Keratinocytes")] <- "Keratinocyte"

combined_data$Cell_Line[combined_data$Cell_Line %in% c("Cardiac (left ventricle)", "heart")] <- "Heart"
combined_data$Cell_Line[combined_data$Cell_Line %in% c("testis", "testes")] <- "Testes"
combined_data$Cell_Line[combined_data$Cell_Line == "mouse eye"] <- "Eye"
combined_data$Cell_Line[combined_data$Cell_Line == "fat"] <- "Fat"
combined_data$Cell_Line[combined_data$Cell_Line == "brain"] <- "Brain"
combined_data$Cell_Line[combined_data$Cell_Line == "embryo"] <- "Embryo"
combined_data$Cell_Line[combined_data$Cell_Line == "lung"] <- "Lung"
combined_data$Cell_Line[combined_data$Cell_Line == "kidney"] <- "Kidney"
combined_data$Cell_Line[combined_data$Cell_Line == "liver"] <- "Liver"
combined_data$Cell_Line[combined_data$Cell_Line == "pancreas"] <- "Pancreas"
combined_data$Cell_Line[combined_data$Cell_Line == "skeletal muscle"] <- "Skeletal muscle"

#simplify names
combined_data[combined_data == "Dorsal section of lumbar spinal cord"] <- "Spinal cord"
combined_data[combined_data == "KRPC-A cells (from gen. eng. murine tumors)"] <- "KRPC-A"
combined_data[combined_data == "skin squamous tumours (skin papilloma)"] <- "Skin papilloma"
combined_data[combined_data == "embryonic fibroblast"] <- "MEF"
combined_data[combined_data == "lymphoid ba/f3 cells"] <- "Lymphoid"

# Capitalize
combined_data[combined_data == "neutrophils"] <- "Neutrophils"
combined_data[combined_data == "spleen"] <- "Spleen"

# Define ORF type
combined_data <- combined_data %>%
  mutate(ORF_Class = ifelse(str_detect(Type, "dORF"), "dORF", "uORF")) %>%
  mutate(Group = paste(Context, ORF_Class)) 

### Themes ###
theme_custom <- theme_minimal() +
  theme(
    panel.grid = element_blank(),                
    panel.border = element_blank(),              
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"), 
    axis.ticks = element_blank(),                
    axis.text = element_text(color = "black", size = 7), 
    axis.title = element_text(size = 8),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=7),
    legend.key.size = unit(0.5, "lines"),
  )

theme_custom_x0 <- theme_minimal() +
  theme(
    panel.grid = element_blank(),                 
    panel.border = element_blank(),               
    axis.line.x = element_blank(),                
    axis.line.y = element_line(color = "black"),  
    axis.ticks = element_blank(),                 
    axis.text = element_text(color = "black", size = 7), 
    axis.title = element_text(size = 7),
    legend.title=element_text(size=8), 
    legend.text=element_text(size=7),
    legend.key.size = unit(0.2, "lines"),
    legend.key.height = unit(0.2, "lines"),
    legend.key.width = unit(0.2, "lines"),
    legend.spacing.y = unit(0.1, "lines")
  ) 

### Dot Plot of Translation Level ###
# Collapse to one value per Gene/Start/Stop/Cell_Line
per_orf <- combined_data %>%
  group_by(Cell_Line, Group, Gene, Start, Stop) %>%
  summarise(reads = mean(reads_CPM, na.rm = TRUE), .groups = "drop")

plot_data <- per_orf %>%
  group_by(Cell_Line, Group) %>%
  summarise(
    Med_Reads = median(reads, na.rm = TRUE),
    Pct_Translated = mean(reads > 0, na.rm = TRUE) * 100,  # percent non-zero
    .groups = "drop"
  )


# Create a data frame for shading
background_rects <- data.frame(
  Group = factor(c("Neural dORF", "Neural uORF", "Embryo dORF", "Embryo uORF"),
                 levels = levels(factor(plot_data$Group))),
  fill_color = c("steelblue", "turquoise", "brown", "tomato")
)

plot_data <- plot_data %>%
  mutate(
    Dataset = case_when(
      str_detect(Group, "Neural") ~ "Neural",
      str_detect(Group, "Embryo") ~ "Embryo"
    ),
    ORF_Class = case_when(
      str_detect(Group, "uORF") ~ "uORF",
      str_detect(Group, "dORF") ~ "dORF"
    )
  )

neural_order <- plot_data %>%
  filter(Dataset == "Neural") %>%
  group_by(Cell_Line) %>%
  summarise(med = median(Med_Reads),
            pct_med = median(Pct_Translated)) %>%
  arrange(desc(med), desc(pct_med)) %>%
  pull(Cell_Line)  

plot_data$Cell_Line <- factor(plot_data$Cell_Line, levels = neural_order)

# Plot
p_dotplot <- ggplot(plot_data, aes(x = Cell_Line, y = Group)) +
  geom_rect(data = background_rects,
            aes(ymin = as.numeric(Group) - 0.5,
                ymax = as.numeric(Group) + 0.5),
            xmin = -Inf, xmax = Inf,
            fill = background_rects$fill_color,
            alpha = 0.05,
            inherit.aes = FALSE) +
  geom_point(aes(size = Pct_Translated, fill = Med_Reads+0.001), shape = 21,
             stroke = 0.2, 
             color = "black") +
  scale_fill_viridis(
    option = "plasma",
    name = "Median reads\n(CPM)",
    trans = "log10",
    labels = scales::label_number(accuracy = 0.001)
  ) +
  scale_size(range = c(1, 3), name = "Percent\ntranslated") +
  labs(x = NULL, y = NULL) +
  theme_custom +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.box.margin = margin(t = 50)
  )




### Embryo violin plot ###
plot_partial_correlation <- function(data, fillColor, orftype, xPos, yOff) {
  
  data <- data %>% filter(!is.na(RNASeq_Exp))
  data$Name = paste0(data$Gene,"_",data$Start)
  data$RNA_CPM = data$RNASeq_Gene/data$RNASeq_Exp * 10^6

  # Run Spearman partial correlation for each ORF
  unique_orfs = unique(data$Name)
  outcorr <- data.frame(Name = unique_orfs, R = NA, p = NA)

  for (orf in unique_orfs) {
    orf_data <- data %>% filter(Name == orf)

    if (sum(orf_data$gene_CPM) == 0) {
      outcorr[outcorr$Name == orf, "R"] <- NA
      outcorr[outcorr$Name == orf, "p"] <- NA
      next
    }

    test <- pcor.test(orf_data$reads_CPM, orf_data$gene_CPM, orf_data$RNA_CPM, method = "spearman")
    outcorr[outcorr$Name == orf, "R"] <- test$estimate
    outcorr[outcorr$Name == orf, "p"] <- test$p.value

  }
  outcorr$padj <- p.adjust(outcorr$p, method = "fdr")
  mean_R <- mean(outcorr$R, na.rm = TRUE)

  # Generate violin plot
 plot_data = outcorr %>% filter(!is.na(R))
 p =  ggplot(plot_data, aes(x = "", y = R)) +
    geom_violin(trim = TRUE, fill = fillColor) +
    geom_hline(yintercept = mean_R, linetype = "dashed", color = "black") +
   annotate("text", x = xPos, y = mean_R + yOff, label = bquote(bar(rho) == .(round(mean_R, 2))), color = "black", size = 2)+
   labs(x = "", y = expression("Spearman " * rho)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
   coord_cartesian(xlim = c(1, 1.3)) +
     theme_custom_x0
  return(p)
}

dataE = combined_data[combined_data$Context == "Embryo", ]
uorfE <- dataE %>% filter(!str_detect(Type, "dORF"))
dorfE <- dataE %>% filter(str_detect(Type, "dORF"))
p_uCorr = plot_partial_correlation(uorfE, "turquoise", "uORF", 1.69, 0.037)
p_dCorr = plot_partial_correlation(dorfE, "steelblue", "dORF", 1.61, 0.025)




### Neural specificity analysis + plot ###
specificity_correlation <- function(data, orftype, fillColor) {
  # cell line average
  avg <- data %>%
    group_by(Cell_Line, Gene, Start, Stop) %>%
    summarise(`ORF Reads` = mean(reads_CPM), .groups = "drop") %>%
    mutate(Name = paste0(Gene, "_", Start))
  
  avg_wide <- avg %>%
    dplyr::select(Name, Cell_Line, `ORF Reads`) %>%
    pivot_wider(names_from = Cell_Line, values_from = `ORF Reads`, values_fill = 0)
  
  # Compute Gini coefficient for each ORF across cell lines
  expr_cols <- setdiff(names(avg_wide), "Name")
  avg_wide$Gini <- apply(avg_wide[, expr_cols], 1, function(x) {
    if (sum(x) == 0) return(NA)
    Gini(x, na.rm = TRUE)
  })
  
  data <- data %>% filter(!is.na(RNASeq_Exp))
  data$Name = paste0(data$Gene,"_",data$Start)
  data$RNA_CPM = data$RNASeq_Gene/data$RNASeq_Exp * 10^6
  
  # Run Spearman partial correlation for each ORF
  unique_orfs = unique(data$Name)
  outcorr <- data.frame(Name = unique_orfs, R = NA, p = NA)
  
  for (orf in unique_orfs) {
    orf_data <- data %>% filter(Name == orf)
    
    if (sum(orf_data$gene_CPM) == 0) {
      outcorr[outcorr$Name == orf, "R"] <- NA
      outcorr[outcorr$Name == orf, "p"] <- NA
      next
    }
    
    test <- pcor.test(orf_data$reads_CPM, orf_data$gene_CPM, orf_data$RNA_CPM, method = "spearman")
    outcorr[outcorr$Name == orf, "R"] <- test$estimate
    outcorr[outcorr$Name == orf, "p"] <- test$p.value
  }
  
  # Spearman for Gini vs. ORF-CDS rho
  merged <- left_join(avg_wide[, c("Name", "Gini")], outcorr, by = "Name") %>%
    filter(!is.na(Gini) & !is.na(R))
  
  included_names <- merged$Name
  all_names <- intersect(avg_wide$Name, outcorr$Name)
  excluded_names <- setdiff(all_names, included_names)
  
  # Print excluded ORFs - usually due to no CDS reads in all datasets
  if (length(excluded_names) > 0) {
    message(length(excluded_names), " ORFs were excluded from the final correlation due to missing Gini or R.")
    print(excluded_names) 
  }
  
  cor_result <- cor.test(merged$Gini, merged$R, method = "spearman")
  print(cor_result)
  
  # Plot Gini vs Spearman rho
  merged$Group <- orftype
  p <- ggplot(merged, aes(x = Gini, y = R, color = Group)) +
    geom_point(alpha = 0.9, size = 1) +
    scale_color_manual(values = setNames(fillColor, orftype)) +
    labs(
      x = "Gini Coefficient",
      y = expression("Partial Spearman " * rho)
    ) +
    annotate(
      "text",
      x = 0.05, y = max(merged$R), 
      label = bquote(rho == .(round(cor_result$estimate, 2))),
      hjust = 0, size = 2.5
    ) +
    theme_custom_x0 +
    theme(
      legend.position = "none"
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) 
  p = ggMarginal(p, type = "density", groupFill = TRUE, margins = "y")
  return(list(plot = p, corr_df = merged, corr_res = cor_result))
}

dataN = combined_data[combined_data$Context == "Neural", ]
uorfN <- dataN %>% filter(!str_detect(Type, "dORF"))
dorfN <- dataN %>% filter(str_detect(Type, "dORF"))
outputU = specificity_correlation(uorfN, "uORF", "tomato")
outputD = specificity_correlation(dorfN, "dORF", "brown")

### Arrangement ###
row2 = plot_grid(p_uCorr, p_dCorr, outputU$plot, outputD$plot, "", labels = c("B", "C", "D", "E", ""), label_size = 10, ncol = 5, rel_widths = c(1,1,1,1,0.05))
fig4 = plot_grid(p_dotplot, row2, ncol = 1, labels = c("A", ""), label_size = 10, rel_heights = c(1.5, 1))
ggsave("fig4.pdf", fig4, width = 6.3, height = 4.3,  units = "in")

