library(tidyverse)

### Cell Line Average Table ###
# organize  data
orf_info <- read_csv("outputFTembryoORF_comp.csv", col_names = TRUE)
ribobase <- read_csv("outputRNASeqORFembryo.csv", col_names = TRUE)
ribobase <- left_join(ribobase, dplyr::select(orf_info, Gene, Start, Stop, True_Stop), 
                    by = c("Gene", "Start", "Stop"))
ribobase$CPM = ribobase$`uORF Reads`/ribobase$`Experiment Reads` * 1000000

ribobase <- ribobase %>%
  group_by(Cell_Line, Gene, Start, True_Stop) %>%
  summarise(`ORF Reads` = mean(CPM))

ribobase <- pivot_wider(ribobase, names_from = Cell_Line, values_from = 'ORF Reads') %>% 
  dplyr::rename(
    Stop = True_Stop
  )

write.csv(ribobase, "embryoCelllineAvg.csv", row.names = FALSE)



### ORF Summary Table ###
preDf = read_csv("candidateORFs.csv", col_names = TRUE)
postDf = read_csv("outputFT.csv", col_names = TRUE)
# Select relevant columns from preDf
preDf <- preDf %>%
  dplyr::select(Transcript, Gene, True_Stop, Start, Start_Codon, Type)

# Add "FT Accept" column: 1 if match found in postDf by Gene, Start, True_Stop
preDf <- preDf %>%
  mutate(`FT Accept` = ifelse(paste(Gene, Start, True_Stop) %in%
                                paste(postDf$Gene, postDf$Start, postDf$True_Stop), 1, 0)) %>%
  dplyr::rename(Stop = True_Stop)

write.csv(preDf, "embryoORFs.csv", row.names = FALSE)
