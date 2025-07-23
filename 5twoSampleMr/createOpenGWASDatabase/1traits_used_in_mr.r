library(TwoSampleMR)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)

# --- Setup ---
outputPath <- './prj_190_viking_somalogic/Scripts/twoSampleMr/startingDataset.tsv'

# --- Step 1: Load available outcomes ---
allOutcomes <- available_outcomes()

# Filter only relevant outcome sets
filteredOutcomes <- allOutcomes %>%
  filter(grepl("ebi", id) | grepl("ieu", id) | grepl("ukb-b", id)) %>%
  rename(id.outcome = id) %>%
  select(id.outcome, trait, population, sample_size, ncase, ncontrol, 
         subcategory, category, ontology, consortium, author, year, pmid)

# --- Step 2: Normalize trait names for harmonization ---
harmonizedOutcomes <- filteredOutcomes %>%
  mutate(outcome_name = tolower(trait)) %>%
  mutate(outcome_name = gsub("high density lipoprotein cholesterol levels", "hdl cholesterol", outcome_name)) %>%
  mutate(outcome_name = gsub("low density lipoprotein cholesterol levels", "ldl cholesterol", outcome_name)) %>%
  mutate(outcome_name = gsub("total cholesterol levels|cholesterol, total", "total cholesterol", outcome_name)) %>%
  mutate(outcome_name = gsub("triglyceride levels", "triglycerides", outcome_name)) %>%
  mutate(outcome_name = gsub("c-reactive protein levels?|c-reactive protein level", "c-reactive protein", outcome_name)) %>%
  mutate(outcome_name = gsub("body fat percentage", "body fat", outcome_name)) %>%
  mutate(outcome_name = gsub("eosinophil counts", "eosinophil cell count", outcome_name)) %>%
  mutate(outcome_name = gsub("monocyte count", "monocyte cell count", outcome_name)) %>%
  mutate(outcome_name = gsub("neutrophil count", "neutrophil cell count", outcome_name)) %>%
  mutate(outcome_name = gsub("neurociticism", "neuroticism", outcome_name)) %>%
  mutate(outcome_name = gsub("peak expiratory flow \\(pef\\)", "peak expiratory flow", outcome_name)) %>%
  mutate(outcome_name = gsub("urate levels", "urate", outcome_name))

# --- Step 3: Deduplicate outcomes by outcome_name ---
# Keep preference: European > Mixed > largest ncase/sample_size
deduplicatedOutcomes <- ddply(harmonizedOutcomes, .(outcome_name), function(group) {
  message("Processing outcome: ", unique(group$outcome_name))

  if (any(group$population == "European", na.rm = TRUE)) {
    group <- group %>% filter(population == "European")
  } else if (any(group$population == "Mixed", na.rm = TRUE)) {
    group <- group %>% filter(population == "Mixed")
  }

  if (any(!is.na(group$ncase))) {
    group %>% slice(which.max(ncase))
  } else {
    group %>% slice(which.max(sample_size))
  }
})

# --- Step 4: Export cleaned list ---
selectedOutcomes <- deduplicatedOutcomes %>%
  select(id.outcome, outcome_name, population, sample_size, 
         ncase, ncontrol, consortium, author, year, pmid) %>%
  arrange(outcome_name)

fwrite(selectedOutcomes, file = outputPath, sep = "\t", quote = FALSE, row.names = FALSE)
