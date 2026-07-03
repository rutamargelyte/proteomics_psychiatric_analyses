#Tractability evidence extraction for the biomarkers

setwd("/dir/proteomewide/supplement")

#Extract the biomarkers with evidence from the supplement file

genes<- read.csv("S_ALL_MR_FLAGS.csv")
genes<- genes[which(genes$Evidence_Tier=="Tier_A"|genes$Evidence_Tier=="Tier_B"|genes$Evidence_Tier=="Tier_C"),]
focus<- genes$ENSEMBL_ID 
focus <- unique(genes$ENSEMBL_ID)


###

library(httr)
library(jsonlite)
library(dplyr)

#Create function to extract tractability data

get_tractability_group <- function(ensembl_id) {
  
  query <- paste0('
  {
    target(ensemblId: "', ensembl_id, '") {
      approvedSymbol
      tractability {
        label
        modality
        value
      }
    }
  }
  ')
  
  res <- POST(
    "https://api.platform.opentargets.org/api/v4/graphql",
    body = list(query = query),
    encode = "json"
  )
  
  parsed <- fromJSON(content(res, "text"))
  
  tract <- parsed$data$target$tractability
  symbol <- parsed$data$target$approvedSymbol
  

  if (is.null(tract) || nrow(tract) == 0) {
    return(data.frame(
      ensembl_id = ensembl_id,
      gene_symbol = symbol,
      group = NA
    ))
  }
  
  # Define groups
  group1_labels <- c("Approved Drug", "Advanced Clinical", "Phase 1 Clinical")
  
  group2_labels_sm <- c(
    "Structure with Ligand",
    "High-Quality Ligand",
    "High-Quality Pocket",
    "Med-Quality Pocket",
    "Druggable Family"
  )
  
  group2_labels_ab <- c(
    "UniProt loc high conf",
    "GO CC high conf"
  )
  
  # Assign group 
 
  group <- case_when(
    any(tract$label %in% group1_labels & tract$value) ~ "Group 1: Strong druggability",
    
    any(
      (tract$modality == "SM" & tract$label %in% group2_labels_sm & tract$value) |
      (tract$modality == "AB" & tract$label %in% group2_labels_ab & tract$value)
    ) ~ "Group 2: Likely druggable",
    
    TRUE ~ "Group 3: Low/unknown druggability"
  )
  
  return(data.frame(
    ensembl_id = ensembl_id,
    gene_symbol = symbol,
    group = group
  ))
}

results_list <- lapply(focus, get_tractability_group)
results_df <- bind_rows(results_list)

write.csv(results_df, "open_targets_tractability_results.csv", row.names = FALSE)

#Extract the trials and associated diseases: 

library(dplyr)
library(purrr)
library(otargen)

# Here only the ones with strong druggability are relevant
strong_genes <- results_df %>%
  filter(group == "Group 1: Strong druggability")

get_drugs_otargen <- function(ensembl_id) {
  
  res <- knownDrugsGeneQuery(ensgId = ensembl_id)
  
  # Always return a dataframe
  if (is.null(res) || nrow(res) == 0) {
    return(tibble())
  }
  
  res %>%
    mutate(ensembl_id = ensembl_id)
}

# Run the function across all with strong druggability
drugs_df <- map_dfr(strong_genes$ensembl_id, get_drugs_otargen)

#Tidy up and save

library(dplyr)
library(purrr)
library(tidyr)

drugs_df_flat <- drugs_df %>%
  mutate(
    diseases = map(diseases, ~ {
      if (is.null(.x)) return(NULL)
      .x
    }),
    clinicalReports = map(clinicalReports, ~ {
      if (is.null(.x)) return(NULL)
      .x
    })
  )
  
diseases_df <- drugs_df %>%
  select(ensembl_id, drug.name, diseases) %>%
  unnest(diseases)

  trials_df <- drugs_df %>%
  select(ensembl_id, drug.name, clinicalReports) %>%
  unnest(clinicalReports)
  
write.csv(diseases_df, "diseases_clean.csv", row.names = FALSE)
write.csv(trials_df, "trials_clean.csv", row.names = FALSE)  