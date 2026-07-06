#Extraction from IMPC database

setwd("/dir/proteomewide/supplement")

#Extract the markers with evidence from the supplement file

all <- read.csv("S_ALL_MR_FLAGS.csv")

genes<- all[which(all$Evidence_Tier=="Tier_A"|all$Evidence_Tier=="Tier_B"|all$Evidence_Tier=="Tier_C"),]

library(homologene)
library(httr)
library(jsonlite)
library(dplyr)

# Map human-mouse orthologs

orthologs <- homologene(genes = genes$HGNC_gene, inTax = 9606, outTax = 10090)


human_mouse <- orthologs %>%
  rename(
    human_symbol = `9606`,
    mouse_symbol = `10090`,
    human_entrez = `9606_ID`,
    mouse_entrez = `10090_ID`
  ) %>%
  distinct()


# Query IMPC Solr core

query_impc_core <- function(core = "genotype-phenotype", params) {
  base_url <- "https://www.ebi.ac.uk/mi/impc/solr"
  full_url <- paste0(base_url, "/", core, "/select")
  
  resp <- GET(full_url, query = params)
  stop_for_status(resp)
  
  # parse JSON
  content(resp, as = "text", encoding = "UTF-8") %>%
    fromJSON(flatten = TRUE)
}

# Fetch phenotypes for one gene- this is useful for subsequent full gene list function
get_gene_phenotypes <- function(gene_symbol, rows = 20000) {
  params <- list(
    q = paste0("marker_symbol:", gene_symbol),
    wt = "json",
    fl = paste(
      "marker_symbol",
      "mp_term_name",
      "mp_id",
      "p_value",
      "zygosity",
      "sex",
      sep = ","
    ),
    rows = rows
  )
  
  solr_resp <- query_impc_core("genotype-phenotype", params)
  
  docs <- solr_resp$response$docs
  if (length(docs) == 0) return(data.frame())
  
  as_tibble(docs)
}

# Fetch for my list of genes
genes <- human_mouse$mouse_symbol 

all_results <- lapply(genes, function(g) {
  df <- get_gene_phenotypes(g)
  if (nrow(df) > 0) df$queried_gene <- g
  df
})

all_df <- bind_rows(all_results)


print(head(all_df))

#Merge with the human mouse df to retain entrez ids

human_mouse <- human_mouse %>%
  rename(marker_symbol = mouse_symbol)

merged_df <- all_df %>%
  left_join(human_mouse, by = "marker_symbol") %>%
  select(
    queried_gene,
    marker_symbol,
    mp_term_name,
    p_value,
    zygosity,
    sex,
    human_symbol,
    human_entrez,
    mouse_entrez
  )

write.csv(merged_df, "IMPC_annotations.csv", row.names = FALSE)
