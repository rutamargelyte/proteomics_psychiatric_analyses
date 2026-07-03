#Extraction of GO and OMIM terms for the biomarkers

setwd("/dir/proteomewide/supplement")

#Extract the biomarkers with evidence

all <- read.csv("S_ALL_MR_FLAGS.csv")

genes<- all[which(all$Evidence_Tier=="Tier_A"|all$Evidence_Tier=="Tier_B"|all$Evidence_Tier=="Tier_C"),]

  
library(biomaRt)
library(dplyr)

# Connect to Ensembl
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

# Initialize lists to store results
omim_list <- list()
go_list <- list()

# Loop over genes
for (gene in genes$HGNC_gene) {
  
  # OMIM
  omim <- getBM(
    attributes = c("hgnc_symbol", "mim_morbid_accession", "mim_morbid_description"),
    filters = "hgnc_symbol",
    values = gene,
    mart = ensembl
  )
  
  if(nrow(omim) == 0) {
    omim <- data.frame(hgnc_symbol = gene, mim_morbid_accession = NA, mim_morbid_description = NA)
  }
  
  omim_list[[gene]] <- omim
  
  # GO
  go <- getBM(
    attributes = c("hgnc_symbol", "go_id", "name_1006", "namespace_1003"),
    filters = "hgnc_symbol",
    values = gene,
    mart = ensembl
  )
  
  if(nrow(go) == 0) {
    go <- data.frame(hgnc_symbol = gene, go_id = NA, name_1006 = NA, namespace_1003 = NA)
  }
  
  go_list[[gene]] <- go
}

# Combine all genes for OMIM
omim_df <- do.call(rbind, omim_list)
colnames(omim_df) <- c("Gene", "OMIM_ID", "OMIM_Name")
omim_df <- unique(omim_df)

# Combine all genes for GO
go_df <- do.call(rbind, go_list)
colnames(go_df) <- c("Gene", "GO_ID", "GO_Term", "GO_Category")
go_df <- unique(go_df)

# Preview
print("OMIM:")
print(omim_df)
print("GO:")
print(go_df)

# Save to separate CSV files
write.csv(omim_df, "Gene_OMIM_annotations.csv", row.names = FALSE)
write.csv(go_df, "Gene_GO_annotations.csv", row.names = FALSE)

