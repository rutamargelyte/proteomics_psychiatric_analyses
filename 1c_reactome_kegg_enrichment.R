
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)


####select UKB cross-sectional and longitudinal  results####

    a <- results_all %>%
      filter(adjustment=="adj4") %>%
      filter(imputation=="imp3") %>%
      filter(outcome_diagnosis=="anxiety") %>%
      filter(outcome_type=="i") %>%
      mutate(pass = if_else(
        (outcome %in% c("anxiety_prevalent", "depression_prevalent", 
                        "anxiety_incident", "depression_incident",
                        "bipolar_prevalent", "psychotic_prevalent") & 
           (Ash_lfsr_Significant_001 == TRUE | FDR_Significant_001 == TRUE)) |
          (outcome %in% c("bipolar_incident", "psychotic_incident") & 
             (Ash_lfsr_Significant_05 == TRUE | FDR_Significant_05 == TRUE)),
        1, 0 )) %>%
      dplyr::select(Protein, outcome_diagnosis, Estimate, P_value, pass)
    
    b <- results_all %>%
      filter(adjustment=="adj4") %>%
      filter(imputation=="imp3") %>%
      filter(outcome_diagnosis=="anxiety") %>%
      filter(outcome_type=="p") %>%
      mutate(pass = if_else(
        (outcome %in% c("anxiety_prevalent", "depression_prevalent", 
                        "anxiety_incident", "depression_incident",
                        "bipolar_prevalent", "psychotic_prevalent") & 
           (Ash_lfsr_Significant_001 == TRUE | FDR_Significant_001 == TRUE)) |
          (outcome %in% c("bipolar_incident", "psychotic_incident") & 
             (Ash_lfsr_Significant_05 == TRUE | FDR_Significant_05 == TRUE)),
        1, 0 )) %>%
      dplyr::select(Protein, outcome_diagnosis, Estimate, P_value, pass)
    
    
    merged <- a %>%
      dplyr::select(Protein, Estimate_inc = Estimate, P_value_inc = P_value, pass_inc = pass) %>%
      inner_join(
        b %>%
          dplyr::select(Protein, Estimate_prev = Estimate, P_value_prev = P_value, pass_prev = pass),
        by = "Protein"
      ) 
    
    merged$Protein<-toupper(merged$Protein)
    
    
    merged <- merged %>%   
      mutate(class = case_when(
        pass_prev == 1 & pass_inc == 0 ~ "prev_only",
        pass_prev == 1 & pass_inc == 1 ~ "prev_inc",
        pass_prev == 0 & pass_inc == 1 ~ "inc_only",
        TRUE ~ "none"
      ))
    
    table(merged$class)
    
    
    
    # #proteins/genes shared across disorders
    # c <- results_all %>%
    #   filter(adjustment=="adj4") %>%
    #   filter(imputation=="imp3") %>%
    #   #filter(outcome_diagnosis=="psychotic") %>%
    #   #filter(outcome_type=="i") %>%
    #   mutate(pass = if_else(
    #     (outcome %in% c("anxiety_prevalent", "depression_prevalent",
    #                     "anxiety_incident", "depression_incident",
    #                     "bipolar_prevalent", "psychotic_prevalent") &
    #        (Ash_lfsr_Significant_001 == TRUE | FDR_Significant_001 == TRUE)) |
    #       (outcome %in% c("bipolar_incident", "psychotic_incident") &
    #          (Ash_lfsr_Significant_05 == TRUE | FDR_Significant_05 == TRUE)),
    #     1, 0 )) %>%
    #   dplyr::select(Protein, outcome_diagnosis, pass)
    # 
    # c$Protein<-toupper(c$Protein)
    # 
    # c <- c %>% filter(pass=="1") %>%
    #   distinct() %>%
    #   group_by(Protein) %>%
    #   mutate(shared = sum(pass, na.rm = TRUE)) %>%
    #   ungroup() %>%
    #   dplyr::select(-outcome_diagnosis) %>%
    #   distinct()
    # 
    # table(c$shared)


####prepare gene sets####
    
  #gene background (full protein list)
  
  
  background <- c(merged$Protein)
  background <- gsub("^HLA_A$", "HLA-A", background)
  background <- gsub("^HLA_DRA$", "HLA-DRA", background)
  background <- gsub("^HLA_E$", "HLA-E", background)
  background <- gsub("^ERVV_1$", "ERVV-1", background)
  background <- ifelse(background == "C19ORF12", "C19orf12",
                       ifelse(background == "C9ORF40", "C9orf40",
                              ifelse(background == "C2ORF69", "C2orf69",
                                     ifelse(background == "C7ORF50", "C7orf50",
                                            background))))
  background <- ifelse(background == "NTPROBNP", "NPPB",background)
  background <- ifelse(background == "BAP18", "BACC1",background)
  background <- ifelse(background == "GBA", "GBA1",background)
  background <- ifelse(background == "DUSP13", "DUSP13A",background)
  background <- unlist(strsplit(background, "_"))
  #back1$background[duplicated(back1$background)]
  background <- unique(background)
  
  
  #genes of interest
  
  #genes <- c(merged$Protein[merged$class=="Yes"])
  genes <- c(merged$Protein[merged$class=="inc_only"])
  #genes <- c(merged$Protein[merged$class=="prev_only"])
  #genes <- c(merged$Protein[merged$class=="prev_inc"])
  #genes <- c(merged$Protein[merged$class=="prev_inc"|merged$class=="inc_only"])
  #genes <- c(merged$Protein[merged$class=="prev_inc"|merged$class=="prev_only"])
  #genes <- c(merged$Protein[merged$class=="prev_inc"|merged$class=="inc_only"|merged$class=="prev_only"])
  #genes <- c(c$Protein[c$shared=="3"|c$shared=="4"])
  #genes <- c(c$Protein[c$shared=="4"])
  genes <- gsub("^HLA_A$", "HLA-A", genes)
  genes <- gsub("^HLA_DRA$", "HLA-DRA", genes)
  genes <- gsub("^HLA_E$", "HLA-E", genes)
  genes <- gsub("^ERVV_1$", "ERVV-1", genes)
  genes <- ifelse(genes == "C19ORF12", "C19orf12",
                  ifelse(genes == "C9ORF40", "C9orf40",
                         ifelse(genes == "C2ORF69", "C2orf69",
                                ifelse(genes == "C7ORF50", "C7orf50",
                                       genes))))
  genes <- ifelse(genes == "NTPROBNP", "NPPB",genes)
  genes <- ifelse(genes == "BAP18", "BACC1",genes)
  genes <- ifelse(genes == "GBA", "GBA1",genes)
  genes <- ifelse(genes == "DUSP13", "DUSP13A",genes)
  genes <- unlist(strsplit(genes, "_"))
  genes <- unique(genes)
  
  # Convert genes
  
  gene_df <- bitr(genes,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
  
  unmapped <- genes[!genes %in% gene_df$SYMBOL]
  unmapped
  
  # stop / skip if nothing unmapped
  if (length(unmapped) > 0) {
    
    res <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unmapped,
      keytype = "ALIAS",
      columns = "ENTREZID"
    )
    
    res <- dplyr::rename(res, SYMBOL = ALIAS)
    
    gene_df <- dplyr::bind_rows(gene_df, res)
    
  }
  
  # Convert background
  
  bg_df <- bitr(background,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
  
  unmapped <- background[!background %in% bg_df$SYMBOL]
  unmapped
  
  res<-AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unmapped, 
    keytype = "ALIAS",
    columns = "ENTREZID")
  res <- res %>% dplyr::rename(SYMBOL = ALIAS)
  
  bg_df<-rbind(bg_df,res)
  
  rm(res)

  
  gene_ids <- gene_df$ENTREZID
  bg_ids   <- bg_df$ENTREZID


  
####REACTOME enrichment analysis####
  
  reactome_res <- enrichPathway(
    gene          = gene_ids,
    organism      = "human",
    universe      = bg_ids,        # <-- custom background
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
  
  reactome_df <- as.data.frame(reactome_res)
  table(reactome_df$pvalue<0.05)
  table(reactome_df$pvalue<0.01) #9po 10io
  table(reactome_df$p.adjust<0.05)
  
  
  write.csv(reactome_df, "react_df_anx_full_inconly.csv", row.names = FALSE)
  

####KEGG enrichment analysis####

    kegg_res <- enrichKEGG(
      gene         = gene_ids,
      organism     = "hsa",        # human
      universe     = bg_ids,       # custom background
      pvalueCutoff = 1,
      qvalueCutoff =  1
    )

    kegg_df <- as.data.frame(kegg_res)
    table(kegg_df$pvalue<0.05)
    table(kegg_df$pvalue<0.01) #9po 10io
    table(kegg_df$p.adjust<0.05)


    # table(kegg_df$category)
    # table(kegg_df$subcategory)
    # table(kegg_df$subcategory,kegg_df$category)

    kegg_df$category<- ifelse(kegg_df$Description=="Virion - Hepatitis viruses", "Genetic Information Processing",  kegg_df$category)
    kegg_df$subcategory<- ifelse(kegg_df$Description=="Virion - Hepatitis viruses", "Information processing in viruses",  kegg_df$subcategory)

    kegg_df$category<- ifelse(kegg_df$Description=="Hormone signaling", "Environmental Information Processing",  kegg_df$category)
    kegg_df$subcategory<- ifelse(kegg_df$Description=="Hormone signaling", "Signaling molecules and interaction",  kegg_df$subcategory)

    kegg_df$category<- ifelse(kegg_df$Description=="Cornified envelope formation", "Organismal Systems",  kegg_df$category)
    kegg_df$subcategory<- ifelse(kegg_df$Description=="Cornified envelope formation", "Development and regeneration",  kegg_df$subcategory)

    kegg_df$category<- ifelse(kegg_df$Description=="Integrin signaling", "Environmental Information Processing",  kegg_df$category)
    kegg_df$subcategory<- ifelse(kegg_df$Description=="Integrin signaling", "Signaling molecules and interaction",  kegg_df$subcategory)

    kegg_df$category<- ifelse(kegg_df$Description=="Neuroactive ligand signaling", "Environmental Information Processing",  kegg_df$category)
    kegg_df$subcategory<- ifelse(kegg_df$Description=="Neuroactive ligand signaling", "Signaling molecules and interaction",  kegg_df$subcategory)

    kegg_df$category<- ifelse(kegg_df$Description=="Cadherin signaling", "Environmental Information Processing",  kegg_df$category)
    kegg_df$subcategory<- ifelse(kegg_df$Description=="Cadherin signaling", "Signaling molecules and interaction",  kegg_df$subcategory)

    kegg_df$category<- ifelse(kegg_df$Description=="IgSF CAM signaling", "Environmental Information Processing",  kegg_df$category)
    kegg_df$subcategory<- ifelse(kegg_df$Description=="IgSF CAM signaling", "Signaling molecules and interaction",  kegg_df$subcategory)


    write.csv(kegg_df, "kegg_df_anx_full_inconly.csv", row.names = FALSE)



