
  library(dplyr)
  library(tidyr)
  library(lspline)
  library(survival)
  
  library(readr)
  library (ashr)
  library(stringr)

  gc() 
  options(max.print = 1000000) 
  options(scipen = 9999) 

####load UKB data####

  load("part_outc.RData") # psychiatric outcomes
  load("part_outc_contr_olink.RData") #outcome control conditions
  load("part_covar.RData") # covariates/exclusions
  load("olink_valuesimp0.RData") #unimputed protein data
  olink_valuesimp3<-readRDS("olink_valuesimp3.rds") #imputed protein data
  protein_splines_assessment_hour <- read.csv("protein_splines_assessment_hour.csv") #transformation parameters for covariates
  protein_splines_fasting_time <- read.csv("protein_splines_fasting_time.csv") #transformation parameters for covariates
  load("imputed_data.RData") #imputed covariate data
  load("part_outc_contr_olink.RData") #control conditions

####analysis_diagnoses prevalent (cross-sectional)####

  #main outcomes
  binary_outcomes <- c("depression_prevalent","anxiety_prevalent", "psychotic_prevalent", "bipolar_prevalent")
  
  ##pos neg controls
  #binary_outcomes <- c( "ra_prevalent","refacc_prevalent","cataract_prevalent") 
  
  results_all_outcomes <- list()
  
  
  # Loop through each binary outcome
  for (outcome_var in binary_outcomes) {
    results_list <- list()
    #}
    
    print(outcome_var)
    
    #cohort exclusions (cross-sectional)
    excl1 <- part_outc_olink[[outcome_var]] == 0 & part_outc_olink[[paste0(outcome_var, "_other")]] == 1 #exclude with other diagnosis from outc==0
    excl2 <- (part_covar_olink$crp > 10| is.na(part_covar_olink$crp)==TRUE) #crp < 10 only and no NAs
    valid_rows <- !(excl1 | excl2)
    #valid_rows <- TRUE # take all records!!
    
    ###datasets - for unadjusted and minimal model
    subset_outcome <- part_outc_olink[valid_rows, ]
    subset_protein_values <- olink_valuesimp3[valid_rows, ]  #rf imputed proteins
    #subset_protein_values <- olink_valuesimp0[valid_rows, ] #unimputed porteins
    subset_covar <- part_covar_olink[valid_rows, ]
    
    
    # ###datasets - for fully adjusted models
    # imputed_data_c <- complete(imputed_data, "long", include = TRUE)  # Convert imputed data to long format
    # imputed_data_c<-data.frame(imputed_data_c[, c(1:10,11,43:45,47:56 )])   #age+sex+socio(tdi,eth,ed)+life(bmi,phys,smok,alc)+prot+outc
    # imputed_data_c <- imputed_data_c[valid_rows, ]
  
    no<-0
    # Loop through each protein column
    for (protein in colnames(subset_protein_values)[-1]) { 
      
      protein<-"gdf15"
      
      ##data - for unadjusted and minimal model
      protein_values <- subset_protein_values[[protein]]
      outcome <- subset_outcome[[outcome_var]]
      age <- subset_covar$age
      sex <- subset_covar$sex
      assessment_season4 <- subset_covar$assessment_season4
      assessment_hour <- subset_covar$assessment_hour
      fasting_time <- subset_covar$fasting_time
      
      
      # ##data - for fully adjusted models
      # imputed_data_cv<-imputed_data_c
      # imputed_data_cv[[outcome_var]] <-rep(subset_outcome[[outcome_var]], 6) # for strict outcome
      # imputed_data_cv$outcome <- imputed_data_cv[[outcome_var]] # Add outcome
      # imputed_data_cv$protein <- rep(subset_protein_values[[protein]], 6) # Add exposure
      # imputed_data_m <- as.mids(imputed_data_cv)  # Convert back to mids object
      
      #select prarameters for covariate spline tranformation
      knot1 <- protein_splines_assessment_hour$Best_Knot[protein_splines_assessment_hour$Protein == protein]
      splines1 <- protein_splines_assessment_hour$splines[protein_splines_assessment_hour$Protein == protein]
      
      knot2 <- protein_splines_fasting_time$Best_Knot[protein_splines_fasting_time$Protein == protein]
      splines2 <- protein_splines_fasting_time$splines[protein_splines_fasting_time$Protein == protein]
      
      no<-no+1
      print(no)
      print(protein)
      
      # Fit the logistic regression model
      
      ####unadjusted
      #model <- glm(outcome ~ protein_values, family = binomial(link = "logit"))
      
      ###minimal: age sex protein related
        if (splines1 == 1&splines2==1) {
          model <- glm(outcome ~ protein_values + age + sex + assessment_season4 + lspline(assessment_hour, knot1) + lspline(fasting_time, knot2), family = binomial(link = "logit"))
        } else if (splines1 == 1&splines2==0){
          model <- glm(outcome ~ protein_values + age + sex + assessment_season4 + lspline(assessment_hour, knot1) + fasting_time, family = binomial(link = "logit"))
        }  else if (splines1 == 0&splines2==1){
          model <- glm(outcome ~ protein_values + age + sex + assessment_season4 + assessment_hour + lspline(fasting_time, knot2), family = binomial(link = "logit"))
        } else {
          model <- glm(outcome ~ protein_values + age + sex + assessment_season4 + assessment_hour + fasting_time, family = binomial(link = "logit"))
        }
      
      # ###full (minimal + demog + lifestyle) ---bit for imputed covariate data
      # if (splines1 == 1&splines2==1) {
      #   model <- with(imputed_data_m, glm(outcome ~ protein + age + sex + tdi + ethnicity + education + bmi + alcoholstatus + smokingstatus + physicalactivity
      #                                     + assessment_season4 + lspline(assessment_hour, knot1) + lspline(fasting_time, knot2) , family = binomial))
      # } else if (splines1 == 1&splines2==0){
      #   model <- with(imputed_data_m, glm(outcome ~ protein + age + sex + tdi + ethnicity + education + bmi + alcoholstatus + smokingstatus + physicalactivity
      #                                     + assessment_season4 + lspline(assessment_hour, knot1) + fasting_time , family = binomial))
      # }  else if (splines1 == 0&splines2==1){
      #   model <- with(imputed_data_m, glm(outcome ~ protein + age + sex + tdi + ethnicity + education + bmi + alcoholstatus + smokingstatus + physicalactivity
      #                                     + assessment_season4 + assessment_hour + lspline(fasting_time, knot2) , family = binomial))
      # } else {
      #   model <- with(imputed_data_m, glm(outcome ~ protein + age + sex + tdi + ethnicity + education + bmi + alcoholstatus + smokingstatus + physicalactivity
      #                                     + assessment_season4 + assessment_hour + fasting_time , family = binomial))
      # }
  
      
      #summary(model)
  
      
      #summary stats - for unadjusted and minimal models
        summary_stats <- summary(model)$coefficients
        protein_result <- data.frame(
          Protein = protein,
          Estimate = summary_stats["protein_values", "Estimate"],
          Std_Error = summary_stats["protein_values", "Std. Error"],
          Z_value = summary_stats["protein_values", "z value"],
          P_value = summary_stats["protein_values", "Pr(>|z|)"],
          row_count = nrow(model$model),
          outcome_count = as.numeric(sum(model$model$outcome))
        )
      
      # ###summary stats - for fully adjusted model
      # summary_stats<-data.frame (summary(pool(model)))
      # summary_stats <- summary_stats[summary_stats$term == "protein", ]
      # protein_result <- data.frame(
      #   Protein = as.character(protein),
      #   Estimate = summary_stats$estimate,
      #   Std_Error = summary_stats$std.error,
      #   Z_value = summary_stats$statistic,
      #   P_value = summary_stats$p.value,
      #   row_count = nrow(imputed_data_m$data),
      #   outcome_count = sum(imputed_data_m$data$outcome)
      # )
      
      # Append results for the protein to the list
      results_list[[protein]] <- protein_result
      
    }
    
    # Combine all protein results into a single data frame for the current outcome
    results_table <- do.call(rbind, results_list)
    
    results_table$P_value_Bonferroni <- 0.05 / (2920 * 4)
    results_table$Bonferroni_Significant <- results_table$P_value < 0.05 / (2920 * 4)
    
    results_table$P_value_FDR <- p.adjust(results_table$P_value, method = "fdr")
    results_table$FDR_Significant <- results_table$P_value_FDR < 0.05
    
    
    # Save the results table to a CSV file for the current outcome
    file_name <- paste0("logistic_results_adj2crp_imp_", outcome_var, ".csv")
    write.csv(results_table, file_name, row.names = FALSE)
    
    # Add the results table to the main list of outcomes
    results_all_outcomes[[outcome_var]] <- results_table
    
    gc()
    
  }

####analysis_diagnoses incident (longitudinal)####
  
  #main outcomes
  binary_outcomes <- c("depression_incident", "psychotic_incident", "bipolar_incident", "anxiety_incident")
  exclusions <- c("depression_prevalent", "psychotic_prevalent", "bipolar_prevalent", "anxiety_prevalent")
  censoring <- c("depression_ttd", "psychotic_ttd", "bipolar_ttd", "anxiety_ttd")
  
  results_all_outcomes <- list()
  
  # Loop through each outcome and corresponding exclusion
  for (i in seq_along(binary_outcomes)) {
    outcome_var <- binary_outcomes[i]
    exclusion_var <- exclusions[i]
    censoring_var <- censoring[i]  
    #}
    
    print(outcome_var)
    
    #cohort exclusions
    excl1 <- part_outc_olink[[exclusion_var]] != 1 &  ceiling(part_outc_olink[[censoring_var]]/30.44)>12  # exclude prevalent case and exclude censored cases within 1 year after baseline
    excl2 <- (part_covar_olink$crp <=10 & is.na(part_covar_olink$crp)==FALSE ) #crp < 10 only and no NAs
    valid_rows <- (excl1 & excl2)
    #valid_rows <- TRUE # take all!!
    
    ###datasets - for unadjusted and minimal models
    subset_outcome <- part_outc_olink[valid_rows, ]
    subset_protein_values <- olink_valuesimp3[valid_rows, ]  # use imputed rf proteins
    #subset_protein_values <- olink_valuesimp0[valid_rows, ] # use unimputed proteins
    subset_covar <- part_covar_olink[valid_rows, ]
    
    
    # ###datasets - for fully adjusted models only
    # imputed_data_c <- complete(imputed_data, "long", include = TRUE)  # Convert imputed data to long format
    # imputed_data_c<-data.frame(imputed_data_c[, c(1:10,11,43:45,47:56 )])   #age+sex+socio(tdi,eth,ed)+life(bmi,phys,smok,alc)+prot+outc
    # imputed_data_c <- imputed_data_c[valid_rows, ]

    
    # Initialize a list to store results for the current outcome
    results_list <- list()
    
    no<-0
    # Loop through each protein column
    for (protein in colnames(subset_protein_values)[-1]) { # Exclude 'eid' column
      
      #protein<-"gdf15"
      
      
      protein_values <- subset_protein_values[[protein]]
      time_to_event <- ceiling(subset_outcome[[censoring_var]]/30.44)-12 #365.25 #/30.44/12 #convert to month to year - now monthly!!!
      event_status <- subset_outcome[[outcome_var]]
     
      ###data - for unadjusted and minimal models
      age <- subset_covar$age
      sex <- subset_covar$sex
      assessment_season4 <- subset_covar$assessment_season4
      assessment_hour <- subset_covar$assessment_hour
      fasting_time <- subset_covar$fasting_time
    
      
      # ###data - for full model
      # imputed_data_cv<-imputed_data_c
      # imputed_data_cv[[outcome_var]] <-rep(subset_outcome[[outcome_var]], 6) # for strict outcome
      # imputed_data_cv$event_status <- imputed_data_cv[[outcome_var]] # Add outcome (event)
      # imputed_data_cv$time_to_event <- rep((ceiling(subset_outcome[[censoring_var]]/30.44)-12), 6) # Add time to event
      # imputed_data_cv$protein <- rep(subset_protein_values[[protein]], 6) # Add exposure
      # imputed_data_m <- as.mids(imputed_data_cv)  # Convert back to mids object
      
      
      knot1 <- protein_splines_assessment_hour$Best_Knot[protein_splines_assessment_hour$Protein == protein]
      splines1 <- protein_splines_assessment_hour$splines[protein_splines_assessment_hour$Protein == protein]
      
      knot2 <- protein_splines_fasting_time$Best_Knot[protein_splines_fasting_time$Protein == protein]
      splines2 <- protein_splines_fasting_time$splines[protein_splines_fasting_time$Protein == protein]
      
      no<-no+1
      print(no)
      print(protein)
      

      # Fit the Cox regression model
      
      ##unadjusted
      #model <- coxph(Surv(time_to_event, event_status) ~ protein_values, data = subset_outcome)
      
      #minimal: age sex protein related
      if (splines1 == 1&splines2==1) {
        model <- coxph(Surv(time_to_event, event_status) ~ protein_values + age + sex + assessment_season4 + lspline(assessment_hour, knot1) + lspline(fasting_time, knot2), data = subset_outcome )
      } else if (splines1 == 1&splines2==0){
        model <- coxph(Surv(time_to_event, event_status) ~ protein_values + age + sex + assessment_season4 + lspline(assessment_hour, knot1) + fasting_time, data = subset_outcome )
      }  else if (splines1 == 0&splines2==1){
        model <- coxph(Surv(time_to_event, event_status) ~ protein_values + age + sex + assessment_season4 + assessment_hour + lspline(fasting_time, knot2), data = subset_outcome )
      } else {
        model <- coxph(Surv(time_to_event, event_status) ~ protein_values + age + sex + assessment_season4 + assessment_hour + fasting_time, data = subset_outcome )
      }
      

      # ###full (minimal + demog + lifestyle) ---bit for imputed covariate data
      # if (splines1 == 1 & splines2 == 1) {
      #   model <- with(imputed_data_m, coxph(Surv(time_to_event, event_status) ~ protein + age + sex + tdi + ethnicity + education + bmi + alcoholstatus + smokingstatus + physicalactivity
      #                                       + assessment_season4 + lspline(assessment_hour, knot1) + lspline(fasting_time, knot2)))
      # } else if (splines1 == 1 & splines2 == 0) {
      #   model <- with(imputed_data_m, coxph(Surv(time_to_event, event_status) ~ protein + age + sex + tdi + ethnicity + education + bmi + alcoholstatus + smokingstatus + physicalactivity
      #                                       + assessment_season4 + lspline(assessment_hour, knot1) + fasting_time))
      # } else if (splines1 == 0 & splines2 == 1) {
      #   model <- with(imputed_data_m, coxph(Surv(time_to_event, event_status) ~ protein + age + sex + tdi + ethnicity + education + bmi + alcoholstatus + smokingstatus + physicalactivity
      #                                       + assessment_season4 + assessment_hour + lspline(fasting_time, knot2)))
      # } else {
      #   model <- with(imputed_data_m, coxph(Surv(time_to_event, event_status) ~ protein + age + sex + tdi + ethnicity + education + bmi + alcoholstatus + smokingstatus + physicalactivity
      #                                       + assessment_season4 + assessment_hour + fasting_time))
      # }
      
      

      # Extract summary statistics
      
      #summary stats - for unadjusted and minimal models
      summary_stats <- summary(model)$coefficients
      protein_result <- data.frame(
        Protein = protein,
        Estimate = summary_stats["protein_values", "coef"],
        Std_Error = summary_stats["protein_values", "se(coef)"],
        Z_value = summary_stats["protein_values", "z"],
        # HR = exp(summary_stats["protein_values", "coef"]),
        # HR_Lower_CI = exp(summary_stats["protein_values", "coef"] - 1.96 * summary_stats["protein_values", "se(coef)"]),
        # HR_Upper_CI = exp(summary_stats["protein_values", "coef"] + 1.96 * summary_stats["protein_values", "se(coef)"]),
        P_value = summary_stats["protein_values", "Pr(>|z|)"],
        row_count = as.numeric(summary(model)$n),
        outcome_count = as.numeric(summary(model)$nevent)
      )
      
      
      # #summary stats - for full model only
      # summary_stats<-data.frame (summary(pool(model)))
      # summary_stats <- summary_stats[summary_stats$term == "protein", ]
      # protein_result <- data.frame(
      #   Protein = as.character(protein),
      #   Estimate = summary_stats$estimate,
      #   Std_Error = summary_stats$std.error,
      #   Z_value = summary_stats$statistic,
      #   P_value = summary_stats$p.value,
      #   row_count = nrow(imputed_data_m$data),
      #   outcome_count = sum(as.numeric(imputed_data_m$data$event_status))
      # )
      
      
      
      # Append results for the protein to the list
      results_list[[protein]] <- protein_result
    }
    
    # Combine all protein results into a single data frame for the current outcome
    results_table <- do.call(rbind, results_list)
    
    
    # Perform multiple testing corrections
    results_table$P_value_Bonferroni <- 0.05 / (2920 * 4)
    results_table$Bonferroni_Significant <- results_table$P_value < 0.05 / (2920 * 4)
    results_table$P_value_FDR <- p.adjust(results_table$P_value, method = "fdr")
    results_table$FDR_Significant <- results_table$P_value_FDR < 0.05
    #results_table$P_value_FDR2 <- results_table$P_value * (2920 * 4) / rank(results_table$P_value, ties.method = "min")
    #results_table$P_value_FDR2 <- pmin(results_table$P_value_FDR2, 1)  # Ensure values do not exceed 1
    #results_table$FDR2_Significant <- results_table$P_value_FDR2 < 0.05
    
    # Save the results table to a CSV file for the current outcome
    file_name <- paste0("cox_results_adj2crp_imp3_", outcome_var, ".csv")
    write.csv(results_table, file_name, row.names = FALSE)
    
    # Add the results table to the main list of outcomes
    results_all_outcomes[[outcome_var]] <- results_table
    
    gc()
    
  }  
  
  
  #cleanup#
  
  #remove values
  rm(list = setdiff(ls(), c(lsf.str(), names(which(sapply(mget(ls()), is.data.frame))))))
  
  rm(subset_outcome)
  rm(subset_covar)
  rm(subset_protein_values)
  rm(protein_result)
  
####merge files for further analyses, derive FDR/LFSR, derive significance flags#### 

  
  ####path and files#
  
  #path <- "/proj_proteomics_psychiatry/rproj/backup_assoc/crp_excl"
  path <- "C:/Users/rm15367/OneDrive - University of Bristol/Desktop/GENEPI/proj_proteomics_psychiatry/rproj/backup_assoc/crp_excl"
  
  file_list <- list.files(path = path, pattern = "results_.*\\.csv$", full.names = TRUE)
  
  file_list <- file_list[grepl("_adj0crp_|_adj2crp_|_adj4crp_", file_list)] 
  
  file_list <- file_list[grepl("_imp3|_imp0", file_list)]
  
  file_list <- file_list[grepl("_anxiety_|_depression_|_psychotic_|_bipolar_", file_list)]
  
  file_list
  
  
  ####import table#
  
  data_list <- lapply(file_list, function(file) {
    df <- read.csv(file)  # Read the CSV file
    
    if ("row_count....summary.model..n" %in% colnames(df)) {
      df <- df %>% rename(row_count = 'row_count....summary.model..n')
    }
    
    if ("outcome_count....summary.model..nevent" %in% colnames(df)) {
      df <- df %>% rename(outcome_count = 'outcome_count....summary.model..nevent')
    }
    
    # Add a column with the source file name
    df$source_file <- basename(file)
    
    return(df)
  })
  
  ####Combine all data frames into one big table and derive descriptive columns#
  
  results_all <- bind_rows(data_list)
  
  # Create 'outcome' column based on conditions
  results_all <- results_all %>%
    dplyr::mutate(outcome = dplyr::case_when(
      str_detect(source_file, "imp0_") ~ str_extract(source_file, "(?<=imp0_).*?(?=\\.csv)"),
      str_detect(source_file, "imp1_") ~ str_extract(source_file, "(?<=imp1_).*?(?=\\.csv)"),
      str_detect(source_file, "imp2_") ~ str_extract(source_file, "(?<=imp2_).*?(?=\\.csv)"),
      str_detect(source_file, "imp3_") ~ str_extract(source_file, "(?<=imp3_).*?(?=\\.csv)"),
      TRUE ~ NA_character_  # Default if none of the conditions are met
    ))
  
  # Create 'outcome_type' column based on conditions
  results_all <- results_all %>%
    mutate(outcome_type = case_when(
      grepl("_prevalent", outcome) ~ "p",
      grepl("_incident", outcome) ~ "i",
      TRUE ~ NA_character_  # Assign NA if none of the conditions are met
    ))
  
  # Create 'outcome_diagnosis' column based on conditions
  results_all <- results_all %>%
    mutate(outcome_diagnosis = case_when(
      grepl("psychotic", outcome) ~ "psychotic",
      grepl("anxiety", outcome) ~ "anxiety",
      grepl("bipolar", outcome) ~ "bipolar",
      grepl("depression", outcome) ~ "depression",
      TRUE ~ NA_character_  # Assign NA if none of the conditions are met
    ))
  
  
  # Create 'imputation' column based on conditions
  results_all <- results_all %>%
    mutate(imputation = case_when(
      grepl("_imp0_", source_file) ~ "imp0",     #none
      grepl("_imp1_", source_file) ~ "imp1",     #minsampe
      grepl("_imp2_", source_file) ~ "imp2",     #pmm
      grepl("_imp3_", source_file) ~ "imp3",     #rf
      TRUE ~ NA_character_  # Assign NA if none of the conditions are met
    ))
  
  
  # Create 'adjustment' column based on conditions
  results_all <- results_all %>%
    mutate(adjustment = case_when(
      grepl("_adj0", source_file) ~ "adj0",    #none
      grepl("_adj1", source_file) ~ "adj1",    #protein related
      grepl("_adj2_", source_file)~ "adj2",    #protein related + age/sex 
      grepl("_adj2a", source_file)~ "adj2a",   #protein related + age/sex + psycomorb
      grepl("_adj3", source_file) ~ "adj3",    #protein related + age/sex + sociodemo
      grepl("_adj4_", source_file) ~ "adj4",   #protein related + age/sex + sociodemo + lifestyle
      grepl("_adj4a_", source_file) ~ "adj4a", #protein related + age/sex + sociodemo + lifestyle +psycomorb
      grepl("_adj5", source_file) ~ "adj5", #protein related + age/sex + sociodemo + lifestyle + autoimmune
      grepl("_adj0crp", source_file) ~ "adj0", #crp>10 excluded
      grepl("_adj2crp", source_file) ~ "adj2", #crp>10 excluded
      grepl("_adj4crp", source_file) ~ "adj4", #crp>10 excluded
      TRUE ~ NA_character_  # Assign NA if none of the conditions are met
    ))
  
  
  #fdr sign thresholds#

  results_all$FDR_Significant_001<-results_all$P_value_FDR <0.001
  colnames(results_all)[colnames(results_all) == "FDR_Significant"] <- "FDR_Significant_05"

  
  #ash (lfsr) fit by outcome#
  
  
  results_all$split<-paste(results_all$outcome_type,results_all$imputation, sep = "_")
  table(results_all$outcome, results_all$split)
  
  results_all$split<-paste(results_all$split,results_all$adjustment, sep = "_")
  table(results_all$outcome, results_all$split)
  
  results_all$outcome2<-ifelse(results_all$outcome_diagnosis=="anxiety", "1", ifelse(results_all$outcome_diagnosis=="depression","1","2"))
  table(results_all$outcome, results_all$outcome2)
  
  results_all$outcome2<-paste(results_all$outcome2,results_all$split, sep = "_")
  table(results_all$outcome2, results_all$outcome)
  
  # Add ash columns to results_all with NA to begin with
  results_all <- results_all %>%
    mutate(
      Estimate_shrunk = NA_real_,
      Std_Error_shrunk = NA_real_,
      lfdr = NA_real_,
      lfsr = NA_real_,
      svalue = NA_real_
    )
  
  # Loop through each outcome
  for (o in unique(results_all$outcome2)) {
    message("Running ash for outcome: ", o)
    
    # Subset data for this outcome
    subset_idx <- which(results_all$outcome2 == o)
    dat <- results_all[subset_idx, ]
    
    # Run ash
    ash_fit <- ash(betahat = dat$Estimate,
                   sebetahat = dat$Std_Error,
                   #mixcompdist = "normal",
                   method = "fdr") #shrink fdr
    
    # Add ash results
    results_all$Estimate_shrunk[subset_idx] <- get_pm(ash_fit)
    results_all$Std_Error_shrunk[subset_idx]   <- get_psd(ash_fit)
    results_all$lfdr[subset_idx]           <- get_lfdr(ash_fit) 
    results_all$lfsr[subset_idx]           <- get_lfsr(ash_fit) 
    results_all$svalue[subset_idx]         <- get_svalue(ash_fit) 
  }
  
  #lfdr - local false discovery rate - posterior probability that the true effect is zero 
  # - Bayesian analogue of the p-value “What’s the chance this association is truly non-zero?”
  #lfsr - local false sign rate - posterior probability that the sign (direction) of the effect is wrong 
  # - More stringent than lfdr — even if the effect is non-zero, is it in the right direction? “Can I trust the direction of this association?”
  #svalue - Empirical Bayes analogue of the q-value, but for false sign rate - Adjusted version of lfsr for multiple testing.
  
  
  results_all$Ash_lfsr_Significant_001 <- results_all$lfsr <0.001
  results_all$Ash_lfsr_Significant_05 <- results_all$lfsr <0.05
  
  results_all <- results_all %>% relocate(Bonferroni_Significant, .before = Ash_lfsr_Significant_001)
  results_all <- results_all %>% relocate(FDR_Significant_05, .before = Ash_lfsr_Significant_001)
  results_all <- results_all %>% relocate(FDR_Significant_001, .before = Ash_lfsr_Significant_001)
  
  rm(dat)
  rm(ash_fit)
  rm(data_list)
  rm(file_list)
  rm(path)
  rm(subset_idx)
  rm(o)
  
  table (results_all$source_file, results_all$Ash_lfsr_Significant_001, useNA = "always")
  
  
  
  
  