

library(missForest)
library(doParallel)
library(mice)

####protein random forest imputation####

    setwd("/user/work/rmusername")

    
    #random forest imputation - per panel + age and sex 
    
    #list_datasets <- ls(pattern = "^olink_values_")
    #list_datasets <- c("olink_values_cardiometabolic","olink_values_cardiometabolic_ii",
    #                    "olink_values_inflammation","olink_values_inflammation_ii",
    #                    "olink_values_neurology","olink_values_neurology_ii",     
    #                    "olink_values_oncology","olink_values_oncology_ii" 
    
    # Load an RDS files
    
    args <- commandArgs(trailingOnly=TRUE)
    dataset_name <- args[1]
    
    
    file_path <- paste0("/user/work/rmusername/", dataset_name, ".rds")
    dataset <- readRDS(file_path)
    olink_part_covar<-readRDS("/user/work/rmusername/olink_part_covar.rds")
    
    
    #run imputation
    
    dataset<- merge(dataset, subset(olink_part_covar, select = c(eid, age, sex)), by = "eid")
    row.names(dataset) <- dataset$eid 
    
    gc()
    doParallel::registerDoParallel(cores = 24)
    m.forest <- missForest(dataset, verbose = T, replace = F, parallelize = "forests", ntree = 50)
    gc()
    
    
    new_name <- sub("^olink_values_", "olink_valuesimp3_", dataset_name)
    save(m.forest, file = paste0(new_name,".RData"))
    m.forest$OOBerror
    
    
#####covariate imputation mice####
    
    
    data.frame(Column_Number = seq_along(part_covar_olink), Column_Name = names(part_covar_olink))
    
    a<-colSums(is.na(part_covar_olink))
    a
    
    a <- colSums(is.na(part_covar_olink[, c(1,3:7,9:14,18:37,40,43,46,49,51,55,58,61, 62,64,66,68,69,70)]))
    a
    
    
    data.frame(Column_Number = seq_along(part_outc_olink), Column_Name = names(part_outc_olink))
    
    a<-colSums(is.na(part_outc_olink))
    a
    
    a <- colSums(is.na(part_outc_olink[, c(1,11,12,18,19,25,26,32,33)]))
    a
    
    
    a<-data.frame(part_covar_olink[, c(1,3:7,9:14,18:37,40,43,46,49,52,55,58,61, 62,64,66,68,69,70)],part_outc_olink[, c(11,12,18,19,25,26,32,33)])
    a$sex<-as.factor(a$sex)
    a$ethnicity<-as.factor(a$ethnicity)
    a$alcoholstatus<-as.factor(a$alcoholstatus)
    a$smokingstatus<-as.factor(a$smokingstatus)
    a$education<-as.factor(a$education)
    a$employment<-as.factor(a$employment)
    a$batch<-as.factor(a$batch)
    a$region<-as.factor(a$region)
    a$assessment_season4<-as.factor(a$assessment_season4)
    a$smokingstatus<-as.factor(a$smokingstatus)
    a$smokingstatus<-as.factor(a$smokingstatus)
    
    a[paste0("p22009_a", 1:20)] <- lapply(a[paste0("p22009_a", 1:20)], as.numeric)
    colSums(is.na(a))
    
    pred_matrix <- make.predictorMatrix(a)
    col <- which(colnames(a) == "eid")
    pred_matrix[, col] <- 0
    pred_matrix[col, ] <- 0
    
    init <- mice(a, predictorMatrix = pred_matrix, maxit = 0, dryrun = TRUE)  # Perform a dry run to initialize
    meth <- init$method         # Extract the method vector
    meth
    
    #meth["age"] <- "pmm"       # Predictive Mean Matching for 'age'
    #meth["sex"] <- "logreg" # Logistic Regression for 'gender'
    #meth["eid"] <- ""  # Exclude 'eid' from imputation
    
    imputed_data <- mice(a, predictorMatrix = pred_matrix, method = meth, m = 10,  maxit = 10, seed = 123)
    save(imputed_data, file = "imputed_data.RData")
    
    imputed_data_c <- complete(imputed_data, "long", include = TRUE)  # Convert imputed data to long format
    
    imputed_data_c<-data.frame(imputed_data_c[, c(1:5,11,43:45,47:56 )]) #age+sex+socio (tdi,eth,ed)+prot+outc
    imputed_data_c$outcome <- imputed_data_c$depression_prevalent
    imputed_data_m <- as.mids(imputed_data_c)  # Convert back to mids object
