This repository includes code accompanying manuscript: 

## **Convergent proteogenomic evidence prioritises five causal proteins and new drug targets for major psychiatric disorders** 

### **Authors:**  
Ruta Margelyte and Christina Dardani

### **Contents:**
- 1_1_prot_covar_imputation.R [link](1_1_prot_covar_imputation.R)
- 1_2_xs_longit_analyses.R [link](1_2_xs_longit_analyses.R)
- 1_3_react_kegg_enrichment.R [link](1_3_react_kegg_enrichment.R)
- 2_1_genetic_instrument_region_extraction_colocalisation.R [link](2_1_genetic_instrument_region_extraction_colocalisation.R)
- 2_2_MR.R [link](2_2_MR.R)
- 2_3_bidirectional_MR.R [link](2_3_bidirectional_MR.R)
- 3_1_IPMC_extraction.R [link](3_1_IPMC_extraction.R)
- 3_2_OMIM_extraction.R [link](3_2_OMIM_extraction.R)
- 3_3_tractability_extraction.R [link](3_3_tractability_extraction.R)

### **Data availability:** 
- Data supporting the results of the present study are available from the UKB (https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access) to researchers with UKB approval.
- UKB blood pQTL data: http://ukb-ppp.gwas.eu
- Brain pQTL data at Synapse portal: https://www.synapse.org/Synapse:syn23627957
- Blood eQTL data: https://www.eqtlgen.org/phase1.html
- Brain cortex eQTL data at MetaBrain platform: https://www.metabrain.nl/
- GWAS data on anxiety, bipolar, and schizophrenia: https://pgc.unc.edu/for-researchers/download-results/
- GWAS data on depression: https://ipsych.dk/en/research/downloads/
- ExWAS data (psychiatric disorders): https://www.ebi.ac.uk/gwas/publications/34662886
- Therapeutic tractability and clinical trials data at the Open Targets Platform: https://platform.opentargets.org/
- GO annotations: https://geneontology.org/docs/go-annotations/
- Tissue expression at The Human Protein Atlas: https://www.proteinatlas.org/
- KEGG pathways: https://www.genome.jp/kegg/
- REACTOME pathways: https://reactome.org/
- IMPC annotations: https://www.mousephenotype.org/
- OMIM annotations: https://www.omim.org/

### **Software:**
- Cross-sectional and longitudinal cohort analyses used R (v4.4.2) packages stats (v4.4.2), survival (v4.4.2), lspline (1.0-0), and ashr (v2.2-63). 
- Missing data imputation was done with R packages mice (v3.18.0) [84] and missForest (v1.5) [42]. 
- Enrichment analyses were performed with R packages ReactomePA (v1.50.0) and clusterProfiler (v4.14.6) [85]. 
- Genetic analyses: Blood plasma pQTL data and blood cell derived eQTL data were extracted and processed using R package gwasvcf (v1.0) (https://github.com/MRCIEU/gwasvcf) [86]. The summary data from MetaBrain were lifted over from GRCh38 to GRCh37 using the UCSC liftover tool to match the build of the rest of the data [87]. Two-sample MR, Steiger filtering, and bi-directional MR analyses were conducted using functions from R packages TwoSampleMR (v0.5.6) (https://github.com/MRCIEU/TwoSampleMR) and mrpipeline (v1.0) (https://github.com/jwr-git/mrpipeline) [88]. The PWCoCo algorithm was implemented using the Pair-Wise Conditional analysis and Colocalisation analysis package v1.0 (https://github.com/jwr-git/pwcoco) [67]. 
- Data on the therapeutic tractability, known drugs, and clinical trials were extracted using functions from R package otargen (v2.0.1) [89]. 
- Gene Ontology (GO) and OMIM annotations were retrieved using R package biomaRt (v2.62.1).  
- Mouse knockout phenotype data were obtained from the IMPC database via its public API. Human-to-mouse orthologs were mapped using R package homologene (v1.48.6).  
- Data was visualised using R package ggplot2 (v4.0.1). 
- Genetic analyses and protein imputation were carried out using the computational facilities of the Advanced Computing Research Centre of the University of Bristol (http://www.bris.ac.uk/acrc/).  
