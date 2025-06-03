#metadata prep.R
#This file prepares the metadata for GSE146639 (microglia) and GSE174367 (brain) for use in project.
#Resulting files are also uploaded to the GitHub

library(tidyverse)
library(readxl)
library(writexl)

#replace with your base path here
base_path <- "path/to/project/folder"

#GSE146639 metadata
GSE146639_matrix <- read_excel(
  file.path(base_path, "GSE146639_series_matrix.xlsx"), 
  col_names = FALSE, skip=32)

#format and transpose matrix
colnames(GSE146639_matrix) <- GSE146639_matrix[[1]]
GSE146639_transposed<- t(GSE146639_matrix)
GSE146639_transposed <- as.data.frame(GSE146639_transposed)
colnames(GSE146639_transposed) <- GSE146639_transposed[1, ]  
GSE146639_transposed <- GSE146639_transposed[-1, ]
colnames(GSE146639_transposed) <- make.names(colnames(GSE146639_transposed), unique = TRUE)

#Select and format relevant columns
GSE146639_metadata <- GSE146639_transposed %>% 
  select(
    `X.Sample_title`,
    `X.Sample_geo_accession`,
    `X.Sample_characteristics_ch1.1`,
    `X.Sample_characteristics_ch1.2`,
    `X.Sample_characteristics_ch1.4`,
    `X.Sample_characteristics_ch1.5`,
    `X.Sample_characteristics_ch1.6`,
    `X.Sample_characteristics_ch1.7`) %>%
  rename(
    Sample_id = `X.Sample_geo_accession`,
    donor_group = `X.Sample_characteristics_ch1.1`,
    donor_id = `X.Sample_characteristics_ch1.2`,
    brain_region = `X.Sample_characteristics_ch1.4`,
    age = `X.Sample_characteristics_ch1.6`,
    gender = `X.Sample_characteristics_ch1.7`,
    bulk_single = `X.Sample_title`,
    seq_pool = `X.Sample_characteristics_ch1.5`)%>%
  mutate(across(everything(), ~ gsub("^[^:]+: ", "", .)),
         donor_group = if_else(donor_group == "CTR+", "CTR", donor_group))%>%
  filter(grepl("bulk", bulk_single, ignore.case = TRUE))

rownames(GSE146639_metadata) <- NULL
#save metadata
write_xlsx(GSE146639_metadata, file.path(base_path, "GSE146639_metadata.xlsx"))

#GSE174367 metadata
GSE174367_matrix <- read_excel(
  file.path(base_path, "GSE174367_series_matrix.xlsx"),
  col_names = FALSE, skip = 34)

#format and transpose matrix
GSE174367_transposed<- as.data.frame(t(GSE174367_matrix))
colnames(GSE174367_transposed) <- GSE174367_matrix[[1]]
GSE174367_transposed <- GSE174367_transposed[-1, ] 
colnames(GSE174367_transposed) <- make.names(colnames(GSE174367_transposed), unique = TRUE)

#select and format relevant columns
GSE174367_metadata <- GSE174367_transposed %>% 
  select(
    `X.Sample_title`,
    `X.Sample_geo_accession`,
    `X.Sample_characteristics_ch1`,
    `X.Sample_characteristics_ch1.1`,
    `X.Sample_characteristics_ch1.3`,
    `X.Sample_characteristics_ch1.4`) %>%
  rename(
    sample = `X.Sample_geo_accession`,
    donor_group = `X.Sample_characteristics_ch1.3`,
    brain_region = `X.Sample_characteristics_ch1.4`,
    age = `X.Sample_characteristics_ch1`,
    gender = `X.Sample_characteristics_ch1.1`,
    bulk_single = `X.Sample_title`) %>%
  mutate(across(everything(), ~ gsub("^[^:]+: ", "", .)),
         donor_group = if_else(donor_group == "Control", "CTR", donor_group))%>%
  extract(bulk_single, into = c("Sample_id", "bulk_single"), 
          regex = "^(.*?)\\s*\\((.*?)\\)$")%>%
  filter(grepl("bulk", bulk_single, ignore.case = TRUE))

rownames(GSE174367_metadata) <- NULL
#save metadata
write_xlsx(GSE174367_metadata, file.path(base_path, "GSE174367_metadata.xlsx"))
