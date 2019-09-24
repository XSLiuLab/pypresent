library(readxl)
library(tidyverse)


van2015 = read_excel("data/Science2015_Clinical.xlsx")
van2015 = van2015 %>% 
  mutate(HLA = paste(
    hla.a1, hla.a2,
    hla.b1, hla.b2,
    hla.c1, hla.c2, 
    sep = ";"
  )) %>% 
  select(patient, HLA)

van2015.mut = read_excel("data/Science2015_TMB_List_AllPatients.xlsx")  

van2015.mut = van2015.mut %>% 
  filter(Variant_Type == "SNP") %>% 
  select(patient, Hugo_Symbol, Chromosome, 
         Start_position, End_position, 
         Reference_Allele, Tumor_Seq_Allele2) %>% 
  filter(!is.na(Hugo_Symbol))

write_csv(van2015, "van2015_sample.csv")
write_csv(van2015.mut, "van2015_mutation.csv")


snyder2017 = read_csv("data/Snyder2017_clinical.csv")
snyder2017 = snyder2017 %>% 
  select(patient_id, hla_allele_list) %>% 
  rename(patient = patient_id, 
         HLA = hla_allele_list) %>% 
  mutate(
    HLA = HLA %>% 
      str_remove_all("\\[") %>% 
      str_remove_all("\\]") %>% 
      str_remove_all("\\*") %>% 
      str_remove_all(" ") %>% 
      str_replace_all(",", ";") %>% 
      str_replace_all("([A-Z])", "HLA-\\1")
  )

snyder2017.mut = read_csv("data/Snyder2017_variants.csv")
snyder2017.mut = snyder2017.mut %>% 
  mutate(
    patient = patient_id,
    Hugo_Symbol = gene_name,
    Chromosome = chr, 
    Start_position = start,
    End_position = start,
    Reference_Allele = ref, 
    Tumor_Seq_Allele2 = alt
  ) %>% 
  select(patient, Hugo_Symbol, Chromosome, 
         Start_position, End_position, 
         Reference_Allele, Tumor_Seq_Allele2) %>% 
  filter(!is.na(Hugo_Symbol))

write_csv(snyder2017, "snyder2017_sample.csv")
write_csv(snyder2017.mut, "snyder2017_mutation.csv")
