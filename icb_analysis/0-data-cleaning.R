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


ph = bind_rows(van2015 %>% 
                 mutate(study="Van2015"),
               snyder2017 %>% 
                 mutate(study="snyder2017"))

write_csv(ph, path = "phenotype.csv")

mut = bind_rows(van2015.mut, snyder2017.mut %>% 
                  mutate(patient = as.character(patient),
                         Chromosome = as.character(Chromosome))) %>% 
  filter(Chromosome %in% c(1:22, "X", "Y"))

library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens

vep_mut = mut %>% 
  mutate_if(is.numeric, as.integer)

ref = getSeq(Hsapiens, 
             paste0("chr", vep_mut$Chromosome), 
             start = vep_mut$Start_position,
             end = vep_mut$End_position)

ref_seq = as.character(ref)
vep_mut = vep_mut %>% 
  mutate(ref = ref_seq)

all(vep_mut$Reference_Allele == vep_mut$ref)
# So all mutation are based on strand "+"

vep_mut = vep_mut %>% 
  mutate(
    chromosome = Chromosome,
    start = Start_position,
    end = End_position,
    allele = paste(Reference_Allele, Tumor_Seq_Allele2, sep = "/"),
    strand = "+",
    identifier = paste(patient, Chromosome, Start_position, sep = "-")
  ) %>% 
  select(chromosome, start, end, allele, strand, identifier)

vep_mut = vep_mut %>% 
  arrange(chromosome, start)

write_tsv(vep_mut, path = "vep_input.tsv", col_names = FALSE)
# Submit this data to VEP web
