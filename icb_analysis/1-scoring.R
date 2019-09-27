library(reticulate)
use_condaenv(condaenv = "pypresent", required = TRUE)
# py_config()
bt = import_builtins()
sys = import("sys")
sys$path = c(sys$path, '/Users/wsx/Documents/GitHub/pypresent/')

mutation = import("mutation")
allele = import("allele")
patient = import("patient")

Mutation = mutation$Mutation
Allele = allele$Allele
Patient = patient$Patient

# Check
mut = Mutation(20L, 'H', from_file = TRUE, 
               gene_fasta_file = "../data/protein.fa",
               id = 'OR5AC1_X20K', native_aa = "H")
mut
aI = Allele('HLA-A01:01', mhc_class='I')
aI$allele_score(mut)

allelesI = c('HLA-A11:01', 'HLA-A31:35', 'HLA-B07:02', 'HLA-B40:30', 'HLA-C06:19', 'HLA-C07:07')
p = Patient(allelesI, list())
p$patient_score(mut, mhc_class='I')

# It works


# Scoring -----------------------------------------------------------------
vep_out = data.table::fread("vep_out.txt.gz")
# Keep only variants with AA change
vep_out = vep_out %>% 
  filter(ENSP != "-", Protein_position != "-")

vep_mutation = vep_out %>% 
  select(c(colnames(vep_out)[1:11], 
           c("cDNA_position", "CDS_position","Protein_position", "Amino_acids", "ENSP")))
colnames(vep_mutation)[1] = "variant_id"
rm(vep_out)


phenotype = read_csv("phenotype.csv")

proteins = unique(vep_mutation$ENSP)

get_seqs = function(x, proteomes="proteomes/Homo_sapiens.GRCh37.75.pep.all.fa") {
  cmds = paste(
    "sed",
    "-n",
    paste0("'/^>", x, "/,/>/p'"),
    proteomes, 
    "|",
    "grep -v '>'"
  )
  #message(cmds)
  Seqs = system(cmds, intern = TRUE)
  Seqs = paste0(Seqs, collapse = "")
  Seqs
}

library(furrr)
plan(multiprocess)

protein_seqs = future_map_chr(proteins, get_seqs, .progress = TRUE)
Proteins = dplyr::tibble(
  ENSP = proteins, 
  Seqs = protein_seqs
)
save(protein_seqs, file = "protein_seqs.RData")

RES = vep_mutation %>% 
  mutate(patient = stringr::str_split(variant_id, "-", simplify = T)[, 1]) %>% 
  left_join(phenotype, by = "patient") %>% 
  left_join(Proteins, by = "ENSP")

save(RES, file = "RES.RData")


# Scan all variants and score them
# elf, residue, aa, from_file=True, gene_fasta_file='', gene_sequence='',
# id='mutationID', native_aa=None, native=False

scoring = function(AA_position, AA, Seqs, id, HLA) {
  # AA_postion is an integer
  AA = ifelse(nchar(AA) == 1, AA, substr(AA, 3, 3))
  #sed -n '/^>ENSP00000317992/,/>/p' Homo_sapiens.GRCh37.75.pep.all.fa | grep -v '>'
  
  #message(Seqs)
  #message("Protein length:", nchar(Seqs), " vs AA:", AA_position)
  
  # patient_id = stringr::str_split(id, "-", simplify = T)[1, 1]
  # HLA = phenotype %>% 
  #   dplyr::filter(patient == patient_id) %>% 
  #   dplyr::pull(HLA)
  HLA = str_split(HLA, ";", simplify = T) %>% as.character()
  #print(HLA)
  
  score = tryCatch(
    {
      mut = Mutation(AA_position, AA, from_file=FALSE, gene_sequence=Seqs, id=id)
      p = Patient(HLA, list())
      p$patient_score(mut, mhc_class='I')
    }, error = function(e) {
      cat(AA_position, "vs", nchar(Seqs), "\n", sep = " ", file = "run.log", append = TRUE)
      NA
    }
  )
  cat(id, paste0(score, "\n"), file = "run.tsv", sep = "\t", append = TRUE)
  score
}

# 1
system.time(
  RES %>% 
    head(10) %>% 
    rowwise() %>% 
    do(PHBR = scoring(AA_position = as.integer(.data$Protein_position), 
                      AA = .data$Amino_acids, 
                      Seqs = .data$Seqs, 
                      id = .data$variant_id,
                      HLA = .data$HLA)) -> zz1
)



# 2
library(furrr)
plan(multiprocess)

system.time(
  {
    nbs = 10
    pmap_dbl(list(AA_position = as.integer(RES$Protein_position)[1:nbs],
                  AA = RES$Amino_acids[1:nbs],
                  Seqs = RES$Seqs[1:nbs],
                  id = RES$variant_id[1:nbs],
                  HLA =RES$HLA[1:nbs]), scoring) -> zz2
  }
)

# 3
# Use this
library(parallel)
system.time(
  {
    #nbs = nrow(RES)
    nbs = 10
    PHBR = mcmapply(scoring, AA_position = as.integer(RES$Protein_position)[1:nbs],
                    AA = RES$Amino_acids[1:nbs],
                    Seqs = RES$Seqs[1:nbs],
                    id = RES$variant_id[1:nbs],
                    HLA =RES$HLA[1:nbs], mc.cores = 6)
  }
)

#save(PHBR, file = "PHBR.RData")

