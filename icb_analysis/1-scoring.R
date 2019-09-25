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
# It works


# Scoring -----------------------------------------------------------------



# mut = Mutation(21, 'K', from_file=False, gene_sequence=sequence, native_aa='P', id='OR5AC1_P21K')
# 
# 
# allelesI = ['HLA-A11:01', 'HLA-A31:35', 'HLA-B07:02', 'HLA-B40:30', 'HLA-C06:19', 'HLA-C07:07']
# allelesII = ['DRB1_1114', 'DRB1_1301', 'DRB1_1114', 'DRB1_1301', \
#              'HLA-DQA10511-DQB10308', 'HLA-DQA10509-DQB10628', \
#              'HLA-DQA10511-DQB10628', 'HLA-DQA10509-DQB10308', \
#              'HLA-DPA10301-DPB19201', 'HLA-DPA10103-DPB13401', \
#              'HLA-DPA10301-DPB13401', 'HLA-DPA10103-DPB19201']
# p = Patient(allelesI, allelesII)
