library(biomaRt)

listEnsemblArchives()
listMarts(host = 'http://grch37.ensembl.org')
ensembl = useMart('http://grch37.ensembl.org', 
                  biomart = "ENSEMBL_MART_ENSEMBL", 
                  dataset = "hsapiens_gene_ensembl")


protein = getSequence(id=c("CDC20"),
                      type="hgnc_symbol",
                      seqType="peptide", 
                      mart=ensembl, 
                      verbose = T)
protein

getBM(c("hgnc_symbol", "chromosome_name",
        "start_position", "end_position"),
      "hgnc_symbol", "TP53", ensembl)


# library(biomartr)
# HS.proteome.refseq <- getProteome( db       = "refseq", 
#                                    organism = "Homo sapiens",
#                                    path     = "proteomes", 
#                                    gunzip = TRUE)
# 
# 
# download.file("ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.75.pep.all.fa.gz",
#               destfile = "proteomes/Homo_sapiens.GRCh37.75.pep.all.fa.gz")
# system("gunzip proteomes/Homo_sapiens.GRCh37.75.pep.all.fa.gz")
