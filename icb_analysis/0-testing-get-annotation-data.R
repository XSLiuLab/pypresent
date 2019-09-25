library(biomaRt)

listEnsemblArchives()
listMarts(host = 'http://grch37.ensembl.org')
ensembl = useMart('http://grch37.ensembl.org', 
                  biomart = "ENSEMBL_MART_ENSEMBL", 
                  dataset = "hsapiens_gene_ensembl")


protein = getSequence(id=c("TP53"),
                      type="hgnc_symbol",
                      seqType="peptide", 
                      mart=ensembl, 
                      verbose = T)
protein

getSequence(chromosome = 21,
            start = 43702476,
            end = 43702478,
            type="hgnc_symbol",
            seqType="peptide", 
            mart=ensembl)

getBM(c("hgnc_symbol", "chromosome_name", "strand",
        "start_position", "end_position", 
        "transcript_start", "transcript_end", 
        "transcription_start_site", "peptide"),
      "hgnc_symbol", "TP53", ensembl) -> zz


# library(biomartr)
# HS.proteome.refseq <- getProteome( db       = "refseq", 
#                                    organism = "Homo sapiens",
#                                    path     = "proteomes", 
#                                    gunzip = TRUE)
# 
# 
download.file("ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.75.pep.all.fa.gz",
              destfile = "proteomes/Homo_sapiens.GRCh37.75.pep.all.fa.gz")
system("gunzip proteomes/Homo_sapiens.GRCh37.75.pep.all.fa.gz")
