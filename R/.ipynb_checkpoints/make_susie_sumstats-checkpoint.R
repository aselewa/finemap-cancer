
source('R/GWAS_funs_libs.R')

args <- commandArgs(trailingOnly = T)

sumstats.fname <- args[1]
prior.fname <- args[2]
loci_qval <- args[3]
output <- args[4]

sumstats <- vroom::vroom(sumstats.fname, delim = '\t', col_names = T)

pip <- as_tibble(data.table::fread(prior.fname, sep=' ', header=F)) # vroom cant parse this file...
colnames(pip) <- c('snp','torus_pip')
sumstats <- inner_join(sumstats, pip, by='snp')

# keep loci at 10% FDR
chunk.fdr <- read.delim(loci_qval, sep='', header=F, stringsAsFactors = F) # or this one..
chunks <- chunk.fdr$V2[chunk.fdr$V3 < 0.1]

# choose chunks that pass FDR
sumstats <- sumstats[sumstats$locus %in% chunks, ]

# ready for fine mapping
vroom::vroom_write(sumstats, 
                   delim = '\t', 
                   col_names = T, 
                   path = output)


