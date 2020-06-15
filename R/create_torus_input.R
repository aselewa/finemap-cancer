source('R/GWAS_funs_libs.R')

args <- commandArgs(trailingOnly = TRUE)

## Start here with cleaned summary statistics
cleaned.gwas <- vroom::vroom(args[1], delim = '\t', col_names = T)

# Add annotations
annotations <- list.files(path = args[2], pattern = '*.bed', full.names = T)

cleaned.gwas.annots <- annotator(cleaned.gwas, annotations = annotations)

# vroom write seg faults..
data.table::fwrite(x = cleaned.gwas.annots[,-c(1:6,8:12)], file = args[3], quote = F, sep = '\t', col.names = T, row.names = F)
data.table::fwrite(x = cleaned.gwas[,c('snp','locus','zscore')], file = args[4], quote = F, sep = '\t', col.names = T, row.names = F)


