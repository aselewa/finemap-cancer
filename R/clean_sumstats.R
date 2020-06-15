
source('R/GWAS_funs_libs.R')

print('Cleaning summary statistics...')

args <- commandArgs(trailingOnly=T)
sumstats <- vroom::vroom(args[1], col_names = T)
cols.to.keep <- unlist(strsplit(as.character(args[2]), split=','))
bigsnp.1kg <- args[3]
ldBlocks_bed <- args[4]
outName <- args[5]

cleaned_sumstats <- clean_sumstats(sumstats, cols.to.keep)

# match to reference panel
bigsnp.1kg <- snp_attach(rdsfile = bigsnp.1kg)
cleaned_sumstats <- merge.bigsnp.gwas(cleaned_sumstats, bigsnp.1kg)

# assign each snp to an LD block
cleaned_sumstats <- assign.locus.snp(cleaned.sumstats = cleaned_sumstats, ldBed = ldBlocks_bed)

print('Writing output..')
vroom::vroom_write(cleaned_sumstats, outName)

print('Done!')
