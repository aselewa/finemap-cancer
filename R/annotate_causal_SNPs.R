
source('R/GWAS_funs_libs.R')

### Annotate SNPs

trait = 'OV'
eqtl.files <- c('Breast_Mammary_var_genes.txt.gz','Colon_var_genes.txt.gz','Prostate_var_genes.txt.gz','Ovary_var_genes.txt.gz')
names(eqtl.files) <- c('BRCA','COAD','PRAD','OV')
curr.eqtl <- eqtl.files[trait]

# load SuSiE results and normalized summary statistics
sumstats <- vroom::vroom(paste0('susie_finemapping/',trait,'_susie_L1.txt.gz'), delim = '\t', col_names = T)
top.snps <- sumstats[sumstats$susie_pip > 0.05 & sumstats$CS == 1, ]

# Which annotation did each SNP come from?
annots <- c('H3K27ac.bed','H3K4me1.bed','H3K4me3.bed','H3K9ac.bed','ATAC.bed',paste0(trait,'_NON_immune.bed'))
annots <- paste0('torus_bed_files/', annots)
top.snps <- annotator_merged(gwas = top.snps, annots)
# keep immune only
top.snps <- top.snps[top.snps$annots != "", ]
non.immune <- paste0(trait,'_NON_immune')
top.snps <- top.snps[!grepl(pattern = non.immune, x = top.snps$annots), ]

# add HiC Links
hic_genes_file <- vroom::vroom('/project2/xinhe/alan/Cancer/GERMLINE/HiC/pcHiC_combined_immune.txt') %>% drop_na()
top.snps <- compute_hic_gene_pip(top.snps, hic_genes_file)
vroom::vroom_write(top.snps, path = paste0(trait,'_gene_pip.txt'))

# Which gene is each SNP inside of?
genes <- suppressMessages(vroom::vroom('/project2/xinhe/alan/Cancer/refGenome/b37/GENES_COORDs_b37.txt', col_names = T))
top.snps <- add.nearby.gene(top.snps, genes, dist=0, gene_type='inside_gene')

# Which gene's exon/UTR is it in?
exon_utr <- suppressMessages(vroom::vroom('/project2/xinhe/alan/Cancer/refGenome/b37/EXON_UTR_CORDS_b37.txt', col_names = T))
top.snps <- add.nearby.gene(top.snps, exon_utr, dist=0, gene_type='inside_exon_utr')

# Which promoter is each SNP in?
tss.df <- as_tibble(data.frame(chr=genes$chr, start=genes$start-2000, end=genes$start+1000, gene=genes$gene))
top.snps <- add.nearby.gene(top.snps, tss.df, dist = 0, gene_type='in_promoter')

# Which cancer driver is near each SNP?
drivers <- readLines('/project2/xinhe/alan/Cancer/GERMLINE/ANNOTATIONS/AllDrivers_FDR_0.10.txt')
top.snps <- add.nearby.gene(top.snps, genes[genes$gene %in% drivers, ], dist = 500000, gene_type='nearby_driver')

# Which immune gene is near each SNP?
immune.genes <- suppressMessages(vroom::vroom('/project2/xinhe/alan/Cancer/refGenome/b37/immune_gene_cords.txt', col_names = T))
top.snps <- add.nearby.gene(top.snps, immune.genes, dist = 500000, gene_type='nearby_immune')

# Which SNPs are eQTLs in breast and whole blood?
gtex <- suppressMessages(vroom::vroom(paste0('/project2/xinhe/alan/Cancer/GERMLINE/GTEx_Analysis_v7_eQTL/',curr.eqtl), col_names = T, delim = '\t'))
top.snps <- add.gtex.annotation(top.snps, gtex, tissue_type = 'cancer_gtex_egenes')
gtex <- suppressMessages(vroom::vroom('/project2/xinhe/alan/Cancer/GERMLINE/GTEx_Analysis_v7_eQTL/Whole_Blood_var_genes.txt.gz', col_names = T, delim = '\t'))
top.snps <- add.gtex.annotation(top.snps, gtex, tissue_type = 'whole_blood_gtex_eGenes')

# What trait is found for each SNP in GWAS Catalog?
gwas <- suppressMessages(vroom::vroom('/project2/xinhe/alan/Cancer/GERMLINE/GWAS_catalog/gwas_catalog_v1.0-associations_e98_r2020-02-08.tsv',col_names = T))
top.snps <- add.gwas_catalog(top.snps, gwas)

# What region is each SNP in?
top.snps <- add.region(top.snps)
top.snps$chrPos <- paste0('chr',top.snps$chr,':',top.snps$pos)

# do some formatting..
top.snps <- top.snps %>% replace(. == 'NA' | . == "", NA)
top.snps <- top.snps[order(top.snps$region), ]
top.snps$susie_pip <- round(top.snps$susie_pip, 2)
top.snps$zscore <- round(top.snps$zscore, 2)

region.pip <- top.snps %>% group_by(region) %>% summarise(Region_PIP=sum(susie_pip))
top.snps <- left_join(top.snps, region.pip, by="region")

# add HiC Links
hic_genes_file <- vroom::vroom('/project2/xinhe/alan/Cancer/GERMLINE/HiC/pcHiC_combined_immune.txt') %>% drop_na()
top.snps <- compute_hic_gene_pip(top.snps, hic_genes_file)

# select relevant columns
all.top.snps.out <- top.snps[ , c('region','snp','chrPos','inside_gene','inside_exon_utr','in_promoter','nearby_driver','nearby_immune','cancer_gtex_egenes','whole_blood_gtex_eGenes','annots','hic_genes','gwas_catalog','zscore','susie_pip','Region_PIP')]
colnames(all.top.snps.out) <- c('Region','SNP','Position','Gene Inside','Exon-UTR Inside','Promoter Inside','Nearby Drivers (500kb)','Nearby Immune Gene (500kb)',paste0(trait,' eGenes'),'Whole-Blood eGenes','Immune Mark','HiC Target','GWAS Catalog','ZScore','PIP','Region PIP')

# write results
vroom::vroom_write(all.top.snps.out, path = paste0('susie_finemapping/',trait,'_Immune_finemapped_L1_pip0.05_annotated.txt'), delim = '\t', col_names = T)


