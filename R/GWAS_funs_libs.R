suppressMessages(library(tidyverse))
suppressMessages(library(bigsnpr))
suppressMessages(library(susieR))
suppressMessages(library(ggplot2))
suppressMessages(library(GenomicRanges))

make_ranges <- function(seqname, start, end){
  return(GRanges(seqnames = seqname, ranges = IRanges(start = start, end = end)))
}


locus.lookup <- function(locus){
  euro_ld <- read.delim('../ANNOTATIONS/Euro_LD_Chunks.bed',sep = '\t', stringsAsFactors = F, header = F)
  currLocus <- euro_ld[euro_ld$V4 == locus, ]
  chrom <- currLocus$V1
  start <- currLocus$V2
  end <- currLocus$V3
  return(list(locus=locus, chrom=chrom, start=start,end=end))
}

# Cleans summary statistics and normalizes them for downstream analysis
clean_sumstats <- function(sumstats, cols.to.keep){
  
  stopifnot(!is.null(sumstats))
  stopifnot(length(cols.to.keep) == 8)
  
  chr <- cols.to.keep[1]
  pos <- cols.to.keep[2]
  beta <- cols.to.keep[3]
  se <- cols.to.keep[4]
  a0 <- cols.to.keep[5]
  a1 <- cols.to.keep[6]
  rs <- cols.to.keep[7]
  pval <- cols.to.keep[8]
  
  # keep SNPs in 1kg
  #sumstats <- inner_join(sumstats, snps.to.keep, by=rs)
  # Extract relevant columns
  clean.sumstats <- sumstats[ ,c(chr, pos, beta, se, a0, a1, rs, pval)]
  colnames(clean.sumstats) <- c('chr','pos','beta','se','a0','a1','snp', 'pval')
  
  # drop XY chromosomes
  clean.sumstats <- clean.sumstats[!(clean.sumstats$chr %in% c("X","Y")), ]
  
  # make chromosomes integers
  clean.sumstats$chr <- as.integer(clean.sumstats$chr)
  
  # Compute Zscores
  zscore <- clean.sumstats$beta/clean.sumstats$se
  clean.sumstats['zscore'] <- zscore
  clean.sumstats <- clean.sumstats[!is.na(zscore),]
  
  # Keep SNPs only, no indels
  nucs <- c('A','C','T','G')
  bol <- (clean.sumstats$a0 %in% nucs) & (clean.sumstats$a1 %in% nucs)
  clean.sumstats <- clean.sumstats[bol,]
  
  # sort by chromosome and position
  clean.sumstats <- clean.sumstats[order(clean.sumstats$chr, clean.sumstats$pos), ]
  
  # drop duplicate SNPs
  chrpos <- paste0(clean.sumstats$chr, '_', clean.sumstats$pos)
  clean.sumstats <- clean.sumstats[!duplicated(chrpos), ]
  
  return(clean.sumstats)
}

# Assigns each SNP to one ld-block
assign.locus.snp <- function(cleaned.sumstats, ldBed){
  
  ld <- vroom::vroom(ldBed, col_names = F)
  ldRanges <- make_ranges(ld$X1, ld$X2, ld$X3)
  ldRanges <- plyranges::mutate(ldRanges, locus=ld$X4)
  
  snpRanges <- GRanges(seqnames = cleaned.sumstats$chr, 
                       ranges   = IRanges(start = cleaned.sumstats$pos, 
                                          end   = cleaned.sumstats$pos,
                                          names = cleaned.sumstats$snp))
  
  snpRanges <- plyranges::mutate(snpRanges, snp=names(snpRanges))
  
  snp.ld.overlap <- plyranges::join_overlap_inner(snpRanges, ldRanges)
  snp.ld.block <- as_tibble(snp.ld.overlap@elementMetadata)
  snp.ld.block <- snp.ld.block[!duplicated(snp.ld.block$snp), ] # some SNPs are in multiple ld-blocks due to edge of ld blocks
  cleaned.annot.sumstats <- inner_join(cleaned.sumstats, snp.ld.block, 'snp')
  
  return(cleaned.annot.sumstats)
}


# Each annotation gets assigned SNPs based on overlap
annotator <- function(gwas, annotations){
  
  snpRanges <- make_ranges(gwas$chr, gwas$pos, gwas$pos)
  snpRanges <- plyranges::mutate(snpRanges, snp=gwas$snp)
  
  for(f in annotations){
    
    name <- paste0(basename(f),'_d')
    curr <- rtracklayer::import(f, format='bed')
    subdf <- subsetByOverlaps(snpRanges, curr)
    snpsIn <- unique(subdf$snp)
    
    gwas <- mutate(gwas, !!name := ifelse(snp %in% snpsIn,1,0))
  }
  return(gwas)
}

# Annotations for causal SNPs (apply these after fine-mapping!)
annotator_merged <- function(gwas, annotations){
  
  snpRanges <- make_ranges(gwas$chr, gwas$pos, gwas$pos)
  snpRanges <- plyranges::mutate(snpRanges, snp=gwas$snp)
  gwas['annots'] <- ''
  
  for(f in annotations){
    
    curr <- rtracklayer::import(f, format='bed')
    subdf <- subsetByOverlaps(snpRanges, curr)
    snpsIn <- unique(subdf$snp)
    
    if(length(snpsIn)>0){
      curr <- gwas %>% pull(annots)
      curr <- curr[gwas$snp %in% snpsIn]
      delims <- rep(';', length(curr))
      delims[which(curr == '')] <- ''
      gwas[gwas$snp %in% snpsIn,"annots"] <- paste0(curr,delims,gsub(pattern = '.bed',replacement = '', x = basename(f)))
    }
  }
  return(gwas)
}

merge.bigsnp.gwas <- function(gwas, bigSNP){
  
  map <- bigSNP$map
  snp_info <- map[,c('chromosome','physical.pos','allele1','allele2')]
  colnames(snp_info) <- c('chr','pos','a0','a1')
  
  matched.gwas <- as_tibble(bigsnpr::snp_match(gwas, 
                                               snp_info, 
                                               strand_flip = T, 
                                               match.min.prop = 1)) %>% dplyr::rename(og_index = `_NUM_ID_.ss`) %>% dplyr::rename(bigSNP_index = `_NUM_ID_`) %>% mutate(zscore = beta/se)
  
  return(matched.gwas)
}

prune.regions <- function(sumstats, ref_panel){
  
  df_list <- list()
  r2_thresh <- 0.1
  for(l in unique(sumstats$locus)){
    sub.sumstats <- sumstats[sumstats$locus == l, ]
    
    topSnpPos <- which.max(sub.sumstats$susie_pip)
    topSnp <- sub.sumstats$bigSNP_index[topSnpPos]
    topSnpG <- ref_panel$genotypes[ ,topSnp]
    G <- ref_panel$genotypes[ , sub.sumstats$bigSNP_index]

    r2 <- as.vector(cor(topSnpG, G, method = 'pearson'))
    subRegions <- sub.sumstats[r2 > r2_thresh, ]
    subRegions$r2 <- r2[r2>r2_thresh]
    df_list[[l]] <- subRegions
  }
  return(Reduce(rbind, df_list))
}

# SUSIE related functions

run.susie <- function(sumstats, ref_panel, ldchunk, L, prior){
  
  sub.sumstats <- sumstats[sumstats$locus == ldchunk, ]
  if(nrow(sub.sumstats) > 1){
    X <- ref_panel$genotypes[ ,sub.sumstats$bigSNP_index]
    X <- scale(X, center = T, scale = T)
    zhat <- sub.sumstats$zscore
    R <- cov2cor((crossprod(X) + tcrossprod(zhat))/nrow(X))
    if(prior){
      res <- susie_rss(z = zhat, prior_weights = sub.sumstats$torus_pip, R = R, L = L)
    }
    else{
      res <- susie_rss(z = zhat, R = R, L = L)
    }
    return(res)
  }
}
# merges susie results with original summary statistics data frame
# ASSUMES L = 1! ONLY ONE CREDIBLE SET PER LOCUS! 
merge_susie_sumstats <- function(susie_results, sumstats){
  
  sumstats$susie_pip <- 0
  sumstats$CS <- 0
  loci <- names(susie_results)
  
  for(l in loci){
    n.snps <- length(susie_results[[l]]$pip)
    sumstats[sumstats$locus == as.numeric(l), "susie_pip"] <- susie_results[[l]]$pip
    
    snps.in.cs <- rep(0, n.snps)
    if(!is.null(susie_results[[l]]$sets$cs)){
      snps.in.cs[unlist(susie_results[[l]]$sets$cs$L1)] <- 1
    }
    sumstats[sumstats$locus == as.numeric(l), "CS"] <- snps.in.cs
  }
  return(sumstats)
}

# Annotations for causal SNPs (apply these after fine-mapping!)

add.gtex.annotation <- function(sumstats, gtex, tissue_type=''){
  stopifnot('chr' %in% colnames(gtex))
  stopifnot('pos' %in% colnames(gtex))
  stopifnot('name' %in% colnames(gtex))

  # load GTEx v7 significant snp-gene pairs
  gtex$var_id <- paste0(gtex$chr, '_', gtex$pos)
  gtex <- gtex[, c('var_id', 'name')]
  sumstats$var_id <- paste0(sumstats$chr, '_', sumstats$pos)
  
  # first left join, to get annotations of fine-mapped snps
  annot.sumstats <- left_join(x = sumstats, y = gtex, by='var_id')
  # collapse annotation into one
  annots <- annot.sumstats %>% group_by(var_id) %>% summarise(!!tissue_type := paste(unique(name), collapse = ';'))
  # join again to get original df with new annotation
  annot.sumstats <- inner_join(sumstats, annots, by = 'var_id') %>% dplyr::select(-var_id)
  
  return(annot.sumstats)
}


add.nearby.gene <- function(sumstats, gene.df, dist = 10000, gene_type){
  
  gene.ranges <- make_ranges(gene.df$chr, gene.df$start, gene.df$end)
  gene.ranges <- plyranges::mutate(gene.ranges, name=gene.df$gene)
  
  snp.ranges <- GRanges(seqnames = sumstats$chr, 
                               ranges = IRanges(start = sumstats$pos, 
                                                end   = sumstats$pos))
  snp.ranges <- plyranges::mutate(snp.ranges, snp = sumstats$snp)
  
  nearby.genes <- plyranges::join_overlap_inner(gene.ranges, snp.ranges, maxgap=dist) %>%
    as_tibble() %>%
    group_by(snp) %>%
    summarise(!!gene_type := paste(unique(name), collapse=';'))
  
  annot.sumstats <- left_join(sumstats, nearby.genes, by='snp')
  return(annot.sumstats)
}

add.consequence <- function(sumstats){
  
  stopifnot('snp' %in% colnames(sumstats))
  
  snp.list <- sumstats$snp
  snp.mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_SNP",
                               dataset = "hsapiens_snp")
  
  result <- suppressMessages(biomaRt::getBM(attributes = c("refsnp_id","consequence_type_tv"),
                                            filters    = "snp_filter", 
                                            values     = snp.list, 
                                            mart       = snp.mart) %>% as_tibble %>% group_by(refsnp_id) %>% summarise(consequence = paste0(consequence_type_tv, collapse = ';'))
  )
  
  colnames(result) <- c('snp', 'consequence')
  annot.sumstats <- left_join(sumstats, result, by='snp')
  return(annot.sumstats)
}


add.gwas_catalog <- function(sumstats, gwas){

  gwas <- gwas[ , c('SNPS','DISEASE/TRAIT')]
  colnames(gwas) <- c('snp','other_gwas_trait')

  annot.sumstats <- left_join(sumstats, gwas, by='snp') %>% 
    group_by(snp) %>% 
    summarise(gwas_catalog = paste0(unique(other_gwas_trait), collapse=';'))
  
  sumstats <- left_join(sumstats, annot.sumstats, by='snp')
  return(sumstats)
}

add.region <- function(sumstats){
  
  loci <- sumstats$locus
  L <- length(loci)
  region <- rep(0, L)
  for(l in 1:L){
    res <- locus.lookup(loci[l])
    region[l] <- paste0("chr",res$chr,":",round(res$start/10^6, 1),'-',round(res$end/10^6, 1),'M')
  }
  sumstats$region <- region
  return(sumstats)
}

# link SNPs to genes
compute_hic_gene_pip <- function(sumstats, hic_genes_file)
{
  snpRanges <- make_ranges(sumstats$chr, sumstats$pos, sumstats$pos)
  snpRanges <- plyranges::mutate(snpRanges, snp=sumstats$snp)

  enhancerRanges <- make_ranges(hic_genes_file$enhancer_chr, hic_genes_file$enhancer_start, hic_genes_file$enhancer_end)
  enhancerRanges <- plyranges::mutate(enhancerRanges, genes=hic_genes_file$gene)
  
  join_snps_distal <- plyranges::join_overlap_inner(enhancerRanges,snpRanges) %>% 
    as_tibble() %>% 
    select(c(snp, genes)) 
  #%>%  
   # group_by(snp) %>% 
   # summarise(hic_genes=paste0(unique(genes), collapse = ';'))
  
  annot.sumstats <- left_join(sumstats, join_snps_distal, on='snp') %>% group_by(genes) %>% summarise(Gene_PIP=sum(susie_pip)) %>% arrange(desc(Gene_PIP))

  return(annot.sumstats)
}


