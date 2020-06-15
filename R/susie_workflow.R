
# Command line script for running SuSiE
# Run this on RCC

source('R/GWAS_funs_libs.R')

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

# main
args <- commandArgs(trailingOnly = T)
sumstats <- args[1]
bigSNP <- args[2]
out.prefix <- args[3]
priortype <- args[4]
chrom <- args[5]

susie.df <- vroom::vroom(sumstats,delim = '\t', col_names = T)
susie.df <- susie.df[susie.df$chr == chrom, ]

bigsnp.1kg <- snp_attach(rdsfile = bigSNP)
chunks <- unique(susie.df$locus)

if(length(chunks) == 0){
  if(priortype=="torus"){
    system(paste0('touch ',out.prefix,'_susie_L1.txt.gz'))
  }
  else if(priortype== "uniform"){
    system(paste0('touch ',out.prefix,'_susie_L1_UNIFORM.txt.gz'))
  }
} else{
  if(priortype == "torus"){
    susie_res <- list()
    for(i in 1:length(chunks)){
      print(paste0(i,' out of ', length(chunks)))
      susie_res[[as.character(chunks[i])]] <- run.susie(susie.df, bigsnp.1kg, chunks[i], L = 1, prior = T)
    }
    
    susie.df.merged <- merge_susie_sumstats(susie_res, susie.df)
    fname <- paste0(out.prefix,'_susie_L1.txt.gz')
    
  } else if(priortype == "uniform"){
    susie_res <- list()
    for(i in 1:length(chunks)){
      print(paste0(i,' out of ', length(chunks)))
      susie_res[[as.character(chunks[i])]] <- run.susie(susie.df, bigsnp.1kg, chunks[i], L = 1, prior = F)
    }
    susie.df.merged <- merge_susie_sumstats(susie_res, susie.df)
    fname <- paste0(out.prefix,'_susie_L1_UNIFORM.txt.gz')
  } else{
    stop('prior type not recognized')
  }
  vroom::vroom_write(susie.df.merged, path = fname, col_names = T)
}






