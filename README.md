## Introduction

### Step 0. Setting up

Most software is already installed on our directory. The only things you will need to install are the R packages.

On midway2 do:
```
module load R/3.5.1

install.packages('tidyverse')
install.packages('bigsnpr')
install.packages("vroom")
devtools::install_github("stephenslab/susieR@0.9.0")
BiocManager::install(c("GenomicRanges", "plyranges","rtracklayer"))
```

**Step 1.** 

Obtain summary statistics and place them here
`/project2/xinhe/CANCER_GERMLINE/SUMMARY_STATISTICS/`

You will need summary statistics with the following 8 columns:
* chromosome (just the number, no "chr", hg19!)
* position (base pair position, hg19!)
* beta (if you have Odds Ratio, you will need to transform it to log(Odds Ratio))
* standard error (SE)
* reference allele (A,C,T,G)
* association/effect allele (A,C,T,G)
* some kind of SNP ID (rsID, or chr:position:a1:a2)
* p-value

Write down the column names somewhere for you will need them next.

**Step 2.** 

Munge/clean summary statistics

`R/clean_sumstats
