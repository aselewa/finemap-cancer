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

### Step 1. Summary statistics 

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

### Step 2. Pre-processing

Create a working directory and clone this repo into it
```
mkdir my_project_directory
cd my_project_directory
git clone git@github.com:aselewa/finemap-cancer.git .
```

Run the preprocessing step with THREE arguments:
```
./preprocess.sh SUMSTATS COLUMNS prefix
```
* SUMSTATS : full path to your summary statistics
* COLUMNS: a comma(,) delimited string with the columns we described above! Make sure they are in the order described above.
* prefix: a meaningful name for the cancer/data type. i.e. BrCa for breast cancer

If this worked, you should see a folder called `torus_ready_sumstats` and inside there is `prefix_sumstats.txt.gz`

```

```

