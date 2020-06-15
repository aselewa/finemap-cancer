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

### Step 2. Clean summary statistics

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

Example:
```
./preprocess.sh /project2/xinhe/alan/Cancer/GERMLINE/SUMMARY_STATISTICS/colorectal_raw_sumstats.txt.gz \
                chr,position_b37,ocac_overall_OR,ocac_overall_se,a0,a1,1KG_Phase3_ID,ocac_overall_pvalue \
                "ColCa"
```

This can take about 10 minutes depending on the size of your GWAS. If this worked, you should see a folder called `cleaned_sumstats` and inside there is `{prefix}_cleaned_sumstats.txt.gz` inside. The first few lines of the new file should look like:

```
zcat ColCa_cleaned_sumstats.txt.gz | head -n 5
chr     pos     a0      a1      beta    se      snp     pval    zscore  og_index        bigSNP_index    locus
1       13110   A       G       -0.02066        0.03231 1       0.52245 -0.6394305168678428     3       4       1
1       13116   G       T       -0.01255        0.01895 rs201725126     0.50772 -0.662269129287599      4       5       1
1       13118   G       A       -0.01255        0.01895 rs200579949     0.50772 -0.662269129287599      5       6       1
1       13550   A       G       0.13184 0.11236 1       0.24066 1.173371306514774       7       11      1
```

IMPORTANT: you can have multiple summary statistics in the `cleaned_sumstats/` folder. Just run the `preprocess.sh` recipe for each summary statistics trait. MAKE SURE THEY HAVE DIFFERENT PREFIX! 

### Step 3. Annotations

Annotations should be a bed-file with only the first 3 columns: chromosome, start, end. Chromosome should just be the number, no "chr" and should be in hg19/b37 coordinates. In your project directory, you should see a folder called `annotations/` Go into this. You should now see a folder called `{prefix}/` where `prefix` is what you choose in step 2. Place your annotation file (`.bed` extension in here). 

It may be useful to look at existing annotations that I placed at `/project2/xinhe/CANCER_GERMLINE/ANNOTATIONS/`

Your annotation should look like this:
```
head mybed.bed

1       564999  566395
1       568591  569959
1       569105  570623
1       713254  715133
1       762120  763780
```

IMPORTANT: you can place as many bed-files/annotations in the `prefix` directory as you want! Just make sure they have different names. These are the annotations that will be used for this specific trait.

### Step 4. Running finemapping

We are now ready to run the main pipeline for finemapping.

With your project directory (`my_project_directory`) from earlier, you should now see `cleaned_sumstats` and `annotations/`. 
To run the pipeline, execute the following:

```
sbatch submit_snake.sbatch
```

Make sure the job is running:

```
qstat -u CNETID
```

Once it is over, you should see a directory called `susie_finemapping` which contains the finemapped results for each chromosome. 
