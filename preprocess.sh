
module load R/3.5.1

SUMSTATS=$1
COLUMNS=$2
PREFIX=$3
bigsnp_1kg="/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds"
ldBlocks_bed="/project2/xinhe/CANCER_GERMLINE/ANNOTATIONS/Euro_LD_Chunks.bed"

mkdir -p cleaned_sumstats

Rscript R/clean_sumstats.R $SUMSTATS \
                           $COLUMNS \
                           $bigsnp_1kg \
                           $ldBlocks_bed \
                           cleaned_sumstats/"$PREFIX"_cleaned_sumstats.txt.gz

# set-up corresponding annotation directory
mkdir -p annotations/"$PREFIX"
