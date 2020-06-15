
module load R/3.5.1

SUMSTATS="/project2/xinhe/alan/Cancer/GERMLINE/SUMMARY_STATISTICS"
bigsnp_1kg="/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds"
ldBlocks_bed="/project2/xinhe/alan/Cancer/GERMLINE/ANNOTATIONS/Euro_LD_Chunks.bed"

mkdir -p torus_ready_sumstats

# colorectal
Rscript R/clean_sumstats.R \ 
        "$SUMSTATS"/colorectal_raw_sumstats.txt.gz \
        chromosome,base_pair_location,beta,standard_error,OTHER_ALLELE,EFFECT_ALLELE,variant_id,p_value \
        $bigsnp_1kg \
        $ldBlocks_bed \
        torus_ready_sumstats/COAD_sumstats.txt.gz

# ovarian
Rscript R/clean_sumstats.R \
        "$SUMSTATS"/ovca_overall_raw_sumstats.txt.gz \
        chr,position_b37,ocac_overall_OR,ocac_overall_se,a0,a1,1KG_Phase3_ID,ocac_overall_pvalue \
        $bigsnp_1kg \
        $ldBlocks_bed \
        torus_ready_sumstats/OV_sumstats.txt.gz

# prostate
Rscript R/clean_sumstats.R \
        "$SUMSTATS"/prad_raw_sumstats.txt.gz \
        $bigsnp_1kg \
        $ldBlocks_bed \
        chr,Position,OR,SE,Baseline,Effect,SNP,Pnorm \
        torus_ready_sumstats/PRAD_sumstats.txt.gz
        

# breast
Rscript R/clean_sumstats.R \
        "$SUMSTATS"/brca_raw_onco2_sumstats.txt.gz \
        chr,position_b37,bcac_onco2_beta,bcac_onco2_se,a0,a1,phase3_1kg_id,bcac_onco2_P1df_Wald \
        $bigsnp_1kg \
        $ldBlocks_bed \        
        torus_ready_sumstats/BRCA_sumstats.txt.gz
        
# breast er-pos
Rscript R/clean_sumstats.R \
"$SUMSTATS"/brca_raw_sumstats.txt.gz \
chr,position_b37,bcac_onco2_erpos_beta,bcac_onco2_erpos_se,a0,a1,phase3_1kg_id,bcac_onco2_erpos_P1df_Wald \
$bigsnp_1kg \
$ldBlocks_bed \
torus_ready_sumstats/BRCA_erpos_cleaned_sumstats.txt.gz

# breast er-neg
Rscript R/clean_sumstats.R \
"$SUMSTATS"/brca_raw_sumstats.txt.gz \
chr,position_b37,bcac_onco2_erneg_beta,bcac_onco2_erneg_se,a0,a1,phase3_1kg_id,bcac_onco2_erneg_P1df_Wald \
$bigsnp_1kg \
$ldBlocks_bed \
torus_ready_sumstats/BRCA_erneg_cleaned_sumstats.txt.gz


