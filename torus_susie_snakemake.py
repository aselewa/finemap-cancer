# A snakemake protocol for performing enrichment analysis and functionally-informed finemapping
# Uses PolyFun

import os
import glob
import pandas as pd

logdir = "log/"
if not os.path.isdir(logdir):
  os.mkdir(logdir)

# REQUIRED
pd = 'PROJECT DIRECTORY HERE'

# Inputs (summary stats, a directory of annotations) REQUIRED
cleaned_sumstats = pd + 'torus_ready_sumstats/'
bed_dir = pd + 'torus_bed_files/'

# Outputs created (DO NOT EDIT)
torus_output = pd + 'torus_output/'
annotations = pd + 'torus_annotations/'
finemapping = pd + 'susie_finemapping/'

# software/metadata (DO NOT EDIT)
torus = '/project2/xinhe/software/dap/torus_src/torus'
bigsnp_1kg = '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds'

# Global wildcards
traits = glob_wildcards(cleaned_sumstats + "{traits}_cleaned_sumstats.txt.gz").traits
chrom = list(map(str, range(1,23)))

localrules: concat_finemapping_results, concat_unif_finemapping_results

rule all:
    input:
        expand(annotations + '{t}_torus_annotation.txt.gz', t=traits),
        expand(annotations + '{t}_torus_sumstats.txt.gz', t=traits),
        expand(torus_output + '{t}/{t}_enrichment.txt', t=traits),
        expand(torus_output + '{t}/{t}_SNP_PIP.txt', t=traits),
        expand(torus_output + '{t}/{t}_loci_qval.txt', t=traits),
        expand(finemapping + '{t}/{t}_chr{c}_susie_L1.txt.gz',t=traits,c=chrom),
        expand(finemapping + '{t}/{t}_chr{c}_susie_L1_UNIFORM.txt.gz',t=traits,c=chrom)
        
rule create_torus_input:
    input:
        cleaned_sumstats + '{traits}_cleaned_sumstats.txt.gz',
        bed_dir + '{traits}'
    output:
        annotations + '{traits}_torus_annotation.txt.gz',
        annotations + '{traits}_torus_sumstats.txt.gz'
    shell:
        "Rscript R/create_torus_input.R {input[0]} {input[1]} {output[0]} {output[1]}"
        
rule run_torus:
     input:
        annotations + '{traits}_torus_sumstats.txt.gz',
        annotations + '{traits}_torus_annotation.txt.gz'
     output:
        torus_output + '{traits}/{traits}_enrichment.txt',
        torus_output + '{traits}/{traits}_SNP_PIP.txt'
     shell:
        "{torus} -d {input[0]} -annot {input[1]} --load_zval -dump_prior {wildcards.traits}_curr_prior > {output[0]} && cat {wildcards.traits}_curr_prior/* > {output[1]} && rm -rf {wildcards.traits}_curr_prior"
        
        
rule run_torus_qtl_mode:
     input:
        annotations + '{traits}_torus_sumstats.txt.gz',
        annotations + '{traits}_torus_annotation.txt.gz'
     output:
        torus_output + '{traits}/{traits}_loci_qval.txt'
     shell:
        "{torus} -d {input[0]} -annot {input[1]} --load_zval -qtl > {output}"
          
        
rule merge_sumstats_torus:
     input:
         cleaned_sumstats + '{traits}_cleaned_sumstats.txt.gz',
         torus_output + '{traits}/{traits}_SNP_PIP.txt',
         torus_output + '{traits}/{traits}_loci_qval.txt'
     output:
         finemapping + '{traits}_SuSiE_ready_sumstats.txt.gz'
     shell:
         "Rscript R/make_susie_sumstats.R {input[0]} {input[1]} {input[2]} {output}"

rule run_finemapping:
      input:
          finemapping + '{traits}_SuSiE_ready_sumstats.txt.gz',
          bigsnp_1kg
      output:
          finemapping + '{traits}/{traits}_chr{chrom}_susie_L1.txt.gz'
      params:
          output_prefix= finemapping + '{traits}/{traits}_chr{chrom}'
      shell:
          "Rscript R/susie_workflow.R {input[0]} {input[1]} {params.output_prefix} torus {wildcards.chrom}"
          
rule run_finemapping_uniform:
      input:
          finemapping + '{traits}_SuSiE_ready_sumstats.txt.gz',
          bigsnp_1kg
      output:
          finemapping + '{traits}/{traits}_chr{chrom}_susie_L1_UNIFORM.txt.gz'
      params:
          output_prefix = finemapping + '{traits}/{traits}_chr{chrom}'
      shell:
          "Rscript R/susie_workflow.R {input[0]} {input[1]} {params.output_prefix} uniform {wildcards.chrom}"
          
