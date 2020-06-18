#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --output=snakelog.out
#SBATCH --time=36:00:00
#SBATCH --partition=xinhe
#SBATCH --mem=4G
#SBATCH --tasks-per-node=4

module load R/3.5.1
SNAKEMAKE="/project2/xinhe/software/miniconda3/bin/snakemake"

$SNAKEMAKE \
    -kp \
    --ri \
    -j 150 \
    --latency-wait 50 \
    --cluster-config cluster_config.json \
    -c "sbatch \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --tasks-per-node={cluster.tasks} \
        --partition={cluster.partition} \
        --job-name={cluster.name} \
        --output={cluster.logfile}" \
    -s torus_susie_snakemake.py
