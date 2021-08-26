#!/bin/bash

set -euo pipefail

printf "Running snakemake...\n"

#snakemake --forceall --dag | dot -Tpdf > dag.pdf

mkdir -p logs

# Step 1 would be the alignment
# We just need to get an alignment, then take it for manual curation with RDP
# Also, let the Lineages stuff run here.


# This runs the non-lineage defined version of everything.
#snakemake \
#      -s Snakefile_Recombinants \
#      --cluster-config cluster.yaml \
#      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime=240:00:00 -e logs -o logs" \
#      --jobs 20 all \
#      --keep-going \
#      --reason \
#      --latency-wait 120 \
#      --use-conda

# Get lineages
#python3 scripts/FindLineages.py data/REV3L/REV3L_orthologs.csv REV3L
# Actually use the LineageAnnotation.py for this, pass in the treefile for recombinants, and csv, outputs the .clade files.

# Do analysis with lineages
#bash AssignLineages.sh

# take in the lineage file, and the recombinantion free file, do lineage assignment, RELAX and CFEL for now.
snakemake \
      -s Snakefile_Recombinants \
      --cluster-config cluster.yaml \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime=240:00:00 -e logs -o logs" \
      --jobs 20 all \
      --keep-going \
      --reason \
      --latency-wait 120 \
      --use-conda




exit 0
#End of file





 
