"""``snakemake`` file that runs the Gene Investigator analysis.
Written by Alexander G Lucaci.
"""

import itertools
import os
import sys
import csv
import json
from pathlib import Path
from snakemake.utils import min_version
#min_version('6.3.0')

#----------------------------------------------------------------------------
# Configuration
#----------------------------------------------------------------------------
configfile: 'config.yaml'
Nucleotide_file = config["Nucleotide"]
Protein_file = config["Protein"]
Label = config["Label"]
#HYPHY = config["HyPhy"]
PREMSA = config["pre-msa"]
#MAFFT = config["MAFFT"]
POSTMSA = config["post-msa"]
#IQTREE = config["IQTREE"]
FMM = config["FMM"]
BUSTEDS_MH = config["BUSTEDSMH"]
MSS = config["MSS"]
BUSTEDMSS = config["BUSTEDMSS"]
CODONSTSV = config["CODONSTSV"]

# Set output directory
BASEDIR = os.getcwd()
OUTDIR = os.path.join(BASEDIR, "results/" + Label)

# Create output dir.
Path(OUTDIR).mkdir(parents=True, exist_ok=True)

#----------------------------------------------------------------------------
# Helper functions
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Rule all
#----------------------------------------------------------------------------
rule all:
    input:
        os.path.join(OUTDIR, Label),
        os.path.join(OUTDIR, Label + "_protein.fas"),
        os.path.join(OUTDIR, Label + "_nuc.fas"),
        os.path.join(OUTDIR, Label + "_protein.aln"),
        os.path.join(OUTDIR, Label + "_codons.fasta"),
        os.path.join(OUTDIR, Label + "_codons_duplicates.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.treefile"),
        os.path.join(OUTDIR, Label + "_codons.fasta.FEL.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.FUBAR.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.BUSTEDS.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.MEME.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.ABSREL.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.SLAC.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.BGM.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.PRIME.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.ABSREL-MH.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.BUSTEDS-MH.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.MSS.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.BUSTED-MSS.json"),
        os.path.join(OUTDIR, Label + "_codons.fasta.FMM.json")
#end rule all

#----------------------------------------------------------------------------
# Rules
#----------------------------------------------------------------------------
rule get_codons:
    output:
        codons = os.path.join(OUTDIR, Label)
    params:
        Nuc = Nucleotide_file,
        Prot = Protein_file,
        Out = os.path.join(OUTDIR, Label)
    conda: 'environment.yaml'
    script:
        "scripts/codons.py"
#end rule get_codons

rule pre_msa:
    input: 
        codons = rules.get_codons.output.codons
    output: 
        protein_fas = os.path.join(OUTDIR, Label + "_protein.fas"),
        nucleotide_fas = os.path.join(OUTDIR, Label + "_nuc.fas")
    conda: 'environment.yaml'
    shell: 
        "hyphy {PREMSA} --input {input.codons}"
#end rule pre_msa

rule mafft:
    input:
        protein = rules.pre_msa.output.protein_fas
    output:
        protein_aln = os.path.join(OUTDIR, Label + "_protein.aln")
    conda: 'environment.yaml'
    shell:
        "mafft --auto {input.protein} > {output.protein_aln}"
#end rule mafft

rule post_msa:
    input: 
        protein_aln = rules.mafft.output.protein_aln,
        nucleotide_seqs = rules.pre_msa.output.nucleotide_fas      
    output: 
        codons_fas = os.path.join(OUTDIR, Label + "_codons.fasta"),
        duplicates_json = os.path.join(OUTDIR, Label + "_codons_duplicates.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy {POSTMSA} --protein-msa {input.protein_aln} --nucleotide-sequences {input.nucleotide_seqs} --output {output.codons_fas} --duplicates {output.duplicates_json}"
#end rule pre_msa

#----------------------------------------------------------------------------
# AlignmentProfiler
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# TN93
#----------------------------------------------------------------------------

rule iqtree:
    input:
        codons_fas = rules.post_msa.output.codons_fas
    output:
        tree = os.path.join(OUTDIR, Label + "_codons.fasta.treefile")
    conda: 'environment.yaml'
    shell:
        "iqtree -s {input.codons_fas}"
#end rule iqtree

#----------------------------------------------------------------------------
# FADE, need to root on something.
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Annotate tree for taxonomy.
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Selection Analyses
#----------------------------------------------------------------------------

rule FEL:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.FEL.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy FEL --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule FEL

rule FUBAR:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.FUBAR.json")
    conda: 'environment.yaml'    
    shell: 
        "hyphy FUBAR --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule FUBAR

rule BUSTEDS:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.BUSTEDS.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy BUSTED --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule BUSTEDS

rule MEME:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.MEME.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy MEME --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule MEME

rule ABSREL:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.ABSREL.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy ABSREL --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule ABSREL

rule SLAC:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.SLAC.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy SLAC --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule 

rule BGM:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.BGM.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy BGM --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule 

rule PRIME:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.PRIME.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy PRIME --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --impute-states Yes"
#end rule 

rule ABSRELMH:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.ABSREL-MH.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy ABSREL --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --multiple-hits Double+Triple"
#end rule ABSRELMH

rule BUSTEDSMH:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.BUSTEDS-MH.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy {BUSTEDS_MH} --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule BUSTEDSMH

rule MSS:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.MSS.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy {MSS} --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --neutral NEUTRAL --classes {CODONSTSV} --type global --frequencies CF3x4 --ci Yes --lrt Yes"
#end rule MSS

rule BUSTEDMSS:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.BUSTED-MSS.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy {BUSTEDMSS} --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --classes {CODONSTSV} --neutral NEUTRAL"
#end rule BUSTEDMSS

rule FMM:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.fasta.FMM.json")
    conda: 'environment.yaml'
    shell: 
        "hyphy {FMM} --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --triple-islands Yes"
#end rule FMM


#----------------------------------------------------------------------------
# End of file
#----------------------------------------------------------------------------
