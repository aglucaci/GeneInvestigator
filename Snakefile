"""``snakemake`` file that runs the Gene Investigator analysis.
Written by Alexander G Lucaci.
"""

import itertools
import os
import sys
import csv
import json
from snakemake.utils import min_version
#min_version('6.3.0')

#----------------------------------------------------------------------------
# Configuration
#----------------------------------------------------------------------------
configfile: 'config.yaml'

Nucleotide_file = config["Nucleotide"]
Protein_file = config["Protein"]
Label = config["Label"]
HYPHY = config["HyPhy"]
PREMSA = config["pre-msa"]
MAFFT = config["MAFFT"]
POSTMSA = config["post-msa"]
IQTREE = config["IQTREE"]

FMM = config["FMM"]
BUSTEDS_MH = config["BUSTEDSMH"]

#----------------------------------------------------------------------------
# Helper functions
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Rule all
#----------------------------------------------------------------------------
rule all:
    input:
        os.path.join("results", Label),
        os.path.join("results", Label + "_protein.fas"),
        os.path.join("results", Label + "_nuc.fas"),
        os.path.join("results", Label + "_protein.aln"),
        os.path.join("results", Label + "_codons.fasta"),
        os.path.join("results", Label + "_codons_duplicates.json"),
        os.path.join("results", Label + "_codons.fasta.treefile"),
        os.path.join("results", Label + "_codons.fasta.FEL.json"),
        os.path.join("results", Label + "_codons.fasta.FUBAR.json"),
        os.path.join("results", Label + "_codons.fasta.BUSTEDS.json"),
        os.path.join("results", Label + "_codons.fasta.MEME.json"),
        os.path.join("results", Label + "_codons.fasta.ABSREL.json"),
        os.path.join("results", Label + "_codons.fasta.SLAC.json"),
        os.path.join("results", Label + "_codons.fasta.BGM.json"),
        os.path.join("results", Label + "_codons.fasta.PRIME.json"),
        os.path.join("results", Label + "_codons.fasta.ABSREL-MH.json")


#----------------------------------------------------------------------------
# Rules
#----------------------------------------------------------------------------
rule get_codons:
    output:
        codons = "results/" + Label
    params:
        Nuc = Nucleotide_file,
        Prot = Protein_file,
        Out = os.path.join("results", Label)
    script:
        "scripts/codons.py"
#end rule get_codons

rule pre_msa:
    input: 
        codons = rules.get_codons.output.codons
    output: 
        protein_fas = os.path.join("results", Label + "_protein.fas"),
        nucleotide_fas = os.path.join("results", Label + "_nuc.fas")
    shell: 
        "{HYPHY} {PREMSA} --input {input.codons}"
#end rule pre_msa

rule mafft:
    input:
        protein = rules.pre_msa.output.protein_fas
    output:
        protein_aln = os.path.join("results", Label + "_protein.aln")
    shell:
        "{MAFFT} --auto {input.protein} > {output.protein_aln}"
#end rule mafft

rule post_msa:
    input: 
        protein_aln = rules.mafft.output.protein_aln,
        nucleotide_seqs = rules.pre_msa.output.nucleotide_fas      
    output: 
        codons_fas = os.path.join("results", Label + "_codons.fasta"),
        duplicates_json = os.path.join("results", Label + "_codons_duplicates.json")
    shell: 
        "{HYPHY} {POSTMSA} --protein-msa {input.protein_aln} --nucleotide-sequences {input.nucleotide_seqs} --output {output.codons_fas} --duplicates {output.duplicates_json}"
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
        tree = os.path.join("results", Label + "_codons.fasta.treefile")
    shell:
        "{IQTREE} -s {input.codons_fas}"
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
        results = os.path.join("results", Label + "_codons.fasta.FEL.json")
    shell: 
        "{HYPHY} FEL --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule FEL

rule FUBAR:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join("results", Label + "_codons.fasta.FUBAR.json")
    shell: 
        "{HYPHY} FUBAR --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule FUBAR

rule BUSTEDS:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join("results", Label + "_codons.fasta.BUSTEDS.json")
    shell: 
        "{HYPHY} BUSTED --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule BUSTEDS

rule MEME:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join("results", Label + "_codons.fasta.MEME.json")
    shell: 
        "{HYPHY} MEME --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule MEME

rule ABSREL:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join("results", Label + "_codons.fasta.ABSREL.json")
    shell: 
        "{HYPHY} ABSREL --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule ABSREL

rule SLAC:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join("results", Label + "_codons.fasta.SLAC.json")
    shell: 
        "{HYPHY} SLAC --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule 

rule BGM:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join("results", Label + "_codons.fasta.BGM.json")
    shell: 
        "{HYPHY} BGM --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule 

rule PRIME:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join("results", Label + "_codons.fasta.PRIME.json")
    shell: 
        "{HYPHY} PRIME --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --impute-states Yes"
#end rule 

rule ABSRELMH:
    input: 
        codon_aln = rules.post_msa.output.codons_fas,
        tree = rules.iqtree.output.tree      
    output: 
        results = os.path.join("results", Label + "_codons.fasta.ABSREL-MH.json")
    shell: 
        "{HYPHY} ABSREL --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --multiple-hits Double+Triple"
#end rule ABSREL

#----------------------------------------------------------------------------
# End of file
#----------------------------------------------------------------------------