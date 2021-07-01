#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 13:04:15 2020
@author: alexander g. lucaci
The idea for this is that:
    
    Given a protein sequence and a transcript sequence
    I find the codons by stepping over the transcript sequence until the translated sequence matches the protein sequence
    that way, I have only the codons and not the additional sequences from the transcript
    (Which may be useful later)
    I will also create two output files
        One with the STOP codon stripped (this makes it hyphy compatible.)
        One with the STOP codons (may be useful later, codon bias?)
"""

# =============================================================================
# Imports
# =============================================================================
from Bio import SeqIO
#from Bio.Alphabet import generic_rna
import sys
import argparse

# =============================================================================
# Declares
# =============================================================================
#PROTEIN_FASTA = sys.argv[1]
#TRANSCRIPTS_FASTA = sys.argv[2]
#OUTPUT = sys.argv[3]

PROTEIN_FASTA = snakemake.params.Prot
TRANSCRIPTS_FASTA = snakemake.params.Nuc
OUTPUT = snakemake.params.Out

results = []
no_match = []

#arguments = argparse.ArgumentParser(description='Combine alignments into a single file, adding a reference sequence as well')
#arguments.add_argument('-i', '--input',            help = 'FASTA file to process',                         required = True, type = str )
#arguments.add_argument('-o', '--output',           help = 'Directory for output and working files',        required = True, type = str)
#settings = arguments.parse_args()

# =============================================================================
# Helper functions
# =============================================================================
#turn into class
def check_against_protein_fasta(test_protein_seq, protein_file, species):
    # If the proteins are not the same length, exit
    pass

def already_in_results(transcript_desc):
    global results
    Found = False
    for record in results: # results stores transcript records that passed
        if transcript_desc == record.description: # already exists?
            #print("# Already found this, move on.")
            #start += 1
            #continue
            Found = True
            break
        #end if
    #end for
    return Found
#end method

def Process(protein_desc, protein_seq, TRANSCRIPTS_FASTA, species): #protein species
    global results, no_match
    #Loop over TRANSCRIPT_SEQ
    start = 0
    NT_SEQ_LENGTH = len(protein_seq) * 3
    # loop over all of the TRANSCRIPTS_Fasta seqs
    with open(TRANSCRIPTS_FASTA, "r") as transcript_handle:
        for m, transcript_record in enumerate(SeqIO.parse(transcript_handle, "fasta")):
            DONE = False
            # Grab Transcript Data
            transcript_id = transcript_record.id
            transcript_desc = transcript_record.description
            transcript_seq = transcript_record.seq
 
            if species not in transcript_desc: 
                #print("# Mismatch between species") # move on to the next one
                continue # only look at sequences from your species, not something similar.
            #end if
            #print("TX DESC:", transcript_desc)
            # can be a separate subroutine.
            start = 0
            NT_SEQ_LENGTH = len(protein_seq) * 3
            #print(NT_SEQ_LENGTH)
            while start < len(str(transcript_seq)):
                coding_dna = ""
                try:
                    coding_seq = transcript_seq[start: start + NT_SEQ_LENGTH]
                    #print(coding_seq)
                    coding_dna = coding_seq.translate() #translated, universal code
                except:
                    pass
                    #print("#ERROR", coding_dna)
                    #report ERROR
                #end try
                
                exists = already_in_results(transcript_desc)
                if coding_dna == str(protein_seq) and exists == False: # Exit upon first match, may be useful to see how many matches.
                  DONE = True
                  #print(codon_dna)
                  break
                else:
                  start += 1
                #end if
            #end while
            if DONE == True:  break
        #end for
    #end with
    if DONE == True:
        #return transcript_id, transcript_desc, coding_seq
        transcript_record.seq = coding_seq
        
        # If it fails, return a known
        return transcript_record
    else:
        return "NO_MATCH"
    #end if
#end method

# =============================================================================
# Main subroutine.
# =============================================================================
def progressBar(value, endvalue, bar_length=50):
    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()
#end method

def main(PROTEIN, TRANSCRIPTS): # Really to verify things.
    print("# TRANSCRIPT INPUT FILE:", TRANSCRIPTS)
    print("# PROTEIN INPUT FILE:", PROTEIN)
    protein_list = []
    transcript_list = []
    with open(TRANSCRIPTS, "r") as handle:
        #x = SeqIO.parse(handle, "fasta")
        #print(len(x))
        trans_count = 0 
        for record in SeqIO.parse(handle, "fasta"):
            trans_count +=1
            transcript_list.append(record.description)
        print("# Transcripts:", trans_count)    
    handle.close()
    
    with open(PROTEIN, "r") as handle:
        #x = SeqIO.parse(handle, "fasta")
        #print(len(x))
        prot_count = 0 
        for record in SeqIO.parse(handle, "fasta"):
            prot_count +=1
            protein_list.append(record.description)
        print("# Proteins:", prot_count)
    handle.close()
    
    #assert(trans_count == prot_count, "Counts do not match") # Check to make sure we have the same number of proteins and transcripts
    return trans_count, prot_count
    #for n, item in enumerate(transcript_list):
    #    print(item, [protein_list[n]])
#end method

# =============================================================================
# Main
# =============================================================================
#Verify files exist
print("# =============================================================================")
print("# Processing... ")
print("# =============================================================================")
trans_count, prot_count = main(PROTEIN_FASTA, TRANSCRIPTS_FASTA)
#Looks like species all match up in transcript and protein fasta.
#This is exceptional, will need to look for species name (from protein desc.) in transcript desc.

# DEBUG
#sys.exit(1)
      
# Create empty output file.
with open(OUTPUT, "w") as fh:
    fh.write("")
fh.close()

# Main program
successful_count = 0 
num_errors = 0
errors_IDs = []

# Grab the mrna transcript
# iterate over the possible proteins from it.
# Does any of them match a protein within the protein file? if so, have a dict store the transcript, protein pair.
# output the fasta
#{} = {"NT ID -- Short form i.e. the id": {"NT ID -- Full version i.e. the description": NT_ID, 
  #             "Codon Sequence": SEQ, 
  #             "Protein ID": PROTID, 
  #             "PROTEIN_SEQ": Protein_seq}}


# Grab the protein, can I find a match in any of the mRNA transcript? if so how many?
with open(PROTEIN_FASTA, "r") as prot_handle:
    for n, record in enumerate(SeqIO.parse(prot_handle, "fasta")):
        # Grab protein data.
        protein_id = record.id 
        protein_desc = record.description
        protein_seq = record.seq
        
        print("# Processing:", record.description)

        species = ""

        if "[" in record.description:
            species = record.description.split("[")[1].replace("]", "")
 
        #if "PREDICTED:" in record.description:
        #    species = " ".join(record.description.split(" ")[2:4])

        #if species == "":
        #    species = " ".join(record.description.split(" ")[1:3])
        print("# Species:", species)

        #print([species])
        #continue
        #print(n+1, protein_desc)
        #transcript_id, transcript_desc, coding_seq = Process (protein_desc, protein_seq, TRANSCRIPTS_FASTA)
        tx_record = Process (protein_desc, protein_seq, TRANSCRIPTS_FASTA, species)
        #results[tx_record.id] = {"TX_DESC": tx_record.desc, "CodonSequence": tx_record.seq}

        if type(tx_record) != str:
            print("# Match:", tx_record.description, "\n")
            results.append(tx_record)
        else:
            # No match...
            print("# -- NO Match -- \n")
            no_match.append(protein_desc)
            pass
    #end for
#end with


# Write out records
print("# Writing to:", OUTPUT)
SeqIO.write(results, OUTPUT, "fasta")

# Report on no matches
print("--- The following had no matches")
print("Total:", len(no_match))
for item in no_match:
    print(item)
#sys.exit(0)
# =============================================================================
# End of file    
# =============================================================================
