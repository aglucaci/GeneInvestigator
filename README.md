# Gene Investigator
A simple application to interrogate the molecular evolution of a single gene

Installation and dependencies
This application is currently designed to run in an HPC environment.

There is an assumption that the freely available Anaconda software is installed on your machine.

You will also need to download the standalone hyphy-analyses repository (https://github.com/veg/hyphy-analyses). Make sure to modify the config.yaml to point to the correct directory on your system

To install -- Steps necessary to complete before running
git clone https://github.com/aglucaci/GeneInvestigator.git
conda env create -f environment.yaml. This will create a virtual environment called (GeneInvestigator) with the necessary dependencies.
At this point, run conda activate GeneInvestigator and your environment will be ready to go.


## Data retrival via NCBI Orthologs
For example, if we are interested in the TP53 gene: https://www.ncbi.nlm.nih.gov/gene/7157/ortholog/?scope=117570&term=TP53

Download all information: Tabular data, RefSeq Transcripts, and RefSeq Protein. (Typically one gene per species, but all transcripts per species is also available) 
