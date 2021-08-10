# Gene Investigator
A simple application to interrogate the molecular evolution of a single gene

## Installation and dependencies
This application is currently designed to run in an HPC environment.

There is an assumption that the freely available Anaconda software is installed on your machine.

You will also need to download the standalone hyphy-analyses repository (https://github.com/veg/hyphy-analyses). Make sure to modify the config.yaml to point to the correct directory on your system

### To install -- Steps necessary to complete before running
1. `git clone https://github.com/aglucaci/GeneInvestigator.git`
2. `conda env create -f environment.yaml`. This will create a virtual environment called (GeneInvestigator) with the necessary dependencies.
3. At this point, run `conda activate GeneInvestigator` and your environment will be ready to go.

## Data retrival via NCBI Orthologs
For example, if we are interested in the TP53 gene: https://www.ncbi.nlm.nih.gov/gene/7157/ortholog/?scope=117570&term=TP53

Download all information: Tabular data, RefSeq Transcripts, and RefSeq Protein. (Typically one gene per species, but all transcripts per species is also available) 

## Results
The following are JSON files produced by HyPhy analyses. These can be visualized by the appropriate module from HyPhy Vision (http://vision.hyphy.org/). Analysis file names contain the method used (SLAC, FEL, PRIME, FADE, MEME, CFEL, etc), and if appropriate -- the set of branches to which the analysis was applied.

```
── results/TP53
│   ├── TP53.FEL.json
│   ├── TP53.FUBAR.json
│   ├── TP53.BUSTEDS.json
│   ├── TP53.MEME.json
│   ├── TP53.ABSREL.json
│   ├── TP53.SLAC.json
│   ├── TP53.BGM.json
│   ├── TP53.PRIME.json
│   ├── TP53.FMM.json
│   ├── TP53.ABSREL-MH.json
│   ├── TP53.BUSTEDS-MH.json
│   ├── TP53.MSS.json
│   ├── TP53.BUSTEDS-MSS.json
```


Alternative name: Analysis of Orthologous Collections (AOC)