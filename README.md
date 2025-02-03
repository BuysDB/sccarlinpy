

This repository contains code to perform single cell lineage tracing preprocessing steps.
## Overview of steps of the analysis
- Demultiplexing of reads
- Adapter trimming
- Alignment using Minimap2 to filter off-target reads 
- Smith waterman alignment to determine scar 
- Merging data from multiple runs [optional]
- Filtering of doublet cells
- Export to table (cell to scar)


## Execution
Install the package using pip. 
Copy the snakefile to the directory containing the FastQ files.
Edit config.json to set the reference parameters.
Run the snakefile using snakemake.
The output of the snakemake is a *.scstats.tsv file which contains the scar sequences found in each cell.
Run the notebook `1_Merge_scars.ipynb` to merge the scar sequences from multiple runs.
Then run notebook `2_Filter scars and flag doublets.ipynb` to filter scar seq. (This requires cell-meta data from Seurat)
Finally run `3_Expor.ipynb` to export all data to one file called `full_dataset_calls.csv` which can be joined to the cell meta in Seurat.
