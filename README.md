# EGFR Inhibition Promotes Enteroendocrine Cell Differentiation 

<p align="left">
  <img src="docs/eec_subtypes_v3.png" alt="D10 Lapatinib EEC subtypes" width="600">
</p>


## What's this?

This repository is a copy of the code used in the single cell analyses performed in the 'EGFR INHIBITION PROMOTES ENTEROENDOCRINE CELL DIFFERENTIATION CONTRIBUTING TO TREATMENT-ASSOCIATED DIARRHEA' (Biorxiv 2026). 

The code is provided here in the interest of transparency and reproducibility, and displays all steps for data processing downstream of CellRanger v10.0 called on the raw FASTQ, which can be accessed on GEO.

## Downloading the data.

The final GEO submission will be available at GSE329493. The data are available as raw FASTQ files, as well as an intermediate matrix, and metadata sheet. Separate libraries were prepared from duodenal organoids sequenced at Day 2, 4, and 10 of the breault lab differentiation protocol, in both a control and a Lapatinib treated condition for each, excluding D10 (Lapatinib only).

## Running the code.

The code can be run with the conda environment specified in conda/, and also packaged as a binary targeting MacOS ARM architecture. Each notebook folder corresponds with a dataset that can run as a sequential set of notebooks, drawing on the src for standardised steps between each dataset.
