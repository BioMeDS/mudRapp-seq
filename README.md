# mudRapp-seq

This repository contains the code that accompanies the paper introducing 
"*In-situ sequencing with multiple direct RNA-assisted padlock probing (mudRapp-seq)*":

> Ahmad S, Schaust J, Ambil UB, Fischer SC, Ankenbrand MJ, Smyth RP. *Localization of influenza viral RNA in cells by direct RNA padlock probing and in-situ sequencing.* (in preparation)

## Data

Raw data is archived independently in the [BioImage Archive (BIA)](https://www.ebi.ac.uk/bioimage-archive/) with accession number **TODO**.
In order to reproduce our analyses, download the raw data from BIA and put them into the `data/raw` folder.

## Analyses

### Python environment

To re-create the python environment with [`mamba`](https://github.com/mamba-org/mamba) run:

```bash
mamba env create -f environment.yml
mamba activate mudRapp-seq
```

### Data formatting

As [`starfish`](https://github.com/spacetx/starfish) is used, the raw data needs to be restructured in SpaceTx format.

In order to create the formatted data in `data/spacetx` run these steps in the root of the `mudRapp-seq` repo and in the `mudRapp-seq` environment:

```bash
python code/data_formatting/cDNA_vRNA.py
```